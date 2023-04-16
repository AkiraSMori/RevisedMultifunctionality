#### Results are stored as CSV files ####
essentials <- c("tidyverse", "reshape2", "broom.mixed", "data.table", "multifunc",
                "glmmTMB", "effects", "picante", "vegan", "RColorBrewer")
lapply(essentials, require, character.only=TRUE)
cols1<-rev(brewer.pal(9, "Spectral"))
cols2<-rev(brewer.pal(11, "Spectral"))
cols3<-rev(brewer.pal(10, "Spectral"))
#####################################
#### Modifications for functions ####
#####################################
getFuncMaxed<-function (adf, vars=NA, thresh=0.05, proportion=F, prepend="Diversity", maxN=1)
{
  if(is.na(vars)[1])
    stop("You need to specify some response variable names")
  vars<-whichVars(adf, vars)
  getMaxValue<-function(x){
    l<-length(x)
    mean(sort(x, na.last=F)[l:(l-maxN+1)], na.rm=T)
  }
  funcMaxed<-rowSums(colwise(function(x) x >= thresh *
                               getMaxValue(x))(adf[,which(names(adf) %in% vars)]), na.rm=T)
  ret<-data.frame(cbind(adf[, which(names(adf) %in% prepend)], funcMaxed))
  names(ret)<-c(names(adf)[which(names(adf) %in% prepend)],
                "funcMaxed")
  ret$nFunc<-rowSums((adf[,which(names(adf) %in% vars)]!="NA")==TRUE, na.rm=T)
  if(proportion)
    ret$funcMaxed<-ret$funcMaxed/ret$nFunc
  ret
} #na.rm=T allows to have different number of functions considered.
## For output summary from glmmTMB ##
getCoefTab2<-function(eqn, fun=glmmTMB, data, groupVar="thresholds", coefVar=NULL,
                      ...)
{
  ret<-plyr::ddply(data, .variables=groupVar, function(adf)
  {
    options(warn = 0)
    aFit<-try(fun(eqn, data=adf, offset=log(nFunc))) #offset term
    options(warn = 0)
    if("try-error" %in% class(aFit))
      return(rep(NA, 4))
    if("glm" %in% class(aFit)) {
      if (!aFit$converged)
        return(rep(NA, 4))
    }
    coefInfo<-summary(aFit)$coef
    coefInfo<-cbind(coef_name = rownames(coefInfo$cond), as.data.frame(coefInfo$cond))
    return(coefInfo)
  })
  if(!is.null(coefVar))
    ret<-subset(ret, ret$coef_name==coefVar)
  return(ret)
}
## MF species importace function ##
SpImportSES<-function(dat=dat, vars=vars, prepend=c("plot","Diversity"),
                      threshmin=0.05, threshmax=0.95, maxN=length(vars),
                      gamma=length(splist), splist=splist,
                      Nsim=999)
{
  fitThresh<-getFuncsMaxed(dat, vars, threshmin=0.05, threshmax=0.95,
                           prepend=c("plot","Diversity"), maxN=maxN)
  datComm<-dat %>% dplyr::select(plot, all_of(splist))
  ## Difference metric ##
  datCommPivot<-datComm %>% tidyr::pivot_longer(col=-plot, names_to="Species",
                                                values_to="Presence")
  fitThreshComm<-full_join(fitThresh, datCommPivot, by="plot")
  resDiffMt<-fitThreshComm %>%
    group_by(thresholds, Species, Presence) %>%
    dplyr::summarise(meanMF=mean(funcMaxed), .groups = "drop") %>%
    group_by(thresholds, Species) %>%
    dplyr::summarise(DiffMt=meanMF[Presence==1]-meanMF[Presence==0], .groups = "drop")
  #Presence - Absence
  ## Simulation ##
  obsPivot<-fitThreshComm %>% group_by(thresholds, Species)
  Nplot<-length(unique(obsPivot$plot))
  randMatrix<-matrix(nrow=nrow(resDiffMt), ncol=Nsim)
  for(n in 1:Nsim)
  {
    randPivot<-obsPivot %>% group_by(thresholds, Species) %>% slice_sample(., n=Nplot,
                                                                           replace=FALSE) %>%
      select(thresholds, Species, funcMaxed)
    obsPivot$funcMaxed<-randPivot$funcMaxed
    tmpRand<-obsPivot %>%
      group_by(thresholds, Species, Presence) %>%
      dplyr::summarise(meanMF=mean(funcMaxed), .groups = "drop") %>%
      group_by(thresholds, Species) %>%
      dplyr::summarise(DiffMt=meanMF[Presence==1]-meanMF[Presence==0], .groups = "drop")
    #Presence - Absence
    randMatrix[,n]<-tmpRand$DiffMt
  }
  resDiffMt$meanNullDiff<-apply(randMatrix, 1, mean)
  resDiffMt$sdNullDiff<-apply(randMatrix, 1, sd)
  resDiffMt <- resDiffMt %>% mutate(., SESDiffMt=((DiffMt-meanNullDiff)/sdNullDiff) %>%
                                      unlist()) %>%
    as.data.frame()
  return(resDiffMt)
}
## Randomize MF function ##
ThreshSES<-function(dat, vars, pred="Diversity", Nsim=999, threshmin=0.05, threshmax=0.95,
                    maxN=1)
{
  res<-getFuncsMaxed(dat, vars, prepend=pred, threshmin=threshmin, threshmax=threshmax,
                     maxN=maxN)
  null<-matrix(0,nrow=nrow(res),ncol=Nsim)
  for(i in 1:Nsim)
  {
    use<-dat[,vars]
    tmp<-as.data.frame(apply(use,2,function(x)as.numeric(x[sample(1:length(x))])))
    null[,i]<-getFuncsMaxed(tmp, vars, threshmin=threshmin,
                            threshmax=threshmax,maxN=maxN)$funcMaxed
  }
  res$meanNull<-apply(null,1,mean)
  res$sdNull<-apply(null,1,sd)
  res$sesFuncMaxed<-(res$funcMaxed-res$meanNull)/res$sdNull
  return(res)
}
#####################################
#####################################
###########################################################
#### Biodepth Germnay data from Byrnes et al. 2014 MEE ####
###########################################################
#### Read in data to run sample analyses on the biodepth data ####
data(all_biodepth)
allVars<-qw(biomassY3, root3, N.g.m2, light3, N.Soil, wood3, cotton3)
varIdx<-which(names(all_biodepth) %in% allVars)
Dat<-subset(all_biodepth, all_biodepth$location=="Germany")
DatPlot<-Dat[,c("plot","Diversity")]
UsedVars<-whichVars(Dat, allVars)
DatMF<-Dat[,UsedVars]
DatComm<-Dat[,relevantSp(Dat,26:ncol(Dat))]
gamma<-ncol(DatComm)
Freq<-colSums(DatComm)
Spp<-names(Freq) ## Species<-relevantSp(Dat,26:ncol(Dat))
Dat<-cbind(DatPlot, DatMF, DatComm) #Updated Dat
###########################################################
#####################################
######## Species importance ########
#####################################
#### Calculations for all possible combinations of functions ####
AllCombSESMatrix<-rep(NA, 8)
for (k in 2:length(UsedVars))
{
  Comb<-combn(length(UsedVars),k)
  CombNo<-ncol(Comb) #Number of combinations
  CombSESMatrix<-rep(NA, 8)
  for(l in 1:CombNo)
  {
    SubVars<-UsedVars[Comb[,l]]
    OmitVars<-UsedVars[-Comb[,l]] #Vars that will be omitted
    if(length(SubVars)<length(UsedVars))
    {
      TempDat<-Dat[,-which(colnames(Dat) %in% OmitVars)] #Remove OmitVars
    } else {
      TempDat<-Dat #In case, no OmitVars (i.e., all functions)
    }
    Vars<-whichVars(TempDat, SubVars)
    SESMatrix<-SpImportSES(dat=TempDat, vars=Vars,
                           threshmin=0.05, threshmax=0.95, prepend=c("plot","Diversity"),
                           maxN=length(Vars), gamma=gamma, Nsim=999, splist=Spp)
    SESMatrix[,c("combID", "nFunc")]<-cbind.data.frame(l,k)
    names(CombSESMatrix)<-names(SESMatrix)
    CombSESMatrix<-dplyr::bind_rows(CombSESMatrix, SESMatrix)
  }
  names(AllCombSESMatrix)<-names(CombSESMatrix)
  AllCombSESMatrix<-dplyr::bind_rows(AllCombSESMatrix, CombSESMatrix)
}
AllCombSESMatrix<-dplyr::filter(AllCombSESMatrix, !is.na(thresholds))
AllCombSESMatrix<-AllCombSESMatrix %>% mutate(UniqCombID = group_indices(., nFunc,
                                                                         combID)) #Assign unique id
fwrite(AllCombSESMatrix, "ResSpImporSES.csv")
AllCombSESMatrix<-fread("ResSpImporSES.csv")
#### Show results for all individual combinations and means of each species at a given number of functions ####
ResSum<-AllCombSESMatrix %>%
  group_by(nFunc, thresholds, Species) %>%
  dplyr::summarise(mean(SESDiffMt))
colnames(ResSum)<-c("nFunc", "thresholds", "Species", "SES")
fig.1<-ggplot() +
  geom_boxplot(data=AllCombSESMatrix, aes(y=SESDiffMt, x=thresholds*100, group=thresholds),
               col="light grey") +
  geom_boxplot(data=ResSum, aes(x=thresholds*100, y=SES, group=thresholds), col="black") +
  geom_hline(yintercept=c(-1.96,1.96), linetype="dashed") +
  xlab("Threshold (%)") + ylab("Standardized effect size (SES)") +
  facet_wrap(~nFunc,ncol=length(UsedVars)-1) +
  theme_bw(base_size=14)
windows()
fig.1
SumSig<-ResSum %>% dplyr::filter(SES>=1.96) %>%
  group_by(nFunc, thresholds) %>% dplyr::summarise(n()) %>% data.frame()
colnames(SumSig)<-c("nFunc", "thresholds", "no.spp")
SumSig$direction<-c("positive")
SumNonSig<-ResSum %>% dplyr::filter(SES<=-1.96) %>%
  group_by(nFunc, thresholds) %>% dplyr::summarise(n()*-1) %>% data.frame()
colnames(SumNonSig)<-c("nFunc", "thresholds", "no.spp")
SumNonSig$direction<-c("negative")
ResSumSig<-bind_rows(SumSig, SumNonSig)
fig.2<-ggplot() +
  geom_area(data=filter(ResSumSig, direction=="positive"), aes(x=thresholds*100, y=no.spp,
                                                               group=-nFunc, fill=nFunc), col="white") +
  geom_area(data=filter(ResSumSig, direction=="negative"), aes(x=thresholds*100, y=no.spp,
                                                               group=-nFunc, fill=nFunc), col="white") +
  xlim(0,100) +
  scale_fill_gradientn(colours = cols1[2:5]) +
  xlab("Threshold (%)") + ylab("Number of species")
windows()
fig.2
#####################################
#####################################
##################################
######## Offset approach ########
##################################
#### Calculations for all possible combinations of functions ####
AllCombFitThresh<-rep(NA, 6)
AllCombFitLinearSlopes<-rep(NA, 8)
for (k in 2:length(UsedVars))
{
  Comb<-combn(length(UsedVars),k)
  CombNo<-ncol(Comb) #Number of combinations
  CombFitLinearSlopes<-rep(NA, 8)
  for(l in 1:CombNo)
  {
    SubVars<-UsedVars[Comb[,l]]
    OmitVars<-UsedVars[-Comb[,l]] #Vars that will be omitted
    if(length(SubVars)<length(UsedVars))
    {
      TempDat<-Dat[,-which(colnames(Dat) %in% OmitVars)] #Remove OmitVars
    } else {
      TempDat<-Dat #In case, no OmitVars (i.e., all functions)
    }
    Vars<-whichVars(TempDat, SubVars)
    FitThreshTemp<-getFuncsMaxed(TempDat, Vars, threshmin=0.05, threshmax=0.95,
                                 proportion=F,
                                 prepend=c("plot","Diversity"), maxN=length(SubVars))
    SampFitLinearSlopes<-NA
    SampFitLinearSlopes<-try(getCoefTab(funcMaxed~Diversity+offset(log(nFunc)),
                                        data=FitThreshTemp,
                                        groupVar="thresholds", coefVar="Diversity",
                                        fun=glm, family="poisson"), silent=F)
    #offset: This allows some observations for functions are missing.
    #The offset term is not neccesary for the Germany data because no function data is lacking.
    #The function "getCoefTab" can be thus used.
    ifelse(class(SampFitLinearSlopes)=="try-error", SampFitLinearSlopes,
           SampFitLinearSlopes[,c("combID", "nFunc")]<-cbind.data.frame(l,k))
    FitThreshTemp <- FitThreshTemp %>%
      select(thresholds, plot, Diversity, funcMaxed, nFunc)
    FitThreshTemp$combID<-l
    names(AllCombFitThresh)<-names(FitThreshTemp)
    AllCombFitThresh<-dplyr::bind_rows(AllCombFitThresh, FitThreshTemp)
    CombFitLinearSlopes<-rbind(CombFitLinearSlopes, SampFitLinearSlopes)
  }
  names(AllCombFitLinearSlopes)<-names(CombFitLinearSlopes)
  AllCombFitLinearSlopes<-dplyr::bind_rows(AllCombFitLinearSlopes, CombFitLinearSlopes)
}
AllCombFitLinearSlopes<-dplyr::filter(AllCombFitLinearSlopes, !is.na(thresholds))
AllCombFitThresh<-dplyr::filter(AllCombFitThresh, !is.na(thresholds))
AllCombFitThresh<-AllCombFitThresh %>% mutate(UniqCombID = group_by(., nFunc, combID))
#Assign unique id
fwrite(AllCombFitLinearSlopes, "AllCombFitLinearSlopes.csv")
AllCombSESMatrix<-fread("AllCombFitLinearSlopes.csv")
fwrite(AllCombFitThresh, "AllCombFitThresh.csv")
AllCombFitThresh<-fread("AllCombFitThresh.csv")
#### Show results for all individual combinations and means at a given toal number of functions considered ####
ResSum2<-AllCombFitLinearSlopes %>% group_by(nFunc, thresholds) %>%
  dplyr::summarise(mean(estimate))
colnames(ResSum2)<-c("nFunc", "thresholds", "estimate")
## Y-axis converted to proportional changes (%) ##
fig.3<-ggplot() +
  geom_smooth(data=AllCombFitLinearSlopes,
              aes(x=thresholds*100, y=100*estimate, group=combID), se=FALSE, col="light grey", size=0.5) +
  geom_abline(intercept=0, slope=0, lwd=0.5, linetype=2) +
  geom_smooth(data=ResSum2,
              aes(x=thresholds*100, y=100*estimate, group=nFunc), se=FALSE, col="dark blue",
              size=1) +
  ylab("Change in number of functions per addition of 1 species (%)") + xlab("Threshold (%)") +
  facet_wrap(~nFunc,ncol=length(UsedVars)-1) +
  theme_bw(base_size=14)
windows()
fig.3
#### Get common estimates by accounting for differences in total number of functions (offset slopes; Fig. 4) ####
AllCombFitThresh$UniqCombID<-as.factor(AllCombFitThresh$UniqCombID)
PlotX<-data.frame(
  Diversity=rep(seq(1, max(Dat$Diversity),length=max(Dat$Diversity)), 4),
  nFunc=c(rep(2,length=max(Dat$Diversity)), rep(3,length=max(Dat$Diversity)),
          rep(4,length=max(Dat$Diversity)), rep(5,length=max(Dat$Diversity))),
  UniqCombID=NA
)
threshmin<-5
threshmax<-95
AllCombOffsetSlopes<-rep(NA, 10)
AllCombOffsetFit<-rep(NA, 5)
for(i in threshmin:threshmax) {
  mod.offset<-try(glmmTMB(data=dplyr::filter(AllCombFitThresh, thresholds==i/100),
                          funcMaxed~Diversity+(1|UniqCombID), offset=log(nFunc),
                          #nFunc(offset) / UniqCombID(random factor)
                          family="poisson"), silent=F)
  if(class(mod.offset)!="try-error")
  {
    CommonEstimate<-cbind(tidy(mod.offset, effects="fixed", conf.int = TRUE)[2,],
                          thresholds=i/100)
    names(AllCombOffsetSlopes)<-names(CommonEstimate)
    AllCombOffsetSlopes<-dplyr::bind_rows(AllCombOffsetSlopes, CommonEstimate)
    PredY<-predict(mod.offset, newdata=PlotX, re.form=NULL, interval = c("confidence"))
    TempRes<-data.frame(PlotX, PredY, thresholds=i/100)
  }
  names(AllCombOffsetFit)<-names(TempRes)
  AllCombOffsetFit<-dplyr::bind_rows(AllCombOffsetFit, TempRes)
  #"PredY" is neccesary to show different slopes for each number of functions; Fig. 4a).
  #Function "getCoefTab2 (for glmmTMB)" can be used if "PredY" (estimated based on "mod.offset") is not required.
}
AllCombOffsetSlopes<-dplyr::filter(AllCombOffsetSlopes, !is.na(thresholds))
AllCombOffsetFit<-dplyr::filter(AllCombOffsetFit, !is.na(thresholds))
fwrite(AllCombOffsetSlopes, "AllCombOffsetSlopes.csv")
AllCombOffsetSlopes<-fread("AllCombOffsetSlopes.csv")
fwrite(AllCombOffsetFit, "AllCombOffsetFit.csv")
AllCombOffsetFit<-fread("AllCombOffsetFit.csv")
## Slopes for every 10% of thresholds ##
AllCombOffsetFitSubset<-AllCombOffsetFit %>%
  dplyr::filter(thresholds==0.1|thresholds==0.2|thresholds==0.3|thresholds==0.4|
                  thresholds==0.5|thresholds==0.6|thresholds==0.7|thresholds==0.8|thresholds==0.9)
AllCombOffsetFitSubset$percent<-paste(100*AllCombOffsetFitSubset$thresholds, "%", sep="")
fig.4a<-ggplot(AllCombOffsetFitSubset,aes(x=Diversity, y=PredY)) +
  geom_line(aes(color=as.factor(nFunc)), size=0.5) +
  scale_color_manual(values=cols2[3:9]) +
  ylab(expression("Number of functions >= Threshold")) + xlab("Species richness") +
  facet_wrap(.~percent, ncol=9) +
  theme_bw(base_size=14)
windows()
fig.4a
## Y-axis converted to proporetional changes (%) ##
fig.4d<-ggplot(AllCombOffsetSlopes, aes(x=thresholds)) +
  geom_ribbon(fill="light blue", alpha=0.5, aes(x=thresholds*100, ymin=100*conf.low,
                                                ymax=100*conf.high)) +
  geom_point(aes(x=thresholds*100, y=100*estimate)) +
  geom_abline(intercept=0, slope=0, lwd=0.5, linetype=2) +
  ylab("Change in number of functions per addition of 1 species (%)") + xlab("Threshold (%)") +
  theme_bw(base_size=14)
windows()
fig.4d
##################################
##################################
#######################################
######## Randomizing approach ########
#######################################
DatThresh<-getFuncsMaxed(Dat, UsedVars, threshmin=0.05, threshmax=0.95,
                         prepend=c("plot","Diversity"), maxN=length(UsedVars))
gcPlot<-subset(DatThresh, DatThresh$thresholds %in% qw(0.25, 0.5, 0.75)) #Using qw as %in% is a string comparison operator
gcPlot$percent<-paste(100*gcPlot$thresholds, "%", sep="")
qplot(Diversity, funcMaxed, data=gcPlot, facets=~percent, geom="jitter") +
  stat_smooth(method="glm", method.args=list(family=quasipoisson(link="identity")),
              colour="red", lwd=1.2) +
  ylab(expression("Number of functions >= Threshold")) +
  xlab("Diversity") + theme_bw(base_size=14)
StandFit<-ThreshSES(dat=Dat, UsedVars, pred=c("Diversity"), Nsim=999, threshmin=0.05,
                    threshmax=0.95,
                    maxN=length(UsedVars))
StandSlopes<-getCoefTab(sesFuncMaxed~Diversity, data=StandFit, coefVar="Diversity", fun=lm)
fwrite(StandSlopes, "StandSlopes.csv")
StandSlopes<-fread("StandSlopes.csv")
fig.5a.inset<-ggplot(StandSlopes, aes(x=thresholds)) +
  geom_ribbon(aes(x=thresholds*100, ymin=estimate-1.96*.data[["std.error"]],
                  ymax=estimate+1.96*.data[["std.error"]]), fill="light grey", alpha=0.5)+
  geom_point(aes(x=thresholds*100, y=estimate)) +
  geom_hline(yintercept=0, lwd=1) +
  geom_vline(xintercept=StandSlopes$thresholds[StandSlopes$estimate==max(StandSlopes$estimate)]*100, linetype="dashed") +
  ylab("Change in Standardized effect size (SES) for increase of 1 species") + xlab("Threshold (%)") +
  theme_bw(base_size=14)
windows()
fig.5a.inset
StandFit$percent<-100*StandFit$thresholds
dplyr::inner_join(StandFit, StandSlopes, by="thresholds") -> StandFit2
StandFit2$significnace<-StandFit2[,"p.value"]<=0.05
fig.5a<-ggplot(data=StandFit2, aes(x=Diversity, y=sesFuncMaxed, group=percent)) +
  stat_smooth(method="lm", lwd=0.5, fill=NA, aes(color=percent, linetype=significnace)) +
  scale_color_gradientn(name="Percent of ¥nMaximum", colours=cols3) +
  scale_linetype_manual(values=c("dotted", "solid")) +
  ylab(expression("Standardized effect size (SES)")) + xlab("Species richness") +
  theme_bw(base_size=14)
windows()
fig.5a
#######################################
#######################################