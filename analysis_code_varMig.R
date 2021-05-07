setwd("~/Dropbox/varmodel2/results")
library(RSQLite)
library(reshape2)
library("vegan")
library(tidyverse)
library(igraph)
library(cowplot)
library(scales)
library(adegenet)
library("factoextra")

prefix<-"Aug19/AugNew"
allParaTable<- paste( prefix, "_param.csv",sep="")
sumTable<- paste(prefix, "_AllSumTab.txt",sep="")
netSameSize<-paste(prefix, "_allNetCal.txt",sep="")
PTSquantile<-paste(prefix, "_allPTS.txt",sep="")
DURquantile<-paste(prefix, "_allDurTab.txt",sep="")

#import all the parameters varied for each scenario
allPara<-read.table(allParaTable,sep=",",header=T)
allPara<-as_tibble(allPara)
allPara<-allPara%>%select(NO, BITING_RATE_MEAN,N_GENES_INITIAL)

#import summary stats file
sumTableNames<-c("time","pop_id","n_infections","n_infected","n_infected_bites","n_total_bites","n_circulating_strains","n_circulating_genes","num","selMode","IRS","r","EIR","Prevalence","MOI","pool_size","n_circulating_alleles.x","n_circulating_alleles.y")

summaryTable<-read.table(sumTable,sep="\t",header=F)
colnames(summaryTable)<-sumTableNames
summaryTable<-as_tibble(summaryTable)

summaryTable<-left_join(summaryTable, allPara,by=c("num"="NO"))

#import PTS quantile file
ptsTableNames<-c("time","q25","q50","q75","q100","q125","q150","q175","q200","q225","q250","q275","q300","q325","q350","q375","q400","q425","q450","q475","q500","q525","q550","q575","q600","q625","q650","q675","q700","q725","q750","q775","q800","q825","q850","q875","q900","q925","q950","q975","q1000", "num","selMode","IRS","r")

ptsTable<-read.table(PTSquantile,sep="\t",header=F)
colnames(ptsTable)<-ptsTableNames
ptsTable<-as_tibble(ptsTable)
summaryTable<-summaryTable%>%
  left_join(ptsTable%>%select(time,q50,num, selMode, IRS, r), by=c("time", "num", "selMode", "IRS", "r"))


###Fig 1, Intervention--------
# import biting Rates
foi<-read.csv("~/Dropbox/varmodel2/ms/Figures/bitingRates.csv", header=T)
foiBr<-rbind(foi,foi,foi)
foiBr<-foiBr%>%mutate(mbr = rep(c(0.0001,0.00025,0.0005),each = 1080))
dat_rect <- data.frame(
  xmin = 361/360,
  xmax = 720/360,
  ymin = -Inf,
  ymax = Inf
)
ggplot()+geom_rect(data=dat_rect, aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),fill="pink", alpha=0.5)+
  geom_line(data = foiBr, aes(Day/360, Br*mbr, lty =as.factor(mbr) ))+
  annotate("text", x = c(100/360, 550/360,820/360), y = c(1.25, 1.25,1.25), 
           label = c("pre-IRS", "during IRS","post-IRS"), size = 4,
           fontface = "bold")+
  ylab("Biting rate per host per day")+xlab("Year")+
  scale_linetype(name = "",labels = c("low","medium","high"), guide = guide_legend(direction="horizontal"))+
  theme_light()+
  theme(legend.position = c(0.5,0.92),
        legend.background = element_rect(colour = 'black', fill = 'white', linetype='solid'))
ggsave("~/Dropbox/varmodel2/ms/Figures/Fig1C.pdf", width=8, height=4)

#Fig 2 extinction prob+prevalence after IRS------

extStats<-summaryTable%>%filter(IRS>0, time==39960)%>%group_by(num, IRS,BITING_RATE_MEAN, N_GENES_INITIAL, selMode)%>%summarise(count = n(), extTotal = sum(n_infections == 0), extProb = extTotal/count, sdExtProb =sd(n_infections == 0)/count)

bitesTotal<-summaryTable%>%filter(IRS==0, selMode =="S")%>%group_by(num)%>%summarise(count = n(), bites = floor(sum(n_total_bites)/count*12/10000))

transCat<-data.frame(bites = c(44, 110, 221), transmission = c("low","medium","high"))
transCat$transmission <- factor(transCat$transmission, levels= c("low","medium","high"))


extStats<-extStats%>%left_join(bitesTotal, by = "num")%>%left_join(transCat, by = "bites")


extStats$transmission <- factor(extStats$transmission, levels= c("low","medium","high"))

p1<-ggplot(extStats,aes(N_GENES_INITIAL,extProb, col = as.factor(selMode)))+geom_point(alpha=0.4)+geom_line(lty=3)+facet_grid(transmission~IRS,labeller = label_both)+ 
  scale_colour_manual(name = "Selection",values = c("blue","red")) + theme_light()+ylab("Extinction probability")+xlab("Gene pool size")+
  geom_errorbar(aes(ymin=extProb-sdExtProb, ymax=extProb+sdExtProb), width=.2,position=position_dodge(0.05))
p1

p1t<-ggplot(extStats%>%filter(transmission=="high"),aes(N_GENES_INITIAL,extProb, col = as.factor(selMode)))+geom_point(alpha=0.4)+geom_line(lty=3)+facet_grid(transmission~IRS,labeller = label_both)+ 
  scale_colour_manual(name = "Selection",values = c("blue","red")) + theme_light()+ylab("Extinction probability")+xlab("Gene pool size")+
  geom_errorbar(aes(ymin=extProb-sdExtProb, ymax=extProb+sdExtProb), width=.2,position=position_dodge(0.05))
p1t

extCount<-summaryTable%>%filter(IRS>0, time==39960)%>%group_by(num, IRS, selMode, r)%>%summarise(ext = (n_infections == 0))

summaryTable<-summaryTable%>%left_join(bitesTotal, by = "num")
summaryTable<-summaryTable%>%left_join(extCount, by = c("num", "IRS", "selMode", "r"))
summaryTable<-summaryTable%>%left_join(transCat, by = "bites")

p2<-summaryTable%>%filter(time > 28800,time <= 28800+(IRS)*360) %>% group_by(num, IRS,BITING_RATE_MEAN, N_GENES_INITIAL, selMode,transmission,r)%>%summarise(mPrev = mean(Prevalence))%>%
  group_by(num, IRS,BITING_RATE_MEAN, N_GENES_INITIAL, selMode,transmission)%>%summarize(count = n(),mPrevalence = mean(mPrev), prevSd = sd(mPrev))%>%  
  ggplot(aes(N_GENES_INITIAL,mPrevalence, col = as.factor(selMode)))+geom_point(size = 2,alpha = 0.2)+geom_line(lty=3)+
  geom_errorbar(aes(ymin=mPrevalence-prevSd, ymax=mPrevalence+prevSd), width=.2,position=position_dodge(0.05))+
  facet_grid(transmission~IRS, scale = "free_y" ,labeller = label_both)+ 
  scale_colour_manual(name = "Selection",values = c("blue","red")) + theme_light()+ylab("Prevalence")+xlab("Gene pool size")
p2

p2t<-summaryTable%>%filter(time > 28800,time <= 28800+(IRS)*360) %>% group_by(num, IRS,BITING_RATE_MEAN, N_GENES_INITIAL, selMode,transmission,r)%>%summarise(mPrev = mean(Prevalence))%>%
  group_by(num, IRS,BITING_RATE_MEAN, N_GENES_INITIAL, selMode,transmission)%>%summarize(count = n(),mPrevalence = mean(mPrev), prevSd = sd(mPrev))%>%
  filter(transmission=="high")%>%
  ggplot(aes(N_GENES_INITIAL,mPrevalence, col = as.factor(selMode)))+geom_point(size = 2,alpha = 0.2)+geom_line(lty=3)+
  geom_errorbar(aes(ymin=mPrevalence-prevSd, ymax=mPrevalence+prevSd), width=.2,position=position_dodge(0.05))+
  facet_grid(transmission~IRS, scale = "free_y" ,labeller = label_both)+ 
  scale_colour_manual(name = "Selection",values = c("blue","red")) + theme_light()+ylab("Prevalence")+xlab("Gene pool size")
p2t

plot2by2<-plot_grid(p1,p2, labels = c("(A)","(B)"), ncol =1)
save_plot(paste("~/Dropbox/varmodel2/ms/Figures/","SupplFig2_ext.pdf", sep=""), plot2by2,
          ncol = 1, # we're saving a grid plot of 2 columns
          nrow = 2, # and 2 rows
          # each individual subplot should have an aspect ratio of 1.3
          base_aspect_ratio = 1.8,
          base_height = 4.5
)

save_plot(paste("~/Dropbox/varmodel2/ms/Figures/","SupplFig2_ext.png", sep=""), plot2by2,
          ncol = 1, # we're saving a grid plot of 2 columns
          nrow = 2, # and 2 rows
          # each individual subplot should have an aspect ratio of 1.3
          base_aspect_ratio = 1.8,
          base_height = 4.5
)


plot2by2<-plot_grid(p1t,p2t, labels = c("(A)","(B)"), ncol =1)
save_plot(paste("~/Dropbox/varmodel2/ms/Figures/","Fig2_ext.pdf", sep=""), plot2by2,
          ncol = 1, # we're saving a grid plot of 2 columns
          nrow = 2, # and 2 rows
          # each individual subplot should have an aspect ratio of 1.3
          base_aspect_ratio = 2.5,
          base_height=2.5
)

save_plot(paste("~/Dropbox/varmodel2/ms/Figures/","Fig2_ext.png", sep=""), plot2by2,
          ncol = 1, # we're saving a grid plot of 2 columns
          nrow = 2, # and 2 rows
          # each individual subplot should have an aspect ratio of 1.3
          base_aspect_ratio = 2.5,
          base_height=2.5
)

###Fig S1------
colnames(summaryTable)[20]<-"genePool"

p3<-summaryTable%>%filter(time > 28800+(IRS+2)*360, ext==FALSE) %>% group_by(num, IRS,BITING_RATE_MEAN, genePool, selMode,transmission,r)%>%summarise(mPrev = mean(Prevalence))%>%
  group_by(num, IRS,BITING_RATE_MEAN, genePool, selMode,transmission)%>%summarise(count = n(),mPrevalence = mean(mPrev), prevSd = sd(mPrev))%>%  
  ggplot(aes(genePool,mPrevalence, col = as.factor(selMode)))+geom_point(size = 2,alpha = 0.2)+geom_line(lty=3)+
  geom_errorbar(aes(ymin=mPrevalence-prevSd, ymax=mPrevalence+prevSd), width=.2,position=position_dodge(0.05))+
  facet_grid(transmission~IRS, scale = "free_y" ,labeller = label_both)+ 
  scale_colour_manual(name = "Selection",values = c("blue","red")) + theme_light()+ylab("Prevalance")+xlab("Gene pool size")

ggsave(paste("~/Dropbox/varmodel2/ms/Figures/","FigS1_prevAfterIRS.pdf", sep=""),p3, width = 8.1, height = 4.5)


###Fig 3 DOI+PTS--------
allDur<-read.table(DURquantile, header=F, sep="\t")
colnames(allDur)<-c("year", "q50","q75","q80","q85","q90","q95","mdur","sdDur", "num","selMode","IRS","r")
allDur<-allDur%>%left_join(allPara, by = c("num"="NO"))
allDur<-allDur%>%left_join(bitesTotal, by = "num")
allDur<-allDur%>%left_join(transCat, by = "bites")

allDur<-allDur%>%mutate(duringIRS = ((year>80)*(year<=80+IRS)))
allDur$yearWhole<-as.factor(allDur$year)

summaryTable<-summaryTable%>%mutate(year = ceiling(time/360))

summaryTable<-summaryTable%>%mutate(yearWhole = ceiling(time/360))
summaryTable$yearWhole<-as.factor(summaryTable$yearWhole)
summaryTable<-summaryTable%>%
  mutate(yearCut = cut_width(year, width=1, boundary=0) )
df <- data.frame(
  yearWhole = rep(unique(summaryTable$yearWhole)[71],3),
  upper = unique(summaryTable$yearWhole)[c(73,76,81)],
  y = c(0,0,0),
  IRS = c(2,5,10)
)

p4<-ggplot()+ geom_rect(data = df, aes(xmin = yearWhole, xmax = upper, ymin = y, ymax = y + 400), fill="pink", alpha=0.5)+geom_boxplot(data=allDur%>%filter(year>=79, year<=99, transmission =="high"), aes(yearWhole, q95, fill = selMode),outlier.shape = NA, lwd=0.2)+
  facet_grid(.~IRS,labeller = label_both)+
  scale_fill_manual(name = "Selection",values = c("blue","red"))+
  xlab("Year")+ylab("95th percentile value of\nduration of infection")+ theme_cowplot()+scale_x_discrete(labels=c("-2","","0","","2","","","5","","","","","10","","","","","15","","",""))

p5<-ggplot()+ geom_rect(data = df, aes(xmin = yearWhole, xmax = upper, ymin = y, ymax = y + 0.2), fill="pink", alpha=0.5)+
  geom_boxplot(data = summaryTable%>%filter( year>=79, year<=99,transmission=="high",ext==FALSE, genePool %in% c(14400)), aes(yearWhole,q50, fill = selMode), outlier.shape = NA, lwd=0.2)+
  facet_grid(genePool~IRS,labeller = label_both, scale = "free_y")+
  scale_fill_manual(name = "Selection",values = c("blue","red"))+
  xlab("Year")+ylab("5th percentile of PTS")+ theme_cowplot()+scale_x_discrete(labels=c("-2","","0","","2","","","5","","","","","10","","","","","15","","",""))


plot2by2<-plot_grid(p4,p5, labels = c("(A)","(B)"), ncol =1, rel_heights=c(1,1))
save_plot(paste("~/Dropbox/varmodel2/ms/Figures/","Fig3_DOI_PTS.pdf", sep=""), plot2by2,
          ncol = 1, # we're saving a grid plot of 2 columns
          nrow = 2, # and 2 rows
          # each individual subplot should have an aspect ratio of 1.3
          base_aspect_ratio = 2.2,
          base_height = 4
)

save_plot(paste("~/Dropbox/varmodel2/ms/Figures/","Fig3_DOI_PTS.png", sep=""), plot2by2,
          ncol = 1, # we're saving a grid plot of 2 columns
          nrow = 2, # and 2 rows
          # each individual subplot should have an aspect ratio of 1.3
          base_aspect_ratio = 2.2,
          base_height = 4
)

###Fig S2+S3-----
p6<-ggplot()+geom_rect(data = df, aes(xmin = yearWhole, xmax = upper, ymin = y, ymax = Inf), fill="pink", alpha=0.5)+geom_boxplot(data = allDur%>%filter( year>=79, year<=99), aes(as.factor(year), q95, fill = selMode),outlier.shape = NA, lwd=0.2)+
  facet_grid(transmission~IRS,labeller = label_both)+
  scale_fill_manual(name = "Selection",values = c("blue","red"))+
  xlab("Year")+ylab("95% quantile value of\nduration of infection") + theme_cowplot()+scale_x_discrete(labels=c("-2","","0","","2","","","5","","","","","10","","","","","15","","",""))

ggsave(paste("~/Dropbox/varmodel2/ms/Figures/","FigS2_95DOI.pdf",sep=""),p6, width=8.8, height=4)


p7<-ggplot()+geom_rect(data = df, aes(xmin = yearWhole, xmax = upper, ymin = y, ymax = Inf), fill="pink", alpha=0.5)+geom_boxplot(data = allDur%>%filter( year>=79, year<=99), aes(as.factor(year), mdur, fill = selMode),outlier.shape = NA, lwd=0.2)+
  facet_grid(transmission~IRS,labeller = label_both)+ 
  scale_fill_manual(name = "Selection",values = c("blue","red"))+
  scale_color_manual(name = "during IRS",values = c("grey","black"))+
  xlab("Year")+ylab("Mean duration of infection")+ theme_cowplot()+scale_x_discrete(labels=c("-2","","0","","2","","","5","","","","","10","","","","","15","","",""))

ggsave(paste("~/Dropbox/varmodel2/ms/Figures/","FigS3_meanDOI.pdf", sep=""), p7, width=9, height=7)

p8<-ggplot()+ geom_rect(data = df, aes(xmin = yearWhole, xmax = upper, ymin = y, ymax = Inf), fill="pink", alpha=0.5)+
  geom_boxplot(data = summaryTable%>%filter( year>=79, year<=99,transmission=="high",ext==FALSE, genePool>=7200), aes(yearWhole,q50, fill = selMode), outlier.shape = NA, lwd=0.2)+
  facet_grid(genePool~IRS,labeller = label_both, scale = "free_y")+
  scale_fill_manual(name = "Selection",values = c("blue","red"))+
  xlab("Year")+ylab("5th percentile of PTS")+ theme_cowplot()+scale_x_discrete(labels=c("-2","","0","","2","","","5","","","","","10","","","","","15","","",""))
ggsave(paste("~/Dropbox/varmodel2/ms/Figures/","FigS4_PTS5perc.pdf", sep=""), p8, width=9, height=7)


###Fig4, new, get individual PTS distribution---------
ptsCat<-data.frame(variable = colnames(ptsTable[2:41]), PTS = seq(0.025,1,0.025))


ptsMelt<-melt(ptsTable%>%filter(num==27),
              id.vars = c("time","num","selMode","IRS","r"))
ptsMelt<-ptsMelt%>%left_join(ptsCat, by = "variable")%>%mutate(year = time/360, timeMark1=(year>80), timeMark2 = (year>80+IRS), group = timeMark1+timeMark2)
#create annotation
#2yr
ann_text <- data.frame(label = c("", "Gene pool size = 14400\nIRS = 2 yrs"),selMode   = c("G","S"),group = 0)

pts2y<-ggplot(ptsMelt%>%filter(IRS==2), aes(as.factor(PTS),value, fill=as.factor(group)))+geom_boxplot(outlier.shape = NA, lwd=0.2)+facet_grid(selMode~.)+
  scale_fill_manual(name = "IRS",values = c("orange","purple","darkgreen"), labels = c("before", "during","after")) +ylab("Density")+xlab("PTS")+
  geom_text(data = ann_text, aes(x =c(30,30),y = c(5,5),label=label))+
  theme_cowplot()+theme(legend.justification=c(0.9,0.95), legend.position=c(0.9,0.95))+scale_x_discrete(labels=c("","","","0.1","","","","0.2","","","","0.3","","","","0.4","","","","0.5","","","","0.6","","","","0.7","","","","0.8","","","","0.9","","","","1"))

#5yr
ann_text <- data.frame(label = c("", "Gene pool size = 14400\nIRS = 5 yrs"),selMode   = c("G","S"),group = 0)

pts5y<-ggplot(ptsMelt%>%filter(IRS==5), aes(as.factor(PTS),value, fill=as.factor(group)))+geom_boxplot(outlier.shape = NA, lwd=0.2)+facet_grid(selMode~.)+
  scale_fill_manual(name = "IRS",values = c("orange","purple","darkgreen"), labels = c("before", "during","after")) +ylab("Density")+xlab("PTS")+
  geom_text(data = ann_text, aes(x =c(30,30),y = c(5,5),label=label))+
  theme_cowplot()+theme(legend.justification=c(0.9,0.95), legend.position=c(0.9,0.95))+scale_x_discrete(labels=c("","","","0.1","","","","0.2","","","","0.3","","","","0.4","","","","0.5","","","","0.6","","","","0.7","","","","0.8","","","","0.9","","","","1"))

#10yr
ann_text <- data.frame(label = c("", "Gene pool size = 14400\nIRS = 10 yrs"),selMode   = c("G","S"),group = 0)

pts10y<-ggplot(ptsMelt%>%filter(IRS==10), aes(as.factor(PTS),value, fill=as.factor(group)))+geom_boxplot(outlier.shape = NA, lwd=0.2)+facet_grid(selMode~.)+
  scale_fill_manual(name = "IRS",values = c("orange","purple","darkgreen"), labels = c("before", "during","after")) +ylab("Density")+xlab("PTS")+
  geom_text(data = ann_text, aes(x =c(30,30),y = c(5,5),label=label))+
  theme_cowplot()+theme(legend.justification=c(0.9,0.95), legend.position=c(0.9,0.95))+scale_x_discrete(labels=c("","","","0.1","","","","0.2","","","","0.3","","","","0.4","","","","0.5","","","","0.6","","","","0.7","","","","0.8","","","","0.9","","","","1"))

#include in the main text 5yr intervention
ggsave("~/Dropbox/varmodel2/ms/Figures/Fig4_ptsDist.pdf",pts5y, width=6, height=4)
ggsave("~/Dropbox/varmodel2/ms/Figures/Fig4_ptsDist.png",pts5y, width=6, height=4)

#5 and 10yr include in suppl
plot2by2<-plot_grid(pts2y,pts10y, labels = c("(A)","(B)"), ncol =1)
save_plot(paste("~/Dropbox/varmodel2/ms/Figures/","SupplFig3_PTS.pdf", sep=""), plot2by2,
          ncol = 1, # we're saving a grid plot of 2 columns
          nrow = 2, # and 2 rows
          # each individual subplot should have an aspect ratio of 1.3
          base_aspect_ratio = 1.8,
          base_height = 4.5
)

save_plot(paste("~/Dropbox/varmodel2/ms/Figures/","SupplFig3_PTS.png", sep=""), plot2by2,
          ncol = 1, # we're saving a grid plot of 2 columns
          nrow = 2, # and 2 rows
          # each individual subplot should have an aspect ratio of 1.3
          base_aspect_ratio = 1.8,
          base_height = 4.5
)





###Fig 5 gene diversity-----
###network features
netTable<-read.table(netSameSize,sep="\t",header=F)
allTableNames<-c("time", "nicheDiv","alleleShannon","alleleSimpson","geneShannon","geneSimpson","f_01_averageLocalClusteringCoeff","f_02_averageLocalClusteringCoeffWeighted","f_03_globalClusteringCoeff","f_04_gdensity","f_05_proportionSingletons","f_06_proportionEndpoints","f_07_meanDegree","f_08_assortativityDegree","f_09_meanStrength","f_10_straightness","f_11_entropyDegreeDistribution","f_12_ratioComponents","f_13_giantConnectedRatio","f_14_evennessComponentSize","f_15_CentralPointDominance","f_16_meanEccentricity","f_17_gdiameter","f_18_meanDiameterComponents","f_19_globalEfficiency","f_20_averageClosenessCentrality","f_21_mot1","f_22_mot2","f_23_mot3","f_24_mot4","f_25_mot5","f_26_mot6","f_27_mot7","f_28_mot8","f_29_mot9","f_30_mot10","f_31_mot11","f_32_mot12","f_33_reciprocity","f_34_inoutCorrelation","f_35_communitySizeEvenness","f_36_numberCommonCommunities","f_37_CommunityRatio","f_38_modularity","f_39_meanFST","f_40_maxFST","f_41_minFST","f_42_giniFST","num","selMode","IRS","r",)
colnames(netTable)<-allTableNames

netTable<-as_tibble(netTable)

netTable<-left_join(netTable, allPara,by=c("num"="NO"))
netTable<-left_join(netTable, extCount,by=c("num","IRS","selMode","r"))
netTable<-netTable%>%mutate(year = time/360)

vlineDat<-data.frame(IRS = c(2,2,5,5,10,10), xint = c(80,82,80,85,80,90))

##calculate mean diversity and variation---
genDiv<-netTable%>%filter(IRS>0, BITING_RATE_MEAN==5E-04,  N_GENES_INITIAL==14400, time>25000, year<95, ext==F)%>%
  group_by(IRS, selMode,year)%>%summarize(mGS = mean(geneShannon),
                                     vGS = sd(geneShannon))

p8<-genDiv%>%
  ggplot(aes(year, mGS, col = selMode))+geom_line()+
  geom_errorbar(aes(ymin=mGS-vGS, ymax=mGS+vGS), width=.2,position=position_dodge(0.05),alpha=0.5)+
  theme_cowplot()+
  facet_grid(.~IRS,labeller = label_both )+ 
  scale_colour_manual(name = "Selection",values = c("blue","red"))+
  geom_vline(data=vlineDat, aes(xintercept=xint), linetype="dashed")+
  ylab("shannon diversity of genes")


summaryTable<-summaryTable%>%mutate(year = time/360)
summaryTable<-left_join(summaryTable, extCount,by=c("num","IRS","selMode","r"))

prevAv<-summaryTable%>%filter(IRS>0, 0.00010 == 0.0001, N_GENES_INITIAL==14400, year>78, year<95, ext==F)%>%
  group_by(IRS, selMode, year)%>%summarize(mP = mean(Prevalence),
                                           vP = sd(Prevalence))


p9<-prevAv%>%
  ggplot(aes(year, mP, col = selMode))+geom_line()+theme_cowplot()+  geom_errorbar(aes(ymin=mP-vP, ymax=mP+vP), width=.2,position=position_dodge(0.05),alpha=0.5)+
  facet_grid(.~IRS, scale = "free_y",labeller = label_both )+
  geom_vline(data=vlineDat, aes(xintercept=xint), linetype="dashed")+ 
  scale_colour_manual(name = "Selection",values = c("blue","red"))


plot2by2<-plot_grid(p8,p9, labels = c("(A)","(B)"), ncol =1)
save_plot(paste("~/Dropbox/varmodel2/ms/Figures/","Fig5_geneDiv.pdf", sep=""), plot2by2,
          ncol = 1, # we're saving a grid plot of 2 columns
          nrow = 2, # and 2 rows
          # each individual subplot should have an aspect ratio of 1.3
          base_aspect_ratio = 3,
          base_height = 3.5
)

save_plot(paste("~/Dropbox/varmodel2/ms/Figures/","Fig5_geneDiv.png", sep=""), plot2by2,
          ncol = 1, # we're saving a grid plot of 2 columns
          nrow = 2, # and 2 rows
          # each individual subplot should have an aspect ratio of 1.3
          base_aspect_ratio = 3,
          base_height = 3.5
)

###Fig 6 network trajectory--------
library(adegenet)
netNoCutTable<-netTable%>%mutate(year = ceiling(time/360))

###one example 14400
#for (gn in c(4800, 7200, 9600, 14400, 19200)) {
#  for (i in c(2,5,10)){
    for (gn in c(14400)) {
      for (i in c(5)){
    yearStart <- c(77, 80, 80+i)
    dur<-c(3,i,20-i)
        subTable<-netNoCutTable%>%filter(N_GENES_INITIAL==gn,
                                BITING_RATE_MEAN==5e-04, IRS==i, 
                                ext==FALSE)%>%mutate(timeMark1=(year>80), timeMark2 = (year>80+IRS), group = timeMark1+timeMark2)

        subTable<-netNoCutTable%>%filter(N_GENES_INITIAL==gn,
                                         BITING_RATE_MEAN==5e-04, IRS==i)%>%mutate(timeMark1=(year>80), timeMark2 = (year>80+IRS), group = timeMark1+timeMark2)
        
   subTable$yearStart <-yearStart[subTable$group+1]
   subTable$dur <- dur[subTable$group+1]
   
    subTable$group[subTable$group==1]<-"2_during IRS"
    subTable$group[subTable$group==2]<-"3_after IRS"
    subTable$group[subTable$group==0]<-"1_before IRS"
    testpca<-dudi.pca(df = subTable[, c(24:64)]%>%drop_na(), scale = T, scannf = FALSE, nf = 4)
    temp<-cbind(subTable%>%drop_na()%>%select(year, group, yearStart, dur, selMode), testpca$li)

    temp2<-temp%>%filter(group=="1_before IRS")
    temp2<-rbind(temp2,temp2,temp2)
    temp2$group<-rep(c("1_before IRS","2_during IRS","3_after IRS"), each = nrow(temp2)/3)
    
    
    p10<-ggplot(temp, aes(Axis1, Axis2, color = (year-yearStart)/dur))+
      geom_vline(xintercept = 0)+geom_hline(yintercept = 0)+
      geom_point(alpha=0.5, size=0.8)+stat_ellipse(data = temp2, linetype = 2)+
      facet_grid(selMode~group)+ 
      theme_light()+
      scale_color_gradientn(colors = c("orange","red","blue"), name = "time")+xlab("PC1")+ylab("PC2")+labs(title = paste("pool size =", gn, "IRS=", i, "years"))

        ggsave(paste("~/Dropbox/varmodel2/ms/Figures/network2/","gene_",gn,"_IRS",i, "_networkNoCut.pdf",sep = ""),width=6,height=4)
    
    p11<-fviz_contrib(testpca, choice = "var", axes = 1, top = 10)
    ggsave(paste("~/Dropbox/varmodel2/ms/Figures/network2/","gene_",gn,"_IRS",i, "_contrAx1.pdf",sep = ""),p,width=8,height=4)
    
    p12<-fviz_contrib(testpca, choice = "var", axes = 2, top = 10)
    ggsave(paste("~/Dropbox/varmodel2/ms/Figures/network2/","gene_",gn,"_IRS",i, "_contrAx2.pdf",sep = ""),p,width=8,height=4)
    
  }
}

#get the calculation of contributions of PCs
a1<-p11$data
a2<-p12$data
a1<-a1%>%arrange(desc(contrib))%>%mutate(PC = 1)
a2<-a2%>%arrange(desc(contrib))%>%mutate(PC = 2)
newDat<-rbind(a1[1:10,],a2[1:10,])

networkCode<-read.csv(paste("~/Dropbox/varmodel2/ms/Figures/","networkCode.csv", sep=""),header=T)
newDat<-newDat%>%left_join(networkCode, by = "name")

p11a<-ggplot(newDat%>%filter(PC==1), aes(x=reorder(Network.Properties,contrib), y=contrib,fill=Category)) +
  geom_bar(stat='identity') +
  coord_flip()+theme_cowplot()+xlab("")+
  geom_hline(yintercept = 100/41, linetype="dashed", color = "red")
p12a<-ggplot(newDat%>%filter(PC==2), aes(x=reorder(Network.Properties,contrib), y=contrib,fill=Category)) +
  geom_bar(stat='identity') +
  coord_flip()+theme_cowplot()+xlab("")+
  geom_hline(yintercept = 100/41, linetype="dashed", color = "red")

plot2by2<-plot_grid(p11a,p12a, labels = c("(A) PC1","(B) PC2"), ncol =1)
save_plot(paste("~/Dropbox/varmodel2/ms/Figures/","SupplFig_PCContrib.pdf", sep=""), plot2by2,
          ncol = 1, # we're saving a grid plot of 2 columns
          nrow = 2, # and 2 rows
          # each individual subplot should have an aspect ratio of 1.3
          base_aspect_ratio = 3.5
)

save_plot(paste("~/Dropbox/varmodel2/ms/Figures/","SupplFig_PCContrib.png", sep=""), plot2by2,
          ncol = 1, # we're saving a grid plot of 2 columns
          nrow = 2, # and 2 rows
          # each individual subplot should have an aspect ratio of 1.3
          base_aspect_ratio = 3.5
)



