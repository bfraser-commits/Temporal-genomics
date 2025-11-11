library(ggplot2)
library(grid)
library(plyr)
library(ggpubr)
library(gridExtra)
library(dplyr)
library(data.table)
library(viridis)
library(cowplot)

### This script will process the ROH and create the plot seen in fig 1C
### Read in the data
ghp_hom<-read.table("Documents/FIBR_temporal/ROH_plink/GHP_500K_1het.hom", header=TRUE)
C13_hom<-read.table("Documents/FIBR_temporal/ROH_plink/CA13_500K_1het.hom", header=T)
C18_hom<-read.table("Documents/FIBR_temporal/ROH_plink/CA18_500K_1het.hom", header=T)
T13_hom<-read.table("Documents/FIBR_temporal/ROH_plink/TA13_500K_1het.hom", header=TRUE)
T18_hom<-read.table("Documents/FIBR_temporal/ROH_plink/TA18_500K_1het.hom", header=TRUE)
LL13_hom<-read.table("Documents/FIBR_temporal/ROH_plink/LL13_500K_1het.hom", header=TRUE)
LL18_hom<-read.table("Documents/FIBR_temporal/ROH_plink/LL18_500K_1het.hom", header=TRUE)
UL13_hom<-read.table("Documents/FIBR_temporal/ROH_plink/UL13_500K_1het.hom", header=TRUE)
UL18_hom<-read.table("Documents/FIBR_temporal/ROH_plink/UL18_500K_1het.hom", header=TRUE)

# remove the sex chromosome
ghp_hom<-ghp_hom[!(ghp_hom$CHR==12),]
C13_hom<-C13_hom[!(C13_hom$CHR==12),]
C18_hom<-C18_hom[!(C18_hom$CHR==12),]
T13_hom<-T13_hom[!(T13_hom$CHR==12),]
T18_hom<-T18_hom[!(T18_hom$CHR==12),]
LL13_hom<-LL13_hom[!(LL13_hom$CHR==12),]
LL18_hom<-LL18_hom[!(LL18_hom$CHR==12),]
UL13_hom<-UL13_hom[!(UL13_hom$CHR==12),]
UL18_hom<-UL18_hom[!(UL18_hom$CHR==12),]

### Add population group
ghp_hom$POP<-as.factor(rep('GHP',nrow(ghp_hom)))
C13_hom$POP<-as.factor(rep('CA13',nrow(C13_hom)))
C18_hom$POP<-as.factor(rep('CA18',nrow(C18_hom)))
T13_hom$POP<-as.factor(rep('TA13',nrow(T13_hom)))
T18_hom$POP<-as.factor(rep('TA18',nrow(T18_hom)))
LL13_hom$POP<-as.factor(rep('LL13',nrow(LL13_hom)))
LL18_hom$POP<-as.factor(rep('LL18',nrow(LL18_hom)))
UL13_hom$POP<-as.factor(rep('UL13',nrow(UL13_hom)))
UL18_hom$POP<-as.factor(rep('UL18',nrow(UL18_hom)))

## Add MB column
ghp_hom$MB<-ghp_hom$KB/1000
C13_hom$MB<-C13_hom$KB/1000
C18_hom$MB<-C18_hom$KB/1000
T13_hom$MB<-T13_hom$KB/1000
T18_hom$MB<-T18_hom$KB/1000
LL13_hom$MB<-LL13_hom$KB/1000
LL18_hom$MB<-LL18_hom$KB/1000
UL13_hom$MB<-UL13_hom$KB/1000
UL18_hom$MB<-UL18_hom$KB/1000

####prep data###

#GHP
### Bin the data in 4 bins, and summarise: 1) the total length of ROHs in a bin and 2) the number of ROH's in a bin 
ghp_hom$GROUP<-cut(ghp_hom$MB, breaks=c(0.5,0.75,1,1.5,2,2.5,Inf),labels=c("0.50-0.75","0.75-1.0","1.0-1.5","1.5-2.0","2.0-2.5",">2.5"))
# sum per group per individual
ghp_sumind<-aggregate(list(MB=ghp_hom$MB), by=list(group=ghp_hom$GROUP,ind=ghp_hom$FID), FUN=sum)
# average length per group
ghp_sum<-aggregate(list(sum=ghp_sumind$MB), by=list(Category=ghp_sumind$group), FUN=mean)
ghp_sum$POP<-rep('GHP',nrow(ghp_sum))
ghp_sum$prop<-ghp_sum$sum/sum(ghp_sum$sum)
ghp_sum$propBUCK<-ghp_sum$sum/718
# counts per group per individual
ghp_freq2<- ddply(ghp_hom, .(ghp_hom$GROUP, ghp_hom$FID), nrow)
# average count per group
ghp_freq<- aggregate(list(mean_count=ghp_freq2$V1), by=list(Category=ghp_freq2$`ghp_hom$GROUP`), FUN=mean)
ghp_freq$POP<-rep('GHP',nrow(ghp_freq))
ghp_freq$prop<-ghp_freq$mean_count/sum(ghp_freq$mean_count)


#C13
C13_hom$GROUP<-cut(C13_hom$MB, breaks=c(0.5,0.75,1,1.5,2,2.5,Inf),labels=c("0.50-0.75","0.75-1.0","1.0-1.5","1.5-2.0","2.0-2.5",">2.5"))
C13_sumind<-aggregate(list(Mb=C13_hom$MB), by=list(group=C13_hom$GROUP,ind=C13_hom$FID), FUN=sum)
C13_sum<-aggregate(list(sum=C13_sumind$Mb), by=list(Category=C13_sumind$group), FUN=mean)
C13_sum$POP<-rep('CA13',nrow(C13_sum))
C13_sum$prop<-C13_sum$sum/sum(C13_sum$sum)
C13_sum$propBUCK<-C13_sum$sum/718
C13_freq2<- ddply(C13_hom, .(C13_hom$GROUP, C13_hom$FID), nrow)
C13_freq<- aggregate(list(mean_count=C13_freq2$V1), by=list(Category=C13_freq2$`C13_hom$GROUP`), FUN=mean)
C13_freq$POP<-rep('CA13',nrow(C13_freq))
C13_freq$prop<-C13_freq$mean_count/sum(C13_freq$mean_count)

#CA18
C18_hom$GROUP<-cut(C18_hom$MB, breaks=c(0.5,0.75,1,1.5,2,2.5,Inf),labels=c("0.50-0.75","0.75-1.0","1.0-1.5","1.5-2.0","2.0-2.5",">2.5"))
C18_sumind<-aggregate(list(Mb=C18_hom$MB), by=list(group=C18_hom$GROUP,ind=C18_hom$FID), FUN=sum)
C18_sum<-aggregate(list(sum=C18_sumind$Mb), by=list(Category=C18_sumind$group), FUN=mean)
C18_sum$POP<-rep('CA18',nrow(C18_sum))
C18_sum$prop<-C18_sum$sum/sum(C18_sum$sum)
C18_sum$propBUCK<-C18_sum$sum/718
C18_freq2<- ddply(C18_hom, .(C18_hom$GROUP, C18_hom$FID), nrow)
C18_freq<- aggregate(list(mean_count=C18_freq2$V1), by=list(Category=C18_freq2$`C18_hom$GROUP`), FUN=mean)
C18_freq$POP<-rep('CA18',nrow(C18_freq))
C18_freq$prop<-C18_freq$mean_count/sum(C18_freq$mean_count)


#TA13
T13_hom$GROUP<-cut(T13_hom$MB, breaks=c(0.5,0.75,1,1.5,2,2.5,Inf),labels=c("0.50-0.75","0.75-1.0","1.0-1.5","1.5-2.0","2.0-2.5",">2.5"))
T13_sumind<-aggregate(list(MB=T13_hom$MB), by=list(group=T13_hom$GROUP,ind=T13_hom$FID), FUN=sum)
T13_sum<-aggregate(list(sum=T13_sumind$MB), by=list(Category=T13_sumind$group), FUN=mean)
T13_sum$POP<-rep('TA13',nrow(T13_sum))
T13_sum$prop<-T13_sum$sum/sum(T13_sum$sum)
T13_sum$propBUCK<-T13_sum$sum/718
T13_freq2<- ddply(T13_hom, .(T13_hom$GROUP, T13_hom$FID), nrow)
T13_freq<- aggregate(list(mean_count=T13_freq2$V1), by=list(Category=T13_freq2$`T13_hom$GROUP`), FUN=mean)
T13_freq$POP<-rep('TA13',nrow(T13_freq))
T13_freq$prop<-T13_freq$mean_count/sum(T13_freq$mean_count)

#TA18
T18_hom$GROUP<-cut(T18_hom$MB, breaks=c(0.5,0.75,1,1.5,2,2.5,Inf),labels=c("0.50-0.75","0.75-1.0","1.0-1.5","1.5-2.0","2.0-2.5",">2.5"))
T18_sumind<-aggregate(list(MB=T18_hom$MB), by=list(group=T18_hom$GROUP,ind=T18_hom$FID), FUN=sum)
T18_sum<-aggregate(list(sum=T18_sumind$MB), by=list(Category=T18_sumind$group), FUN=mean)
T18_sum$POP<-rep('TA18',nrow(T18_sum))
T18_sum$prop<-T18_sum$sum/sum(T18_sum$sum)
T18_sum$propBUCK<-T18_sum$sum/718
T18_freq2<- ddply(T18_hom, .(T18_hom$GROUP, T18_hom$FID), nrow)
T18_freq<- aggregate(list(mean_count=T18_freq2$V1), by=list(Category=T18_freq2$`T18_hom$GROUP`), FUN=mean)
T18_freq$POP<-rep('TA18',nrow(T18_freq))
T18_freq$prop<-T18_freq$mean_count/sum(T18_freq$mean_count)

#LL13
LL13_hom$GROUP<-cut(LL13_hom$MB, breaks=c(0.5,0.75,1,1.5,2,2.5,Inf),labels=c("0.50-0.75","0.75-1.0","1.0-1.5","1.5-2.0","2.0-2.5",">2.5"))
LL13_sumind<-aggregate(list(MB=LL13_hom$MB), by=list(group=LL13_hom$GROUP,ind=LL13_hom$FID), FUN=sum)
LL13_sum<-aggregate(list(sum=LL13_sumind$MB), by=list(Category=LL13_sumind$group), FUN=mean)
LL13_sum$POP<-rep('LL13',nrow(LL13_sum))
LL13_sum$prop<-LL13_sum$sum/sum(LL13_sum$sum)
LL13_sum$propBUCK<-LL13_sum$sum/718
LL13_freq2<- ddply(LL13_hom, .(LL13_hom$GROUP, LL13_hom$FID), nrow)
LL13_freq<- aggregate(list(mean_count=LL13_freq2$V1), by=list(Category=LL13_freq2$`LL13_hom$GROUP`), FUN=mean)
LL13_freq$POP<-rep('LL13',nrow(LL13_freq))
LL13_freq$prop<-LL13_freq$mean_count/sum(LL13_freq$mean_count)

#LL18
LL18_hom$GROUP<-cut(LL18_hom$MB, breaks=c(0.5,0.75,1,1.5,2,2.5,Inf),labels=c("0.50-0.75","0.75-1.0","1.0-1.5","1.5-2.0","2.0-2.5",">2.5"))
LL18_sumind<-aggregate(list(MB=LL18_hom$MB), by=list(group=LL18_hom$GROUP,ind=LL18_hom$FID), FUN=sum)
LL18_sum<-aggregate(list(sum=LL18_sumind$MB), by=list(Category=LL18_sumind$group), FUN=mean)
LL18_sum$POP<-rep('LL18',nrow(LL18_sum))
LL18_sum$prop<-LL18_sum$sum/sum(LL18_sum$sum)
LL18_sum$propBUCK<-LL18_sum$sum/718
LL18_freq2<- ddply(LL18_hom, .(LL18_hom$GROUP, LL18_hom$FID), nrow)
LL18_freq<- aggregate(list(mean_count=LL18_freq2$V1), by=list(Category=LL18_freq2$`LL18_hom$GROUP`), FUN=mean)
LL18_freq$POP<-rep('LL18',nrow(LL18_freq))
LL18_freq$prop<-LL18_freq$mean_count/sum(LL18_freq$mean_count)

#UL13
UL13_hom$GROUP<-cut(UL13_hom$MB, breaks=c(0.5,0.75,1,1.5,2,2.5,Inf),labels=c("0.50-0.75","0.75-1.0","1.0-1.5","1.5-2.0","2.0-2.5",">2.5"))
UL13_sumind<-aggregate(list(MB=UL13_hom$MB), by=list(group=UL13_hom$GROUP,ind=UL13_hom$FID), FUN=sum)
UL13_sum<-aggregate(list(sum=UL13_sumind$MB), by=list(Category=UL13_sumind$group), FUN=mean)
UL13_sum$POP<-rep('UL13',nrow(UL13_sum))
UL13_sum$prop<-UL13_sum$sum/sum(UL13_sum$sum)
UL13_sum$propBUCK<-UL13_sum$sum/718
UL13_freq2<- ddply(UL13_hom, .(UL13_hom$GROUP, UL13_hom$FID), nrow)
UL13_freq<- aggregate(list(mean_count=UL13_freq2$V1), by=list(Category=UL13_freq2$`UL13_hom$GROUP`), FUN=mean)
UL13_freq$POP<-rep('UL13',nrow(UL13_freq))
UL13_freq$prop<-UL13_freq$mean_count/sum(UL13_freq$mean_count)

#UL18
UL18_hom$GROUP<-cut(UL18_hom$MB, breaks=c(0.5,0.75,1,1.5,2,2.5,Inf),labels=c("0.50-0.75","0.75-1.0","1.0-1.5","1.5-2.0","2.0-2.5",">2.5"))
UL18_sumind<-aggregate(list(MB=UL18_hom$MB), by=list(group=UL18_hom$GROUP,ind=UL18_hom$FID), FUN=sum)
UL18_sum<-aggregate(list(sum=UL18_sumind$MB), by=list(Category=UL18_sumind$group), FUN=mean)
UL18_sum$POP<-rep('UL18',nrow(UL18_sum))
UL18_sum$prop<-UL18_sum$sum/sum(UL18_sum$sum)
UL18_sum$propBUCK<-UL18_sum$sum/718
UL18_freq2<- ddply(UL18_hom, .(UL18_hom$GROUP, UL18_hom$FID), nrow)
UL18_freq<- aggregate(list(mean_count=UL18_freq2$V1), by=list(Category=UL18_freq2$`UL18_hom$GROUP`), FUN=mean)
UL18_freq$POP<-rep('UL18',nrow(UL18_freq))
UL18_freq$prop<-UL18_freq$mean_count/sum(UL18_freq$mean_count)


freq_all<-rbind(ghp_freq,T13_freq,T18_freq, C13_freq, C18_freq, UL13_freq, UL18_freq, LL13_freq, LL18_freq)
freq_all$POP<-factor(freq_all$POP, levels = c('GHP','CA13', 'CA18', 'TA13', 'TA18', 'LL13', 'LL18', 'UL13', 'UL18'))

sum_all<-rbind(ghp_sum,T13_sum, T18_sum, C13_sum, C18_sum, UL13_sum, UL18_sum, LL13_sum, LL18_sum)
sum_all$POP<-factor(sum_all$POP, levels = c('GHP','CA13', 'CA18', 'TA13', 'TA18', 'LL13', 'LL18', 'UL13', 'UL18'))

write.table(freq_all,"Documents/FIBR_temporal/ROH_plink/freq_of_ROH_auto_500K.txt", quote = F, sep = '\t', row.names = F)
write.table(sum_all,"Documents/FIBR_temporal/ROH_plink/sum_of_ROH_auto_500K.txt", quote = F, sep = '\t', row.names = F)

### plot bars. On x is the size bins. Y = the sum of the ROH's for that bin in each population
sums<-ggplot(data=sum_all, aes(x=POP,y=sum,fill=Category))+
  geom_bar(stat="identity",width = 0.9,position=position_dodge2(width = 2, preserve = "single"),color='#238443')+
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text = element_text(size=20),
        axis.title=element_text(size=20),
        plot.title = element_text(hjust = 0.5,size=20),
        legend.text = element_text(size = 16),
        legend.title = element_text(size=16),
        plot.margin=unit(c(0.5,0.5,0.5,0.5),"cm"),
        legend.position = c(0.60, 0.85), 
        legend.direction = "horizontal",  
      legend.box = "horizontal")+ 
    scale_fill_manual(values = c('#f7fcb9','#addd8e','#41ab5d','#005a32'))+
  scale_y_continuous(name='Sum of ROHs (Mb)',expand = c(0,0))+
  xlab('Population')+
  labs(fill = "Size class (Mb)")+
  coord_cartesian(ylim=c(0,15))+
  guides(fill=guide_legend(ncol=4))


png('Documents/FIBR_temporal/ROH_plink/sums_figure.png',width = 10,height = 5, units='in', res=300)
sums
dev.off()

# plot the nr of ROH
freqs<-ggplot(data=freq_all, aes(x=POP,y=mean_count,fill=Category))+
  geom_bar(stat="identity",width = 0.9,position=position_dodge2(width = 2, preserve = "single"),color='#238443')+
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text = element_text(size=20),
        axis.title = element_text(size=20),
        axis.title.x=element_blank(),
        axis.title.y=element_text(size=20),
        plot.title = element_text(hjust = 0.5,size=20),
        legend.text = element_text(size = 16),
        legend.title = element_text(size=16),
        plot.margin=unit(c(0.5,0.5,0.5,0.5),"cm"),
        legend.position = c(0.60, 0.85), 
        legend.direction = "horizontal",  
        legend.box = "horizontal")+ 
  scale_fill_manual(values =  c('#f7fcb9','#addd8e','#41ab5d','#005a32'))+
  scale_y_continuous(name='Count of ROHs (Mb)',expand = c(0,0))+
  xlab('Population')+
  labs(fill = "Size class (Mb)")+
  coord_cartesian(ylim=c(0,20))

png('Documents/FIBR_temporal/ROH_plink/freqs_of_ROH_auto_500K_bars_line.png',width = 800,height = 300)
freqs
dev.off()


# mean Mb in ROH, min length + max length
ghp_length <- aggregate(ghp_hom$MB, by = list(ind=ghp_hom$FID), FUN= sum)
ghp_mean<-mean(ghp_length$x)
ghp_mean
min(ghp_length$x)
max(ghp_length$x)

C13_length <- aggregate(C13_hom$MB, by = list(ind=C13_hom$FID), FUN= sum)
C13_mean<-mean(C13_length$x)
C13_mean
min(C13_length$x)
max(C13_length$x)

C18_length <- aggregate(C18_hom$MB, by = list(ind=C18_hom$FID), FUN= sum)
C18_mean<-mean(C18_length$x)
C18_mean
min(C18_length$x)
max(C18_length$x)

T13_length <- aggregate(T13_hom$MB, by = list(ind=T13_hom$FID), FUN= sum)
T13_mean<-mean(T13_length$x)
T13_mean
min(T13_length$x)
max(T13_length$x)

T18_length <- aggregate(T18_hom$MB, by = list(ind=T18_hom$FID), FUN= sum)
T18_mean<-mean(T18_length$x)
T18_mean
min(T18_length$x)
max(T18_length$x)

LL13_length <- aggregate(LL13_hom$MB, by = list(ind=LL13_hom$FID), FUN= sum)
LL13_mean<-mean(LL13_length$x)
LL13_mean
min(LL13_length$x)
max(LL13_length$x)

LL18_length <- aggregate(LL18_hom$MB, by = list(ind=LL18_hom$FID), FUN= sum)
LL18_mean<-mean(LL18_length$x)
LL18_mean
min(LL18_length$x)
max(LL18_length$x)

UL13_length <- aggregate(UL13_hom$MB, by = list(ind=UL13_hom$FID), FUN= sum)
UL13_mean<-mean(UL13_length$x)
UL13_mean
min(UL13_length$x)
max(UL13_length$x)

UL18_length <- aggregate(UL18_hom$MB, by = list(ind=UL18_hom$FID), FUN= sum)
UL18_mean<-mean(UL18_length$x)
UL18_mean
min(UL18_length$x)
max(UL18_length$x)


# calculate the inbreeding coefficient
## Froh = SUM(LRroh)/Lauto
ghp_froh<-aggregate(list(Mb=ghp_hom$MB), by=list(ind=ghp_hom$FID), FUN=sum)
ghp_froh$FROH<-ghp_froh$Mb/718
mean(ghp_froh$FROH)
min(ghp_froh$FROH)
max(ghp_froh$FROH)
ghp_froh$POP <- "GH"
ghp_froh

C13_froh<-aggregate(list(Mb=C13_hom$MB), by=list(ind=C13_hom$FID), FUN=sum)
C13_froh$FROH<-C13_froh$Mb/718
mean(C13_froh$FROH)
min(C13_froh$FROH)
max(C13_froh$FROH)
C13_froh$POP <- "C13"
C13_froh

C18_froh<-aggregate(list(Mb=C18_hom$MB), by=list(ind=C18_hom$FID), FUN=sum)
C18_froh$FROH<-C18_froh$Mb/718
mean(C18_froh$FROH)
min(C18_froh$FROH)
max(C18_froh$FROH)
C18_froh$POP <- "C18"

T13_froh<-aggregate(list(Mb=T13_hom$MB), by=list(ind=T13_hom$FID), FUN=sum)
T13_froh$FROH<-T13_froh$Mb/718
mean(T13_froh$FROH)
min(T13_froh$FROH)
max(T13_froh$FROH)
T13_froh$POP <- "T13"

T18_froh<-aggregate(list(Mb=T18_hom$MB), by=list(ind=T18_hom$FID), FUN=sum)
T18_froh$FROH<-T18_froh$Mb/718
mean(T18_froh$FROH)
min(T18_froh$FROH)
max(T18_froh$FROH)
T18_froh$POP <- "T18"

LL13_froh<-aggregate(list(Mb=LL13_hom$MB), by=list(ind=LL13_hom$FID), FUN=sum)
LL13_froh$FROH<-LL13_froh$Mb/718
mean(LL13_froh$FROH)
min(LL13_froh$FROH)
max(LL13_froh$FROH)
LL13_froh$POP <- "LL13"

LL18_froh<-aggregate(list(Mb=LL18_hom$MB), by=list(ind=LL18_hom$FID), FUN=sum)
LL18_froh$FROH<-LL18_froh$Mb/718
mean(LL18_froh$FROH)
min(LL18_froh$FROH)
max(LL18_froh$FROH)
LL18_froh$POP <- "LL18"

UL13_froh<-aggregate(list(Mb=UL13_hom$MB), by=list(ind=UL13_hom$FID), FUN=sum)
UL13_froh$FROH<-UL13_froh$Mb/718
mean(UL13_froh$FROH)
min(UL13_froh$FROH)
max(UL13_froh$FROH)
UL13_froh$POP <- "UL13"

UL18_froh<-aggregate(list(Mb=UL18_hom$MB), by=list(ind=UL18_hom$FID), FUN=sum)
UL18_froh$FROH<-UL18_froh$Mb/718
mean(UL18_froh$FROH)
min(UL18_froh$FROH)
max(UL18_froh$FROH)
UL18_froh$POP <- "UL18"


#do a MWU on GH vs each pop
GH_C13 <- data.frame(rbind(ghp_froh, C13_froh))
wilcox.test(FROH ~ POP, data = GH_C13)

GH_C18 <- data.frame(rbind(ghp_froh, C18_froh))
wilcox.test(FROH ~ POP, data = GH_C18)

GH_T13 <- data.frame(rbind(ghp_froh, T13_froh))
wilcox.test(FROH ~ POP, data = GH_T13)

GH_T18 <- data.frame(rbind(ghp_froh, T18_froh))
wilcox.test(FROH ~ POP, data = GH_T18)

GH_LL13 <- data.frame(rbind(ghp_froh, LL13_froh))
wilcox.test(FROH ~ POP, data = GH_LL13)

GH_LL18 <- data.frame(rbind(ghp_froh, LL18_froh))
wilcox.test(FROH ~ POP, data = GH_LL18)

GH_UL13 <- data.frame(rbind(ghp_froh, UL13_froh))
wilcox.test(FROH ~ POP, data = GH_UL13)

GH_UL18 <- data.frame(rbind(ghp_froh, UL18_froh))
wilcox.test(FROH ~ POP, data = GH_UL18)


#do a MWU between 13 and 18 of each pop
C13_C18 <- data.frame(rbind(C18_froh, C13_froh))
wilcox.test(FROH ~ POP, data = C13_C18)

T13_T18 <- data.frame(rbind(T18_froh, T13_froh))
wilcox.test(FROH ~ POP, data = T13_T18)


LL13_LL18 <- data.frame(rbind(LL18_froh, LL13_froh))
wilcox.test(FROH ~ POP, data = LL13_LL18)


UL13_UL18 <- data.frame(rbind(UL18_froh, UL13_froh))
wilcox.test(FROH ~ POP, data = UL13_UL18)

