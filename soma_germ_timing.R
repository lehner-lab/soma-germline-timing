library(ggplot2)
library(gplots)
library(binom)
library(MASS)
library(drc)


# set workind directory
wd <- '~/workdir/'
setwd(wd)
#read the data
data<-read.table('data.txt',sep='\t',header=T,stringsAsFactors = T)

data_filt<-na.omit(data)

#define the transitions (molting, first embryo)

mt<-(data_filt[,'transition'])
#mt<-factor(mt, levels=unique(mt))

#define the experimental factors

#I now use e description field as a quick way to include all relevant experimental info this makes things easier when you run the for loop to create the t50 table
description<-data_filt[,'Series']


#this for loop produces the table of T50 across all experimental factors (using the description field you only need to include one for loop for all experimental factors)
#it first fit a binomial glm then finds the t at which some fraction of the population have transitioned through a developmental milestone

# choose the estimated probability (i.e. 0.5 = t50, 0.9 = t90)
p<-0.5

# choose the confidence interval (0.95 = 95% confidence)
CI<-0.95

t50_table_tot<-data.frame(matrix(nrow=0,ncol=4))
for(d in unique(description)){
  #for(a in unique(starve)){
  #for(h in unique(hatch_day)){
  #  for(a in unique(mother_age)){
  for(m in levels(mt)){
    data_try<-data_filt[mt==m&description==d,]
    if(dim(data_try)[1]>1){
      
      #this function fits a glm the y is binomial format (before and after the transition) and x time
      fit<-glm(as.matrix(data_try[,2:1])~data_try[,'time.after.plate'],family=binomial('logit'))
      
      #this function finds the time at which p = proportion of the population have undergone the transition
      xp<-dose.p(fit,p=p)
      
      #this finds the 95% confidence intervals
      xp.ci <- xp + attr(xp, "SE") %*% matrix(qnorm(1-(1-CI)/2)*c(-1,1), nrow=1)
      
      t50 <- (cbind(xp, attr(xp, "SE"), xp.ci[,1], xp.ci[,2]))
      dimnames(t50)[[2]] <- c("time", "SE", "LCL","UCL")
      dimnames(t50)[[1]]<-paste(m,d)
      t50_table_tot<-rbind(t50_table_tot,t50)
    }    
  }
}
# }
# }


#number of transitions
ncol=2
t50_table_tot[,'transition']<-matrix(unlist(strsplit(rownames(t50_table_tot),' ')),ncol=ncol,byrow=T)[,1]


#add experimental factors to t50 table

#here I use the description field to summarize all experimental factors 

t50_table_tot[,'description']<-matrix(unlist(strsplit(rownames(t50_table_tot),' ')),ncol=ncol,byrow=T)[,2]

t50_table_tot[,'description']<-factor(t50_table_tot[,'description'],levels = sort(unique(t50_table_tot[,'description'])))


#prepare table for plotting

# subsets T50 table according to transition
L4_A<-t50_table_tot[t50_table_tot[,'transition']=='L4/A',]
embryo<-t50_table_tot[t50_table_tot[,'transition']=='embryo',]

#create a table with time difference between transitions to plot the shifts
data_diff<-embryo[1]-L4_A[1]
data_diff[,'SE']<-sqrt(embryo[,'SE']^2+L4_A[,'SE']^2)
data_diff[,'UCL']<-data_diff[,'time']+data_diff[,'SE']*1.96
data_diff[,'LCL']<-data_diff[,'time']-data_diff[,'SE']*1.96

#add experimental factors to data_diff table
data_diff[,'description']<-L4_A[,'description']

#write data diff table 
#write.table(data_diff,file='~/Dropbox (CRG)/work/expr_variance/ultimate_experiment/soma_germ_uncoupling/microscope/vit_RNAi_5_6_12_5_2017_data_diff.txt',sep='\t',col.names=NA)

#plot t50 data (absolute timings) using ggplot

#choose your x

x<-'description'

q<-ggplot(data = t50_table_tot, aes_string(x = x, y = 'time'
                           
                           , ymin = 'LCL'
                           
                           , ymax = 'UCL'
                           
                           ,col= 'transition'
)) +
  geom_point() +
  geom_errorbar( width = 0.5) +
  #facet_grid(~genotype)+
  #geom_line()+
  #coord_flip() +
  scale_colour_manual(values = c('dark red','blue','dark green')) +
  # theme_bw() +
  #ylim(c(3.3,5.5))+
  theme_classic() # + ylim(0,52)
q+theme(axis.text.x = element_text(angle = 90, hjust = 1))

#plot data diff table (shifts)
q<-ggplot(data = data_diff, aes_string(x = x, y = 'time'
                                
                                , ymin = 'LCL'
                                
                                , ymax = 'UCL'
                                
                                # ,col= juice
)) +
  geom_point() +
  geom_errorbar( width = 0.5) +
  #facet_grid(~genotype)+
  #geom_line()+
  #coord_flip() +
  #scale_colour_manual(values = c('dark red','blue','dark green')) +
  # theme_bw() +
  #ylim(c(2,6))+
  theme_classic()#+ ylim(43,54.5)
q+theme(axis.text.x = element_text(angle = 90, hjust = 1))


# Plot binomial glm fit


for(d in unique(description)){
  #for(a in unique(starve)){
  #for(h in unique(hatch_day)){
  #  for(a in unique(mother_age)){
  plot(c(39,53),c(0,1),type='n',xlab='time',ylab='proportion',main='d')
  for(m in levels(mt)){
    data_try<-data_filt[mt==m&description==d,]
    if(dim(data_try)[1]>1){
BA<-as.matrix(data_try[,2:1])
time<-data_try[,'time.after.plate']
fit<-glm(BA~time,family=binomial('logit'))
est<-predict(fit,type='resp',se.fit = T)
#points(time, est$fit)

if(m=='embryo'){
  sel<-t50_table_tot[,'transition']==m&t50_table_tot[,'description']==d
  t50<-t50_table_tot[sel,]
  time_inter<-seq(t50[1,1]-1.6,t50[1,1]+1.6,0.1)
  plotCI(time,est$fit,uiw= est$se.fit*1.96, add = T,err='y',col="red", barcol="red", lwd=1,gap=0.3)
  plotCI(t50[1,1],0.5,uiw= t50[1,4]-t50[1,1],liw= t50[1,1]-t50[1,3], add = T,err='x',col="red", barcol="red", lwd=1,gap=0.3)
  pred<-predict(fit,data.frame(time=time_inter),type='resp',se.fit = T)
  lines(time_inter,pred$fit,col='red',lwd=2)
  lines(time_inter,pred$fit-pred$se.fit,lty=2,col='red')
  lines(time_inter,pred$fit+pred$se.fit,lty=2,col='red')
} else {
  sel<-t50_table_tot[,'transition']==m&t50_table_tot[,'description']==d
  t50<-t50_table_tot[sel,]
  plotCI(time,est$fit,uiw= est$se.fit*1.96, add = T,err='y',col="dark green", barcol="dark green", lwd=1,gap=0.3)
  plotCI(t50[1,1],0.5,uiw= t50[1,4]-t50[1,1],liw=t50[1,1]-t50[1,3], add = T,err='x',col="dark green", barcol="dark green", lwd=1,gap=0.3)
  time_inter<-seq(t50[1,1]-1.6,t50[1,1]+1.6,0.1)
  pred<-predict(fit,data.frame(time=time_inter),type='resp',se.fit = T)
  lines(time_inter,pred$fit,col='dark green',lwd=2)
  lines(time_inter,pred$fit-pred$se.fit,lty=2,col='dark green')
  lines(time_inter,pred$fit+pred$se.fit,lty=2,col='dark green')
}

    }
  }
}

  

