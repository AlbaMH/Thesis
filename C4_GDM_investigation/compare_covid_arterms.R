# Compare AR versions of COVID model 




################################
## Nowcasting and Forecasting  #
## COVID-19 and Other Diseases #
################################

# Code structure..

############
# Preamble #
############

# Set your working directory.

setwd("~/OneDrive/PhD/Covid delay Biometrics revisions/AR COVID")

# The following packages will need to be installed and loaded.

library(tidyverse)
library(reshape2)
library(mgcv)
library(coda)
library(ggfan)
library(gridExtra)
library(viridis)
library(scales)
library(abind)
library(nimble)
library(doParallel)
library(lubridate)
library(openxlsx)

# Read in some necessary functions.
source('Functions.R')

data_file <- 'COVID_Data.xlsx'

# Load in the data for each region.
library(readxl)
EOE <- openxlsx::read.xlsx(data_file,sheet=2,rowNames=TRUE,detectDates=TRUE)
LDN <- openxlsx::read.xlsx(data_file,sheet=3,rowNames=TRUE,detectDates=TRUE)
MID <- openxlsx::read.xlsx(data_file,sheet=4,rowNames=TRUE,detectDates=TRUE)
NEY <- openxlsx::read.xlsx(data_file,sheet=5,rowNames=TRUE,detectDates=TRUE)
NW <- openxlsx::read.xlsx(data_file,sheet=6,rowNames=TRUE,detectDates=TRUE)
SE <- openxlsx::read.xlsx(data_file,sheet=7,rowNames=TRUE,detectDates=TRUE)
SW <- openxlsx::read.xlsx(data_file,sheet=8,rowNames=TRUE,detectDates=TRUE)

# The raw data is arranged by date of death and date reported. We need
# it to be arranged by date of death and days of delay since death.

# Begin processing data into the format for modelling.
regions_raw <- abind(EOE,LDN,MID,NEY,NW,SE,SW,along=3)[-(1:32),,]

L <- dim(regions_raw)[1] # Number of rows in the raw data.
C <- 33 # Length of time series up to the present day in the data.
N=C+7 # Length of the time series including 7 days of forecasting.
S=7 # Number of regions (excluding England).

for(s in 1:S){
  for(i in 2:L){
    regions_raw[i,1:(i-1),s] <- NA # Partial reports which haven't happened yet are NA.
  }
}
regions <- array(NA,dim=c(dim(regions_raw)[1]+N-C,dim(regions_raw)[2],S))

for(s in 1:S){
  for(i in 1:L){
    regions[i,,s] <- as.numeric(c(regions_raw[i,!is.na(regions_raw[i,,s]),s],
                                  rep(NA,sum(is.na(regions_raw[i,,s])))))
  }
}
# Now the data is in a suitable format for modelling.

D_max <- 14 # Cut-off / maximum delay.
D=6 # Number of delays to explicitly model (see Stoner and Economou 2019 (Biometrics) for more details).
n_chains <- 4 # Number of MCMC chains to run. Reduce if you don't have more than 4 CPU cores!
n_knots <- c(8,8) # Number of knots of the temporal and seasonal splines, respectively.

dates <- seq(as.Date("2020/4/2"),by=1,length=L) # Dates corresponding to the data.
days <- weekdays(dates)
days <- days%>%
  factor(levels=c('Monday','Tuesday','Wednesday','Thursday','Friday','Saturday','Sunday'))%>%
  as.numeric()

# Put the raw data into a data frame for plotting.
raw_data=tibble(t=dates,day=factor(c('Monday','Tuesday','Wednesday','Thursday','Friday','Saturday','Sunday')[days],
                                   levels=c('Monday','Tuesday','Wednesday','Thursday','Friday','Saturday','Sunday')),
                act=apply(regions_raw,1,sum,na.rm=T),ann=apply(regions_raw,2,sum,na.rm=T))
raw_data=rbind(mutate(raw_data,type='New Deaths'),mutate(raw_data,ann=cumsum(ann),act=cumsum(act),type='Cumulative'))%>%
  mutate(type=factor(type,levels=c('New Deaths','Cumulative')))

## Figure 2: actual deaths versus announced deaths.
new_deaths <- ggplot(filter(raw_data,type=='New Deaths'))+
  geom_smooth(aes(x=t,y=act,linetype='Trend'),se=FALSE,colour='gray40',size=0.5)+
  geom_line(aes(x=t,y=act,linetype='Actual'))+geom_line(aes(x=t,y=ann,linetype='Announced'))+
  geom_point(aes(x=t,y=ann,colour=day,shape=day),size=2)+
  labs(x='Date',y='Deaths',title='Hospital Deaths from COVID-19 in England')+
  theme_minimal()+
  scale_color_brewer(name=NULL,palette='Dark2')+
  scale_shape_manual(name=NULL,values=1:7)+
  scale_linetype_manual(values=c(3,2,1),name=NULL)+scale_x_date(limits=as.Date(c('2020-04-07','2020-05-23')))+
  theme(text = element_text(size=12), plot.title = element_text(size=16),plot.subtitle=element_text(size=14),
        legend.title = element_text(size=14),strip.text = element_text(size=14))
#ggsave(new_deaths,filename='Plots/new_deaths.pdf',width=9,height=3.25)

# Vector containing region names.
region_names <- c('East of England','London','Midlands',
                  'North East and Yorkshire',
                  'North West','South East','South West')

# Store counts in arrays (y_full/z_full is uncensored and y/z is censored).
z_full <- aperm(regions,c(1,3,2))

z_england <- apply(z_full,c(1,3),sum)
z_data_england <- cbind(z_england[,1:4],apply(z_england[,5:dim(z_england)[2]],1,sum,na.rm=T))%>%melt(varnames=c('t','d'))%>%
  mutate(d=as.character(d))
z_data_england$d[z_data_england$d=='5']='5+'
z_data_england <- mutate(z_data_england,d=factor(z_data_england$d,levels=rev(c(as.character(1:7),'5+'))),t=dates[t])

# Total reported daily deaths after 14 days of dealy.
y_full <- apply(regions[,1:D_max,],c(1,3),sum)

z <- z_full[1:N,,]
for(s in 1:S){
  z[,s,][outer(1:dim(z[,s,])[1], 0:(dim(z[,s,])[2]-1), FUN = "+") > C] <- NA
}

y <- apply(z[,,1:D_max],c(1,2),sum)

y_so_far <- apply(z[,,1:D_max],c(1,2),sum,na.rm=T)

# % reported within D_max
t_max <- sum(!is.na(apply(y_full,1,sum)))
sum(z_full[1:t_max,,1:7],na.rm=T)/sum(z_full[1:t_max,,1:28],na.rm=T)
sum(z_full[1:t_max,,1:14],na.rm=T)/sum(z_full[1:t_max,,1:28],na.rm=T)

## Figure 1: Example delayed reporting bar plot.
n_time=5
n_delay=5
tick_labels=numeric(n_time)
for(i in 1:(n_time-1)){
  tick_labels[i]=paste('t',as.character(-n_time+i),sep='')
}
tick_labels[n_time]='t'
bar_data=data.frame(time=rep(1:n_time,n_delay),delay=as.factor(sort(rep(1:n_delay,n_time))),
                    count=as.numeric(z[(C-n_time+1):C,1,1:n_delay]))
bar_plot=ggplot()+geom_col(data=data.frame(time=1:n_time,total=y_full[(C-n_time+1):C,1]),aes(x=time,y=total))+
  geom_col(data=bar_data,aes(x=time,y=count,fill=delay),position = position_stack(reverse = TRUE))+
  labs(x='Day of Death',y='Reported Deaths',title='Delay Structure of COVID-19 Data',subtitle='Hospital Deaths in East of England') +
  theme_minimal()+scale_fill_brewer(palette='Dark2',name='Delay')+scale_x_reverse(breaks=1:n_time,labels=tick_labels)+coord_flip()
#ggsave(bar_plot,file='Plots/bar_plot_covid.pdf',width=4.5,height=3)

# Total cases data into long format.
y_data <- y_full[1:N,]%>%melt(varnames=c('t','s'),value.name='y')%>%mutate(r=region_names[s])
Y_data <- tibble(t=1:N,y=apply(y_full[1:N,],1,sum))%>%mutate(s=S+1,r='England')
y_data <- rbind(y_data,Y_data)%>%mutate(r=factor(r,levels=c(region_names,'England')))

# plot residuals 
load("~/media/alba/Disk 1/OneDrive/PhD/Thesis code/Chapter 5 COVID-19 England/AR COVID/covid_samples_noAR.RData")
load("~/media/alba/Disk 1/OneDrive/PhD/Thesis code/Chapter 5 COVID-19 England/AR COVID/covid_samples_both.RData")
load("~/media/alba/Disk 1/OneDrive/PhD/Thesis code/Chapter 5 COVID-19 England/AR COVID/covid_samples_justAR.RData")


# covid model with no AR term.
# Combine all MCMC chains.
covid_combined_samples_noAR <- as_tibble(do.call('rbind',covid_samples_noAR))

# Number of MCMC samples.
n_sim <- dim(covid_combined_samples_noAR)[1]
# Negative-Binomial dispersion parameters.
theta_noAR<-which(str_detect(colnames(covid_combined_samples_noAR),c('theta')))
covid_theta_noAR<-covid_combined_samples_noAR[,theta_noAR]%>%as.matrix()%>%array(dim=c(n_sim,S))
# Negative-Binomial means.
lambda_noAR<-which(str_detect(colnames(covid_combined_samples_noAR),c('lambda')))
covid_lambda_noAR<-covid_combined_samples_noAR[,lambda_noAR]%>%as.matrix()%>%array(dim=c(n_sim,N,S))
# Now-casting and forecasting samples of total covid deaths.
covid_y_noAR <- array(dim=c(n_sim,N,S))
y_noAR<-which(str_detect(colnames(covid_combined_samples_noAR),c('y')))
covid_y_noAR[,1:C,]<-covid_combined_samples_noAR[,y_noAR]%>%as.matrix()%>%array(dim=c(n_sim,C,S))
for(s in 1:S){
  covid_y_noAR[,(C+1):N,s] <- rnbinom(n_sim*(N-C),mu=covid_lambda_noAR[,(C+1):N,s],size=covid_theta_noAR[,s])
}
# Samples of the total for England.
covid_Y_noAR <- apply(covid_y_noAR,c(1,2),sum)
# Total cases for each region.
covid_y_noAR <- covid_y_noAR%>%apply(c(2,3),quantile,c(0.025,0.1,0.175,0.25,0.5,0.75,0.825,0.9,0.975))%>%melt(varnames=c('quantile','t','s'),value.name='y')%>%
  spread(quantile,y)%>%mutate(r=region_names[s])
covid_Y_noAR <- covid_Y_noAR%>%apply(2,quantile,c(0.025,0.1,0.175,0.25,0.5,0.75,0.825,0.9,0.975))%>%melt(varnames=c('quantile','t'),value.name='y')%>%
  spread(quantile,y)%>%mutate(s=S+1,r='England')
covid_y_all_noAR <- rbind(covid_y_noAR,covid_Y_noAR)%>%mutate(r=factor(r,levels=c(region_names,'England')))

# Check Residuals.
lambda_medians_noAR <- covid_lambda_noAR%>%apply(c(2,3),quantile,c(0.025,0.1,0.175,0.25,0.5,0.75,0.825,0.9,0.975))%>%melt(varnames=c('quantile','t','s'),value.name='y')%>%
  spread(quantile,y)%>%mutate(r=region_names[s])
#res_s<-1
scaled_residuals_noAR<-matrix(nrow=N,ncol=S)
for(res_s in 1:S){
  all_data<-inner_join(covid_y_all_noAR,y_data, by=c("r","t","s"))
  predicted_data<-filter(all_data, s==res_s)
  mu<-filter(lambda_medians_noAR, s==res_s)
  var_fn<-mu$`50%`+mu$`50%`^2/median(covid_theta_noAR[,res_s])
  residuals<-(predicted_data$y-mu$`50%`)
  scaled_residuals_noAR[,res_s]<-(predicted_data$y-mu$`50%`)/sqrt(var_fn)
  
  # Residuals over time.
  #plot(x=1:length(residuals), y=residuals)
  #plot(x=1:length(residuals), y=scaled_residuals)
  plot(x=1:length(residuals), y=scaled_residuals_noAR[,res_s], type="l")
  
  # Residual ACF
  # calculate the ACF for lags between 1 and 21 (inclusive).
#  autocorrelation <- acf(scaled_residuals)
 # partial_autocorrelation <- pacf(scaled_residuals)
  
  # k_check<- gam(predicted_data$y~s(predicted_data$t, bs='cs' ,k =27))
  # gam.check(k_check)
}


# covid model with just an AR term.
# Combine all MCMC chains.
covid_combined_samples_justAR <- as_tibble(do.call('rbind',covid_samples_ar))

# Number of MCMC samples.
n_sim <- dim(covid_combined_samples_justAR)[1]
# Negative-Binomial dispersion parameters.
theta_justAR<-which(str_detect(colnames(covid_combined_samples_justAR),c('theta')))
covid_theta_justAR<-covid_combined_samples_justAR[,theta_justAR]%>%as.matrix()%>%array(dim=c(n_sim,S))
# Negative-Binomial means.
lambda_justAR<-which(str_detect(colnames(covid_combined_samples_justAR),c('lambda')))
covid_lambda_justAR<-covid_combined_samples_justAR[,lambda_justAR]%>%as.matrix()%>%array(dim=c(n_sim,N,S))
# Now-casting and forecasting samples of total covid deaths.
covid_y_justAR <- array(dim=c(n_sim,N,S))
y_justAR<-which(str_detect(colnames(covid_combined_samples_justAR),c('y')))
covid_y_justAR[,1:C,]<-covid_combined_samples_justAR[,y_justAR]%>%as.matrix()%>%array(dim=c(n_sim,C,S))
for(s in 1:S){
  covid_y_justAR[,(C+1):N,s] <- rnbinom(n_sim*(N-C),mu=covid_lambda_justAR[,(C+1):N,s],size=covid_theta_justAR[,s])
}
# Samples of the total for England.
covid_Y_justAR <- apply(covid_y_justAR,c(1,2),sum)
# Total cases for each region.
covid_y_justAR <- covid_y_justAR%>%apply(c(2,3),quantile,c(0.025,0.1,0.175,0.25,0.5,0.75,0.825,0.9,0.975))%>%melt(varnames=c('quantile','t','s'),value.name='y')%>%
  spread(quantile,y)%>%mutate(r=region_names[s])
covid_Y_justAR <- covid_Y_justAR%>%apply(2,quantile,c(0.025,0.1,0.175,0.25,0.5,0.75,0.825,0.9,0.975))%>%melt(varnames=c('quantile','t'),value.name='y')%>%
  spread(quantile,y)%>%mutate(s=S+1,r='England')
covid_y_all_justAR <- rbind(covid_y_justAR,covid_Y_justAR)%>%mutate(r=factor(r,levels=c(region_names,'England')))

# Check Residuals.
lambda_medians_justAR <- covid_lambda_justAR%>%apply(c(2,3),quantile,c(0.025,0.1,0.175,0.25,0.5,0.75,0.825,0.9,0.975))%>%melt(varnames=c('quantile','t','s'),value.name='y')%>%
  spread(quantile,y)%>%mutate(r=region_names[s])
#res_s<-1
scaled_residuals_justAR<-matrix(nrow=N,ncol=S)
for(res_s in 1:S){
  all_data<-inner_join(covid_y_all_justAR,y_data, by=c("r","t","s"))
  predicted_data<-filter(all_data, s==res_s)
  mu<-filter(lambda_medians_justAR, s==res_s)
  var_fn<-mu$`50%`+mu$`50%`^2/median(covid_theta_justAR[,res_s])
  residuals<-(predicted_data$y-mu$`50%`)
  scaled_residuals_justAR[,res_s]<-(predicted_data$y-mu$`50%`)/sqrt(var_fn)
  
  # Residuals over time.
  #plot(x=1:length(residuals), y=residuals)
  #plot(x=1:length(residuals), y=scaled_residuals)
  plot(x=1:length(residuals), y=scaled_residuals_justAR[,res_s], type="l")
  
  # Residual ACF
  # calculate the ACF for lags between 1 and 21 (inclusive).
  #  autocorrelation <- acf(scaled_residuals)
  # partial_autocorrelation <- pacf(scaled_residuals)
  
  # k_check<- gam(predicted_data$y~s(predicted_data$t, bs='cs' ,k =27))
  # gam.check(k_check)
}



# covid model with an AR term and a spline.
# Combine all MCMC chains.
covid_combined_samples_both <- as_tibble(do.call('rbind',covid_samples_both))

# Number of MCMC samples.
n_sim <- dim(covid_combined_samples_both)[1]
# Negative-Binomial dispersion parameters.
theta_both<-which(str_detect(colnames(covid_combined_samples_both),c('theta')))
covid_theta_both<-covid_combined_samples_both[,theta_both]%>%as.matrix()%>%array(dim=c(n_sim,S))
# Negative-Binomial means.
lambda_both<-which(str_detect(colnames(covid_combined_samples_both),c('lambda')))
covid_lambda_both<-covid_combined_samples_both[,lambda_both]%>%as.matrix()%>%array(dim=c(n_sim,N,S))
# Now-casting and forecasting samples of total covid deaths.
covid_y_both <- array(dim=c(n_sim,N,S))
y_both<-which(str_detect(colnames(covid_combined_samples_both),c('y')))
covid_y_both[,1:C,]<-covid_combined_samples_both[,y_both]%>%as.matrix()%>%array(dim=c(n_sim,C,S))
for(s in 1:S){
  covid_y_both[,(C+1):N,s] <- rnbinom(n_sim*(N-C),mu=covid_lambda_both[,(C+1):N,s],size=covid_theta_both[,s])
}
# Samples of the total for England.
covid_Y_both <- apply(covid_y_both,c(1,2),sum)
# Total cases for each region.
covid_y_both <- covid_y_both%>%apply(c(2,3),quantile,c(0.025,0.1,0.175,0.25,0.5,0.75,0.825,0.9,0.975))%>%melt(varnames=c('quantile','t','s'),value.name='y')%>%
  spread(quantile,y)%>%mutate(r=region_names[s])
covid_Y_both <- covid_Y_both%>%apply(2,quantile,c(0.025,0.1,0.175,0.25,0.5,0.75,0.825,0.9,0.975))%>%melt(varnames=c('quantile','t'),value.name='y')%>%
  spread(quantile,y)%>%mutate(s=S+1,r='England')
covid_y_all_both <- rbind(covid_y_both,covid_Y_both)%>%mutate(r=factor(r,levels=c(region_names,'England')))

# Check Residuals.
lambda_medians_both <- covid_lambda_both%>%apply(c(2,3),quantile,c(0.025,0.1,0.175,0.25,0.5,0.75,0.825,0.9,0.975))%>%melt(varnames=c('quantile','t','s'),value.name='y')%>%
  spread(quantile,y)%>%mutate(r=region_names[s])
#res_s<-1
scaled_residuals_both<-matrix(nrow=N,ncol=S)
for(res_s in 1:S){
  all_data<-inner_join(covid_y_all_both,y_data, by=c("r","t","s"))
  predicted_data<-filter(all_data, s==res_s)
  mu<-filter(lambda_medians_both, s==res_s)
  var_fn<-mu$`50%`+mu$`50%`^2/median(covid_theta_both[,res_s])
  residuals<-(predicted_data$y-mu$`50%`)
  scaled_residuals_both[,res_s]<-(predicted_data$y-mu$`50%`)/sqrt(var_fn)
  
  # Residuals over time.
  #plot(x=1:length(residuals), y=residuals)
  #plot(x=1:length(residuals), y=scaled_residuals)
  plot(x=1:length(residuals), y=scaled_residuals_both[,res_s], type="l")
  
  # Residual ACF
  # calculate the ACF for lags between 1 and 21 (inclusive).
   # autocorrelation <- acf(scaled_residuals)
   # partial_autocorrelation <- pacf(scaled_residuals)
  
  # k_check<- gam(predicted_data$y~s(predicted_data$t, bs='cs' ,k =27))
  # gam.check(k_check)
}


All_Residuals<-rbind(melt(scaled_residuals_noAR, value.name = 'scaled_residuals', varnames=c('t','s'))%>%mutate(model='Cubic spline'),
                     melt(scaled_residuals_justAR, value.name = 'scaled_residuals', varnames=c('t','s'))%>%mutate(model='Autoregressive term'),
                     melt(scaled_residuals_both, value.name = 'scaled_residuals', varnames=c('t','s'))%>%mutate(model='Cubic spline & autoregressive term'))

residual_plot_AR<-ggplot(data=All_Residuals)+geom_line(aes(x=dates[t],y=scaled_residuals,colour=model),alpha=0.7)+
  theme_minimal()+ geom_vline(xintercept = dates[C+1], linetype="dashed")+
  scale_color_discrete(name=NULL)+
  scale_x_date(date_breaks = "2 week", date_labels =  "%d %b %Y")+
  facet_wrap(~region_names[s], scales='free_x')+
  labs(title = 'Scaled residuals of the expected mean GDM model predictions', subtitle = "English regional COVID-19 deaths",
       y="Scaled residuals",x="Date")+geom_hline(yintercept=0, colour="darkgrey")+
  theme(text = element_text(size=12), plot.title = element_text(size=16),plot.subtitle=element_text(size=14),
        legend.title = element_text(size=12),strip.text = element_text(size=14), legend.position = "bottom")
ggsave(residual_plot_AR,file='Plots/residual_plot_AR.pdf',width=10,height=9)       









