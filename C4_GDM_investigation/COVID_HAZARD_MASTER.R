################################
## Nowcasting and Forecasting  #
## COVID-19 and Other Diseases #
################################

# Code structure..

############
# Preamble #
############

# Set your working directory.
setwd("~/Hazard_code")

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
library(INLA)
library(rjags)
library(NobBS)
library(gganimate)

# Read in some necessary functions.
source('Functions.R')

# Run to get latest data: 
#source('data_code.R')
# Run to use last saved latest data:
load("regions_raw.Rdata")

# Number of Regions.
S<-dim(regions_raw)[3]-1 # Number of regions (excluding England).

# Create matrix with Date of Death (rows) and Date the Death was Announced (columns).
Len<-dim(regions_raw)[1]
regions_raw <- regions_raw[1:Len,1:Len,2:(S+1)]

L <- dim(regions_raw)[1] # Number of rows in the raw data.
N<- L-14 # Length of the time series including 7 days of forecasting.
C <- N-7 # Length of time series up to the present day in the data.

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

dates <- seq(as.Date("2020/4/2"),by=1,length=L) # Dates corresponding to the data.
days <- weekdays(dates)
days <- days%>%
  factor(levels=c('Monday','Tuesday','Wednesday','Thursday','Friday','Saturday','Sunday'))%>%
  as.numeric()
save(dates, file='dates.Rdata')

# Put the raw data into a data frame for plotting.
raw_data=tibble(t=dates,day=factor(c('Monday','Tuesday','Wednesday','Thursday','Friday','Saturday','Sunday')[days],
                                   levels=c('Monday','Tuesday','Wednesday','Thursday','Friday','Saturday','Sunday')),
                act=apply(regions_raw,1,sum,na.rm=T),ann=apply(regions_raw,2,sum,na.rm=T))
raw_data=rbind(mutate(raw_data,type='New Deaths'),mutate(raw_data,ann=cumsum(ann),act=cumsum(act),type='Cumulative'))%>%
  mutate(type=factor(type,levels=c('New Deaths','Cumulative')))

## Figure 2: actual deaths versus announced deaths.
new_deaths <- ggplot(filter(raw_data,type=='New Deaths'))+
  geom_smooth(aes(x=t,y=act,linetype='Trend'),method="gam",n=100,se=FALSE,colour='gray40',size=0.5)+
  geom_line(aes(x=t,y=act,linetype='Actual'))+
  geom_line(aes(x=(t+1),y=ann,linetype='Announced'))+
  geom_point(aes(x=t,y=ann,colour=day,shape=day),size=2)+
  labs(x='Date',y='Deaths',title='Hospital Deaths from COVID-19 in England')+
  theme_minimal()+
  scale_color_brewer(name=NULL,palette='Dark2')+
  scale_shape_manual(name=NULL,values=1:7)+
  scale_linetype_manual(values=c(3,2,1),name=NULL)+scale_x_date(limits=as.Date(c('2020-04-07',dates[N])))+
  theme(text = element_text(size=12), plot.title = element_text(size=16),plot.subtitle=element_text(size=14),
        legend.title = element_text(size=14),strip.text = element_text(size=14))
ggsave(new_deaths,filename='Plots/new_deaths_long.pdf',width=9,height=3.25)

# Vector containing region names.
region_names <- c('East of England','London','Midlands',
                  'North East and Yorkshire',
                  'North West','South East','South West')
save(region_names, file='region_names.Rdata')

# Store counts in arrays (y_full/z_full is uncensored and y/z is censored).
z_full <- aperm(regions,c(1,3,2))

z_england <- apply(z_full,c(1,3),sum)
z_data_england <- cbind(z_england[,1:4],apply(z_england[,5:dim(z_england)[2]],1,sum,na.rm=T))%>%melt(varnames=c('t','d'))%>%
  mutate(d=as.character(d))
z_data_england$d[z_data_england$d=='5']='5+'
z_data_england <- mutate(z_data_england,d=factor(z_data_england$d,levels=rev(c(as.character(1:7),'5+'))),t=dates[t])

save(z_full, file='z_full.Rdata')

# Total reported daily deaths after 14 days of delay.
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
ggsave(bar_plot,file='Plots/bar_plot_covid_long.pdf',width=4.5,height=3)

# Observed index for each region at each delay.
obs_index<-matrix(NA, nrow=S, ncol=D)
for(s in 1:S){
  for(d in 1:D){
    obs_index[s,d]<-which(is.na(z[,s,d])==TRUE)[1]-1
  }
}

# Rolling Prediction Experiment 
n_cores<-16 # Number of cores to run cluster's on. 
start_date<- as.Date("2020-09-03") # Date to start the rolling predictions experiments.
end_date<- as.Date("2021-03-22") # Date to end the rolling predictions experiments.
start_t<-as.numeric(difftime(start_date,min(dates)-1,units=c("days")))
end_t<-as.numeric(difftime(end_date,min(dates)-1,units=c("days")))
n_nowcasts <- 24 # Number of days to preform nowcasts in the rolling predictions experiments.
# Rolling day indices for data censoring and model fitting.
C_list <- floor(seq(from=start_t, to=end_t, length=n_nowcasts))

#######################
#       Part One:     #
#    GDM Model for    #
# the 4th of May 2020 #
#######################

# Set up the splines using jagam.
n_knots <- c(30,8) # Number of knots of the temporal and seasonal splines, respectively.
blank_data=data_frame(y=rnorm(N,0,1),t=1:N,w=(1:N)%%7)
blank_jagam=jagam(y~ti(t,bs='cs',k=n_knots[1])+ti(w,bs='cc',k=n_knots[2])+
                    ti(t,w, bs=c('cs','cc'), k=n_knots[c(1,2)]),data=blank_data,file='blank.jags',
                  knots=list(t=seq(1,C,length=n_knots[1]),w=seq(0,7,length=n_knots[2])))

# Register the Beta-Binomial as a distribution for NIMBLE (see Functions.R for dbetabin).
registerDistributions(list(dbetabin=list(
  BUGSdist='dbetabin(mu,phi,size)',discrete=TRUE)))


covid_code<- nimbleCode({
  for(s in 1:S){
    for(t in 1:N){
      # Mean total deaths.
      log(lambda[t,s]) <- iota[s] + delta[t,s] 
      # Relative proportions (Beta-Binomial means).
      logit(mu[t,s,1]) <-  gamma[t,s]
      for(d in 2:D){
        logit(mu[t,s,d]) <-  beta[t,s,d-1] + gamma[t,s]
      }
    }
    for(t in 1:C){
      # Model for total counts
      y[t,s] ~ dnegbin( theta[s]/(theta[s]+lambda[t,s]),theta[s])
      # Model for delayed counts
      z[t,s,1] ~ dbetabin(mu[t,s,1],phi[s,1],y[t,s])
    }
    for(d in 2:D){
      for(t in 1:obs_index[s,d]){
        z[t,s,d] ~ dbetabin(mu[t,s,d],phi[s,d],y[t,s]-sum(z[t,s,1:(d-1)]))
      }
    }
  }
  ## Overall Temporal spline effect
  Omega_alpha[1:K_t,1:K_t] <- S_t[1:K_t,1:K_t]*tau_alpha
  kappa_alpha[1:K_t] ~ dmnorm(zeros[1:K_t],Omega_alpha[1:K_t,1:K_t])
  alpha[1:N] <- X_t[1:N,1:K_t]%*%kappa_alpha[1:K_t]
  
  ## Regional effects
  for(s in 1:S){
    Omega_delta[1:K_t,1:K_t,s] <- S_t[1:K_t,1:K_t]*tau_delta[s]
    kappa_delta[1:K_t,s] ~ dmnorm(kappa_alpha[1:K_t],Omega_delta[1:K_t,1:K_t,s])
    delta[1:N,s] <- X_t[1:N,1:K_t]%*%kappa_delta[1:K_t,s]
    
    # Time and week day splines
    Omega_gamma_time[1:K_t,1:K_t,s] <- S_t[1:K_t,1:K_t]*tau_gamma_t[s]
    Omega_gamma_week[1:K_w,1:K_w,s] <- S_w[1:K_w,1:K_w]*tau_gamma_w[s]
    Omega_gamma_inter[1:K_tw,1:K_tw,s] <- S_tw[1:K_tw,1:K_tw]*tau_gamma_tw[s,1]+S_tw[1:K_tw,(K_tw+1):(2*K_tw)]*tau_gamma_tw[s,2]
    
    kappa_gamma[1:K_t,s] ~ dmnorm(zeros[1:K_t],Omega_gamma_time[1:K_t,1:K_t,s])
    kappa_gamma[(K_t+1):(K_t+K_w),s] ~ dmnorm(zeros[1:K_w],Omega_gamma_week[1:K_w,1:K_w,s])
    kappa_gamma[(K_t+K_w+1):(K_t+K_w+K_tw),s] ~ dmnorm(zeros[1:K_tw],Omega_gamma_inter[1:K_tw,1:K_tw,s])
    
    gamma_time[1:(N),s] <-  gamma0[s] + X_tw[1:(N),1:(K_t)]%*%kappa_gamma[1:(K_t),s]
    gamma_week[1:(N),s] <- X_tw[1:(N),(K_t+1):(K_t+K_w)]%*%kappa_gamma[(K_t+1):(K_t+K_w),s]
    gamma_inter[1:(N),s] <-   X_tw[1:(N),(K_t+K_w+1):(K_t+K_w+K_tw)]%*%kappa_gamma[(K_t+K_w+1):(K_t+K_w+K_tw),s]
    
    gamma[1:N,s]<-matrix(gamma_time[1:(N),s]+gamma_week[1:(N),s]+gamma_inter[1:(N),s])
    
    # Delay splines
    for(d in 1:(D-1)){
      Omega_beta[1:K_t,1:K_t,s,d] <- S_t[1:K_t,1:K_t]*tau_beta[s,d]
      kappa_beta[1:K_t,s,d] ~ dmnorm(zeros[1:K_t],Omega_beta[1:K_t,1:K_t,s,d])
      beta[1:N,s,d] <- rep(psi[s,d],N) + X_t[1:N,1:K_t]%*%kappa_beta[1:K_t,s,d]
    }
  }
  
  # Smoothing parameter priors.
  tau_alpha ~ dinvgamma(0.5,0.5) # Equivalent to Half-Normal(0,1) on 1/sqrt(tau).
  for(s in 1:S){
    tau_delta[s] ~ dinvgamma(0.5,0.5)
    tau_gamma_t[s] ~ dinvgamma(0.5,0.5)
    tau_gamma_w[s] ~ dinvgamma(0.5,0.5)
    for(j in 1:2){
      tau_gamma_tw[s,j] ~ dinvgamma(0.5,0.5)
    }
    gamma0[s] ~ dnorm(0,sd=10) # Spatial intercept for time-week.
    for(d in 1:(D-1)){
      tau_beta[s,d]~ dinvgamma(0.5,0.5)
      psi[s,d] ~ dnorm(logit(1/(D+2-d)),sd=2)
    }
    for(d in 1:D){
      phi[s,d] ~ dgamma(2,0.02) # Beta-Binomial dispersion parameters.
    }
    iota[s] ~ dnorm(0,sd=10) # Spatial intercept.
    theta[s] ~ dgamma(2,0.02) # Negative-Binomial dispersion parameters.
  }
})


# Constants (e.g. number of days to model, number of knots) for NIMBLE.
covid_constants <- list(N=N,C=C,S=S,D=D,obs_index=obs_index,
                        K_t=n_knots[1]-1,
                        K_w=dim(blank_jagam$jags.data$S2)[1],
                        K_tw=dim(blank_jagam$jags.data$S3)[1])

# Data (e.g. total deaths, partial counts, spline model matrix) for NIMBLE.
covid_data <- list(z=z[1:C,,1:D],y=y[1:C,],
                   X_t=blank_jagam$jags.data$X[,2:(covid_constants$K_t+1)],
                   X_tw=blank_jagam$jags.data$X[,2:(covid_constants$K_t+covid_constants$K_w+covid_constants$K_tw+1)],
                   S_t=blank_jagam$jags.data$S1,
                   S_w=blank_jagam$jags.data$S2,
                   S_tw=blank_jagam$jags.data$S3,
                   zeros=rep(0,max(covid_constants$K_tw,covid_constants$K_t,covid_constants$K_w)))

# Set up the NIMBLE model for each chain.
covid_inits <- covid_model <- covid_compiled_model <- covid_mcmc_config <- covid_mcmc <- covid_compiled_mcmc <- list()
compile_time<-system.time({
  for(i in 1:n_chains){
    # Generate random initial values.
    covid_inits[[i]] <- list(kappa_alpha=rnorm(covid_constants$K_t,0,0.1),
                             kappa_delta=matrix(rnorm(S*covid_constants$K_t,0,0.1),ncol=S),
                             kappa_gamma=matrix(rnorm(S*(covid_constants$K_t+covid_constants$K_w+covid_constants$K_tw),0,0.1),ncol=S),
                             kappa_beta=array(rnorm(S*D*(covid_constants$K_t),0,0.01),dim=c(covid_constants$K_t,S,D)),
                             iota=rnorm(S,0,1),
                             theta=rexp(S,0.01),
                             gamma0=rnorm(S,-10,1),
                             psi=matrix(rnorm(S*D,0,1), nrow=S),
                             tau_alpha=rinvgamma(1,0.5,0.5),
                             tau_delta=rinvgamma(S,0.5,0.5),
                             tau_gamma_t=rinvgamma(S,5,0.5),
                             tau_gamma_w=rinvgamma(S,5,0.5),
                             tau_gamma_tw=matrix(rinvgamma(2*S,5,0.5), nrow=S),
                             tau_beta=matrix(rinvgamma(D*S,5,0.5), nrow=S),
                             phi=matrix(rexp(D*S,0.01),nrow=S,ncol=D),
                             y=matrix(NA,nrow=C,ncol=S),
                             z=array(NA,dim=c(C,S,D)))
    for(s in 1:S){
      for(t in 1:C){
        for(d in 1:D){
          if(is.na(covid_data$z[t,s,d])) covid_inits[[i]]$z[t,s,d]=rpois(1,median(covid_data$z[,s,d],na.rm=T))
        }
        if(is.na(covid_data$y[t,s])) covid_inits[[i]]$y[t,s]=sum(c(covid_inits[[i]]$z[t,s,],covid_data$z[t,s,],rpois(1,median(covid_data$y[,s]-rowSums(covid_data$z[,s,]),na.rm=T))),na.rm=TRUE)
      }
    }
    # Build the model.
    covid_model[[i]] <- nimbleModel(covid_code,covid_constants,covid_data,covid_inits[[i]])
    # Compile the model.
    covid_compiled_model[[i]] <- compileNimble(covid_model[[i]])
    # Set up the MCMC.
    covid_mcmc_config[[i]] <- configureMCMC(covid_model[[i]], #control = list(scale=0.1),
                                            monitors=c('alpha','iota','gamma0',
                                                       'psi', 'theta','phi',
                                                       'kappa_gamma','kappa_beta',
                                                       'delta','y'))
    for(s in 1:S){
      for(d in 1:(D-1)){
        covid_mcmc_config[[i]]$removeSamplers(paste('kappa_beta[1:',(covid_constants$K_t),', ',s,', ',d,']', sep=''))
        covid_mcmc_config[[i]]$addSampler(target=paste('kappa_beta[1:',(covid_constants$K_t),', ',s,', ',d,']', sep=''),type='ess') 
      }
      # covid_mcmc_config[[i]]$removeSamplers(paste('kappa_gamma[1:',(covid_constants$K_t),', ',s,']', sep=''))
      # covid_mcmc_config[[i]]$addSampler(target=paste('kappa_gamma[1:',(covid_constants$K_t),', ',s,']', sep=''),type='ess') 
      covid_mcmc_config[[i]]$removeSamplers(paste('kappa_gamma[',(covid_constants$K_t+1),':',(covid_constants$K_t+covid_constants$K_w),', ',s,']', sep=''))
      covid_mcmc_config[[i]]$addSampler(target=paste('kappa_gamma[',(covid_constants$K_t+1),':',(covid_constants$K_t+covid_constants$K_w),', ',s,']', sep=''),type='ess') 
      covid_mcmc_config[[i]]$removeSamplers(paste('kappa_gamma[',(covid_constants$K_t+covid_constants$K_w+1),':',(covid_constants$K_t+covid_constants$K_w+covid_constants$K_tw),', ',s,']', sep=''))
      covid_mcmc_config[[i]]$addSampler(target=paste('kappa_gamma[',(covid_constants$K_t+covid_constants$K_w+1),':',(covid_constants$K_t+covid_constants$K_w+covid_constants$K_tw),', ',s,']', sep=''),type='ess') 
      
    }
    covid_mcmc[[i]] <- buildMCMC(covid_mcmc_config[[i]])
    # Compile the MCMC.
    covid_compiled_mcmc[[i]] <- compileNimble(covid_mcmc[[i]],project=covid_model[[i]])
  }
})

# Run the model in parallel. Change 'dopar' to 'do' to run the model on only one core.
registerDoParallel(cores=n_chains)
run_time<-system.time({
  covid_samples <- as.mcmc.list(foreach(i=1:n_chains)%dopar%{
    runMCMC(covid_compiled_mcmc[[i]],niter=100000,nburnin=50000,inits=covid_inits[[i]],nchains=1,samplesAsCodaMCMC = TRUE,thin=10)
  })
})

save(compile_time,run_time,file="covid_time.RData")
save(covid_samples,file='covid_samples_hazard.RData')

# Combine all MCMC chains.
covid_combined_samples <- as_tibble(do.call('rbind',covid_samples))

# Number of MCMC samples.
n_sim <- dim(covid_combined_samples)[1]

## Check MCMC convergence for delta, and theta.
delta_index <- which(dimnames(covid_combined_samples)[[2]]=='delta[1, 1]'):which(dimnames(covid_combined_samples)[[2]]==paste('delta[',N,', ',S,']',sep=''))
delta_psrf <- sapply(delta_index,function(x)gelman.diag(covid_samples[,x],autoburnin = FALSE,transform = TRUE,multivariate = FALSE)$psrf[,1])

theta_index <- which(dimnames(covid_combined_samples)[[2]]=='theta[1]'):which(dimnames(covid_combined_samples)[[2]]==paste('theta[',S,']',sep=''))
theta_psrf <- sapply(theta_index,function(x)gelman.diag(covid_samples[,x],autoburnin = FALSE,transform = TRUE,multivariate = FALSE)$psrf[,1])

# Proportion of PSRFs below 1.05.
mean(delta_psrf<1.05)
mean(theta_psrf<1.05)

# Proportion of PSRFs below 1.2.
mean(delta_psrf<1.10)

# Overall temporal effect on covid fatality.
covid_alpha <- select(covid_combined_samples,contains('alpha'))%>%as.matrix()%>%apply(2,median)%>%melt(value.name='y')%>%mutate(x=1:N)

# Overall temporal effect on cumulative proportion reported.
covid_iota <- select(covid_combined_samples,starts_with('iota'))%>%as.matrix()%>%apply(2,median)
cbind(region_names,covid_iota)

# Regional means of total covid incidence.
covid_iota <- select(covid_combined_samples,starts_with('iota'))%>%as.matrix()%>%array(dim=c(n_sim,S))
covid_delta <- select(covid_combined_samples,starts_with('delta'))%>%as.matrix()%>%array(dim=c(n_sim,N,S))
covid_lambda<-array(dim=c(n_sim,N,S))
for(t in 1:N){
  for(s in 1:S){
    covid_lambda[,t,s]<-exp(covid_iota[,s]+covid_delta[,t,s])
  }
}

# Regional temporal effect of covid incidence.
covid_delta <- select(covid_combined_samples,starts_with('delta'))%>%as.matrix()%>%array(dim=c(n_sim,N,S))%>%
  apply(c(2,3),quantile,c(0.025,0.5,0.975))%>%melt(varnames=c('quantile','t','s'),value.name='y')%>%
  spread(quantile,y)%>%mutate(s=factor(region_names[s],levels=region_names))

# Time and week day splines
covid_kappagamma <- select(covid_combined_samples,starts_with('kappa_gamma'))%>%as.matrix()%>%array(dim=c(n_sim,covid_constants$K_t+covid_constants$K_w+covid_constants$K_tw,S))
covid_gamma0<- select(covid_combined_samples,starts_with('gamma0'))%>%as.matrix()%>%array(dim=c(n_sim,S))
gamma_time<-gamma_week<-gamma_inter<-array(dim=c(n_sim,N,S))
for(i in 1:n_sim){
  for(s in 1:S){
    # Time spline
    gamma_time[i,1:(N),s] <-  covid_gamma0[i,s] + covid_data$X_tw[1:(N),1:(covid_constants$K_t)]%*%covid_kappagamma[i,1:(covid_constants$K_t),s]
    # Day of week spline
    gamma_week[i,1:(N),s] <- covid_data$X_tw[1:(N),(covid_constants$K_t+1):(covid_constants$K_t+covid_constants$K_w)]%*%covid_kappagamma[i,(covid_constants$K_t+1):(covid_constants$K_t+covid_constants$K_w),s]
    # Time and day of week interaction
    gamma_inter[i,1:(N),s] <-   covid_data$X_tw[1:(N),(covid_constants$K_t+covid_constants$K_w+1):(covid_constants$K_t+covid_constants$K_w+covid_constants$K_tw)]%*%covid_kappagamma[i,(covid_constants$K_t+covid_constants$K_w+1):(covid_constants$K_t+covid_constants$K_w+covid_constants$K_tw),s]
  }}
# Overall time and day of week spline
gamma<-gamma_time+gamma_week+gamma_inter

covid_gamma <- gamma%>%
  apply(c(2,3),quantile,c(0.025,0.5,0.975))%>%melt(varnames=c('quantile','t','s'),value.name='y')%>%
  spread(quantile,y)%>%mutate(s=factor(region_names[s],levels=region_names),w=t%%7, w_d=weekdays(dates[t]))

covid_gamma_time <- gamma_time%>% 
  apply(c(2,3),quantile,c(0.025,0.5,0.975))%>%melt(varnames=c('quantile','t','s'),value.name='y')%>%
  spread(quantile,y)%>%mutate(s=factor(region_names[s],levels=region_names),w=t%%7,  w_d=weekdays(dates[t]))

covid_gamma_week <- gamma_week%>% 
  apply(c(2,3),quantile,c(0.025,0.5,0.975))%>%melt(varnames=c('quantile','t','s'),value.name='y')%>%
  spread(quantile,y)%>%mutate(s=factor(region_names[s],levels=region_names),w=t%%7,  w_d=weekdays(dates[t]))

covid_gamma_inter <- gamma_inter%>% 
  apply(c(2,3),quantile,c(0.025,0.5,0.975))%>%melt(varnames=c('quantile','t','s'),value.name='y')%>%
  spread(quantile,y)%>%mutate(s=factor(region_names[s],levels=region_names),w=t%%7, w_d=weekdays(dates[t]))

##Figure 3: Time effects Plots.
delta_plot <- ggplot(covid_delta)+
  geom_line(aes(colour=as.factor(s),x=dates[t],y=`50%`))+geom_line(data=covid_alpha,aes(x=dates[x],y=y),linetype=2)+theme_minimal()+
  labs(x=NULL,y=NULL,title='Temporal Trends',subtitle='Daily Death Rate')+
  scale_color_brewer(name=NULL,palette='Dark2',guide=guide_legend(nrow=7))+
  theme(axis.text=element_text(size=10),legend.text = element_text(size=10),legend.key.size = unit(0.75,'lines'),
        plot.title = element_text(size=16),plot.subtitle = element_text(size=14))
ggsave(delta_plot,file='Plots/delta_plot_hazard_long.pdf',width=9,height=3,scale=1)  

# Overall time and day of week effect plot
gamma_plot1 <- ggplot(covid_gamma)+
  geom_line(aes(colour=as.factor(s),x=dates[t],y=`50%`))+
  geom_point(aes(colour=as.factor(s),x=dates[t],y=`50%`,shape=as.factor(w)))+
  labs(x='Time (days)',y=NULL,title="Time and Week day effect on proportions of deaths reported")+
  theme_minimal()+
  facet_wrap(~as.factor(s),nrow=2,scales='free')+
  scale_color_brewer(name=NULL,palette='Dark2')+
  scale_fill_brewer(name=NULL,palette='Dark2')+
  scale_shape_manual(name=NULL,values=1:S)+
  theme(axis.text=element_text(size=10),legend.text = element_text(size=10),
        legend.key.size = unit(0.75,'lines'),
        plot.title = element_text(size=16),plot.subtitle = element_text(size=14))+
  scale_x_date(date_labels = '%d %b', limits = c(dates[N-50],dates[N]))
ggsave(gamma_plot1,file='Plots/gamma_plot.pdf',width=9,height=3)  

# Time effect plot
gamma_plot <- ggplot(covid_gamma_time)+
  geom_line(aes(colour=as.factor(s),x=dates[t],y=`50%`))+
  labs(x=NULL,y=NULL,title='',subtitle='Time effect on proportions of deaths reported')+
  guides(colour=FALSE)+theme_minimal()+
  scale_color_brewer(palette='Dark2',guide=FALSE)+
  theme(axis.text=element_text(size=10),plot.title = element_text(size=16),plot.subtitle = element_text(size=14))
ggsave(gamma_plot,file='Plots/gamma_time_plot.pdf',width=9,height=3,scale=1)  
# Day of week effect plot
gamma_plot_week <- ggplot(filter(covid_gamma_week, t%in%c(5:11)))+
  geom_line(aes(colour=as.factor(s),x=w,y=`50%`))+
  labs(x=NULL,y=NULL,title='',subtitle='Week day effect on proportions of deaths reported')+
  guides(colour=FALSE)+theme_minimal()+
  facet_wrap(~as.factor(s),nrow=2,scales='free')+
  scale_color_brewer(palette='Dark2',guide=FALSE)+
  scale_x_continuous(labels=c('Mon','Tue','Wed','Thu','Fri','Sat','Sun'),breaks=0:6)+
  theme(axis.text=element_text(size=10),plot.title = element_text(size=16),plot.subtitle = element_text(size=14))
ggsave(gamma_plot_week,file='Plots/gamma_week_plot.pdf',width=9,height=3)  
# Interaction of time and day of week effect plot
gamma_plot_inter <- ggplot(covid_gamma_inter)+
  geom_line(aes(colour=as.factor(s),x=dates[t],y=`50%`))+
  labs(x=NULL,y=NULL,title='',subtitle='Time and Week day interaction effect on proportions of deaths reported')+
  guides(colour=FALSE)+theme_minimal()+
  facet_wrap(~as.factor(s),nrow=2,scales='free')+
  scale_color_brewer(palette='Dark2',guide=FALSE)+
  theme(axis.text=element_text(size=10),plot.title = element_text(size=16),plot.subtitle = element_text(size=14))+
ggsave(gamma_plot_inter,file='Plots/gamma_inter_plot.pdf',width=20,height=3)  

# Animation of overall time and day of week effects by week 
start<-covid_gamma$t[which(covid_gamma$w_d=="Monday")[1]]
sundays_end<-which(covid_gamma$w_d=="Sunday")
end<-covid_gamma$t[sundays_end[length(sundays_end)]]
animate_gamma<-covid_gamma%>%filter(t>(start-1),t<(end+1))%>%mutate(t=t-start+1)%>%mutate(week=ceiling(t/7))

gamma_animate <- ggplot(animate_gamma)+
  geom_point(aes(colour=as.factor(s),x=w,y=`50%`))+
  geom_line(aes(colour=as.factor(s),x=w,y=`50%`))+
  theme_minimal()+
  geom_ribbon(aes(x=w,ymin=`2.5%`,ymax=`97.5%`,fill=s),alpha=0.25)+
  facet_wrap(~as.factor(s),nrow=2,scales='free')+
  scale_color_brewer(name=NULL,palette='Dark2')+
  scale_fill_brewer(name=NULL,palette='Dark2')+
  scale_shape_manual(name=NULL,values=1:S)+
  theme(axis.text=element_text(size=10),legend.text = element_text(size=10),
        legend.key.size = unit(0.75,'lines'),
        plot.title = element_text(size=16),plot.subtitle = element_text(size=14))+
  scale_x_continuous(labels=c('Mon','Tue','Wed','Thu','Fri','Sat','Sun'),breaks=0:6)
plot_animate<-gamma_animate + transition_time(week) +
  labs(title = "Week: {frame_time}", y=NULL)
gganim<-animate(plot_animate, height=365,width=1000,fps=30,duration=45,end_pause = 50,res=100,rewind=F,
                                renderer = gifski_renderer())
anim_save("Plots/gamma_animate.gif",animation=gganim)

# Animation of day of week effect and time and day of week interaction by week 
animate_gamma_wi<-array(dim=c(n_sim,N,S))
for(t in 1:N){
  animate_gamma_wi[,t,]<-gamma_time[,1,]+gamma_week[,t,]+gamma_inter[,t,]
}
gamma_wi<-animate_gamma_wi%>% 
  apply(c(2,3),quantile,c(0.025,0.5,0.975))%>%melt(varnames=c('quantile','t','s'),value.name='y')%>%
  spread(quantile,y)%>%mutate(s=factor(region_names[s],levels=region_names),w=t%%7, w_d=weekdays(dates[t]))%>%
  filter(t>(start-1),t<(end+1))%>%mutate(t=t-start+1)%>%mutate(week=ceiling(t/7))

gamma_animate_wi <- ggplot(gamma_wi)+
  geom_point(aes(colour=as.factor(s),x=w,y=`50%`))+
  geom_line(aes(colour=as.factor(s),x=w,y=`50%`))+
  theme_minimal()+
  geom_ribbon(aes(x=w,ymin=`2.5%`,ymax=`97.5%`,fill=s),alpha=0.25)+
  facet_wrap(~as.factor(s),nrow=2,scales='free')+
  scale_color_brewer(name=NULL,palette='Dark2')+
  scale_fill_brewer(name=NULL,palette='Dark2')+
  scale_shape_manual(name=NULL,values=1:S)+
  theme(axis.text=element_text(size=10),legend.text = element_text(size=10),
        legend.key.size = unit(0.75,'lines'),
        plot.title = element_text(size=16),plot.subtitle = element_text(size=14))+
  scale_x_continuous(labels=c('Mon','Tue','Wed','Thu','Fri','Sat','Sun'),breaks=0:6)
plot_animate_wi<-gamma_animate_wi + transition_time(week) +
  labs(title = "Week: {frame_time}", y=NULL)
gganim_wi<-animate(plot_animate_wi, height=365,width=1000,fps=30,duration=45,end_pause = 50,res=100,rewind=F,
                  renderer = gifski_renderer())
anim_save("Plots/gamma_week_inter_animate.gif",animation=gganim_wi)

# Effect of time and delay on relative proportions (Beta-Binomial means).
covid_psi<-select(covid_combined_samples,starts_with('psi'))%>%as.matrix()%>%array(dim=c(n_sim,S,D-1))
covid_kappabeta<-select(covid_combined_samples,starts_with('kappa_beta'))%>%as.matrix()%>%array(dim=c(n_sim,covid_constants$K_t,S,D-1))
beta<-array(dim=c(n_sim,N,S,D-1))
for(i in 1:n_sim){
  for(d in 1:(D-1)){
    for(s in 1:S){
      beta[i,1:N,s,d] <- rep(covid_psi[i,s,d],N) + covid_data$X_t[1:N,1:covid_constants$K_t]%*%covid_kappabeta[i,1:covid_constants$K_t,s,d]
    }
  }
}
covid_beta<-beta%>% 
  apply(c(2,3,4),quantile,c(0.025,0.5,0.975))%>%melt(varnames=c('quantile','t','s','d'),value.name='y')%>%
  spread(quantile,y)%>%mutate(s=factor(region_names[s],levels=region_names),w=t%%7, w_d=weekdays(dates[t]))

# Delay effect plot for given times.
beta_plot <- ggplot(filter(covid_beta, t%in%c(10,200,300,400)))+
  geom_line(aes(colour=as.factor(dates[t]),x=d,y=`50%`))+
  geom_point(aes(colour=as.factor(dates[t]),x=d,y=`50%`))+
  labs(x='Time (days)',title='Time and delay effects on proportion of deaths reported', y=NULL)+
  theme_minimal()+
  facet_wrap(~as.factor(s),nrow=2,scales='free')+
  scale_color_brewer(name=NULL,palette='Dark2')+
  scale_fill_brewer(name=NULL,palette='Dark2')+
  scale_shape_manual(name=NULL,values=1:S)+
  theme(axis.text=element_text(size=10),legend.text = element_text(size=10),
        legend.key.size = unit(0.75,'lines'),
        plot.title = element_text(size=16),plot.subtitle = element_text(size=14))
ggsave(beta_plot,file='Plots/beta_plot.pdf',width=9.5,height=3)  

# Animation of effect of delay over time (days).
beta_animate <- ggplot(covid_beta)+
  geom_line(aes(colour=as.factor(s),x=d,y=(`50%`)))+
  geom_point(aes(colour=as.factor(s),x=d,y=(`50%`)))+
  labs(x='Time (days)',title='Time and delay effects on proportion of deaths reported', y=NULL)+
  theme_minimal()+
  facet_wrap(~as.factor(s),nrow=2,scales='free')+
  scale_color_brewer(name=NULL,palette='Dark2')+
  scale_fill_brewer(name=NULL,palette='Dark2')+
  scale_shape_manual(name=NULL,values=1:S)+
  theme(axis.text=element_text(size=10),legend.text = element_text(size=10),
        legend.key.size = unit(0.75,'lines'),
        plot.title = element_text(size=16),plot.subtitle = element_text(size=14))

plot_beta_animate<-beta_animate + transition_time(t) +
  labs(title = "Time (days): {frame_time}", y=NULL)
gganim_beta<-animate(plot_beta_animate, height=365,width=1000,fps=30,duration=50,end_pause = 50,res=100,rewind=F,
                    renderer = gifski_renderer())
anim_save("Plots/beta_animate.gif",animation=gganim_beta)

# Regional time and delay effect on cumulative proportion reported.
expit_beta<-expit(beta)
covid_c_beta<-array(dim=c(n_sim,N,S,D-1))
covid_c_beta[,,,1]<-expit_beta[,,,1]
for(d in 2:(D-1)){
  covid_c_beta[,,,d]<-covid_c_beta[,,,d-1]+expit_beta[,,,d]*(1-covid_c_beta[,,,d-1])
}

covid_c_beta_plot<-covid_c_beta%>% 
  apply(c(2,3,4),quantile,c(0.025,0.5,0.975))%>%melt(varnames=c('quantile','t','s','d'),value.name='y')%>%
  spread(quantile,y)%>%mutate(s=factor(region_names[s],levels=region_names),w=t%%7, w_d=weekdays(dates[t]))

# Plot of delay effects on cumulative proportions reported. 
c_beta_plot <- ggplot(filter(covid_c_beta_plot, t%in%c(10,200,300,400)))+
  geom_line(aes(colour=as.factor(dates[t]),x=d,y=`50%`))+
  geom_point(aes(colour=as.factor(dates[t]),x=d,y=`50%`))+
  labs(x='Time (days)',title='Time and delay effects on cumulative proportion of deaths reported', y=NULL)+
  theme_minimal()+
  facet_wrap(~as.factor(s),nrow=2,scales='free')+
  scale_color_brewer(name=NULL,palette='Dark2')+
  scale_fill_brewer(name=NULL,palette='Dark2')+
  scale_shape_manual(name=NULL,values=1:S)+
  theme(axis.text=element_text(size=10),legend.text = element_text(size=10),
        legend.key.size = unit(0.75,'lines'),
        plot.title = element_text(size=16),plot.subtitle = element_text(size=14))
ggsave(c_beta_plot,file='Plots/cumulative_beta_plot.pdf',width=9.5,height=3)  

# Animation of delay effect on cumulatie proportions reported over time (days).
c_beta_animate <- ggplot(covid_c_beta_plot)+
  geom_line(aes(colour=as.factor(s),x=d,y=(`50%`)))+
  geom_point(aes(colour=as.factor(s),x=d,y=(`50%`)))+
  labs(x='Time (days)',title='Time and delay effects on cumulative proportion of deaths reported',y=NULL)+
  theme_minimal()+
  facet_wrap(~as.factor(s),nrow=2,scales='free')+
  scale_color_brewer(name=NULL,palette='Dark2')+
  scale_fill_brewer(name=NULL,palette='Dark2')+
  scale_shape_manual(name=NULL,values=1:S)+
  theme(axis.text=element_text(size=10),legend.text = element_text(size=10),
        legend.key.size = unit(0.75,'lines'),
        plot.title = element_text(size=16),plot.subtitle = element_text(size=14))

c_plot_beta_animate<-c_beta_animate + transition_time(dates[t]) +
  labs(title = "Time (days): {frame_time}", y=NULL)
c_gganim_beta<-animate(c_plot_beta_animate, height=365,width=1000,fps=30,duration=45,end_pause = 50,res=100,rewind=F,
                       renderer = gifski_renderer())
anim_save("Plots/cumulative_beta_animate.gif",animation=c_gganim_beta)

# Negative-Binomial dispersion parameters.
covid_theta <- select(covid_combined_samples,starts_with('theta'))%>%as.matrix()

# Now-casting and forecasting samples of total covid deaths.
covid_y <- array(dim=c(n_sim,N,S))
covid_y[,1:C,] <- select(covid_combined_samples,starts_with('y'))%>%as.matrix()%>%array(dim=c(n_sim,C,S))
for(s in 1:S){
  covid_y[,(C+1):N,s] <- rnbinom(n_sim*(N-C),mu=covid_lambda[,(C+1):N,s],size=covid_theta[,s])
}

# Samples of the total for England.
covid_Y <- apply(covid_y,c(1,2),sum)

# Total cases for each region.
covid_y <- covid_y%>%apply(c(2,3),quantile,c(0.025,0.1,0.175,0.25,0.5,0.75,0.825,0.9,0.975))%>%melt(varnames=c('quantile','t','s'),value.name='y')%>%
  spread(quantile,y)%>%mutate(r=region_names[s])

covid_Y <- covid_Y%>%apply(2,quantile,c(0.025,0.1,0.175,0.25,0.5,0.75,0.825,0.9,0.975))%>%melt(varnames=c('quantile','t'),value.name='y')%>%
  spread(quantile,y)%>%mutate(s=S+1,r='England')

covid_y <- rbind(covid_y,covid_Y)%>%mutate(r=factor(r,levels=c(region_names,'England')))

# Total cases data into long format.
y_data <- y_full[1:N,]%>%melt(varnames=c('t','s'),value.name='y')%>%mutate(r=region_names[s])
Y_data <- tibble(t=1:N,y=apply(y_full[1:N,],1,sum))%>%mutate(s=S+1,r='England')

y_data <- rbind(y_data,Y_data)%>%mutate(r=factor(r,levels=c(region_names,'England')))

y_data_so_far <- y_so_far[1:N,]%>%melt(varnames=c('t','s'),value.name='y')%>%mutate(r=region_names[s])
Y_data_so_far <- tibble(t=1:N,y=apply(y_so_far[1:N,],1,sum))%>%mutate(s=S+1,r='England')

y_data_so_far <- rbind(y_data_so_far,Y_data_so_far)%>%mutate(r=factor(r,levels=c(region_names,'England')))

# Figure 6: Plots of nowcasting and forecasting predictions for each region.
forecast_plot<- ggplot(filter(covid_y, t>470))+
  geom_ribbon(aes(x=dates[t],ymin=`2.5%`,ymax=`97.5%`,fill=as.factor(s)),alpha=0.2)+
  geom_ribbon(aes(x=dates[t],ymin=`10%`,ymax=`90%`,fill=as.factor(s)),alpha=0.2)+
  geom_ribbon(aes(x=dates[t],ymin=`17.5%`,ymax=`82.5%`,fill=as.factor(s)),alpha=0.2)+
  geom_ribbon(aes(x=dates[t],ymin=`25%`,ymax=`75%`,fill=as.factor(s)),alpha=0.2)+
  geom_line(aes(x=dates[t],y=`50%`,colour=as.factor(s)))+
  geom_point(data=filter(y_data, t>470),aes(x=dates[t],y=y,colour=as.factor(s)), size=1.5)+
  facet_wrap(~r,nrow=2,scales='free')+
  labs(x='Date of Death',y=NULL,title='Predicted Daily Hospital Deaths from COVID-19')+
  guides(colour=FALSE,fill=FALSE)+
  coord_cartesian()+
  geom_vline(xintercept=dates[C+1], alpha=0.6, linetype="longdash")+
  theme_minimal()+
  scale_fill_brewer(type = "qual",palette=2)+
  scale_colour_brewer(type = "qual",palette=2)+
  theme(strip.text = element_text(size=12),
        plot.caption = element_text(size=11, hjust=0),plot.title=element_text(size=16),
        axis.text = element_text(size=10),axis.title = element_text(size=12))+
  scale_x_date(date_labels = '%d %b ', limits=c(as.Date.character("2021-07-16"),as.Date.character("2021-09-02")))
ggsave(forecast_plot,file='Plots/forecast_hazard.pdf',width=11,height=4.5)

# Model Checking 
# Regional temporal effect of Covid incidence.
covid_iota <- select(covid_combined_samples,starts_with('iota'))%>%as.matrix()%>%array(dim=c(n_sim,S))
covid_delta <- select(covid_combined_samples,starts_with('delta'))%>%as.matrix()%>%array(dim=c(n_sim,N,S))
covid_lambda<-array(dim=c(n_sim,N,S))
for(t in 1:N){
  covid_lambda[,t,]<-exp(covid_iota[]+covid_delta[,t,])
}
covid_theta <- select(covid_combined_samples,starts_with('theta'))%>%as.matrix()

# Monte Carlo simulations of total covid deaths (y).
sim_y<-array(dim=c(n_sim,C,S))
for(t in 1:C){
  for(s in 1:S){
    sim_y[,t,s] <- rnbinom(n_sim,mu=covid_lambda[,t,s],size=covid_theta[,s])
  }
}
q_d<-91 # Average days in yearly quarter.
K=floor((which(is.na(y))[1]-1)/(q_d)) # Number of yearly quarters observed.
quarter_names<-c('Quarter 1')
for(k in 2:K){
  quarter_names<-c(quarter_names, paste("Quater ", k))
}
quarter_names<-c(quarter_names, 'Full time series')

# Simulated variance/mean.
var_mean_sim<-array(dim=c(n_sim,K+1,S))
for(i in 1:n_sim){
  for(s in 1:S){
    for(k in 1:K){
      var_mean_sim[i,k,s]<-var(sim_y[i,((k-1)*q_d+1):(k*q_d),s])/mean(sim_y[i,((k-1)*q_d+1):(k*q_d),s])
    }
    var_mean_sim[i,K+1,s]<-var(sim_y[i,1:(K*q_d),s])/mean(sim_y[i,1:(K*q_d),s])
  }
}
var_mean_melt_sim<-melt(var_mean_sim, varnames = c('iter','quarter','region'), value.name = "sim")%>%
  mutate(all=quarter==(K+1))

# Observed variance/mean.
var_mean_data<-array(dim=c(K+1,S))
for(s in 1:S){
  for(k in 1:K){
    var_mean_data[k,s]<-var(y[((k-1)*q_d+1):(k*q_d),s])/mean(y[((k-1)*q_d+1):(k*q_d),s])
  }
  var_mean_data[K+1,s]<-var(y[1:(K*q_d),s])/mean(y[1:(K*q_d),s])
}
var_mean_melt_data<-melt(var_mean_data, varnames = c('quarter','region'), value.name = "data")%>%
  mutate(all=quarter==(K+1))

var_mean<-full_join(var_mean_melt_data, var_mean_melt_sim, by=c('quarter','region'))
var_mean_plot<-var_mean%>%group_by(region,quarter)%>%summarise(median=median(sim),
                                                               lower=quantile(sim,0.025), 
                                                               upper=quantile(sim,0.975), 
                                                               data=median(data))
# Plot of model compared to data Variance/Mean.
sample_var_mean<-ggplot(var_mean_plot)+
  geom_point(aes(x=data, y=median, col=region_names[region], shape=region_names[region]))+
  facet_wrap(~quarter_names[quarter], scale="free")+
  theme_minimal()+geom_abline()+
  geom_errorbar(aes(x=data, ymin=lower, ymax=upper, col=region_names[region]))+
  labs(x='Sample Variance/Mean',y="Simulated Variance/Mean",title=NULL)+
  scale_color_discrete(name = "Region")+
  scale_shape_manual(values=1:7,name = "Region")
ggsave(sample_var_mean,file='Plots/variance_mean_check.pdf',width=11,height=4.5)

# Simulated Mean.
mean_sim<-array(dim=c(n_sim,K+1,S))
for(i in 1:n_sim){
  for(s in 1:S){
    for(k in 1:K){
      mean_sim[i,k,s]<-mean(sim_y[i,((k-1)*q_d+1):(k*q_d),s])
    }
    mean_sim[i,K+1,s]<-mean(sim_y[i,1:(K*q_d),s])
  }
}
mean_melt_sim<-melt(mean_sim, varnames = c('iter','quarter','region'), value.name = "sim")%>%
  mutate(all=quarter==(K+1))

# Observed Mean.
mean_data<-array(dim=c(K+1,S))
for(s in 1:S){
  for(k in 1:K){
    mean_data[k,s]<-mean(y[((k-1)*q_d+1):(k*q_d),s])
  }
  mean_data[K+1,s]<-mean(y[1:(K*q_d),s])
}
mean_melt_data<-melt(mean_data, varnames = c('quarter','region'), value.name = "data")%>%
  mutate(all=quarter==(K+1))

mean<-full_join(mean_melt_data, mean_melt_sim, by=c('quarter','region'))
mean_plot<-mean%>%group_by(region,quarter)%>%summarise(median=median(sim),
                                                       lower=quantile(sim,0.025), 
                                                       upper=quantile(sim,0.975), 
                                                       data=median(data))
# Plot of model compared to data Mean.
sample_mean<-ggplot(mean_plot)+
  geom_point(aes(x=data, y=median, col=region_names[region], shape=region_names[region]))+
  facet_wrap(~quarter_names[quarter], scale="free")+
  theme_minimal()+geom_abline()+
  geom_errorbar(aes(x=data, ymin=lower, ymax=upper, col=region_names[region]))+
  labs(x='Sample Mean',y="Simulated Mean",title=NULL)+
  scale_color_discrete(name = "Region")+
  scale_shape_manual(values=1:7,name = "Region")
ggsave(sample_mean,file='Plots/mean_check.pdf',width=11,height=4.5)


# Simulated Variance.
var_sim<-array(dim=c(n_sim,K+1,S))
for(i in 1:n_sim){
  for(s in 1:S){
    for(k in 1:K){
      var_sim[i,k,s]<-var(sim_y[i,((k-1)*q_d+1):(k*q_d),s])
    }
    var_sim[i,K+1,s]<-var(sim_y[i,1:(K*q_d),s])
  }
}
var_melt_sim<-melt(var_sim, varnames = c('iter','quarter','region'), value.name = "sim")%>%
  mutate(all=quarter==(K+1))

# Observed Variance.
var_data<-array(dim=c(K+1,S))
for(s in 1:S){
  for(k in 1:K){
    var_data[k,s]<-var(y[((k-1)*q_d+1):(k*q_d),s])
  }
  var_data[K+1,s]<-var(y[1:(K*q_d),s])
}
var_melt_data<-melt(var_data, varnames = c('quarter','region'), value.name = "data")%>%
  mutate(all=quarter==(K+1))

var<-full_join(var_melt_data, var_melt_sim, by=c('quarter','region'))
var_plot<-var%>%group_by(region,quarter)%>%summarise(median=median(sim),
                                                     lower=quantile(sim,0.025), 
                                                     upper=quantile(sim,0.975), 
                                                     data=median(data))
# Plot of model compared to data Variance.
sample_var<-ggplot(var_plot)+
  geom_point(aes(x=data, y=median, col=region_names[region], shape=region_names[region]))+
  facet_wrap(~quarter_names[quarter], scale="free")+
  theme_minimal()+geom_abline()+
  geom_errorbar(aes(x=data, ymin=lower, ymax=upper, col=region_names[region]))+
  labs(x='Sample Variance',y="Simulated Variance",title=NULL)+
  scale_color_discrete(name = "Region")+
  scale_shape_manual(values=1:7,name = "Region")
ggsave(sample_var,file='Plots/variance_check.pdf',width=11,height=4.5)

# Simulated IQR/Mean.
IQR_mean_sim<-array(dim=c(n_sim,K+1,S))
for(i in 1:n_sim){
  for(s in 1:S){
    for(k in 1:K){
      IQR_mean_sim[i,k,s]<-IQR(sim_y[i,((k-1)*q_d+1):(k*q_d),s])/mean(sim_y[i,((k-1)*q_d+1):(k*q_d),s])
    }
    IQR_mean_sim[i,K+1,s]<-IQR(sim_y[i,1:(K*q_d),s])/mean(sim_y[i,1:(K*q_d),s])
  }
}
IQR_mean_melt_sim<-melt(IQR_mean_sim, varnames = c('iter','quarter','region'), value.name = "sim")%>%
  mutate(all=quarter==(K+1))

# Observed IQR/Mean.
IQR_mean_data<-array(dim=c(K+1,S))
for(s in 1:S){
  for(k in 1:K){
    IQR_mean_data[k,s]<-IQR(y[((k-1)*q_d+1):(k*q_d),s])/mean(y[((k-1)*q_d+1):(k*q_d),s])
  }
  IQR_mean_data[K+1,s]<-IQR(y[1:(K*q_d),s])/mean(y[1:(K*q_d),s])
}
IQR_mean_melt_data<-melt(IQR_mean_data, varnames = c('quarter','region'), value.name = "data")%>%
  mutate(all=quarter==(K+1))

IQR_mean<-full_join(IQR_mean_melt_data, IQR_mean_melt_sim, by=c('quarter','region'))
IQR_mean_plot<-IQR_mean%>%group_by(region,quarter)%>%summarise(median=median(sim),
                                                               lower=quantile(sim,0.025), 
                                                               upper=quantile(sim,0.975), 
                                                               data=median(data))
# Plot of model compared to data IQR/Mean.
sample_IQR_mean<-ggplot(IQR_mean_plot)+
  geom_point(aes(x=data, y=median, col=region_names[region], shape=region_names[region]))+
  facet_wrap(~quarter_names[quarter], scale="free")+
  theme_minimal()+geom_abline()+
  geom_errorbar(aes(x=data, ymin=lower, ymax=upper, col=region_names[region]))+
  labs(x='Sample IQR/Mean',y="Simulated IQR/Mean",title=NULL)+
  scale_color_discrete(name = "Region")+
  scale_shape_manual(values=1:7,name = "Region")
ggsave(sample_IQR_mean,file='Plots/IQR_mean_check.pdf',width=11,height=4.5)



######################
#    Part Two:       #
# Rolling Prediction # 
#    Experiment      #
######################

# Make Cluster.
this_cluster<- makeCluster(n_cores)
N<-end_t+7 #forecast end date for rolling prediction experiment 

### cluster function
run_MCMC_allcode <- function( seed, C_list, n_knots, N, S, D, D_max, z) {
  
  #load libraries
  library(nimble)
  library(tidyverse)
  library(mgcv)
  
  # Register the Beta-Binomial as a distribution for NIMBLE (see Functions.R for dbetabin).
  dbetabin=nimbleFunction(run=function(x=double(0),mu=double(0),phi=double(0),size=double(0),log=integer(0)){
    returnType(double(0))
    if(x>=0&x<=size){
      return(lgamma(size+1)+lgamma(x+mu*phi)+lgamma(size-x+(1-mu)*phi)+lgamma(phi)-
               lgamma(size+phi)-lgamma(mu*phi)-lgamma((1-mu)*phi)-lgamma(size-x+1)-lgamma(x+1))
    }else{
      return(-Inf)
    }
  })
  
  rbetabin=nimbleFunction(run=function(n=integer(0),mu=double(0),phi=double(0),size=double(0)){
    pi=rbeta(1,mu*phi,(1-mu)*phi)
    returnType(double(0))
    return(rbinom(1,size,pi))
  })
  
  assign('dbetabin', dbetabin, envir = .GlobalEnv)
  assign('rbetabin', rbetabin, envir = .GlobalEnv)
  
  # Code for independent time series GDM models.
  covid_rolling_gdm_code <-  nimbleCode({
    for(s in 1:S){
      for(t in 1:N){
        # Mean total deaths.
        log(lambda[t,s]) <- iota[s] + delta[t,s] 
        # Relative proportions (Beta-Binomial means).
        logit(mu[t,s,1]) <-  gamma[t,s]
        for(d in 2:D){
          logit(mu[t,s,d]) <-  beta[t,s,d-1] + gamma[t,s]
        }
      }
      for(t in 1:C){
        # Model for total counts
        y[t,s] ~ dnegbin( theta[s]/(theta[s]+lambda[t,s]),theta[s])
        # Model for delayed counts
        z[t,s,1] ~ dbetabin(mu[t,s,1],phi[s,1],y[t,s])
      }
      for(d in 2:D){
        for(t in 1:obs_index[s,d]){
          z[t,s,d] ~ dbetabin(mu[t,s,d],phi[s,d],y[t,s]-sum(z[t,s,1:(d-1)]))
        }
      }
    }
    ## Regional effects
    for(s in 1:S){
      Omega_delta[1:K_t,1:K_t,s] <- S_t[1:K_t,1:K_t]*tau_delta[s]
      kappa_delta[1:K_t,s] ~ dmnorm(zeros[1:K_t],Omega_delta[1:K_t,1:K_t,s])
      delta[1:N,s] <- X_t[1:N,1:K_t]%*%kappa_delta[1:K_t,s]
      
      # Time and week day splines
      Omega_gamma_time[1:K_t,1:K_t,s] <- S_t[1:K_t,1:K_t]*tau_gamma_t[s]
      Omega_gamma_week[1:K_w,1:K_w,s] <- S_w[1:K_w,1:K_w]*tau_gamma_w[s]
      Omega_gamma_inter[1:K_tw,1:K_tw,s] <- S_tw[1:K_tw,1:K_tw]*tau_gamma_tw[s,1]+S_tw[1:K_tw,(K_tw+1):(2*K_tw)]*tau_gamma_tw[s,2]
      
      kappa_gamma[1:K_t,s] ~ dmnorm(zeros[1:K_t],Omega_gamma_time[1:K_t,1:K_t,s])
      kappa_gamma[(K_t+1):(K_t+K_w),s] ~ dmnorm(zeros[1:K_w],Omega_gamma_week[1:K_w,1:K_w,s])
      kappa_gamma[(K_t+K_w+1):(K_t+K_w+K_tw),s] ~ dmnorm(zeros[1:K_tw],Omega_gamma_inter[1:K_tw,1:K_tw,s])
      
      gamma_time[1:(N),s] <-  gamma0[s] + X_tw[1:(N),1:(K_t)]%*%kappa_gamma[1:(K_t),s]
      gamma_week[1:(N),s] <- X_tw[1:(N),(K_t+1):(K_t+K_w)]%*%kappa_gamma[(K_t+1):(K_t+K_w),s]
      gamma_inter[1:(N),s] <-   X_tw[1:(N),(K_t+K_w+1):(K_t+K_w+K_tw)]%*%kappa_gamma[(K_t+K_w+1):(K_t+K_w+K_tw),s]
      
      gamma[1:N,s]<-matrix(gamma_time[1:(N),s]+gamma_week[1:(N),s]+gamma_inter[1:(N),s])
      
      # Delay splines
      for(d in 1:(D-1)){
        Omega_beta[1:K_t,1:K_t,s,d] <- S_t[1:K_t,1:K_t]*tau_beta[s,d]
        kappa_beta[1:K_t,s,d] ~ dmnorm(zeros[1:K_t],Omega_beta[1:K_t,1:K_t,s,d])
        beta[1:N,s,d] <- rep(psi[s,d],N) + X_t[1:N,1:K_t]%*%kappa_beta[1:K_t,s,d]
      }
    }
    
    # Smoothing parameter priors.
    for(s in 1:S){
      tau_delta[s] ~ dinvgamma(0.5,0.5)
      tau_gamma_t[s] ~ dinvgamma(0.5,0.5)
      tau_gamma_w[s] ~ dinvgamma(0.5,0.5)
      for(j in 1:2){
        tau_gamma_tw[s,j] ~ dinvgamma(0.5,0.5)
      }
      gamma0[s] ~ dnorm(0,sd=10) # Spatial intercept for time-week.
      for(d in 1:(D-1)){
        tau_beta[s,d]~ dinvgamma(0.5,0.5)
        psi[s,d] ~ dnorm(logit(1/(D+2-d)),sd=2)
      }
      for(d in 1:D){
        phi[s,d] ~ dgamma(2,0.02) # Beta-Binomial dispersion parameters.
      }
      iota[s] ~ dnorm(0,sd=10) # Spatial intercept.
      theta[s] ~ dgamma(2,0.02) # Negative-Binomial dispersion parameters.
    }
  })
  
  covid_rolling_gdm_constants <- covid_rolling_gdm_data <- covid_rolling_gdm_inits <- covid_rolling_gdm_model <- covid_rolling_gdm_compiled_model <- covid_rolling_gdm_mcmc_config <- covid_rolling_gdm_mcmc <- covid_rolling_gdm_compiled_mcmc <- list()
  blank_data=data_frame(y=rnorm(N,0,1),t=1:N,w=(1:N)%%7)
  blank_jagam=jagam(y~ti(t,bs='cs',k=n_knots[1])+ti(w,bs='cc',k=n_knots[2])+
                      ti(t,w, bs=c('cs','cc'), k=n_knots[c(1,2)]),data=blank_data,file='blank.jags',
                    knots=list(t=seq(1,N,length=n_knots[1]),w=seq(0,7,length=n_knots[2])))
  
  
  
  # Censor the delayed counts according to
  # the data which would be available on day i.
  censored_z <- z
  for(s in 1:S){
    censored_z[,s,][outer(1:dim(z[,s,])[1], 0:(dim(z[,s,])[2]-1), FUN = "+") > C_list[seed]] <- NA
  }
  # Observed index for each region at each delay.
  obs_index<-matrix(NA, nrow=S, ncol=D)
  for(s in 1:S){
    for(d in 1:D){
      obs_index[s,d]<-which(is.na(censored_z[,s,d])==TRUE)[1]-1
    }
  }
  
  # Censor the totals.
  censored_y <- apply(censored_z[,,1:D_max],c(1,2),sum)
  
  # Constants (e.g. number of weeks to model, number of knots) for NIMBLE.
  covid_rolling_gdm_constants <-  list(N=N,C=C_list[seed],S=S,D=D,obs_index=obs_index,
                                       K_t=n_knots[1]-1,
                                       K_w=dim(blank_jagam$jags.data$S2)[1],
                                       K_tw=dim(blank_jagam$jags.data$S3)[1])
  
  # Data (e.g. total cases, partial counts, spline model matrix) for NIMBLE.
  covid_rolling_gdm_data <- list(z=censored_z[1:C_list[seed],,1:D],y=censored_y[1:C_list[seed],],
                                 X_t=blank_jagam$jags.data$X[,2:(covid_rolling_gdm_constants$K_t+1)],
                                 X_tw=blank_jagam$jags.data$X[,2:(covid_rolling_gdm_constants$K_t+covid_rolling_gdm_constants$K_w+covid_rolling_gdm_constants$K_tw+1)],
                                 S_t=blank_jagam$jags.data$S1,
                                 S_w=blank_jagam$jags.data$S2,
                                 S_tw=blank_jagam$jags.data$S3,
                                 zeros=rep(0,max(covid_rolling_gdm_constants$K_tw,covid_rolling_gdm_constants$K_t,covid_rolling_gdm_constants$K_w)))
  
  # Generate random initial values.
  covid_rolling_gdm_inits <- list(kappa_delta=matrix(rnorm(S*covid_rolling_gdm_constants$K_t,0,0.1),ncol=S),
                                  kappa_gamma=matrix(rnorm(S*(covid_rolling_gdm_constants$K_t+covid_rolling_gdm_constants$K_w+covid_rolling_gdm_constants$K_tw),0,0.1),ncol=S),
                                  kappa_beta=array(rnorm(S*D*(covid_rolling_gdm_constants$K_t),0,0.01),dim=c(covid_rolling_gdm_constants$K_t,S,D)),
                                  iota=rnorm(S,0,1),
                                  theta=rexp(S,0.01),
                                  gamma0=rnorm(S,-10,1),
                                  psi=matrix(rnorm(S*D,0,1), nrow=S),
                                  tau_delta=rinvgamma(S,0.5,0.5),
                                  tau_gamma_t=rinvgamma(S,5,0.5),
                                  tau_gamma_w=rinvgamma(S,5,0.5),
                                  tau_gamma_tw=matrix(rinvgamma(2*S,5,0.5), nrow=S),
                                  tau_beta=matrix(rinvgamma(D*S,5,0.5), nrow=S),
                                  phi=matrix(rexp(D*S,0.01),nrow=S,ncol=D),
                                  y=matrix(NA,nrow=C_list[seed],ncol=S),
                                  z=array(NA,dim=c(C_list[seed],S,D)))
  for(s in 1:S){
    for(t in 1:C_list[seed]){
      for(d in 1:D){
        if(is.na(covid_rolling_gdm_data$z[t,s,d])) covid_rolling_gdm_inits$z[t,s,d]=rpois(1,median(covid_rolling_gdm_data$z[,s,d],na.rm=T))
      }
      if(is.na(covid_rolling_gdm_data$y[t,s])) covid_rolling_gdm_inits$y[t,s]=sum(c(covid_rolling_gdm_inits$z[t,s,],
                                                                                    covid_rolling_gdm_data$z[t,s,],rpois(1,median(covid_rolling_gdm_data$y[,s]-
                                                                                                                                    rowSums(covid_rolling_gdm_data$z[,s,]),na.rm=T))),na.rm=TRUE)
    }
  }
  # Build the model.
  covid_rolling_gdm_model <- nimbleModel(covid_rolling_gdm_code,covid_rolling_gdm_constants,covid_rolling_gdm_data,covid_rolling_gdm_inits)
  # Compile the model.
  covid_rolling_gdm_compiled_model <- compileNimble(covid_rolling_gdm_model)
  # Set up the MCMC.
  covid_rolling_gdm_mcmc_config <- configureMCMC(covid_rolling_gdm_model,monitors=c('iota','gamma0','lambda','y',
                                                                                    'psi', 'theta','phi',
                                                                                    'kappa_gamma','kappa_beta',
                                                                                    'delta'))
  for(s in 1:S){
    for(d in 1:(D-1)){
      covid_rolling_gdm_mcmc_config$removeSamplers(paste('kappa_beta[1:',(covid_rolling_gdm_constants$K_t),', ',s,', ',d,']', sep=''))
      covid_rolling_gdm_mcmc_config$addSampler(target=paste('kappa_beta[1:',(covid_rolling_gdm_constants$K_t),', ',s,', ',d,']', sep=''),type='ess') 
    }
    # covid_rolling_gdm_mcmc_config$removeSamplers(paste('kappa_gamma[1:',(covid_rolling_gdm_constants$K_t),', ',s,']', sep=''))
    # covid_rolling_gdm_mcmc_config$addSampler(target=paste('kappa_gamma[1:',(covid_rolling_gdm_constants$K_t),', ',s,']', sep=''),type='ess') 
    covid_rolling_gdm_mcmc_config$removeSamplers(paste('kappa_gamma[',(covid_rolling_gdm_constants$K_t+1),':',(covid_rolling_gdm_constants$K_t+covid_rolling_gdm_constants$K_w),', ',s,']', sep=''))
    covid_rolling_gdm_mcmc_config$addSampler(target=paste('kappa_gamma[',(covid_rolling_gdm_constants$K_t+1),':',(covid_rolling_gdm_constants$K_t+covid_rolling_gdm_constants$K_w),', ',s,']', sep=''),type='ess') 
    covid_rolling_gdm_mcmc_config$removeSamplers(paste('kappa_gamma[',(covid_rolling_gdm_constants$K_t+covid_rolling_gdm_constants$K_w+1),':',(covid_rolling_gdm_constants$K_t+covid_rolling_gdm_constants$K_w+covid_rolling_gdm_constants$K_tw),', ',s,']', sep=''))
    covid_rolling_gdm_mcmc_config$addSampler(target=paste('kappa_gamma[',(covid_rolling_gdm_constants$K_t+covid_rolling_gdm_constants$K_w+1),':',(covid_rolling_gdm_constants$K_t+covid_rolling_gdm_constants$K_w+covid_rolling_gdm_constants$K_tw),', ',s,']', sep=''),type='ess') 
    
  }
  covid_rolling_gdm_mcmc <- buildMCMC(covid_rolling_gdm_mcmc_config)
  # Compile the MCMC.
  covid_rolling_gdm_compiled_mcmc <- compileNimble(covid_rolling_gdm_mcmc,project=covid_rolling_gdm_model)
  
  covid_samples_cluster <- runMCMC(covid_rolling_gdm_compiled_mcmc,niter=150000,nburnin=50000,inits=covid_rolling_gdm_inits,nchains=1,samplesAsCodaMCMC = TRUE,thin=10)
  
  return(covid_samples_cluster)
  
}


rolling_time<-system.time({
  covid_rolling_gdm_samples  <- parLapply(cl = this_cluster, X = 1:length(C_list), 
                                          fun = run_MCMC_allcode, 
                                          C_list=C_list, n_knots=n_knots, N=N, S=S, 
                                          D_max=D_max, D=D, z=z_full[1:N,,])
})

stopCluster(this_cluster)

save(rolling_time,file='rolling_time.RData' )
save(covid_rolling_gdm_samples,file='covid_rolling_gdm_samples.RData')

n_sim_censored <- dim(covid_rolling_gdm_samples[[1]])[1]

# Nowcasting and forecasting samples.
covid_rolling_gdm_y <- lapply(covid_rolling_gdm_samples,function(x){
  # Negative-Binomial dispersion parameters.
  theta <- select(as.data.frame(x),starts_with('theta'))%>%as.matrix()
  # Negative-Binomial means.
  lambda <- select(as.data.frame(x),starts_with('lambda'))%>%as.matrix()%>%array(dim=c(n_sim_censored,N,S))
  
  # Now-casting and forecasting samples of total covid cases.
  y <- array(dim=c(n_sim_censored,N,S))
  y_samples <- select(as.data.frame(x),starts_with('y'))%>%as.matrix()
  C <- dim(y_samples)[2]/S
  y[,1:C,] <- array(y_samples,dim=c(n_sim_censored,C,S))
  for(s in 1:S){
    y[,(C+1):N,s] <- rnbinom(n_sim_censored*(N-C),mu=lambda[,(C+1):N,s],size=theta[,s])
  }
  y <- abind(y,apply(y,c(1,2),sum),along=3)
  return(y)
})

covid_rolling_gdm_y <- do.call('abind',list(covid_rolling_gdm_y,along=4))

covid_rolling_gdm_y_cumulative <- apply(covid_rolling_gdm_y,c(1,3,4),cumsum)%>%aperm(c(2,1,3,4))%>%apply(c(2,3,4),quantile,c(0.025,0.5,0.975))%>%melt(varnames=c('quantile','t','s','c'),value.name='y')%>%
  spread(quantile,y)%>%mutate(r=c(region_names,'England')[s],c=dates[C_list[c]],delta=dates[t]-c)

covid_rolling_gdm_y <- apply(covid_rolling_gdm_y,c(2,3,4),quantile,c(0.025,0.5,0.975))%>%melt(varnames=c('quantile','t','s','c'),value.name='y')%>%
  spread(quantile,y)%>%mutate(r=c(region_names,'England')[s],c=dates[C_list[c]],delta=dates[t]-c)


# Put the censored total counts into long format.
censored_y_data <- array(dim=c(N,S+1,length(C_list)))
for(i in 1:length(C_list)){
  censored_z <- z[1:N,,]
  for(s in 1:S){
    censored_z[,s,][outer(1:dim(z[1:N,s,])[1], 0:(dim(z[1:N,s,])[2]-1), FUN = "+") > C_list[i]] <- NA
  }
  censored_y_data[,1:S,i] <- apply(censored_z,c(1,2),sum,na.rm=TRUE)
  censored_y_data[,S+1,i] <- apply(censored_y_data[,1:S,i],1,sum)
}  
censored_y_data <- censored_y_data%>%
  melt(varnames=c('t','s','c'),value.name='y')%>%
  mutate(r=factor(c(region_names,'England')[s],
                  levels=c(region_names,'England')),c=dates[C_list[c]])

# Total cases data into long format.
y_data <- y_full[1:N,]%>%melt(varnames=c('t','s'),value.name='y')%>%mutate(r=region_names[s])
Y_data <- tibble(t=1:N,y=apply(y_full[1:N,],1,sum))%>%mutate(s=S+1,r='England')
y_data <- rbind(y_data,Y_data)%>%mutate(r=factor(r,levels=c(region_names,'England')))

# Arrange predictions and data by 
# days relative to model fit date (delta).
covid_rolling_gdm_y <- mutate(covid_rolling_gdm_y,delta=dates[t]-c)
covid_rolling_gdm_y <- left_join(covid_rolling_gdm_y,y_data)%>%mutate(r=factor(r,levels=c(region_names,'England')))

delta_seq <- (-3):6
n_delta <- length(delta_seq)
coverage <- numeric(n_delta)
for(i in 1:n_delta){
  filtered <- filter(covid_rolling_gdm_y,delta==delta_seq[i])
  coverage[i] <- mean(filtered$y>=filtered$`2.5%`&filtered$y<=filtered$`97.5%`,na.rm=T)
}
coverage_data <- tibble(delta=delta_seq,coverage=coverage)


# Figure 7: Predictions by day relative to model fit date.  
coverage_plot <- ggplot(filter(covid_rolling_gdm_y,delta>=-3&delta<=6))+
  geom_abline(slope=1,intercept=0)+
  geom_errorbar(aes(x=y,ymin=`2.5%`,ymax=`97.5%`,colour=r))+
  geom_point(aes(x=y,y=`50%`,colour=r,shape=r))+
  geom_text(data=coverage_data,aes(x=10,y=316,label=round(coverage,2)))+
  facet_wrap(~factor(sapply(delta,function(x)paste(x,'Days')),levels=unique(sapply(delta,function(x)paste(x,'Days')))),nrow=2)+
  scale_fill_brewer(name=NULL,palette='Dark2')+scale_colour_brewer(name=NULL,palette='Dark2')+
  scale_shape_manual(name=NULL,values=1:8)+
  scale_x_continuous(trans='log10',breaks=c(10,100,1000))+
  scale_y_continuous(trans='log10',breaks=c(10,100,1000))+
  coord_cartesian(ylim=c(3.16,1000),xlim=c(3.16,1778))+
  theme_minimal()+labs(y='Predicted Deaths',x='Observed Deaths',title='Predictions of Hospital Deaths from COVID-19 in England',
                       subtitle='by Time from Day Model Fitted')+
  theme(strip.text = element_text(size=12),
        plot.subtitle=element_text(size=14),
        plot.caption = element_text(size=10),plot.title=element_text(size=16),
        axis.text = element_text(size=10),axis.title = element_text(size=12),
        legend.text = element_text(size=10))
ggsave(coverage_plot,width=9,height=4,filename='Plots/coverage_plot.pdf')

# Multiple Nowcasts plotted by region.
C_days<-C_list[floor(seq(from=1, to=length(unique(covid_rolling_gdm_y$c)), length=8))] 
reg_all_GDM<- ggplot()+
  geom_point(data=covid_rolling_gdm_y, aes(x = dates[t], y = y))+
  geom_ribbon(data=filter(covid_rolling_gdm_y,c==dates[C_days[8]], t<(C_days[8]+1)), 
              aes(x=dates[t],ymin=`2.5%`,ymax=`97.5%`,colour=as.factor(dates[C_days[8]]),fill=as.factor(dates[C_days[8]])),alpha=0.35)+
  geom_ribbon(data=filter(covid_rolling_gdm_y,c==dates[C_days[7]], t<(C_days[7]+1)), 
              aes(x=dates[t],ymin=`2.5%`,ymax=`97.5%`,colour=as.factor(dates[C_days[7]]),fill=as.factor(dates[C_days[7]])),alpha=0.35)+
  geom_ribbon(data=filter(covid_rolling_gdm_y,c==dates[C_days[6]], t<(C_days[6]+1)), 
              aes(x=dates[t],ymin=`2.5%`,ymax=`97.5%`,colour=as.factor(dates[C_days[6]]),fill=as.factor(dates[C_days[6]])),alpha=0.35)+
  geom_ribbon(data=filter(covid_rolling_gdm_y,c==dates[C_days[5]], t<(C_days[5]+1)), 
              aes(x=dates[t],ymin=`2.5%`,ymax=`97.5%`,colour=as.factor(dates[C_days[5]]),fill=as.factor(dates[C_days[5]])),alpha=0.35)+
  geom_ribbon(data=filter(covid_rolling_gdm_y,c==dates[C_days[4]], t<(C_days[4]+1)), 
              aes(x=dates[t],ymin=`2.5%`,ymax=`97.5%`,colour=as.factor(dates[C_days[4]]),fill=as.factor(dates[C_days[4]])),alpha=0.35)+
  geom_ribbon(data=filter(covid_rolling_gdm_y,c==dates[C_days[3]], t<(C_days[3]+1)), 
              aes(x=dates[t],ymin=`2.5%`,ymax=`97.5%`,colour=as.factor(dates[C_days[3]]),fill=as.factor(dates[C_days[3]])),alpha=0.35)+
  geom_ribbon(data=filter(covid_rolling_gdm_y,c==dates[C_days[2]], t<(C_days[2]+1)), 
              aes(x=dates[t],ymin=`2.5%`,ymax=`97.5%`,colour=as.factor(dates[C_days[2]]),fill=as.factor(dates[C_days[2]])),alpha=0.35)+
  geom_ribbon(data=filter(covid_rolling_gdm_y,c==dates[C_days[1]], t<(C_days[1]+1)), 
              aes(x=dates[t],ymin=`2.5%`,ymax=`97.5%`,colour=as.factor(dates[C_days[1]]),fill=as.factor(dates[C_days[1]])),alpha=0.35)+
  facet_wrap(~r,nrow=4,scales='free')+
  theme_minimal()+
  scale_color_brewer(name='Model Fit Date',palette='Spectral')+
  scale_fill_brewer(name='Model Fit Date',palette='Spectral')+
  scale_x_date(date_labels = '%b%y',limits=c(dates[start_t-14], dates[end_t+7]))+
  ylab("Number of  cases") + xlab("Time (days after 2020-04-02)")
reg_all_GDM
ggsave(reg_all_GDM,width=17,height=17,filename='GDM_rolling.pdf')



#####################
#    Part Three:    #
# Marginal NB Model #
#####################

## Rolling prediction experiment using the marginal model.
# Make Cluster.
this_cluster<- makeCluster(n_cores)
N<-end_t+7 #forecast end date for rolling prediction experiment 

run_MCMC_allcode <- function( seed, C_list, n_knots, N, S, D, D_max, z) {
  
  #load libraries
  library(nimble)
  library(tidyverse)
  library(mgcv)
  library(abind)
  
  # Register the Beta-Binomial as a distribution for NIMBLE (see Functions.R for dbetabin).
  dbetabin=nimbleFunction(run=function(x=double(0),mu=double(0),phi=double(0),size=double(0),log=integer(0)){
    returnType(double(0))
    if(x>=0&x<=size){
      return(lgamma(size+1)+lgamma(x+mu*phi)+lgamma(size-x+(1-mu)*phi)+lgamma(phi)-
               lgamma(size+phi)-lgamma(mu*phi)-lgamma((1-mu)*phi)-lgamma(size-x+1)-lgamma(x+1))
    }else{
      return(-Inf)
    }
  })
  
  rbetabin=nimbleFunction(run=function(n=integer(0),mu=double(0),phi=double(0),size=double(0)){
    pi=rbeta(1,mu*phi,(1-mu)*phi)
    returnType(double(0))
    return(rbinom(1,size,pi))
  })
  
  assign('dbetabin', dbetabin, envir = .GlobalEnv)
  assign('rbetabin', rbetabin, envir = .GlobalEnv)
  
  # Marginal NB model code. 
  covid_nb_code <- nimbleCode({
    for(s in 1:S){
      for(t in 1:N){
        # Mean total deaths.
        log(lambda[t,s]) <- iota[s] + delta[t,s] 
        # Relative proportions (Beta-Binomial means).
        logit(nu[t,s,1]) <-  gamma[t,s]
        for(d in 2:(D+1)){
          logit(nu[t,s,d]) <-  beta[t,s,d-1]
        }
        
        # Absolute proportions.
        mu[t,s,1] <- nu[t,s,1]
        # Cumulative proportions.
        p[t,s,1] <- nu[t,s,1]
        for(d in 2:(D+1)){
          p[t,s,d]<-p[t,s,d-1]+nu[t,s,d]*(1-p[t,s,d-1])
          mu[t,s,d] <- nu[t,s,d]*(1-p[t,s,d-1])
        }
        for(d in 1:(D+1)){
          mu_lambda[t,s,d] <- mu[t,s,d]*lambda[t,s]
        }
      }
      for(t in 1:C){
        # Model for partial reports.
        for(d in 1:(D+1)){
          z[t,s,d] ~ dnegbin(theta[s]/(theta[s]+mu_lambda[t,s,d]),theta[s])
        }
      }
    }
    ## Regional effects
    for(s in 1:S){
      Omega_delta[1:K_t,1:K_t,s] <- S_t[1:K_t,1:K_t]*tau_delta[s]
      kappa_delta[1:K_t,s] ~ dmnorm(zeros[1:K_t],Omega_delta[1:K_t,1:K_t,s])
      delta[1:N,s] <- X_t[1:N,1:K_t]%*%kappa_delta[1:K_t,s]
      
      # Time and week day splines
      Omega_gamma_time[1:K_t,1:K_t,s] <- S_t[1:K_t,1:K_t]*tau_gamma_t[s]
      Omega_gamma_week[1:K_w,1:K_w,s] <- S_w[1:K_w,1:K_w]*tau_gamma_w[s]
      Omega_gamma_inter[1:K_tw,1:K_tw,s] <- S_tw[1:K_tw,1:K_tw]*tau_gamma_tw[s,1]+S_tw[1:K_tw,(K_tw+1):(2*K_tw)]*tau_gamma_tw[s,2]
      
      kappa_gamma[1:K_t,s] ~ dmnorm(zeros[1:K_t],Omega_gamma_time[1:K_t,1:K_t,s])
      kappa_gamma[(K_t+1):(K_t+K_w),s] ~ dmnorm(zeros[1:K_w],Omega_gamma_week[1:K_w,1:K_w,s])
      kappa_gamma[(K_t+K_w+1):(K_t+K_w+K_tw),s] ~ dmnorm(zeros[1:K_tw],Omega_gamma_inter[1:K_tw,1:K_tw,s])
      
      gamma_time[1:(N),s] <-  gamma0[s] + X_tw[1:(N),1:(K_t)]%*%kappa_gamma[1:(K_t),s]
      gamma_week[1:(N),s] <- X_tw[1:(N),(K_t+1):(K_t+K_w)]%*%kappa_gamma[(K_t+1):(K_t+K_w),s]
      gamma_inter[1:(N),s] <-   X_tw[1:(N),(K_t+K_w+1):(K_t+K_w+K_tw)]%*%kappa_gamma[(K_t+K_w+1):(K_t+K_w+K_tw),s]
      
      gamma[1:N,s]<-matrix(gamma_time[1:(N),s]+gamma_week[1:(N),s]+gamma_inter[1:(N),s])
      
      # Delay splines
      for(d in 1:(D)){
        Omega_beta[1:K_t,1:K_t,s,d] <- S_t[1:K_t,1:K_t]*tau_beta[s,d]
        kappa_beta[1:K_t,s,d] ~ dmnorm(zeros[1:K_t],Omega_beta[1:K_t,1:K_t,s,d])
        beta[1:N,s,d] <- rep(psi[s,d],N) + X_t[1:N,1:K_t]%*%kappa_beta[1:K_t,s,d]
      }
    }
    
    # Smoothing parameter priors.
    for(s in 1:S){
      tau_delta[s] ~ dinvgamma(0.5,0.5)
      tau_gamma_t[s] ~ dinvgamma(0.5,0.5)
      tau_gamma_w[s] ~ dinvgamma(0.5,0.5)
      for(j in 1:2){
        tau_gamma_tw[s,j] ~ dinvgamma(0.5,0.5)
      }
      gamma0[s] ~ dnorm(0,sd=10) # Spatial intercept for time-week.
      for(d in 1:(D-1)){
        tau_beta[s,d]~ dinvgamma(0.5,0.5)
        psi[s,d] ~ dnorm(logit(1/(D+2-d)),sd=2)
      }
      for(d in 1:D){
        phi[s,d] ~ dgamma(2,0.02) # Beta-Binomial dispersion parameters.
      }
      iota[s] ~ dnorm(0,sd=10) # Spatial intercept.
      theta[s] ~ dgamma(2,0.02) # Negative-Binomial dispersion parameters.
    }
  })
  blank_data=data_frame(y=rnorm(N,0,1),t=1:N,w=(1:N)%%7)
  blank_jagam_c=jagam(y~ti(t,bs='cs',k=n_knots[1])+ti(w,bs='cc',k=n_knots[2])+
                        ti(t,w, bs=c('cs','cc'), k=n_knots[c(1,2)]),data=blank_data,file='blank.jags',
                      knots=list(t=seq(1,N,length=n_knots[1]),w=seq(0,7,length=n_knots[2])))
  
  
  
  # Censor data which wouldn't be available.
  censored_z <- z
  for(s in 1:S){
    censored_z[,s,][outer(1:dim(z[,s,])[1], 0:(dim(z[,s,])[2]-1), FUN = "+") > C_list[seed]] <- NA
  }
  
  # Observed index for each region at each delay.
  obs_index<-matrix(NA, nrow=S, ncol=D)
  for(s in 1:S){
    for(d in 1:D){
      obs_index[s,d]<-which(is.na(censored_z[,s,d])==TRUE)[1]-1
    }
  }
  
  censored_y <- apply(censored_z[,,1:D_max],c(1,2),sum)
  
  # Constants (e.g. number of weeks to model, number of knots) for NIMBLE.
  covid_rolling_nb_constants <-list(N=N,C=C_list[seed],D=D,S=S,obs_index=obs_index,
                                    K_t=n_knots[1]-1,
                                    K_w=dim(blank_jagam_c$jags.data$S2)[1],
                                    K_tw=dim(blank_jagam_c$jags.data$S3)[1])
  # Data (e.g. total cases, partial counts, spline model matrix) for NIMBLE.
  covid_rolling_nb_data <- list(z=abind(censored_z[1:C_list[seed],,1:D],
                                        apply(censored_z[1:C_list[seed],,1:D_max],c(1,2),sum)-
                                          apply(censored_z[1:C_list[seed],,1:D],c(1,2),sum),along=3),
                                X_t=blank_jagam_c$jags.data$X[,2:(covid_rolling_nb_constants$K_t+1)],
                                X_tw=blank_jagam_c$jags.data$X[,2:(covid_rolling_nb_constants$K_t+covid_rolling_nb_constants$K_w+covid_rolling_nb_constants$K_tw+1)],
                                S_t=blank_jagam_c$jags.data$S1,
                                S_w=blank_jagam_c$jags.data$S2,
                                S_tw=blank_jagam_c$jags.data$S3,
                                zeros=rep(0,max(covid_rolling_nb_constants$K_tw,covid_rolling_nb_constants$K_t,covid_rolling_nb_constants$K_w)))
  
  # Generate random initial values.
  covid_rolling_nb_inits <- list(kappa_delta=matrix(rnorm(S*covid_rolling_nb_constants$K_t,0,0.1),ncol=S),
                                 kappa_gamma=matrix(rnorm(S*(covid_rolling_nb_constants$K_t+covid_rolling_nb_constants$K_w+covid_rolling_nb_constants$K_tw),0,0.1),ncol=S),
                                 kappa_beta=array(rnorm(S*D*(covid_rolling_nb_constants$K_t),0,0.01),dim=c(covid_rolling_nb_constants$K_t,S,D)),
                                 iota=rnorm(S,0,1),
                                 theta=rexp(S,0.01),
                                 gamma0=rnorm(S,-10,1),
                                 psi=matrix(rnorm(S*D,0,1), nrow=S),
                                 tau_delta=rinvgamma(S,0.5,0.5),
                                 tau_gamma_t=rinvgamma(S,5,0.5),
                                 tau_gamma_w=rinvgamma(S,5,0.5),
                                 tau_gamma_tw=matrix(rinvgamma(2*S,5,0.5), nrow=S),
                                 tau_beta=matrix(rinvgamma(D*S,5,0.5), nrow=S),
                                 phi=matrix(rexp(D*S,0.01),nrow=S,ncol=D),
                                 z=array(NA,dim=c(C_list[seed],S,D+1)))
  for(s in 1:S){
    for(t in 1:C_list[seed]){
      for(d in 1:(D+1)){
        if(is.na(covid_rolling_nb_data$z[t,s,d])) covid_rolling_nb_inits$z[t,s,d]=rpois(1,median(covid_rolling_nb_data$z[,s,d],na.rm=T))
      }
    }
  }
  # Build the model.
  covid_rolling_nb_model <- nimbleModel(covid_nb_code,covid_rolling_nb_constants,covid_rolling_nb_data,covid_rolling_nb_inits)
  # Compile the model.
  covid_rolling_nb_compiled_model <- compileNimble(covid_rolling_nb_model)
  # Set up the MCMC.
  covid_rolling_nb_mcmc_config <- configureMCMC(covid_rolling_nb_model,monitors=c('iota','gamma0','mu_lambda', 'z',
                                                                                  'psi', 'theta','phi','lambda',
                                                                                  'kappa_gamma','kappa_beta',
                                                                                  'delta'))
  for(s in 1:S){
    for(d in 1:(D)){
      covid_rolling_nb_mcmc_config$removeSamplers(paste('kappa_beta[1:',(covid_rolling_nb_constants$K_t),', ',s,', ',d,']', sep=''))
      covid_rolling_nb_mcmc_config$addSampler(target=paste('kappa_beta[1:',(covid_rolling_nb_constants$K_t),', ',s,', ',d,']', sep=''),type='ess') 
    }
    # covid_rolling_nb_mcmc_config$removeSamplers(paste('kappa_gamma[1:',(covid_rolling_nb_constants$K_t),', ',s,']', sep=''))
    # covid_rolling_nb_mcmc_config$addSampler(target=paste('kappa_gamma[1:',(covid_rolling_nb_constants$K_t),', ',s,']', sep=''),type='ess') 
    covid_rolling_nb_mcmc_config$removeSamplers(paste('kappa_gamma[',(covid_rolling_nb_constants$K_t+1),':',(covid_rolling_nb_constants$K_t+covid_rolling_nb_constants$K_w),', ',s,']', sep=''))
    covid_rolling_nb_mcmc_config$addSampler(target=paste('kappa_gamma[',(covid_rolling_nb_constants$K_t+1),':',(covid_rolling_nb_constants$K_t+covid_rolling_nb_constants$K_w),', ',s,']', sep=''),type='ess') 
    covid_rolling_nb_mcmc_config$removeSamplers(paste('kappa_gamma[',(covid_rolling_nb_constants$K_t+covid_rolling_nb_constants$K_w+1),':',(covid_rolling_nb_constants$K_t+covid_rolling_nb_constants$K_w+covid_rolling_nb_constants$K_tw),', ',s,']', sep=''))
    covid_rolling_nb_mcmc_config$addSampler(target=paste('kappa_gamma[',(covid_rolling_nb_constants$K_t+covid_rolling_nb_constants$K_w+1),':',(covid_rolling_nb_constants$K_t+covid_rolling_nb_constants$K_w+covid_rolling_nb_constants$K_tw),', ',s,']', sep=''),type='ess') 
    
  }
  covid_rolling_nb_mcmc <- buildMCMC(covid_rolling_nb_mcmc_config)
  # Compile the MCMC.
  covid_rolling_nb_compiled_mcmc <- compileNimble(covid_rolling_nb_mcmc,project=covid_rolling_nb_model)
  
  covid_samples_cluster <- runMCMC(covid_rolling_nb_compiled_mcmc,niter=150000,nburnin=50000,inits=covid_rolling_nb_inits,nchains=1,samplesAsCodaMCMC = TRUE,thin=10)
  
  return(covid_samples_cluster)
  
}

nb_time<-system.time({
  covid_rolling_nb_samples  <- parLapply(cl = this_cluster, X = 1:length(C_list), 
                                          fun = run_MCMC_allcode, 
                                          C_list=C_list, n_knots=n_knots, N=N, S=S, 
                                          D_max=D_max, D=D, z=z)
})
stopCluster(this_cluster)

save(nb_time,file="nb_times_rollling.RData")
save(covid_rolling_nb_samples,file='covid_rolling_nb_samples.RData')

n_sim_censored <- dim(covid_rolling_nb_samples[[1]])[1]
# Simulate nowcasting and forecasting samples.
covid_rolling_nb_y <- lapply(covid_rolling_nb_samples,function(x){
  # Negative-Binomial dispersion parameters.
  theta <- select(as.data.frame(x),starts_with('theta'))%>%as.matrix()
  # Negative-Binomial means.
  mu_lambda <- select(as.data.frame(x),starts_with('mu_lambda'))%>%as.matrix()%>%array(dim=c(n_sim_censored,N,S,D+1))
  
  # Now-casting and forecasting samples of total covid cases.
  z <- array(dim=c(n_sim_censored,N,S,D+1))
  z_samples <- select(as.data.frame(x),starts_with('z'))%>%as.matrix()
  C <- dim(z_samples)[2]/(S*(D+1))
  z[,1:C,,] <- array(z_samples,dim=c(n_sim_censored,C,S,D+1))
  for(t in (C+1):N){
    for(d in 1:(D+1)){
      z[,t,,d] <- rnbinom(n_sim_censored*S,mu=mu_lambda[,t,,d],size=theta)
    }
  }
  y <- apply(z,c(1,2,3),sum)
  y <- abind(y,apply(y,c(1,2),sum),along=3)
  return(y)
})
covid_rolling_nb_y <- do.call('abind',list(covid_rolling_nb_y,along=4))

# Put into long format and join with observed values.
covid_rolling_nb_y <- apply(covid_rolling_nb_y,c(2,3,4),quantile,c(0.025,0.5,0.975),na.rm=T)%>%melt(varnames=c('quantile','t','s','c'),value.name='y')%>%
  spread(quantile,y)%>%mutate(r=c(region_names,'England')[s],c=dates[C_list[c]],delta=dates[t]-c)

covid_rolling_nb_y <- mutate(covid_rolling_nb_y,delta=dates[t]-c)
covid_rolling_nb_y <- left_join(covid_rolling_nb_y,y_data)%>%mutate(r=factor(r,levels=c(region_names,'England')))

#####################
#    Part Four:    #
#   Other Models    #
#  and Comparison   #
#####################

# Run the INLA model.
source('INLA_Model.R') # Takes a long time (hours).
# Run the NobBS models.
source('NobBS_Models.R')

# Put nowcasting quantiles from competing models into the same format.
inla_censored_y <- INLA_nowcasts%>%mutate(delta=t-now_day,median=Median,c=dates[now_day],model='INLA')%>%select(-Mode,-Median,-now_day)%>%
  left_join(y_data)%>%mutate(r=factor(r,levels=c(region_names)))
nobbs_censored_y <- NobBS_nowcasts%>%mutate(delta=t-now_day,t=as.numeric(t),median=Median,c=dates[now_day],model='NobBS')%>%
  select(-Mode,-Median,-now_day,-now_date,-n.reported,-onset_date)%>%
  left_join(y_data)%>%mutate(r=factor(r,levels=c(region_names)))
nobbs_14_censored_y <- NobBS_nowcasts_window%>%mutate(delta=t-now_day,t=as.numeric(t),median=Median,c=dates[now_day],model='NobBS-14')%>%
  select(-Mode,-Median,-now_day,-now_date,-n.reported,-onset_date)%>%
  left_join(y_data)%>%mutate(r=factor(r,levels=c(region_names)))

# Join the nowcasts together for comparison.
nowcasts_compare <- rbind(mutate(filter(covid_rolling_gdm_y,delta==0,c!="2020-04-14",s!=8),model='GDM',
                                 lower=`2.5%`,median=`50%`,upper=`97.5%`,`2.5%`=NULL,`50%`=NULL,`97.5%`=NULL),
                          mutate(filter(covid_rolling_nb_y,delta==0,s!=8),model='NB',
                                 lower=`2.5%`,median=`50%`,upper=`97.5%`,`2.5%`=NULL,`50%`=NULL,`97.5%`=NULL),
                          filter(nobbs_censored_y,delta==0),filter(inla_censored_y,delta==0),
                          filter(nobbs_14_censored_y,delta==0))


# Figure 8: Nowcast scatter plot.
scatter_compare <- ggplot(filter(nowcasts_compare,s!=8))+
  geom_abline(intercept = 0,slope=1)+ geom_point(aes(x=median,y=y,colour=model))+
  scale_colour_brewer(name=NULL,palette='Dark2')+
  facet_wrap(~factor(model,levels=c('GDM','NB','INLA','NobBS','NobBS-14')),nrow=1)+
  theme_minimal()+
  guides(colour=FALSE)+
  labs(x='Predicted Deaths',y='Observed Deaths',title='Same-Day Nowcasts of COVID-19 Deaths')+
  theme(strip.text = element_text(size=12),legend.position = 'bottom',
        plot.subtitle=element_text(size=14),
        plot.caption = element_text(size=10),plot.title=element_text(size=16),
        axis.text = element_text(size=10),axis.title = element_text(size=12),
        legend.text = element_text(size=10))
ggsave(scatter_compare,filename = 'Plots/scatter_compare.pdf',width=9,height=2.5)

# Compute nowcast metrics writing out to a spreadsheet.
preds_compare <- rbind(mutate(filter(covid_rolling_gdm_y,c!="2020-04-14",s!=8),model='GDM',
                              lower=`2.5%`,median=`50%`,upper=`97.5%`,`2.5%`=NULL,`50%`=NULL,`97.5%`=NULL),
                       mutate(filter(covid_rolling_nb_y,s!=8),model='NB',
                              lower=`2.5%`,median=`50%`,upper=`97.5%`,`2.5%`=NULL,`50%`=NULL,`97.5%`=NULL),
                       filter(nobbs_censored_y),filter(inla_censored_y),
                       filter(nobbs_14_censored_y))

# Compute metrics.
compare_models <- nowcasts_compare%>%group_by(model,r)%>%
  summarise(coverage=round(mean(y>=lower&y<=upper),2),rmse=round(sqrt(mean((y-median)^2))),piw=round(mean(upper-lower)),
            bias=round(mean(median-y)))

compare_models%>%group_by(model)%>%summarise_all(mean)%>%select(-r)

compare_models_out <- (compare_models%>%
                         myspread(model,c(coverage,bias,rmse,piw)))[,c('r','GDM_rmse','NB_rmse','INLA_rmse','NobBS_rmse','NobBS-14_rmse',
                                                                       'GDM_bias','NB_bias','INLA_bias','NobBS_bias','NobBS-14_bias',
                                                                       'GDM_piw','NB_piw','INLA_piw','NobBS_piw','NobBS-14_piw',
                                                                       'GDM_coverage','NB_coverage','INLA_coverage','NobBS_coverage','NobBS-14_coverage')]
compare_models_out <- rbind(compare_models_out,nowcasts_compare%>%group_by(model)%>%
                              summarise(coverage=round(mean(y>=lower&y<=upper),2),bias=round(mean(median-y)),rmse=round(sqrt(mean((y-median)^2))),piw=round(mean(upper-lower)))%>%
                              mutate(r='Overall')%>%myspread(model,c(coverage,rmse,piw,bias))
)
# Write to a spreadsheet.
write.xlsx(compare_models_out,file='compare_models_rounded.xlsx')

# Figure 9: Prediction interval width plots.
piw_plot <- ggplot()+geom_point(data=nowcasts_compare,aes(x=dates[t],y=upper-lower,colour=model,shape=model))+
  facet_wrap(~r,scales='free_y',nrow=2)+
  labs(x=NULL,y='95% Prediction Interval Width',title='Same-Day Nowcasts of COVID-19 Deaths')+
  scale_shape_discrete(name=NULL,guide=guide_legend(ncol=1))+
  scale_colour_brewer(name=NULL,palette='Dark2',guide=guide_legend(ncol=1))+
  theme_minimal()+
  theme(strip.text = element_text(size=12),legend.position = c(1-1/8,0.25),
        plot.subtitle=element_text(size=14),
        plot.caption = element_text(size=10),plot.title=element_text(size=16),
        axis.text = element_text(size=10),axis.title = element_text(size=12),
        legend.text = element_text(size=10))
ggsave(piw_plot,filename = 'Plots/piw.pdf',width=9,height=4,scale=1)


