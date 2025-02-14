############################################
#  Simulation Experiment: DATA GENERATION  #
############################################

library(ggplot2)
library(tidyverse)
library(reshape2)
library(mgcv)
library(coda)
library(abind)
library(nimble)
library(doParallel)
library(lubridate)
library(formattable)
library(compositions)
library(gridExtra)
library(MASS)
library(dplyr)

# Read in some necessary functions.
set.seed(345)

# Beta-Binomial (rbetabin) Function
rbetabin=nimbleFunction(run=function(n=integer(0),mu=double(0),phi=double(0),size=double(0)){
  pi=rbeta(n,mu*phi,(1-mu)*phi)
  returnType(double(0))
  return(rbinom(n,size,pi))
})

# Number of weeks in the time series.
N <- 100
# Number of delays.
D <- 7
# Number of regions.
S <- 27
# Linear temporal trend in total cases.

##### FIT GDM TO DATA #####
# Create function for cluster to carry out MCMC model:
Simulation_data<- function(seed, S, N, D){
  # Load libraries.
  library(nimble)
  library(tidyverse)
  library(mgcv)
  
  # Define the beta-Binomial as a distribution for NIMBLE.
  dbetabin=nimbleFunction(run=function(x=double(0),mu=double(0),phi=double(0),size=double(0),log=integer(0)){
    returnType(double(0))
    phi <- min(phi,1e+04) # Hard upper limit on phi for computational stability.
    if(x>=0&x<=size){
      return(lgamma(size+1)+lgamma(x+mu*phi)+lgamma(size-x+(1-mu)*phi)+lgamma(phi)-
               lgamma(size+phi)-lgamma(mu*phi)-lgamma((1-mu)*phi)-lgamma(size-x+1)-lgamma(x+1))
    }else{
      return(-Inf)
    }
  })
  
  rbetabin=nimbleFunction(run=function(n=integer(0),mu=double(0),phi=double(0),size=double(0)){
    phi <- min(phi,1e+04) # Hard upper limit on phi for computational stability.
    pi=rbeta(1,mu*phi,(1-mu)*phi)
    returnType(double(0))
    return(rbinom(1,size,pi))
  })
  
  
  assign('dbetabin', dbetabin, envir = .GlobalEnv)
  assign('rbetabin', rbetabin, envir = .GlobalEnv)
  
  registerDistributions(list(dbetabin=list(
    BUGSdist='dbetabin(mu,phi,size)',discrete=TRUE)))
  
  # set seed 
  set.seed(seed)
  # Generate spline bases functions for modelling and simulation.
  blank_data=tibble(y=rnorm(N,0,1),t=1:N)
  blank_jagam=jagam(y~s(t,k=floor(N/5), bs='tp'),
                    data=blank_data,file='blank.jags',
                    knots=list(t=seq(1,N,length=floor(N/5))))
  
  
  # Simulation data:
  # Simulate random zero mean total cases temporal tend from thin plate spline.
  K_t=dim(blank_jagam$jags.data$S1)[1]
  S_t=blank_jagam$jags.data$S1
  X_t=blank_jagam$jags.data$X[,2:(K_t+1)]
  
  # Set spline precision parameter.
  sigma_alpha_sim <- c(5,5)
  # Generate thin plate splines 
  Omega_alpha_sim <- S_t[1:(K_t),1:(K_t)]/sigma_alpha_sim[1]^2 + S_t[1:(K_t),(K_t+1):(2*K_t)]/sigma_alpha_sim[2]^2
  kappa_alpha_sim <- rmnorm_chol(n=1,mean=rep(0,K_t), chol(Omega_alpha_sim[1:(K_t),1:(K_t)]))
  sim_alpha_raw<- X_t[1:N,1:(K_t)]%*%kappa_alpha_sim[1:(K_t)]
  # Calculate linear trend in simulated spline
  alpha_linear_trend<-round(as.numeric(lm(data=tibble(y=sim_alpha_raw,x=(1:N-mean(1:N))))$coefficients[2]),3)
  # Take linear trend out of spline
  sim_alpha<-sim_alpha_raw - (1:N-mean(1:N))*alpha_linear_trend
  
  zeta_value<-sort(rep(c(-0.1,0,0.1),9))
  
  # Spatial intercept of mean reported cases - iota 
  sim_iota<- abs(rnorm(S,5,sd=0.25))
  # Overall temporal trend in total cases for all zeta values
  # Simulate mean number of total cases - lambda
  sim_lambda <- array(dim=c(N,S))
  sign<-(rbernoulli(1,0.5)-1)
  sign[sign==0]=1
  #sign
  for(s in 1:S){
    for(t in 1:N){
      sim_lambda[t,s]<-exp(sim_iota[s]+sim_alpha[t]*sign+zeta_value[s]*((t-mean(1:N))/sd(1:N)))
    }
  }
  # Simulate Negative Binomial dispersion parameter.
  sim_theta<-rgamma(S,2,0.02) 
  # Negative Binomial distribution for total counts (y).
  sim_y<-array(dim=c(N,S))
  for(s in 1:S){
    sim_y[,s]<-rnbinom(N,mu=sim_lambda[,s],size=sim_theta[s])
  }
  
  # GDM simulation model for partial reports (z).
  # Independent delay effects linear improvement over time (psi). 
  sim_psi<-matrix(rep((rnorm(D,0.25,0.5)),S), nrow=S, ncol=D, byrow = TRUE)
  for(s in 1:S){
    sim_psi[s,]<-sort(sim_psi[s,])
  }
  
  # Linear temporal trend in reporting performance over time (eta)
  sim_eta<-matrix(nrow=N,ncol=S)
  eta_slopes<-rep(rep(seq(from=0.2, to=-0.2, length=3),3),3)
  for(s in 1:S){
    for(t in 1:N){
      sim_eta[t,s]<-eta_slopes[s]*((t-mean(1:N))/sd(1:N))
    }
  }
  
  # Linear relationship between mean cases (lambda) and delay length.
  delta_slopes <- rep(c(rep(-0.2,3),rep(0,3),rep(0.2,3)),3)
  sim_delta<-array(dim=c(N,S))
  for(s in 1:S){
    for(t in 1:N){
      sim_delta[t,s]<-delta_slopes[s]*(sim_y[t,s]-mean(sim_y[,s]))/sd(sim_y[,s])
    }
  }
  sim_p<-array(dim=c(N,S,D))
  for(d in 1:D){
    for(s in 1:S){
      sim_p[,s,d] <-  iprobit(sim_eta[,s]+sim_psi[s,d]+sim_delta[,s])
    }}
  
  # Relative proportions (eta-Binomial means - nu).
  sim_nu<-array(dim=c(N,S,D))
  for(s in 1:S){
    sim_nu[,s,1] <- sim_p[,s,1]
    for(d in 2:D){
      sim_nu[,s,d] <- (sim_p[,s,d]-sim_p[,s,d-1])/(1-sim_p[,s,d-1])
    }}
  
  # eta-Binomial dispersion parameters (phi).
  # sim_phi<-matrix(rep(rgamma(D,2,0.02),S), nrow=S, ncol=D, byrow=TRUE) 
  sim_phi<-matrix((rgamma(D*S,2,0.02)), nrow=S, ncol=D, byrow=TRUE) 
  sim_z<-array(dim=c(N,S,D))
  for(t in 1:N){
    for(s in 1:S){
      sim_z[t,s,1] <- rbetabin(1,sim_nu[t,s,1],sim_phi[s,1],sim_y[t,s])
      for(d in 2:(D)){
        sim_z[t,s,d] <- rbetabin(1,sim_nu[t,s,d],sim_phi[s,d],sim_y[t,s]-sum(sim_z[t,s,1:(d-1)]))
      }
    }
  }
  
  simulation_data <- list(y=sim_y, z=sim_z[,1:S,1:D], p=sim_p, lambda=sim_lambda, zeta=zeta_value,eta=eta_slopes,delta=delta_slopes)
  
  
  return(simulation_data)
}


# Number of simulations.
sims=100
# Set maximum COVID-19 lab report delay.
n_cores<-min(detectCores(),sims)

# Make Cluster for MCMC. 
this_cluster_simulation_data<- makeCluster(n_cores)

# Run Cluster.
time_data_gen <- system.time({
  simulation_data_list <-  parLapply(cl = this_cluster_simulation_data, X = 1:sims, 
                                     fun = Simulation_data,
                                     S=S,
                                     N=N,
                                     D=D)
})

# Stop Cluster.
stopCluster(this_cluster_simulation_data)

save(simulation_data_list,time_data_gen,file="Simulation_data_y.RData")

# Simulated data from simulation_data.R script:
sim_y<-list()
sim_lambda<-list()
sim_p<-list()
for(i in 1:sims){
  sim_y[[i]]<-matrix(simulation_data_list[[i]]$y,nrow=N,ncol=S)%>%melt(varnames=c('t','s'),value.name='y')%>%
    mutate(sim=i)
  sim_lambda[[i]]<-matrix(simulation_data_list[[i]]$lambda,nrow=N,ncol=S)%>%melt(varnames=c('t','s'),value.name='lambda')%>%
    mutate(sim=i)
  sim_p[[i]]<-array(simulation_data_list[[i]]$p,dim=dim(simulation_data_list[[i]]$p))%>%melt(varnames=c('t','s','d'),value.name='p')%>%
    mutate(sim=i)
}
# Table of simulated total counts (y)
sim_y_long<-unlist(sim_y)%>%array(dim=c(dim(sim_y[[1]]),sims))%>%apply(2,c)
colnames(sim_y_long)<- colnames(sim_y[[1]])
sim_y_long<-as_tibble(sim_y_long)%>%
  mutate(eta=simulation_data_list[[1]]$eta[s],zeta=simulation_data_list[[1]]$zeta[s],delta=simulation_data_list[[1]]$delta[s])%>%
  mutate(s=as.factor(s),sim=as.factor(sim),zeta=factor(zeta,levels=unique(zeta)),
         eta=factor(eta,levels=unique(eta)),delta=factor(delta,levels=unique(delta)))
# Table of simulated mean total counts (lambda).
sim_lambda_long<-unlist(sim_lambda)%>%array(dim=c(dim(sim_lambda[[1]]),sims))%>%apply(2,c)
colnames(sim_lambda_long)<- colnames(sim_lambda[[1]])
sim_lambda_long<-as_tibble(sim_lambda_long)%>%
  mutate(eta=simulation_data_list[[1]]$eta[s],zeta=simulation_data_list[[1]]$zeta[s],delta=simulation_data_list[[1]]$delta[s])%>%
  mutate(s=as.factor(s),sim=as.factor(sim),zeta=factor(zeta,levels=unique(zeta)),
         eta=factor(eta,levels=unique(eta)),delta=factor(delta,levels=unique(delta)))
# Table of simulated cumulative proportions (p).
sim_p_long<-unlist(sim_p)%>%array(dim=c(dim(sim_p[[1]]),sims))%>%apply(2,c)
colnames(sim_p_long)<- colnames(sim_p[[1]])
sim_p_long<-as_tibble(sim_p_long)%>%mutate(eta=simulation_data_list[[1]]$eta[s],zeta=simulation_data_list[[1]]$zeta[s],delta=simulation_data_list[[1]]$delta[s])%>%
  mutate(s=as.factor(s),sim=as.factor(sim),zeta=factor(zeta,levels=unique(zeta)),
         eta=factor(eta,levels=unique(eta)),delta=factor(delta,levels=unique(delta)))

# Plot of cumulative proportions
p_sim_plot <- ggplot(filter(sim_p_long, sim%in%c(1),s%in%c(10:18)))+geom_line(aes(x=t,y=p, colour=as.factor(d-1)))+
  theme_minimal()+
  facet_grid(paste('Case load effect:',delta)~paste('Delay temporal trend:',eta))+
  labs(title = "Simulated cumulative proportions reported",
       subtitle="For no linear temporal trend in the total cases (zeta=0)" ,x="Time (t)",
       y=expression(S['t,d']))+
  scale_colour_viridis_d(name='Delay (d)',begin=0.1,end=0.8)+
  scale_linetype(name='Simulation')+
  theme_minimal()+
  theme(text = element_text(size=14), axis.title = element_text(size=16), plot.title = element_text(size=20),plot.subtitle=element_text(size=18),
        legend.title = element_text(size=16),strip.text = element_text(size=14),
        axis.text = element_text(size=9),legend.position="bottom")+guides(colour=guide_legend(nrow=1))
p_sim_plot
ggsave(p_sim_plot,file='Plots/p_sim_plot_y.pdf',width=9,height=8)       

# plot simulated link between case load and delay 
data_sim_long<-full_join(full_join(sim_p_long,sim_y_long, by=c('t','s','eta','zeta','delta','sim')),sim_lambda_long, by=c('t','s','eta','zeta','delta','sim'))


# Plot of mean reported cases:
lambda_sim_plot <- ggplot(filter(sim_lambda_long, sim%in%c(1:4),s%in%c(1,10,19)))+geom_line(aes(x=t,y=lambda, colour=zeta, linetype=(sim)))+
  theme_minimal()+
  facet_wrap(~paste('Totals temporal trend:',zeta), scales='free')+
  labs(title = "Simulation examples of the mean total cases temporal trend",
       subtitle="For three different overall linear temporal trends in the total cases", x="Time (t)",y=expression(lambda['t,s']))+
  scale_colour_viridis_d(name='Linear temporal trend in total cases (zeta)',begin=0.1,end=0.8)+
  scale_linetype(name='Simulation')+
  theme_minimal()+
  theme(text = element_text(size=14), axis.title = element_text(size=16), plot.title = element_text(size=20),plot.subtitle=element_text(size=18),
        legend.title = element_text(size=16),strip.text = element_text(size=14),
        axis.text = element_text(size=9),legend.position="bottom")
lambda_sim_plot
ggsave(lambda_sim_plot,file='Plots/lambda_sim_plot_y.pdf',width=10,height=4)       

p_y_sim_plot <- ggplot(filter(data_sim_long, sim%in%c(1)))+
  geom_point(aes(x=log(y),y=probit(p), colour=as.factor(d)))+
  geom_smooth(aes(x=log(y),y=probit(p), colour=as.factor(d)))+
  theme_minimal()+
  facet_grid(paste('Case load effect:',delta)~paste('Delay temporal trend:',eta))+
  labs(title = "Simulated cumulative proportions reported",
       subtitle="For three values of case load effect and three values of temporal trend in the delay" ,x=expression(y['t']),
       y=expression(S['t,d']))+
  scale_colour_viridis_d(name='Delay',begin=0.1,end=0.8)+
  scale_linetype(name='Simulation')+
  theme_minimal()+
  theme(text = element_text(size=12), plot.title = element_text(size=16),plot.subtitle=element_text(size=14),
        legend.title = element_text(size=12),strip.text = element_text(size=14), legend.position = "bottom")
p_y_sim_plot

p_lambda_sim_plot <- ggplot(filter(data_sim_long, sim%in%c(3),d==1))+
  geom_point(aes(x=log(lambda),y=probit(p), colour=as.factor(d)))+
  geom_smooth(aes(x=log(lambda),y=probit(p), colour=as.factor(d)))+
  theme_minimal()+
  facet_grid(paste('Case load effect:',delta)~paste('Delay temporal trend:',eta))+
  labs(title = "Simulated cumulative proportions reported",
       subtitle="For three values of case load effect and three values of temporal trend in the delay" ,x=expression(lambda['t']),
       y=expression(S['t,d']))+
  scale_colour_viridis_d(name='Delay',begin=0.1,end=0.8)+
  scale_linetype(name='Simulation')+
  theme_minimal()+
  theme(text = element_text(size=12), plot.title = element_text(size=16),plot.subtitle=element_text(size=14),
        legend.title = element_text(size=12),strip.text = element_text(size=14), legend.position = "bottom")
p_lambda_sim_plot
