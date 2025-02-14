
## SARI GDM with no missing values, no positive values, fitted using a cluster

library(tidyverse)
library(nimble)
library(reshape2)
library(mgcv)
library(doParallel)
library(coda)
library(ggfan)
library(gridExtra)
library(viridis)
library(scales)
nimbleOptions(oldConjugacyChecking = FALSE)
nimbleOptions(useNewConfigureMCMC = TRUE)

# Set your working directory.

source('Functions.R') # Load in the script containing various functions.
source('GDM_Survivor.R') # Load in the script containing the GDM survivor model.
load('SARI_Data.RData') # Load in the SARI data.

# About the objects in sari_data.RData:
# sari is a three-dimensional array: dimension 1 is time; dimension 2
# is region; dimension 3 is delay.
# population contains the populations for each health region from the 2010 census.

S=22 # Number of regions to model.
D=2 # Number of delays to explicit
n_chains <- 4 # Number of MCMC chains to run. Reduce if you don't have more than 8 CPU cores!
n_knots <- c(18,9) # Number of knots of the temporal and seasonal splines, respectively.
N <- 230 # Number of weeks to consider.
C <- 224 # Present-day week.

dates <- seq(as.Date("2013/1/7"),by=7,length=243) # Dates are approximate and note from here onwards we assume exactly 52 weeks in the year.

sari_total <- apply(sari,c(1,2),sum) # Calculate the total reported cases after 27 weeks.
sari_total_data <- melt(sari_total,varnames=c('t','r'),value.name='y')%>%mutate(r=as.factor(r),w=t%%52,pop=population[r],rate=100000*y/pop)

sari_state <- apply(sari,c(1,3),sum)
sari_state_total <- apply(sari_state,1,sum)


# Plot the total cases by region.
data_plot <- ggplot(sari_total_data)+geom_area(aes(x=dates[t],fill=as.factor(r),y=rate),position='identity',alpha=0.5)+
  labs(x=NULL,y=NULL,title='Recorded Severe Acute Respiratory Infection (SARI) Cases',subtitle='Weekly Total per 100,000 People, State of Paraná, Brazil')+guides(fill=FALSE)+
  scale_fill_viridis_d(end=0.8)+theme_minimal()+coord_cartesian(ylim=c(0,6))+geom_vline(xintercept=dates[C])
ggsave(data_plot,filename = 'Plots/sari_data_plot.pdf',width=9,height=2.5)

# Register the Beta-Binomial as a distribution for NIMBLE (see Functions.R for dbetabin).
registerDistributions(list(dbetabin=list(
  BUGSdist='dbetabin(mu,phi,size)',discrete=TRUE)))

week_index=rep(1:52,10) # Vector for index which week in the year each week is.

reduced_sari=sari[1:N,1:S,] # Reduce the data to the time period we want to consider.

reduced_y=rowSums(reduced_sari,dims=2) # Total cases of this time period.

# Censor the sari data according to the chosen present-day week (C).
censored_sari <- reduced_sari
for(s in 1:S){
  censored_sari[,s,][outer(1:dim(reduced_sari[,s,])[1], 0:(dim(reduced_sari[,s,])[2]-1), FUN = "+") > C] <- NA
}

censored_sari_state <- apply(censored_sari,c(1,3),sum)

# Example delayed reporting bar plot.
n_time=5
n_delay=5
tick_labels=numeric(n_time)
for(i in 1:(n_time-1)){
  tick_labels[i]=paste('t',as.character(-n_time+i),sep='')
}
tick_labels[n_time]='t'
bar_data=data.frame(time=rep(1:n_time,n_delay),delay=as.factor(sort(rep(1:n_delay,n_time))),
                    count=as.numeric(censored_sari_state[(C-n_time+1):C,1:n_delay]))
bar_plot=ggplot()+geom_col(data=data.frame(time=1:n_time,total=sari_state_total[(C-n_time+1):C]),aes(x=time,y=total))+
  geom_col(data=bar_data,aes(x=time,y=count,fill=delay),position = position_stack(reverse = TRUE))+
  labs(x='Week',y='Reported Cases',title='Delay Structure of SARI Data',subtitle='State of Paraná, Brazil') +
  theme_minimal()+scale_fill_viridis(name='Delay',discrete = TRUE,begin=0.2,end=0.8)+scale_x_reverse(breaks=1:n_time,labels=tick_labels)+coord_flip()
ggsave(bar_plot,file='Plots/sari_bar_plot.pdf',width=4.5,height=3)


# Set up the splines using jagam.
blank_data=data_frame(y=rnorm(N,0,1),t=1:N,w=week_index[1:N])
blank_jagam=jagam(y~s(t,bs='cs',k=n_knots[1])+s(w,bs='cc',k=n_knots[2]),data=blank_data,file='blank.jags',
                  knots=list(t=seq(1,C,length=n_knots[1]),w=seq(0,52,length=n_knots[2])))

obs_index<-matrix(NA, nrow=S, ncol=D)
for(s in 1:S){
  for(d in 1:D){
    if(sum(is.na(censored_sari[,s,d]))>0){obs_index[s,d]<-which(is.na(censored_sari[,s,d])==TRUE)[1]-1}else{
      obs_index[s,d]<-C  }
}}

# Constants (e.g. number of weeks to model, number of knots) for NIMBLE.
sari_constants <- list(N=N,C=C,S=S,K_t=n_knots[1]-1,K_w=n_knots[2]-2, D=D,obs_index=obs_index,w=week_index[1:N])
# Data (e.g. total cases, partial counts, spline model matrix) for NIMBLE.
sari_data <- list(z=censored_sari[1:C,,1:D],y=apply(censored_sari,c(1,2),sum)[1:C,],X_t=blank_jagam$jags.data$X[,2:(n_knots[1])],S_t=blank_jagam$jags.data$S1,
                  X_w=blank_jagam$jags.data$X[,(n_knots[1]+1):(n_knots[1]+n_knots[2]-2)],S_w=blank_jagam$jags.data$S2,
                  zeros=rep(0,max(sari_constants$K_t,sari_constants$K_w)),population=population)


library(parallel)
sari_nona_cluster<- makeCluster(4)

### cluster function 

run_sari_nona_cluster<- function( seed, sari_constants, sari_data ) {
  
  
  
  #load libraries
  library(nimble)
  library(tidyverse)
  library(trustOptim)
  library(numDeriv)
  
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
  
  #constants
  S<-sari_constants$S
  C<-sari_constants$C
  N<-sari_constants$N
  D<-sari_constants$D
  
  # NIMBLE model code.
  sari_code <- nimbleCode({
    for(s in 1:S){
      for(t in 1:N){
        # Mean total cases.
        log(lambda[t,s]) <-  log(population[s]) + iota[s] + delta[t,s] + xi[w[t],s]
        for(d in 1:D){
          # Expected cumulative proportions.
          probit(p[t,s,d]) <- beta[s,d]+gamma[t,s]
        }
        # Relative proportions (Beta-Binomial means).
        nu[t,s,1] <- p[t,s,1]
        for(d in 2:D){
          nu[t,s,d] <- (p[t,s,d]-p[t,s,d-1])/(1-p[t,s,d-1])
        }
      }
      for(t in 1:C){
        # Model for total cases.
        y[t,s] ~ dnegbin(theta[s]/(theta[s]+lambda[t,s]),theta[s])
        # Model for partial counts.
        z[t,s,1] ~ dbetabin(nu[t,s,1],phi[s,1],y[t,s])
      }
      for(d in 2:D){
        for(t in 1:obs_index[s,d]){
          z[t,s,d] ~ dbetabin(nu[t,s,d],phi[s,d],y[t,s]-sum(z[t,s,1:(d-1)]))
        }
      }
    }
    
    ## Overall effects
    # Temporal spline
    Omega_alpha[1:K_t,1:K_t] <- S_t[1:K_t,1:K_t]*exp(log_tau_alpha)
    kappa_alpha[1:K_t] ~ dmnorm(zeros[1:K_t],Omega_alpha[1:K_t,1:K_t])
    alpha[1:N] <- X_t[1:N,1:K_t]%*%kappa_alpha[1:K_t]
    # Seasonal spline
    Omega_eta[1:K_w,1:K_w] <- S_w[1:K_w,1:K_w]*exp(log_tau_eta)
    kappa_eta[1:K_w] ~ dmnorm(zeros[1:K_w],Omega_eta[1:K_w,1:K_w])
    eta[1:52] <- X_w[1:52,1:K_w]%*%kappa_eta[1:K_w]
    # Delay spline
    Omega_psi[1:K_t,1:K_t] <- S_t[1:K_t,1:K_t]*exp(log_tau_psi)
    kappa_psi[1:K_t] ~ dmnorm(zeros[1:K_t],Omega_psi[1:K_t,1:K_t])
    psi[1:N] <- X_t[1:N,1:K_t]%*%kappa_psi[1:K_t]
    ## Regional effects
    for(s in 1:S){
      Omega_delta[1:K_t,1:K_t,s] <- S_t[1:K_t,1:K_t]*exp(log_tau_delta[s])
      kappa_delta[1:K_t,s] ~ dmnorm(kappa_alpha[1:K_t],Omega_delta[1:K_t,1:K_t,s])
      delta[1:N,s] <- X_t[1:N,1:K_t]%*%kappa_delta[1:K_t,s]
      # Seasonal spline
      Omega_xi[1:K_w,1:K_w,s] <- S_w[1:K_w,1:K_w]*exp(log_tau_xi[s])
      kappa_xi[1:K_w,s] ~ dmnorm(kappa_eta[1:K_w],Omega_xi[1:K_w,1:K_w,s])
      xi[1:52,s] <- X_w[1:52,1:K_w]%*%kappa_xi[1:K_w,s]
      # Delay spline
      Omega_gamma[1:K_t,1:K_t,s] <- S_t[1:K_t,1:K_t]*exp(log_tau_gamma[s])
      kappa_gamma[1:K_t,s] ~ dmnorm(kappa_psi[1:K_t],Omega_gamma[1:K_t,1:K_t,s])
      gamma[1:N,s] <- X_t[1:N,1:K_t]%*%kappa_gamma[1:K_t,s]
    }
    # Smoothing parameter priors.
    # tau_alpha ~ dinvgamma(0.5,0.5) # Equivalent to Half-Normal(0,1) on 1/sqrt(tau).
    #  tau_eta ~ dinvgamma(0.5,0.5)
    # tau_psi ~ dinvgamma(0.5,0.5)
    log_tau_alpha ~ dnorm(1.26,sd=2.23)
    log_tau_psi ~ dnorm(1.26,sd=2.23)
    log_tau_eta ~ dnorm(1.26,sd=2.23)
    
    for(s in 1:S){
      #  tau_delta[s] ~ dinvgamma(0.5,0.5)
      #  tau_xi[s] ~ dinvgamma(0.5,0.5)
      # tau_gamma[s] ~ dinvgamma(0.5,0.5)
      log_tau_delta[s] ~ dnorm(1.26,sd=2.23)
      log_tau_gamma[s] ~ dnorm(1.26,sd=2.23)
      log_tau_xi[s] ~ dnorm(1.26,sd=2.23)
      for(d in 1:D){
        # phi[s,d] ~ dgamma(2,0.02) # Beta-Binomial dispersion parameters.
        log_phi[s,d] ~ dnorm(4.3,sd=0.8)
        phi[s,d]<-exp(log_phi[s,d])
      }
      
      omega_beta[s,1] ~ dnorm(0,sd=10)
      beta[s,1]<-omega_beta[s,1]
      for(d in 2:D){
        omega_beta[s,d] ~ dnorm(0,sd=2)
        beta[s,d]<-beta[s,d-1]+exp(omega_beta[s,d])  # Independent delay effects
      }
      
      iota[s] ~ dnorm(0,sd=10) # Spatial intercept.
      #theta[s] ~ dgamma(2,0.02) # Negative-Binomial dispersion parameters.
      log_theta[s] ~ dnorm(4.3,sd=0.8)
      theta[s]<-exp(log_theta[s])
    }
  })
  
    # Generate random initial values.
    sari_inits <- list(kappa_alpha=rnorm(sari_constants$K_t,0,0.1),
                            kappa_psi=rnorm(sari_constants$K_t,0,0.1),kappa_eta=rnorm(sari_constants$K_w,0,0.1),
                            kappa_delta=matrix(rnorm(S*sari_constants$K_t,0,0.1),ncol=S),
                            kappa_gamma=matrix(rnorm(S*sari_constants$K_t,0,0.1),ncol=S),
                            kappa_xi=matrix(rnorm(S*sari_constants$K_w,0,0.1),ncol=S),
                            iota=rnorm(S,0,1),
                            log_theta=rnorm(S,4,1.3),
                            log_tau_gamma=rnorm(S,1.26,sd=2.23),
                            log_tau_delta=rnorm(S,1.26,sd=2.23),
                            log_tau_xi=rnorm(S,1.26,sd=2.23),
                            log_phi=matrix(rnorm(S*D,4,1.3),nrow=S,ncol=D),
                            log_tau_psi=rnorm(1,1.26,sd=2.23),
                            log_tau_alpha=rnorm(1,1.26,sd=2.23),
                            log_tau_eta=rnorm(1,1.26,sd=2.23),
                            omega_beta=matrix(rnorm(D*S,0,0.1),nrow=S,ncol=D),
                            y=matrix(NA,nrow=C,ncol=S))
    for(s in 1:S){
      for(t in 1:C){
        if(is.na(sari_data$y[t,s])) sari_inits$y[t,s]=sum(c(sari_inits$z[t,s,],sari_data$z[t,s,],rpois(1,median(sari_data$y[,s]-rowSums(sari_data$z[,s,]),na.rm=T))),na.rm=TRUE)
      }
    }
    # Build the model.
    sari_model <- nimbleModel(sari_code,sari_constants,sari_data,sari_inits)
    # Compile the model.
    sari_compiled_model<- compileNimble(sari_model)
    # Set up the MCMC.
    sari_mcmc_config <- configureMCMC(sari_model,monitors=c('alpha','eta','iota','beta','psi',
                                                                      'theta','phi','gamma','delta','xi','y','lambda'))
    sari_mcmc <- buildMCMC(sari_mcmc_config)
    # Compile the MCMC.
    sari_compiled_mcmc <- compileNimble(sari_mcmc,project=sari_model)
    sari_samples <- runMCMC(sari_compiled_mcmc,niter=1200000,nburnin=200000,inits=sari_inits,nchains=1,samplesAsCodaMCMC = TRUE,thin=1000)
    
    return(sari_samples)
}


time_nona_cluster<-system.time({
  sari_samples <- parLapply(cl = sari_nona_cluster, X = 1:4, 
                            fun = run_sari_nona_cluster, 
                            sari_constants = sari_constants, sari_data=sari_data)
})

stopCluster(sari_nona_cluster)


#PLOTS

sari_samples<-as.mcmc.list(sari_samples)
# Combine all MCMC chains.
sari_combined_samples <- as_tibble(do.call('rbind',sari_samples))

time_nona_cluster[3]/60/60

## Check convergence of lambda, unobserved y, and theta.
lambda_index <- which(dimnames(sari_combined_samples)[[2]]=='lambda[1, 1]'):which(dimnames(sari_combined_samples)[[2]]==paste('lambda[',N,', ',S,']',sep=''))
lambda_psrf <- sapply(lambda_index,function(x)gelman.diag(sari_samples[,x],autoburnin = FALSE,transform = TRUE,multivariate = FALSE)$psrf[,1])

y_index <- which(dimnames(sari_combined_samples)[[2]]=='y[1, 1]'):which(dimnames(sari_combined_samples)[[2]]==paste('y[',C,', ',S,']',sep=''))
y_index_matrix <- matrix(y_index,nrow=C,ncol=S)
y_unob_index <- as.numeric(y_index_matrix[(C-27+2):C,])

y_unob_psrf <- sapply(y_unob_index,function(x)gelman.diag(sari_samples[,x],autoburnin = FALSE,transform = TRUE,multivariate = FALSE)$psrf[,1])

theta_index <- which(dimnames(sari_combined_samples)[[2]]=='theta[1]'):which(dimnames(sari_combined_samples)[[2]]==paste('theta[',S,']',sep=''))
theta_psrf <- sapply(theta_index,function(x)gelman.diag(sari_samples[,x],autoburnin = FALSE,transform = TRUE,multivariate = FALSE)$psrf[,1])

# Proportion of PSRFs below 1.05.
mean(lambda_psrf<1.05)
mean(na.omit(y_unob_psrf)<1.05)
mean(theta_psrf<1.05)

# Proportion of PSRFs below 1.2.
mean(lambda_psrf<1.20)
mean(na.omit(y_unob_psrf<1.20))

#ESS
ESS_SARI_all_opt<-effectiveSize(sari_samples)
#effective sample size
print('summary of effective sample size of posterior samples ')
summary(ESS_SARI_all_opt)
#ESS lambda
print('summary of ESS of lambda samples')
ESS_SARI_lambda_opt<-effectiveSize(sari_samples[,lambda_index])
summary(ESS_SARI_lambda_opt)
#ESS unobserved y
print('summary of ESS of unobserved y samples')
ESS_SARI_y_unob_opt<-effectiveSize(sari_samples[,y_unob_index])
summary(ESS_SARI_y_unob_opt)

# Number of MCMC samples.
n_sim <- dim(sari_combined_samples)[1]

# Overall temporal effect of SARI incidence.
sari_alpha <- select(sari_combined_samples,contains('alpha'))%>%as.matrix()%>%apply(2,median)%>%melt(value.name='y')%>%mutate(x=1:N)
alpha_plot <- ggplot(sari_alpha)+geom_line(aes(x=x,y=y))+labs(x='Week',y=NULL,title='Temporal Effect',subtitle='on Dengue Incidence')

# Overall seasonal effect on SARI incidence.
sari_eta <- select(sari_combined_samples,starts_with('eta'))%>%as.matrix()%>%apply(2,median)%>%melt(value.name='y')%>%mutate(x=1:52)
eta_plot <- ggplot(sari_eta)+geom_line(aes(x=x,y=y))+labs(x='Week',y=NULL,title='Seasonal Effect',subtitle='on Dengue Incidence')

# Overall temporal effect on cumulative proportion reported.
sari_psi <- select(sari_combined_samples,starts_with('psi'))%>%as.matrix()%>%apply(2,median)%>%melt(value.name='y')%>%mutate(x=1:N)
psi_plot <- ggplot(sari_psi)+geom_line(aes(x=x,y=y))+labs(x='Week',y=NULL,title='Temporal Effect',subtitle='on Cumulative Proportion Reported')

# Regional temporal effect of SARI incidence.
sari_delta <- select(sari_combined_samples,starts_with('delta'))%>%as.matrix()%>%array(dim=c(n_sim,N,S))%>%
  apply(c(2,3),quantile,c(0.025,0.5,0.975))%>%melt(varnames=c('quantile','t','s'),value.name='y')%>%
  spread(quantile,y)

# Regional seasonal effect of SARI incidence.
sari_xi <- select(sari_combined_samples,starts_with('xi'))%>%as.matrix()%>%array(dim=c(n_sim,52,S))%>%
  apply(c(2,3),quantile,c(0.025,0.5,0.975))%>%melt(varnames=c('quantile','t','s'),value.name='y')%>%
  spread(quantile,y)

# Regional temporal effect on cumulative proportion reported.
sari_gamma <- select(sari_combined_samples,starts_with('gamma'))%>%as.matrix()%>%array(dim=c(n_sim,N,S))%>%
  apply(c(2,3),quantile,c(0.025,0.5,0.975))%>%melt(varnames=c('quantile','t','s'),value.name='y')%>%
  spread(quantile,y)


# Plot the temporal and seasonal effects.
delta_plot <- ggplot(sari_delta)+
  geom_line(aes(colour=as.factor(s),x=dates[t],y=`50%`))+geom_line(data=sari_alpha,aes(x=dates[x],y=y),linetype=2)+theme_minimal()+
  labs(x=NULL,y=NULL,title='Temporal Effect',subtitle='on SARI Incidence')+guides(colour=FALSE)+scale_colour_viridis_d(end=0.8)
xi_plot <- ggplot(sari_xi)+
  geom_line(aes(colour=as.factor(s),x=t,y=`50%`))+geom_line(data=sari_eta,aes(x=x,y=y),linetype=2)+scale_colour_viridis_d(end=0.8)+
  labs(x=NULL,y=NULL,title='Seasonal Effect',subtitle='on SARI Incidence')+guides(colour=FALSE)+theme_minimal()+
  scale_x_continuous(breaks=c(52/24+0.5+(0:5)*52/6),labels=c('Jan','Mar','May','Jul','Sep','Nov'))
gamma_plot <- ggplot(sari_gamma)+
  geom_line(aes(colour=as.factor(s),x=dates[t],y=`50%`))+geom_line(data=sari_psi,aes(x=dates[x],y=y),linetype=2)+scale_colour_viridis_d(end=0.8)+
  labs(x=NULL,y=NULL,title='Temporal Effect',subtitle='on Cumulative Proportion Reported')+guides(colour=FALSE)+theme_minimal()

effect_plots <- arrangeGrob(delta_plot,xi_plot,gamma_plot,nrow=1)
#ggsave(effect_plots,file='Plots/sari_effect_plots.pdf',width=9,height=2.5)  


# Negative-Binomial dispersion parameters.
sari_theta <- select(sari_combined_samples,starts_with('theta'))%>%as.matrix()
theta_optim<-apply(sari_theta,2,median)
# Beta-Binomial dispersion 
sari_phi <- select(sari_combined_samples,starts_with('phi'))%>%as.matrix()%>%array(dim=c(n_sim,S,D))
phi_optim<-apply(sari_phi,c(2,3),median)

# Negative-Binomial means.
sari_lambda <- select(sari_combined_samples,starts_with('lambda'))%>%as.matrix()%>%array(dim=c(n_sim,N,S))

# Now-casting and forecasting samples of total SARI cases.
sari_y <- array(dim=c(n_sim,N,S))
sari_y[,1:C,] <- select(sari_combined_samples,starts_with('y'))%>%as.matrix()%>%array(dim=c(n_sim,C,S))
for(s in 1:S){
  sari_y[,(C+1):N,s] <- rnbinom(n_sim*(N-C),mu=sari_lambda[,(C+1):N,s],size=sari_theta[,s])
}

# Vector containing health region names.
region_names <- rep(NA,22)
region_names[c(2,17,15)] <- c('Curitiba','Londrina','Maringá')

# Total cases for each region.
sari_y <- sari_y%>%apply(c(2,3),quantile,c(0.025,0.1,0.175,0.25,0.5,0.75,0.825,0.9,0.975))%>%melt(varnames=c('quantile','t','s'),value.name='y')%>%
  spread(quantile,y)%>%mutate(r=region_names[s])

# Total cases data into long format.
y_data <- reduced_y%>%melt(varnames=c('t','s'),value.name='y')%>%mutate(r=region_names[s])
Y_data <- apply(reduced_y,1,sum)%>%melt(varnames=c('t'),value.name='y')%>%mutate(t=1:N)

# Predictions for the three most populous regions.
forecast_plot <- ggplot(filter(sari_y,s%in%c(2,17,15)))+
  geom_ribbon(aes(x=dates[t],ymin=`2.5%`,ymax=`97.5%`,fill=as.factor(s)),alpha=0.25)+
  geom_ribbon(aes(x=dates[t],ymin=`10%`,ymax=`90%`,fill=as.factor(s)),alpha=0.25)+
  geom_ribbon(aes(x=dates[t],ymin=`17.5%`,ymax=`82.5%`,fill=as.factor(s)),alpha=0.25)+
  geom_ribbon(aes(x=dates[t],ymin=`25%`,ymax=`75%`,fill=as.factor(s)),alpha=0.25)+
  geom_line(aes(x=dates[t],y=`50%`,colour=as.factor(s)))+
  geom_point(data=filter(y_data,s%in%c(2,17,15)),aes(x=dates[t],y=y,colour=as.factor(s)))+
  facet_wrap(~r,nrow=1,scales='free')+
  labs(x=NULL,y=NULL,title='Predicted Weekly Total SARI Cases')+guides(colour=FALSE,fill=FALSE)+
  coord_cartesian(xlim=c(dates[200],dates[230]))+geom_vline(xintercept=dates[C])+theme_minimal()+
  scale_fill_viridis_d(begin=0.2,end=0.8)+scale_colour_viridis_d(begin=0.2,end=0.8)+
  theme(strip.text = element_text(size=12))+scale_x_date(labels = date_format("%b"),date_breaks = '2 months')
ggsave(forecast_plot,file='Plots/sari_nona_cluster_forecast.pdf',width=9,height=2.5)
forecast_plot




#TRACE PLOTS
#trace plots
# lambda
plot(sari_samples[,'lambda[1, 1]'], main="optim lambda[1,1]")
plot(sari_samples[,'lambda[36, 7]'], main="optim lambda[36,7]")
plot(sari_samples[,'lambda[21, 5]'], main="optim lambda[21,5]")
#beta
plot(sari_samples[,'beta[1, 1]'], main="optim beta[1,1]")
plot(sari_samples[,'beta[2, 2]'], main="optim beta[2,4]")
#delta
plot(sari_samples[,'delta[1, 2]'], main="optim delta[1, 2]")
plot(sari_samples[,'delta[2, 3]'], main="optim delta[2, 3]")
#gamma
plot(sari_samples[,'gamma[15, 4]'], main="optim gamma[15, 4]")
plot(sari_samples[,'gamma[30, 20]'], main="optim gamma[30, 20]")
#xi
plot(sari_samples[,'xi[13, 4]'], main="optim xi[13, 4]")
plot(sari_samples[,'xi[50, 18]'], main="optim xi[50, 18]")
#theta
plot(sari_samples[,'theta[4]'], main="optim theta[4]")
plot(sari_samples[,'theta[5]'], main="optim theta[5]")
#eta 
plot(sari_samples[,'eta[2]'], main="optim eta[2]")
plot(sari_samples[,'eta[4]'], main="optim eta[4]")
plot(sari_samples[,'eta[9]'], main="optim eta[9]")
plot(sari_samples[,'eta[12]'], main="optim eta[12]")
#alpha
plot(sari_samples[,'alpha[3]'], main="optim alpha[3]")
plot(sari_samples[,'alpha[5]'], main="optim alpha[5]")
plot(sari_samples[,'alpha[13]'], main="optim alpha[13]")
plot(sari_samples[,'alpha[20]'], main="optim alpha[20]")
#psi
plot(sari_samples[,'psi[1]'], main="optim psi[1]")
plot(sari_samples[,'psi[11]'], main="optim psi[11]")
plot(sari_samples[,'psi[19]'], main="optim psi[19]")
plot(sari_samples[,'psi[22]'], main="optim psi[22]")
