### original GDM for SARI data

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
    }
}

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
  Omega_alpha[1:K_t,1:K_t] <- S_t[1:K_t,1:K_t]*tau_alpha
  kappa_alpha[1:K_t] ~ dmnorm(zeros[1:K_t],Omega_alpha[1:K_t,1:K_t])
  alpha[1:N] <- X_t[1:N,1:K_t]%*%kappa_alpha[1:K_t]
  # Seasonal spline
  Omega_eta[1:K_w,1:K_w] <- S_w[1:K_w,1:K_w]*tau_eta
  kappa_eta[1:K_w] ~ dmnorm(zeros[1:K_w],Omega_eta[1:K_w,1:K_w])
  eta[1:52] <- X_w[1:52,1:K_w]%*%kappa_eta[1:K_w]
  # Delay spline
  Omega_psi[1:K_t,1:K_t] <- S_t[1:K_t,1:K_t]*tau_psi
  kappa_psi[1:K_t] ~ dmnorm(zeros[1:K_t],Omega_psi[1:K_t,1:K_t])
  psi[1:N] <- X_t[1:N,1:K_t]%*%kappa_psi[1:K_t]
  ## Regional effects
  for(s in 1:S){
    Omega_delta[1:K_t,1:K_t,s] <- S_t[1:K_t,1:K_t]*tau_delta[s]
    kappa_delta[1:K_t,s] ~ dmnorm(kappa_alpha[1:K_t],Omega_delta[1:K_t,1:K_t,s])
    delta[1:N,s] <- X_t[1:N,1:K_t]%*%kappa_delta[1:K_t,s]
    # Seasonal spline
    Omega_xi[1:K_w,1:K_w,s] <- S_w[1:K_w,1:K_w]*tau_xi[s]
    kappa_xi[1:K_w,s] ~ dmnorm(kappa_eta[1:K_w],Omega_xi[1:K_w,1:K_w,s])
    xi[1:52,s] <- X_w[1:52,1:K_w]%*%kappa_xi[1:K_w,s]
    # Delay spline
    Omega_gamma[1:K_t,1:K_t,s] <- S_t[1:K_t,1:K_t]*tau_gamma[s]
    kappa_gamma[1:K_t,s] ~ dmnorm(kappa_psi[1:K_t],Omega_gamma[1:K_t,1:K_t,s])
    gamma[1:N,s] <- X_t[1:N,1:K_t]%*%kappa_gamma[1:K_t,s]
  }
  # Smoothing parameter priors.
  tau_alpha ~ dinvgamma(0.5,0.5) # Equivalent to Half-Normal(0,1) on 1/sqrt(tau).
  tau_eta ~ dinvgamma(0.5,0.5)
  tau_psi ~ dinvgamma(0.5,0.5)
  for(s in 1:S){
    tau_delta[s] ~ dinvgamma(0.5,0.5)
    tau_xi[s] ~ dinvgamma(0.5,0.5)
    tau_gamma[s] ~ dinvgamma(0.5,0.5)
    for(d in 1:D){
      phi[s,d] ~ dgamma(2,0.02) # Beta-Binomial dispersion parameters.
    }
    beta[s,1] ~ dnorm(0,sd=10) # Independent delay effects.
    for(d in 2:D){
      beta[s,d] ~ T(dnorm(beta[s,d-1],sd=10),beta[s,d-1],)
    }
    iota[s] ~ dnorm(0,sd=10) # Spatial intercept.
    theta[s] ~ dgamma(2,0.02) # Negative-Binomial dispersion parameters.
    
  }
})

# Constants (e.g. number of weeks to model, number of knots) for NIMBLE.
sari_constants <- list(N=N,C=C,S=S,K_t=n_knots[1]-1,K_w=n_knots[2]-2,w=week_index[1:N], obs_index=obs_index)
# Data (e.g. total cases, partial counts, spline model matrix) for NIMBLE.
sari_data <- list(z=censored_sari[1:C,,1:D],y=apply(censored_sari,c(1,2),sum)[1:C,],X_t=blank_jagam$jags.data$X[,2:(n_knots[1])],S_t=blank_jagam$jags.data$S1,
                  X_w=blank_jagam$jags.data$X[,(n_knots[1]+1):(n_knots[1]+n_knots[2]-2)],S_w=blank_jagam$jags.data$S2,
                  zeros=rep(0,max(sari_constants$K_t,sari_constants$K_w)),population=population)

# Set up the NIMBLE model for each chain.
sari_inits <- sari_model <- sari_compiled_model <- sari_mcmc_config <- sari_mcmc <- sari_compiled_mcmc <- list()
for(i in 1:n_chains){
  # Generate random initial values.
  sari_inits[[i]] <- list(kappa_alpha=rnorm(sari_constants$K_t,0,0.1),kappa_eta=rnorm(sari_constants$K_w,0,0.1),
                          kappa_psi=rnorm(sari_constants$K_t,0,0.1),kappa_delta=matrix(rnorm(S*sari_constants$K_t,0,0.1),ncol=S),
                          kappa_xi=matrix(rnorm(S*sari_constants$K_w,0,0.1),ncol=S),kappa_gamma=matrix(rnorm(S*sari_constants$K_t,0,0.1),ncol=S),
                          iota=rnorm(S,0,1),theta=rexp(S,0.01),tau_gamma=rinvgamma(S,0.5,0.5),
                          tau_delta=rinvgamma(S,0.5,0.5),tau_xi= rinvgamma(S,0.5,0.5),
                          phi=matrix(rexp(D*S,0.01),nrow=S,ncol=D),tau_psi=rinvgamma(1,0.5,0.5),
                          tau_alpha=rinvgamma(1,0.5,0.5),tau_eta= rinvgamma(1,0.5,0.5),
                          beta=t(apply(matrix(rnorm(D*S,0,0.1),nrow=S,ncol=D),1,sort)),y=matrix(NA,nrow=C,ncol=S)  )
  for(s in 1:S){
    for(t in 1:C){
      if(is.na(sari_data$y[t,s])) sari_inits[[i]]$y[t,s]=sum(c(sari_inits[[i]]$z[t,s,],sari_data$z[t,s,],rpois(1,median(sari_data$y[,s]-rowSums(sari_data$z[,s,]),na.rm=T))),na.rm=TRUE)
    }
  }
}
SARI_time<-system.time({
  for(i in 1:n_chains){
    # Build the model.
    sari_model[[i]] <- nimbleModel(sari_code,sari_constants,sari_data,sari_inits[[i]])
    # Compile the model.
    sari_compiled_model[[i]] <- compileNimble(sari_model[[i]])
    # Set up the MCMC.
    sari_mcmc_config[[i]] <- configureMCMC(sari_model[[i]],monitors=c('alpha','eta','iota','beta','psi',
                                                                      'theta','phi','gamma','delta','xi','y','lambda'))
    sari_mcmc[[i]] <- buildMCMC(sari_mcmc_config[[i]])
    # Compile the MCMC.
    sari_compiled_mcmc[[i]] <- compileNimble(sari_mcmc[[i]],project=sari_model[[i]])
  }
  
  # Run the model in parallel. Change 'dopar' to 'do' to run the model on only one core.
  registerDoParallel(cores=n_chains)
  sari_samples <- as.mcmc.list(foreach(i=1:n_chains)%dopar%{
    runMCMC(sari_compiled_mcmc[[i]],niter=2600000,nburnin=1000000,inits=sari_inits[[i]],nchains=1,samplesAsCodaMCMC = TRUE,thin=1000)
  })
  
})

# Combine all MCMC chains.
sari_combined_samples <- as_tibble(do.call('rbind',sari_samples))

SARI_time

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
mean(na.omit(y_unob_psrf)<1.20)

#ESS
ESS_SARI_all<-effectiveSize(sari_samples)
#effective sample size
print('summary of effective sample size of posterior samples ')
summary(ESS_SARI_all)
#ESS lambda
print('summary of ESS of lambda samples')
ESS_SARI_lambda<-effectiveSize(sari_samples[,lambda_index])
summary(ESS_SARI_lambda)
#ESS unobserved y
print('summary of ESS of unobserved y samples')
ESS_SARI_y_unob<-effectiveSize(sari_samples[,y_unob_index])
summary(ESS_SARI_y_unob)

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
effect_plots

# Negative-Binomial dispersion parameters.
sari_theta <- select(sari_combined_samples,starts_with('theta'))%>%as.matrix()
theta_og<-apply(sari_theta,2,median)

sari_phi <- select(sari_combined_samples,starts_with('phi'))%>%as.matrix()%>%array(dim=c(n_sim,S,D))
phi_og<-apply(sari_phi,c(2,3),median)

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
forecast_plot

ggsave(forecast_plot,file='Plots/sari_forecast_nona_pos.pdf',width=9,height=2.5)

# Calculate back-casting and forecasting coverage.
backcast_coverage <- forecast_coverage <- numeric(S)
for(i in 1:S){
  backcast_coverage[i] <- mean(sari_total[(C-27+2):(C-1),i]>=filter(sari_y,s==i,t%in%((C-27+2):(C-1)))$`2.5%`&sari_total[(C-27+2):(C-1),i]<=filter(sari_y,s==i,t%in%((C-27+2):(C-1)))$`97.5%`)
  forecast_coverage[i] <- mean(sari_total[(C+1):N,i]>=filter(sari_y,s==i,t%in%((C+1):N))$`2.5%`&sari_total[(C+1):N,i]<=filter(sari_y,s==i,t%in%((C+1):N))$`97.5%`)
}
mean(backcast_coverage)
mean(forecast_coverage)

# Calculate now-casting coverage.
mean(sari_total[C,]>=filter(sari_y,t==C)$`2.5%`&sari_total[C,]<=filter(sari_y,t==C)$`97.5%`)

D_max <- 27 # Maximum delay before y all cases assumed reported.

# Simulate posterior replicates of the total cases for each region.
y_replicates <- array(dim=c(n_sim,C-D_max+1,S))
for(i in 1:S){
  y_replicates[,,i] <- rnbinom(n_sim*(C-D_max+1),mu=sari_lambda[,1:(C-D_max+1),i],size=sari_theta[,i])
}

# Replicates of the total for the whole state.
Y_replicates <- apply(y_replicates,c(1,2),sum)



#save(SARI_time,ESS_SARI_all,ESS_SARI_lambda, ESS_SARI_y_unob, file = 'SARI_original.RData')
save(sari_samples,file='sari_samples_nona_pos.RData')



#trace plots
# lambda
plot(sari_samples[,'lambda[1, 1]'], main="long lambda[1,1]")
plot(sari_samples[,'lambda[36, 7]'], main="long lambda[36,7]")
plot(sari_samples[,'lambda[21, 5]'], main="long lambda[21,5]")
#beta
plot(sari_samples[,'beta[1, 1]'], main="long beta[1,1]")
plot(sari_samples[,'beta[2, 2]'], main="long beta[2,4]")
#theta
plot(sari_samples[,'theta[4]'], main="long theta[4]")
plot(sari_samples[,'theta[5]'], main="long theta[5]")
#delta
plot(sari_samples[,'delta[1, 2]'], main="long delta[1, 2]")
plot(sari_samples[,'delta[2, 3]'], main="long delta[2, 3]")
#alpha
plot(sari_samples[,'alpha[200]'], main="long alpha[200]")
plot(sari_samples[,'alpha[112]'], main="long alpha[112]")
#gamma
plot(sari_samples[,'gamma[15, 4]'], main="long gamma[15, 4]")
plot(sari_samples[,'gamma[150, 4]'], main="long gamma[150, 4]")
#psi
plot(sari_samples[,'psi[145]'], main="long psi[145]")
plot(sari_samples[,'psi[220]'], main="long psi[220]")
#xi
plot(sari_samples[,'xi[13, 4]'], main="long xi[13, 4]")
plot(sari_samples[,'xi[30, 6]'], main="long xi[30, 6]")
#eta
plot(sari_samples[,'eta[32]'], main="long eta[32]")
plot(sari_samples[,'eta[46]'], main="long eta[46]")


