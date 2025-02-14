
library(tidyverse)
library(nimble)
library(reshape2)
library(tidyr)
library(dplyr)
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
#source('GDM_Survivor.R') # Load in the script containing the GDM survivor model.
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
ggsave(data_plot,filename = 'Sari opt plots/sari_data_plot.pdf',width=9,height=2.5)

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


# Total cases data into long format.
y_data <- reduced_y%>%melt(varnames=c('t','s'),value.name='y')%>%mutate(r=region_names[s])
Y_data <- apply(reduced_y,1,sum)%>%melt(varnames=c('t'),value.name='y')%>%mutate(t=1:N)

#load MCMC samples from GDM model
#load("~/Dropbox/GDM optimisation/sari_samples.RData")

# # Combine all MCMC chains.
# sari_combined_samples <- as_tibble(do.call('rbind',sari_samples))
# 
# 
# ## Check convergence of lambda, unobserved y, and theta.
# lambda_index <- which(dimnames(sari_combined_samples)[[2]]=='lambda[1, 1]'):which(dimnames(sari_combined_samples)[[2]]==paste('lambda[',N,', ',S,']',sep=''))
# lambda_psrf <- sapply(lambda_index,function(x)gelman.diag(sari_samples[,x],autoburnin = FALSE,transform = TRUE,multivariate = FALSE)$psrf[,1])
# 
# y_index <- which(dimnames(sari_combined_samples)[[2]]=='y[1, 1]'):which(dimnames(sari_combined_samples)[[2]]==paste('y[',C,', ',S,']',sep=''))
# y_index_matrix <- matrix(y_index,nrow=C,ncol=S)
# y_unob_index <- as.numeric(y_index_matrix[(C-27+2):C,])
# 
# y_unob_psrf <- sapply(y_unob_index,function(x)gelman.diag(sari_samples[,x],autoburnin = FALSE,transform = TRUE,multivariate = FALSE)$psrf[,1])
# 
# theta_index <- which(dimnames(sari_combined_samples)[[2]]=='theta[1]'):which(dimnames(sari_combined_samples)[[2]]==paste('theta[',S,']',sep=''))
# theta_psrf <- sapply(theta_index,function(x)gelman.diag(sari_samples[,x],autoburnin = FALSE,transform = TRUE,multivariate = FALSE)$psrf[,1])
# 
# # Proportion of PSRFs below 1.05.
# mean(lambda_psrf<1.05)
# mean(na.omit(y_unob_psrf)<1.05)
# mean(theta_psrf<1.05)
# 
# # Proportion of PSRFs below 1.2.
# mean(lambda_psrf<1.20)
# mean(na.omit(y_unob_psrf)<1.20)
# 
# #ESS
# ESS_SARI_all<-effectiveSize(sari_samples)
# #effective sample size
# print('summary of effective sample size of posterior samples ')
# summary(ESS_SARI_all)
# #ESS lambda
# print('summary of ESS of lambda samples')
# ESS_SARI_lambda<-effectiveSize(sari_samples[,lambda_index])
# summary(ESS_SARI_lambda)
# #ESS unobserved y
# print('summary of ESS of unobserved y samples')
# ESS_SARI_y_unob<-effectiveSize(sari_samples[,y_unob_index])
# summary(ESS_SARI_y_unob)

# # Number of MCMC samples.
# n_sim <- dim(sari_combined_samples)[1]
# 
# # Overall temporal effect of SARI incidence.
# sari_alpha <- select(sari_combined_samples,contains('alpha'))%>%as.matrix()%>%apply(2,median)%>%melt(value.name='y')%>%mutate(x=1:N)
# alpha_plot <- ggplot(sari_alpha)+geom_line(aes(x=x,y=y))+labs(x='Week',y=NULL,title='Temporal Effect',subtitle='on Dengue Incidence')
# 
# # Overall seasonal effect on SARI incidence.
# sari_eta <- select(sari_combined_samples,starts_with('eta'))%>%as.matrix()%>%apply(2,median)%>%melt(value.name='y')%>%mutate(x=1:52)
# eta_plot <- ggplot(sari_eta)+geom_line(aes(x=x,y=y))+labs(x='Week',y=NULL,title='Seasonal Effect',subtitle='on Dengue Incidence')
# 
# # Overall temporal effect on cumulative proportion reported.
# sari_psi <- select(sari_combined_samples,starts_with('psi'))%>%as.matrix()%>%apply(2,median)%>%melt(value.name='y')%>%mutate(x=1:N)
# psi_plot <- ggplot(sari_psi)+geom_line(aes(x=x,y=y))+labs(x='Week',y=NULL,title='Temporal Effect',subtitle='on Cumulative Proportion Reported')
# 
# # Regional temporal effect of SARI incidence.
# sari_delta <- select(sari_combined_samples,starts_with('delta'))%>%as.matrix()%>%array(dim=c(n_sim,N,S))%>%
#   apply(c(2,3),quantile,c(0.025,0.5,0.975))%>%melt(varnames=c('quantile','t','s'),value.name='y')%>%
#   spread(quantile,y)
# 
# # Regional seasonal effect of SARI incidence.
# sari_xi <- select(sari_combined_samples,starts_with('xi'))%>%as.matrix()%>%array(dim=c(n_sim,52,S))%>%
#   apply(c(2,3),quantile,c(0.025,0.5,0.975))%>%melt(varnames=c('quantile','t','s'),value.name='y')%>%
#   spread(quantile,y)
# 
# # Regional temporal effect on cumulative proportion reported.
# sari_gamma <- select(sari_combined_samples,starts_with('gamma'))%>%as.matrix()%>%array(dim=c(n_sim,N,S))%>%
#   apply(c(2,3),quantile,c(0.025,0.5,0.975))%>%melt(varnames=c('quantile','t','s'),value.name='y')%>%
#   spread(quantile,y)
# 
# 
# # Plot the temporal and seasonal effects.
# delta_plot <- ggplot(sari_delta)+
#   geom_line(aes(colour=as.factor(s),x=dates[t],y=`50%`))+geom_line(data=sari_alpha,aes(x=dates[x],y=y),linetype=2)+theme_minimal()+
#   labs(x=NULL,y=NULL,title='Temporal Effect',subtitle='on SARI Incidence')+guides(colour=FALSE)+scale_colour_viridis_d(end=0.8)
# xi_plot <- ggplot(sari_xi)+
#   geom_line(aes(colour=as.factor(s),x=t,y=`50%`))+geom_line(data=sari_eta,aes(x=x,y=y),linetype=2)+scale_colour_viridis_d(end=0.8)+
#   labs(x=NULL,y=NULL,title='Seasonal Effect',subtitle='on SARI Incidence')+guides(colour=FALSE)+theme_minimal()+
#   scale_x_continuous(breaks=c(52/24+0.5+(0:5)*52/6),labels=c('Jan','Mar','May','Jul','Sep','Nov'))
# gamma_plot <- ggplot(sari_gamma)+
#   geom_line(aes(colour=as.factor(s),x=dates[t],y=`50%`))+geom_line(data=sari_psi,aes(x=dates[x],y=y),linetype=2)+scale_colour_viridis_d(end=0.8)+
#   labs(x=NULL,y=NULL,title='Temporal Effect',subtitle='on Cumulative Proportion Reported')+guides(colour=FALSE)+theme_minimal()
# 
# effect_plots <- arrangeGrob(delta_plot,xi_plot,gamma_plot,nrow=1)
# #ggsave(effect_plots,file='Sari opt plots/sari_effect_plots.pdf',width=9,height=2.5)  
# effect_plots
# 
# # Negative-Binomial dispersion parameters.
# sari_theta <- select(sari_combined_samples,starts_with('theta'))%>%as.matrix()
# theta_og<-apply(sari_theta,2,median)
# 
# sari_phi <- select(sari_combined_samples,starts_with('phi'))%>%as.matrix()%>%array(dim=c(n_sim,S,D))
# phi_og<-apply(sari_phi,c(2,3),median)
# 
# # Negative-Binomial means.
# sari_lambda <- select(sari_combined_samples,starts_with('lambda'))%>%as.matrix()%>%array(dim=c(n_sim,N,S))
# 
# # Now-casting and forecasting samples of total SARI cases.
# sari_y <- array(dim=c(n_sim,N,S))
# sari_y[,1:C,] <- select(sari_combined_samples,starts_with('y'))%>%as.matrix()%>%array(dim=c(n_sim,C,S))
# for(s in 1:S){
#   sari_y[,(C+1):N,s] <- rnbinom(n_sim*(N-C),mu=sari_lambda[,(C+1):N,s],size=sari_theta[,s])
# }
# 
# # Vector containing health region names.
# region_names <- rep(NA,22)
# region_names[c(2,17,15)] <- c('Curitiba','Londrina','Maringá')
# 
# # Total cases for each region.
# sari_y <- sari_y%>%apply(c(2,3),quantile,c(0.025,0.1,0.175,0.25,0.5,0.75,0.825,0.9,0.975))%>%melt(varnames=c('quantile','t','s'),value.name='y')%>%
#   spread(quantile,y)%>%mutate(r=region_names[s])
# 
# # Total cases data into long format.
# y_data <- reduced_y%>%melt(varnames=c('t','s'),value.name='y')%>%mutate(r=region_names[s])
# Y_data <- apply(reduced_y,1,sum)%>%melt(varnames=c('t'),value.name='y')%>%mutate(t=1:N)
# 
# # Predictions for the three most populous regions.
# forecast_plot <- ggplot(filter(sari_y,s%in%c(2,17,15)))+
#   geom_ribbon(aes(x=dates[t],ymin=`2.5%`,ymax=`97.5%`,fill=as.factor(s)),alpha=0.25)+
#   geom_ribbon(aes(x=dates[t],ymin=`10%`,ymax=`90%`,fill=as.factor(s)),alpha=0.25)+
#   geom_ribbon(aes(x=dates[t],ymin=`17.5%`,ymax=`82.5%`,fill=as.factor(s)),alpha=0.25)+
#   geom_ribbon(aes(x=dates[t],ymin=`25%`,ymax=`75%`,fill=as.factor(s)),alpha=0.25)+
#   geom_line(aes(x=dates[t],y=`50%`,colour=as.factor(s)))+
#   geom_point(data=filter(y_data,s%in%c(2,17,15)),aes(x=dates[t],y=y,colour=as.factor(s)))+
#   facet_wrap(~r,nrow=1,scales='free')+
#   labs(x=NULL,y=NULL,title='Predicted Weekly Total SARI Cases')+guides(colour=FALSE,fill=FALSE)+
#   coord_cartesian(xlim=c(dates[200],dates[230]))+geom_vline(xintercept=dates[C])+theme_minimal()+
#   scale_fill_viridis_d(begin=0.2,end=0.8)+scale_colour_viridis_d(begin=0.2,end=0.8)+
#   theme(strip.text = element_text(size=12))+scale_x_date(labels = date_format("%b"),date_breaks = '2 months')
# ggsave(forecast_plot,file='Sari opt plots/sari_forecast1.pdf',width=9,height=2.5)
# forecast_plot

##### OPTIMISATION OF SARI MODEL

# Set up the splines using jagam.
blank_data=data_frame(y=rnorm(N,0,1),t=1:N,w=week_index[1:N])
blank_jagam=jagam(y~s(t,bs='cs',k=n_knots[1])+s(w,bs='cc',k=n_knots[2]),data=blank_data,file='blank.jags',
                  knots=list(t=seq(1,C,length=n_knots[1]),w=seq(0,52,length=n_knots[2])))
#Maximum Delay
D_max<-dim(censored_sari)[3]
#total counts
y<-apply(censored_sari,c(1,2),sum)[1:C,]
#Delay's with remainder
D_rem<-D+1
S_all<-S

f<-function(i, censored_sari,N,C,S,D, population){
  library(tidyverse)
  library(nimble)
  library(reshape2)
  library(mgcv)
  library(doParallel)
  library(coda)
  library(gridExtra)
  library(abind)
  library(viridis)
  library(scales)
  library(trustOptim)
  library(numDeriv)
  nimbleOptions(oldConjugacyChecking = FALSE)
  nimbleOptions(useNewConfigureMCMC = TRUE)
  #population
  #log_pop<-log(population[i])
  
  # Set up the splines using jagam.
  n_knots <- c(18,9)
  week_index=rep(1:52,10) # Vector for index which week in the year each week is.
  blank_data=data_frame(y=rnorm(N,0,1),t=1:N,w=week_index[1:N])
  blank_jagam=jagam(y~s(t,bs='cs',k=n_knots[1])+s(w,bs='cc',k=n_knots[2]),data=blank_data,file='blank.jags',
                    knots=list(t=seq(1,C,length=n_knots[1]),w=seq(0,52,length=n_knots[2])))
  #partial counts 
  D_max<-dim(censored_sari)[3]
  z<-array(censored_sari[,i,],dim=c(N,S,D_max))
  D_rem<-D+1
  
  #unobserved y
  y<-apply(censored_sari,c(1,2),sum)[1:C,]
  unobs_y<-which(is.na(y[,i])) 
  
  #remainder term for z's
  remainder=matrix((apply(z,2,rowSums, na.rm=TRUE)-apply(array(z[,,1:D],dim=c(N,S,D)),2,rowSums, na.rm=TRUE))[1:C,], ncol=S)
  remainder[unobs_y[1]:C,]<-NA # reminder unknown when y is unknown
  z=array(abind(array(z[1:C,,1:D],dim=c(C,S,D)), remainder)[1:C,,1:(D+1)], dim=c(C,S,D+1))
  
  obs_index<-matrix(NA, nrow=S, ncol=D_rem)
  for(s in 1:S){
    for(d in 1:D_rem){
      if(sum(is.na(z[,s,d]))>0){obs_index[s,d]<-which(is.na(z[,s,d])==TRUE)[1]-1}else{
        obs_index[s,d]<-C  }
    }
  }
  
  # Constants (e.g. number of weeks to model, number of knots) for NIMBLE.
  sari_constants_NBz <- list(N=N,C=C,S=S,K_t=n_knots[1]-1,K_w=n_knots[2]-2, D=D_rem,obs_index=obs_index,w=week_index[1:N])
  # Data (e.g. total cases, partial counts, spline model matrix) for NIMBLE.
  sari_data_NBz <- list(z=array(z[1:C,,1:D_rem], dim=c(C,S,D_rem)),X_t=blank_jagam$jags.data$X[,2:(n_knots[1])],S_t=blank_jagam$jags.data$S1,
                        X_w=blank_jagam$jags.data$X[,(n_knots[1]+1):(n_knots[1]+n_knots[2]-2)],S_w=blank_jagam$jags.data$S2,
                        zeros=rep(0,max(sari_constants_NBz$K_t,sari_constants_NBz$K_w)))
  
  # Optimisation model 
  sari_gamma_code <- nimbleCode({
    for(s in 1:S){
      for(t in 1:N){
        # Mean total cases.
        log(lambda[t,s]) <-  iota[s] + delta[t,s] + xi[w[t],s]
        for(d in 1:(D-1)){
          # Expected cumulative proportions.
          probit(p[t,s,d]) <- beta[s,d]+gamma[t,s]
        }
        p[t,s,D]<-1
        # Relative proportions (Negative-Binomial proportions).
        mu[t,s,1]<- p[t,s,1]
        for(d in 2:D){
          mu[t,s,d] <- p[t,s,d] - p[t,s,(d-1)]
        }
      }
    }
    
    # Model for partial reports.
    for(s in 1:S){
      for(d in 1:D){
        for(t in 1:obs_index[s,d]){
          prob[t,s,d]<-(phi[s,d])/(phi[s,d]+lambda[t,s]*mu[t,s,d])
          z[t,s,d] ~ dnegbin( prob[t,s,d], phi[s,d])  
        }
      }
    }
    
    ## Regional effects
    for(s in 1:S){
      Omega_delta[1:K_t,1:K_t,s] <- S_t[1:K_t,1:K_t]*exp(log_tau_delta[s])
      kappa_delta[1:K_t,s] ~ dmnorm(zeros[1:K_t],Omega_delta[1:K_t,1:K_t,s])
      delta[1:N,s] <- X_t[1:N,1:K_t]%*%kappa_delta[1:K_t,s]
      # Seasonal spline
      Omega_xi[1:K_w,1:K_w,s] <- S_w[1:K_w,1:K_w]*exp(log_tau_xi[s])
      kappa_xi[1:K_w,s] ~ dmnorm(zeros[1:K_w],Omega_xi[1:K_w,1:K_w,s])
      xi[1:52,s] <- X_w[1:52,1:K_w]%*%kappa_xi[1:K_w,s]
      # Delay spline
      Omega_gamma[1:K_t,1:K_t,s] <- S_t[1:K_t,1:K_t]*exp(log_tau_gamma[s])
      kappa_gamma[1:K_t,s] ~ dmnorm(zeros[1:K_t],Omega_gamma[1:K_t,1:K_t,s])
      gamma[1:N,s] <- X_t[1:N,1:K_t]%*%kappa_gamma[1:K_t,s]
    }
    
    for(s in 1:S){
      # Smoothing parameter priors.
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
      for(d in 2:(D-1)){
        omega_beta[s,d] ~ dnorm(0,sd=2)
        beta[s,d]<-beta[s,d-1]+exp(omega_beta[s,d])  # Independent delay effects
      }
      
      iota[s] ~ dnorm(0,sd=10) # Spatial intercept.
      #theta[s] ~ dgamma(2,0.02) # Negative-Binomial dispersion parameters.
      #log_theta[s] ~ dnorm(4.3,sd=0.8)
      #theta[s]<-exp(log_theta[s])
    }
  })
  
  
  #initial values for gamma parameter (mean of y's)
  z_inits<-array(NA, dim=c(N,S,D_rem))
  z_inits[1:C,,1:D_rem]<-z
  for(t in 1:N){
    for(s in 1:S){
      for(d in 1:D_rem){
        if(is.na(z_inits[t,s,d])){
          z_inits[t,s,d]<-rpois(1,median(z_inits[(t-4):(t-1),s,d],na.rm=T))
        }
      }
      
    }
  }
  gamma_inits<-apply(z_inits,2,rowSums)
  
  # Generate random initial values.
  sari_inits_NBz <- list( #kappa_psi=rnorm(sari_constants_NBz$K_t,0,0.1),
    #kappa_eta=rnorm(sari_constants_NBz$K_w,0,0.1),
    kappa_gamma=matrix(rnorm(S*sari_constants_NBz$K_t,0,0.1),ncol=S),
    kappa_xi=matrix(rnorm(S*sari_constants_NBz$K_w,0,0.1),ncol=S),
    #kappa_alpha=rnorm(sari_constants_NBz$K_t,0,0.1),#
    kappa_delta=matrix(rnorm(S*sari_constants_NBz$K_t,0,0.1),ncol=S),#
    #log_tau_alpha=rnorm(1,1.26,sd=2.23),#
    log_tau_delta=rnorm(S,1.26,sd=2.23),#
    iota=rnorm(S,0,1),#
    # log_theta=rnorm(S,4,1.3),
    log_tau_gamma=rnorm(S,1.26,sd=2.23),
    log_tau_xi=rnorm(S,1.26,sd=2.23),
    log_phi=matrix(rnorm(S*D_rem,4,1.3),nrow=S,ncol=D_rem),
    #log_tau_eta=rnorm(1,1.26,sd=2.23),
    #log_tau_psi=rnorm(1,1.26,sd=2.23),
    omega_beta=matrix(rnorm((D_rem-1)*S,0,0.1),nrow=S,ncol=(D_rem-1)))
  
  
  # Build the model.
  sari_model_gamma <- nimbleModel(sari_gamma_code,sari_constants_NBz,sari_data_NBz,sari_inits_NBz)
  # Compile the model.
  sari_compiled_model_gamma <- compileNimble(sari_model_gamma)
  # Set up the MCMC.
  sari_mcmc_config_gamma<- configureMCMC(sari_model_gamma, monitors=c('beta','xi','gamma','delta','iota', 'lambda', 'phi'),useConjugacy = FALSE)
  
  
  # objective function
  objective <- nimbleFunction(
    setup = function(model, target) {
      calcNodes <- model$getDependencies(target)     
      n_par<-length(values(model, target))
      
    },
    run = function(par = double(1)) { 
      returnType(double(0))
      values(model, target) <<- par
      ans <- model$calculate(calcNodes) 
      return(ans)
    }
  )
  
  
  
  # optimisation model run
  
  #paramers to optimise
  target_all <- c(sari_mcmc_config_gamma$getSamplers()%>%lapply(function(x)x$target)%>%unlist())
  n_par <- length(values(sari_model_gamma,target_all))
  
  #Initial vlaues for the optimisation
  opt_inits <- rnorm(n_par, 0, 0.01) #rep(0.01, n_par)
  #opt_inits[gamma_index]<-c(log(as.numeric(gamma_inits+0.05))) #log and add 0.05 as 0 values can't be logged
  
  rObjective <- objective(sari_model_gamma, target_all)
  cObjective <- compileNimble(rObjective, project = sari_model_gamma)
  
  #function for estimating gradient 
  grad.fun <- function(x){
    grad(func=cObjective$run, x=x, method="simple") #"Richardson" simple
  }
  
  #maximum iterations of optimisation 
  max_it<-5000
  
  # run and time optimisation
  
  optC_gamma<-trust.optim(fn = cObjective$run, gr= grad.fun,  x = opt_inits, control = list(function.scale.factor = -1, maxit=max_it), method = c("BFGS")) #other method: "SR1" "BFGS"
  
  # maxit=100 by default
  
  target<-target_all #don't need to run through all p_gammas
  par_z<-data.frame()
  par_z[1,1:length(target)]<-NA
  colnames(par_z)<-target
  model_values<-rep(0, n_par)
  for(j in 1:length(target)){   
    values(sari_model_gamma, target_all[j])<-NA
    par_index<-which(is.na(values(sari_model_gamma, target_all)))
    par_vector<-optC_gamma$solution[par_index]
    par_z[1:length(par_vector),j]<-par_vector
    values(sari_model_gamma, target_all)<-model_values
  }
  
  
  
  # get parameter estimates 
  
  # calculate beta
  beta_opt<-matrix(NA,nrow=S,ncol=(D_rem-1))
  beta_start1<-which(target_all=='omega_beta[1, 1]')
  beta_opt[,1]<-as.numeric(par_z[1,beta_start1:(beta_start1+S-1)])
  beta_start2<-beta_start1+S
  beta_opt[,2:(D_rem-1)]<-matrix(as.numeric(par_z[1,beta_start2:(beta_start2+(D_rem-2)*S-1)]), ncol=c(D_rem-2), byrow=TRUE)
  omega_beta_opt<-beta_opt
  for(s in 1:S){ 
    for(d in 2:(D_rem-1)){
      beta_opt[s,d]<-beta_opt[s,d-1]+exp(omega_beta_opt[s,d])
    }}
  
  
  #calculate gamma
  is_kappa_gamma<-sapply(target_all,function(x)substr(x,1,11))=='kappa_gamma'
  which_kappa_gamma<-as.numeric(which(is_kappa_gamma))
  kappa_gamma<- matrix(NA, nrow=sari_constants_NBz$K_t, ncol=S)
  for(s in 1:S){
    kappa_gamma[,s] <-as.numeric(na.omit(par_z[,which_kappa_gamma[s]]))
  }
  gamma_opt<-matrix(NA, nrow=N, ncol=S)
  for(s in 1:S){
    gamma_opt[,s] <- sari_data_NBz$X_t[1:N,1:sari_constants_NBz$K_t]%*%kappa_gamma[,s]
  }
  
  #calculate xi
  is_kappa_xi<-sapply(target_all,function(x)substr(x,1,8))=='kappa_xi'
  which_kappa_xi<-as.numeric(which(is_kappa_xi))
  kappa_xi<- matrix(NA, nrow=sari_constants_NBz$K_w, ncol=S)
  for(s in 1:S){
    kappa_xi[,s] <-as.numeric(na.omit(par_z[,which_kappa_xi[s]]))
  }
  xi_opt<-matrix(NA, nrow=52, ncol=S)
  for(s in 1:S){
    xi_opt[,s] <- sari_data_NBz$X_w[1:52,1:sari_constants_NBz$K_w]%*%kappa_xi[,s]
  }          
  
  # calculate phi
  phi_opt<-matrix(NA,nrow=S,ncol=D_rem)
  phi_start<-which(target_all=='log_phi[1, 1]')
  phi_opt<-matrix(exp(as.numeric(par_z[1,phi_start:(phi_start+D_rem*S-1)])), byrow=TRUE, ncol=D_rem, nrow=S)
  
  #calculate  delta's
  is_kappa_delta<-sapply(target_all,function(x)substr(x,1,11))=='kappa_delta'
  which_kappa_delta<-as.numeric(which(is_kappa_delta))
  kappa_delta<- matrix(NA, nrow=sari_constants_NBz$K_t, ncol=S)
  for(s in 1:S){
    kappa_delta[,s] <-as.numeric(na.omit(par_z[,which_kappa_delta[s]]))
  }
  delta_opt<-matrix(NA, nrow=N, ncol=S)
  for(s in 1:S){
    delta_opt[,s] <- sari_data_NBz$X_t[1:N,1:sari_constants_NBz$K_t]%*%kappa_delta[,s]
  }   
  
  #iota
  log_pop<-log(population[i])
  iota_start<-which(target_all=='iota[1]')
  iota_opt<-as.numeric(par_z[1,iota_start:(iota_start+S-1)])
  iota_opt_offset<-iota_opt-log_pop
  
  p_opt<-array(NA, dim=c(N,S,D_rem))
  mu_opt<-array(NA, dim=c(N,S,D_rem))
  lambda_opt<-array(NA, dim=c(N,S))
  for(t in 1:N){
    for(s in 1:S){
      #log(population[s]) + iota[s] + delta[t,s] + xi[w[t],s] #log(population[i])+
      lambda_opt[t,s]<-exp(iota_opt[s]+delta_opt[t,s]+xi_opt[week_index[t],s])
      for(d in 1:(D_rem-1)){
        p_opt[t,s,d]<-pnorm(beta_opt[s,d] + gamma_opt[t,s]) 
      }
      p_opt[t,s,D_rem]<-1
    }
  }
  for(s in 1:S){
    for(t in 1:N){
      mu_opt[t,s,1]<- p_opt[t,s,1]
      for(d in 2:(D_rem-1)){
        mu_opt[t,s,d] <- p_opt[t,s,d] - p_opt[t,s,(d-1)]
      }
      mu_opt[t,s,D_rem] <- 1 - p_opt[t,s,(D_rem-1)]
    }
  }
  
  #p_gamma 
  #  order<-numeric(N)
  #  pgam_name <- matrix(target_all[which_gamma],ncol=S)
  #  for(t in 1:N){
  #    p_ind<-paste('p_gamma[',t,',',' 1]', sep='')
  #   order[t]<-which(pgam_name[,1]==p_ind)
  #  }
  #  pgam_opt_unordered<-matrix(as.numeric(exp(optC_gamma$solution[gamma_index])),ncol=S)
  #  pgam_opt<-matrix(NA,ncol=S, nrow=N)
  #  for(t in 1:N){
  #    pgam_opt[t,]<-pgam_opt_unordered[order[t],]
  #  }
  
  
  z_opt<-z
  for(t in 1:C){
    for(s in 1:S){
      for(d in 1:D_rem){
        if(is.na(z_opt[t,s,d])){
          z_opt[t,s,d]<-lambda_opt[t,s]*mu_opt[t,s,d]  
        }
      }
    }
  }
  
  # start_theta<-which(target_all=='log_theta[1]')
  lambda_opt_gamma<-ceiling(lambda_opt)
  y_init_gamma<-rbind(ceiling(apply(z_opt,2,rowSums)),matrix(lambda_opt_gamma[(C+1):N,],ncol=S))
  # theta_opt_gamma<- exp(as.numeric(par_z[1,start_theta:c(start_theta+S-1)]))
  phi_opt_gamma<-phi_opt
  
  #other parameters
  start_tau_gamma<-which(target_all=='log_tau_gamma[1]')
  start_tau_delta<-which(target_all=='log_tau_delta[1]')
  start_tau_xi<-which(target_all=='log_tau_xi[1]')
  
  
  
  y_na<-y_init_gamma[1:C,1:S] #use optimised y
  not_na<-is.na(apply(censored_sari,c(1,2),sum)[1:C,i])==0
  y_na[not_na]<-rep(NA,sum(not_na)) #NA's for observed y's
  
  sari_inits <- list(#kappa_alpha=kappa_alpha,
    #kappa_psi=kappa_psi,
    #kappa_eta=kappa_eta,
    kappa_delta=kappa_delta,
    kappa_gamma=kappa_gamma,
    kappa_xi=kappa_xi,
    iota=iota_opt_offset,
    #log_theta=log(theta_opt_gamma),#
    log_tau_gamma=as.numeric(par_z[1,start_tau_gamma:c(start_tau_gamma+S-1)]),
    log_tau_delta=as.numeric(par_z[1,start_tau_delta:c(start_tau_delta+S-1)]),
    log_tau_xi=as.numeric(par_z[1,start_tau_xi:c(start_tau_xi+S-1)]),
    log_phi=log(phi_opt[1:S,1:D]),
    #log_tau_psi=log_tau_psi,
    #log_tau_alpha=log_tau_alpha,
    #log_tau_eta=log_tau_eta,
    omega_beta=omega_beta_opt[1:S,1:D],
    y=y_na)
  list(sari_inits=sari_inits, y_init_gamma=y_init_gamma,lambda_opt=lambda_opt, p_opt=p_opt, mu_opt=mu_opt)
  
}

time_optim<-system.time({
  cl <- makeCluster(4)
  save_nogam <- parLapply(cl, X=1:S_all, fun=f, N=N, C=C,D=2, S=1, censored_sari=censored_sari, population=population)
  stopCluster(cl)
})

#Plot optimisation 
stat_t=1
for(i in 1:S_all){
  # png(paste('Sari opt Sari opt plots/R',i,'sari_no_gam.png'), width = 700, height = 350)
  plot(stat_t:N, reduced_y[stat_t:N,i])
  title(paste(i, '- No gamma phi '))
  lines(stat_t:N, save_nogam[[i]]$lambda_opt[stat_t:N,], col='blue') #lambda optimised
  abline(v=C)
  abline(v=199)
  lines(stat_t:N, save_nogam[[i]]$y_init_gamma[stat_t:N,], col='green') #z's sum
  #dev.off() 
}

## SET INTIAL VALUES 
kappa_delta<-matrix(NA, nrow=nrow(save_nogam[[1]]$sari_inits$kappa_delta), ncol=S_all)
kappa_gamma<-matrix(NA, nrow=nrow(save_nogam[[1]]$sari_init$kappa_gamma), ncol=S_all)
kappa_xi<-matrix(NA, nrow=nrow(save_nogam[[1]]$sari_inits$kappa_xi), ncol=S_all)
iota<-rep(NA,S_all)
#log_theta<-rep(NA,S_all)
log_tau_gamma<-rep(NA,S_all)
log_tau_delta<-rep(NA,S_all)
log_tau_xi<-rep(NA,S_all)
log_phi<-matrix(NA, nrow=S_all, ncol=D)
omega_beta<-matrix(NA, nrow=S_all, ncol=D)
y<-matrix(NA, nrow=C, ncol=S_all)

for(i in 1:S_all){
  kappa_delta[,i]<-save_nogam[[i]]$sari_inits$kappa_delta
  kappa_gamma[,i]<-save_nogam[[i]]$sari_inits$kappa_gamma
  kappa_xi[,i]<-save_nogam[[i]]$sari_inits$kappa_xi
  iota[i]<-save_nogam[[i]]$sari_inits$iota
  #log_theta[i]<-save_nogam[[i]]$sari_inits$log_theta
  log_tau_gamma[i]<-save_nogam[[i]]$sari_inits$log_tau_gamma
  log_tau_delta[i]<-save_nogam[[i]]$sari_inits$log_tau_delta
  log_tau_xi[i]<-save_nogam[[i]]$sari_inits$log_tau_xi
  log_phi[i,]<-save_nogam[[i]]$sari_inits$log_phi
  omega_beta[i,]<-save_nogam[[i]]$sari_inits$omega_beta
  y[,i]<-save_nogam[[i]]$sari_inits$y
}

kappa_eta <-apply(kappa_xi, 1, mean)
kappa_alpha <-apply(kappa_delta, 1, mean)
kappa_psi <-apply(kappa_gamma, 1, mean)
log_tau_eta<-mean(log_tau_xi)
log_tau_alpha<-mean(log_tau_delta)
log_tau_psi<-mean(log_tau_gamma)

sari_inits <- list(kappa_alpha=kappa_alpha,
                   kappa_psi=kappa_psi,
                   kappa_eta=kappa_eta,
                   kappa_delta=kappa_delta,
                   kappa_gamma=kappa_gamma,
                   kappa_xi=kappa_xi,
                   iota=iota,
                   log_theta=rnorm(S,4,1.3),#
                   log_tau_gamma=log_tau_gamma,
                   log_tau_delta=log_tau_delta,
                   log_tau_xi=log_tau_xi,
                   log_phi=log_phi,
                   log_tau_psi=log_tau_psi,
                   log_tau_alpha=log_tau_alpha,
                   log_tau_eta=log_tau_eta,
                   omega_beta=omega_beta,
                   y=y)

###### GDM MODEL? ###########

### FIT GDM
### set up cluster

S=22 # Number of regions to model.
D=2 # Number of delays to explicit
n_chains <- 4 # Number of MCMC chains to run. Reduce if you don't have more than 8 CPU cores!
n_knots <- c(18,9) # Number of knots of the temporal and seasonal splines, respectively.
N <- 230 # Number of weeks to consider.
C <- 224 # Present-day week.

week_index=rep(1:52,10) # Vector for index which week in the year each week is.

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
sari_nona_cluster<- makeCluster(n_chains)

### cluster function 

run_sari_nona_cluster<- function( seed, sari_constants, sari_data, sari_inits ) {
  
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
  
  # Objective function for optimisation
  objective <- nimbleFunction(
    setup = function(model, target) {
      calcNodes <- model$getDependencies(target)     
      n_par<-length(values(model, target))
      
    },
    run = function(par = double(1)) { 
      returnType(double(0))
      values(model, target) <<- par
      ans <- model$calculate(calcNodes) 
      return(ans)
    }
  )
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
  
  # Build the model.
  sari_model <- nimbleModel(sari_code,sari_constants,sari_data,sari_inits)
  # Compile the model.
  sari_compiled_model<- compileNimble(sari_model)
  # Set up the MCMC.
  sari_mcmc_config <- configureMCMC(sari_model,monitors=c('alpha','eta','iota','beta','psi',
                                                          'theta','phi','gamma','delta','xi','y','lambda'))
  
  
  ## OPTIMISE DISPERSION PARAMETERS 
  target_full <- c(sari_mcmc_config$getSamplers()%>%lapply(function(x)x$target)%>%unlist())
  is_theta <- sapply(target_full,function(x)substr(x,1,7))=='log_the'
  is_phi <- sapply(target_full,function(x)substr(x,1,7))=='log_phi'
  which_theta<-which(is_theta)
  which_phi<-which(is_phi)
  n_phi<-length(which_phi)
  target_disp<-target_full[sort(c(which_theta,which_phi))]
  n_par2 <- length(values(sari_model,target_disp))
  
  #initial values for the optimisation
  opt_inits2 <- rnorm(n_par2,0.01,0.1)
  
  rObjective <- objective(sari_model, target_disp)
  cObjective <- compileNimble(rObjective, project = sari_model)
  
  #function for estimating gradient 
  grad.fun <- function(x){
    grad(func=cObjective$run, x=x, method="simple") #"Richardson" simple
  }
  #optimisation function
  optC_disp<-trust.optim(fn = cObjective$run, gr= grad.fun,  x = opt_inits2,
                         control = list(function.scale.factor = -1, maxit=500),
                         method = c("BFGS")) #hs= #only needed for sparse #other method: "SR1" "BFGS"
  
  
  #set dispersion parameters 
  theta_opt_gamma2<- as.numeric(exp(optC_disp$solution[(n_phi+1):n_par2]))
  phi_opt_gamma2<-as.numeric(exp(optC_disp$solution[1:n_phi]))%>%matrix(nrow=S, ncol=D, byrow=TRUE)
  
  sari_inits$log_theta<-as.numeric(optC_disp$solution[(n_phi+1):n_par2])
  sari_inits$log_phi<-log(phi_opt_gamma2)#
  sari_model$setInits(list(log_theta=as.numeric(optC_disp$solution[(n_phi+1):n_par2]), log_phi=log(phi_opt_gamma2)))
  
  
  sari_mcmc <- buildMCMC(sari_mcmc_config)
  # Compile the MCMC.
  sari_compiled_mcmc <- compileNimble(sari_mcmc,project=sari_model)
  sari_samples <- runMCMC(sari_compiled_mcmc,niter=750000,nburnin=50000,inits=sari_inits,nchains=1,samplesAsCodaMCMC = TRUE,thin=1000)
  
  list(sari_samples=sari_samples,theta_opt=theta_opt_gamma2, phi_opt=phi_opt_gamma2)
}


time_nona_cluster<-system.time({
  sari_samples <- parLapply(cl = sari_nona_cluster, X = 1:n_chains, 
                            fun = run_sari_nona_cluster, sari_inits=sari_inits,
                            sari_constants = sari_constants, sari_data=sari_data)
})

stopCluster(sari_nona_cluster)

save(sari_samples, file='SARI_samples_optim.RData')

sari_sampes_mcmc<-list()
for(j in 1:n_chains){
  # Survivor age model
  sari_sampes_mcmc[[j]]<-unlist(sari_samples[[j]]$sari_samples)
  #covid_combined_samples_mcmc[[j]] <- as_tibble(do.call('rbind',SARI_Survivor_linear_x[[j]]))
}
sari_sampes_mcmc<-as.mcmc.list(sari_sampes_mcmc)
sari_combined_samples <- as_tibble(do.call('rbind',sari_sampes_mcmc))


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
#ggsave(effect_plots,file='Sari opt plots/sari_effect_plots.pdf',width=9,height=2.5)  
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

##### OPTIMISATION PLOTS 

lambda_NB<-matrix(NA, nrow = N, ncol=S_all)
y_init_region<-matrix(NA, nrow = N, ncol=S_all)
for(i in 1:S_all){
  y_init_region[,i]<-save_nogam[[i]]$y_init_gamma
  lambda_NB[,i]<-save_nogam[[i]]$lambda_opt
}

opt_df<-cbind(melt(y_init_region),melt(lambda_NB)[3])
colnames(opt_df)<-c("t","s","y_gamma","lambda_NB")

colours <- c("MCMC" = "orange", "NB Optimisation" = "blue", "GDM Optimisation" = "dark green")
fill.col <- c("95%" = "orange")

opt_df$r<-region_names[opt_df$s]
y_data$r<-region_names[y_data$s]
y_data<-y_data[-which(y_data$s==8),]
sari_y<-sari_y[-which(sari_y$s==8),]



region_all <- ggplot(filter(sari_y,s%in%c(2,17,15)))+
  geom_ribbon(aes(x=dates[t],ymin=`2.5%`,ymax=`97.5%`,fill="95%"),alpha=0.25)+
  geom_line(aes(x=dates[t],y=`50%`,colour=c("MCMC")))+
  geom_line(data=filter(opt_df,s%in%c(2,17,15)), aes(x=dates[t],y=`y_gamma`,colour=c("NB Optimisation")))+
  geom_line(data=filter(opt_df,s%in%c(2,17,15)), aes(x=dates[t],y=`lambda_NB`,colour=c("NB Optimisation")), linetype="dashed")+
  geom_point(data=filter(y_data,s%in%c(2,17,15)),aes(x=dates[t],y=y))+
  labs(x=NULL,y=NULL,title='Predicted weekly SARI cases',
       subtitle='For three regional health disricts in Paraná, Brazil',
       caption=NULL, colour="Model", fill="Prediction Interval")+
  coord_cartesian(xlim=c(dates[200],dates[230]))+geom_vline(xintercept=dates[C])+
  theme_minimal()+
  facet_wrap(~r,ncol=3,scales='free')+
  scale_fill_manual(values=fill.col)+
  scale_color_manual(values=colours)+
  theme(text = element_text(size=14), axis.title = element_text(size=16), plot.title = element_text(size=20),plot.subtitle=element_text(size=18),
        legend.title = element_text(size=16),strip.text = element_text(size=14),
        axis.text = element_text(size=9),legend.position="bottom")+guides(colour=guide_legend(nrow=1))+
  scale_x_date(labels = date_format("%b %Y"),date_breaks = '2 months')

region_all
ggsave(region_all,filename = 'Sari opt plots/sari_region_all_5000.pdf',width=10,height=4)


#ggsave(region_all,filename = 'Sari opt plots/Covid_SARI_5000_l.png',width=12,height=3)



# THETA

disp_df<-data.frame(sari_samples[[1]]$theta_opt)
colnames(disp_df)<-c("theta_opt_gamma2")
disp_df$theta_median<-apply(sari_theta,2,quantile, 0.5)
disp_df$theta_025<-apply(sari_theta,2,quantile, 0.025)
disp_df$theta_975<-apply(sari_theta,2,quantile, 0.975)
disp_df$s<-c(1:S_all)

theta_opt<-ggplot(data=disp_df)+
  geom_line(aes(x=c(1:S),y=`theta_opt_gamma2`,colour=c("GDM Optimisation")))+
  #geom_line(aes(x=c(1:7),y=`theta_opt_gamma_og`,colour=c("NB Optimisation")))+
  geom_line(aes(x=c(1:S),y=`theta_median`,colour=c("MCMC")))+
  geom_ribbon(aes(x=c(1:S), ymin=`theta_025`, ymax=`theta_975`, fill=c("95%")),alpha=0.25)+
  labs(x='Regions',y=NULL,title='Theta comparison',
       caption=NULL, colour="Model", fill="Prediction Interval")+ theme_minimal()+
  scale_fill_manual(values=fill.col)+
  scale_color_manual(values=colours)+
  theme(strip.text = element_text(size=12),
        plot.caption = element_text(size=10),plot.title=element_text(size=16),
        axis.text = element_text(size=10),axis.title = element_text(size=12))

ggsave(theta_opt,filename = 'Sari opt plots/SARI_theta_opt_GDM_5000.png',width=8,height=4)
theta_opt

#phi

phi_opt_gamma2<-matrix(sari_samples[[1]]$phi_opt, nrow=S_all, ncol=D, byrow=TRUE)

phi_df<-data.frame(cbind(melt(phi_opt_gamma2)))
colnames(phi_df)<-c("r","d","phi_opt_gamma")
phi_mcmc<-cbind(melt(apply(sari_phi,c(2,3),quantile, 0.5)),
                melt(apply(sari_phi,c(2,3),quantile, 0.025))[,3],
                melt(apply(sari_phi,c(2,3),quantile, 0.975))[,3])
colnames(phi_mcmc)<-c("r","d","phi_median","phi_025","phi_975")
phi_df<-full_join(phi_df,phi_mcmc,by=c("r","d"))%>%mutate(r=factor(r))

phi_opt<-ggplot(data=filter(phi_df))+ 
  geom_line(aes(x=d,y=`phi_opt_gamma`,colour=c("GDM Optimisation")))+
 # geom_line(aes(x=d,y=`phi_opt_gamma_og`,colour=c("NB Optimisation")))+
  geom_line(aes(x=d,y=`phi_median`,colour=c("MCMC")))+
  geom_ribbon(aes(x=d, ymin=`phi_025`, ymax=`phi_975`, fill=c("95%")),alpha=0.25)+
  labs(x='delays',y=NULL,title='Phi comparison',
       caption=NULL, colour="Model", fill="Prediction Interval")+ theme_minimal()+
  scale_fill_manual(values=fill.col)+
  scale_color_manual(values=colours)+
  facet_wrap(~r,nrow=2,scales='free')+
  theme(strip.text = element_text(size=12),
        plot.caption = element_text(size=10),plot.title=element_text(size=16),
        axis.text = element_text(size=10),axis.title = element_text(size=12))
ggsave(phi_opt,filename = 'Sari opt plots/SARI_phi_opt_5000.png',width=21,height=4)
phi_opt

