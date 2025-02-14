############################################
#  Simulation Experiment: DATA GENERATION  #
############################################

library(nimble)
library(compositions)
library(tidyverse)
library(mgcv)

# Beta-Binomial (rbetabin) Function
rbetabin=nimbleFunction(run=function(n=integer(0),mu=double(0),phi=double(0),size=double(0)){
  pi=rbeta(n,mu*phi,(1-mu)*phi)
  returnType(double(0))
  return(rbinom(n,size,pi))
})


##### FIT GDM TO DATA #####
# Create function for cluster to carry out MCMC model:
Simulation_data<- function(seed, S, N, D){
  # Load libraries.
  library(nimble)
  library(tidyverse)
  library(mgcv)
  library(compositions)
  
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
  #set.seed(seed)
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
  
  # Spatial intercept of mean reported cases - iota 
  sim_iota<- rep(rnorm(1,5,sd=0.25),S)
  # Overall temporal trend in total cases for all zeta values
  # Simulate mean number of total cases - lambda
  sim_lambda <- array(dim=c(N,S))
  sign<-(rbernoulli(1,0.5)-1)
  sign[sign==0]=1
  #sign
  for(s in 1:S){
    for(t in 1:N){
      sim_lambda[t,s]<-exp(sim_iota[s]+sim_alpha[t]*sign)
                           
    }
  }
  # Simulate Negative Binomial dispersion parameter.
  sim_theta<-rep(rgamma(1,2,0.02),S)
  # Negative Binomial distribution for total counts (y).
  sim_y<-array(dim=c(N,S))
  for(s in 1:S){
    sim_y[,s]<-rnbinom(N,mu=sim_lambda[,s],size=sim_theta[s])
  }
  
  # GDM simulation model for partial reports (z).
  # Independent delay effects (psi). 
  sim_psi<-matrix(rep(rnorm(D,0,sd=0.25),S), nrow=S, ncol=D)
  # for(s in 1:S){
  #   sim_psi[s,]<-sort(sim_psi[s,])
  # }
  
  # Linear temporal trend in reporting performance over time (eta)
  sim_eta<-matrix(nrow=N,ncol=S)
  # Magnitude of linear trend in the delay distribution:
  eta_slopes<-rnorm(1,0.2,sd=0.1)
  for(s in 1:S){
    for(t in 1:N){
      sim_eta[t,s]<-eta_slopes*((t-mean(1:N))/sd(1:N))
    }
  }
  
  # Linear relationship between mean cases (lambda) and delay length.
  # delta_slopes <- rep(c(rep(-0.2,3),rep(0,3),rep(0.2,3)),3)
  # sim_delta<-array(dim=c(N,S))
  # for(s in 1:S){
  #   for(t in 1:N){
  #     sim_delta[t,s]<-delta_slopes[s]*(sim_y[t,s]-mean(sim_y[,s]))/sd(sim_y[,s])
  #   }
  # }
  # sim_p<-array(dim=c(N,S,D))
  # for(d in 1:D){
  #   for(s in 1:S){
  #     sim_p[,s,d] <-  iprobit(sim_eta[,s]+sim_psi[s,d])#+sim_delta[,s])
  #   }}
  
  # Linear predictor (Beta-Binomial means).
   sim_u<-array(dim=c(N,D+1,S))
   sim_nu<-array(dim=c(N,S,D))
  for(s in 1:S){
    for(t in 1:N){
  for(d in 1:D){
    sim_u[t,d,s]<- sim_eta[t,s]+sim_psi[s,d]
  }
  sim_u[t,D+1,s]<- -sum(sim_u[t,1:D,s])
  sim_ap<-clrInv(sim_u[t,1:D,s])
  # Relative proportions (Beta-Binomial means).
  for(d in 1:D){
    sim_nu[t,s,d] <- exp(sim_u[t,d,s])/(sum(exp(sim_u[t,d:(D+1),s])))
  }
  }}
   
  # Relative proportions (eta-Binomial means - nu).
  # for(s in 1:S){
  #   sim_nu[,s,1] <- sim_p[,s,1]
  #   for(d in 2:D){
  #     sim_nu[,s,d] <- (sim_p[,s,d]-sim_p[,s,d-1])/(1-sim_p[,s,d-1])
  #   }}
  
  # eta-Binomial dispersion parameters (phi).
  # sim_phi<-matrix(rep(rgamma(D,2,0.02),S), nrow=S, ncol=D, byrow=TRUE) 
  sim_phi<-seq(from=2, to=10, length=S)
  sim_z<-array(dim=c(N,S,D))
  for(t in 1:N){
    for(s in 1:S){
      sim_z[t,s,1] <- rbetabin(1,sim_nu[t,s,1],sim_phi[s],sim_y[t,s])
      for(d in 2:(D)){
        sim_z[t,s,d] <- rbetabin(1,sim_nu[t,s,d],sim_phi[s],sim_y[t,s]-sum(sim_z[t,s,1:(d-1)]))
      }
    }
  }
  
  simulation_data <- list(y=sim_y, z=array(sim_z[,1:S,1:D],dim=c(N,S,D)), p=sim_nu, lambda=sim_lambda, phi=sim_phi,eta=eta_slopes)
  
  
  return(simulation_data)
}
