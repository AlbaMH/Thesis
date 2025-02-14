################################################
#  Simulation Experiment: PREDICTION EXPERIMENT  #
################################################

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
library(ggh4x)
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


##### FIT GDM TO SIMULATED DATA #####

# Create function for cluster to carry out MCMC model:
Simulation_survivor<- function(X, sim_data, S, N, D, n_iter, n_burn, n_thin, n_chains){
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
  
  delta_converter=nimbleFunction(run=function(x=double(0),grid=double(1),grid_length=double(0)){
    dif=abs(x-grid)
    min_dif=min(dif)
    min_index=dif==min_dif
    index_delta=((1:grid_length)[min_index])
    
    returnType(double(1))
    return(index_delta)
  })
  
  
  assign('dbetabin', dbetabin, envir = .GlobalEnv)
  assign('rbetabin', rbetabin, envir = .GlobalEnv)
  assign('delta_converter', delta_converter, envir = .GlobalEnv)
  
  registerDistributions(list(dbetabin=list(
    BUGSdist='dbetabin(mu,phi,size)',discrete=TRUE)))
  
  # NIMBLE model code. 
  simulation_code <- nimbleCode({
    for(s in 1:S){
      for(t in 1:N){
        # Negative Binomial Model for total simulation cases.
        y[t,s] ~ dnegbin(prob=theta[s]/(theta[s]+lambda[t,s]),size=theta[s])  
        
        # Model for partial delayed simulation counts (first delay).
        z[t,s,1] ~ dbetabin(nu[t,1,s],phi[1,s],y[t,s])
      }
      for(d in 2:D){
        for(t in 1:obs_index[s,d]){
          # Model for partial delayed simulation counts.
          z[t,s,d] ~ dbetabin(nu[t,d,s],phi[d,s],y[t,s]-sum(z[t,s,1:(d-1)]))
        }
      }
      for(t in 1:N){
        # Linear predictor (beta-Binomial means).
        for(d in 1:D){
          # Expected cumulative proportions.
          probit(p[t,d,s]) <- psi[d,s] + eta_slope[s]*((t-mean(1:N))/sd(1:N)) + delta_slope[s]*((y[t,s]-mean(y[,s]))/sd(y[,s]))
        }
        # Relative proportions (beta-Binomial means).
        nu[t,1,s] <- p[t,1,s]
        for(d in 2:D){
          nu[t,d,s] <- (p[t,d,s]-p[t,d-1,s])/(1-p[t,d-1,s])
        }
        # Mean for total simulation cases .
        log(lambda[t,s]) <- iota[s] + alpha_spline[t,s] 
        
      }
      
      
      # Linear coefficients 
    #  delta_slope[s] ~ dnorm(0,sd=10) 
      eta_slope[s] ~ dnorm(0,sd=10)
      
      #  Temporal spline for expected mean cases.
      Omega_alpha[1:(K_t),1:(K_t),s] <- S_t[1:(K_t),1:(K_t)]/sigma_alpha[s,1]^2+S_t[1:(K_t),(K_t+1):(2*K_t)]/sigma_alpha[s,2]^2
      kappa_alpha[1:(K_t),s] ~ dmnorm(zeros[1:(K_t)],Omega_alpha[1:(K_t),1:(K_t),s])
      alpha_spline[1:N,s] <- X_t[1:N,1:(K_t)]%*%kappa_alpha[1:(K_t),s]
      
      sigma_alpha[s,1] ~ T(dnorm(0,sd=10),0,)
      sigma_alpha[s,2] ~ T(dnorm(0,sd=10),0,)
      
      
      # Cumulative delay curves.
      psi[1,s] ~ dnorm(-1.5,sd=5)
      for(d in 2:D){
        psi[d,s] ~ T(dnorm(psi[d-1,s],sd=5),psi[d-1,s],)
      }
      
      # beta-Binomial dispersion parameters.
      for(d in 1:D){
        phi[d,s] ~ dgamma(2,0.02) 
      }
      theta[s] ~ dgamma(2,0.02) # Negative-Binomial dispersion parameters.
      
      # Intercepts
      iota[s]  ~ dnorm(5,sd=10)
    }
    
  })
  
  # Generate spline bases functions for modelling and simulation.
  blank_data=tibble(y=rnorm(N,0,1),t=1:N)
  blank_jagam=jagam(y~s(t,k=floor(N/5), bs='tp'),
                    data=blank_data,file='blank.jags',
                    knots=list(t=seq(1,N,length=floor(N/5))))
  
  # Simulated data
  sim_y<-sim_data[[X]]$y
  sim_z<-sim_data[[X]]$z
  # Constants from simulated data:
  N<-dim(sim_y)[1]
  S<-dim(sim_y)[2]
  D<-dim(sim_z)[3]
  
  z<- sim_z[,1:S,1:D]
  
  # Reported partial counts:
  censored_z <- z
  
  for(s in 1:S){
    censored_z[,s,][outer(1:dim(z[,s,])[1], 0:(dim(z[,s,])[2]-1), FUN = "+") > N] <- NA
  }
  
  # Maximum observed index for each region at each delay.
  obs_index<-array(NA, dim=c(S,D))
  for(s in 1:S){
    for(d in 1:D){
      obs_index[s,d]<-which(is.na(censored_z[,s,d])==TRUE)[1]-1
    }
  }
  
  
  # Censor the totals.
  censored_y <- sim_y
  censored_y[(N-D+1):N,]<-NA
  # Model Data
  simulation_data <- list(y=censored_y[1:N,1:S], z=censored_z[1:N,1:S,1:D])
  
  # Constants
  simulation_constants <-list(N=N,  D=D, S=S, obs_index=obs_index,
                              delta_slope=sim_data[[X]]$delta,
                              K_t=dim(blank_jagam$jags.data$S1)[1])
  simulation_constants$S_t <- blank_jagam$jags.data$S1
  simulation_constants$X_t <- blank_jagam$jags.data$X[,2:(simulation_constants$K_t+1)]
  simulation_constants$zeros <- rep(0,max(simulation_constants$K_t))
  
  #simulation_inits<-simulation_model<-simulation_compiled_model<- simulation_mcmc_config<-simulation_mcmc <- simulation_compiled_mcmc<-list()
  # Generate random initial values.
  simulation_inits <- list(
  #  delta_slope=rnorm(S,0,sd=0.1),
    kappa_alpha=matrix(rnorm(S*(simulation_constants$K_t),0,0.01),ncol=S),
    sigma_alpha=matrix(runif(S*2,0,1),ncol=2),
    # zeta_slope=rnorm(S,0,sd=0.1),
    #kappa_zeta=matrix(rnorm(S*(simulation_constants$K_t),0,0.01),ncol=S),
    eta_slope=rnorm(S,0,sd=0.1),
    psi=matrix(rep(sort(rnorm(D,0,1)),S),ncol=S, byrow = FALSE),
    theta=abs(rnorm(S,50,sd=10)),
    iota=rnorm(S,5,1),
    phi=matrix(abs(rnorm(D*S,50,sd=10)),ncol=S),
    y=array(dim=c(N,S)) )
  # Set initial values for total counts y.
  for(t in 1:N){
    for(s in 1:S){
      if(is.na(censored_y[t,s])) simulation_inits$y[t,s]=sum(simulation_data$z[t,s,],rpois(1,median(simulation_data$y[,s]-rowSums(simulation_data$z[,s,]),na.rm=T)),na.rm=TRUE)
    }
  }
  
  
  # Build the model.
  simulation_model <- nimbleModel(simulation_code,simulation_constants,simulation_data,simulation_inits)
  # Compile the model.
  simulation_compiled_model <- compileNimble(simulation_model)
  
  # Set up the MCMC.
  simulation_mcmc_config <- configureMCMC(simulation_model,monitors=c("lambda","iota","alpha_spline",
                                                                      "psi","theta","y",
                                                                      "eta_slope","kappa_alpha","sigma_alpha"
  ),
  useConjugacy = FALSE)
  
  simulation_mcmc_config$removeSamplers('eta_slope','psi','iota','kappa_alpha','sigma_alpha')#
  
  for(s in 1:S){
    simulation_mcmc_config$addSampler(target=c(#paste('delta_slope[',s,']', sep=''), 
      paste('eta_slope[',s,']',sep=''), 
      paste('psi[1:',D,', ',s,']', sep=''),
      paste('sigma_alpha[',s,', 1:2]', sep=''),
      paste('kappa_alpha[1:',(simulation_constants$K_t),', ',s,']',sep=''),
      paste('iota[',s,']',sep='')),type='AF_slice')
    
  }
  
  # Build MCMC. 
  simulation_mcmc<- buildMCMC(simulation_mcmc_config)
  # Compile the MCMC.
  simulation_compiled_mcmc <- compileNimble(simulation_mcmc,project=simulation_model)
  # Run MCMC. 
  simulation_cluster<- runMCMC(simulation_compiled_mcmc,niter=n_iter,nburnin=n_burn,thin=n_thin,inits=simulation_inits,nchains=n_chains,samplesAsCodaMCMC = TRUE)
  
  return(simulation_cluster)
}


# MCMC parameters.
# Set chains.
n_chains=2
# Set iterations.
n_iter=15000
# Set burn-in.
n_burn=10000
# Set thinning.
n_thin=5
# Set number of knots for model splines
#n_knots<-# knots for time, time-delay interaction and delay respectively 
# Load simulated data sets:
load("~/media/alba/Disk 1/OneDrive/PhD/Case load simulations/Simulation studies/new data/Simulation_data_y.RData")

# Number of simulations.
sims=100
# Set maximum COVID-19 lab report delay.
n_cores<-min(detectCores(),sims,25)

# Make Cluster for MCMC. 
this_cluster_simulation_fixed<- makeCluster(n_cores)

# Run Cluster.
time_surv_sim_fixed <- system.time({
  simulation_surv_sim_fixed <-  parLapply(cl = this_cluster_simulation_fixed, X = 1:sims, 
                                          fun = Simulation_survivor,
                                          sim_data=simulation_data_list,
                                          S=S,
                                          N=N,
                                          D=D,
                                          n_iter=n_iter,
                                          n_burn=n_burn,
                                          n_thin=n_thin,
                                          n_chains=n_chains)
})

# Stop Cluster.
stopCluster(this_cluster_simulation_fixed)

# Save survivor model output. 
save(simulation_surv_sim_fixed,time_surv_sim_fixed, file='simulation_parameter_fixed_delta.RData')

# Check PSRF of model parameters.
parameter_group=c("lambda","iota","alpha_spline",
                  "psi","theta",#"delta_slope",
                  "eta_slope","kappa_alpha","sigma_alpha")
sari_psrf<-list()
sari_mean_psrf<-sari_max_psrf<-sari_leq_105 <-sari_leq_110<-array(NA,dim=c(length(parameter_group),sims))
for(k in 1:length(parameter_group)){
  sari_psrf[[k]]<-list()  
  for(j in 1:sims){
    sari_parameter_names <- colnames((as.matrix(simulation_surv_sim_fixed[[j]])))
    sari_index <- which(startsWith(sari_parameter_names,parameter_group[k]))
    sari_psrf[[k]][[j]]<- sapply(sari_index,function(v)gelman.diag(simulation_surv_sim_fixed[[j]][,v],autoburnin = FALSE,transform = TRUE,multivariate = FALSE)$psrf[,1])
    sari_mean_psrf[k,j] <- mean(sari_psrf[[k]][[j]])
    sari_max_psrf[k,j] <- max(sari_psrf[[k]][[j]])
    sari_leq_105[k,j] <- mean(sari_psrf[[k]][[j]]<=1.05)
    sari_leq_110[k,j] <- mean(sari_psrf[[k]][[j]][j]<=1.1)
  }
}
sari_MAX_psrf<-cbind(parameter_group,apply(sari_max_psrf,1,max))
sari_MEAN_psrf<-cbind(parameter_group,apply(sari_mean_psrf,1,mean))



plot(simulation_surv_sim_fixed[[1]][,"delta_slope[1]"],main="delta_slope")
plot(simulation_surv_sim_fixed[[2]][,"delta_slope[2]"],main="delta_slope")
plot(simulation_surv_sim_fixed[[1]][,"delta_slope[3]"],main="delta_slope")
plot(simulation_surv_sim_fixed[[1]][,"delta_slope[4]"],main="delta_slope")
plot(simulation_surv_sim_fixed[[2]][,"delta_slope[5]"],main="delta_slope")
plot(simulation_surv_sim_fixed[[1]][,"delta_slope[6]"],main="delta_slope")
plot(simulation_surv_sim_fixed[[2]][,"delta_slope[7]"],main="delta_slope")
plot(simulation_surv_sim_fixed[[1]][,"delta_slope[8]"],main="delta_slope")
plot(simulation_surv_sim_fixed[[2]][,"delta_slope[9]"],main="delta_slope")
# plot(simulation_surv_sim_fixed[[3]][,"delta_slope[7]"],main="delta_slope")
# plot(simulation_surv_sim_fixed[[5]][,"delta_slope[7]"],main="delta_slope")
# plot(simulation_surv_sim_fixed[[6]][,"delta_slope[7]"],main="delta_slope")
# plot(simulation_surv_sim_fixed[[7]][,"delta_slope[7]"],main="delta_slope")
# plot(simulation_surv_sim_fixed[[9]][,"delta_slope[7]"],main="delta_slope")
# plot(simulation_surv_sim_fixed[[10]][,"delta_slope[7]"],main="delta_slope")
# 
# plot(simulation_surv_sim_fixed[[8]][,"delta_slope[7]"],main="delta_slope")
# plot(simulation_surv_sim_fixed[[4]][,"delta_slope[7]"],main="delta_slope")
# 
# plot(simulation_surv_sim_fixed[[1]][,"alpha[1, 1]"],main="zeta[1] 1")
# plot(simulation_surv_sim_fixed[[2]][,"zeta[50]"],main="zeta[50] 2")
# plot(simulation_surv_sim_fixed[[3]][,"zeta[38]"],main="zeta[38] 3")
# 
# plot(simulation_surv_sim_fixed[[1]][,"eta_slope[1]"],main="eta_slope 1")
# plot(simulation_surv_sim_fixed[[2]][,"eta_slope"],main="eta_slope 2")
# plot(simulation_surv_sim_fixed[[3]][,"eta_slope"],main="eta_slope 3")

simulation_mcmc<-list()
simulation_combined<-list()
for(j in 1:sims){
  simulation_mcmc[[j]]<-as.mcmc.list(simulation_surv_sim_fixed[[j]])
  # Combine all MCMC chains.
  simulation_combined[[j]] <- as_tibble(do.call('rbind',simulation_mcmc[[j]]))
}

# Model coefficient output
n_sim<-dim(simulation_combined[[1]])[1]
delta_slope<-zeta_slope<-eta_slope<-list()
eta_ind<-rep(c('Improving','Not changing','Worsening'),S/3)
delta_ind<-c(rep('Negative effect',S/3),rep('No effect',S/3),rep('Positive effect',S/3))
zeta_ind<-c(rep('Improving',S/3),rep('Not changing',S/3),rep('Worsening',S/3))

# Simulated parameter values:
eta_value<-simulation_data_list[[1]]$eta
delta_value<-simulation_data_list[[1]]$delta
zeta_value<-simulation_data_list[[1]]$zeta
actual_delta <- tibble(delta_slopes=as.numeric(delta_value),s=(1:S))
actual_zeta <- tibble(zeta_slopes=zeta_value,s=1:S)
actual_eta <-tibble(eta_slopes=as.numeric(eta_value),s=(1:S))



# y_fixed[1,100,1]
# simulation_combined[[j]][,7264]
# y_quant_fixed[[j]]<-simulation_combined[[j]][,y_cols]%>%as.matrix()%>%array(dim=c(n_sim,N,S))
# y_quant_fixed[[j]][,100,10]


alpha<-kappa_alpha<-lambda<-alpha_spline<-y_quant_fixed<-list()
for(j in 1:sims){
  y_cols<-which(str_detect(colnames(simulation_combined[[j]]),c('y')))
  y_quant_fixed[[j]]<-simulation_combined[[j]][,y_cols]%>%as.matrix()%>%array(dim=c(n_sim,N,S))%>%
    apply(c(2,3),quantile,c(0.025,0.5,0.975))%>%melt(varnames=c('quantile','t','s'),value.name='y')%>%#mutate(quantile=c(0.025,0.5,0.975))%>%
    spread(quantile,y)%>%mutate(zeta_ind=zeta_ind[s],eta_ind=eta_ind[s],delta_ind=delta_ind[s],
                                zeta=zeta_value[s],eta=eta_value[s],delta=delta_value[s],
                                sim=j)
  # y_quant[[j]]<- full_join(y_slope[[j]],simulation_data_list[[j]]$y, by='s')
  
  # delta_cols<-which(str_detect(colnames(simulation_combined[[j]]),c('delta_slope')))
  # delta_slope[[j]]<-simulation_combined[[j]][,delta_cols]%>%as.matrix()%>%array(dim=c(n_sim,S))%>%
  #   apply(c(2),quantile,c(0.025,0.5,0.975))%>%melt(varnames=c('quantile','s'),value.name='y')%>%#mutate(quantile=c(0.025,0.5,0.975))%>%
  #   spread(quantile,y)%>%mutate(zeta_ind=zeta_ind[s],eta_ind=eta_ind[s],delta_ind=delta_ind[s],
  #                               zeta=zeta_value[s],eta=eta_value[s],delta=delta_value[s],
  #                               sim=j)
  # delta_slope[[j]]<- full_join(delta_slope[[j]],actual_delta, by='s')
  # 
  # 
  eta_cols<-which(str_detect(colnames(simulation_combined[[j]]),c('eta_slope')))
  eta_slope[[j]]<-simulation_combined[[j]][,eta_cols]%>%as.matrix()%>%array(dim=c(n_sim,S))%>%
    apply(c(2),quantile,c(0.025,0.5,0.975))%>%melt(varnames=c('quantile','s'),value.name='y')%>%#mutate(quantile=c(0.025,0.5,0.975))%>%
    spread(quantile,y)%>%mutate(zeta_ind=zeta_ind[s],eta_ind=eta_ind[s],delta_ind=delta_ind[s],
                                zeta=zeta_value[s],eta=eta_value[s],delta=delta_value[s],
                                sim=j)
  eta_slope[[j]]<- full_join(eta_slope[[j]],actual_eta, by='s')

  alpha_cols<-which(str_detect(colnames(simulation_combined[[j]]),c('alpha_spline')))
  alpha_spline[[j]]<-simulation_combined[[j]][,alpha_cols]%>%as.matrix()%>%array(dim=c(n_sim,N,S))

  
}
zeta_slope<-list()
for(j in 1:sims){
  zeta_slope[[j]]<-matrix(nrow=n_sim,ncol=S)    
  for(s in 1:S){
    for(i in 1:n_sim){
      zeta_slope[[j]][i,s]<-round(as.numeric(lm(data=tibble(y=alpha_spline[[j]][i,,s],x=((1:N)-mean(1:N))/sd(1:N)))$coefficients[2]),3)
    }
  }
}


zeta_slope_quantiles<-list()
for(j in 1:sims){
  zeta_slope_quantiles[[j]]<- zeta_slope[[j]]%>%
    apply(c(2),quantile,c(0.025,0.5,0.975))%>%melt(varnames=c('quantile','s'),value.name='y')%>%#mutate(quantile=c(0.025,0.5,0.975))%>%
    spread(quantile,y)%>%mutate(zeta_ind='Up',eta_ind=eta_ind[s],delta_ind=delta_ind[s],
                                zeta=zeta_value[s],eta=eta_value[s],delta=delta_value[s],
                                sim=j)
  zeta_slope_quantiles[[j]]<- full_join(zeta_slope_quantiles[[j]],actual_zeta, by='s')
  
}

zeta_slope_all_fixed<-unlist(zeta_slope_quantiles)%>%array(c(S,12,sims))
colnames(zeta_slope_all_fixed)<-c('s','2.5%','50%','97.5%','zeta_ind','eta_ind','delta_ind','zeta','eta','delta','sim','actual')
zeta_slope_all_fixed<-apply(zeta_slope_all_fixed,2,c)%>%as_tibble()%>%mutate(`2.5%`=as.numeric(`2.5%`),
                                                                             `50%`=as.numeric(`50%`),
                                                                             `97.5%`=as.numeric(`97.5%`),
                                                                             `actual`=as.numeric(`actual`))

# delta_slope_all_fixed<-unlist(delta_slope)%>%array(c(S,12,sims))
# colnames(delta_slope_all_fixed)<-c('s','2.5%','50%','97.5%','zeta_ind','eta_ind','delta_ind','zeta','eta','delta','sim','actual')
# delta_slope_all_fixed<-apply(delta_slope_all_fixed,2,c)%>%as_tibble()%>%mutate(`2.5%`=as.numeric(`2.5%`),
#                                                                                `50%`=as.numeric(`50%`),
#                                                                                `97.5%`=as.numeric(`97.5%`),
#                                                                                `actual`=as.numeric(`actual`))


eta_slope_all_fixed<-unlist(eta_slope)%>%array(c(S,12,sims))
colnames(eta_slope_all_fixed)<-c('s','2.5%','50%','97.5%','zeta_ind','eta_ind','delta_ind','zeta','eta','delta','sim','actual')
eta_slope_all_fixed<-apply(eta_slope_all_fixed,2,c)%>%as_tibble()%>%mutate(`2.5%`=as.numeric(`2.5%`),
                                                                           `50%`=as.numeric(`50%`),
                                                                           `97.5%`=as.numeric(`97.5%`),
                                                                           `actual`=as.numeric(`actual`))


save(y_quant_fixed,eta_slope_all_fixed,zeta_slope_all_fixed, file='matrix_simulations_prediction_fixed_delta.RData')

y_quant_fixed[[1]][2700,1:10]


# width=930, height=920
red <- "#E6194B"
green <- "#2AA198"
blue <- "#4363D8"
yellow <- "#FFE119"
gray <- "#A9A9A9"

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
  mutate(sim=as.factor(sim),zeta=factor(zeta,levels=unique(zeta)),
         eta=factor(eta,levels=unique(eta)),delta=factor(delta,levels=unique(delta)))






# Nowcasts from Survivor model.
n_sim<-dim(y_quant_fixed[[1]])[1]/S
y_fixed<-y_quant_fixed%>%unlist()%>%array(dim=c(dim(y_quant_fixed[[1]]),sims))
colnames(y_fixed)<- colnames(y_quant_fixed[[1]])
y_fixed<-apply(y_fixed,2,c)%>%as_tibble()%>%mutate(`2.5%`=as.numeric(`2.5%`),
                                                           `50%`=as.numeric(`50%`), 
                                                           `97.5%`=as.numeric(`97.5%`),
                                                           t=as.numeric(t),
                                                           s=as.numeric(s))%>%mutate(model='Survivor linear link (delta fixed)')

# Nowcast from Survivor linear model.
y_linear<-y_quant_linear_log[c(1:10)]%>%unlist()%>%array(dim=c(dim(y_quant_linear_log[[1]]),sims))
colnames(y_linear)<- colnames(y_quant_linear_log[[1]])
y_linear<-apply(y_linear,2,c)%>%as_tibble()%>%mutate(`2.5%`=as.numeric(`2.5%`),
                                                     `50%`=as.numeric(`50%`), 
                                                     `97.5%`=as.numeric(`97.5%`),t=as.numeric(t),
                                                     s=as.numeric(s))%>%mutate(model='Survivor linear link')

# Nowcasts from both models:
y_all<-rbind(y_fixed,y_linear)
# Add true y values:
y_all_sim<-right_join(sim_y_long,y_all, by=c('t','s','eta','zeta','delta','sim'))

ggplot()+geom_line(data=filter(y_all_sim,sim==5,delta==c(-0.2)),aes(x=t,y=y))+
  geom_line(data=filter(y_all_sim,sim==5,delta==c(-0.2)),aes(x=t,y=`50%`,colour=model))+
  facet_wrap(zeta~eta)+
  geom_ribbon(data=filter(y_all_sim,sim==5,delta==c(-0.2)),aes(x=t,ymin=`2.5%`,ymax=`97.5%`,fill=model),alpha=0.5)

# Filter for dates that are censored:
y_all_d<-y_all_sim%>%mutate(d=t-N)%>%filter(d>-D)

y_all_summary<-y_all_d%>%mutate(model=factor(model,levels=c('Survivor linear link (delta fixed)','Survivor linear link'
)),y=as.numeric(y),
d=as.numeric(d),sim=as.numeric(sim))

# compare bias for positive and negative delta
y_all_summary_metrics<-y_all_summary%>%group_by(model,d,s,delta)%>%
  summarise(Coverage=mean(y>=`2.5%`&y<=`97.5%`),
            `Prediction interval width`=mean(`97.5%`-`2.5%`),
            `Mean absolute error`=mean(abs(`50%`-y)),
            `Bias`=mean(`50%`-y))%>%
  pivot_longer(cols=c('Coverage','Prediction interval width','Mean absolute error','Bias'),names_to = 'metric')%>%
  mutate(model=factor(model,levels=c('Survivor linear link (delta fixed)','Survivor linear link'
  )),
  metric=factor(metric,levels=c('Mean absolute error','Prediction interval width','Coverage','Bias')))

y_all_summary_average<-y_all_summary_metrics%>%group_by(model,d,metric,delta)%>%
  summarise(value.mean=mean(as.numeric(value)), value.median=median(as.numeric(value)))

y_all_summary_average$value.mean[y_all_summary_average$metric%in%c("Prediction interval width","Mean absolute error","Bias")]<-NA
y_all_summary_average$value.median[y_all_summary_average$metric%in%c("Coverage")]<-NA
y_all_summary_average<-y_all_summary_average%>%group_by(model,d,metric,delta)%>%mutate(value=sum(value.mean,value.median,na.rm=TRUE))

# Plot of prediction performance metrics.
compare_delta_all_metric <- ggplot(y_all_summary_average)+
  geom_line(aes(x=d,y=value,colour=model),alpha=0.6)+
  geom_point(aes(x=d,y=value,colour=model,shape=model),alpha=0.6)+
  facet_grid2(delta~metric, scales = c("free"),independent = c("all"))+
  scale_colour_manual(name=NULL,values = c("#2AA198", "#D55E00"))+
  scale_shape_manual(name=NULL,values=c(15,3))+
  scale_x_continuous()+
  labs(x='Prediction time difference',y=NULL,
       title="Simulation experiment",
       subtitle='Performance metrics from the prediction experiment')+
  theme_minimal()+
  theme(legend.position = "bottom",text = element_text(size=12), plot.title = element_text(size=16),
        plot.subtitle=element_text(size=14), legend.title = element_text(size=14),
        strip.text = element_text(size=14))
compare_delta_all_metric



# Nowcasts from Survivor model.
n_sim<-dim(y_quant_fixed[[1]])[1]/S
y_fixed<-y_quant_fixed%>%unlist()%>%array(dim=c(dim(y_quant_fixed[[1]]),sims))
colnames(y_fixed)<- colnames(y_quant_fixed[[1]])
y_fixed<-apply(y_fixed,2,c)%>%as_tibble()%>%mutate(`2.5%`=as.numeric(`2.5%`),
                                                           `50%`=as.numeric(`50%`), 
                                                           `97.5%`=as.numeric(`97.5%`),
                                                           t=as.numeric(t),
                                                           s=as.numeric(s))%>%mutate(model='Survivor linear link (delta fixed)')
# Nowcast from Survivor linear model.
y_nonlinear<-y_quant[c(1:10)]%>%unlist()%>%array(dim=c(dim(y_quant[[1]]),sims))
colnames(y_nonlinear)<- colnames(y_quant[[1]])
y_nonlinear<-apply(y_nonlinear,2,c)%>%as_tibble()%>%mutate(`2.5%`=as.numeric(`2.5%`),
                                                           `50%`=as.numeric(`50%`), 
                                                           `97.5%`=as.numeric(`97.5%`),t=as.numeric(t),
                                                           s=as.numeric(s))%>%mutate(model='Survivor')

# Nowcasts from both models:
y_all<-rbind(y_fixed,y_nonlinear)
# Add true y values:
y_all_sim<-right_join(sim_y_long,y_all, by=c('t','s','eta','zeta','delta','sim'))

ggplot()+geom_line(data=filter(y_all_sim,sim==10,delta==c(-0.2)),aes(x=t,y=y))+
  geom_line(data=filter(y_all_sim,sim==10,delta==c(-0.2)),aes(x=t,y=`50%`,colour=model))+
  facet_grid2(paste('Delay temporal trend:',eta)~paste('Totals temporal trend:',zeta))+
  geom_ribbon(data=filter(y_all_sim,sim==10,delta==c(-0.2)),aes(x=t,ymin=`2.5%`,ymax=`97.5%`,fill=model),alpha=0.25)

ggplot()+geom_line(data=filter(y_all_sim,sim==10,delta==c(0.2)),aes(x=t,y=y))+
  geom_line(data=filter(y_all_sim,sim==10,delta==c(0.2)),aes(x=t,y=`50%`,colour=model))+
  facet_grid2(paste('Delay temporal trend:',eta)~paste('Totals temporal trend:',zeta))+
  geom_ribbon(data=filter(y_all_sim,sim==10,delta==c(0.2)),aes(x=t,ymin=`2.5%`,ymax=`97.5%`,fill=model),alpha=0.25)

# Filter for dates that are censored:
y_all_d<-y_all_sim%>%mutate(d=t-N)%>%filter(d>-D)

y_all_summary<-y_all_d%>%mutate(model=factor(model,levels=c('Survivor linear link (delta fixed)','Survivor'
)),y=as.numeric(y),
d=as.numeric(d),sim=as.numeric(sim))
# compare bias for positive and negative delta
y_all_summary_metrics<-y_all_summary%>%group_by(model,d,s,delta)%>%
  summarise(Coverage=mean(y>=`2.5%`&y<=`97.5%`),
            `Prediction interval width`=mean(`97.5%`-`2.5%`),
            `Mean absolute error`=mean(abs(`50%`-y)),
            `Bias`=mean(`50%`-y))%>%
  pivot_longer(cols=c('Coverage','Prediction interval width','Mean absolute error','Bias'),names_to = 'metric')%>%
  mutate(model=factor(model,levels=c('Survivor linear link (delta fixed)','Survivor'
  )),
  metric=factor(metric,levels=c('Mean absolute error','Prediction interval width','Coverage','Bias')))


y_all_summary_average_no_zero<-filter(y_all_summary_metrics,delta!='0')
y_all_summary_average_no_zero<-y_all_summary_metrics%>%group_by(model,d,metric)%>%
  summarise(value.mean=mean(as.numeric(value)), value.median=median(as.numeric(value)))

y_all_summary_average_no_zero$value.mean[y_all_summary_average_no_zero$metric%in%c("Prediction Interval Width","Mean Absolute Error")]<-NA
y_all_summary_average_no_zero$value.median[y_all_summary_average_no_zero$metric%in%c("Coverage")]<-NA
y_all_summary_average_no_zero<-y_all_summary_average_no_zero%>%group_by(model,d,metric)%>%mutate(value=sum(value.mean,value.median,na.rm=TRUE))
# Plot of prediction performance metrics.
compare_survivor_non_zero <- ggplot(y_all_summary_average_no_zero)+
  geom_line(aes(x=d,y=value,colour=model),alpha=0.6)+
  geom_point(aes(x=d,y=value,colour=model,shape=model),alpha=0.6)+
  facet_grid2(~metric, scales = c("free"),independent = c("all"))+
  scale_colour_manual(name=NULL,values = c( "#2AA198", "#4363D8"))+
  scale_shape_manual(name=NULL,values=c(15,3))+
  scale_x_continuous()+
  labs(x='Prediction time difference',y=NULL,
       title="Simulation experiment",
       subtitle='Performance metrics from the prediction experiment')+
  theme_minimal()+
  theme(legend.position = "bottom",text = element_text(size=12), plot.title = element_text(size=16),
        plot.subtitle=element_text(size=14), legend.title = element_text(size=14),
        strip.text = element_text(size=14))
compare_survivor_non_zero
ggsave(compare_survivor_non_zero, file="Plots/Survivor_fixed_log_compare.pdf", width=12, height=4)
# 

y_all_summary_average<-y_all_summary_metrics%>%group_by(model,d,metric,delta)%>%
  summarise(value.mean=mean(as.numeric(value)), value.median=median(as.numeric(value)))

y_all_summary_average$value.mean[y_all_summary_average$metric%in%c("Prediction interval width","Mean absolute error","Bias")]<-NA
y_all_summary_average$value.median[y_all_summary_average$metric%in%c("Coverage")]<-NA
y_all_summary_average<-y_all_summary_average%>%group_by(model,d,metric,delta)%>%mutate(value=sum(value.mean,value.median,na.rm=TRUE))
# Plot of prediction performance metrics.
compare_delta_all_metric <- ggplot(y_all_summary_average)+
  geom_line(aes(x=d,y=value,colour=model),alpha=0.6)+
  geom_point(aes(x=d,y=value,colour=model,shape=model),alpha=0.6)+
  facet_grid2(delta~metric, scales = c("free"),independent = c("all"))+
  scale_colour_manual(name=NULL,values = c("#2AA198","#4363D8"))+
  scale_shape_manual(name=NULL,values=c(15,3))+
  scale_x_continuous()+
  labs(x='Prediction time difference',y=NULL,
       title="Simulation experiment",
       subtitle='Performance metrics from the prediction experiment')+
  theme_minimal()+
  theme(legend.position = "bottom",text = element_text(size=12), plot.title = element_text(size=16),
        plot.subtitle=element_text(size=14), legend.title = element_text(size=14),
        strip.text = element_text(size=14))
compare_delta_all_metric


compare_delta_mae <- ggplot(filter(y_all_summary_average,metric=="Mean absolute error"))+
  geom_line(aes(x=d,y=value,colour=model),alpha=0.6)+
  geom_point(aes(x=d,y=value,colour=model,shape=model),alpha=0.6)+
  facet_grid2(~paste('Case load effect:',delta),switch=c('y'))+
  scale_colour_manual(name=NULL,values = c("#2AA198", "#4363D8"))+
  scale_shape_manual(name=NULL,values=c(15,3))+
  scale_x_continuous()+
  labs(x=NULL,y="Mean absolute error",
       title="Simulation experiment",
       subtitle='Performance metrics from the prediction experiment')+
  theme_minimal()+guides(colour="none", shape='none')+
  theme(legend.position = "bottom",text = element_text(size=14), axis.title = element_text(size=16), plot.title = element_text(size=20),plot.subtitle=element_text(size=18),
        legend.title = element_text(size=16),strip.text = element_text(size=14))

compare_delta_piw <- ggplot(filter(y_all_summary_average,metric=="Prediction interval width"))+
  geom_line(aes(x=d,y=value,colour=model),alpha=0.6)+
  geom_point(aes(x=d,y=value,colour=model,shape=model),alpha=0.6)+
  facet_grid2(~paste('Case load effect:',delta))+
  scale_colour_manual(name=NULL,values =c("#2AA198", "#4363D8"))+
  scale_shape_manual(name=NULL,values=c(15,3))+
  scale_x_continuous()+
  labs(x=NULL,y="Prediction interval width",
       title=NULL,
       subtitle=NULL)+
  theme_minimal()+guides(colour="none", shape='none')+
  theme(legend.position = "bottom",text = element_text(size=14), axis.title = element_text(size=16), plot.title = element_text(size=20),plot.subtitle=element_text(size=18),
        legend.title = element_text(size=16),strip.text = element_text(size=14))

compare_delta_cov <- ggplot(filter(y_all_summary_average,metric=="Coverage"))+
  geom_line(aes(x=d,y=value,colour=model),alpha=0.6)+
  geom_point(aes(x=d,y=value,colour=model,shape=model),alpha=0.6)+
  facet_grid2(~paste('Case load effect:',delta))+
  scale_colour_manual(name=NULL,values = c("#2AA198", "#4363D8"))+
  scale_shape_manual(name=NULL,values=c(15,3))+
  scale_x_continuous()+
  labs(x='Prediction time difference',y="Coverage",
       title=NULL,
       subtitle=NULL)+
  theme_minimal()+ 
  theme(legend.position = "bottom",text = element_text(size=14), axis.title = element_text(size=16), plot.title = element_text(size=20),plot.subtitle=element_text(size=18),
        legend.title = element_text(size=16),strip.text = element_text(size=14))

compare_delta_sep<-grid.arrange(compare_delta_mae,compare_delta_piw,compare_delta_cov,nrow=3,heights=c(4.1,3.2,4.1))


ggsave(compare_delta_sep, file="Plots/fixed_delta_pred_sep.pdf", width=9, height=11)



y_all_summary_metrics_neg<-filter(y_all_summary,delta==(-0.2))%>%group_by(model,d,s,eta,zeta)%>%
  summarise(Coverage=mean(y>=`2.5%`&y<=`97.5%`),
            `Prediction interval width`=mean(`97.5%`-`2.5%`),
            `Mean absolute error`=mean(abs(`50%`-y)),
            `Bias`=mean(`50%`-y))%>%
  pivot_longer(cols=c('Coverage','Prediction interval width','Mean absolute error','Bias'),names_to = 'metric')%>%
  mutate(model=factor(model,levels=c('Survivor linear link (delta fixed)','Survivor'
  )),
  metric=factor(metric,levels=c('Mean absolute error','Prediction interval width','Coverage','Bias')))

y_all_summary_average_neg<-y_all_summary_metrics_neg%>%group_by(model,d,metric,eta,zeta)%>%
  summarise(value.mean=mean(as.numeric(value)), value.median=median(as.numeric(value)))

y_all_summary_average_neg$value.mean[y_all_summary_average_neg$metric%in%c("Prediction interval width","Mean absolute error","Bias")]<-NA
y_all_summary_average_neg$value.median[y_all_summary_average_neg$metric%in%c("Coverage")]<-NA
y_all_summary_average_neg<-y_all_summary_average_neg%>%group_by(model,d,metric,eta,zeta)%>%mutate(value=sum(value.mean,value.median,na.rm=TRUE))
compare_delta_mae <- ggplot(filter(y_all_summary_average_neg,metric=="Mean absolute error"))+
  geom_line(aes(x=d,y=value,colour=model),alpha=0.6)+
  geom_point(aes(x=d,y=value,colour=model,shape=model),alpha=0.6)+
  facet_grid2(paste('Delay temporal trend:',eta)~paste('Totals temporal trend:',zeta), scales = c("free"),independent = c("all"))+
  scale_colour_manual(name=NULL,values = c( "#2AA198","#4363D8"))+
  scale_shape_manual(name=NULL,values=c(15,3))+
  scale_x_continuous()+
  labs(x='Prediction time difference',y=NULL,
       title="Simulation experiment",
       subtitle='Mean absolute error: delta=-0.2')+
  theme_minimal()+
  theme(legend.position = "bottom",text = element_text(size=12), plot.title = element_text(size=16),
        plot.subtitle=element_text(size=14), legend.title = element_text(size=14),
        strip.text = element_text(size=14))
compare_delta_mae

compare_delta_bias <- ggplot(filter(y_all_summary_average_neg,metric=="Bias"))+
  geom_line(aes(x=d,y=value,colour=model),alpha=0.6)+
  geom_point(aes(x=d,y=value,colour=model,shape=model),alpha=0.6)+
  facet_grid2(paste('Delay temporal trend:',eta)~paste('Totals temporal trend:',zeta), scales = c("free"),independent = c("all"))+
  scale_colour_manual(name=NULL,values = c( "#2AA198", "#4363D8"))+
  scale_shape_manual(name=NULL,values=c(15,3))+
  scale_x_continuous()+
  labs(x='Prediction time difference',y=NULL,
       title="Simulation experiment",
       subtitle='Bias: delta=-0.2')+
  theme_minimal()+
  theme(legend.position = "bottom",text = element_text(size=12), plot.title = element_text(size=16),
        plot.subtitle=element_text(size=14), legend.title = element_text(size=14),
        strip.text = element_text(size=14))
compare_delta_bias

compare_delta_piw <- ggplot(filter(y_all_summary_average_neg,metric=="Prediction interval width"))+
  geom_line(aes(x=d,y=value,colour=model),alpha=0.6)+
  geom_point(aes(x=d,y=value,colour=model,shape=model),alpha=0.6)+
  facet_grid2(paste('Delay temporal trend:',eta)~paste('Totals temporal trend:',zeta), scales = c("free"),independent = c("all"))+
  scale_colour_manual(name=NULL,values = c( "#2AA198", "#4363D8"))+
  scale_shape_manual(name=NULL,values=c(15,3))+
  scale_x_continuous()+
  labs(x='Prediction time difference',y=NULL,
       title="Simulation experiment",
       subtitle='Prediction interval width: delta=-0.2')+
  theme_minimal()+
  theme(legend.position = "bottom",text = element_text(size=12), plot.title = element_text(size=16),
        plot.subtitle=element_text(size=14), legend.title = element_text(size=14),
        strip.text = element_text(size=14))
compare_delta_piw
ggsave(compare_delta_piw, file="Plots/FIXED_neg_compare_delta_piw.pdf", width=9, height=10)

y_all_summary_metrics_pos<-filter(y_all_summary,delta==(0.2))%>%group_by(model,d,s,eta,zeta)%>%
  summarise(Coverage=mean(y>=`2.5%`&y<=`97.5%`),
            `Prediction interval width`=mean(`97.5%`-`2.5%`),
            `Mean absolute error`=mean(abs(`50%`-y)),
            `Bias`=mean(`50%`-y))%>%
  pivot_longer(cols=c('Coverage','Prediction interval width','Mean absolute error','Bias'),names_to = 'metric')%>%
  mutate(model=factor(model,levels=c('Survivor linear link (delta fixed)','Survivor'
  )),
  metric=factor(metric,levels=c('Mean absolute error','Prediction interval width','Coverage','Bias')))

y_all_summary_average_pos<-y_all_summary_metrics_pos%>%group_by(model,d,metric,eta,zeta)%>%
  summarise(value.mean=mean(as.numeric(value)), value.median=median(as.numeric(value)))

y_all_summary_average_pos$value.mean[y_all_summary_average_pos$metric%in%c("Prediction interval width","Mean absolute error","Bias")]<-NA
y_all_summary_average_pos$value.median[y_all_summary_average_pos$metric%in%c("Coverage")]<-NA
y_all_summary_average_pos<-y_all_summary_average_pos%>%group_by(model,d,metric,eta,zeta)%>%mutate(value=sum(value.mean,value.median,na.rm=TRUE))
compare_delta_mae_pos <- ggplot(filter(y_all_summary_average_pos,metric=="Mean absolute error"))+
  geom_line(aes(x=d,y=value,colour=model),alpha=0.6)+
  geom_point(aes(x=d,y=value,colour=model,shape=model),alpha=0.6)+
  facet_grid2(paste('Delay temporal trend:',eta)~paste('Totals temporal trend:',zeta), scales = c("free"),independent = c("all"))+
  scale_colour_manual(name=NULL,values = c( "#2AA198", "#4363D8"))+
  scale_shape_manual(name=NULL,values=c(15,3))+
  scale_x_continuous()+
  labs(x='Prediction time difference',y=NULL,
       title="Simulation experiment",
       subtitle='Mean absolute error: delta= 0.2')+
  theme_minimal()+
  theme(legend.position = "bottom",text = element_text(size=12), plot.title = element_text(size=16),
        plot.subtitle=element_text(size=14), legend.title = element_text(size=14),
        strip.text = element_text(size=14))
compare_delta_mae_pos

compare_delta_bias_pos <- ggplot(filter(y_all_summary_average_pos,metric=="Bias"))+
  geom_line(aes(x=d,y=value,colour=model),alpha=0.6)+
  geom_point(aes(x=d,y=value,colour=model,shape=model),alpha=0.6)+
  facet_grid2(paste('Delay temporal trend:',eta)~paste('Totals temporal trend:',zeta), scales = c("free"),independent = c("all"))+
  scale_colour_manual(name=NULL,values = c("#2AA198", "#4363D8"))+
  scale_shape_manual(name=NULL,values=c(15,3))+
  scale_x_continuous()+
  labs(x='Prediction time difference',y=NULL,
       title="Simulation experiment",
       subtitle='Bias: delta= 0.2')+
  theme_minimal()+
  theme(legend.position = "bottom",text = element_text(size=12), plot.title = element_text(size=16),
        plot.subtitle=element_text(size=14), legend.title = element_text(size=14),
        strip.text = element_text(size=14))
compare_delta_bias_pos


compare_delta_piw_pos <- ggplot(filter(y_all_summary_average_pos,metric=="Prediction interval width"))+
  geom_line(aes(x=d,y=value,colour=model),alpha=0.6)+
  geom_point(aes(x=d,y=value,colour=model,shape=model),alpha=0.6)+
  facet_grid2(paste('Delay temporal trend:',eta)~paste('Totals temporal trend:',zeta), scales = c("free"),independent = c("all"))+
  scale_colour_manual(name=NULL,values = c("#2AA198", "#4363D8"))+
  scale_shape_manual(name=NULL,values=c(15,3))+
  scale_x_continuous()+
  labs(x='Prediction time difference',y=NULL,
       title="Simulation experiment",
       subtitle='Prediction interval width: delta= 0.2')+
  theme_minimal()+
  theme(legend.position = "bottom",text = element_text(size=12), plot.title = element_text(size=16),
        plot.subtitle=element_text(size=14), legend.title = element_text(size=14),
        strip.text = element_text(size=14))
compare_delta_piw_pos
ggsave(compare_delta_piw_pos, file="Plots/FIXED_pos_compare_delta_piw.pdf", width=9, height=10)



# 
# # Nowcasts from Survivor model.
# n_sim<-dim(y_quant_fixed[[1]])[1]/S
# y_fixed<-y_quant_fixed%>%unlist()%>%array(dim=c(dim(y_quant_fixed[[1]]),sims))
# colnames(y_fixed)<- colnames(y_quant_fixed[[1]])
# y_fixed<-apply(y_fixed,2,c)%>%as_tibble()%>%mutate(`2.5%`=as.numeric(`2.5%`),
#                                                            `50%`=as.numeric(`50%`), 
#                                                            `97.5%`=as.numeric(`97.5%`),
#                                                            t=as.numeric(t),
#                                                            s=as.numeric(s))%>%mutate(model='Survivor linear link (delta fixed)')
# # Nowcast from Survivor linear model.
# y_fixed<-y_quant_fixed[c(1:10)]%>%unlist()%>%array(dim=c(dim(y_quant_fixed[[1]]),sims))
# colnames(y_fixed)<- colnames(y_quant_fixed[[1]])
# y_fixed<-apply(y_fixed,2,c)%>%as_tibble()%>%mutate(`2.5%`=as.numeric(`2.5%`),
#                                                      `50%`=as.numeric(`50%`), 
#                                                      `97.5%`=as.numeric(`97.5%`),t=as.numeric(t),
#                                                      s=as.numeric(s))%>%mutate(model='Survivor linear link log fixed')
# # 
# # # Nowcasts from both models:
# y_all<-rbind(y_fixed,y_fixed)
# # Add true y values:
# y_all_sim<-right_join(sim_y_long,y_all, by=c('t','s','eta','zeta','delta','sim'))
# 
# ggplot()+geom_line(data=filter(y_all_sim,sim==5,delta==c(-0.2)),aes(x=t,y=y))+
#   geom_line(data=filter(y_all_sim,sim==5,delta==c(-0.2)),aes(x=t,y=`50%`,colour=model))+
#   facet_wrap(zeta~eta)+
#   geom_ribbon(data=filter(y_all_sim,sim==5,delta==c(-0.2)),aes(x=t,ymin=`2.5%`,ymax=`97.5%`,fill=model),alpha=0.5)
# 
# # Filter for dates that are censored:
# y_all_d<-y_all_sim%>%mutate(d=t-N)%>%filter(d>-D)
# 
# y_all_summary<-y_all_d%>%mutate(model=factor(model,levels=c('Survivor linear link (delta fixed)','Survivor linear link log fixed'
# )),y=as.numeric(y),
# d=as.numeric(d),sim=as.numeric(sim))
# # compare bias for positive and negative delta
# y_all_summary_metrics<-y_all_summary%>%group_by(model,d,s,delta)%>%
#   summarise(Coverage=mean(y>=`2.5%`&y<=`97.5%`),
#             `Prediction interval width`=mean(`97.5%`-`2.5%`),
#             `Mean absolute error`=mean(abs(`50%`-y)),
#             `Bias`=mean(`50%`-y))%>%
#   pivot_longer(cols=c('Coverage','Prediction interval width','Mean absolute error','Bias'),names_to = 'metric')%>%
#   mutate(model=factor(model,levels=c('Survivor linear link (delta fixed)','Survivor linear link log fixed'
#   )),
#   metric=factor(metric,levels=c('Mean absolute error','Prediction interval width','Coverage','Bias')))
# 
# y_all_summary_average<-y_all_summary_metrics%>%group_by(model,d,metric,delta)%>%
#   summarise(value.mean=mean(as.numeric(value)), value.median=median(as.numeric(value)))
# 
# y_all_summary_average$value.mean[y_all_summary_average$metric%in%c("Prediction interval width","Mean absolute error","Bias")]<-NA
# y_all_summary_average$value.median[y_all_summary_average$metric%in%c("Coverage")]<-NA
# y_all_summary_average<-y_all_summary_average%>%group_by(model,d,metric,delta)%>%mutate(value=sum(value.mean,value.median,na.rm=TRUE))
# # Plot of prediction performance metrics.
# compare_delta_all_metric <- ggplot(y_all_summary_average)+
#   geom_line(aes(x=d,y=value,colour=model),alpha=0.6)+
#   geom_point(aes(x=d,y=value,colour=model,shape=model),alpha=0.6)+
#   facet_grid2(delta~metric, scales = c("free"),independent = c("all"))+
#   scale_colour_manual(name=NULL,values = c("#E6194B", "#D55E00"))+
#   scale_shape_manual(name=NULL,values=c(15,3))+
#   scale_x_continuous()+
#   labs(x='Prediction time difference',y=NULL,
#        title="Simulation experiment",
#        subtitle='Performance metrics from the prediction experiment')+
#   theme_minimal()+
#   theme(legend.position = "bottom",text = element_text(size=12), plot.title = element_text(size=16),
#         plot.subtitle=element_text(size=14), legend.title = element_text(size=14),
#         strip.text = element_text(size=14))
# compare_delta_all_metric
# 
# 

# plot data vs estimates 
# total trends parameter
ggplot(data=zeta_slope_all_fixed)+
  geom_errorbar(aes(x=delta,  ymin=`2.5%`,ymax=`97.5%`,colour=as.factor(eta)),alpha=0.5)+
  facet_wrap(~paste('Totals temporal trend:',zeta))+
  geom_point(aes(x=delta, y=`50%`,colour=as.factor(eta)),shape=3,size=3)+
  geom_point(aes(x=delta, y=actual),size=3)+theme_minimal()+
  theme(text = element_text(size=12), plot.title = element_text(size=16),plot.subtitle=element_text(size=14),
        legend.title = element_text(size=14),strip.text = element_text(size=14),legend.position="bottom")+
  scale_colour_viridis_d(name='Temporal trend of reporting (eta)',begin=0.1,end=0.9)+
  scale_shape(name='Effect of case load on delay (delta)')+
  labs(title='Fixed delta model estimates against simulated parameter values',
       subtitle = 'For linear temporal trend in total counts (zeta)',
       x='Effect of case load on delay (delta)',y='Value for temporal trend in total cases (zeta)')

ggplot(data=zeta_slope_all_linear_log)+
  geom_errorbar(aes(x=delta,  ymin=`2.5%`,ymax=`97.5%`,colour=as.factor(eta)),alpha=0.5)+
  facet_wrap(~paste('Totals temporal trend:',zeta))+
  geom_point(aes(x=delta, y=`50%`,colour=as.factor(eta)),shape=3,size=3)+
  geom_point(aes(x=delta, y=actual),size=3)+theme_minimal()+
  theme(text = element_text(size=12), plot.title = element_text(size=16),plot.subtitle=element_text(size=14),
        legend.title = element_text(size=14),strip.text = element_text(size=14),legend.position="bottom")+
  scale_colour_viridis_d(name='Temporal trend of reporting (eta)',begin=0.1,end=0.9)+
  scale_shape(name='Effect of case load on delay (delta)')+
  labs(title='Linear link model estimates against simulated parameter values',
       subtitle = 'For linear temporal trend in total counts (zeta)',
       x='Effect of case load on delay (delta)',y='Value for temporal trend in total cases (zeta)')

# reporting trend parameter
ggplot(data=eta_slope_all_fixed)+
  geom_errorbar(aes(x=eta,  ymin=`2.5%`,ymax=`97.5%`,colour=as.factor(delta)),alpha=0.5)+
  facet_wrap(~paste('Totals temporal trend:',zeta))+
  geom_point(aes(x=eta, y=`50%`,colour=as.factor(delta)),shape=3,size=3)+
  geom_point(aes(x=eta, y=actual),size=3)+theme_minimal()+
  theme(text = element_text(size=12), plot.title = element_text(size=16),plot.subtitle=element_text(size=14),
        legend.title = element_text(size=14),strip.text = element_text(size=14),legend.position="bottom")+
  scale_colour_viridis_d(name='Effect of case load on delay (delta)',begin=0.1,end=0.9)+
  scale_shape(name='Effect of case load on delay (delta)')+
  labs(title='Fixed delta model estimates against simulated parameter values',
       subtitle = 'For linear temporal trend in reporting delays (eta)',
       x='Temporal trend of reporting (eta)',y='Temporal trend of reporting (eta)')


ggplot(data=eta_slope_all_linear_log)+
  geom_errorbar(aes(x=eta,  ymin=`2.5%`,ymax=`97.5%`,colour=as.factor(delta)),alpha=0.5)+
  facet_wrap(~paste('Totals temporal trend:',zeta))+
  geom_point(aes(x=eta, y=`50%`,colour=as.factor(delta)),shape=3,size=3)+
  geom_point(aes(x=eta, y=actual),size=3)+theme_minimal()+
  theme(text = element_text(size=12), plot.title = element_text(size=16),plot.subtitle=element_text(size=14),
        legend.title = element_text(size=14),strip.text = element_text(size=14),legend.position="bottom")+
  scale_colour_viridis_d(name='Effect of case load on delay (delta)',begin=0.1,end=0.9)+
  scale_shape(name='Effect of case load on delay (delta)')+
  labs(title='Linear link model estimates against simulated parameter values',
       subtitle = 'For linear temporal trend in reporting delays (eta)',
       x='Temporal trend of reporting (eta)',y='Temporal trend of reporting (eta)')

