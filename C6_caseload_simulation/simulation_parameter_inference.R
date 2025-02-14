################################################
#  Simulation Experiment: PARAMETER INFERENCE  #
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


##### FIT GDM TO SIMULATED DATA #####

# Create function for cluster to carry out MCMC model:
Simulation_survivor<- function(seed, sim_data, S, N, D, n_iter, n_burn, n_thin, n_chains){
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
        for(t in 1:N){
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
      delta_slope[s] ~ dnorm(0,sd=10) 
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
  sim_y<-sim_data[[seed]]$y
  sim_z<-sim_data[[seed]]$z
  # Constants from simulated data:
  N<-dim(sim_y)[1]
  S<-dim(sim_y)[2]
  D<-dim(sim_z)[3]
  
  # Model Data
  simulation_data <- list(y=sim_y[1:N,1:S], z=sim_z[1:N,1:S,1:D])
  
  # Constants
  simulation_constants <-list(N=N,  D=D, S=S, 
                              K_t=dim(blank_jagam$jags.data$S1)[1])
  simulation_constants$S_t <- blank_jagam$jags.data$S1
  simulation_constants$X_t <- blank_jagam$jags.data$X[,2:(simulation_constants$K_t+1)]
  simulation_constants$zeros <- rep(0,max(simulation_constants$K_t))
  
  #simulation_inits<-simulation_model<-simulation_compiled_model<- simulation_mcmc_config<-simulation_mcmc <- simulation_compiled_mcmc<-list()
  # Generate random initial values.
  simulation_inits <- list(
    delta_slope=rnorm(S,0,sd=0.1),
    kappa_alpha=matrix(rnorm(S*(simulation_constants$K_t),0,0.01),ncol=S),
    sigma_alpha=matrix(runif(S*2,0,1),ncol=2),
    # zeta_slope=rnorm(S,0,sd=0.1),
    #kappa_zeta=matrix(rnorm(S*(simulation_constants$K_t),0,0.01),ncol=S),
    eta_slope=rnorm(S,0,sd=0.1),
    psi=matrix(rep(sort(rnorm(D,0,1)),S),ncol=S, byrow = FALSE),
    theta=abs(rnorm(S,50,sd=10)),
    iota=rnorm(S,5,1),
    phi=matrix(abs(rnorm(D*S,50,sd=10)),ncol=S))
  
  # Build the model.
  simulation_model <- nimbleModel(simulation_code,simulation_constants,simulation_data,simulation_inits)
  # Compile the model.
  simulation_compiled_model <- compileNimble(simulation_model)
  
  # Set up the MCMC.
  simulation_mcmc_config <- configureMCMC(simulation_model,monitors=c("lambda","iota","alpha_spline",
                                                                      "psi","theta","delta_slope",
                                                                      "eta_slope","kappa_alpha","sigma_alpha"
  ),
  useConjugacy = FALSE)
  
  simulation_mcmc_config$removeSamplers('eta_slope','psi','iota','kappa_alpha','sigma_alpha','delta_slope')#
  
  for(s in 1:S){
    simulation_mcmc_config$ (target=c(paste('delta_slope[',s,']', sep=''), 
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
n_iter=20000
# Set burn-in.
n_burn=15000
# Set thinning.
n_thin=5
# Set number of knots for model splines
#n_knots<-# knots for time, time-delay interaction and delay respectively 
# Load simulated data sets:
#load("~/media/alba/Disk 1/OneDrive/PhD/Case load simulations/Simulation studies/Simulation_data_y.RData")
# Number of simulations.
sims=100
# Set maximum COVID-19 lab report delay.
n_cores<-min(detectCores(),sims,25)

# Make Cluster for MCMC. 
this_cluster_simulation<- makeCluster(n_cores)

# Run Cluster.
time_surv_sim_inference <- system.time({
  simulation_surv_sim_inference <-  parLapply(cl = this_cluster_simulation, X = 1:sims, 
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
stopCluster(this_cluster_simulation)

# Save survivor model output. 
save(simulation_surv_sim_inference,time_surv_sim_inference, file='simulation_parameter_inference.RData')

# Check PSRF of model parameters.
parameter_group=c("lambda","iota","alpha_spline",
                  "psi","theta","delta_slope",
                  "eta_slope","kappa_alpha","sigma_alpha")
sari_psrf<-list()
sari_mean_psrf<-sari_max_psrf<-sari_leq_105 <-sari_leq_110<-array(NA,dim=c(length(parameter_group),sims))
for(k in 1:length(parameter_group)){
  sari_psrf[[k]]<-list()  
  for(j in 1:sims){
    sari_parameter_names <- colnames((as.matrix(simulation_surv_sim_inference[[j]])))
    sari_index <- which(startsWith(sari_parameter_names,parameter_group[k]))
    sari_psrf[[k]][[j]]<- sapply(sari_index,function(v)gelman.diag(simulation_surv_sim_inference[[j]][,v],autoburnin = FALSE,transform = TRUE,multivariate = FALSE)$psrf[,1])
    sari_mean_psrf[k,j] <- mean(sari_psrf[[k]][[j]])
    sari_max_psrf[k,j] <- max(sari_psrf[[k]][[j]])
    sari_leq_105[k,j] <- mean(sari_psrf[[k]][[j]]<=1.05)
    sari_leq_110[k,j] <- mean(sari_psrf[[k]][[j]][j]<=1.1)
  }
}
sari_MAX_psrf<-cbind(parameter_group,apply(sari_max_psrf,1,max))
sari_MEAN_psrf<-cbind(parameter_group,apply(sari_mean_psrf,1,mean))



plot(simulation_surv_sim_inference[[1]][,"delta_slope[1]"],main="delta_slope")
plot(simulation_surv_sim_inference[[2]][,"delta_slope[2]"],main="delta_slope")
plot(simulation_surv_sim_inference[[1]][,"delta_slope[3]"],main="delta_slope")
plot(simulation_surv_sim_inference[[1]][,"delta_slope[4]"],main="delta_slope")
plot(simulation_surv_sim_inference[[2]][,"delta_slope[5]"],main="delta_slope")
plot(simulation_surv_sim_inference[[1]][,"delta_slope[6]"],main="delta_slope")
plot(simulation_surv_sim_inference[[2]][,"delta_slope[7]"],main="delta_slope")
plot(simulation_surv_sim_inference[[1]][,"delta_slope[8]"],main="delta_slope")
plot(simulation_surv_sim_inference[[2]][,"delta_slope[9]"],main="delta_slope")
# plot(simulation_surv_sim_inference[[3]][,"delta_slope[7]"],main="delta_slope")
# plot(simulation_surv_sim_inference[[5]][,"delta_slope[7]"],main="delta_slope")
# plot(simulation_surv_sim_inference[[6]][,"delta_slope[7]"],main="delta_slope")
# plot(simulation_surv_sim_inference[[7]][,"delta_slope[7]"],main="delta_slope")
# plot(simulation_surv_sim_inference[[9]][,"delta_slope[7]"],main="delta_slope")
# plot(simulation_surv_sim_inference[[10]][,"delta_slope[7]"],main="delta_slope")
# 
# plot(simulation_surv_sim_inference[[8]][,"delta_slope[7]"],main="delta_slope")
# plot(simulation_surv_sim_inference[[4]][,"delta_slope[7]"],main="delta_slope")
# 
# plot(simulation_surv_sim_inference[[1]][,"alpha[1, 1]"],main="zeta[1] 1")
# plot(simulation_surv_sim_inference[[2]][,"zeta[50]"],main="zeta[50] 2")
# plot(simulation_surv_sim_inference[[3]][,"zeta[38]"],main="zeta[38] 3")
# 
# plot(simulation_surv_sim_inference[[1]][,"eta_slope[1]"],main="eta_slope 1")
# plot(simulation_surv_sim_inference[[2]][,"eta_slope"],main="eta_slope 2")
# plot(simulation_surv_sim_inference[[3]][,"eta_slope"],main="eta_slope 3")

simulation_mcmc<-list()
simulation_combined<-list()
for(j in 1:sims){
  simulation_mcmc[[j]]<-as.mcmc.list(simulation_surv_sim_inference[[j]])
  # Combine all MCMC chains.
  simulation_combined[[j]] <- as_tibble(do.call('rbind',simulation_mcmc[[j]]))
}

# Model coefficient output
n_sim<-dim(simulation_combined[[1]])[1]
delta_slope<-zeta_slope<-eta_slope<-list()

# Values from simulated data
eta_value<-simulation_data_list[[1]]$eta
delta_value<-simulation_data_list[[1]]$delta
zeta_value<-simulation_data_list[[1]]$zeta
# What each value translates to:
# Temporal trend in reporting is...
eta_ind<-rep(c('Improving','Not changing','Worsening'),S/3)
# Case load effect is...
delta_ind<-rep(c(rep('Negative effect',S/9),rep('No effect',S/9),rep('Positive effect',S/9)),S/9)
# Temporal trend in total cases is...
zeta_ind<-c(rep('Downward',S/3),rep('Stable',S/3),rep('Upward',S/3))

# Add simulated parameter values to table:
actual_delta <- tibble(delta_slopes=as.numeric(delta_value),s=(1:S))
actual_zeta <- tibble(zeta_slopes=as.numeric(zeta_value),s=1:S)
actual_eta <-tibble(eta_slopes=as.numeric(eta_value),s=(1:S))
alpha<-kappa_alpha<-lambda<-alpha_spline<-list()


for(j in 1:sims){
  delta_cols<-which(str_detect(colnames(simulation_combined[[j]]),c('delta_slope')))
  delta_slope[[j]]<-simulation_combined[[j]][,delta_cols]%>%as.matrix()%>%array(dim=c(n_sim,S))%>%
    apply(c(2),quantile,c(0.025,0.5,0.975))%>%melt(varnames=c('quantile','s'),value.name='y')%>%#mutate(quantile=c(0.025,0.5,0.975))%>%
    spread(quantile,y)%>%mutate(zeta_ind=zeta_ind[s],eta_ind=eta_ind[s],delta_ind=delta_ind[s],
                                zeta=zeta_value[s],eta=eta_value[s],delta=delta_value[s],
                                sim=j)
  delta_slope[[j]]<- full_join(delta_slope[[j]],actual_delta, by='s')
  
  
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
    spread(quantile,y)%>%mutate(zeta_ind=zeta_ind[s],eta_ind=eta_ind[s],delta_ind=delta_ind[s],
                                zeta=zeta_value[s],eta=eta_value[s],delta=delta_value[s],
                                sim=j)
  zeta_slope_quantiles[[j]]<- full_join(zeta_slope_quantiles[[j]],actual_zeta, by='s')
  
}

zeta_slope_all<-unlist(zeta_slope_quantiles)%>%array(c(S,12,sims))
colnames(zeta_slope_all)<-c('s','2.5%','50%','97.5%','zeta_ind','eta_ind','delta_ind','zeta','eta','delta','sim','actual')
zeta_slope_all<-apply(zeta_slope_all,2,c)%>%as_tibble()%>%mutate(`2.5%`=as.numeric(`2.5%`),
                                                                 `50%`=as.numeric(`50%`),
                                                                 `97.5%`=as.numeric(`97.5%`),
                                                                 `actual`=as.numeric(`actual`))

delta_slope_all<-unlist(delta_slope)%>%array(c(S,12,sims))
colnames(delta_slope_all)<-c('s','2.5%','50%','97.5%','zeta_ind','eta_ind','delta_ind','zeta','eta','delta','sim','actual')
delta_slope_all<-apply(delta_slope_all,2,c)%>%as_tibble()%>%mutate(`2.5%`=as.numeric(`2.5%`),
                                                                   `50%`=as.numeric(`50%`),
                                                                   `97.5%`=as.numeric(`97.5%`),
                                                                   `actual`=as.numeric(`actual`))


eta_slope_all<-unlist(eta_slope)%>%array(c(S,12,sims))
colnames(eta_slope_all)<-c('s','2.5%','50%','97.5%','zeta_ind','eta_ind','delta_ind','zeta','eta','delta','sim','actual')
eta_slope_all<-apply(eta_slope_all,2,c)%>%as_tibble()%>%mutate(`2.5%`=as.numeric(`2.5%`),
                                                               `50%`=as.numeric(`50%`),
                                                               `97.5%`=as.numeric(`97.5%`),
                                                               `actual`=as.numeric(`actual`))

save(delta_slope_all,eta_slope_all,zeta_slope_all, file='matrix_simulations_parameter_inference.RData')

# Plots: simulated data vs model estimates 
plot_zeta_all_sims<-ggplot(data=zeta_slope_all)+
  geom_errorbar(aes(x=delta,  ymin=`2.5%`,ymax=`97.5%`,colour=as.factor(eta)),alpha=0.5)+
  facet_wrap(~paste('Totals temporal trend:',zeta))+
  geom_point(aes(x=delta, y=`50%`,colour=as.factor(eta)),shape=3,size=3)+
  geom_point(aes(x=delta, y=actual),size=3)+theme_minimal()+
  theme(text = element_text(size=12), plot.title = element_text(size=16),plot.subtitle=element_text(size=14),
        legend.title = element_text(size=14),strip.text = element_text(size=14),legend.position="bottom")+
  scale_colour_viridis_d(name='Temporal trend of reporting (eta)',begin=0.1,end=0.9)+
  scale_shape(name='Effect of case load on delay (delta)')+
  labs(title='Model estimates against simulated parameter values',
       subtitle = 'For linear temporal trend in total cases (zeta)',
       x='Effect of case load on delay (delta)',y='Linear trend in total cases (zeta)')
ggsave(plot_zeta_all_sims,file='Plots/plot_zeta_all_sims.pdf',width=9,height=4)       

plot_eta_all_sims<-ggplot(data=eta_slope_all)+
  geom_errorbar(aes(x=eta,  ymin=`2.5%`,ymax=`97.5%`,colour=as.factor(delta)),alpha=0.5)+
  facet_wrap(~paste('Totals temporal trend:',zeta))+
  geom_point(aes(x=eta, y=`50%`,colour=as.factor(delta)),shape=3,size=3)+
  geom_point(aes(x=eta, y=actual),size=3)+theme_minimal()+
  theme(text = element_text(size=12), plot.title = element_text(size=16),plot.subtitle=element_text(size=14),
        legend.title = element_text(size=14),strip.text = element_text(size=14),legend.position="bottom")+
  scale_colour_viridis_d(name='Effect of case load on delay (delta)',begin=0.1,end=0.9)+
  scale_shape(name='Effect of case load on delay (delta)')+
  labs(title='Model estimates against simulated parameter values',
       subtitle = 'For linear temporal trend in reporting delays (eta)',
       x='Temporal trend of reporting (eta)',y='Temporal trend of reporting (eta)')
ggsave(plot_eta_all_sims,file='Plots/plot_eta_all_sims.pdf',width=9,height=4)       

plot_delta_all_sims<-ggplot(data=delta_slope_all)+
  geom_errorbar(aes(x=delta,  ymin=`2.5%`,ymax=`97.5%`,colour=as.factor(eta)),alpha=0.5)+
  facet_wrap(~paste('Totals temporal trend:',zeta))+
  geom_point(aes(x=delta, y=`50%`,colour=as.factor(eta)),shape=3,size=3)+
  geom_point(aes(x=delta, y=actual),size=3)+theme_minimal()+
  theme(text = element_text(size=12), plot.title = element_text(size=16),plot.subtitle=element_text(size=14),
        legend.title = element_text(size=14),strip.text = element_text(size=14),legend.position="bottom")+
  scale_colour_viridis_d(name='Temporal trend of reporting (eta)',begin=0.1,end=0.9)+
  scale_shape(name='Effect of case load on delay (delta)')+
  labs(title='Model estimates against simulated parameter values',
       subtitle = 'For effect of case load on reporting delay (delta)',
       x='Effect of case load on delay (delta)',y='Effect of case load on delay (delta)')
ggsave(plot_delta_all_sims,file='Plots/plot_delta_all_sims.pdf',width=9,height=4)       



# delta_slope_combined<-rbind(delta_slope_all_up, delta_slope_all_stable, delta_slope_all_down)
# eta_slope_combined<-rbind(eta_slope_all_up, eta_slope_all_stable, eta_slope_all_down)
# zeta_slope_combined<-rbind(zeta_slope_all_up, zeta_slope_all_stable, zeta_slope_all_down)

# Simulation performance metrics 
delta_metric<-delta_slope_all%>%group_by(s,zeta,eta,delta)%>%
  summarise(Coverage=mean(actual>=`2.5%`&actual<=`97.5%`),
            `Prediction interval width`=mean(`97.5%`-`2.5%`),
            `Mean absolute percentage error`=mean(abs((`50%`-actual)/actual)),
            `Bias`=mean(`50%`-actual))%>%
  pivot_longer(cols=c('Coverage','Mean absolute percentage error','Prediction interval width','Bias'),names_to = 'metric')

eta_metric<-eta_slope_all%>%group_by(s,zeta,eta,delta)%>%
  summarise(Coverage=mean(actual>=`2.5%`&actual<=`97.5%`),
            `Prediction interval width`=mean(`97.5%`-`2.5%`),
            `Mean absolute percentage error`=mean(abs((`50%`-actual)/actual)),
            `Bias`=mean(`50%`-actual))%>%
  pivot_longer(cols=c('Coverage','Mean absolute percentage error','Prediction interval width','Bias'),names_to = 'metric')

zeta_metric<-zeta_slope_all%>%group_by(s,zeta,eta,delta)%>%
  summarise(Coverage=mean(actual>=`2.5%`&actual<=`97.5%`),
            `Prediction interval width`=mean(`97.5%`-`2.5%`),
            `Mean absolute percentage error`=mean(abs((`50%`-actual)/actual)),
            `Bias`=mean(`50%`-actual))%>%
  pivot_longer(cols=c('Coverage','Mean absolute percentage error','Prediction interval width','Bias'),names_to = 'metric')



# Simulation performance metrics 
delta_table<-delta_slope_all%>%group_by(zeta,eta,delta)%>%
  summarise(Coverage=mean(actual>=`2.5%`&actual<=`97.5%`),
            `Prediction interval width`=round(mean(`97.5%`-`2.5%`),3),
            `Bias`=round(mean(`50%`-actual),3))

colnames(delta_table)<-c('Cases temporal trend (zeta)','Delay temporal trend (eta)','Effect of cases on delay (delta)',
                         'Coverage (delta)','PIW (delta)','Bias (delta)')


eta_table<-eta_slope_all%>%group_by(zeta,eta,delta)%>%
  summarise(Coverage=mean(actual>=`2.5%`&actual<=`97.5%`),
            `Prediction interval width`=round(mean(`97.5%`-`2.5%`),3),
            `Bias`=round(mean(`50%`-actual),3))
colnames(eta_table)<-c('Cases temporal trend (zeta)','Delay temporal trend (eta)','Effect of cases on delay (delta)',
                       'Coverage (eta)','PIW (eta)','Bias (eta)')



zeta_table<-zeta_slope_all%>%group_by(zeta,eta,delta)%>%
  summarise(Coverage=mean(actual>=`2.5%`&actual<=`97.5%`),
            `Prediction interval width`=round(mean(`97.5%`-`2.5%`),3),
            `Bias`=round(mean(`50%`-actual),3))

colnames(zeta_table)<-c('Cases temporal trend (zeta)','Delay temporal trend (eta)','Effect of cases on delay (delta)',
                        'Coverage (zeta)','PIW (zeta)','Bias (zeta)')

# Custom colours for tables.
lightGreen = "#9BE2B7"
lightRed = "#F19999"
lightBlue="#9CD4EF"
darkGreen = "#32B767"
darkRed = "#DD5A5A"
darkBlue="#1491CD"



# Formats for text and fill of tables.
colour_formatter <- 
  formatter("span", 
            style = x ~ style(
              #font.weight = "bold",
              color = ifelse(x > 0, darkGreen, ifelse(x < 0, darkRed, darkBlue))))
fill_formatter <- 
  formatter("span", 
            style = x ~ style( display = "block", 
                               padding = "0 4px", 
                               `color` = "black", 
                               `border-radius` = "4px", 
                               `background-color` = ifelse(x>0, lightGreen, ifelse(x<0, lightRed, lightBlue))))


# 
# delta_formatter <- 
#   formatter("span", 
#             style = x ~ style( display = "block", 
#                                padding = "0 4px", 
#                                `color` = "white", 
#                                `border-radius` = "4px", 
#                                `background-color` = ifelse(x=="Positive effect", customGreen, ifelse(x=="Negative effect", customRed, customBlue))))

# Delta table:
formattable(delta_table, list(`Bias (delta)`=colour_formatter, 
                              `Delay temporal trend (eta)`=fill_formatter,
                              `Effect of cases on delay (delta)`=fill_formatter,
                              `Cases temporal trend (zeta)`=fill_formatter ))
# Eta table:
formattable(eta_table, list(`Bias (eta)`=colour_formatter, 
                            `Delay temporal trend (eta)`=fill_formatter,
                            `Effect of cases on delay (delta)`=fill_formatter,
                            `Cases temporal trend (zeta)`=fill_formatter ))
# Zeta table:
formattable(zeta_table, list(`Bias (zeta)`=colour_formatter, 
                             `Delay temporal trend (eta)`=fill_formatter,
                             `Effect of cases on delay (delta)`=fill_formatter,
                             `Cases temporal trend (zeta)`=fill_formatter ))
# width=930, height=920

# Colour blind friendly palatte:
red <- "#E6194B"
green <- "#2AA198"
blue <- "#4363D8"
yellow <- "#FFE119"
gray <- "#A9A9A9"

# Combine all tables:
ALL_table<-full_join(full_join(zeta_table,eta_table, by=colnames(delta_table)[1:3]),delta_table, by=colnames(delta_table)[1:3])
colnames(ALL_table)[1:3]<-c("Simulated zeta","Simulated eta", "Simulated delta")
formattable(ALL_table, list(`Bias (zeta)`=colour_formatter, 
                            `Bias (eta)`=colour_formatter, 
                            `Bias (delta)`=colour_formatter, 
                             `Simulated eta`=fill_formatter,
                             `Simulated delta`=fill_formatter,
                             `Simulated zeta`=fill_formatter ))

# Performance metric plots
delta_metric_plot<-ggplot(data=filter(delta_metric, metric!='Mean absolute percentage error'))+
  geom_point(aes(x=zeta,y=value,colour=delta,shape=eta),size=3,alpha=0.7)+facet_wrap(~metric, scales = 'free')+
  theme_light()+
  theme(text = element_text(size=14), axis.title = element_text(size=16), plot.title = element_text(size=20),plot.subtitle=element_text(size=18),
        legend.title = element_text(size=16),strip.text = element_text(size=16),
        axis.text = element_text(size=12),legend.position="bottom")+guides(colour=guide_legend(nrow=1))+
  scale_colour_manual(name='Temporal trend of reporting (eta)',values=c(red,blue,green))+
  scale_shape_manual(name='Effect of case load on delay (delta)',values=c(0, 8, 17))+
  labs(x='Linear temporal trend in total cases (zeta)',title='Simulation performance metrics',
       subtitle='Values for effect of case load on delay (delta)',y=NULL)
ggsave(delta_metric_plot,file='Plots/delta_metric_plot.pdf',width=12,height=5)

eta_metric_plot<-ggplot(data=filter(eta_metric, metric!='Mean absolute percentage error'))+
  geom_point(aes(x=zeta,y=value,colour=delta,shape=eta),size=3,alpha=0.7)+facet_wrap(~metric, scales = 'free')+
  theme_light()+
  theme(text = element_text(size=14), axis.title = element_text(size=16), plot.title = element_text(size=20),plot.subtitle=element_text(size=18),
        legend.title = element_text(size=16),strip.text = element_text(size=16),
        axis.text = element_text(size=12),legend.position="bottom")+guides(colour=guide_legend(nrow=1))+
  scale_colour_manual(name='Temporal trend of reporting (eta)',values=c(red,blue,green))+
  scale_shape_manual(name='Effect of case load on delay (delta)',values=c(0, 8, 17))+
  labs(x='Linear temporal trend in total cases (zeta)',title='Simulation performance metrics',
       subtitle='Values for linear temporal trend of reporting (eta)',y=NULL)
ggsave(eta_metric_plot,file='Plots/eta_metric_plot.pdf',width=12,height=5)

zeta_metric_plot<-ggplot(data=filter(zeta_metric, metric!='Mean absolute percentage error'))+
  geom_point(aes(x=zeta,y=value,colour=delta,shape=eta),size=3,alpha=0.7)+facet_wrap(~metric, scales = 'free')+
  theme_light()+
  theme(text = element_text(size=14), axis.title = element_text(size=16), plot.title = element_text(size=20),plot.subtitle=element_text(size=18),
        legend.title = element_text(size=16),strip.text = element_text(size=16),
        axis.text = element_text(size=12),legend.position="bottom")+guides(colour=guide_legend(nrow=1))+
  scale_colour_manual(name='Temporal trend of reporting (eta)',values=c(red,blue,green))+
  scale_shape_manual(name='Effect of case load on delay (delta)',values=c(0, 8, 17))+
  labs(x='Linear temporal trend in total cases (zeta)',title='Simulation performance metrics',
       subtitle ='Values for linear temporal trend in total cases (zeta)',y=NULL)
ggsave(zeta_metric_plot,file='Plots/zeta_metric_plot.pdf',width=12,height=5)


bias_metric<-rbind(filter(delta_metric,metric=='Bias')%>%mutate(parameter='delta'),filter(eta_metric,metric=='Bias')%>%mutate(parameter='eta'),filter(zeta_metric,metric=='Bias')%>%mutate(parameter='zeta'))
bias_metric_plot<-ggplot(data=bias_metric)+
  geom_point(aes(x=zeta,y=value,colour=delta,shape=eta),size=3,alpha=0.7)+
  facet_wrap(~parameter)+theme_light()+
  geom_hline(yintercept = 0, colour="darkgrey")+
  theme(text = element_text(size=14), axis.title = element_text(size=16), plot.title = element_text(size=20),plot.subtitle=element_text(size=18),
        legend.title = element_text(size=16),strip.text = element_text(size=16),
        axis.text = element_text(size=12),legend.position="bottom")+guides(colour=guide_legend(nrow=1))+
  scale_colour_manual(name='Temporal trend of reporting (eta)',values=c(red,blue,green))+
  scale_shape_manual(name='Effect of case load on delay (delta)',values=c(6, 8, 17))+
  labs(x='Linear temporal trend in total cases (zeta)',title='Bias performance metric',
       subtitle ='Values for the three paramerers of interest',y=NULL)
bias_metric_plot
ggsave(bias_metric_plot,file='Plots/bias_metric_plot.pdf',width=12,height=5)


Coverage_metric<-rbind(filter(delta_metric,metric=='Coverage')%>%mutate(parameter='delta'),filter(eta_metric,metric=='Coverage')%>%mutate(parameter='eta'),filter(zeta_metric,metric=='Coverage')%>%mutate(parameter='zeta'))
Coverage_metric_plot<-ggplot(data=Coverage_metric)+
  geom_point(aes(x=zeta,y=value,colour=delta,shape=eta),size=3,alpha=0.7)+
  facet_wrap(~parameter, scales = 'free')+theme_light()+
  theme(text = element_text(size=14), axis.title = element_text(size=16), plot.title = element_text(size=20),plot.subtitle=element_text(size=18),
        legend.title = element_text(size=16),strip.text = element_text(size=16),
        axis.text = element_text(size=12),legend.position="bottom")+guides(colour=guide_legend(nrow=1))+
  scale_colour_manual(name='Temporal trend of reporting (eta)',values=c(red,blue,green))+
  scale_shape_manual(name='Effect of case load on delay (delta)',values=c(6, 8, 17))+
  labs(x='Linear temporal trend in total cases (zeta)',title='Coverage performance metric',
       subtitle ='Values for the three paramerers of interest',y=NULL)
ggsave(Coverage_metric_plot,file='Plots/Coverage_metric_plot.pdf',width=12,height=5)


piw_metric<-rbind(filter(delta_metric,metric=='Prediction interval width')%>%mutate(parameter='delta'),filter(eta_metric,metric=='Prediction interval width')%>%mutate(parameter='eta'),filter(zeta_metric,metric=='Prediction interval width')%>%mutate(parameter='zeta'))
piw_metric_plot<-ggplot(data=piw_metric)+
  geom_point(aes(x=zeta,y=value,colour=delta,shape=eta),size=3,alpha=0.7)+
  facet_wrap(~parameter, scales = 'free')+theme_light()+
  theme(text = element_text(size=14), axis.title = element_text(size=16), plot.title = element_text(size=20),plot.subtitle=element_text(size=18),
        legend.title = element_text(size=16),strip.text = element_text(size=16),
        axis.text = element_text(size=12),legend.position="bottom")+guides(colour=guide_legend(nrow=1))+
  scale_colour_manual(name='Temporal trend of reporting (eta)',values=c(red,blue,green))+
  scale_shape_manual(name='Effect of case load on delay (delta)',values=c(6, 8, 17))+
  labs(x='Linear temporal trend in total cases (zeta)',title='Prediction interval width performance metric',
       subtitle ='Values for the three paramerers of interest',y=NULL)
ggsave(piw_metric_plot,file='Plots/piw_metric_plot.pdf',width=12,height=5)

piw_metric_box_plot<-ggplot(data=piw_metric)+
  geom_boxplot(aes(x=zeta,y=value,colour=delta),alpha=0.1)+
  geom_point(aes(x=zeta,y=value,colour=delta,shape=eta),size=3,alpha=0.7)+
  facet_wrap(~parameter, scales = 'free')+theme_light()+
  theme(text = element_text(size=14), axis.title = element_text(size=16), plot.title = element_text(size=20),plot.subtitle=element_text(size=18),
        legend.title = element_text(size=16),strip.text = element_text(size=16),
        axis.text = element_text(size=12),legend.position="bottom")+guides(colour=guide_legend(nrow=1))+
  scale_colour_manual(name='Temporal trend of reporting (eta)',values=c(red,blue,green))+
  scale_shape_manual(name='Effect of case load on delay (delta)',values=c(6, 8, 17))+
  labs(x='Linear temporal trend in total cases (zeta)',title='Prediction interval width performance metric',
       subtitle ='Values for the three paramerers of interest',y=NULL)
ggsave(piw_metric_box_plot,file='Plots/piw_metric_box_plot.pdf',width=12,height=5)

