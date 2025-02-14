# The following packages will need to be installed and loaded.

library(tidyverse)
library(reshape2)
library(mgcv)
library(coda)
library(gridExtra)
library(viridis)
library(scales)
library(abind)
library(nimble)
library(doParallel)
library(lubridate)
library(openxlsx)
#library(INLA)
library(readr)

# Read in some necessary functions.
#source('Functions.R')

# Regional populations in England:
# Reading in the table from Wikipedia
page = read_html("https://en.wikipedia.org/wiki/Regions_of_England")
# Obtain the piece of the web page that corresponds to the "wikitable" node
my.table = html_node(page, ".wikitable")
# Convert the html table eleument into a data frame
Engalnd_region_metrics= html_table(my.table, fill = TRUE)[-1,]


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

# Vector containing region names.
region_names <- c('East of England','London','Midlands',
                  'North East and Yorkshire',
                  'North West','South East','South West')

population_table<-tibble(Region=region_names, population=as.numeric(0), s=1:7)
Engalnd_region_metrics[,6]<-as.numeric(str_replace_all(as.character(unlist(Engalnd_region_metrics[,6])),",",""))
#East of England
population_table[1,2]<-Engalnd_region_metrics[6,6]
#London
population_table[2,2]<-Engalnd_region_metrics[7,6]
#Midlands = East Midlands + West Midlands
population_table[3,2]<-Engalnd_region_metrics[4,6]+Engalnd_region_metrics[5,6]
#North East and Yorkshire = North East + Yorkshire and the Humber
population_table[4,2]<-Engalnd_region_metrics[1,6]+Engalnd_region_metrics[3,6]
#North West
population_table[5,2]<-Engalnd_region_metrics[2,6]
#South East
population_table[6,2]<-Engalnd_region_metrics[8,6]
#South West
population_table[7,2]<-Engalnd_region_metrics[9,6]

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


## OPTIMISATION BY REGION CLUSTER
#NB marginal model 

# Set up the splines using jagam.
blank_data=data_frame(y=rnorm(N,0,1),t=1:N,w=(1:N)%%7)
blank_jagam=jagam(y~s(t,bs='cs',k=n_knots[1])+s(w,bs='cc',k=n_knots[2]),data=blank_data,file='blank.jags',
                  knots=list(t=seq(1,C,length=n_knots[1]),w=seq(0,7,length=n_knots[2])))


#PARTIAL COUNTS
z <- z_full[1:N,,]
for(s in 1:S){
  z[,s,][outer(1:dim(z[,s,])[1], 0:(dim(z[,s,])[2]-1), FUN = "+") > C] <- NA
}

# Create Optimisation Cluster
library(parallel)


# Function that optimises each region seperately 
run_optim_allcode <- function( i , N, S, C, D, n_knots, partial_counts, population_table) {
  
  #load libraries
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
  
  # Set up the splines using jagam.
  blank_data=data_frame(y=rnorm(N,0,1),t=1:N,w=(1:N)%%7)
  blank_jagam=jagam(y~s(t,bs='cs',k=n_knots[1])+s(w,bs='cc',k=n_knots[2]),data=blank_data,file='blank.jags',
                    knots=list(t=seq(1,C,length=n_knots[1]),w=seq(0,7,length=n_knots[2])))
  
  
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
  #registerDistributions(list(dbetabin=list(
  #BUGSdist='dbetabin(mu,phi,size)',discrete=TRUE)))
  
  #partial counts 
  D_max<-dim(partial_counts)[3]
  z<-array(partial_counts[,i,],dim=c(N,S,D_max))
  D_rem<-D+1
  
  #unobserved y
  y<-apply(partial_counts,c(1,2),sum)[1:C,]
  unobs_y<-which(is.na(y[,i])) 
  
  #remainder term for z's
  remainder=matrix((apply(z,2,rowSums, na.rm=TRUE)-apply(array(z[,,1:D],dim=c(N,S,D)),2,rowSums, na.rm=TRUE))[1:C,], ncol=S)
  
  if(is.na(unobs_y[1])==FALSE){remainder[unobs_y[1]:C,]<-NA}# reminder unknown when y is unknown
  z=array(abind(array(z[1:C,,1:D],dim=c(C,S,D)), remainder)[1:C,,1:(D+1)], dim=c(C,S,D+1))
  
  obs_index<-matrix(NA, nrow=S, ncol=D_rem)
  for(s in 1:S){
    for(d in 1:D_rem){
      if(sum(is.na(z[,s,d]))>0){obs_index[s,d]<-which(is.na(z[,s,d])==TRUE)[1]-1}else{
        obs_index[s,d]<-C
      }
    }
  }
  
  #NB model for partial counts
  covid_gamma_code <- nimbleCode({
    for(s in 1:S){
      for(t in 1:N){
        # Mean total deaths.
        log(lambda[t,s]) <- log(Pop[s]) + iota[s] + delta[t,s]
        for(d in 1:(D-1)){
          # Expected cumulative proportions.
          probit(p[t,s,d]) <- beta[s,d] + gamma[t,s] + xi[t,s]
        }
        p[t,s,D]<-1
        # Relative proportions (Negative-Binomial proportions).
        mu[t,s,1]<- p[t,s,1]
        for(d in 2:(D)){
          mu[t,s,d] <- p[t,s,d] - p[t,s,(d-1)]
        }
      }
    }
    for(s in 1:S){
      for(d in 1:D)
        for(t in 1:obs_index[s,d]){
          prob[t,s,d]<-(phi[s,d])/(phi[s,d]+lambda[t,s]*mu[t,s,d])
          z[t,s,d] ~ dnegbin( prob[t,s,d], phi[s,d])  
        }
    }
    
    
    ## Overall effects
    ## Regional effects
    for(s in 1:S){
      Omega_delta[1:K_t,1:K_t,s] <- S_t[1:K_t,1:K_t]*exp(log_tau_delta[s])
      kappa_delta[1:K_t,s] ~ dmnorm(zeros[1:K_t],Omega_delta[1:K_t,1:K_t,s])
      delta[1:N,s] <- X_t[1:N,1:K_t]%*%kappa_delta[1:K_t,s]
      # Delay spline
      Omega_gamma[1:K_t,1:K_t,s] <- S_t[1:K_t,1:K_t]*exp(log_tau_gamma[s])
      kappa_gamma[1:K_t,s] ~ dmnorm(zeros[1:K_t],Omega_gamma[1:K_t,1:K_t,s])
      gamma[1:N,s] <- X_t[1:N,1:K_t]%*%kappa_gamma[1:K_t,s]
      # Time of week spline
      Omega_xi[1:K_w,1:K_w,s] <- S_w[1:K_w,1:K_w]*exp(log_tau_xi[s])
      kappa_xi[1:K_w,s] ~ dmnorm(zeros[1:K_w],Omega_xi[1:K_w,1:K_w,s])
      xi[1:N,s] <- X_w[1:N,1:K_w]%*%kappa_xi[1:K_w,s]
    }
    
    # Now need to provide initial values for log_tau_alpha etc not tau_alpha
    for(s in 1:S){
      log_tau_delta[s] ~ dnorm(1.26,sd=2.23)
      log_tau_gamma[s] ~ dnorm(1.26,sd=2.23)
      log_tau_xi[s] ~ dnorm(1.26,sd=2.23)
      for(d in 1:D){
        log_phi[s,d] ~ dnorm(4.3,sd=0.8)
        phi[s,d]<-exp(log_phi[s,d])
        
      }
      omega_beta[s,1] ~ dnorm(0,sd=10)
      beta[s,1]<-omega_beta[s,1]
      for(d in 2:(D-1)){
        omega_beta[s,d] ~ dnorm(0,sd=2)
        beta[s,d]<-beta[s,d-1]+exp(omega_beta[s,d])  # Independent delay effects
      }
      iota[s] ~ dnorm(-10,sd=10) # Spatial intercept.
      
    }
  })
  
  
  # Constants (e.g. number of days to model, number of knots) for NIMBLE.
  covid_constants_NBz <- list(N=N,C=C,S=S,D=D_rem,K_t=n_knots[1]-1,K_w=n_knots[2]-2, obs_index=obs_index, Pop=as.numeric(unlist(population_table[,2])))
  
  # Data (e.g. total deaths, partial counts, spline model matrix) for NIMBLE.
  covid_data_NBz <- list(z=array(z[1:C,1:S,1:D_rem], dim=c(C,S,D_rem)),
                         X_t=blank_jagam$jags.data$X[,2:(n_knots[1])],
                         X_w=blank_jagam$jags.data$X[,(n_knots[1]+1):(n_knots[1]+n_knots[2]-2)],
                         S_t=blank_jagam$jags.data$S1,
                         S_w=blank_jagam$jags.data$S2,
                         zeros=rep(0,max(covid_constants_NBz$K_t,covid_constants_NBz$K_w)))
  
  
  
  # Generate random initial values.
  covid_inits_NBz <- list( kappa_gamma=matrix(rnorm(S*covid_constants_NBz$K_t,0,0.1),ncol=S),
                           kappa_xi=matrix(rnorm(S*covid_constants_NBz$K_w,0,0.1),ncol=S),
                           kappa_delta=matrix(rnorm(S*covid_constants_NBz$K_t,0,0.1),ncol=S),#
                           log_tau_delta=rnorm(S,1.26,sd=2.23),#
                           iota=rnorm(S,-10,1),#
                           # log_theta=rnorm(S,4,1.3),
                           log_tau_gamma=rnorm(S,1.26,sd=2.23),
                           log_tau_xi=rnorm(S,1.26,sd=2.23),
                           log_phi=matrix(rnorm(S*D_rem,4,1.3),nrow=S,ncol=D_rem),
                           omega_beta=matrix(rnorm((D_rem-1)*S,0,0.1),nrow=S,ncol=(D_rem-1)))
  
  
  # Build the model.
  covid_model_gamma <- nimbleModel(covid_gamma_code,covid_constants_NBz,covid_data_NBz,covid_inits_NBz)
  # Compile the model.
  covid_compiled_model_gamma <- compileNimble(covid_model_gamma)
  # Set up the MCMC.
  covid_mcmc_config_gamma<- configureMCMC(covid_model_gamma, monitors=c('beta','xi','phi','gamma','delta','iota', 'lambda'),useConjugacy = FALSE)
  
  ## gamma objective function
  
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
  
  ## OPTIMISATION 
  #paramers to optimise
  target_all <- c(covid_mcmc_config_gamma$getSamplers()%>%lapply(function(x)x$target)%>%unlist())
  n_par <- length(values(covid_model_gamma,target_all))
  
  # Initial vlaues for the optimisation
  opt_inits <- rep(0.01, n_par)
  
  rObjective <- objective(covid_model_gamma, target_all)
  cObjective <- compileNimble(rObjective, project = covid_model_gamma)
  #function for estimating gradient 
  grad.fun <- function(x){
    grad(func=cObjective$run, x=x, method="simple") #"Richardson"\simple
  }
  max_it<-5000
  # run and time optimisation
  
  optC_gamma<-trust.optim(fn = cObjective$run, gr= grad.fun,  x = opt_inits, control = list(function.scale.factor = -1, maxit=max_it), method = c("BFGS")) #hs= #only needed for sparse #other method: "SR1" "BFGS"
  
  
  #model constant
  D_rem<-covid_constants_NBz$D
  D<-covid_constants_NBz$D
  S<-covid_constants_NBz$S
  N<-covid_constants_NBz$N
  C<-covid_constants_NBz$C
  
  target<-target_all
  par_z<-data.frame()
  par_z[1,1:length(target)]<-NA  
  colnames(par_z)<-target # data frame for each target entry
  model_values<-rep(0, n_par)
  for(j in 1:length(target)){   
    values(covid_model_gamma, target_all[j])<-NA #identify parameter j values
    par_index<-which(is.na(values(covid_model_gamma, target_all))) # index for parameter j values
    par_vector<-optC_gamma$solution[par_index] # just parameter j values 
    par_z[1:length(par_vector),j]<-par_vector
    values(covid_model_gamma, target_all)<-model_values
  }
  # get parameter estimates 
  
  #NOTE: look at par_z header to find out if parameters are in nice order 
  
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
  kappa_gamma<- matrix(NA, nrow=covid_constants_NBz$K_t, ncol=S)
  for(s in 1:S){
    kappa_gamma[,s] <-as.numeric(na.omit(par_z[,which_kappa_gamma[s]]))
  }
  gamma_opt<-matrix(NA, nrow=N, ncol=S)
  for(s in 1:S){
    gamma_opt[,s] <- covid_data_NBz$X_t[1:N,1:covid_constants_NBz$K_t]%*%kappa_gamma[,s]
  }
  
  #calculate xi
  is_kappa_xi<-sapply(target_all,function(x)substr(x,1,8))=='kappa_xi'
  which_kappa_xi<-as.numeric(which(is_kappa_xi))
  kappa_xi<- matrix(NA, nrow=covid_constants_NBz$K_w, ncol=S)
  for(s in 1:S){
    kappa_xi[,s] <-as.numeric(na.omit(par_z[,which_kappa_xi[s]]))
  }
  xi_opt<-matrix(NA, nrow=N, ncol=S)
  for(s in 1:S){
    xi_opt[,s] <- covid_data_NBz$X_w[1:N,1:covid_constants_NBz$K_w]%*%kappa_xi[,s]
  }          
  
  # calculate phi
  phi_opt<-matrix(NA,nrow=S,ncol=D_rem)
  phi_start<-which(target_all=='log_phi[1, 1]')
  phi_opt<-matrix(exp(as.numeric(par_z[1,phi_start:(phi_start+D_rem*S-1)])), byrow=TRUE, ncol=D_rem, nrow=S)
  
  #calculate  delta's
  is_kappa_delta<-sapply(target_all,function(x)substr(x,1,11))=='kappa_delta'
  which_kappa_delta<-as.numeric(which(is_kappa_delta))
  kappa_delta<- matrix(NA, nrow=covid_constants_NBz$K_t, ncol=S)
  for(s in 1:S){
    kappa_delta[,s] <-as.numeric(na.omit(par_z[,which_kappa_delta[s]]))
  }
  delta_opt<-matrix(NA, nrow=N, ncol=S)
  for(s in 1:S){
    delta_opt[,s] <- covid_data_NBz$X_t[1:N,1:covid_constants_NBz$K_t]%*%kappa_delta[,s]
  }   
  
  #iota
  iota_start<-which(target_all=='iota[1]')
  iota_opt<-as.numeric(par_z[1,iota_start:(iota_start+S-1)])
  
  
  p_opt<-array(NA, dim=c(N,S,D_rem))
  mu_opt<-array(NA, dim=c(N,S,D_rem))
  lambda_opt<-array(NA, dim=c(N,S))
  for(t in 1:N){
    for(s in 1:S){
      #log(population[s]) + iota[s] + delta[t,s] + xi[w[t],s]
      lambda_opt[t,s]<-exp(log(covid_constants_NBz$Pop[s])+iota_opt[s]+delta_opt[t,s])
      for(d in 1:(D_rem-1)){
        p_opt[t,s,d]<-pnorm(beta_opt[s,d] + gamma_opt[t,s]+xi_opt[t,s]) 
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
  
  lambda_opt_gamma<-ceiling(lambda_opt)
  y_init_gamma<-rbind(matrix(ceiling(apply(z_opt,2,rowSums)),nrow=C),matrix(lambda_opt_gamma[(C+1):N,], ncol=S))
  
  
  phi_opt_gamma<-phi_opt
  p_opt_gamma<-p_opt
  
  start_tau_gamma<-which(target_all=='log_tau_gamma[1]')
  start_tau_delta<-which(target_all=='log_tau_delta[1]')
  start_tau_xi<-which(target_all=='log_tau_xi[1]')
  
  y_na<-y_init_gamma[1:C,1:S] #use optimised y
  not_na<-is.na(apply(partial_counts,c(1,2),sum)[1:C,i])==0#observed y
  y_na[not_na]<-rep(NA,sum(not_na)) #NA's for observed y's
  
  #set inits
  covid_inits_gam <- list(kappa_delta=kappa_delta,
                          kappa_gamma=kappa_gamma,
                          kappa_xi=kappa_xi,
                          iota=as.numeric(par_z[1,iota_start:(iota_start+S-1)]),
                          log_tau_gamma=as.numeric(par_z[1,start_tau_gamma:c(start_tau_gamma+S-1)]),
                          log_tau_delta=as.numeric(par_z[1,start_tau_delta:c(start_tau_delta+S-1)]),
                          log_tau_xi=as.numeric(par_z[1,start_tau_xi:c(start_tau_xi+S-1)]),
                          log_phi=log(phi_opt[1:S,1:(D_rem-1)]),
                          omega_beta=omega_beta_opt[1:S,1:(D_rem-1)],
                          y=y_na)
  #function return:
  list(optim_inits=covid_inits_gam,y_init_gamma=y_init_gamma,lambda_opt=lambda_opt)
}

this_cluster_optim<- makeCluster(n_chains)

S_all<-7
D_max<-14
time_optimcluster2<-system.time({
  optim_output  <- parLapply(cl = this_cluster_optim, X = 1:S_all, fun = run_optim_allcode, 
                             N=N, S=1, D=6, C=C, n_knots=n_knots, 
                             partial_counts=z[1:N,1:S_all,1:D_max],
                             population_table=population_table)
})

stopCluster(this_cluster_optim)

#Plot Optimisation 
stat_t=1
for(i in 1:S_all){
  plot(stat_t:N, y_full[stat_t:N,i])
  lines(stat_t:N, optim_output[[i]]$lambda_opt[stat_t:N,], col='blue') #lambda optimised
  title(paste('Region',i,': optimised y and mean'))
  abline(v=C) #Start of forecast 
  lines(stat_t:N, optim_output[[i]]$y_init_gamma[stat_t:N,], col='green') #z's sum
}

## SET INTIAL VALUES 
kappa_delta<-matrix(NA, nrow=nrow(optim_output[[1]]$optim_inits$kappa_delta), ncol=S)
kappa_gamma<-matrix(NA, nrow=nrow(optim_output[[1]]$optim_inits$kappa_gamma), ncol=S)
kappa_xi<-matrix(NA, nrow=nrow(optim_output[[1]]$optim_inits$kappa_xi), ncol=S)
iota<-rep(NA,S)
#log_theta<-rep(NA,S)
log_tau_gamma<-rep(NA,S)
log_tau_delta<-rep(NA,S)
log_tau_xi<-rep(NA,S)
log_phi<-matrix(NA, nrow=S, ncol=(D))
omega_beta<-matrix(NA, nrow=S, ncol=D)
y_init<-matrix(NA, nrow=C, ncol=S)

# Combine regions
for(i in 1:S){
  kappa_delta[,i]<-optim_output[[i]]$optim_inits$kappa_delta
  kappa_gamma[,i]<-optim_output[[i]]$optim_inits$kappa_gamma
  kappa_xi[,i]<-optim_output[[i]]$optim_inits$kappa_xi
  iota[i]<-optim_output[[i]]$optim_inits$iota
  #log_theta[i]<-optim_output[[i]]$optim_inits$log_theta
  log_tau_gamma[i]<-optim_output[[i]]$optim_inits$log_tau_gamma
  log_tau_delta[i]<-optim_output[[i]]$optim_inits$log_tau_delta
  log_tau_xi[i]<-optim_output[[i]]$optim_inits$log_tau_xi
  log_phi[i,]<-optim_output[[i]]$optim_inits$log_phi
  omega_beta[i,]<-optim_output[[i]]$optim_inits$omega_beta
  y_init[,i]<-optim_output[[i]]$optim_inits$y
}

kappa_eta <-apply(kappa_xi, 1, mean)
kappa_alpha <-apply(kappa_delta, 1, mean)
kappa_psi <-apply(kappa_gamma, 1, mean)
log_tau_eta<-mean(log_tau_xi)
log_tau_alpha<-mean(log_tau_delta)
log_tau_psi<-mean(log_tau_gamma)

## Intial values for MCMC:
covid_inits_gam <- list(kappa_alpha=kappa_alpha,
                        kappa_psi=kappa_psi,
                        kappa_eta=kappa_eta,
                        kappa_delta=kappa_delta,
                        kappa_gamma=kappa_gamma,
                        kappa_xi=kappa_xi,
                        iota=iota,             
                        log_theta=rnorm(S,4,1.3),
                        log_tau_gamma=log_tau_gamma,
                        log_tau_delta=log_tau_delta,
                        log_tau_xi=log_tau_xi,
                        log_phi=log_phi[,1:D],
                        log_tau_psi=log_tau_psi,
                        log_tau_alpha=log_tau_alpha,
                        log_tau_eta=log_tau_eta,
                        omega_beta=omega_beta,
                        y=y_init)



# GDM MODEL
D=6 
z=z[1:C,,1:D]
y=y[1:C,]

#observed index for each region at each delay
obs_index<-matrix(NA, nrow=S, ncol=D)
for(s in 1:S){
  for(d in 1:D){
    obs_index[s,d]<-which(is.na(z[,s,d])==TRUE)[1]-1
  }
}


# Constants (e.g. number of days to model, number of knots) for NIMBLE.
covid_constants <- list(N=N,C=C,S=S,D=D,K_t=n_knots[1]-1,K_w=n_knots[2]-2, obs_index=obs_index, Pop=as.numeric(unlist(population_table[,2])) )
# Data (e.g. total deaths, partial counts, spline model matrix) for NIMBLE.
covid_data <- list(z=z[1:C,,1:D],y=y[1:C,],X_t=blank_jagam$jags.data$X[,2:(n_knots[1])],
                   X_w=blank_jagam$jags.data$X[,(n_knots[1]+1):(n_knots[1]+n_knots[2]-2)],
                   S_t=blank_jagam$jags.data$S1,S_w=blank_jagam$jags.data$S2,
                   zeros=rep(0,max(covid_constants$K_t,covid_constants$K_w)))


#MAKE CLUSTER FOR MCMC
this_cluster_gam<- makeCluster(n_chains)

### cluster function
run_MCMC_allcode <- function( seed, covid_constants, covid_inits, covid_data) {
  
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
  #registerDistributions(list(dbetabin=list(
  #BUGSdist='dbetabin(mu,phi,size)',discrete=TRUE)))
  
  
  
  N<-covid_constants$N
  C<-covid_constants$C
  S<-covid_constants$S
  D<-covid_constants$D
  n_chains <- 4 # Number of MCMC chains to run. Reduce if you don't have more than 4 CPU cores!
  n_knots <- c(8,8) # Number of knots of the temporal and seasonal splines, respectively.
  
  
  # NIMBLE model code. 
  covid_code_cluster <- nimbleCode({
    for(s in 1:S){
      for(t in 1:N){
        # Mean total deaths.
        log(lambda[t,s]) <- log(Pop[s]) + iota[s] + delta[t,s]
        for(d in 1:D){
          # Expected cumulative proportions.
          probit(p[t,s,d]) <- beta[s,d] + gamma[t,s] + xi[t,s]
        }
        # Relative proportions (Beta-Binomial means).
        nu[t,s,1] <- p[t,s,1]
        for(d in 2:D){
          nu[t,s,d] <- (p[t,s,d]-p[t,s,d-1])/(1-p[t,s,d-1])
        }
      }
      for(t in 1:C){
        # Model for total deaths.
        y[t,s] ~ dnegbin(theta[s]/(theta[s]+lambda[t,s]),theta[s])
        # Model for partial reports.
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
    Omega_alpha[1:K_t,1:K_t] <- S_t[1:K_t,1:K_t]*exp(log_tau_alpha)  #tau_alpha
    kappa_alpha[1:K_t] ~ dmnorm(zeros[1:K_t],Omega_alpha[1:K_t,1:K_t])
    alpha[1:N] <- X_t[1:N,1:K_t]%*%kappa_alpha[1:K_t]
    # Delay spline
    Omega_psi[1:K_t,1:K_t] <- S_t[1:K_t,1:K_t]*exp(log_tau_psi)  #tau_psi
    kappa_psi[1:K_t] ~ dmnorm(zeros[1:K_t],Omega_psi[1:K_t,1:K_t])
    psi[1:N] <- X_t[1:N,1:K_t]%*%kappa_psi[1:K_t]
    # Time of week spline
    Omega_eta[1:K_w,1:K_w] <- S_w[1:K_w,1:K_w]*exp(log_tau_eta)  #tau_eta
    kappa_eta[1:K_w] ~ dmnorm(zeros[1:K_w],Omega_eta[1:K_w,1:K_w])
    eta[1:N] <- X_w[1:N,1:K_w]%*%kappa_eta[1:K_w]
    ## Regional effects
    for(s in 1:S){
      Omega_delta[1:K_t,1:K_t,s] <- S_t[1:K_t,1:K_t]*exp(log_tau_delta[s]) #tau_delta[s]#
      kappa_delta[1:K_t,s] ~ dmnorm(kappa_alpha[1:K_t],Omega_delta[1:K_t,1:K_t,s])
      delta[1:N,s] <- X_t[1:N,1:K_t]%*%kappa_delta[1:K_t,s]
      # Delay spline
      Omega_gamma[1:K_t,1:K_t,s] <- S_t[1:K_t,1:K_t]*exp(log_tau_gamma[s])  #tau_gamma[s]#
      kappa_gamma[1:K_t,s] ~ dmnorm(kappa_psi[1:K_t],Omega_gamma[1:K_t,1:K_t,s])
      gamma[1:N,s] <- X_t[1:N,1:K_t]%*%kappa_gamma[1:K_t,s]
      # Time of week spline
      Omega_xi[1:K_w,1:K_w,s] <- S_w[1:K_w,1:K_w]*exp(log_tau_xi[s])  #tau_xi[s]#
      kappa_xi[1:K_w,s] ~ dmnorm(kappa_eta[1:K_w],Omega_xi[1:K_w,1:K_w,s])
      xi[1:N,s] <- X_w[1:N,1:K_w]%*%kappa_xi[1:K_w,s]
    }
    # Smoothing parameter priors.
    # tau_alpha ~ dinvgamma(0.5,0.5) # Equivalent to Half-Normal(0,1) on 1/sqrt(tau).
    #tau_psi ~ dinvgamma(0.5,0.5)
    #tau_eta ~ dinvgamma(0.5,0.5)
    log_tau_alpha ~ dnorm(1.26,sd=2.23)
    log_tau_psi ~ dnorm(1.26,sd=2.23)
    log_tau_eta ~ dnorm(1.26,sd=2.23)
    for(s in 1:S){
      #tau_delta[s] ~ dinvgamma(0.5,0.5)
      #tau_gamma[s] ~ dinvgamma(0.5,0.5)
      #tau_xi[s] ~ dinvgamma(0.5,0.5)
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
      for(d in 2:(D)){
        omega_beta[s,d] ~ dnorm(0,sd=2)
        beta[s,d]<-beta[s,d-1]+exp(omega_beta[s,d])  # Independent delay effects
      }
      # beta[s,1] ~ dnorm(0,sd=10) # Independent delay effects.
      #  for(d in 2:D){
      #   beta[s,d] ~ T(dnorm(beta[s,d-1],sd=10),beta[s,d-1],)
      # }
      iota[s] ~ dnorm(0,sd=10) # Spatial intercept.
      #theta[s] ~ dgamma(2,0.02) # Negative-Binomial dispersion parameters.
      log_theta[s] ~ dnorm(4.3,sd=0.8)
      theta[s]<-exp(log_theta[s])
    }
  })
  
  
  
  
  # Set up the NIMBLE model for each chain.
  
  
  
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
  #
  # Build the model.
  covid_model_cluster <-nimbleModel(covid_code_cluster,covid_constants,covid_data,covid_inits)
  # Compile the model.
  covid_compiled_model_cluster <- compileNimble(covid_model_cluster)
  
  # Set up the MCMC.
  covid_mcmc_config_cluster <-configureMCMC(covid_model_cluster,
                                            monitors=c('alpha','iota','beta','xi','psi','eta','y',
                                                       'theta','phi','gamma','delta','lambda'))
  
  ## OPTIMISE DISPERSION PARAMETERS 
  target_full <- c(covid_mcmc_config_cluster$getSamplers()%>%lapply(function(x)x$target)%>%unlist())
  is_theta <- sapply(target_full,function(x)substr(x,1,7))=='log_the'
  is_phi <- sapply(target_full,function(x)substr(x,1,7))=='log_phi'
  which_theta<-which(is_theta)
  which_phi<-which(is_phi)
  n_phi<-length(which_phi)
  target_disp<-target_full[sort(c(which_theta,which_phi))]
  n_par2 <- length(values(covid_model_cluster,target_disp))
  
  #initial values for the optimisation
  opt_inits2 <- rnorm(n_par2,0.01,0.1)
  
  rObjective <- objective(covid_model_cluster, target_disp)
  cObjective <- compileNimble(rObjective, project = covid_model_cluster)
  
  #function for estimating gradient 
  grad.fun <- function(x){
    grad(func=cObjective$run, x=x, method="simple") #"Richardson" simple
  }
  
  
  optC_disp<-trust.optim(fn = cObjective$run, gr= grad.fun,  x = opt_inits2,
                         control = list(function.scale.factor = -1, maxit=500),
                         method = c("BFGS")) #hs= #only needed for sparse #other method: "SR1" "BFGS"
  
  
  #set dispersion parameters 
  theta_opt_gamma2<- as.numeric(exp(optC_disp$solution[(n_phi+1):n_par2]))
  phi_opt_gamma2<-as.numeric(exp(optC_disp$solution[1:n_phi]))%>%matrix(nrow=S, ncol=D, byrow=TRUE)
  
  covid_inits$log_theta<-as.numeric(optC_disp$solution[(n_phi+1):n_par2])
  covid_inits$log_phi<-log(phi_opt_gamma2)#
  covid_model_cluster$setInits(list(log_theta=as.numeric(optC_disp$solution[(n_phi+1):n_par2]), log_phi=log(phi_opt_gamma2)))
  
  #run MCMC
  covid_mcmc_cluster <- buildMCMC(covid_mcmc_config_cluster)
  # Compile the MCMC.
  covid_compiled_mcmc_cluster <- compileNimble(covid_mcmc_cluster,project=covid_model_cluster)
  
  
  # Run the model 
  
  covid_samples_cluster2 <- runMCMC(covid_compiled_mcmc_cluster,niter=90000,nburnin=20000,inits=covid_inits,nchains=1,samplesAsCodaMCMC = TRUE,thin=10)
  
  return(list(covid_samples_cluster2=covid_samples_cluster2, covid_inits=covid_inits))
  
#  return(covid_samples_cluster2)
  
}

time_optimcluster_gam2<-system.time({
  covid_samples_output2  <- parLapply(cl = this_cluster_gam, X = 1:n_chains, 
                                      fun = run_MCMC_allcode, 
                                      covid_constants = covid_constants, covid_data=covid_data, covid_inits=covid_inits_gam)
})

stopCluster(this_cluster_gam)

covid_samples<-lapply(covid_samples_output2, getElement, 'covid_samples_cluster2')
covid_inits_MCMC<-lapply(covid_samples_output2, getElement, 'covid_inits')

covid_samples_gam<-as.mcmc.list(covid_samples)

#RUN time
time_optimcluster2
time_optimcluster_gam2
time_optimcluster2[3]+time_optimcluster_gam2[3]

# Combine all MCMC chains.
covid_combined_gam<- as_tibble(do.call('rbind',covid_samples_gam))
# Effective sample size
gam_ESS2<-effectiveSize(covid_samples_gam)

## Check MCMC convergence for lambda, unobserved y, and theta.
lambda_index <- which(dimnames(covid_combined_gam)[[2]]=='lambda[1, 1]'):which(dimnames(covid_combined_gam)[[2]]==paste('lambda[',N,', ',S,']',sep=''))

y_unobserved_index <- which(dimnames(covid_combined_gam)[[2]]=='y[21, 1]'):which(dimnames(covid_combined_gam)[[2]]==paste('y[',C,', ',S,']',sep=''))

#effective sample size
#ESS lambda
gam_ESS_lambda2<-effectiveSize(covid_samples_gam[,lambda_index])

#ESS unobserved y
gam_ESS_y_unob2<-effectiveSize(covid_samples_gam[,y_unobserved_index])

#ESS
summary(gam_ESS2)
summary(gam_ESS_lambda2)
summary(gam_ESS_y_unob2)

## Check convergence of lambda, unobserved y, and theta.
lambda_psrf <- sapply(lambda_index,function(x)gelman.diag(covid_samples_gam[,x],autoburnin = FALSE,transform = TRUE,multivariate = FALSE)$psrf[,1])

y_index <- which(dimnames(covid_combined_gam)[[2]]=='y[1, 1]'):which(dimnames(covid_combined_gam)[[2]]==paste('y[',C,', ',S,']',sep=''))
y_index_matrix <- matrix(y_index,nrow=C,ncol=S)
y_unob_index <- as.numeric(y_index_matrix[(C-27+2):C,])

y_unob_psrf <- sapply(y_unob_index,function(x)gelman.diag(covid_samples_gam[,x],autoburnin = FALSE,transform = TRUE,multivariate = FALSE)$psrf[,1])

theta_index <- which(dimnames(covid_combined_gam)[[2]]=='theta[1]'):which(dimnames(covid_combined_gam)[[2]]==paste('theta[',S,']',sep=''))
theta_psrf <- sapply(theta_index,function(x)gelman.diag(covid_samples_gam[,x],autoburnin = FALSE,transform = TRUE,multivariate = FALSE)$psrf[,1])

# Proportion of PSRFs below 1.05.
mean(lambda_psrf<1.05)
mean(na.omit(y_unob_psrf)<1.05)
mean(theta_psrf<1.05)

# Proportion of PSRFs below 1.1.
mean(lambda_psrf<1.10)
mean(na.omit(y_unob_psrf)<1.10)
mean(theta_psrf<1.10)

#PLOTS
# Number of MCMC samples.
n_sim <- dim(covid_combined_gam)[1]

# Overall temporal effect on covid fatality.
covid_alpha <- select(covid_combined_gam,contains('alpha'))%>%as.matrix()%>%apply(2,median)%>%melt(value.name='y')%>%mutate(x=1:N)

# Overall temporal effect on cumulative proportion reported.
covid_psi <- select(covid_combined_gam,starts_with('psi'))%>%as.matrix()%>%apply(2,median)%>%melt(value.name='y')%>%mutate(x=1:N)

# Overall temporal effect on cumulative proportion reported.
covid_iota <- select(covid_combined_gam,starts_with('iota'))%>%as.matrix()%>%apply(2,median)
cbind(region_names,covid_iota)

# Overall temporal effect on cumulative proportion reported.
covid_eta <- select(covid_combined_gam,starts_with('eta'))%>%as.matrix()%>%apply(2,median)%>%melt(value.name='y')%>%mutate(x=1:N,x=days[x])

# Regional temporal effect of covid incidence.
covid_beta <- select(covid_combined_gam,starts_with('beta'))%>%as.matrix()%>%array(dim=c(n_sim,S,D))%>%
  apply(c(2,3),quantile,c(0.025,0.5,0.975))%>%melt(varnames=c('quantile','s','d'),value.name='y')%>%
  spread(quantile,y)%>%mutate(s=factor(region_names[s],levels=region_names))

# Regional temporal effect of covid incidence.
covid_delta <- select(covid_combined_gam,starts_with('delta'))%>%as.matrix()%>%array(dim=c(n_sim,N,S))%>%
  apply(c(2,3),quantile,c(0.025,0.5,0.975))%>%melt(varnames=c('quantile','t','s'),value.name='y')%>%
  spread(quantile,y)%>%mutate(s=factor(region_names[s],levels=region_names))

# Regional temporal effect on cumulative proportion reported.
covid_gamma <- select(covid_combined_gam,starts_with('gamma'))%>%as.matrix()%>%array(dim=c(n_sim,N,S))%>%
  apply(c(2,3),quantile,c(0.025,0.5,0.975))%>%melt(varnames=c('quantile','t','s'),value.name='y')%>%
  spread(quantile,y)%>%mutate(s=factor(region_names[s],levels=region_names))

# Figure 3: Time effects plots.
delta_plot <- ggplot(covid_delta)+
  geom_line(aes(colour=as.factor(s),x=dates[t],y=`50%`))+geom_line(data=covid_alpha,aes(x=dates[x],y=y),linetype=2)+theme_minimal()+
  labs(x=NULL,y=NULL,title='Temporal Trends',subtitle='Daily Death Rate')+
  scale_color_brewer(name=NULL,palette='Dark2',guide=guide_legend(nrow=7))+
  theme(axis.text=element_text(size=10),legend.position = c(0.35,0.35),legend.text = element_text(size=10),legend.key.size = unit(0.75,'lines'),
        plot.title = element_text(size=16),plot.subtitle = element_text(size=14))
delta_plot
gamma_plot <- ggplot(covid_gamma)+
  geom_line(aes(colour=as.factor(s),x=dates[t],y=`50%`))+geom_line(data=covid_psi,aes(x=dates[x],y=y),linetype=2)+
  labs(x=NULL,y=NULL,title='',subtitle='Cumulative Proportion Reported')+guides(colour=FALSE)+theme_minimal()+
  scale_color_brewer(palette='Dark2',guide=FALSE)+
  theme(axis.text=element_text(size=10),plot.title = element_text(size=16),plot.subtitle = element_text(size=14))
effect_plots <- arrangeGrob(delta_plot,gamma_plot,nrow=1)
gamma_plot


# Figure 4: Delay effects plot.
beta_plot <- ggplot(covid_beta)+
  geom_line(aes(colour=as.factor(s),x=d,y=iprobit(`50%`)))+
  geom_point(aes(colour=as.factor(s),x=d,y=iprobit(`50%`),shape=as.factor(s)))+
  labs(x='Delay (days)',y='Cumulative Proportion Reported',title='Mean Delay Distribution')+
  theme_minimal()+
  scale_color_brewer(name=NULL,palette='Dark2')+
  scale_fill_brewer(name=NULL,palette='Dark2')+
  scale_shape_manual(name=NULL,values=1:S)+
  scale_y_continuous(trans='probit',labels=scales::label_percent())+
  theme(axis.text=element_text(size=10),legend.text = element_text(size=10),
        legend.key.size = unit(0.75,'lines'),legend.position = c(0.75,0.35),
        plot.title = element_text(size=16),plot.subtitle = element_text(size=14))
beta_plot

# Regional time of week effect on cumulative proportion reported.
covid_xi <- select(covid_combined_gam,starts_with('xi'))%>%as.matrix()%>%array(dim=c(n_sim,N,S))%>%
  apply(c(2,3),quantile,c(0.025,0.5,0.975))%>%melt(varnames=c('quantile','t','s'),value.name='y')%>%
  spread(quantile,y)%>%mutate(s=factor(region_names[s],levels=region_names),t=days[t])

covid_xi <- select(covid_combined_gam,starts_with('xi'))%>%as.matrix()%>%array(dim=c(n_sim,N,S))
covid_gamma <- select(covid_combined_gam,starts_with('gamma'))%>%as.matrix()%>%array(dim=c(n_sim,N,S))
covid_beta <- select(covid_combined_gam,starts_with('beta'))%>%as.matrix()%>%array(dim=c(n_sim,S,D))

# Negative-Binomial dispersion parameters.
covid_theta <- select(covid_combined_gam,starts_with('theta'))%>%as.matrix()

# Negative-Binomial means.
covid_lambda <- select(covid_combined_gam,starts_with('lambda'))%>%as.matrix()%>%array(dim=c(n_sim,N,S))


# Now-casting and forecasting samples of total covid deaths.
covid_y <- array(dim=c(n_sim,N,S))
covid_y[,1:C,] <- select(covid_combined_gam,starts_with('y'))%>%as.matrix()%>%array(dim=c(n_sim,C,S))
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
forecast_plot <- ggplot(filter(covid_y))+
  geom_ribbon(aes(x=dates[t],ymin=`2.5%`,ymax=`97.5%`,fill=as.factor(s)),alpha=0.25)+
  geom_ribbon(aes(x=dates[t],ymin=`10%`,ymax=`90%`,fill=as.factor(s)),alpha=0.25)+
  geom_ribbon(aes(x=dates[t],ymin=`17.5%`,ymax=`82.5%`,fill=as.factor(s)),alpha=0.25)+
  geom_ribbon(aes(x=dates[t],ymin=`25%`,ymax=`75%`,fill=as.factor(s)),alpha=0.25)+
  geom_line(aes(x=dates[t],y=`50%`,colour=as.factor(s)))+
  geom_point(data=y_data,aes(x=dates[t],y=y,colour=as.factor(s)))+
  facet_wrap(~r,nrow=2,scales='free')+
  labs(x='Date of Death',y=NULL,title='Predicted Daily Hospital Deaths from COVID-19',
       caption=NULL)+guides(colour=FALSE,fill=FALSE)+
  coord_cartesian()+geom_vline(xintercept=dates[C+1])+theme_minimal()+
  scale_fill_brewer(palette='Dark2')+scale_colour_brewer(palette='Dark2')+
  theme(strip.text = element_text(size=12),
        plot.caption = element_text(size=10),plot.title=element_text(size=16),
        axis.text = element_text(size=10),axis.title = element_text(size=12))+
  scale_x_date(date_labels = '%d %b',breaks=dates[c(5,19,33)])
forecast_plot

#trace plots
# # lambda
# plot(covid_samples_gam[,'lambda[36, 7]'], main="parallel lambda[36,7]")
# plot(covid_samples_gam[,'lambda[21, 5]'], main="parallel lambda[21,5]")
# #beta
# plot(covid_samples_gam[,'beta[1, 1]'], main="parallel beta[1,1]")
# plot(covid_samples_gam[,'beta[2, 4]'], main="parallel beta[2,4]")
# #delta
# plot(covid_samples_gam[,'delta[1, 2]'], main="parallel delta[1, 2]")
# plot(covid_samples_gam[,'delta[2, 3]'], main="parallel delta[2, 3]")
# #gamma
# plot(covid_samples_gam[,'gamma[15, 4]'], main="parallel gamma[15, 4]")
# plot(covid_samples_gam[,'gamma[15, 4]'], main="parallel gamma[34, 6]")
# #xi
# plot(covid_samples_gam[,'xi[13, 4]'], main="parallel xi[29, 7]")
# plot(covid_samples_gam[,'xi[13, 4]'], main="parallel xi[13, 4]")
# #theta
# plot(covid_samples_gam[,'theta[4]'], main="parallel theta[4]")
# plot(covid_samples_gam[,'theta[5]'], main="parallel theta[5]")
# #phi
# plot(covid_samples_gam[,'phi[2, 6]'], main="parallel delta[2, 6]")
# plot(covid_samples_gam[,'phi[7, 1]'], main="parallel delta[7,1]")
# 
# #eta
# plot(covid_samples_gam[,'eta[1]'], main="parallel eta[1]")
# plot(covid_samples_gam[,'eta[30]'], main="parallel eta[30]")
# #alpha
# plot(covid_samples_gam[,'alpha[2]'], main="parallel alpha[2]")
# plot(covid_samples_gam[,'alpha[23]'], main="parallel alpha[23]")
# #psi
# plot(covid_samples_gam[,'psi[15]'], main="parallel psi[15]")
# plot(covid_samples_gam[,'psi[35]'], main="parallel psi[35]")



lambda_NB<-matrix(NA, nrow = N, ncol=S)
y_init_region<-matrix(NA, nrow = N, ncol=S)
for(i in 1:S){
  y_init_region[,i]<-optim_output[[i]]$y_init_gamma
  lambda_NB[,i]<-optim_output[[i]]$lambda_opt
}

opt_df<-cbind(melt(y_init_region),melt(lambda_NB)[3])
colnames(opt_df)<-c("t","s","y_gamma","lambda_NB")

colours <- c("MCMC" = "orange", "NB Optimisation" = "blue", "GDM Optimisation" = "dark green")
fill.col <- c("95%" = "orange")

opt_df$r<-region_names[opt_df$s]
y_data$r<-region_names[y_data$s]
y_data<-y_data[-which(y_data$s==8),]
covid_y<-covid_y[-which(covid_y$s==8),]


region_all <- ggplot(filter(covid_y))+
  geom_ribbon(aes(x=dates[t],ymin=`2.5%`,ymax=`97.5%`,fill="95%"),alpha=0.25)+
  #geom_ribbon(aes(x=dates[t],ymin=`10%`,ymax=`90%`,fill=c("MCMC")),alpha=0.25)+
  #geom_ribbon(aes(x=dates[t],ymin=`17.5%`,ymax=`82.5%`,fill=c("MCMC")),alpha=0.25)+
  #geom_ribbon(aes(x=dates[t],ymin=`25%`,ymax=`75%`,fill=c("MCMC")),alpha=0.25)+
  geom_line(aes(x=dates[t],y=`50%`,colour=c("MCMC")))+
  #geom_line(data=filter(opt_df, s==1), aes(x=dates[t],y=`y_NB`,colour=c("Opt NB")))+
  geom_line(data=opt_df, aes(x=dates[t],y=`y_gamma`,colour=c("NB Optimisation")))+
  geom_line(data=opt_df, aes(x=dates[t],y=`lambda_NB`,colour=c("NB Optimisation")), linetype="dashed")+
  #geom_line(data=filter(opt_df, s==1), aes(x=dates[t],y=`y_cond`,colour=c("Opt conditional")))+
  geom_point(data=filter(y_data, s%in%c(1:7)),aes(x=dates[t],y=y))+
  labs(x='Date of Death',y=NULL,title='Predicted Daily Hospital Deaths from COVID-19: East of England',
       caption=NULL, colour="Model", fill="Prediction interval")+
  coord_cartesian()+geom_vline(xintercept=dates[C+1])+geom_vline(xintercept=dates[21])+
  theme_minimal()+
  facet_wrap(~r,nrow=2,scales='free')+
  scale_fill_manual(values=fill.col)+
  scale_color_manual(values=colours)+
  theme(strip.text = element_text(size=12),
        plot.caption = element_text(size=10),plot.title=element_text(size=16),
        axis.text = element_text(size=10),axis.title = element_text(size=12))+
  scale_x_date(date_labels = '%d %b',breaks=dates[c(5,19,33)])
region_all

#ggsave(region_all,filename = 'Plots/region_all_5000.png',width=13,height=5)


region1_y <- ggplot(filter(covid_y, r==region_names[1]))+
  geom_ribbon(aes(x=dates[t],ymin=`2.5%`,ymax=`97.5%`,fill="95%"),alpha=0.25)+
  #geom_ribbon(aes(x=dates[t],ymin=`10%`,ymax=`90%`,fill=c("MCMC")),alpha=0.25)+
  #geom_ribbon(aes(x=dates[t],ymin=`17.5%`,ymax=`82.5%`,fill=c("MCMC")),alpha=0.25)+
  #geom_ribbon(aes(x=dates[t],ymin=`25%`,ymax=`75%`,fill=c("MCMC")),alpha=0.25)+
  geom_line(aes(x=dates[t],y=`50%`,colour=c("MCMC")))+
  #geom_line(data=filter(opt_df, s==1), aes(x=dates[t],y=`y_NB`,colour=c("Opt NB")))+
  geom_line(data=filter(opt_df, s==1), aes(x=dates[t],y=`y_gamma`,colour=c("NB Optimisation")))+
  #geom_line(data=filter(opt_df, s==1), aes(x=dates[t],y=`y_cond`,colour=c("Opt conditional")))+
  geom_point(data=filter(y_data, r==region_names[1]),aes(x=dates[t],y=y))+
  labs(x='Date of Death',y=NULL,title='Predicted Daily Hospital Deaths from COVID-19: East of England',
       caption=NULL, colour="Model", fill="Prediction interval")+
  coord_cartesian()+geom_vline(xintercept=dates[C+1])+geom_vline(xintercept=dates[21])+
  theme_minimal()+
  scale_fill_manual(values=fill.col)+
  scale_color_manual(values=colours)+
  theme(strip.text = element_text(size=12),
        plot.caption = element_text(size=10),plot.title=element_text(size=16),
        axis.text = element_text(size=10),axis.title = element_text(size=12))+
  scale_x_date(date_labels = '%d %b',breaks=dates[c(5,19,33)])
region1_y

region2_y <- ggplot(filter(covid_y, r==region_names[2]))+
  geom_ribbon(aes(x=dates[t],ymin=`2.5%`,ymax=`97.5%`,fill=c("95%")),alpha=0.25)+
  geom_line(aes(x=dates[t],y=`50%`,colour=c("MCMC")))+
  #geom_line(data=filter(opt_df, s==2), aes(x=dates[t],y=`y_NB`,colour=c("Opt NB")))+
  geom_line(data=filter(opt_df, s==2), aes(x=dates[t],y=`y_gamma`,colour=c("NB Optimisation")))+
  #geom_line(data=filter(opt_df, s==2), aes(x=dates[t],y=`y_cond`,colour=c("Opt conditional")))+
  geom_point(data=filter(y_data, r==region_names[2]),aes(x=dates[t],y=y))+
  labs(x='Date of Death',y=NULL,title=c('Predicted Daily Hospital Deaths from COVID-19: London'),
       caption=NULL, colour="Model", fill="Prediction interval")+
  coord_cartesian()+geom_vline(xintercept=dates[C+1])+geom_vline(xintercept=dates[21])+
  theme_minimal()+
  scale_fill_manual(values=fill.col)+
  scale_color_manual(values=colours)+
  theme(strip.text = element_text(size=12),
        plot.caption = element_text(size=10),plot.title=element_text(size=16),
        axis.text = element_text(size=10),axis.title = element_text(size=12))+
  scale_x_date(date_labels = '%d %b',breaks=dates[c(5,19,33)])
region2_y

region3_y <- ggplot(filter(covid_y, r==region_names[3]))+
  geom_ribbon(aes(x=dates[t],ymin=`2.5%`,ymax=`97.5%`,fill=c("95%")),alpha=0.25)+
  geom_line(aes(x=dates[t],y=`50%`,colour=c("MCMC")))+
  #geom_line(data=filter(opt_df, s==3), aes(x=dates[t],y=`y_NB`,colour=c("Opt NB")))+
  geom_line(data=filter(opt_df, s==3), aes(x=dates[t],y=`y_gamma`,colour=c("NB Optimisation")))+
  #geom_line(data=filter(opt_df, s==3), aes(x=dates[t],y=`y_cond`,colour=c("Opt conditional")))+
  geom_point(data=filter(y_data, r==region_names[3]),aes(x=dates[t],y=y))+
  labs(x='Date of Death',y=NULL,title=c('Predicted Daily Hospital Deaths from COVID-19: Midlands'),
       caption=NULL, colour="Model", fill="Prediction interval")+
  coord_cartesian()+geom_vline(xintercept=dates[C+1])+geom_vline(xintercept=dates[21])+
  theme_minimal()+
  scale_fill_manual(values=fill.col)+
  scale_color_manual(values=colours)+
  theme(strip.text = element_text(size=12),
        plot.caption = element_text(size=10),plot.title=element_text(size=16),
        axis.text = element_text(size=10),axis.title = element_text(size=12))+
  scale_x_date(date_labels = '%d %b',breaks=dates[c(5,19,33)])
region3_y


region4_y <- ggplot(filter(covid_y, r==region_names[4]))+
  geom_ribbon(aes(x=dates[t],ymin=`2.5%`,ymax=`97.5%`,fill=c("95%")),alpha=0.25)+
  geom_line(aes(x=dates[t],y=`50%`,colour=c("MCMC")))+
  #geom_line(data=filter(opt_df, s==4), aes(x=dates[t],y=`y_NB`,colour=c("Opt NB")))+
  geom_line(data=filter(opt_df, s==4), aes(x=dates[t],y=`y_gamma`,colour=c("NB Optimisation")))+
  #geom_line(data=filter(opt_df, s==4), aes(x=dates[t],y=`y_cond`,colour=c("Opt conditional")))+
  geom_point(data=filter(y_data, r==region_names[4]),aes(x=dates[t],y=y))+
  labs(x='Date of Death',y=NULL,title=c('Predicted Daily Hospital Deaths from COVID-19: North East and Yorkshire'),
       caption=NULL, colour="Model", fill="Prediction interval")+
  coord_cartesian()+geom_vline(xintercept=dates[C+1])+geom_vline(xintercept=dates[21])+
  theme_minimal()+
  scale_fill_manual(values=fill.col)+
  scale_color_manual(values=colours)+
  theme(strip.text = element_text(size=12),
        plot.caption = element_text(size=10),plot.title=element_text(size=16),
        axis.text = element_text(size=10),axis.title = element_text(size=12))+
  scale_x_date(date_labels = '%d %b',breaks=dates[c(5,19,33)])
region4_y

region5_y <- ggplot(filter(covid_y, r==region_names[5]))+
  geom_ribbon(aes(x=dates[t],ymin=`2.5%`,ymax=`97.5%`,fill=c("95%")),alpha=0.25)+
  geom_line(aes(x=dates[t],y=`50%`,colour=c("MCMC")))+
  #geom_line(data=filter(opt_df, s==5), aes(x=dates[t],y=`y_NB`,colour=c("Opt NB")))+
  geom_line(data=filter(opt_df, s==5), aes(x=dates[t],y=`y_gamma`,colour=c("NB Optimisation")))+
  #geom_line(data=filter(opt_df, s==5), aes(x=dates[t],y=`y_cond`,colour=c("Opt conditional")))+
  geom_point(data=filter(y_data, r==region_names[5]),aes(x=dates[t],y=y))+
  labs(x='Date of Death',y=NULL,title=c('Predicted Daily Hospital Deaths from COVID-19: North West'),
       caption=NULL, colour="Model", fill="Prediction interval")+
  coord_cartesian()+geom_vline(xintercept=dates[C+1])+geom_vline(xintercept=dates[21])+
  theme_minimal()+
  scale_fill_manual(values=fill.col)+
  scale_color_manual(values=colours)+
  theme(strip.text = element_text(size=12),
        plot.caption = element_text(size=10),plot.title=element_text(size=16),
        axis.text = element_text(size=10),axis.title = element_text(size=12))+
  scale_x_date(date_labels = '%d %b',breaks=dates[c(5,19,33)])
region5_y

region6_y <- ggplot(filter(covid_y, r==region_names[6]))+
  geom_ribbon(aes(x=dates[t],ymin=`2.5%`,ymax=`97.5%`,fill=c("95%")),alpha=0.25)+
  geom_line(aes(x=dates[t],y=`50%`,colour=c("MCMC")))+
  #geom_line(data=filter(opt_df, s==6), aes(x=dates[t],y=`y_NB`,colour=c("Opt NB")))+
  geom_line(data=filter(opt_df, s==6), aes(x=dates[t],y=`y_gamma`,colour=c("NB Optimisation")))+
  #geom_line(data=filter(opt_df, s==6), aes(x=dates[t],y=`y_cond`,colour=c("Opt conditional")))+
  geom_point(data=filter(y_data, r==region_names[6]),aes(x=dates[t],y=y))+
  labs(x='Date of Death',y=NULL,title=c('Predicted Daily Hospital Deaths from COVID-19: South East'),
       caption=NULL, colour="Model", fill="Prediction interval")+
  coord_cartesian()+geom_vline(xintercept=dates[C+1])+geom_vline(xintercept=dates[21])+
  theme_minimal()+
  scale_fill_manual(values=fill.col)+
  scale_color_manual(values=colours)+
  theme(strip.text = element_text(size=12),
        plot.caption = element_text(size=10),plot.title=element_text(size=16),
        axis.text = element_text(size=10),axis.title = element_text(size=12))+
  scale_x_date(date_labels = '%d %b',breaks=dates[c(5,19,33)])
region6_y

region7_y <- ggplot(filter(covid_y, r==region_names[7]))+
  geom_ribbon(aes(x=dates[t],ymin=`2.5%`,ymax=`97.5%`,fill=c("95%")),alpha=0.25)+
  geom_line(aes(x=dates[t],y=`50%`,colour=c("MCMC")))+
  #geom_line(data=filter(opt_df, s==7), aes(x=dates[t],y=`y_NB`,colour=c("Opt NB")))+
  geom_line(data=filter(opt_df, s==7), aes(x=dates[t],y=`y_gamma`,colour=c("NB Optimisation")))+
  #geom_line(data=filter(opt_df, s==7), aes(x=dates[t],y=`y_cond`,colour=c("Opt conditional")))+
  geom_point(data=filter(y_data, r==region_names[7]),aes(x=dates[t],y=y))+
  labs(x='Date of Death',y=NULL,title=c('Predicted Daily Hospital Deaths from COVID-19: South West'),
       caption=NULL, colour="Model", fill="Prediction interval")+
  coord_cartesian()+geom_vline(xintercept=dates[C+1])+geom_vline(xintercept=dates[21])+
  theme_minimal()+
  scale_fill_manual(values=fill.col)+
  scale_color_manual(values=colours)+
  theme(strip.text = element_text(size=12),
        plot.caption = element_text(size=10),plot.title=element_text(size=16),
        axis.text = element_text(size=10),axis.title = element_text(size=12))+
  scale_x_date(date_labels = '%d %b',breaks=dates[c(5,19,33)])
region7_y


# THETA

disp_df<-data.frame( exp(covid_inits_gam$log_theta), exp(covid_inits_MCMC[[1]]$log_theta))
colnames(disp_df)<-c("theta_opt_gamma_og", "theta_opt_2")
disp_df$theta_median<-apply(covid_theta,2,quantile, 0.5)
disp_df$theta_025<-apply(covid_theta,2,quantile, 0.025)
disp_df$theta_975<-apply(covid_theta,2,quantile, 0.975)
disp_df$s<-c(1:S)

theta_opt<-ggplot(data=disp_df)+
  # geom_line(aes(x=c(1:7),y=`theta_opt_gamma2`,colour=c("GDM Optimisation")))+
  geom_line(aes(x=c(1:7),y=`theta_opt_gamma_og`,colour=c("NB Optimisation")))+
  geom_line(aes(x=c(1:7),y=`theta_opt_2`,colour=c("GDM Optimisation")))+
  geom_line(aes(x=c(1:7),y=`theta_median`,colour=c("MCMC")))+
  geom_ribbon(aes(x=c(1:7), ymin=`theta_025`, ymax=`theta_975`, fill=c("95%")),alpha=0.25)+
  labs(x='regions',y=NULL,title='Theta comparison',
       caption=NULL, colour="Model", fill="Prediction interval")+ theme_minimal()+
  scale_fill_manual(values=fill.col)+
  scale_color_manual(values=colours)+
  theme(strip.text = element_text(size=12),
        plot.caption = element_text(size=10),plot.title=element_text(size=16),
        axis.text = element_text(size=10),axis.title = element_text(size=12))

#ggsave(theta_opt,filename = 'Plots/theta_opt_5000.png',width=8,height=4)
theta_opt

covid_phi <- select(covid_combined_gam,starts_with('phi'))%>%as.matrix()%>%array(dim=c(n_sim,S,D))%>%
  apply(c(2,3),quantile,c(0.025,0.5,0.975))%>%melt(varnames=c('quantile','s','d'),value.name='y')%>%
  spread(quantile,y)%>%mutate(s=factor(region_names[s],levels=region_names))

phi_opt_gamma2<-matrix(exp(covid_inits_MCMC[[1]]$log_phi), nrow=S, ncol=D, byrow=TRUE)

phi_df<-data.frame(cbind(melt(phi_opt_gamma2),melt(exp(covid_inits_gam$log_phi))[3]))
colnames(phi_df)<-c("r","d","phi_opt_gamma","phi_opt_gamma_og")
phi_df<-phi_df%>%mutate(r=factor(r))
phi_mcmc<-covid_phi
colnames(phi_mcmc)<-c("r","d","phi_025","phi_median","phi_975")
phi_df<-full_join(phi_df,phi_mcmc,by=c("r","d"))%>%mutate(r=factor(r))
phi_df$r<-region_names[phi_df$r]

phi_opt<-ggplot(data=filter(phi_df))+ 
  geom_line(aes(x=d,y=`phi_opt_gamma`,colour=c("GDM Optimisation")))+
  geom_line(aes(x=d,y=`phi_opt_gamma_og`,colour=c("NB Optimisation")))+
  geom_line(aes(x=d,y=`phi_median`,colour=c("MCMC")))+
  geom_ribbon(aes(x=d, ymin=`phi_025`, ymax=`phi_975`, fill=c("95%")),alpha=0.25)+
  labs(x='delays',y=NULL,title='Phi comparison',
       caption=NULL, colour="Model", fill="Confidence Interval")+ theme_minimal()+
  scale_fill_manual(values=fill.col)+
  scale_color_manual(values=colours)+
  facet_wrap(~r,nrow=2,scales='free')+
  theme(strip.text = element_text(size=12),
        plot.caption = element_text(size=10),plot.title=element_text(size=16),
        axis.text = element_text(size=10),axis.title = element_text(size=12))
#ggsave(phi_opt,filename = 'Plots/phi_opt_5000.png',width=12,height=4)
phi_opt

