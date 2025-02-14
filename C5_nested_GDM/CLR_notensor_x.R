## CLR NO TENSOR MODEL ##
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
library(readr)
library(dplyr)           
library(doParallel)
library(lubridate)
library(formattable)
library(viridis)

# Set working directory.
setwd("~/JASA_code")
# Set seed.
set.seed(60209)

# Define beta binomial distribution.
  rbetabin=nimbleFunction(run=function(n=integer(0),mu=double(0),phi=double(0),size=double(0)){
    pi=rbeta(n,mu*phi,(1-mu)*phi)
    returnType(double(1))
    return(rbinom(n,size,pi))
  })
  
  
#### Code to create data matrices from raw downloaded data ####
  # #Read in data..
  # # for 2021. 
  # INFLUD21 <- read_delim("INFLUD21-30-01-2023.csv", delim = ";", escape_double = FALSE, trim_ws = TRUE)
  # # for 2022. 
  # INFLUD22 <- read_delim("INFLUD22-30-01-2023.csv", 
  #                        delim = ";", escape_double = FALSE, 
  #                        col_types = cols(DT_NOTIFIC = col_date(format = "%d/%m/%Y"), 
  #                                         DT_SIN_PRI = col_date(format = "%d/%m/%Y")), 
  #                        trim_ws = TRUE)
  # # Select attributes of interest.
  # col_ant<-which(colnames(INFLUD21)=="AN_SARS2")
  # col_pcr<-which(colnames(INFLUD21)=="PCR_SARS2")
  # col_results<-which(colnames(INFLUD21)=="PCR_RESUL")
  # col_hospital<-which(colnames(INFLUD21)=="DT_INTERNA")
  # col_collection<-which(colnames(INFLUD21)=="DT_COLETA")
  # col_an_result<-which(colnames(INFLUD21)=="DT_RES_AN")
  # col_pcr_result<-which(colnames(INFLUD21)=="DT_PCR")
  # col_notif_date<-which(colnames(INFLUD21)=="DT_DIGITA")
  # 
  # # Create data frame for 2021 SARI cases.
  # SARI_2021<- INFLUD21[,c(1:16, col_hospital, col_collection, col_results, col_pcr_result, col_notif_date, col_pcr, col_an_result, col_ant)]%>%
  #   mutate(DT_NOTIFIC=as.Date((DT_NOTIFIC),format="%d/%m/%Y"),DT_DIGITA=as.Date((DT_DIGITA),format="%d/%m/%Y"),DT_SIN_PRI=as.Date((DT_SIN_PRI),format="%d/%m/%Y"))
  # dates_2021<-sort(as.Date(unique(SARI_2021$DT_SIN_PRI),format="%d/%m/%Y"))
  # # Create data frame for 2022 SARI cases.
  # SARI_2022<- INFLUD22[,c(1:16, col_hospital, col_collection, col_results, col_pcr_result, col_notif_date, col_pcr, col_an_result, col_ant)]%>%
  #   mutate(DT_NOTIFIC=as.Date((DT_NOTIFIC),format="%d/%m/%Y"),DT_DIGITA=as.Date((DT_DIGITA),format="%d/%m/%Y"),DT_SIN_PRI=as.Date((DT_SIN_PRI),format="%d/%m/%Y"))
  # dates_2022<-sort(as.Date(unique(SARI_2022$DT_SIN_PRI),format="%d/%m/%Y"))
  # # Combine 2021 and 2022 SARI cases. 
  # SARI_2122<-rbind(SARI_2021,mutate(SARI_2022, SEM_PRI=SEM_PRI+52,SEM_NOT=SEM_NOT+52))%>%
  #   mutate(delay=difftime(DT_DIGITA,DT_SIN_PRI, unit=c("days")),delay_noti=difftime(DT_NOTIFIC,DT_SIN_PRI, unit=c("days")))%>%
  #   mutate(DT_COLETA=as.Date((DT_COLETA),format="%d/%m/%Y"),DT_PCR=as.Date((DT_PCR),format="%d/%m/%Y"))%>%
  #   mutate(pcr_delay=difftime(DT_PCR,DT_COLETA, unit=c("days")))%>%mutate(delay_diffweeks=as.numeric(difftime(DT_DIGITA,DT_SIN_PRI, unit=c("weeks"))))
  # 
  # # Delay of zero happened in first week and one happened in second week.
  # SARI_2122$delay_week<-ceiling(SARI_2122$delay_diffweeks+1/7)
  # 
  # # Make matrix with onset week as rows and delay weeks as columns.
  # SARI_delay_long<-filter(SARI_2122, delay_week>0)%>%group_by(SEM_PRI,delay_week,SG_UF_NOT)%>%summarise(cases=n())
  # SARI_delay<-SARI_delay_long%>%ungroup()%>%
  #   pivot_wider(id_cols=c(SEM_PRI,SG_UF_NOT),names_from = delay_week,values_from = cases)
  # # Any missing entries mean no (zero) cases were reported.
  # SARI_delay[is.na(SARI_delay)]<-0
  # 
  # # Latest week we have observed cases.
  # N_raw<-max(SARI_delay$SEM_PRI)
  # # Remove latest week as only 3 regions have observations.
  # SARI_delay<-filter(SARI_delay, SEM_PRI<N_raw)
  # # Number of weeks to model.
  # N<-max(SARI_delay$SEM_PRI)
  # # Maximum delay week observed.
  # D_N<-dim(SARI_delay)[2]-2
  # 
  # # Order regions and weeks in data.
  # SARI_delay_ordered<-as.matrix(SARI_delay)
  # SARI_delay_ordered<-SARI_delay_ordered[order(SARI_delay_ordered[,2]),]
  # SARI_delay_ordered<-SARI_delay_ordered[order(SARI_delay_ordered[,1]),]
  # 
  # # Save spatial constants (federal units).
  # federal_units<-sort(unique(SARI_delay_ordered[,"SG_UF_NOT"]))
  # # Total number of federal units. 
  # S<-length(federal_units)
  # 
  # # Make an array including delay:
  # # SARI_array[REGION, TIME, DELAY]
  # SARI_array<-as.numeric(SARI_delay_ordered[,-(1:2)])%>%as.matrix()%>%array(dim=c(S,N,D_N))
  # 
  # # Total SARI cases for each region and week.
  # total_SARI<-apply(SARI_array,c(1,2),sum)
  # SARI_delay_totals<-melt(total_SARI)
  # colnames(SARI_delay_totals)<-c('SG_UF_NOT',"SARI","totals")
  # 
  # 
  # # Set D_max - maximum delays we will consider. 
  # D_max<-20
  # # Set D - number of delays we explicitly model. 
  # D<-8
  # 
  # # Array for partial SARI counts by region, week and delay.
  # SARI_z<-SARI_array[,,1:(D_max)]
  # # Matrix for total SARI counts for region and week.
  # SARI_y<-apply(SARI_z,c(1,2),sum)
  # 
  # 
  # # Create data frame for COVID cases:
  # # for 2021. 
  # SARI_2021<- INFLUD21[,c(1:16,col_results,col_pcr, col_ant)]
  # dates_2021<-sort(as.Date(unique(SARI_2021$DT_SIN_PRI),format="%d/%m/%Y"))
  # 
  # SARI_counts_21<-SARI_2021%>%group_by(SEM_PRI,SG_UF_NOT)%>%
  #   summarise(SARI=n(), PCR=sum(PCR_SARS2 == 1,na.rm = TRUE ),
  #             AN=sum(AN_SARS2 == 1,na.rm = TRUE ),
  #             SARS2=sum((PCR_SARS2==1)|(AN_SARS2==1),na.rm=TRUE))%>%
  #   mutate(nu=SARS2/SARI)
  # SARI_counts_21<-SARI_counts_21%>%mutate( WEEK=SEM_PRI)
  # 
  # # for 2022.
  # SARI_2022<- INFLUD22[,c(1:16,col_results,col_pcr, col_ant)]
  # dates_2022<-sort(as.Date(unique(SARI_2022$DT_SIN_PRI),format="%d/%m/%Y"))
  # 
  # SARI_counts_22<-SARI_2022%>%group_by(SEM_PRI,SG_UF_NOT)%>%summarise(SARI=n(), PCR=sum(PCR_SARS2 == 1,na.rm = TRUE ),
  #                                                                     AN=sum(AN_SARS2 == 1,na.rm = TRUE ),
  #                                                                     SARS2=sum((PCR_SARS2==1)|(AN_SARS2==1),na.rm = TRUE))%>%
  #   mutate(nu=SARS2/SARI)
  # SARI_counts_22<-SARI_counts_22%>%mutate( WEEK=SEM_PRI+52)
  # # Combine 2021 and 2022 COVID cases.
  # SARI_counts<-rbind(SARI_counts_21,SARI_counts_22)
  # # Earliest date we have cases reported for. 
  # first_date<-min(dates_2021)
  # 
  # # Create a martix for COVID cases.
  # SARI_counts_wide <- SARI_counts %>% 
  #   select(SG_UF_NOT,SARS2,WEEK) %>%
  #   pivot_wider(names_from=SG_UF_NOT, values_from=SARS2) 
  # # Matrix of COVID counts by week and region.
  # SARI_x<-SARI_counts_wide[-N_raw,-(1:2)]
  # 
  # 
  #### Load Pre-Saved Data Matrices ####
  # Commented code above generates following data matrices from raw data files:
  # Load matrix of SARI totals.
  load("~/SARI_age_y.RData")
  # Load matrix of partial SARI counts.
  load("~/SARI_age_z.RData")
  # Load matrix of COVID-positive counts.
  load("~/SARI_age_x.RData")
  
  SARI_x<-t(apply(SARI_age_x,c(1,2),sum))
  SARI_z<-apply(SARI_age_z,c(1,2,3),sum)
  SARI_y<-t(apply(SARI_age_y,c(1,2),sum))
  
  # Set D_max - maximum delays we will consider. 
  D_max<-20
  # Set D - number of delays we explicitly model. 
  D<-8
  # Set moving window size.
  W=60  
  # Set nowcast dates for rolling prediction experiment. 
  Nowcast_list<-floor(seq(from=W, to=N-D_max, length=16))
  Nowcast_dates<-as.Date(first_date)+(Nowcast_list-1)*7
  # Set maximum covid lab report delay. 
  R<-30 
  # MCMC parameters
  n_chains=2 # number of chains 
  n_iter=15000 # number of iterations 
  n_burn=10000 # amount of burn-in
  n_thin=5 # chain thinning 
  # Set number of knots for splines 
  n_knots<-c(15,6,5) # knots for time, time-delay interaction and delay respectively 
  
  # Data frames:
  # for COVID cases. 
  x_all<-melt(as.matrix(SARI_x))
  colnames(x_all)<-c("t","s","x")
  # for SARI cases.
  y_all<-melt(t(SARI_y))
  colnames(y_all)<-c("t","s","y")
  y_all$s<-federal_units[y_all$s]
  
  # Simulate observed COVID (x) values for each nowcasting date.
  censored_x<-list()
  for( i in 1:length(Nowcast_list)){
    censored_x[[i]]<-SARI_x[(1+Nowcast_list[i]-W):Nowcast_list[i],1:S]
    censored_x[[i]][(W-R+1):W,]<-floor(censored_x[[i]][(W-R+1):W,]*seq(from=0.95, to=0.25, length=R))
    
  }
  
  # Create function to fit MCMC model using cluster:
  SARI_CLR_region <- function(seed, SARI_x, SARI_z, Nowcast_list, W, S, R, D, D_max, n_knots, n_iter, n_burn, n_thin, n_chains){
    # Load libraries.
    library(nimble)
    library(tidyverse)
    library(mgcv)
    
    # Define the Beta-Binomial as a distribution for NIMBLE.
    dbetabin=nimbleFunction(run=function(x=double(0),mu=double(0),phi=double(0),size=double(0),log=integer(0)){
      phi <- min(phi,1e+04) # Hard upper limit on phi for computational stability.
      returnType(double(0))
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
    
    # Register the Beta-Binomial as a distribution.
    registerDistributions(list(dbetabin=list(
      BUGSdist='dbetabin(mu,phi,size)',discrete=TRUE)))
    
    # NIMBLE model code. 
    SARI_code <- nimbleCode({
      for(s in 1:S){
        for(t in 1:N){
          # Beta-Binomial model for COVID cases given SARI cases.
          x[t,s] ~ dbetabin(mu[t,s],chi[s],y[s,t])
          # Reporting rate for COVID.
          pi[t,s] <- exp(log_pi[t,s])
          
        }
        for(t in 1:M){
          log_pi[t,s] <- 0
        }
        for(t in (M+1):N){
          log_pi[t,s] <- (-omega[s])*(t-M)      
          # Beta-Binomial model for censored COVID cases given eventual COVID cases.
          c[t,s] ~ dbetabin(pi[t,s],upsilon[s],x[t,s])
          
        }
        # Slope coefficient for reporting rate of COVID cases.
        omega[s] ~ dgamma(shape=10,rate=200)
        for(t in 1:N){
          # Negative Binomial Model for total SARI cases
          y[s,t] ~ dnegbin(prob=theta[s]/(theta[s]+lambda[t,s]),size=theta[s])  
          
          # Model for delayed SARI counts
          z[s,t,1] ~ dbetabin(nu[t,1,s],phi[s,1],y[s,t])
        }
        for(d in 2:D){
          for(t in 1:obs_index[s,d]){
            z[s,t,d] ~ dbetabin(nu[t,d,s],phi[s,d],y[s,t]-sum(z[s,t,1:(d-1)]))
          }
        }
        
        for(t in 1:N){
          # Linear predictor (Beta-Binomial means).
          for(d in 1:D){
            u[t,d,s]<-  eta[t,d,s] + psi[d,s] 
          }
          u[t,D+1,s]<- -sum(u[t,1:D,s])
          # Relative proportions (Beta-Binomial means).
          for(d in 1:D){
            nu[t,d,s] <- exp(u[t,d,s])/(sum(exp(u[t,d:(D+1),s])))
          }
          # Mean for total SARI cases 
          log(lambda[t,s]) <-  zeta0[s] + zeta1[t,s] 
          # Linear predictor for Binomial probability.
          logit(mu[t,s]) <- beta0[s] + beta1[t,s] + zeta1[t,s]*delta[s] 
          
        }
        
        
        
        # Splines of time in cumulative proportion reported.
        for(d in 1:D){
          Omega_eta[1:K_t2,1:K_t2,d,s] <- S_t2[1:K_t2,1:K_t2]/sigma_eta[d,s]^2
          kappa_eta[1:K_t2,d,s] ~ dmnorm(zeros[1:K_t2],Omega_eta[1:K_t2,1:K_t2,d,s])
          eta[1:N,d,s] <- X_t2[1:N,1:K_t2]%*%kappa_eta[1:K_t2,d,s]
        }
        
        # Splines of time for proportion of COVID-19/SARI cases.
        Omega_beta[1:K_t3,1:K_t3,s] <- S_t3[1:K_t3,1:K_t3]/sigma_beta[s]^2
        kappa_beta[1:K_t3,s] ~ dmnorm(zeros[1:K_t3],Omega_beta[1:K_t3,1:K_t3,s])
        beta1[1:N,s] <- X_t3[1:N,1:K_t3]%*%kappa_beta[1:K_t3,s]
        
        # Splines of time for SARI cases.
        Omega_zeta[1:K_t2,1:K_t2,s] <- S_t2[1:K_t2,1:K_t2]/sigma_zeta[s]^2
        kappa_zeta[1:K_t2,s] ~ dmnorm(zeros[1:K_t2],Omega_zeta[1:K_t2,1:K_t2,s])
        zeta1[1:N,s] <- X_t2[1:N,1:K_t2]%*%kappa_zeta[1:K_t2,s]
        
        # Smoothing parameter priors.
        sigma_beta[s] ~ T(dnorm(0,1),0,)
        sigma_zeta[s] ~ T(dnorm(0,1),0,)
      
        
        # Beta-Binomial dispersion parameters.
        for(d in 1:D){
          phi[s,d] ~ dgamma(2,0.02) 
          # Smoothing parameter priors.
          sigma_eta[d,s] ~ T(dnorm(0,1),0,) 
          # Independent Delay effect. 
          psi[d,s] ~ dnorm(logit(1/(D+2-d)),sd=2)
        }

        # Intercepts
        beta0[s] ~ dnorm(0,sd=5)
        zeta0[s]  ~ dnorm(mean.sari[s],sd=10)
        
        delta[s] ~ dnorm(0,sd=5) # Effect scale of SARI cases on proportion of COVID-19/SARI.
        chi[s] ~ dgamma(2,0.02) # Beta-Binomial dispersion parameters.
        upsilon[s] ~ dgamma(2,0.02) # Beta-Binomial dispersion parameters.
        theta[s] ~ dgamma(2,0.02) # Negative-Binomial dispersion parameters.
      }
      
    })
    
    # Constants for current nowcast date.
    N_now<-Nowcast_list[seed]
    N<-W
    C<-N-D_max+1
    M<-N-R
    
    # Censoring of partial SARI counts (z).
    z<- SARI_z[1:S,(1+N_now-W):(N_now),]
    censored_z <- z
    for(s in 1:S){
      censored_z[s,,][outer(1:dim(z[s,,])[1], 0:(dim(z[s,,])[2]-1), FUN = "+") > N] <- NA
    }
    
    # Maximum observed index for each region at each delay.
    obs_index<-matrix(NA, nrow=S, ncol=D)
    for(s in 1:S){
      for(d in 1:D){
        obs_index[s,d]<-which(is.na(censored_z[s,,d])==TRUE)[1]-1
      }
    }
    
    
    # Censor the total SARI counts (y).
    censored_y <- apply(censored_z[,,1:D_max],c(1,2),sum)
    
    # Simulate COVID (x) censoring.
    censored_x<-SARI_x[(1+N_now-W):(N_now),1:S]
    censored_x[(N-R+1):N,]<-floor(censored_x[(N-R+1):N,]*seq(from=0.95, to=0.25, length=R))
    
    
    # Set up spline bases.
    # Cubic spline with shrinkage temporal spline for expected mean SARI cases and absolute proportions reported. 
    blank_data2=tibble(y=rnorm(N,0,1),t=1:N) 
    blank_jagam2=jagam(y~s(t,k=n_knots[1], bs='cs',pc=round(C/2)),data=blank_data2,file='blank.jags',
                       knots=list(t=seq(1,N,length=n_knots[1]))) # #
    # Cubic spline with shrinkage temporal spline for proportion of SARI cases that are COVID-positive.
    blank_data3=tibble(y=rnorm(N,0,1),t=1:N) 
    blank_jagam3=jagam(y~s(t,k=n_knots[3], bs='cs', pc=round(M/2)),data=blank_data3,file='blank.jags',
                       knots=list(t=seq(1,N,length=n_knots[3]))) 
    # Constants.
    SARI_constants <- list(N=N,
                           S=S,
                           M=M,
                           D=D,
                           mean.sari=floor(log(apply(censored_y[,1:(N-D_max)],1,mean))),
                           obs_index=(obs_index),
                           K_t3=dim(blank_jagam3$jags.data$S1)[1],
                           S_t3=blank_jagam3$jags.data$S1,
                           K_t2=dim(blank_jagam2$jags.data$S1)[1],
                           S_t2=blank_jagam2$jags.data$S1,
                           x_cen=as.matrix(censored_x[(M+1):N,]))
  
    SARI_constants$X_t3 <- blank_jagam3$jags.data$X[,2:(SARI_constants$K_t3+1)]
    SARI_constants$X_t2 <- blank_jagam2$jags.data$X[,2:(SARI_constants$K_t2+1)]
    SARI_constants$zeros <- rep(0,max(SARI_constants$K_t3,SARI_constants$K_td,SARI_constants$K_t2)) #,SARI_constants$K_cases
    
    # Data.
    SARI_data <- list(y=censored_y[,1:N],
                      x=as.matrix(censored_x[1:N,],ncol=S),
                      c=as.matrix(censored_x[1:N,],ncol=S),
                      z=censored_z[,,1:D]) 
    SARI_data$x[(M+1):N,] <- NA
    
    # Generate random initial values.
    SARI_inits <- list(
      kappa_beta=matrix(rnorm(S*(SARI_constants$K_t3),0,0.1),ncol=S),
      kappa_zeta=matrix(rnorm(S*(SARI_constants$K_t2),0,0.1),ncol=S),
      kappa_eta=array(rnorm(D*S*(SARI_constants$K_t2),0,0.1),dim=c(SARI_constants$K_t2,D,S)),
      psi=matrix(rnorm(S*D,0,1),ncol=S),
      phi=matrix(abs(rnorm(S*D,30,10)),nrow=S,ncol=D),
      sigma_beta=runif(S,0,1),
      sigma_zeta=runif(S,0,1),
      delta=rnorm(S,0,1),
      chi=runif(S,5,10),
      theta=abs(rnorm(S,50,10)),
      beta0=rnorm(S,0,1),
      zeta0=rnorm(S,SARI_constants$mean.sari,1),
      sigma_eta=matrix(runif(D*S,0,1),ncol=S),
      omega=runif(S,0,0.1),
      y=matrix(NA,nrow=S,ncol=N),
      upsilon=runif(S,5,10),
      x=matrix(NA,nrow = N, ncol=S))
    # Set initial values for y which are greater than the as-of-yet reported counts:
    for(s in 1:S){
      for(t in 1:N){
        if(is.na(censored_y[s,t])) SARI_inits$y[s,t]=sum(SARI_data$z[s,t,],rpois(1,median(SARI_data$y[s,]-rowSums(SARI_data$z[s,,]),na.rm=T)),na.rm=TRUE)
      }
      for(t in (SARI_constants$M+1):N){
        if(is.na(SARI_data$x[t,s])) SARI_inits$x[t,s]=mean(c(SARI_data$c[t,s],SARI_inits$y[s,t],censored_y[s,t]),na.rm=T)
      }
    }
    
    
    # Build the model.
    SARI_model <- nimbleModel(SARI_code,SARI_constants,SARI_data,SARI_inits)
    # Compile the model.
    SARI_compiled_model <- compileNimble(SARI_model)
    # Set up the MCMC.
    SARI_mcmc_config <- configureMCMC(SARI_model,monitors=c("mu","lambda","zeta0","beta0","y","x",
                                                            "psi","chi","theta","phi","delta",
                                                            "sigma_eta","sigma_zeta", "sigma_beta",
                                                            'pi','beta1',"zeta1","omega", "upsilon",
                                                            "kappa_zeta","kappa_beta","kappa_eta"),
                                      useConjugacy = FALSE)
    
    SARI_mcmc_config$removeSamplers('kappa_zeta','sigma_zeta','kappa_beta','sigma_beta',
                                    'kappa_eta','sigma_eta','omega','delta','beta0','zeta0','psi')
    for(s in 1:S){
      SARI_mcmc_config$addSampler(target=c(paste('kappa_zeta[1:',(SARI_constants$K_t2),', ',s,']', sep=''),paste('sigma_zeta[',s,']',sep=''),paste('zeta0[',s,']',sep=''),
                                           paste('kappa_beta[1:',(SARI_constants$K_t3),', ',s,']', sep=''),paste('sigma_beta[',s,']',sep=''),paste('beta0[',s,']',sep=''),
                                           paste('delta[',s,']',sep=''),paste('omega[',s,']',sep='')),type='AF_slice') 
      SARI_mcmc_config$addSampler(target=c(paste('kappa_eta[1:',(SARI_constants$K_t2),', 1:',D,', ',s,']', sep=''),paste('sigma_eta[1:',D,', ',s,']',sep=''),paste('psi[1:',D,', ',s,']',sep='')),type='AF_slice')
    }
    
    SARI_mcmc<- buildMCMC(SARI_mcmc_config)
    # Compile the MCMC.
    SARI_compiled_mcmc <- compileNimble(SARI_mcmc,project=SARI_model)
    # Run the MCMC. 
    SARI_output_cluster<- runMCMC(SARI_compiled_mcmc,niter=n_iter,nburnin=n_burn,thin=n_thin,inits=SARI_inits,nchains=n_chains,samplesAsCodaMCMC = TRUE)
    
    return(SARI_output_cluster)
  }
  
  # CLR model - with no tensor product in the delay distribution. 
  # Number of available computer cores.
  n_cores<-min(detectCores(),length(Nowcast_list))
  # Make Cluster for MCMC. 
  this_cluster_region<- makeCluster(n_cores)
  
  # Run cluster. 
  time_clr_notensor<-system.time({
    SARI_clr_notensor <- parLapply(cl = this_cluster_region, X = 1:length(Nowcast_list), 
                                  fun = SARI_CLR_region,
                                  SARI_x=as.matrix(SARI_x),
                                  SARI_z=as.array(SARI_z),
                                  Nowcast_list=Nowcast_list,
                                  W=W,
                                  S=S,
                                  R=R,
                                  D=D,
                                  D_max=D_max,
                                  n_iter=n_iter,
                                  n_burn=n_burn,
                                  n_thin=n_thin,
                                  n_chains=n_chains,
                                  n_knots=n_knots)
  })
  
  
  stopCluster(this_cluster_region)
  
  # Save Samples.
  save(SARI_clr_notensor,time_clr_notensor,file='cluster_notensor_x.RData')
  #
  
  # Check PSRF of model parameters.
  sari_parameter_group=c("mu","lambda","zeta0","beta0",
                         "psi","chi","theta","phi","delta",
                         "sigma_eta","sigma_zeta", "sigma_beta",
                         'pi','beta1',"zeta1","omega",
                         "kappa_zeta","kappa_beta","kappa_eta")
  sari_psrf<-list()
  sari_mean_psrf<-sari_max_psrf<-sari_leq_105 <-sari_leq_110<-array(NA,dim=c(length(sari_parameter_group),length(Nowcast_list)))
  for(k in 1:length(sari_parameter_group)){
    sari_psrf[[k]]<-list()
    for(j in 1:length(Nowcast_list)){
      sari_parameter_names <- colnames((as.matrix(SARI_clr_notensor[[j]])))
      sari_index <- which(startsWith(sari_parameter_names,sari_parameter_group[k]))
      sari_psrf[[k]][[j]]<- sapply(sari_index,function(v)gelman.diag(SARI_clr_notensor[[j]][,v],autoburnin = FALSE,transform = TRUE,multivariate = FALSE)$psrf[,1])
      sari_mean_psrf[k,j] <- mean(sari_psrf[[k]][[j]])
      sari_max_psrf[k,j] <- max(sari_psrf[[k]][[j]])
      sari_leq_105[k,j] <- mean(sari_psrf[[k]][[j]]<=1.05)
      sari_leq_110[k,j] <- mean(sari_psrf[[k]][[j]][j]<=1.1)
    }
  }
  sari_MAX_psrf<-cbind(sari_parameter_group,apply(sari_max_psrf,1,max))
  sari_MEAN_psrf<-cbind(sari_parameter_group,apply(sari_mean_psrf,1,mean))
  