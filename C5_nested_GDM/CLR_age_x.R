## CLR AGE MODEL ##
library(ggplot2)
library(tidyverse)
library(reshape2)
library(mgcv)
library(coda)
library(abind)
library(nimble)
library(readr)
library(dplyr)           
library(doParallel)
library(lubridate)
library(formattable)
library(gt)
library(gtExtras)
library(grDevices)
library(viridis)
library(grid)
library(gridExtra)

# Set working directory.
setwd("~/JASA_code")
# Set seed.
set.seed(60209)

#### Code to create data matrices from raw downloaded data ####
# 
# #Read in data..
# # for 2021.
# INFLUD21 <- read_delim("INFLUD21-30-01-2023.csv", delim = ";", escape_double = FALSE, trim_ws = TRUE)
# # for 2022.
# INFLUD22 <- read_delim("INFLUD22-30-01-2023.csv", 
#                        delim = ";", escape_double = FALSE, 
#                        col_types = cols(DT_NOTIFIC = col_date(format = "%d/%m/%Y"), 
#                                         DT_SIN_PRI = col_date(format = "%d/%m/%Y")), 
#                        trim_ws = TRUE)
# 
# # Identify columns of interest:
# col_ant<-which(colnames(INFLUD21)=="AN_SARS2")
# col_pcr<-which(colnames(INFLUD21)=="PCR_SARS2")
# col_results<-which(colnames(INFLUD21)=="PCR_RESUL")
# col_hospital<-which(colnames(INFLUD21)=="DT_INTERNA")
# col_collection<-which(colnames(INFLUD21)=="DT_COLETA")
# col_an_result<-which(colnames(INFLUD21)=="DT_RES_AN")
# col_pcr_result<-which(colnames(INFLUD21)=="DT_PCR")
# col_notif_date<-which(colnames(INFLUD21)=="DT_DIGITA")
# 
# # Create data frame for 2021.
# SARI_2021<- INFLUD21[,c(1:16, col_hospital, col_collection, col_results, col_pcr_result, col_notif_date, col_pcr, col_an_result, col_ant)]%>%
#   mutate(DT_NOTIFIC=as.Date((DT_NOTIFIC),format="%d/%m/%Y"),DT_DIGITA=as.Date((DT_DIGITA),format="%d/%m/%Y"),DT_SIN_PRI=as.Date((DT_SIN_PRI),format="%d/%m/%Y"))
# dates_2021<-sort(as.Date(unique(SARI_2021$DT_SIN_PRI),format="%d/%m/%Y"))
# # Create data frame for 2022
# SARI_2022<- INFLUD22[,c(1:16, col_hospital, col_collection, col_results, col_pcr_result, col_notif_date, col_pcr, col_an_result, col_ant)]%>%
#   mutate(DT_NOTIFIC=as.Date((DT_NOTIFIC),format="%d/%m/%Y"),DT_DIGITA=as.Date((DT_DIGITA),format="%d/%m/%Y"),DT_SIN_PRI=as.Date((DT_SIN_PRI),format="%d/%m/%Y"))%>%
#   filter(DT_SIN_PRI<as.Date("2023-01-01")) # cut off at end of 2022 to ensure fully reported data
# dates_2022<-sort(as.Date(unique(SARI_2022$DT_SIN_PRI),format="%d/%m/%Y"))
# # Create data frame for 2021 and 2022. 
# SARI_2122<-rbind(SARI_2021,mutate(SARI_2022, SEM_PRI=SEM_PRI+52,SEM_NOT=SEM_NOT+52))%>%
#   mutate(delay=difftime(DT_DIGITA,DT_SIN_PRI, unit=c("days")),delay_noti=difftime(DT_NOTIFIC,DT_SIN_PRI, unit=c("days")))%>%
#   mutate(DT_COLETA=as.Date((DT_COLETA),format="%d/%m/%Y"),DT_PCR=as.Date((DT_PCR),format="%d/%m/%Y"),DT_NASC=as.Date((DT_NASC),format="%d/%m/%Y"))%>%
#   mutate(pcr_delay=difftime(DT_PCR,DT_COLETA, unit=c("days")))%>%
#   mutate(delay_diffweeks=as.numeric(difftime(DT_DIGITA,DT_SIN_PRI, unit=c("weeks")))) #
# # Use difference between digitalization date (DT_DIGITA) and date of onset of symptoms (DT_SIN_PRI) for delay.
# 
# # Earliest date we have data for:
# first_date<-min(SARI_2122$DT_SIN_PRI)
# 
# # Determine delay in terms of weeks not days:
# SARI_2122$delay_week<-ceiling(SARI_2122$delay_diffweeks+1/7)
# 
# # Create three age groups.
# SARI_2122<-SARI_2122%>%mutate(age=round(as.numeric(time_length(difftime(DT_SIN_PRI,DT_NASC), "years"))))%>%
#   mutate(age_group=cut(age, c(0, 18, 60, Inf), c("0-18", "19-60", ">60"), include.lowest=TRUE))
# # COVID cases are defined by a PCR or ANTIGEN test being positive for SARS2.
# SARI_2122<-SARI_2122%>%mutate(SARS2=as.numeric((PCR_SARS2==1)|(AN_SARS2==1)))
# 
# 
# # Make matrix for each age group with onset week as rows and delay weeks as columns:
# # Filter for first age group
# SARI_0_18<-filter(SARI_2122,age_group=="0-18")
# # Make matrix with onset week as rows and delay weeks as columns
# SARI_0_18_long<-filter(SARI_0_18, delay_week>0)%>%group_by(SEM_PRI,delay_week,SG_UF_NOT)%>%summarise(cases=n())
# 
# # remove latest week as only 3 regions have observations 
# N_raw<-max(SARI_0_18_long$SEM_PRI)
# SARI_0_18_long<-filter(SARI_0_18_long, SEM_PRI<N_raw)
# 
# SARI_0_18_delay<-SARI_0_18_long%>%ungroup()%>%
#   pivot_wider(id_cols=c(SEM_PRI,SG_UF_NOT),names_from = delay_week,values_from = cases)
# SARI_0_18_delay[is.na(SARI_0_18_delay)]<-0
# 
# SARI_0_18_totals<-SARI_0_18_delay
# SARI_0_18_totals$total<-rowSums(SARI_0_18_delay[,-(1:2)], na.rm=T)
# 
# 
# # Filter for second age group
# SARI_19_60<-filter(SARI_2122,age_group=="19-60")
# # Make matrix with onset week as rows and delay weeks as columns
# SARI_19_60_long<-filter(SARI_19_60, delay_week>0)%>%group_by(SEM_PRI,delay_week,SG_UF_NOT)%>%summarise(cases=n())
# 
# # remove latest week as only 3 regions have observations 
# N_raw<-max(SARI_19_60_long$SEM_PRI)
# SARI_19_60_long<-filter(SARI_19_60_long, SEM_PRI<N_raw)
# 
# SARI_19_60_delay<-SARI_19_60_long%>%ungroup()%>%
#   pivot_wider(id_cols=c(SEM_PRI,SG_UF_NOT),names_from = delay_week,values_from = cases)
# SARI_19_60_delay[is.na(SARI_19_60_delay)]<-0
# 
# SARI_19_60_totals<-SARI_19_60_delay
# SARI_19_60_totals$total<-rowSums(SARI_19_60_delay[,-(1:2)], na.rm=T)
# 
# # Filter for third age group
# SARI_61<-filter(SARI_2122,age_group==">60")
# # Make matrix with onset week as rows and delay weeks as columns
# SARI_61_long<-filter(SARI_61, delay_week>0)%>%group_by(SEM_PRI,delay_week,SG_UF_NOT)%>%summarise(cases=n())
# 
# # remove latest week as only 3 regions have observations 
# N_raw<-max(SARI_61_long$SEM_PRI)
# SARI_61_long<-filter(SARI_61_long, SEM_PRI<N_raw)
# 
# SARI_61_delay<-SARI_61_long%>%ungroup()%>%
#   pivot_wider(id_cols=c(SEM_PRI,SG_UF_NOT),names_from = delay_week,values_from = cases)
# SARI_61_delay[is.na(SARI_61_delay)]<-0
# 
# SARI_61_totals<-SARI_61_delay
# SARI_61_totals$total<-rowSums(SARI_61_delay[,-(1:2)], na.rm=T)
# 
# 
# 
# # Spatial names (federal_units) and number of regions (S)
# federal_units<-sort(unique(SARI_19_60_delay$SG_UF_NOT))
# S<-length(federal_units)
# N<-N_raw-1
# # number of age groups
# A<-3
# # Maximum Delay:
# D_max<-20
# # Number of delays explicitly modeled in GDM:
# D<-8
# 
# C <- N-D_max+1
# 
# # Order arrays:
# SARI_0_18_delay_ordered<-SARI_0_18_delay[,1:(D_max+2)]%>%mutate(id=paste(SEM_PRI,SG_UF_NOT,sep=''))
# SARI_0_18_delay_ordered<-SARI_0_18_delay_ordered%>%arrange(SEM_PRI,SG_UF_NOT)
# SARI_19_60_delay_ordered<-SARI_19_60_delay[,1:(D_max+2)]%>%mutate(id=paste(SEM_PRI,SG_UF_NOT,sep=''))
# SARI_19_60_delay_ordered<-SARI_19_60_delay_ordered%>%arrange(SEM_PRI,SG_UF_NOT)
# SARI_61_delay_ordered<-SARI_61_delay[,1:(D_max+2)]%>%mutate(id=paste(SEM_PRI,SG_UF_NOT,sep=''))
# SARI_61_delay_ordered<-SARI_61_delay_ordered%>%arrange(SEM_PRI,SG_UF_NOT)
# 
# 
# # Create blank data frame for all regions and time. 
# SARI_0_18_full <- cbind(tibble(SEM_PRI=sort(rep(1:N,27)), SG_UF_NOT=rep(sort(federal_units),N)),
#                         matrix(0, ncol=D_max,nrow=S*N))%>%mutate(check = NA)%>%mutate(id=paste(SEM_PRI,SG_UF_NOT,sep=''))
# SARI_19_60_full <- cbind(tibble(SEM_PRI=sort(rep(1:N,27)), SG_UF_NOT=rep(sort(federal_units),N)),
#                          matrix(0, ncol=D_max,nrow=S*N))%>%mutate(check = NA)%>%mutate(id=paste(SEM_PRI,SG_UF_NOT,sep=''))
# SARI_61_full <- cbind(tibble(SEM_PRI=sort(rep(1:N,27)), SG_UF_NOT=rep(sort(federal_units),N)),
#                       matrix(0, ncol=D_max,nrow=S*N))%>%mutate(check = NA)%>%mutate(id=paste(SEM_PRI,SG_UF_NOT,sep=''))
# 
# # Fill in available data
# SARI_0_18_full[SARI_0_18_full$id %in% SARI_0_18_delay_ordered$id ,3:(D_max+3)] <- SARI_0_18_delay_ordered[,3:(D_max+3)]
# SARI_19_60_full[SARI_19_60_full$id %in% SARI_19_60_delay_ordered$id ,3:(D_max+3)] <- SARI_19_60_delay_ordered[,3:(D_max+3)]
# SARI_61_full[SARI_61_full$id %in% SARI_61_delay_ordered$id ,3:(D_max+3)] <- SARI_61_delay_ordered[,3:(D_max+3)]
# 
# 
# # Make an array including delay:
# # SARI_..._array[REGION, TIME, DELAY]
# SARI_0_18_array<-SARI_0_18_full[,3:(D_max+2)]%>%as.matrix()%>%array(dim=c(S,N,D_max))
# SARI_19_60_array<-SARI_19_60_full[,3:(D_max+2)]%>%as.matrix()%>%array(dim=c(S,N,D_max))
# SARI_61_array<-SARI_61_full[,3:(D_max+2)]%>%as.matrix()%>%array(dim=c(S,N,D_max))
# # SARI_array[REGION, TIME, DELAY, AGE GROUP]
# SARI_array<-array(abind(SARI_0_18_array,SARI_19_60_array,SARI_61_array),dim=c(S,N,D_max,A))
# # total SARI cases for each region
# total_SARI<-apply(SARI_array,c(1,4),sum)
# SARI_delay_totals<-melt(total_SARI)
# colnames(SARI_delay_totals)<-c('SG_UF_NOT',"AGE GROUP","totals")
# 
# 
# 
# SARI_age_z<-SARI_array[,,1:D_max,] #[REGION, TIME, DELAY, AGE GROUP]
# SARI_age_y<-apply(SARI_age_z,c(1,2,4), sum) #[REGION, TIME, AGE GROUP]
# SARI_dates<-as.Date(first_date)+(1:ncol(SARI_age_y)-1)*7
# 
# # Make array for COVID-19 cases. 
# COVID_age<-filter(SARI_2122, delay_week>0, delay_week<21)%>%group_by(SEM_PRI,age_group,SARS2,SG_UF_NOT)%>%summarise(cases=n())%>%drop_na(age_group,SARS2)%>%
#   pivot_wider(id_cols=c(SEM_PRI,SG_UF_NOT),names_from = age_group,values_from = cases)
# COVID_age[is.na(COVID_age)]<-0
# 
# # Remove latest week as only 3 regions have observations.
# COVID_long<-filter(COVID_age, SEM_PRI<N_raw)
# # Add ID and order matrix.
# COVID_long_ordered<-COVID_long%>%mutate(id=paste(SEM_PRI,SG_UF_NOT,sep=''))%>%arrange(SEM_PRI,SG_UF_NOT)
# 
# # Create blank data frame for all regions and time. 
# COVID_full <- tibble(SEM_PRI=sort(rep(1:N,27)), SG_UF_NOT=rep(sort(federal_units),N),
#                      `0-18` = 0, `19-60` = 0, `>60` = 0, check = NA)
# COVID_full<-COVID_full%>%mutate(id=paste(SEM_PRI,SG_UF_NOT,sep=''))
# 
# # Fill in available data
# COVID_full[COVID_full$id %in% COVID_long_ordered$id ,3:6] <- COVID_long_ordered[,3:6]
# COVID_wide<-COVID_full[,-c(1,2,6,7)]%>%as.matrix()%>%array(dim=c(S,N,A))
# SARI_age_x<-COVID_wide # [REGION, TIME, AGE]
# 
# save(SARI_age_x, file="SARI_age_x.RData")
# save(SARI_age_z, file="SARI_age_z.RData")
# save(SARI_age_y, file="SARI_age_y.RData")

#### Load Pre-Saved Data Matrices ####
# Commented code above generates following data matrices from raw data files:
# Load matrix of SARI totals.
load("~/SARI_age_y.RData")
# Load matrix of partial SARI counts.
load("~/SARI_age_z.RData")
# Load matrix of COVID-positive counts.
load("~/SARI_age_x.RData")

x_all<-melt(apply(SARI_age_x,c(1,2),sum))
colnames(x_all)<-c("s","t","x")
x_all$s<-federal_units[x_all$s]

y_all<-melt(apply(SARI_age_y,c(1,2),sum))
colnames(y_all)<-c("s","t","y")
y_all$s<-federal_units[y_all$s]


# MCMC parameters.
n_chains=2
n_iter=15000
n_burn=10000
n_thin=5
n_knots<-c(15,6,8) # knots for time, time-delay interaction and delay respectively 
R<-30 # maximum COVID-19 lab report delay 
W<-60 # moving window size
# Set D_max - maximum delays we will consider. 
D_max<-20
# Set D - number of delays we explicitly model. 
D<-8

# List of dates to preform the nowcasts for rolling predictions experiment.
Nowcast_list<-floor(seq(from=W, to=N-D_max, length=16))
# Dates for the nowcasts.
Nowcast_dates<- as.Date(first_date)+(Nowcast_list-1)*7

# Create function for cluster to carry out MCMC model:
SARI_age_reg <- function(seed, SARI_x, SARI_z, Nowcast_list, W, S, R, D, D_max, n_knots, n_iter, n_burn, n_thin, n_chains){
  # Load libraries.
  library(nimble)
  library(tidyverse)
  library(mgcv)
  
  # Define the Beta-Binomial as a distribution for NIMBLE.
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
  
  # Register the Beta-Binomial as a distribution.
  registerDistributions(list(dbetabin=list(
    BUGSdist='dbetabin(mu,phi,size)',discrete=TRUE)))
  
  # NIMBLE model code. 
  SARI_code <- nimbleCode({
    for(s in 1:S){
      for(t in 1:N){
        # Beta-Binomial model for COVID cases given SARI cases.
        x[t,s] ~ dbetabin(mu[t,s],chi[s],y[s,t])  
        # Reporting rate for COVID-19.
        pi[t,s] <- exp(log_pi[t,s])
      }
      for(t in 1:M){
        log_pi[t,s] <- 0
      }
      for(t in (M+1):N){
        log_pi[t,s] <- (-omega[s])*(t-M)
        c[t,s] ~ dbetabin(pi[t,s],upsilon[s],x[t,s])
        
      }
      # Slope coefficient for reporting rate of COVID cases.
      omega[s] ~ dgamma(shape=10,rate=200)
      for(t in 1:N){
        # Negative Binomial Model for total SARI cases
        y[s,t] ~ dnegbin(prob=theta[s]/(theta[s]+lambda[t,s]),size=theta[s])  
        
        # Model for delayed SARI counts (first delay).
        z[s,t,1] ~ dbetabin(nu[t,1,s],phi[s,1],y[s,t])
      }
      for(d in 2:D){
        for(t in 1:obs_index[s,d]){
          # Model for delayed SARI counts (from second delay).
          z[s,t,d] ~ dbetabin(nu[t,d,s],phi[s,d],y[s,t]-sum(z[s,t,1:(d-1)]))
        }
      }
      
      for(t in 1:N){
        # Linear predictor (Beta-Binomial means).
        for(d in 1:D){
          # Absolute proportions.
          u[t,d,s] <- gamma0[s] + gamma[t,d,s]
        }
        u[t,D+1,s]<- -sum(u[t,1:D,s])
        for(d in 1:D){
          # Relative proportions (Beta-Binomial means).
          nu[t,d,s] <- exp(u[t,d,s])/(sum(exp(u[t,d:(D+1),s])))
        }
        # Mean for total SARI cases. 
        log(lambda[t,s]) <-  zeta0[s] + zeta1[t,s]   
        # Linear predictor for COVID-positive Beta-Binomial probability.
        logit(mu[t,s]) <- beta0[s] + beta1[t,s] + zeta1[t,s]*delta[s] + beta_cov_0to18[1]*age_0to18[t,1,s] + beta_cov_0to18[2]*age_0to18[t,2,s] + beta_cov_0to18[3]*age_0to18[t,3,s] +  
          beta_cov_19to60[1]*age_19to60[t,1,s] + beta_cov_19to60[2]*age_19to60[t,2,s] + beta_cov_19to60[3]*age_19to60[t,3,s]
        
      }
      
      
      # Time- Delay spline
      Omega_gamma[1:K_td,1:K_td,s] <- S_td[1:K_td,1:K_td]/sigma_gamma[1,s]^2+S_td[1:K_td,(K_td+1):(2*K_td)]/sigma_gamma[2,s]^2+S_td[1:K_td,(2*K_td+1):(3*K_td)]/sigma_gamma[3,s]^2
      kappa_gamma[1:K_td,s] ~ dmnorm(zeros[1:K_td],Omega_gamma[1:K_td,1:K_td,s])
      gamma_long[1:(N*D),s] <-  X_td[1:(N*D),1:K_td]%*%kappa_gamma[1:K_td,s]
      gamma[1:N,1:D,s]<-matrix(gamma_long[1:(N*D),s],ncol=D)
      
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
      for(k in 1:3){
        sigma_gamma[k,s] ~ T(dnorm(0,1),0,)
      }
      
      # Beta-Binomial dispersion parameters.
      for(d in 1:D){
        phi[s,d] ~ dgamma(2,0.02) 
      }
      
      # Intercepts
      beta0[s] ~ dnorm(0,sd=1)
      zeta0[s]  ~ dnorm(mean.sari[s],sd=10)
      gamma0[s]  ~ dnorm(0,sd=10)
      
      upsilon[s] ~ dgamma(2,0.02) # Beta-Binomial dispersion parameters.
      delta[s] ~ dnorm(0,sd=5) # Effect scale of SARI cases on proportion of COVID-19/SARI.
      chi[s] ~ dgamma(2,0.02) # Beta-Binomial dispersion parameters.
      theta[s] ~ dgamma(2,0.02) # Negative-Binomial dispersion parameters.
    }
    
    # Age group coefficients.
    for(j in 1:3){
      beta_cov_0to18[j] ~ dnorm(0,sd=2)
      beta_cov_19to60[j] ~ dnorm(0,sd=2)
    }
    
  })
  
  # Data for current nowcast time 
  N_now<-Nowcast_list[seed]
  N<-W
  C<-N-D_max+1
  M<-N-R
  A<-dim(SARI_z)[4]
  
  # Censore the partial SARI counts (z).
  z<- SARI_z[1:S,(1+N_now-W):(N_now),,1:A]
  censored_z_age <- z
  for(a in 1:A){
    for(s in 1:S){
      a_z<-censored_z_age[s,,,a]
      a_z[outer(1:dim(z[s,,,a])[1], 0:(dim(z[s,,,a])[2]-1), FUN = "+") > N] <- NA
      censored_z_age[s,,,a]<-a_z
    }
  }
  
  # Censor the SARI totals (y).
  censored_y<- apply(censored_z_age[1:S,,1:D_max,1:A],c(1,2),sum)
  censored_y_na<- apply(censored_z_age[1:S,,1:D_max,1:A],c(1,2,4),sum,na.rm=TRUE)
  # Calculate age group covariates. 
  age_group_prop<-censored_y_na[,,1:A]/array(rep(apply(censored_y_na,c(1,2),sum),A),dim=c(S,N,A))
  age_group_prop[age_group_prop=="NaN"]<-0
  # Scale age group covariates. 
  age_group_prop_scaled <- age_group_prop%>%apply(3,function(x){(x-mean(x))/sd(x)})%>%
    array(dim=dim(age_group_prop))%>%
    aperm(c(2,1,3))
  # Covaraites for cubic polynomial:
  age_group_cubic<-array(NA, dim=c(N,3,S,2)) #time,cubic terms order,region,age group
  age_group_cubic[,1,,] <- age_group_prop_scaled[,,1:2]
  age_group_cubic[,2,,] <- age_group_prop_scaled[,,1:2]^2
  age_group_cubic[,3,,] <- age_group_prop_scaled[,,1:2]^3
  
  # Calculate censored z over all age groups.
  censored_z <- apply(censored_z_age,c(1,2,3),sum)
  
  # Maximum observed index for each region at each delay.
  obs_index<-matrix(nrow=D, ncol=S)
  for(d in 1:D){
    for(s in 1:S){
      obs_index[d,s]<-which(is.na(censored_z[s,,d])==TRUE)[1]-1
    }
  }
  
  # Simulate COVID censoring.
  censored_x<-SARI_x[1:S,(1+N_now-W):(N_now),1:A]
  for(a in 1:A){
    for(s in 1:S){
      censored_x[s,(N-R+1):N,a]<-floor(censored_x[s,(N-R+1):N,a]*seq(from=0.95, to=0.25, length=R)) 
    }
  }
  censored_x<-t(apply(censored_x,c(1,2),sum))
  
  # Set up spline bases.
  # Tensor product for absolute proportions of SARI counts reported.
  blank_data=tibble(z=rnorm(N*D,0,1), t=rep(1:N,D),
                    d=as.numeric(matrix(rep(1:D,N), ncol=D, byrow=TRUE)))
  blank_jagam=jagam(z~te(t,d, bs=c('tp','tp'), k=n_knots[1:2]), 
                    data=blank_data, file='blank.jags',knots=list(t=seq(1,N,length=n_knots[1]),d=seq(1,D,length=n_knots[2])))
  # Cubic spline with shrinkage for mean SARI cases.
  blank_data2=tibble(y=rnorm(N,0,1),t=1:N) 
  blank_jagam2=jagam(y~s(t,k=n_knots[1], bs='cs',pc=round(C/2)),data=blank_data2,file='blank.jags',
                     knots=list(t=seq(1,N,length=n_knots[1]))) # #
  # Cubic spline with shrinkage for proportion of SARI cases that are COVID-positive.
  blank_data3=tibble(y=rnorm(N,0,1),t=1:N) 
  blank_jagam3=jagam(y~s(t,k=n_knots[3], bs='cs', pc=round(M/2)),data=blank_data3,file='blank.jags',
                     knots=list(t=seq(1,N,length=n_knots[3]))) 
  # Constants.
  SARI_constants <- list(N=N,
                         S=S,
                         M=M,
                         D=D,
                         age_0to18=array(age_group_cubic[,,,1],dim=c(N,3,S)),
                         age_19to60=array(age_group_cubic[,,,2],dim=c(N,3,S)),
                         mean.sari=floor(log(apply(censored_y[,1:(N-D_max)],1,mean))),
                         obs_index=t(obs_index),
                         K_t3=dim(blank_jagam3$jags.data$S1)[1],
                         S_t3=blank_jagam3$jags.data$S1,
                         K_td=dim(blank_jagam$jags.data$S1)[1],
                         S_td=blank_jagam$jags.data$S1,
                         K_t2=dim(blank_jagam2$jags.data$S1)[1],
                         S_t2=blank_jagam2$jags.data$S1)
  SARI_constants$X_t3 <- blank_jagam3$jags.data$X[,2:(SARI_constants$K_t3+1)]
  SARI_constants$X_td <- blank_jagam$jags.data$X[,2:(SARI_constants$K_td+1)]
  SARI_constants$X_t2 <- blank_jagam2$jags.data$X[,2:(SARI_constants$K_t2+1)]
  SARI_constants$zeros <- rep(0,max(SARI_constants$K_t3,SARI_constants$K_td,SARI_constants$K_t2)) #,SARI_constants$K_cases
  
  # Data:
  SARI_data <- list(y=censored_y[,1:N],
                    x=as.matrix(censored_x[1:N,],ncol=S),
                    c=as.matrix(censored_x[1:N,],ncol=S),
                    z=censored_z[,,1:D]) 
  SARI_data$x[(M+1):N,] <- NA
  
  
  # Generate random initial values.
  SARI_inits <- list(
    beta_cov_0to18=(rnorm(3,0,0.1)),
    beta_cov_19to60=(rnorm(3,0,0.1)),
    kappa_beta=matrix(rnorm(S*(SARI_constants$K_t3),0,0.1),ncol=S),
    kappa_zeta=matrix(rnorm(S*(SARI_constants$K_t2),0,0.1),ncol=S),
    kappa_gamma=matrix(rnorm(S*(SARI_constants$K_td),0,0.1),ncol=S),
    phi=matrix(abs(rnorm(S*D,30,10)),nrow=S,ncol=D),
    sigma_beta=runif(S,0,1),
    sigma_zeta=runif(S,0,1),
    sigma_gamma=matrix(runif(3*S,0,1),ncol=S),
    delta=rnorm(S,0,1),
    chi=runif(S,5,10),
    theta=abs(rnorm(S,50,10)),
    beta0=rnorm(S,0,1),
    gamma0=rnorm(S,0,1),
    zeta0=rnorm(S,SARI_constants$mean.sari,1),
    omega=runif(S,0,0.1),
    y=matrix(NA,nrow=S,ncol=N), 
    x=matrix(NA,nrow=N,ncol=S) 
  )
  # Set initial values for y which are greater than the as-of-yet reported counts:
  for(t in 1:N){
    for(s in 1:S){
      if(is.na(censored_y[s,t])) SARI_inits$y[s,t]=sum(SARI_data$z[s,t,],rpois(1,median(SARI_data$y[s,]-rowSums(SARI_data$z[s,,]),na.rm=T)),na.rm=TRUE)
      if(is.na(SARI_data$x[t,s])) SARI_inits$x[t,s]=mean(c(SARI_data$c[t,s],SARI_inits$y[s,t],censored_y[s,t]),na.rm=T)
    }
  }
  
  # Build the model.
  SARI_model <- nimbleModel(SARI_code,SARI_constants,SARI_data,SARI_inits)
  # Compile the model.
  SARI_compiled_model <- compileNimble(SARI_model)
  # Set up the MCMC.
  SARI_mcmc_config <- configureMCMC(SARI_model,monitors=c("mu","lambda","zeta0","beta0","y","x",
                                                          "gamma0","chi","theta","phi","delta",
                                                          "beta_cov_0to18","beta_cov_19to60",
                                                          "sigma_gamma","sigma_zeta","pi","beta1",
                                                          "sigma_beta","kappa_zeta","kappa_beta","kappa_gamma","omega"),
                                    useConjugacy = FALSE)
  
  SARI_mcmc_config$removeSamplers('kappa_zeta','sigma_zeta','kappa_beta','sigma_beta',
                                  'kappa_gamma','sigma_gamma','omega','delta', 
                                  'gamma0','beta0','zeta0')
  for(s in 1:S){
    SARI_mcmc_config$addSampler(target=c(paste('kappa_zeta[1:',(SARI_constants$K_t2),', ',s,']', sep=''),paste('sigma_zeta[',s,']',sep=''),paste('zeta0[',s,']',sep=''),
                                         paste('kappa_beta[1:',(SARI_constants$K_t3),', ',s,']', sep=''),paste('sigma_beta[',s,']',sep=''),paste('beta0[',s,']',sep=''),
                                         paste('delta[',s,']',sep=''),paste('omega[',s,']',sep='')),type='AF_slice') 
    SARI_mcmc_config$addSampler(target=c(paste('kappa_gamma[1:',(SARI_constants$K_td),', ',s,']', sep=''),paste('sigma_gamma[1:3, ',s,']',sep=''),paste('gamma0[',s,']',sep='')),type='AF_slice')
    
  }
  SARI_mcmc_config$removeSamplers(target=c(paste('beta_cov_0to18[1:3]',sep=''),paste('beta_cov_19to60[1:3]',sep='')))
  SARI_mcmc_config$addSampler(target=c(paste('beta_cov_0to18[1:3]',sep=''),paste('beta_cov_19to60[1:3]',sep='')),type='AF_slice')
  # Build MCMC.
  SARI_mcmc<- buildMCMC(SARI_mcmc_config)
  # Compile the MCMC.
  SARI_compiled_mcmc <- compileNimble(SARI_mcmc,project=SARI_model)
  # Run MCMC.
  SARI_intermediate_cluster<- runMCMC(SARI_compiled_mcmc,niter=n_iter,nburnin=n_burn,thin=n_thin,inits=SARI_inits,nchains=n_chains,samplesAsCodaMCMC = TRUE)
  
  return(SARI_intermediate_cluster)
}



# Number of available computer cores 
n_cores<-min(detectCores(),length(Nowcast_list))

# CLR model with age covariates 
# Make Cluster for MCMC. 
this_cluster_reg_age<- makeCluster(n_cores)
# Run the cluster. 
time_age <- system.time({
  SARI_age_output <-  parLapply(cl = this_cluster_reg_age, X = 1:length(Nowcast_list), 
                                fun = SARI_age_reg,
                                SARI_x= as.array(SARI_age_x),
                                SARI_z=as.array(SARI_age_z),
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
# Stop the cluster. 
stopCluster(this_cluster_reg_age)
# Save CLR model output. 
save(SARI_age_output,time_age, file='cluster_clr_age_x.RData')

# Check PSRF of model parameters.
sari_parameter_group=c("mu","lambda","zeta0","beta0",
                       "gamma0","chi","theta","phi","delta",
                       "beta_cov_0to18","beta_cov_19to60",
                       "sigma_gamma","sigma_zeta","pi","beta1",
                       "sigma_beta","kappa_zeta","kappa_beta","kappa_gamma","omega")
sari_psrf<-list()
sari_mean_psrf<-sari_max_psrf<-sari_leq_105 <-sari_leq_110<-array(NA,dim=c(length(sari_parameter_group),length(Nowcast_list)))
for(k in 1:length(sari_parameter_group)){
  sari_psrf[[k]]<-list()
  for(j in 1:length(Nowcast_list)){
    sari_parameter_names <- colnames((as.matrix(SARI_age_output[[j]])))
    sari_index <- which(startsWith(sari_parameter_names,sari_parameter_group[k]))
    sari_psrf[[k]][[j]]<- sapply(sari_index,function(v)gelman.diag(SARI_age_output[[j]][,v],autoburnin = FALSE,transform = TRUE,multivariate = FALSE)$psrf[,1])
    sari_mean_psrf[k,j] <- mean(sari_psrf[[k]][[j]])
    sari_max_psrf[k,j] <- max(sari_psrf[[k]][[j]])
    sari_leq_105[k,j] <- mean(sari_psrf[[k]][[j]]<=1.05)
    sari_leq_110[k,j] <- mean(sari_psrf[[k]][[j]][j]<=1.1)
  }
}
sari_MAX_psrf<-cbind(sari_parameter_group,apply(sari_max_psrf,1,max))
sari_MEAN_psrf<-cbind(sari_parameter_group,apply(sari_mean_psrf,1,mean))


plot(SARI_age_output[[2]][,"mu[10, 26]"])
plot(SARI_age_output[[2]][,"x[81, 26]"])
plot(SARI_age_output[[2]][,"x[55, 2]"])
plot(SARI_age_output[[2]][,"y[3, 55]"])
