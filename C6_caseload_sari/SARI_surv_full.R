## MODEL FOR TOTAL HOSPITALISATIONS ##
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
library(rvest)

# Set working directory.
#setwd("~/OneDrive/PhD/comp 2023 paper code")

# Set seed.
set.seed(60209)

#Read in data..
INFLUD21 <- read_delim("INFLUD21-30-01-2023.csv", delim = ";", escape_double = FALSE, trim_ws = TRUE)

INFLUD22 <- read_delim("INFLUD22-30-01-2023.csv", 
                       delim = ";", escape_double = FALSE, 
                       col_types = cols(DT_NOTIFIC = col_date(format = "%d/%m/%Y"), 
                                        DT_SIN_PRI = col_date(format = "%d/%m/%Y")), 
                       trim_ws = TRUE)

# Identify columns of interest:
col_ant<-which(colnames(INFLUD21)=="AN_SARS2")
col_pcr<-which(colnames(INFLUD21)=="PCR_SARS2")
col_results<-which(colnames(INFLUD21)=="PCR_RESUL")
col_hospital<-which(colnames(INFLUD21)=="DT_INTERNA")
col_collection<-which(colnames(INFLUD21)=="DT_COLETA")
col_an_result<-which(colnames(INFLUD21)=="DT_RES_AN")
col_pcr_result<-which(colnames(INFLUD21)=="DT_PCR")
col_notif_date<-which(colnames(INFLUD21)=="DT_DIGITA")


SARI_2021<- INFLUD21[,c(1:16, col_hospital, col_collection, col_results, col_pcr_result, col_notif_date, col_pcr, col_an_result, col_ant)]%>%
  mutate(DT_NOTIFIC=as.Date((DT_NOTIFIC),format="%d/%m/%Y"),DT_DIGITA=as.Date((DT_DIGITA),format="%d/%m/%Y"),DT_SIN_PRI=as.Date((DT_SIN_PRI),format="%d/%m/%Y"))
dates_2021<-sort(as.Date(unique(SARI_2021$DT_SIN_PRI),format="%d/%m/%Y"))

SARI_2022<- INFLUD22[,c(1:16, col_hospital, col_collection, col_results, col_pcr_result, col_notif_date, col_pcr, col_an_result, col_ant)]%>%
  mutate(DT_NOTIFIC=as.Date((DT_NOTIFIC),format="%d/%m/%Y"),DT_DIGITA=as.Date((DT_DIGITA),format="%d/%m/%Y"),DT_SIN_PRI=as.Date((DT_SIN_PRI),format="%d/%m/%Y"))%>%
  filter(DT_SIN_PRI<as.Date("2023-01-01")) # cut off at end of 2022 to ensure fully reported data
dates_2022<-sort(as.Date(unique(SARI_2022$DT_SIN_PRI),format="%d/%m/%Y"))

SARI_2122<-rbind(SARI_2021,mutate(SARI_2022, SEM_PRI=SEM_PRI+52, SEM_NOT=SEM_NOT+52))%>%
  mutate(delay=difftime(DT_DIGITA,DT_SIN_PRI, unit=c("days")),delay_noti=difftime(DT_NOTIFIC,DT_SIN_PRI, unit=c("days")))%>%
  mutate(DT_COLETA=as.Date((DT_COLETA),format="%d/%m/%Y"),DT_PCR=as.Date((DT_PCR),format="%d/%m/%Y"),DT_NASC=as.Date((DT_NASC),format="%d/%m/%Y"))%>%
  mutate(pcr_delay=difftime(DT_PCR,DT_COLETA, unit=c("days")))%>%
  mutate(delay_diffweeks=as.numeric(difftime(DT_DIGITA,DT_SIN_PRI, unit=c("weeks")))) #


# Earliest date we have data for:
first_date<-min(SARI_2122$DT_SIN_PRI)

# We shall use the difference between digitalization date (DT_DIGITA) and date of onset of symptoms (DT_SIN_PRI).
# Determine delay in terms of weeks not days:
SARI_2122$delay_week<-ceiling(SARI_2122$delay_diffweeks+1/7)

# Create age groups
SARI_2122<-SARI_2122%>%mutate(age=round(as.numeric(time_length(difftime(DT_SIN_PRI,DT_NASC), "years"))))%>%
  mutate(age_group=cut(age, c(0, 18, 60, Inf), c("0-18", "19-60", ">60"), include.lowest=TRUE))

SARI_2122<-SARI_2122%>%mutate(SARS2=as.numeric((PCR_SARS2==1)|(AN_SARS2==1)))


# Make matrix for each age group with onset week as rows and delay weeks as columns:
# Filter for first age group:
SARI_0_18<-filter(SARI_2122,age_group=="0-18")
# Make matrix with onset week as rows and delay weeks as columns.
SARI_0_18_long<-filter(SARI_0_18, delay_week>0)%>%group_by(SEM_PRI,delay_week, SG_UF_NOT)%>%summarise(cases=n())
# Remove latest week (N_raw) as only 3 regions have observations.
N_raw<-max(SARI_0_18_long$SEM_PRI)
SARI_0_18_long<-filter(SARI_0_18_long, SEM_PRI<N_raw)
SARI_0_18_delay<-SARI_0_18_long%>%ungroup()%>%
  pivot_wider(id_cols=c(SEM_PRI, SG_UF_NOT),names_from = delay_week,values_from = cases)
# If entry is missing then zero cases were reported. 
SARI_0_18_delay[is.na(SARI_0_18_delay)]<-0

# Filter for second age group:
SARI_19_60<-filter(SARI_2122,age_group=="19-60")
# Make matrix with onset week as rows and delay weeks as columns
SARI_19_60_long<-filter(SARI_19_60, delay_week>0)%>%group_by(SEM_PRI,delay_week, SG_UF_NOT)%>%summarise(cases=n())
# Remove latest week (N_raw) as only 3 regions have observations.
N_raw<-max(SARI_19_60_long$SEM_PRI)
SARI_19_60_long<-filter(SARI_19_60_long, SEM_PRI<N_raw)
SARI_19_60_delay<-SARI_19_60_long%>%ungroup()%>%
  pivot_wider(id_cols=c(SEM_PRI, SG_UF_NOT),names_from = delay_week,values_from = cases)
# If entry is missing then zero cases were reported. 
SARI_19_60_delay[is.na(SARI_19_60_delay)]<-0

# Filter for third age group:
SARI_61<-filter(SARI_2122,age_group==">60")
# Make matrix with onset week as rows and delay weeks as columns
SARI_61_long<-filter(SARI_61, delay_week>0)%>%group_by(SEM_PRI,delay_week, SG_UF_NOT)%>%summarise(cases=n())
# Remove latest week (N_raw) as only 3 regions have observations.
N_raw<-max(SARI_61_long$SEM_PRI)
SARI_61_long<-filter(SARI_61_long, SEM_PRI<N_raw)
SARI_61_delay<-SARI_61_long%>%ungroup()%>%
  pivot_wider(id_cols=c(SEM_PRI, SG_UF_NOT),names_from = delay_week,values_from = cases)
# If entry is missing then zero cases were reported. 
SARI_61_delay[is.na(SARI_61_delay)]<-0

# Spatial names (federal_units) and number of regions (S).
federal_units<-sort(unique(SARI_19_60_delay$SG_UF_NOT))
S<-length(federal_units)
# Number of weeks in final data frame.
N_max<-N_raw-1
# Number of age groups.
A<-3
# Maximum Delay to consider.
D_max<-20
# Number of delays explicitly modelled in GDM.
D<-8


# Order attributes in arrays:
SARI_0_18_delay_ordered<-SARI_0_18_delay[,1:(D_max+2)]%>%mutate(id=paste(SEM_PRI, SG_UF_NOT, sep=''))
SARI_0_18_delay_ordered<-SARI_0_18_delay_ordered%>%arrange(SEM_PRI, SG_UF_NOT)
SARI_19_60_delay_ordered<-SARI_19_60_delay[,1:(D_max+2)]%>%mutate(id=paste(SEM_PRI, SG_UF_NOT, sep=''))
SARI_19_60_delay_ordered<-SARI_19_60_delay_ordered%>%arrange(SEM_PRI, SG_UF_NOT)
SARI_61_delay_ordered<-SARI_61_delay[,1:(D_max+2)]%>%mutate(id=paste(SEM_PRI, SG_UF_NOT, sep=''))
SARI_61_delay_ordered<-SARI_61_delay_ordered%>%arrange(SEM_PRI, SG_UF_NOT)

# Create blank data frame for all regions and weeks. 
SARI_0_18_full <- cbind(tibble(SEM_PRI=sort(rep(1:N_max,27)), SG_UF_NOT=rep(sort(federal_units),N_max)),
                        matrix(0, ncol=D_max,nrow=S*N_max))%>%mutate(check = NA)%>%mutate(id=paste(SEM_PRI, SG_UF_NOT, sep=''))
SARI_19_60_full <- cbind(tibble(SEM_PRI=sort(rep(1:N_max,27)), SG_UF_NOT=rep(sort(federal_units),N_max)),
                         matrix(0, ncol=D_max,nrow=S*N_max))%>%mutate(check = NA)%>%mutate(id=paste(SEM_PRI, SG_UF_NOT, sep=''))
SARI_61_full <- cbind(tibble(SEM_PRI=sort(rep(1:N_max,27)), SG_UF_NOT=rep(sort(federal_units),N_max)),
                      matrix(0, ncol=D_max,nrow=S*N_max))%>%mutate(check = NA)%>%mutate(id=paste(SEM_PRI, SG_UF_NOT, sep=''))

# Fill in available data:
SARI_0_18_full[SARI_0_18_full$id %in% SARI_0_18_delay_ordered$id ,3:(D_max+3)] <- SARI_0_18_delay_ordered[,3:(D_max+3)]
SARI_19_60_full[SARI_19_60_full$id %in% SARI_19_60_delay_ordered$id ,3:(D_max+3)] <- SARI_19_60_delay_ordered[,3:(D_max+3)]
SARI_61_full[SARI_61_full$id %in% SARI_61_delay_ordered$id ,3:(D_max+3)] <- SARI_61_delay_ordered[,3:(D_max+3)]


# Make an array including delay:
# SARI_..._array[REGION, TIME, DELAY]
SARI_0_18_array<-SARI_0_18_full[,3:(D_max+2)]%>%as.matrix()%>%array(dim=c(S,N_max,D_max))
SARI_19_60_array<-SARI_19_60_full[,3:(D_max+2)]%>%as.matrix()%>%array(dim=c(S,N_max,D_max))
SARI_61_array<-SARI_61_full[,3:(D_max+2)]%>%as.matrix()%>%array(dim=c(S,N_max,D_max))
# SARI_array[REGION, TIME, DELAY, AGE GROUP]
SARI_array<-array(abind(SARI_0_18_array, SARI_19_60_array, SARI_61_array),dim=c(S,N_max,D_max,A))
# Total SARI cases for each region.
total_SARI<-apply(SARI_array,c(1,4), sum)
SARI_delay_totals<-melt(total_SARI)
colnames(SARI_delay_totals)<-c('SG_UF_NOT',"AGE GROUP","totals")


# Matrix for partial SARI cases.
SARI_linear_z<-SARI_array[,,1:D_max,] #[REGION, TIME, DELAY, AGE GROUP]
# Matrix for total SARI cases.
SARI_linear_y<-apply(SARI_linear_z,c(1,2,4), sum) #[REGION, TIME, AGE GROUP]
# Dates of all weeks in data frame. 
SARI_dates<-as.Date(first_date)+(1:ncol(SARI_linear_y)-1)*7

# Make array for COVID-19 cases:
COVID_age<-filter(SARI_2122, delay_week>0, delay_week<21)%>%group_by(SEM_PRI,age_group, SARS2, SG_UF_NOT)%>%summarise(cases=n())%>%drop_na(age_group, SARS2)%>%
  pivot_wider(id_cols=c(SEM_PRI, SG_UF_NOT),names_from = age_group,values_from = cases)
# If entry is missing then zero cases were reported. 
COVID_age[is.na(COVID_age)]<-0

# Remove latest week (N_raw) as only 3 regions have observations.
COVID_long<-filter(COVID_age, SEM_PRI<N_raw)
# Add ID and order matrix.
COVID_long_ordered<-COVID_long%>%mutate(id=paste(SEM_PRI, SG_UF_NOT, sep=''))%>%arrange(SEM_PRI, SG_UF_NOT)

# Create blank data frame for all regions and time. 
COVID_full <- tibble(SEM_PRI=sort(rep(1:N_max,27)), SG_UF_NOT=rep(sort(federal_units),N_max),
                     `0-18` = 0, `19-60` = 0, `>60` = 0, check = NA)
COVID_full<-COVID_full%>%mutate(id=paste(SEM_PRI, SG_UF_NOT, sep=''))

# Fill in available data:
COVID_full[COVID_full$id %in% COVID_long_ordered$id ,3:6] <- COVID_long_ordered[,3:6]
COVID_wide<-COVID_full[,-c(1,2,6,7)]%>%as.matrix()%>%array(dim=c(S,N_max,A))

# Matrix for COVID cases.
SARI_linear_x<-COVID_wide # [REGION, TIME, AGE]

SARI_x<-t(apply(SARI_linear_x,c(1,2), sum))
SARI_z<-apply(SARI_linear_z,c(1,2,3), sum)
SARI_y<-t(apply(SARI_linear_y,c(1,2), sum))

# PERCENT REPORTED AFTER EIGHT DELAYS
median(apply(SARI_z[,,1:8],c(1,2), sum)/apply(SARI_z,c(1,2), sum), na.rm = TRUE)

# Melt arrays/matrix to create data frames:
# Data frame for COVID cases.
x_all<-melt(apply(SARI_linear_x,c(1,2), sum))
colnames(x_all)<-c("s","t","x")
x_all$s<-federal_units[x_all$s]
# Data frame for SARI cases.
y_all<-melt(apply(SARI_linear_y,c(1,2), sum))
colnames(y_all)<-c("s","t","y")
y_all$s<-federal_units[y_all$s]

# Population Data:
# Reading in the table from Wikipedia.
page = read_html("https://en.wikipedia.org/wiki/Federative_units_of_Brazil")
# Obtain the piece of the web page that corresponds to the "wikitable" node.
my.table = html_node(page, ".wikitable")
# Convert the html table element into a data frame.
my.table = html_table(my.table, fill = TRUE)
population<-my.table[,c(1,2,6)]
colnames(population)<-c("s_full","s","population")
population<-population[-1,]
population<-population[order(population$s),]
population[,3]<-as.numeric(str_replace_all(as.character(unlist(population[,3])),",",""))

sari_wide<-melt(SARI_z)
colnames(sari_wide)<-c("s","t","d","z")
Data_wide<-sari_wide%>%ungroup()%>%
  pivot_wider(id_cols=c(t,s),names_from = d,values_from = z)%>%as.matrix()
Data_wide<-Data_wide[order(Data_wide[,2],Data_wide[,1]),]

Data_delay_sum<-as.tibble(Data_wide)
totals<-Data_delay_sum%>%mutate(total=as.numeric(unlist((apply(Data_wide[,3:22],1, sum, na.rm=TRUE)))))%>%mutate(s=paste(federal_units[s]))
Data_delay_sum$prop_0<-as.numeric(unlist((Data_wide[,3])/apply(Data_wide[,3:22],1, sum, na.rm=TRUE)))
Data_delay_sum$prop_1<-as.numeric(unlist((Data_wide[,3]+Data_wide[,4])/apply(Data_wide[,3:22],1, sum, na.rm=TRUE)))
Data_delay_sum$prop_2<-as.numeric(unlist((Data_wide[,3]+Data_wide[,4]+Data_wide[,5])/apply(Data_wide[,3:22],1, sum, na.rm=TRUE)))
Data_delay_sum$prop_3<-as.numeric(unlist((Data_wide[,3]+Data_wide[,4]+Data_wide[,5]+Data_wide[,6])/apply(Data_wide[,3:22],1, sum, na.rm=TRUE)))
Data_delay_sum$prop_4<-as.numeric(unlist((Data_wide[,3]+Data_wide[,4]+Data_wide[,5]+Data_wide[,6]+Data_wide[,7])/apply(Data_wide[,3:22],1, sum, na.rm=TRUE)))
# 
# Data_delay_sum$prop_0<-as.numeric(unlist((Data_wide[,3])/apply(Data_wide[,3:22],1, sum, na.rm=TRUE)))
# Data_delay_sum$prop_1<-as.numeric(unlist((Data_wide[,4])/(apply(Data_wide[,3:22],1, sum, na.rm=TRUE)-Data_wide[,3])))
# Data_delay_sum$prop_2<-as.numeric(unlist((Data_wide[,5])/(apply(Data_wide[,3:22],1, sum, na.rm=TRUE)-apply(Data_wide[,3:4],1, sum, na.rm=TRUE))))
# Data_delay_sum$prop_3<-as.numeric(unlist((Data_wide[,6])/(apply(Data_wide[,3:22],1, sum, na.rm=TRUE)-apply(Data_wide[,3:5],1, sum, na.rm=TRUE))))
# Data_delay_sum$prop_4<-as.numeric(unlist((Data_wide[,7])/(apply(Data_wide[,3:22],1, sum, na.rm=TRUE)-apply(Data_wide[,3:6],1, sum, na.rm=TRUE))))
Data_delay_long<-as.tibble(Data_delay_sum%>%pivot_longer(!c(t,s), names_to = "delay", values_to = "count")%>%mutate(s=paste(federal_units[s])))
Data_delay_pop<-full_join(Data_delay_long,population, by="s")
Data_delay_prop<-Data_delay_pop%>%filter(delay%in%c("prop_0","prop_1","prop_2","prop_3","prop_4"))
Data_delay_prop<-full_join(Data_delay_prop,totals[,c(1,2,dim(totals)[2])],by = join_by(t, s))
Data_delay_prop<-Data_delay_prop%>%mutate(total.per.capita=total/population)

# plot cumulative proportions against totals
ggplot(Data_delay_prop,aes(x=total,y=probit(count), colour=delay))+
  geom_point(alpha=0.1)+facet_wrap(~s, scales="free")+geom_smooth(alpha=0.1)+scale_x_log10()+
  labs(x="Total counts",y="Cumulative proportions")
ggplot(Data_delay_prop,aes(x=total,y=probit(count), colour=s))+
  geom_point(alpha=0.1)+facet_wrap(~delay, scales="free")+geom_smooth(alpha=0.1)+scale_x_log10()+
  labs(x="Total counts",y="Cumulative proportions")

y_long<-melt(SARI_y, value.name = 'y', varnames = c("t","s"))

# Create function for cluster to carry out MCMC model:
SARI_linear_code <- function(seed, SARI_z, SARI_y, population, N, S, D, D_max, n_knots, n_iter, n_burn, n_thin, n_chains){
  # Load libraries.
  library(nimble)
  library(tidyverse)
  library(mgcv)
  
  # Define the Beta-Binomial as a distribution for NIMBLE.
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
  SARI_code <- nimbleCode({
   
      for(t in 1:N){
        # Negative Binomial Model for total SARI cases.
        y[t] ~ dnegbin(prob=theta/(theta+lambda[t]), size=theta)  
        
        # Model for partial delayed SARI counts (first delay).
        z[t,1] ~ dbetabin(nu[t,1],phi[1],y[t])
      }
      for(d in 2:D){
        for(t in 1:N){
          # Model for partial delayed SARI counts.
          z[t,d] ~ dbetabin(nu[t,d],phi[d],y[t]-sum(z[t,1:(d-1)]))
        }
      }
      for(t in 1:N){
        # Linear predictor (Beta-Binomial means).
        for(d in 1:D){
          # Expected cumulative proportions.
          probit(p[t,d]) <- psi[d] + eta[t] #+ delta*log((y[t]+1)/pop)
        }
        # Relative proportions (Beta-Binomial means).
        nu[t,1] <- p[t,1]
        for(d in 2:D){
          nu[t,d] <- (p[t,d]-p[t,d-1])/(1-p[t,d-1])
        }
        # Mean for total SARI cases .
        log(lambda[t]) <- log(pop) + zeta0 + zeta1[t] + xi[weeks[t]] 
        # Linear predictor for COVID Beta-Binomial probability.
        #    logit(mu[t]) <- beta0[s] + beta1[t] + zeta1[t]*omicron[s] + beta_cov_0to18[1]*age_0to18[t,1] + beta_cov_0to18[2]*age_0to18[t,2] + beta_cov_0to18[3]*age_0to18[t,3] +  
        #     beta_cov_19to60[1]*age_19to60[t,1] + beta_cov_19to60[2]*age_19to60[t,2] + beta_cov_19to60[3]*age_19to60[t,3]
        
      }
      
      
      # Splines of time in cumulative proportion reported.
      Omega_eta[1:K_t2,1:K_t2] <- S_t2[1:K_t2,1:K_t2]/sigma_eta^2
      kappa_eta[1:K_t2] ~ dmnorm(zeros[1:K_t2],Omega_eta[1:K_t2,1:K_t2])
      eta[1:N] <- X_t2[1:N,1:K_t2]%*%kappa_eta[1:K_t2]
      
      # Seasonal spline (for influenza) 
      Omega_xi[1:K_w,1:K_w]<-S_w[1:K_w,1:K_w]/sigma_xi^2
      kappa_xi[1:K_w]~dmnorm(zeros[1:K_w],Omega_xi[1:K_w,1:K_w])
      xi[1:52]<-X_w[1:52,1:K_w]%*%kappa_xi[1:K_w]
      
      
      # Splines of time for proportion of COVID-19/SARI cases.
      #  Omega_beta[1:K_t3,1:K_t3] <- S_t3[1:K_t3,1:K_t3]/sigma_beta[s]^2
      #  kappa_beta[1:K_t3] ~ dmnorm(zeros[1:K_t3],Omega_beta[1:K_t3,1:K_t3])
      #  beta1[1:N] <- X_t3[1:N,1:K_t3]%*%kappa_beta[1:K_t3]
      
      # Splines of time for SARI cases.
      Omega_zeta[1:K_t2,1:K_t2] <- S_t2[1:K_t2,1:K_t2]/sigma_zeta^2
      kappa_zeta[1:K_t2] ~ dmnorm(zeros[1:K_t2],Omega_zeta[1:K_t2,1:K_t2])
      zeta1[1:N] <- X_t2[1:N,1:K_t2]%*%kappa_zeta[1:K_t2]
      
      # for(d in 1:D){
     # delta ~ dnorm(0,1)
      # }
      # Smoothing parameter priors.
      # sigma_beta[s] ~ T(dnorm(0,1),0,)
      sigma_zeta ~ T(dnorm(0,1),0,)
      sigma_eta ~ T(dnorm(0,1),0,) 
      sigma_xi ~ T(dnorm(0,1),0,) 
      
      # Cumulative delay curves.
      psi[1] ~ dnorm(0,sd=10)
      for(d in 2:D){
        psi[d] ~ T(dnorm(psi[d-1], sd=10),psi[d-1],)
      }
      
      # Beta-Binomial dispersion parameters.
      for(d in 1:D){
        phi[d] ~ dgamma(2,0.02) 
      }
      
      # Intercepts
      # beta0[s] ~ dnorm(0, sd=5)
      zeta0  ~ dnorm(-mean.sari, sd=10)
      
      
    
    
  })
  
  # Cubic spline with shrinkage for time (mean SARI cases, cumulative SARI cases reported) %
  # yearly seasonal spline (for 52 weeks):
  week_index<-rep(1:52,9) #first date is 01-01-2013 - NEEDS UPDATING 
  blank_data2=tibble(y=rnorm(N,0,1),t=1:N,w=week_index[1:N]) 
  blank_jagam2=jagam(y~s(t,k=floor(N/10), bs='cs')+
                       s(w,bs='cc',k=n_knots[2]),data=blank_data2,file='blank.jags',
                     knots=list(t=seq(1,N,length=floor(N/10)),w=seq(0,52,length=n_knots[2]))) # #
  
  
  # Constants
  SARI_constants <- list(N=N,
                         D=D,
                         weeks=week_index[1:N],
                         pop=as.numeric(population$population[seed]),
                         #  age_0to18=array(age_group_cubic[,,,1],dim=c(N,3)),
                         #  age_19to60=array(age_group_cubic[,,,2],dim=c(N,3)),
                         mean.sari=floor(mean(log(SARI_y[1:(N),seed]))),
                         K_w=dim(blank_jagam2$jags.data$S2)[1],
                         S_w=blank_jagam2$jags.data$S2,
                         K_t2=dim(blank_jagam2$jags.data$S1)[1],
                         S_t2=blank_jagam2$jags.data$S1)
  
  
  SARI_constants$X_t2 <- blank_jagam2$jags.data$X[,2:(SARI_constants$K_t2+1)]
  SARI_constants$X_w <- blank_jagam2$jags.data$X[,(SARI_constants$K_t2+2):(SARI_constants$K_t2+SARI_constants$K_w+1)]
  
  SARI_constants$zeros <- rep(0,max(SARI_constants$K_t3, SARI_constants$K_t2, SARI_constants$K_w))
  
  # Data:
  SARI_data <- list(y=SARI_y[1:N,seed],
                    #  x=as.matrix(censored_x[1:N,],ncol=S),
                    #  c=as.matrix(censored_x[1:N,],ncol=S),
                    z=matrix(SARI_z[seed,1:N,1:D],ncol=D)) 
  #SARI_data$x[(M+1):N,] <- NA
  
  SARI_inits<-SARI_model<-SARI_compiled_model<- SARI_mcmc_config<-SARI_mcmc <- SARI_compiled_mcmc<-list()
  # Generate random initial values.
  SARI_inits <- list(
  #  delta=rnorm(1,0,0.1),#matrix(rnorm(D*S,0,0.1),ncol=D),
    kappa_xi=rnorm((SARI_constants$K_w),0,0.1),
    kappa_zeta=rnorm((SARI_constants$K_t2),0,0.01),
    kappa_eta=rnorm((SARI_constants$K_t2),0,0.1),
    sigma_xi=runif(1,0,1),
    psi=sort(rnorm(D,0,1)),
    sigma_eta=runif(1,0,1),
    sigma_zeta=runif(1,0,1),
    theta=abs(rnorm(1,50,10)),
    zeta0=rnorm(1,-10,1) )
  
  # Build the model.
  SARI_model <- nimbleModel(SARI_code, SARI_constants, SARI_data, SARI_inits)
  # Compile the model.
  SARI_compiled_model <- compileNimble(SARI_model)
  
  # Set up the MCMC.
  SARI_mcmc_config <- configureMCMC(SARI_model,monitors=c("lambda","zeta0","y",#"phi","chi",
                                                          "psi","theta",#"p",
                                                          #"omicron","mu",
                                                          "sigma_eta","sigma_zeta", "xi",
                                                          "zeta1",#"delta",
                                                          #"kappa_beta", "beta_cov_0to18","beta_cov_19to60"
                                                          #'pi','beta1',"beta0","omega","x","sigma_beta",
                                                          "eta",
                                                          "kappa_zeta","kappa_eta"
  ),
  useConjugacy = FALSE)
  
  SARI_mcmc_config$removeSamplers('kappa_zeta','sigma_zeta',#'kappa_beta','sigma_beta',
                                  'kappa_xi','sigma_xi', 'psi',
                                  #'omega','omicron','beta0','delta',
                                  'kappa_eta','sigma_eta','zeta0')
    SARI_mcmc_config$addSampler(target=c(paste('kappa_zeta[1:',(SARI_constants$K_t2),']', sep=''),paste('sigma_zeta',sep=''),paste('zeta0',sep=''),
                                         paste('kappa_xi[1:',(SARI_constants$K_w),']', sep=''),paste('sigma_xi', sep='')
                                         #   paste('kappa_beta[1:',(SARI_constants$K_t3),', ',']', sep=''),paste('sigma_beta[',']'ep=''),paste('beta0[',']'ep=''),
                                         #   paste('omicron[',']'ep=''),paste('omega[',']'ep='')
    ),type='AF_slice') 
    SARI_mcmc_config$addSampler(target=c(# paste('delta', sep=''),
                                          paste('kappa_eta[1:',(SARI_constants$K_t2),']', sep=''),paste('sigma_eta',sep=''),
                                          paste('psi[1:',D,']', sep='')),type='AF_slice')
  
  
  #SARI_mcmc_config$removeSamplers(target=c(paste('beta_cov_0to18[1:3]'ep=''),paste('beta_cov_19to60[1:3]'ep='')))
  # SARI_mcmc_config$addSampler(target=c(paste('beta_cov_0to18[1:3]'ep=''),paste('beta_cov_19to60[1:3]'ep='')),type='AF_slice')
  # Build MCMC. 
  SARI_mcmc<- buildMCMC(SARI_mcmc_config)
  # Compile the MCMC.
  SARI_compiled_mcmc <- compileNimble(SARI_mcmc,project=SARI_model)
  # Run MCMC. 
  SARI_cluster<- runMCMC(SARI_compiled_mcmc,niter=n_iter,nburnin=n_burn,thin=n_thin,inits=SARI_inits,nchains=n_chains, samplesAsCodaMCMC = TRUE)
  
  return(SARI_cluster)
}


# MCMC parameters.
# Set chains.
n_chains=2
# Set iterations.
n_iter=10000
# Set burn-in.
n_burn=5000
# Set thinning.
n_thin=5
# Set number of knots for model splines
n_knots<-c(15,7,8,10) # knots for time, time-delay interaction and delay respectively 
# Set maximum COVID-19 lab report delay.
R<-30  
N<-(N_max-D_max)
y_long<-filter(y_long,t<=N)

# List of weeks to perform the nowcasts for rolling predictions experiment.
Nowcast_list<-floor(seq(from=100, to=N-D_max, length=4))
# Dates for the nowcasts.
Nowcast_dates<- as.Date(first_date)+(Nowcast_list-1)*7

# Number of available computer cores. 
n_cores<-min(detectCores(),S)

# Survivor model with age covariates. 
# Make Cluster for MCMC. 
this_cluster_full<- makeCluster(n_cores)

# Run Cluster.
time_surv_full<- system.time({
  SARI_linear_output_delta <-  parLapply(cl = this_cluster_full, X = 1:S, 
                                       fun = SARI_linear_code,
                                       SARI_z=SARI_z,
                                       SARI_y=SARI_y,
                                       population=population,
                                       N=N,
                                       S=S,
                                       D=D,
                                       D_max=D_max,
                                       n_iter=n_iter,
                                       n_burn=n_burn,
                                       n_thin=n_thin,
                                       n_chains=n_chains,
                                       n_knots=n_knots)
})

# Stop Cluster.
stopCluster(this_cluster_full)

# Save survivor model output. 
save(SARI_linear_output_delta,time_surv_full, file='plots/survivor_FULL_nolink.RData')
# Check PSRF of model parameters.
sari_parameter_group=c("zeta0","xi","lambda",
                       "psi",#"delta",
                       "sigma_eta","sigma_zeta", "zeta1",
                       "kappa_zeta","kappa_eta")
sari_psrf<-list()
sari_mean_psrf<-sari_max_psrf<-sari_leq_105 <-sari_leq_110<-array(NA,dim=c(length(sari_parameter_group),length(Nowcast_list)))
for(k in 1:length(sari_parameter_group)){
  sari_psrf[[k]]<-list()  
  for(j in 1:length(Nowcast_list)){
    sari_parameter_names <- colnames((as.matrix(SARI_linear_output_delta[[j]])))
    sari_index <- which(startsWith(sari_parameter_names,sari_parameter_group[k]))
    sari_psrf[[k]][[j]]<- sapply(sari_index,function(v)gelman.diag(SARI_linear_output_delta[[j]][,v],autoburnin = FALSE,transform = TRUE,multivariate = FALSE)$psrf[,1])
    sari_mean_psrf[k,j] <- mean(sari_psrf[[k]][[j]])
    sari_max_psrf[k,j] <- max(sari_psrf[[k]][[j]])
    sari_leq_105[k,j] <- mean(sari_psrf[[k]][[j]]<=1.05)
    sari_leq_110[k,j] <- mean(sari_psrf[[k]][[j]][j]<=1.1)
  }
}
sari_MAX_psrf<-cbind(sari_parameter_group,apply(sari_max_psrf,1,max))
sari_MEAN_psrf<-cbind(sari_parameter_group,apply(sari_mean_psrf,1,mean))



plot(SARI_linear_output_delta[[1]][,"delta"],main="delta")
plot(SARI_linear_output_delta[[2]][,"delta"],main="delta")

#### Plots #####
library(geobr)
states<-read_state(year=2013)
state_guide<-tibble(s=states$code_state, s_abb=states$abbrev_state,s_full=states$name_state)
geom_guide<-tibble(s=states$abbrev_state,s_full=states$name_state,geom=states$geom)

n_sim<-dim(SARI_linear_output_delta[[1]]$chain1)[1]*2
SARI_Survivor_linear_x<-list()
SARI_combined_Survivor_linear_x<-list()
surv_linear_mu_x<-covid_surv_linear_quant_x<-sari_surv_linear_quant_x<-list()
SARI_pi_samples_surv_age<-SARI_mu_samples_surv_age<-SARI_y_samples_surv_age<-SARI_chi_surv_age<-SARI_x_samples_surv_age<-list()
for(j in 1:S){
  # Survivor age model
  SARI_Survivor_linear_x[[j]]<-as.mcmc.list(SARI_linear_output_delta[[j]])
  SARI_combined_Survivor_linear_x[[j]] <- as_tibble(do.call('rbind',SARI_Survivor_linear_x[[j]]))
}

# Model coefficient output
delta_slope<-zeta_spline<-eta_spline<-list()
for(j in 1:S){
  delta_slope[[j]]<-select(SARI_combined_Survivor_linear_x[[j]],starts_with('delta'))%>%as.matrix()%>%array(dim=c(n_sim))%>%
    quantile(c(0.025,0.5,0.975))%>%melt(varnames=c('quantile'),value.name='y')%>%mutate(s=j,quantile=c(0.025,0.5,0.975))%>%
    spread(quantile,y)

  zeta_spline[[j]]<-select(SARI_combined_Survivor_linear_x[[j]],starts_with('zeta1'))%>%as.matrix()%>%array(dim=c(n_sim,N))%>%
    apply(c(2),quantile,c(0.025,0.5,0.975))%>%melt(varnames=c('quantile','t'),value.name='y')%>%mutate(s=j)%>%
    spread(quantile,y)
  
  eta_spline[[j]]<-select(SARI_combined_Survivor_linear_x[[j]],starts_with('eta'))%>%as.matrix()%>%array(dim=c(n_sim,N))%>%
    apply(c(2),quantile,c(0.025,0.5,0.975))%>%melt(varnames=c('quantile','t'),value.name='y')%>%mutate(s=j)%>%
    spread(quantile,y)
  
}
#delta matrix
delta_matrix<-matrix(unlist(delta_slope), byrow = TRUE, ncol = 4)%>%as_tibble()
colnames(delta_matrix)<-c('s','2.5%','50%','97.5%')
delta_matrix$s_abb<-factor(federal_units[delta_matrix$s],levels=federal_units[delta_matrix$s])
delta_matrix<-full_join(delta_matrix,state_guide[,2:3], by='s_abb')
delta_matrix$s_full<-factor(delta_matrix$s_full, levels=delta_matrix$s_full)
#eta matrix
eta_matrix<-as_tibble(eta_spline[[1]],ncol = 5)
for(s in 2:S){
  eta_matrix<-rbind(eta_matrix, eta_spline[[s]])
  
}
colnames(eta_matrix)<-c('t','s','2.5%','50%','97.5%')
#zeta matrix
zeta_matrix<-as_tibble(zeta_spline[[1]],ncol = 5)
for(s in 2:S){
  zeta_matrix<-rbind(zeta_matrix, zeta_spline[[s]])
  
}
colnames(zeta_matrix)<-c('t','s','2.5%','50%','97.5%')

delta_comp_plot<-ggplot()+geom_hline(yintercept = 0,colour='darkgrey',linetype='dashed')+
  geom_errorbar(data=delta_matrix,
                aes(x=s_abb,  ymin=`2.5%`,ymax=`97.5%`,colour=s_full))+
  geom_point(data=delta_matrix,aes(x=s_abb, y=`50%`,colour=s_full),shape=3,size=5)+
  scale_colour_viridis_d(name='Region',begin=0.1,end=0.9)+
  scale_fill_viridis_d(name='Region',begin=0.1,end=0.9)+
  labs(x="Brazilian federative unit (s)",y=expression(delta['s']),
       title='SARI hospitalisations',
       subtitle='Linear effect of case load on reporting delay', legend=NULL)+
  theme_minimal()+
  theme(text = element_text(size=12), plot.title = element_text(size=16),plot.subtitle=element_text(size=14),
        legend.title = element_text(size=14),strip.text = element_text(size=14), 
        legend.position = "bottom",legend.box="vertical")+guides(colour=guide_legend(nrow=6,byrow=TRUE))
delta_comp_plot
ggsave(delta_comp_plot,file='plots/full_delta_sari.pdf',width=9,height=6)

delta_matrix$s<-as.character(delta_matrix$s)
delta_matrix_geom<-full_join(delta_matrix,geom_guide, by='s_full')

# customize colour palatte
lightGreen = "#9BE2B7"
lightRed = "#F19999"
lightBlue="#9CD4EF"
darkGreen = "#32B767"
darkRed = "#DD5A5A"
darkBlue="#1491CD"
red <- "#E6194B"
green <- "#2AA198"
blue <- "#4363D8"
yellow <- "#FFE119"
gray <- "#A9A9A9"

# plot delta on map of Brazil
max(max(delta_matrix_geom$`50%`),abs(min(delta_matrix_geom$`50%`)))
limit=c(-max(max(delta_matrix_geom$`50%`),abs(min(delta_matrix_geom$`50%`))),max(max(delta_matrix_geom$`50%`),abs(min(delta_matrix_geom$`50%`))))
delta_link_map<-ggplot() +
  geom_sf(data=states,colour="grey20", size=.05) +
  geom_sf(data=delta_matrix_geom,aes(geometry=geom, fill=`50%`),colour=NA, size=.15) +
#  scale_fill_fermenter(palette = "spectral",limit=limit,name=expression(delta['s']),direction = -1)+ 
  scale_fill_gradient2(low = red, 
                       mid = yellow,
                       high = green ,name=expression(delta['s']))+
  labs(title="SARI hospitalisations" ,subtitle="linear case load effect of reporting delay by federative unit", size=8) +
  theme_minimal()
delta_link_map
ggsave(delta_link_map,file='plots/delta_link_map_sari_cb.pdf',width=9,height=7,dpi = 600)


my.table2<-my.table[-1,c(1:5,7,8,10,11)]
my.table2[,5]<-as.numeric(str_replace_all(as.character(unlist(my.table2[,5])),",",""))
my.table2[,6]<-as.numeric(str_replace_all(as.character(unlist(my.table2[,6])),",",""))
my.table2[,7]<-as.numeric(unlist(my.table2[,7]))
my.table2[,8]<-as.numeric(str_replace_all(as.character(unlist(my.table2[,8])),",",""))
my.table2[,9]<-as.numeric(unlist(my.table2[,9]))
colnames(my.table2)<-c('name.full','s_abb','capital','largest.city','area.km2','population','density.perkm2','GDP','HDI')
#load("~/media/alba/Disk 1/OneDrive/PhD/comp 2023 paper code/FInal_code/OMEGA CODE/SARI totals link/health_expenditure.RData")
load("~/media/alba/Disk 1/OneDrive/PhD/Case load case studies/SARI totals link/health_expenditure.RData")

delta_explore<-full_join(delta_matrix,my.table2,by='s_abb')
delta_explore<-full_join(delta_explore,health_expenditure, by='name.full')

delta_explore<-delta_explore%>%mutate(mean.cases=apply(SARI_y,2,mean))

delta_explore_all<-delta_explore[,c(3,5,7,12,13,15,16)]%>%pivot_longer(cols = 4:7)

ALL_sari<-ggplot()+
  geom_smooth(data=delta_explore_all,aes(x=value, y=`50%`), method='lm', alpha=0.2, colour="grey")+
  geom_hline(yintercept = 0,colour='darkgrey',linetype='dashed')+
  facet_wrap(~name, ncol=1, scales = "free_x",  strip.position="top",
             labeller = as_labeller(c(density.perkm2='Population density (per Km squared on log scale):',
                                      expenditure='Goverment health expenditure (billion Brazillian reals on log scale):',
                                      GDP='Gross domestic product (on log scale):',
                                      mean.cases='Mean SARI cases over all years (on log scale):')))+
  geom_label(data=delta_explore_all,aes(x=value, y=`50%`,label=s_abb,colour=name.full),alpha=0.5)+
  scale_colour_viridis_d(name='Region',begin=0.1,end=0.9)+
  scale_fill_viridis_d(name='Region',begin=0.1,end=0.9)+
  labs(x=NULL,y=expression(delta['s']),title='SARI hospitalisations',
       subtitle='Linear effect of case load on cumulative proportions reported', legend=NULL)+
  scale_x_log10()+theme_minimal()+ 
  theme(text = element_text(size=12), plot.title = element_text(size=16),plot.subtitle=element_text(size=14),
        legend.title = element_text(size=12),strip.text = element_text(size=12), legend.position = "right")+guides(colour=guide_legend(ncol=1,byrow=TRUE))
ggsave(ALL_sari,file='plots/delta_sari_all.pdf',width=9,height=11)


#
delta_cases<-ggplot()+geom_hline(yintercept = 0,colour='darkgrey',linetype='dashed')+
  geom_errorbar(data=delta_explore,
                aes(x=mean.cases,  ymin=`2.5%`,ymax=`97.5%`,colour=s_abb))+
  geom_point(data=delta_explore,aes(x=mean.cases, y=`50%`,colour=s_abb),shape=3,size=4)+
  scale_colour_viridis_d(name='Region',begin=0.1,end=0.9)+
  scale_fill_viridis_d(name='Region',begin=0.1,end=0.9)+
  labs(x="Mean total cases",y=expression(delta['s']),title='Linear effect of case load on delay', legend=NULL)+
  theme(text = element_text(size=12), plot.title = element_text(size=16),plot.subtitle=element_text(size=14),
        legend.title = element_text(size=14),strip.text = element_text(size=14))+scale_x_log10()
delta_cases
ggsave(delta_cases,file='plots/delta_cases.pdf',width=9,height=3)

delta_pop<-ggplot()+geom_hline(yintercept = 0,colour='darkgrey',linetype='dashed')+
  geom_errorbar(data=delta_explore,
                aes(x=population,  ymin=`2.5%`,ymax=`97.5%`,colour=s))+
  geom_point(data=delta_explore,aes(x=population, y=`50%`,colour=s),shape=3,size=4)+
  scale_colour_viridis_d(name='Region',begin=0.1,end=0.9)+
  scale_fill_viridis_d(name='Region',begin=0.1,end=0.9)+
  labs(x="Population (log)",y=expression(delta['s']),title='Linear effect of case load on delay', legend=NULL)+
  theme(text = element_text(size=12), plot.title = element_text(size=16),plot.subtitle=element_text(size=14),
        legend.title = element_text(size=14),strip.text = element_text(size=14))+scale_x_log10()
ggsave(delta_pop,file='plots/delta_pop.pdf',width=9,height=3)

delta_area<-ggplot()+geom_hline(yintercept = 0,colour='darkgrey',linetype='dashed')+
  geom_errorbar(data=delta_explore,
                aes(x=area.km2,  ymin=`2.5%`,ymax=`97.5%`,colour=s))+
  geom_point(data=delta_explore,aes(x=area.km2, y=`50%`,colour=s),shape=3,size=4)+
  scale_colour_viridis_d(name='Region',begin=0.1,end=0.9)+
  scale_fill_viridis_d(name='Region',begin=0.1,end=0.9)+
  labs(x="Area (kM squared log scale)",y=expression(delta['s']),title='Linear effect of case load on delay', legend=NULL)+
  theme(text = element_text(size=12), plot.title = element_text(size=16),plot.subtitle=element_text(size=14),
        legend.title = element_text(size=14),strip.text = element_text(size=14))+scale_x_log10()
ggsave(delta_area,file='plots/delta_area.pdf',width=9,height=3)

delta_density<-ggplot()+geom_hline(yintercept = 0,colour='darkgrey',linetype='dashed')+
  geom_errorbar(data=delta_explore,
                aes(x=density.perkm2,  ymin=`2.5%`,ymax=`97.5%`,colour=s))+
  geom_point(data=delta_explore,aes(x=density.perkm2, y=`50%`,colour=s),shape=3,size=4)+
  scale_colour_viridis_d(name='Region',begin=0.1,end=0.9)+
  scale_fill_viridis_d(name='Region',begin=0.1,end=0.9)+
  labs(x="Density (per Km squared log scale)",y=expression(delta['s']),title='Linear effect of case load on delay', legend=NULL)+
  theme(text = element_text(size=12), plot.title = element_text(size=16),plot.subtitle=element_text(size=14),
        legend.title = element_text(size=14),strip.text = element_text(size=14))+scale_x_log10()
ggsave(delta_density,file='plots/delta_density.pdf',width=9,height=3)

delta_GDP<-ggplot()+geom_hline(yintercept = 0,colour='darkgrey',linetype='dashed')+
  geom_errorbar(data=delta_explore,
                aes(x=GDP,  ymin=`2.5%`,ymax=`97.5%`,colour=s))+
  geom_point(data=delta_explore,aes(x=GDP, y=`50%`,colour=s),shape=3,size=4)+
  scale_colour_viridis_d(name='Region',begin=0.1,end=0.9)+
  scale_fill_viridis_d(name='Region',begin=0.1,end=0.9)+
  labs(x="GDP (log)",y=expression(delta['s']),title='Linear effect of case load on delay', legend=NULL)+
  theme(text = element_text(size=12), plot.title = element_text(size=16),plot.subtitle=element_text(size=14),
        legend.title = element_text(size=14),strip.text = element_text(size=14))+scale_x_log10()
ggsave(delta_GDP,file='plots/delta_GDP.pdf',width=9,height=3)

delta_HDI<-ggplot()+geom_hline(yintercept = 0,colour='darkgrey',linetype='dashed')+
  geom_errorbar(data=delta_explore,
                aes(x=HDI,  ymin=`2.5%`,ymax=`97.5%`,colour=s))+
  geom_point(data=delta_explore,aes(x=HDI, y=`50%`,colour=s),shape=3,size=4)+
  scale_colour_viridis_d(name='Region',begin=0.1,end=0.9)+
  scale_fill_viridis_d(name='Region',begin=0.1,end=0.9)+
  labs(x="HDI",y=expression(delta['s']),title='Linear effect of case load on delay', legend=NULL)+
  theme(text = element_text(size=12), plot.title = element_text(size=16),plot.subtitle=element_text(size=14),
        legend.title = element_text(size=14),strip.text = element_text(size=14))
ggsave(delta_HDI,file='plots/delta_GDP.pdf',width=9,height=3)


delta_expend<-ggplot()+geom_hline(yintercept = 0,colour='darkgrey',linetype='dashed')+
  geom_errorbar(data=delta_explore,
                aes(x=expenditure,  ymin=`2.5%`,ymax=`97.5%`,colour=s))+
  geom_point(data=delta_explore,aes(x=expenditure, y=`50%`,colour=s),shape=3,size=4)+
  scale_colour_viridis_d(name='Region',begin=0.1,end=0.9)+
  scale_fill_viridis_d(name='Region',begin=0.1,end=0.9)+
  labs(x="Goverment health expenditure (billion Brazillian reals log scale)",y=expression(delta['s']),title='Linear effect of case load on delay', legend=NULL)+
  theme(text = element_text(size=12), plot.title = element_text(size=16),plot.subtitle=element_text(size=14),
        legend.title = element_text(size=14),strip.text = element_text(size=14))+scale_x_log10()
ggsave(delta_expend,file='plots/delta_expend.pdf',width=9,height=3)
#
#

eta_comp_plot<-ggplot()+
  geom_ribbon(data=eta_matrix, 
              aes(x=t,  ymin=`2.5%`,ymax=`97.5%`,fill=as.factor(federal_units[s])),alpha=0.2)+
  geom_line(data=eta_matrix,aes(x=t, y=`50%`,colour=as.factor(federal_units[s])))+
  facet_wrap(~federal_units[s], scales = "free")+
  scale_colour_viridis_d(name='Region',begin=0.1,end=0.9)+
  scale_fill_viridis_d(name='Region',begin=0.1,end=0.9)+
  labs(x="Time",y=expression(eta['s']),title='Temporal trend in delay', legend=NULL)+
  theme(text = element_text(size=12), plot.title = element_text(size=16),plot.subtitle=element_text(size=14),
        legend.title = element_text(size=14),strip.text = element_text(size=14))
 ggsave(eta_comp_plot,file='plots/full_eta.pdf',width=9,height=3)


zeta_comp_plot<-ggplot()+
  geom_ribbon(data=zeta_matrix, 
              aes(x=t,  ymin=`2.5%`,ymax=`97.5%`,fill=federal_units[s]),alpha=0.2)+
  geom_line(data=zeta_matrix,aes(x=t, y=`50%`,colour=federal_units[s]))+
  facet_wrap(~federal_units[s], scales = "free")+
  scale_colour_viridis_d(name='Region',begin=0.1,end=0.9)+
  scale_fill_viridis_d(name='Region',begin=0.1,end=0.9)+
  labs(x="Time",y=expression(zeta['s']),title='Temporal trend in total counts', legend=NULL)+
  theme(text = element_text(size=12), plot.title = element_text(size=16),plot.subtitle=element_text(size=14),
        legend.title = element_text(size=14),strip.text = element_text(size=14))
 ggsave(zeta_comp_plot,file='plots/full_zeta.pdf',width=9,height=3)


stripe_indicator=function(x=double(1)){
  stripe<-x%in%c(1:20,40:60,80:100,120:140,160:180)
  return(stripe)
}


for(i in 1:9){
  plot_data<-full_join(mutate(filter(y_long, s%in%c((3*(i-1)+1):(3*i))),m="data",stripe=stripe_indicator(t)),
                       mutate(filter(eta_matrix,s%in%c((3*(i-1)+1):(3*i))), m="spline",stripe=stripe_indicator(t)),
                       by=c("t","s","stripe",'m'))
  assign(paste('plot',i,sep=''),ggplot(data=plot_data)+ geom_rect(aes(xmax = t + 1, 
                                                                      xmin = t, 
                                                                      ymin = -Inf, 
                                                                      ymax = Inf, 
                                                                      fill = stripe), alpha = 0.4, show.legend = FALSE) + 
           scale_fill_manual(values = c("white", "grey50")) +
           scale_color_discrete(name='Region')+
           geom_point( aes(x=t,y=y))+facet_wrap(m~federal_units[s],scales="free",nrow=2)+
           geom_ribbon( aes(x=t,  ymin=`2.5%`,ymax=`97.5%`),alpha=0.2)+
           geom_line(aes(x=t, y=`50%`,colour=as.factor(federal_units[s]))))
  
}

get(paste('plot',1,sep=''))
get(paste('plot',2,sep=''))
get(paste('plot',3,sep=''))
get(paste('plot',4,sep=''))
get(paste('plot',5,sep=''))
get(paste('plot',6,sep=''))
get(paste('plot',7,sep=''))
get(paste('plot',8,sep=''))
get(paste('plot',9,sep=''))
