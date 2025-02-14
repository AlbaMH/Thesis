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

SARI_2122<-rbind(SARI_2021,mutate(SARI_2022, SEM_PRI=SEM_PRI+52,SEM_NOT=SEM_NOT+52))%>%
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
SARI_0_18_long<-filter(SARI_0_18, delay_week>0)%>%group_by(SEM_PRI,delay_week,SG_UF_NOT)%>%summarise(cases=n())
# Remove latest week (N_raw) as only 3 regions have observations.
N_raw<-max(SARI_0_18_long$SEM_PRI)
SARI_0_18_long<-filter(SARI_0_18_long, SEM_PRI<N_raw)
SARI_0_18_delay<-SARI_0_18_long%>%ungroup()%>%
  pivot_wider(id_cols=c(SEM_PRI,SG_UF_NOT),names_from = delay_week,values_from = cases)
# If entry is missing then zero cases were reported. 
SARI_0_18_delay[is.na(SARI_0_18_delay)]<-0

# Filter for second age group:
SARI_19_60<-filter(SARI_2122,age_group=="19-60")
# Make matrix with onset week as rows and delay weeks as columns
SARI_19_60_long<-filter(SARI_19_60, delay_week>0)%>%group_by(SEM_PRI,delay_week,SG_UF_NOT)%>%summarise(cases=n())
# Remove latest week (N_raw) as only 3 regions have observations.
N_raw<-max(SARI_19_60_long$SEM_PRI)
SARI_19_60_long<-filter(SARI_19_60_long, SEM_PRI<N_raw)
SARI_19_60_delay<-SARI_19_60_long%>%ungroup()%>%
  pivot_wider(id_cols=c(SEM_PRI,SG_UF_NOT),names_from = delay_week,values_from = cases)
# If entry is missing then zero cases were reported. 
SARI_19_60_delay[is.na(SARI_19_60_delay)]<-0

# Filter for third age group:
SARI_61<-filter(SARI_2122,age_group==">60")
# Make matrix with onset week as rows and delay weeks as columns
SARI_61_long<-filter(SARI_61, delay_week>0)%>%group_by(SEM_PRI,delay_week,SG_UF_NOT)%>%summarise(cases=n())
# Remove latest week (N_raw) as only 3 regions have observations.
N_raw<-max(SARI_61_long$SEM_PRI)
SARI_61_long<-filter(SARI_61_long, SEM_PRI<N_raw)
SARI_61_delay<-SARI_61_long%>%ungroup()%>%
  pivot_wider(id_cols=c(SEM_PRI,SG_UF_NOT),names_from = delay_week,values_from = cases)
# If entry is missing then zero cases were reported. 
SARI_61_delay[is.na(SARI_61_delay)]<-0

# Spatial names (federal_units) and number of regions (S).
federal_units<-sort(unique(SARI_19_60_delay$SG_UF_NOT))
S<-length(federal_units)
# Number of weeks in final data frame.
N<-N_raw-1
# Number of age groups.
A<-3
# Maximum Delay to consider.
D_max<-20
# Number of delays explicitly modelled in GDM.
D<-8


# Order attributes in arrays:
SARI_0_18_delay_ordered<-SARI_0_18_delay[,1:(D_max+2)]%>%mutate(id=paste(SEM_PRI,SG_UF_NOT,sep=''))
SARI_0_18_delay_ordered<-SARI_0_18_delay_ordered%>%arrange(SEM_PRI,SG_UF_NOT)
SARI_19_60_delay_ordered<-SARI_19_60_delay[,1:(D_max+2)]%>%mutate(id=paste(SEM_PRI,SG_UF_NOT,sep=''))
SARI_19_60_delay_ordered<-SARI_19_60_delay_ordered%>%arrange(SEM_PRI,SG_UF_NOT)
SARI_61_delay_ordered<-SARI_61_delay[,1:(D_max+2)]%>%mutate(id=paste(SEM_PRI,SG_UF_NOT,sep=''))
SARI_61_delay_ordered<-SARI_61_delay_ordered%>%arrange(SEM_PRI,SG_UF_NOT)

# Create blank data frame for all regions and weeks. 
SARI_0_18_full <- cbind(tibble(SEM_PRI=sort(rep(1:N,27)), SG_UF_NOT=rep(sort(federal_units),N)),
                        matrix(0, ncol=D_max,nrow=S*N))%>%mutate(check = NA)%>%mutate(id=paste(SEM_PRI,SG_UF_NOT,sep=''))
SARI_19_60_full <- cbind(tibble(SEM_PRI=sort(rep(1:N,27)), SG_UF_NOT=rep(sort(federal_units),N)),
                         matrix(0, ncol=D_max,nrow=S*N))%>%mutate(check = NA)%>%mutate(id=paste(SEM_PRI,SG_UF_NOT,sep=''))
SARI_61_full <- cbind(tibble(SEM_PRI=sort(rep(1:N,27)), SG_UF_NOT=rep(sort(federal_units),N)),
                      matrix(0, ncol=D_max,nrow=S*N))%>%mutate(check = NA)%>%mutate(id=paste(SEM_PRI,SG_UF_NOT,sep=''))

# Fill in available data:
SARI_0_18_full[SARI_0_18_full$id %in% SARI_0_18_delay_ordered$id ,3:(D_max+3)] <- SARI_0_18_delay_ordered[,3:(D_max+3)]
SARI_19_60_full[SARI_19_60_full$id %in% SARI_19_60_delay_ordered$id ,3:(D_max+3)] <- SARI_19_60_delay_ordered[,3:(D_max+3)]
SARI_61_full[SARI_61_full$id %in% SARI_61_delay_ordered$id ,3:(D_max+3)] <- SARI_61_delay_ordered[,3:(D_max+3)]


# Make an array including delay:
# SARI_..._array[REGION, TIME, DELAY]
SARI_0_18_array<-SARI_0_18_full[,3:(D_max+2)]%>%as.matrix()%>%array(dim=c(S,N,D_max))
SARI_19_60_array<-SARI_19_60_full[,3:(D_max+2)]%>%as.matrix()%>%array(dim=c(S,N,D_max))
SARI_61_array<-SARI_61_full[,3:(D_max+2)]%>%as.matrix()%>%array(dim=c(S,N,D_max))
# SARI_array[REGION, TIME, DELAY, AGE GROUP]
SARI_array<-array(abind(SARI_0_18_array,SARI_19_60_array,SARI_61_array),dim=c(S,N,D_max,A))
# Total SARI cases for each region.
total_SARI<-apply(SARI_array,c(1,4),sum)
SARI_delay_totals<-melt(total_SARI)
colnames(SARI_delay_totals)<-c('SG_UF_NOT',"AGE GROUP","totals")


# Matrix for partial SARI cases.
SARI_sea_z<-SARI_array[,,1:D_max,] #[REGION, TIME, DELAY, AGE GROUP]
# Matrix for total SARI cases.
SARI_sea_y<-apply(SARI_sea_z,c(1,2,4), sum) #[REGION, TIME, AGE GROUP]
# Dates of all weeks in data frame. 
SARI_dates<-as.Date(first_date)+(1:ncol(SARI_sea_y)-1)*7

# Make array for COVID-19 cases:
COVID_age<-filter(SARI_2122, delay_week>0, delay_week<21)%>%group_by(SEM_PRI,age_group,SARS2,SG_UF_NOT)%>%summarise(cases=n())%>%drop_na(age_group,SARS2)%>%
  pivot_wider(id_cols=c(SEM_PRI,SG_UF_NOT),names_from = age_group,values_from = cases)
# If entry is missing then zero cases were reported. 
COVID_age[is.na(COVID_age)]<-0

# Remove latest week (N_raw) as only 3 regions have observations.
COVID_long<-filter(COVID_age, SEM_PRI<N_raw)
# Add ID and order matrix.
COVID_long_ordered<-COVID_long%>%mutate(id=paste(SEM_PRI,SG_UF_NOT,sep=''))%>%arrange(SEM_PRI,SG_UF_NOT)

# Create blank data frame for all regions and time. 
COVID_full <- tibble(SEM_PRI=sort(rep(1:N,27)), SG_UF_NOT=rep(sort(federal_units),N),
                     `0-18` = 0, `19-60` = 0, `>60` = 0, check = NA)
COVID_full<-COVID_full%>%mutate(id=paste(SEM_PRI,SG_UF_NOT,sep=''))

# Fill in available data:
COVID_full[COVID_full$id %in% COVID_long_ordered$id ,3:6] <- COVID_long_ordered[,3:6]
COVID_wide<-COVID_full[,-c(1,2,6,7)]%>%as.matrix()%>%array(dim=c(S,N,A))

# Matrix for COVID cases.
SARI_sea_x<-COVID_wide # [REGION, TIME, AGE]

SARI_x<-t(apply(SARI_sea_x,c(1,2),sum))
SARI_z<-apply(SARI_sea_z,c(1,2,3),sum)
SARI_y<-t(apply(SARI_sea_y,c(1,2),sum))

# Melt arrays/matrix to create data frames:
# Data frame for COVID cases.
x_all<-melt(apply(SARI_sea_x,c(1,2),sum))
colnames(x_all)<-c("s","t","x")
x_all$s<-federal_units[x_all$s]
# Data frame for SARI cases.
y_all<-melt(apply(SARI_sea_y,c(1,2),sum))
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
# population<-my.table[,c(1,2,6)]
Data_delay_long<-as.tibble(Data_delay_sum%>%pivot_longer(!c(t,s), names_to = "delay", values_to = "count")%>%mutate(s=federal_units[s]))
Data_delay_pop<-full_join(Data_delay_long,population, by="s")
Data_delay_prop<-Data_delay_pop%>%filter(delay%in%c("prop_0","prop_1","prop_2","prop_3","prop_4"))
Data_delay_prop<-full_join(Data_delay_prop,totals[,c(1,2,dim(totals)[2])],by = join_by(t, s))
Data_delay_prop<-Data_delay_prop%>%mutate(total.per.capita=total/population)
# plot totals against time
totals_time<-ggplot(Data_delay_prop,aes(x=t,y=total))+
  geom_point(alpha=0.2)+facet_wrap(~s,scales="free")
  labs(x="Total counts (log)",y="Cumulative proportions (probit)",title="SARI: cumulative proportions")
# plot cumulative proportions against time
cumulative_time<-ggplot(Data_delay_prop,aes(x=t,y=probit(count), colour=delay))+
  geom_point(alpha=0.2)+facet_wrap(~s,scales="free")+geom_smooth(alpha=0.2)+
  labs(x="Total counts (log)",y="Cumulative proportions (probit)",title="SARI: cumulative proportions")
# plot cumulative proportions against totals
cumulative_regions<-ggplot(Data_delay_prop,aes(x=log(total),y=probit(count), colour=delay))+
  geom_point(alpha=0.2)+facet_wrap(~s,scales="free")+geom_smooth(alpha=0.2)+
  labs(x="Total counts (log)",y="Cumulative proportions (probit)",title="SARI: cumulative proportions")
# plot cumulative proportions against totals per capita
cumulative_delay<-ggplot(Data_delay_prop,aes(x=log(total),y=probit(count),colour=s))+
  geom_point(alpha=0.2)+facet_wrap(~delay,scales="free")+geom_smooth(alpha=0.2)+
  labs(x="Total counts (log)",y="Cumulative proportions (probit)",title="SARI: cumulative proportions")

Data_delay_sum$prop_0<-as.numeric(unlist((Data_wide[,3])/apply(Data_wide[,3:22],1, sum, na.rm=TRUE)))
Data_delay_sum$prop_1<-as.numeric(unlist((Data_wide[,4])/(apply(Data_wide[,3:22],1, sum, na.rm=TRUE)-Data_wide[,3])))
Data_delay_sum$prop_2<-as.numeric(unlist((Data_wide[,5])/(apply(Data_wide[,3:22],1, sum, na.rm=TRUE)-apply(Data_wide[,3:4],1, sum, na.rm=TRUE))))
Data_delay_sum$prop_3<-as.numeric(unlist((Data_wide[,6])/(apply(Data_wide[,3:22],1, sum, na.rm=TRUE)-apply(Data_wide[,3:5],1, sum, na.rm=TRUE))))
Data_delay_sum$prop_4<-as.numeric(unlist((Data_wide[,7])/(apply(Data_wide[,3:22],1, sum, na.rm=TRUE)-apply(Data_wide[,3:6],1, sum, na.rm=TRUE))))

Data_delay_long<-as.tibble(Data_delay_sum%>%pivot_longer(!c(t,s), names_to = "delay", values_to = "count")%>%mutate(s=federal_units[s]))
Data_delay_pop<-full_join(Data_delay_long,population, by="s")
Data_delay_prop<-Data_delay_pop%>%filter(delay%in%c("prop_0","prop_1","prop_2","prop_3","prop_4"))
Data_delay_prop<-full_join(Data_delay_prop,totals[,c(1,2,dim(totals)[2])],by = join_by(t, s))
Data_delay_prop<-Data_delay_prop%>%mutate(total.per.capita=total/population)
# plot cumulative proportions against totals
relative_regions<-ggplot(Data_delay_prop,aes(x=log(total),y=logit(count), colour=delay))+
  geom_point(alpha=0.2)+facet_wrap(~s,scales="free")+geom_smooth(alpha=0.2)+
  labs(x="Total counts (log)",y="Relative proportions (logit)",title="SARI: relative proportions")
# plot cumulative proportions against totals per capita
relative_delay<-ggplot(Data_delay_prop,aes(x=log(total),y=logit((count)), colour=s))+
  geom_point(alpha=0.2)+facet_wrap(~delay,scales="free")+geom_smooth(alpha=0.2)+
  labs(x="Total counts (log)",y="Relative proportions (logit)",title="SARI: relative proportions")


# Create function for cluster to carry out MCMC model:
SARI_sea_reg <- function(seed, SARI_x, SARI_z, Nowcast_list, population, S, R, D, D_max, n_knots, n_iter, n_burn, n_thin, n_chains){
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
  

  
  assign('dbetabin', dbetabin, envir = .GlobalEnv)
  assign('rbetabin', rbetabin, envir = .GlobalEnv)

  registerDistributions(list(dbetabin=list(
    BUGSdist='dbetabin(mu,phi,size)',discrete=TRUE)))
  
  # NIMBLE model code. 
  SARI_code <- nimbleCode({
    for(s in 1:S){
      # for(t in 1:N){
      # Beta-Binomial model for COVID cases given SARI cases.
      #  x[t,s] ~ dbetabin(mu[t,s],chi[s],y[s,t])
      #  }
      # Slope coefficient for reporting rate of COVID cases.
      #  omega[s] ~ dgamma(shape=10,rate=200)
      # Reporting rate for COVID cases.
      #  for(t in 1:M){
      #   pi[t,s] <- 1
      #  }
      #  for(t in (M+1):N){
      #   pi[t,s] <- exp(-omega[s]*(t-M))
      # Beta-Binomial model for censored COVID cases given eventual COVID cases.
      #   c[t,s] ~ dbetabin(pi[t,s],upsilon[s],x[t,s])
      # }
      for(t in 1:N){
        # Negative Binomial Model for total SARI cases.
        y[s,t] ~ dnegbin(prob=theta[s]/(theta[s]+lambda[t,s]),size=theta[s])  
        
        # Model for partial delayed SARI counts (first delay).
        z[s,t,1] ~ dbetabin(nu[t,1,s],phi[s,1],y[s,t])
      }
      for(d in 2:D){
        for(t in 1:obs_index[s,d]){
          # Model for partial delayed SARI counts.
          z[s,t,d] ~ dbetabin(nu[t,d,s],phi[s,d],y[s,t]-sum(z[s,t,1:(d-1)]))
        }
      }
      for(t in 1:N){
        # Linear predictor (Beta-Binomial means).
       for(d in 1:D){
          # Expected cumulative proportions.
          probit(p[t,d,s]) <- psi[d,s] + eta[t,s] #+ delta[s]*log((y[s,t]+1)/pop[s])
        }
        # Relative proportions (Beta-Binomial means).
        nu[t,1,s] <- p[t,1,s]
        for(d in 2:D){
          nu[t,d,s] <- (p[t,d,s]-p[t,d-1,s])/(1-p[t,d-1,s])
        }
        # Mean for total SARI cases .
        log(lambda[t,s]) <-  log(pop[s]) + zeta0[s] + zeta1[t,s] + xi[weeks[t],s] 
        # Linear predictor for COVID Beta-Binomial probability.
        #    logit(mu[t,s]) <- beta0[s] + beta1[t,s] + zeta1[t,s]*omicron[s] + beta_cov_0to18[1]*age_0to18[t,1,s] + beta_cov_0to18[2]*age_0to18[t,2,s] + beta_cov_0to18[3]*age_0to18[t,3,s] +  
        #     beta_cov_19to60[1]*age_19to60[t,1,s] + beta_cov_19to60[2]*age_19to60[t,2,s] + beta_cov_19to60[3]*age_19to60[t,3,s]
        
      }
      
      
      # Splines of time in cumulative proportion reported.
      Omega_eta[1:K_t2,1:K_t2,s] <- S_t2[1:K_t2,1:K_t2]/sigma_eta[s]^2
      kappa_eta[1:K_t2,s] ~ dmnorm(zeros[1:K_t2],Omega_eta[1:K_t2,1:K_t2,s])
      eta[1:N,s] <- X_t2[1:N,1:K_t2]%*%kappa_eta[1:K_t2,s]
      
      # Seasonal spline (for influenza) 
      Omega_xi[1:K_w,1:K_w,s]<-S_w[1:K_w,1:K_w]/sigma_xi[s]^2
      kappa_xi[1:K_w,s]~dmnorm(zeros[1:K_w],Omega_xi[1:K_w,1:K_w,s])
      xi[1:52,s]<-X_w[1:52,1:K_w]%*%kappa_xi[1:K_w,s]
      
      
      # Splines of time for proportion of COVID-19/SARI cases.
      #  Omega_beta[1:K_t3,1:K_t3,s] <- S_t3[1:K_t3,1:K_t3]/sigma_beta[s]^2
      #  kappa_beta[1:K_t3,s] ~ dmnorm(zeros[1:K_t3],Omega_beta[1:K_t3,1:K_t3,s])
      #  beta1[1:N,s] <- X_t3[1:N,1:K_t3]%*%kappa_beta[1:K_t3,s]
      
      # Splines of time for SARI cases.
      Omega_zeta[1:K_t2,1:K_t2,s] <- S_t2[1:K_t2,1:K_t2]/sigma_zeta[s]^2
      kappa_zeta[1:K_t2,s] ~ dmnorm(zeros[1:K_t2],Omega_zeta[1:K_t2,1:K_t2,s])
      zeta1[1:N,s] <- X_t2[1:N,1:K_t2]%*%kappa_zeta[1:K_t2,s]
      

      # Smoothing parameter priors.
      # sigma_beta[s] ~ T(dnorm(0,1),0,)
      sigma_zeta[s] ~ T(dnorm(0,1),0,)
      sigma_eta[s] ~ T(dnorm(0,1),0,) 
      sigma_xi[s] ~ T(dnorm(0,1),0,) 
      
      # Cumulative delay curves.
      psi[1,s] ~ dnorm(0,sd=10)
      for(d in 2:D){
        psi[d,s] ~ T(dnorm(psi[d-1,s],sd=10),psi[d-1,s],)
      }
      
      # Beta-Binomial dispersion parameters.
      for(d in 1:D){
        phi[s,d] ~ dgamma(2,0.02) 
      }
      
      # Intercepts
      # beta0[s] ~ dnorm(0,sd=5)
      zeta0[s]  ~ dnorm(-mean.sari[s],sd=10)
      
    #  omicron[s] ~ dnorm(0,sd=5) # Effect scale of SARI cases on proportion of COVID-19/SARI.
      chi[s] ~ dgamma(2,0.02) # Beta-Binomial dispersion parameters.
      upsilon[s] ~ dgamma(2,0.02) # Beta-Binomial dispersion parameters.
      theta[s] ~ dgamma(2,0.02) # Negative-Binomial dispersion parameters.
    }

  })
  
  # Constants for current nowcast date: 
  N_now<-Nowcast_list[seed] # week to nowcast up to.
  N<-N_now # number of weeks to model
  C<-N-D_max+1 #time up to SARI cases are fully observed
  M<-N-R #time up to COVID cases are fully observed
  A<-dim(SARI_z)[4] # number of age groups
  #  
  # Censor partial SARI counts.
  z<- SARI_z[1:S,(1):(N_now),,1:A]
  censored_z_age <- z
  for(a in 1:A){
    for(s in 1:S){
      a_z<-censored_z_age[s,,,a]
      a_z[outer(1:dim(z[s,,,a])[1], 0:(dim(z[s,,,a])[2]-1), FUN = "+") > N] <- NA
      censored_z_age[s,,,a]<-a_z
    }
  }
  
  # Censor the total SARI counts.
  censored_y<- apply(censored_z_age[1:S,,1:D_max,1:A],c(1,2),sum)
  censored_y_na<- apply(censored_z_age[1:S,,1:D_max,1:A],c(1,2,4),sum,na.rm=TRUE)
  # # Calculate age group covariates. 
  # age_group_prop<-censored_y_na[,,1:A]/array(rep(apply(censored_y_na,c(1,2),sum),A),dim=c(S,N,A))
  # age_group_prop[age_group_prop=="NaN"]<-0
  # # Scale age group covaraites.
  # age_group_prop_scaled <- age_group_prop%>%apply(3,function(x){(x-mean(x))/sd(x)})%>%
  #   array(dim=dim(age_group_prop))%>%
  #   aperm(c(2,1,3))
  # # Covaraites for cubic polynomial:
  # age_group_cubic<-array(NA, dim=c(N,3,S,2)) #time,cubic terms order,region,age group
  # age_group_cubic[,1,,] <- age_group_prop_scaled[,,1:2]
  # age_group_cubic[,2,,] <- age_group_prop_scaled[,,1:2]^2
  # age_group_cubic[,3,,] <- age_group_prop_scaled[,,1:2]^3
  
  
  # Calculate censored z totals.
  censored_z <- apply(censored_z_age,c(1,2,3),sum)
  
  # Maximum observed index for each region at each delay.
  obs_index<-matrix(nrow=D, ncol=S)
  for(d in 1:D){
    for(s in 1:S){
      obs_index[d,s]<-which(is.na(censored_z[s,,d])==TRUE)[1]-1
    }
  }
  
  # # Simulate COVID censoring.
  # censored_x<-SARI_x[1:S,(1):(N_now),1:A]
  # for(a in 1:A){
  #   for(s in 1:S){
  #     censored_x[s,(N-R+1):N,a]<-floor(censored_x[s,(N-R+1):N,a]*seq(from=0.95, to=0.25, length=R)) 
  #   }
  # }
  # # Censored COVID counts.
  # censored_x<-t(apply(censored_x,c(1,2),sum))
 
  
  
  # Cubic spline with shrinkage for time (mean SARI cases, cumulative SARI cases reported) %
  # yearly seasonal spline (for 52 weeks):
  week_index<-rep(1:52,9) #first date is 01-01-2013 - NEEDS UPDATING 
  blank_data2=tibble(y=rnorm(N,0,1),t=1:N,w=week_index[1:N]) 
  blank_jagam2=jagam(y~s(t,k=floor(N/20), bs='cs',pc=round(C/2))+
                      s(w,bs='cc',k=n_knots[2]),data=blank_data2,file='blank.jags',
                    knots=list(t=seq(1,N,length=floor(N/20)),w=seq(0,52,length=n_knots[2]))) # #


  # Constants
  SARI_constants <- list(N=N,
                         S=S,
                         D=D,
                         weeks=week_index[1:N],
                         pop=as.numeric(population$population),
                         #  age_0to18=array(age_group_cubic[,,,1],dim=c(N,3,S)),
                         #  age_19to60=array(age_group_cubic[,,,2],dim=c(N,3,S)),
                         mean.sari=floor(log(apply(censored_y[,1:(N-D_max)],1,mean))),
                         obs_index=t(obs_index),
                         K_w=dim(blank_jagam2$jags.data$S2)[1],
                         S_w=blank_jagam2$jags.data$S2,
                         K_t2=dim(blank_jagam2$jags.data$S1)[1],
                         S_t2=blank_jagam2$jags.data$S1)
  

  SARI_constants$X_t2 <- blank_jagam2$jags.data$X[,2:(SARI_constants$K_t2+1)]
  SARI_constants$X_w <- blank_jagam2$jags.data$X[,(SARI_constants$K_t2+2):(SARI_constants$K_t2+SARI_constants$K_w+1)]
  
  SARI_constants$zeros <- rep(0,max(SARI_constants$K_t3,SARI_constants$K_t2,SARI_constants$K_w))
  
  # Data:
  SARI_data <- list(y=censored_y[,1:N],
                    #  x=as.matrix(censored_x[1:N,],ncol=S),
                    #  c=as.matrix(censored_x[1:N,],ncol=S),
                    z=censored_z[,,1:D]) 
  #SARI_data$x[(M+1):N,] <- NA
  
  SARI_inits<-SARI_model<-SARI_compiled_model<- SARI_mcmc_config<-SARI_mcmc <- SARI_compiled_mcmc<-list()
  # Generate random initial values.
  SARI_inits <- list(
    # beta_cov_0to18=(rnorm(3,0,0.01)),
    # beta_cov_19to60=(rnorm(3,0,0.01)),
    # kappa_beta=matrix(rnorm(S*(SARI_constants$K_t3),0,0.01),ncol=S),
    kappa_xi=matrix(rnorm(S*(SARI_constants$K_w),0,0.1),ncol=S),
    kappa_zeta=matrix(rnorm(S*(SARI_constants$K_t2),0,0.01),ncol=S),
    kappa_eta=matrix(rnorm(S*(SARI_constants$K_t2),0,0.1),ncol=S),
    kappa_delta=matrix(rnorm(S*(SARI_constants$K_c),0,0.1),ncol=S),
    sigma_delta=runif(2,0,1),
    #  phi=matrix(abs(rnorm(S*D,30,10)),nrow=S,ncol=D),
    sigma_xi=runif(S,0,1),
    psi=matrix(sort(rnorm(D*S,0,1)),ncol=S, byrow = TRUE),
    sigma_eta=runif(S,0,1),
    # sigma_beta=runif(S,0,1),
    sigma_zeta=runif(S,0,1),
    #  omicron=rnorm(S,0,1),
    #  chi=runif(S,5,10),
    #  upsilon=runif(S,5,10),
    theta=abs(rnorm(S,50,10)),
    #  beta0=rnorm(S,0,1),
    zeta0=rnorm(S,-10,1),
    #   omega=runif(S,0,0.1),
    #   pi=exp(t(matrix(c(rep(NA,M*S),runif((N-M)*S,-1,0)),nrow=S))),
    #  x=matrix(NA,nrow=N,ncol=S),
    y=matrix(NA,nrow=S,ncol=N) )
  # Set initial values for y which are greater than the as-of-yet reported counts:
  for(t in 1:N){
    for(s in 1:S){
      if(is.na(censored_y[s,t])) SARI_inits$y[s,t]=sum(SARI_data$z[s,t,],rpois(1,median(SARI_data$y[s,]-rowSums(SARI_data$z[s,,]),na.rm=T)),na.rm=TRUE)
      #  if(is.na(SARI_data$x[t,s])) SARI_inits$x[t,s]=mean(c(SARI_data$c[t,s],SARI_inits$y[s,t],censored_y[s,t]),na.rm=T)
    }
  }
  # Build the model.
  SARI_model <- nimbleModel(SARI_code,SARI_constants,SARI_data,SARI_inits)
  # Compile the model.
  SARI_compiled_model <- compileNimble(SARI_model)
  
  # Set up the MCMC.
  SARI_mcmc_config <- configureMCMC(SARI_model,monitors=c("lambda","zeta0","y",#"phi","chi",
                                                          "psi","theta",
                                                          "sigma_eta","sigma_zeta", "xi",
                                                          "zeta1",
                                                          #"kappa_beta", "beta_cov_0to18","beta_cov_19to60"
                                                          #'pi','beta1',"beta0","omega","x","sigma_beta",
                                                          "eta",
                                                          "kappa_zeta","kappa_eta"
  ),
  useConjugacy = FALSE)
  
  SARI_mcmc_config$removeSamplers('kappa_zeta','sigma_zeta',#'kappa_beta','sigma_beta',
                                  'kappa_xi','sigma_xi', 'psi',
                                  #'omega','omicron','beta0',
                                  'kappa_eta','sigma_eta','zeta0')
  for(s in 1:S){
    SARI_mcmc_config$addSampler(target=c(paste('kappa_zeta[1:',(SARI_constants$K_t2),', ',s,']', sep=''),paste('sigma_zeta[',s,']',sep=''),paste('zeta0[',s,']',sep=''),
                                         paste('kappa_xi[1:',(SARI_constants$K_w),', ',s,']', sep=''),paste('sigma_xi[',s,']',sep='')
                                         #   paste('kappa_beta[1:',(SARI_constants$K_t3),', ',s,']', sep=''),paste('sigma_beta[',s,']',sep=''),paste('beta0[',s,']',sep=''),
                                         #   paste('omicron[',s,']',sep=''),paste('omega[',s,']',sep='')
    ),type='AF_slice') 
    SARI_mcmc_config$addSampler(target=c( paste('kappa_eta[1:',(SARI_constants$K_t2),', ',s,']', sep=''),paste('sigma_eta[',s,']',sep=''),
                                         paste('psi[1:',D,', ',s,']', sep='')),type='AF_slice')
  }
  
  #SARI_mcmc_config$removeSamplers(target=c(paste('beta_cov_0to18[1:3]',sep=''),paste('beta_cov_19to60[1:3]',sep='')))
  # SARI_mcmc_config$addSampler(target=c(paste('beta_cov_0to18[1:3]',sep=''),paste('beta_cov_19to60[1:3]',sep='')),type='AF_slice')
  # Build MCMC. 
  SARI_mcmc<- buildMCMC(SARI_mcmc_config)
  # Compile the MCMC.
  SARI_compiled_mcmc <- compileNimble(SARI_mcmc,project=SARI_model)
  # Run MCMC. 
  SARI_cluster<- runMCMC(SARI_compiled_mcmc,niter=n_iter,nburnin=n_burn,thin=n_thin,inits=SARI_inits,nchains=n_chains,samplesAsCodaMCMC = TRUE)
  
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
n_knots<-c(15,6,8,10) # knots for time, time-delay interaction and delay respectively 
# Set maximum COVID-19 lab report delay.
R<-30  


# List of weeks to perform the nowcasts for rolling predictions experiment.
Nowcast_list<-floor(seq(from=100, to=N-D_max, length=4))
# Dates for the nowcasts.
Nowcast_dates<- as.Date(first_date)+(Nowcast_list-1)*7

# Number of available computer cores. 
n_cores<-min(detectCores(),length(Nowcast_list))

# Survivor model with age covariates. 
# Make Cluster for MCMC. 
this_cluster_reg_age<- makeCluster(n_cores)

# Run Cluster.
time_sea_delta <- system.time({
  SARI_sea_output_delta <-  parLapply(cl = this_cluster_reg_age, X = 1:length(Nowcast_list), 
                                      fun = SARI_sea_reg,
                                      SARI_x= as.array(SARI_sea_x[1:S,,]),
                                      SARI_z=as.array(SARI_sea_z[1:S,,,]),
                                      Nowcast_list=Nowcast_list,
                                      population=population,
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

# Stop Cluster.
stopCluster(this_cluster_reg_age)

# Save survivor model output. 
save(SARI_sea_output_delta,time_sea_delta, file='survivor_sea27.RData')

# Check PSRF of model parameters.
sari_parameter_group=c("zeta0","xi","lambda",
                       "psi","theta",
                       "sigma_eta","sigma_zeta", "zeta1",
                       "kappa_zeta","kappa_eta")
sari_psrf<-list()
sari_mean_psrf<-sari_max_psrf<-sari_leq_105 <-sari_leq_110<-array(NA,dim=c(length(sari_parameter_group),length(Nowcast_list)))
for(k in 1:length(sari_parameter_group)){
  sari_psrf[[k]]<-list()  
  for(j in 1:length(Nowcast_list)){
    sari_parameter_names <- colnames((as.matrix(SARI_sea_output_delta[[j]])))
    sari_index <- which(startsWith(sari_parameter_names,sari_parameter_group[k]))
    sari_psrf[[k]][[j]]<- sapply(sari_index,function(v)gelman.diag(SARI_sea_output_delta[[j]][,v],autoburnin = FALSE,transform = TRUE,multivariate = FALSE)$psrf[,1])
    sari_mean_psrf[k,j] <- mean(sari_psrf[[k]][[j]])
    sari_max_psrf[k,j] <- max(sari_psrf[[k]][[j]])
    sari_leq_105[k,j] <- mean(sari_psrf[[k]][[j]]<=1.05)
    sari_leq_110[k,j] <- mean(sari_psrf[[k]][[j]][j]<=1.1)
  }
}
sari_MAX_psrf<-cbind(sari_parameter_group,apply(sari_max_psrf,1,max))
sari_MEAN_psrf<-cbind(sari_parameter_group,apply(sari_mean_psrf,1,mean))


plot(SARI_sea_output_delta[[1]][,"lambda[1, 1]"],main="lambda[1, 1]")
plot(SARI_sea_output_delta[[2]][,"xi[10, 1]"],main="xi[10, 1]")
plot(SARI_sea_output_delta[[2]][,"xi[1, 1]"],main="xi[1, 1]")

#### Plots #####
library(geobr)
states<-read_state(year=2013)
n_sim<-dim(SARI_sea_output_delta[[1]]$chain1)[1]*2
SARI_Survivor_sea_x<-list()
SARI_combined_Survivor_sea_x<-list()
surv_sea_mu_x<-covid_surv_sea_quant_x<-sari_surv_sea_quant_x<-list()
SARI_pi_samples_surv_age<-SARI_mu_samples_surv_age<-SARI_y_samples_surv_age<-SARI_chi_surv_age<-SARI_x_samples_surv_age<-list()


for(j in length(Nowcast_list):1){
  # Survivor age model
  SARI_Survivor_sea_x[[j]]<-as.mcmc.list(SARI_sea_output_delta[[j]])
  SARI_combined_Survivor_sea_x[[j]] <- as_tibble(do.call('rbind',SARI_Survivor_sea_x[[j]]))
}

lambda_plot<-ggplot()+
  geom_point(data=y_all,aes(x=t,y=y))+
  facet_wrap(~s, scales = 'free')+
  theme_minimal()+labs(y="Mean SARI Hospitalisations", x="Time", title="CLR Regional Model")+
  scale_colour_discrete(name='Nowcast Date')+scale_fill_discrete(name='Nowcast Date')

eta_plot<-ggplot()+
  facet_wrap(~s, scales = 'free')+
  theme_minimal()+labs(y="Delay temporal effect", x="Time", title="CLR Regional Model")+
  scale_colour_discrete(name='Nowcast Date')+scale_fill_discrete(name='Nowcast Date')


y_plot<-ggplot()+
  geom_point(data=y_all,aes(x=t,y=y))+ 
  facet_wrap(~s, scales = 'free')+
  theme_minimal()+labs(y="SARI Hospitalisations", x="Time", title="CLR Regional Model")+#+geom_vline(xintercept = 51)
  scale_colour_discrete(name='Nowcast Date')+scale_fill_discrete(name='Nowcast Date')

delta_plot<-ggplot()+
  theme_minimal()+labs(y="Effect size (iprobit scale)", x="lambda (expected mean SARI counts)", title="Effect of SARI incidence on delay",
                       color = "Nowcast date") + theme(legend.position = "bottom")+
  facet_wrap(~s,scales='free')+
  scale_colour_discrete(name='Nowcast Date')+scale_fill_discrete(name='Nowcast Date')

delta_plot2<-ggplot()+
  theme_minimal()+labs(y="Effect size (probit scale)", x="lambda (expected mean SARI counts)", title="Effect of SARI incidence on delay",
                       color = "Nowcast date") + theme(legend.position = "bottom")+
  facet_wrap(~s,scales='free')+
  scale_colour_discrete(name='Nowcast Date')+scale_fill_discrete(name='Nowcast Date')

omega_plot<-ggplot()+
  facet_wrap(~s, scales = 'free')+
  theme_minimal()+labs(y="zeta", x="Time", title="Temporal logspline for mean SARI cases")+
  scale_colour_discrete(name='Nowcast Date')+scale_fill_discrete(name='Nowcast Date')

zeta_plot<-ggplot()+
  facet_wrap(~s, scales = 'free')+
  theme_minimal()+labs(y="zeta", x="Time", title="Temporal logspline for mean SARI cases")

xi_plot<-ggplot()+
  facet_wrap(~s, scales = 'free')+
  theme_minimal()+labs(y="xi", x="Weeks", title="Seasonal logspline for mean SARI cases")

psi_plot<-ggplot()+
  facet_wrap(~s, scales = 'free')+
  theme_minimal()+labs(y="psi", x="delay", title="Delay trend on cumulative proportions")

SARI_y_samples<-list()
lambda<-lambda_raw<-list()
SARI_pred<-list()
SARI_logspline_quantiles<-list()
theta<-delta<-zeta<-xi<-psi <- list()
delta<-delta_cubic<-delta_quantiles<-omega<-list()
index_delta<-list()
eta<-list()
p_raw<-eta_raw<-psi_raw<-list()
# Explore output
n_sim <- dim(SARI_combined_Survivor_sea_x[[1]])[1]
rbetabin=nimbleFunction(run=function(n=integer(0),mu=double(0),phi=double(0),size=double(0)){
  pi=rbeta(n,mu*phi,(1-mu)*phi)
  returnType(double(0))
  return(rbinom(n,size,pi))
})
week_index<-rep(1:52,9)

# plot output
zeta_raw<-xi_raw<-xi_plus_zeta<-list()
y_pop_delta<-y_pop_delta_long<-list()
p<-list()
for(j in length(Nowcast_list):1){
# Plot Model Variables 
  lambda[[j]] <- select(SARI_combined_Survivor_sea_x[[j]],starts_with('lambda'))%>%as.matrix()%>%array(dim=c(n_sim,Nowcast_list[j],S))%>%
    apply(c(2,3),quantile,c(0.025,0.5,0.975))%>%melt(varnames=c('quantile','t','s'),value.name='y')%>%#mutate(t=t+Nowcast_list[j]-W)%>%
    spread(quantile,y)%>%mutate(Nowcast_date=paste(Nowcast_list[j]),s=federal_units[s])
  
  lambda_plot<-lambda_plot+
    geom_line(data=lambda[[j]], aes(x=t, y=`50%`, colour=Nowcast_date))+
    geom_ribbon(data=lambda[[j]], aes(x=t, ymin=`2.5%`, ymax=`97.5%`,fill=Nowcast_date), alpha=0.2)
  
  
  SARI_logspline_quantiles[[j]]<-select(SARI_combined_Survivor_sea_x[[j]],starts_with('y'))%>%as.matrix()%>%array(dim=c(n_sim,S,Nowcast_list[j]))%>%
    apply(c(2,3),quantile,c(0.025,0.5,0.975))%>%melt(varnames=c('quantile','s','t'),value.name='y')%>%#mutate(t=t+Nowcast_list[j]-W)%>%
    spread(quantile,y)%>%mutate(Nowcast_date=paste(Nowcast_list[j]),s=federal_units[s])
  
  y_plot<-y_plot+
    geom_line(data=SARI_logspline_quantiles[[j]], aes(x=t, y=`50%`, colour=Nowcast_date))+
    geom_ribbon(data=SARI_logspline_quantiles[[j]], aes(x=t, ymin=`2.5%`, ymax=`97.5%`,fill=Nowcast_date))
  
  
  # Plot Zeta Variables 
  zeta[[j]] <- select(SARI_combined_Survivor_sea_x[[j]],starts_with('zeta1'))%>%as.matrix()%>%array(dim=c(n_sim,Nowcast_list[j],S))%>%
    apply(c(2,3),quantile,c(0.025,0.5,0.975))%>%melt(varnames=c('quantile','t','s'),value.name='y')%>%#mutate(t=t+Nowcast_list[j]-W)%>%
    spread(quantile,y)%>%mutate(Nowcast_date=paste(Nowcast_list[j]),s=federal_units[s])
  
  zeta_plot<-zeta_plot+
    geom_line(data=zeta[[j]], aes(x=t, y=`50%`, colour=Nowcast_date))+
    geom_ribbon(data=zeta[[j]], aes(x=t, ymin=`2.5%`, ymax=`97.5%`,fill=Nowcast_date), alpha=0.2)
  

  xi[[j]] <- select(SARI_combined_Survivor_sea_x[[j]],starts_with('xi'))%>%as.matrix()%>%array(dim=c(n_sim,52,S))%>%
    apply(c(2,3),quantile,c(0.025,0.5,0.975))%>%melt(varnames=c('quantile','w','s'),value.name='y')%>%#mutate(t=t+Nowcast_list[j]-W)%>%
    spread(quantile,y)%>%mutate(Nowcast_date=paste(Nowcast_list[j]),s=federal_units[s])
  
  xi_plot<-xi_plot+
    geom_line(data=xi[[j]], aes(x=w, y=`50%`, colour=Nowcast_date))+
    geom_ribbon(data=xi[[j]], aes(x=w, ymin=`2.5%`, ymax=`97.5%`,fill=Nowcast_date), alpha=0.2)
  
  #omega[[j]] <- select(SARI_combined_Survivor_sea_x[[j]],starts_with('omega'))%>%as.matrix()%>%array(dim=c(n_sim,W,S))%>%
  #  apply(c(2,3),quantile,c(0.025,0.5,0.975))%>%melt(varnames=c('quantile','t','s'),value.name='y')%>%mutate(t=t+Nowcast_list[j]-W)%>%
  # spread(quantile,y)%>%mutate(Nowcast_date=paste(Nowcast_list[j]),s=federal_units[s])
  #  omega_plot<-omega_plot+
  #    geom_line(data=omega[[j]], aes(x=t, y=`50%`, colour=Nowcast_date))+
  #   geom_ribbon(data=omega[[j]], aes(x=t, ymin=`2.5%`, ymax=`97.5%`,fill=Nowcast_date), alpha=0.2)
  eta[[j]] <- select(SARI_combined_Survivor_sea_x[[j]],starts_with('eta'))%>%as.matrix()%>%array(dim=c(n_sim,Nowcast_list[j],S))%>%
    apply(c(2,3),quantile,c(0.025,0.5,0.975))%>%melt(varnames=c('quantile','t','s'),value.name='y')%>%#mutate(t=t+Nowcast_list[j]-W)%>%
    spread(quantile,y)%>%mutate(Nowcast_date=paste(Nowcast_list[j]),s=federal_units[s])
  
  eta_plot<-eta_plot+
    geom_line(data=eta[[j]], aes(x=t, y=`50%`, colour=Nowcast_date))+
    geom_ribbon(data=eta[[j]], aes(x=t, ymin=`2.5%`, ymax=`97.5%`,fill=Nowcast_date), alpha=0.2)
  
  
  psi[[j]] <- select(SARI_combined_Survivor_sea_x[[j]],starts_with('psi'))%>%as.matrix()%>%array(dim=c(n_sim,D,S))%>%
    apply(c(2,3),quantile,c(0.025,0.5,0.975))%>%melt(varnames=c('quantile','d','s'),value.name='y')%>%#mutate(t=t+Nowcast_list[j]-W)%>%
    spread(quantile,y)%>%mutate(Nowcast_date=paste(Nowcast_list[j]),s=federal_units[s])
  
  psi_plot<-psi_plot+
    geom_line(data=psi[[j]], aes(x=d, y=`50%`, colour=Nowcast_date))+
    geom_ribbon(data=psi[[j]], aes(x=d, ymin=`2.5%`, ymax=`97.5%`,fill=Nowcast_date), alpha=0.2)
  
  eta_raw[[j]] <- select(SARI_combined_Survivor_sea_x[[j]],starts_with('eta'))%>%as.matrix()%>%array(dim=c(n_sim,Nowcast_list[j],S))
  psi_raw[[j]] <- select(SARI_combined_Survivor_sea_x[[j]],starts_with('psi'))%>%as.matrix()%>%array(dim=c(n_sim,D,S))
  p_raw[[j]]<-array(dim=c(n_sim,Nowcast_list[j],D,S))
   for(t in 1:(Nowcast_list[j])){
    for(s in 1:S){
      for(d in 1:D){
      p_raw[[j]][,t,d,s]<-eta_raw[[j]][,t,s]+psi_raw[[j]][,d,s] 
      }
    }
  }
  p[[j]]<-iprobit(p_raw[[j]])%>%array(dim=c(n_sim,Nowcast_list[j],D,S))%>%
    apply(c(2,3,4),quantile,c(0.025,0.5,0.975))%>%melt(varnames=c('quantile','t','d','s'),value.name='y')%>%#mutate(t=t+Nowcast_list[j]-W)%>%
    spread(quantile,y)%>%mutate(Nowcast_date=paste(Nowcast_list[j]),s=federal_units[s])
}

# SAVE PLOTS
lambda_plot
y_plot
eta_plot
zeta_plot
xi_plot
psi_plot
#ggsave(y_plot,file='Plots/pi_surv_SARI.pdf',width=10,height=7)

delay_name<-unique(Data_delay_prop$delay)
p_long<-p[[4]][,c(1,2,3,5)]%>%mutate(delay=delay_name[d])
check_residuals<-inner_join(p_long,Data_delay_prop,by=c('t','s','delay'))%>%mutate(residual=count-`50%`)

ggplot(filter(check_residuals,s%in%federal_units[c(16,19,21,22,23)], delay=='prop_0'),aes(x=t,y=probit(residual), colour=delay))+
  geom_point(alpha=0.2)+facet_wrap(~s,scales="free")+geom_smooth(alpha=0.2)+
  labs(x="Time",y="Cumulative proportions residuals (probit)",title="SARI: cumulative proportions")

ggplot(filter(check_residuals,s%in%federal_units[c(16,19,21,22,23)]),aes(x=log(total),y=probit(residual), colour=delay))+
  geom_point(alpha=0.2)+facet_wrap(~s,scales="free")+geom_smooth(alpha=0.2)+
  labs(x="Total counts (log)",y="Cumulative proportions residuals (probit)",title="SARI: cumulative proportions")


ggplot(filter(check_residuals,s%in%federal_units[c(5,6,7,8,9)], delay=='prop_0'),aes(x=t,y=probit(residual), colour=delay))+
  geom_point(alpha=0.2)+facet_wrap(~s,scales="free")+geom_smooth(alpha=0.2)+
  labs(x="Time",y="Cumulative proportions residuals (probit)",title="SARI: cumulative proportions")

ggplot(filter(check_residuals,s%in%federal_units[c(5,6,7,8,9)]),aes(x=log(total),y=probit(residual), colour=delay))+
  geom_point(alpha=0.2)+facet_wrap(~s,scales="free")+geom_smooth(alpha=0.2)+
  labs(x="Total counts (log)",y="Cumulative proportions residuals (probit)",title="SARI: cumulative proportions")


# plots of delay structure and prediction example:
j=1
N_now<-Nowcast_list[j] # week to nowcast up to.
N<-N_now # number of weeks to model
C<-N-D_max+1 #time up to SARI cases are fully observed
M<-N-R #time up to COVID cases are fully observed
SARI_z=as.array(SARI_sea_z[1:S,,,])
A<-dim(SARI_z)[4] # number of age groups
#  
# Censor partial SARI counts.
z<- SARI_z[1:S,(1):(N_now),,1:A]
censored_z_age <- z
for(a in 1:A){
  for(s in 1:S){
    a_z<-censored_z_age[s,,,a]
    a_z[outer(1:dim(z[s,,,a])[1], 0:(dim(z[s,,,a])[2]-1), FUN = "+") > N] <- NA
    censored_z_age[s,,,a]<-a_z
  }
}

# Censor the total SARI counts.
censored_y_na<- apply(censored_z_age[1:S,,1:D_max,1:A],c(1,2),sum,na.rm=TRUE)
censored_z <- apply(censored_z_age,c(1,2,3),sum)

censored_z_brazil<-apply(censored_z,c(2,3),sum)
censored_z_melt_cum<-censored_z_brazil%>%melt(value.name = 'z', varnames=c('t','d'))%>%mutate(d=as.factor(d))%>% 
  group_by(t) %>% 
  arrange(d) %>% 
  mutate(CumulativeSum = cumsum(z))

SARI_logspline_all<-select(SARI_combined_Survivor_sea_x[[j]],starts_with('y'))%>%as.matrix()%>%array(dim=c(n_sim,S,Nowcast_list[j]))%>%
  apply(c(1,3),sum)%>%apply(c(2),quantile,c(0.025,0.5,0.975))%>%melt(varnames=c('quantile','t'),value.name='y')%>%#mutate(t=t+Nowcast_list[j]-W)%>%
  spread(quantile,y)%>%mutate(Nowcast_date=paste(Nowcast_list[j]))

y_quantiles<-SARI_logspline_all
y_all_brazil<-y_all%>%group_by(t)%>%summarise(y=sum(y))
y_all_brazil_d<-y_all_brazil%>%mutate(d='Total',z=y,y=NULL)
censored_y_reg<-rowSums(t(censored_y_na), na.rm=TRUE)
censored_y_brazil<-tibble(y=censored_y_reg,t=1:N)
censored_z_full<-full_join(censored_z_melt_cum,y_all_brazil_d,by=c('t','z','d'))

# plot of SARI delay structure 

z_plot<-ggplot()+
  geom_point(data=filter(y_all_brazil,t<101),aes(x=as.Date(first_date)+(t)*7,y=y))+ 
  theme_minimal()+labs(y="Cumulative count reported", x=NULL, title="Severe Acute Respiratory Ilness (SARI) hospitalisations",subtitle = "Brazil")+#+geom_vline(xintercept = 51)
  geom_line(data=censored_z_melt_cum, aes(x=as.Date(first_date)+(t)*7, y=CumulativeSum, colour=d))+#, linetype = 'dashed')+
  scale_colour_discrete(name='Delay (weeks)')+
  scale_x_date(limits = c(as.Date("2022-04-15"),as.Date("2022-12-25")),breaks = "1 month", date_labels =  "%b %Y")+
  scale_y_continuous(limits = c(0,15000))+theme( legend.position="bottom")+
  guides(color = guide_legend(nrow = 2, byrow=TRUE))
z_plot
ggsave(z_plot,file='SARI_z_plot.pdf',width=9,height=5)

# plot of SARI prediction example
y_plot<-ggplot()+
  guides(fill="none", colour="none")+
  scale_linetype(name="Counts")+
  geom_point(data=filter(y_all_brazil,t<101),aes(x=as.Date(first_date)+(t)*7,y=y))+ 
  theme_minimal()+labs(y="Hospitalisations", x=NULL, title="Severe Acute Respiratory Ilness (SARI) hospitalisations",subtitle = "Brazil")+#+geom_vline(xintercept = 51)
  geom_line(data=censored_y_brazil, aes(x=as.Date(first_date)+(t)*7, y=y, linetype='Reported'))+#, linetype = 'dashed')+
  scale_colour_discrete(name='Nowcast Date')+scale_fill_discrete(name='Nowcast Date')+
  geom_line(data=y_quantiles, aes(x=as.Date(first_date)+(t)*7, y=`50%`, colour=Nowcast_date, linetype='Predicted'))+
  scale_x_date(limits = c(as.Date("2022-04-15"),as.Date("2022-12-25")),breaks = "1 month", date_labels =  "%b %Y")+
  scale_y_continuous(limits = c(0,15000))+
  geom_ribbon(data=y_quantiles, aes(x=as.Date(first_date)+(t)*7, ymin=`2.5%`, ymax=`97.5%`,fill=Nowcast_date),alpha=0.5)
y_plot
ggsave(y_plot,file='SARI_y_plot_seasonal.pdf',width=9,height=5)


