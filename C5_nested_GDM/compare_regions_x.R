# Compare are models for regional COVID-19 and SARI hospitalizations in Brazil 
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
library(ggh4x)

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
# 
# #Read in data..
# INFLUD21 <- read_delim("INFLUD21-30-01-2023.csv", delim = ";", escape_double = FALSE, trim_ws = TRUE)
# 
# INFLUD22 <- read_delim("INFLUD22-30-01-2023.csv", 
#                        delim = ";", escape_double = FALSE, 
#                        col_types = cols(DT_NOTIFIC = col_date(format = "%d/%m/%Y"), 
#                                         DT_SIN_PRI = col_date(format = "%d/%m/%Y")), 
#                        trim_ws = TRUE)
# 
# ##### DATA CLEANING #####
# # Select variable of interest from data set. 
# col_ant<-which(colnames(INFLUD21)=="AN_SARS2")
# col_pcr<-which(colnames(INFLUD21)=="PCR_SARS2")
# col_results<-which(colnames(INFLUD21)=="PCR_RESUL")
# col_hospital<-which(colnames(INFLUD21)=="DT_INTERNA")
# col_collection<-which(colnames(INFLUD21)=="DT_COLETA")
# col_an_result<-which(colnames(INFLUD21)=="DT_RES_AN")
# col_pcr_result<-which(colnames(INFLUD21)=="DT_PCR")
# col_notif_date<-which(colnames(INFLUD21)=="DT_DIGITA")
# # Create data frame for 2021.
# SARI_2021<- INFLUD21[,c(1:16, col_hospital, col_collection, col_results, col_pcr_result, col_notif_date, col_pcr, col_an_result, col_ant)]%>%
#   mutate(DT_NOTIFIC=as.Date((DT_NOTIFIC),format="%d/%m/%Y"),DT_DIGITA=as.Date((DT_DIGITA),format="%d/%m/%Y"),DT_SIN_PRI=as.Date((DT_SIN_PRI),format="%d/%m/%Y"))
# dates_2021<-sort(as.Date(unique(SARI_2021$DT_SIN_PRI),format="%d/%m/%Y"))
# # Create data frame for 2022.
# SARI_2022<- INFLUD22[,c(1:16, col_hospital, col_collection, col_results, col_pcr_result, col_notif_date, col_pcr, col_an_result, col_ant)]%>%
#   mutate(DT_NOTIFIC=as.Date((DT_NOTIFIC),format="%d/%m/%Y"),DT_DIGITA=as.Date((DT_DIGITA),format="%d/%m/%Y"),DT_SIN_PRI=as.Date((DT_SIN_PRI),format="%d/%m/%Y"))
# dates_2022<-sort(as.Date(unique(SARI_2022$DT_SIN_PRI),format="%d/%m/%Y"))
# # Create data frame for 2021 and 2022. 
# SARI_2122<-rbind(SARI_2021,mutate(SARI_2022, SEM_PRI=SEM_PRI+52,SEM_NOT=SEM_NOT+52))%>%
#   mutate(delay=difftime(DT_DIGITA,DT_SIN_PRI, unit=c("days")),delay_noti=difftime(DT_NOTIFIC,DT_SIN_PRI, unit=c("days")))%>%
#   mutate(DT_COLETA=as.Date((DT_COLETA),format="%d/%m/%Y"),DT_PCR=as.Date((DT_PCR),format="%d/%m/%Y"))%>%
#   mutate(pcr_delay=difftime(DT_PCR,DT_COLETA, unit=c("days")))%>%mutate(delay_diffweeks=as.numeric(difftime(DT_DIGITA,DT_SIN_PRI, unit=c("weeks"))))
# 
# 
# # delay of zero happened in first week and 1 happened in second week
# SARI_2122$delay_week<-ceiling(SARI_2122$delay_diffweeks+1/7)
# 
# # Make matrix with onset week as rows and delay weeks as columns
# SARI_delay_long<-filter(SARI_2122, delay_week>0)%>%group_by(SEM_PRI,delay_week,SG_UF_NOT)%>%summarise(cases=n())
# 
# SARI_delay<-SARI_delay_long%>%ungroup()%>%
#   pivot_wider(id_cols=c(SEM_PRI,SG_UF_NOT),names_from = delay_week,values_from = cases)
# SARI_delay[is.na(SARI_delay)]<-0
# 
# SARI_delay_totals<-SARI_delay
# SARI_delay_totals$total<-rowSums(SARI_delay[,-(1:2)], na.rm=T)
# 
# N_raw<-max(SARI_delay$SEM_PRI)
# # remove latest week as only 3 regions have observations 
# SARI_delay<-filter(SARI_delay, SEM_PRI<N_raw)
# N<-max(SARI_delay$SEM_PRI)
# D_N<-dim(SARI_delay)[2]-2
# 
# 
# SARI_delay_ordered<-SARI_delay
# SARI_delay_ordered<-SARI_delay_ordered[order(SARI_delay_ordered[,2]),]
# SARI_delay_ordered<-SARI_delay_ordered[order(SARI_delay_ordered[,1]),]
# # Save spatial constants
# federal_units<-sort(unique(SARI_delay_ordered$SG_UF_NOT))
# S<-length(federal_units)
# 
# # Make an array including delay:
# # SARI_array[REGION, TIME, DELAY]
# SARI_array<-SARI_delay_ordered[,-(1:2)]%>%as.matrix()%>%array(dim=c(S,N,D_N))
# 
# 
# 
# # total SARI cases for each region
# total_SARI<-apply(SARI_array,c(1,2),sum)
# SARI_delay_totals<-melt(total_SARI)
# colnames(SARI_delay_totals)<-c('SG_UF_NOT',"SARI","totals")
# 
# 
# # set D_max
# D_max<-20
# D<-8
# 
# # Data frame for partial SARI hospitalisations:
# SARI_z<-SARI_array[,,1:(D_max)]
# # Data frame for total SARI hospitalisations:
# SARI_y<-apply(SARI_z,c(1,2),sum)
# 
# 
# # Data frame for COVID hospitalisations:
# 
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
# SARI_counts_total_21<-SARI_2021%>%group_by(SEM_PRI)%>%
#   summarise(SARI=n(), PCR=sum(PCR_SARS2 == 1,na.rm = TRUE ), 
#             AN=sum(AN_SARS2 == 1,na.rm = TRUE),
#             SARS2=sum((PCR_SARS2==1)|(AN_SARS2==1),na.rm = TRUE))%>%
#   mutate(nu=SARS2/SARI)
# SARI_counts_total_21<-SARI_counts_total_21%>%mutate( WEEK=SEM_PRI)
# 
# SARI_2022<- INFLUD22[,c(1:16,col_results,col_pcr, col_ant)]
# dates_2022<-sort(as.Date(unique(SARI_2022$DT_SIN_PRI),format="%d/%m/%Y"))
# 
# SARI_counts_22<-SARI_2022%>%group_by(SEM_PRI,SG_UF_NOT)%>%summarise(SARI=n(), PCR=sum(PCR_SARS2 == 1,na.rm = TRUE ),
#                                                                     AN=sum(AN_SARS2 == 1,na.rm = TRUE ),
#                                                                     SARS2=sum((PCR_SARS2==1)|(AN_SARS2==1),na.rm = TRUE))%>%
#   mutate(nu=SARS2/SARI)
# SARI_counts_22<-SARI_counts_22%>%mutate( WEEK=SEM_PRI+52)
# 
# SARI_counts<-rbind(SARI_counts_21,SARI_counts_22)
# 
# # Earliest date we have data for.
# first_date<-min(dates_2021)
# 
# SARI_counts_wide <- SARI_counts %>% 
#   select(SG_UF_NOT,SARS2,WEEK) %>%
#   pivot_wider(names_from=SG_UF_NOT, values_from=SARS2) 
# 
# # Data frame for total COVID-positive SARI hospitalisations:
# SARI_x<-SARI_counts_wide[-N_raw,-(1:2)]


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
R<-30 # maximum COVID-19 lab report delay 
W<-60 #moving window size

# List of dates to preform the nowcasts for rolling predictions experiment.
Nowcast_list<-floor(seq(from=W, to=N-D_max, length=16))
# Dates for the nowcasts.
Nowcast_dates<-as.Date(first_date)+(Nowcast_list-1)*7

# Observed COVID-19 hospitalisations.
x_all<-melt(as.matrix(SARI_x))
colnames(x_all)<-c("t","s","x")

# Observed SARI hospitalisations.
y_all<-melt(t(SARI_y))
colnames(y_all)<-c("t","s","y")
y_all$s<-federal_units[y_all$s]

# Set up data for plots:
xy_all<-inner_join(y_all,x_all,by=c('t','s') )%>%mutate(not_cov=y-x)%>%filter(s%in%federal_units[1:S])
N_max<-dim(SARI_y)[2]

# Define colour palatte for plots:
colour_palate<-viridis(n=length(Nowcast_list), end=0.9)

# Define subsection of federal units for plots:
fu_abbr=c("AM","DF","RJ","RR","SP","TO")
regions<-which(federal_units%in%c("AM","DF","RJ","RR","SP","TO"))
fu_full=c("Amazonas","Distrito Federal","Rio de Janeiro","Roraima","São Paulo","Tocantins")

# Calculate censored COVID, SARI and COVID/SARI (mu) available at each nowcast date. 
censored_x_samples<-censored_x_melt<-list()
censored_mu_samples<-list()
censored_z<-list()
censored_y<-censored_y_melt<-list()
for(j in 1:length(Nowcast_list)){
  N_now<-Nowcast_list[j]
  N<-W
  censored_x_samples[[j]]<-SARI_x[(1+N_now-W):(N_now),1:S]
  censored_x_samples[[j]][(N-R+1):N,]<-floor(censored_x_samples[[j]][(N-R+1):N,]*seq(from=0.95, to=0.25, length=R))
  censored_z[[j]]<-SARI_z[,(1+N_now-W):(N_now),]
  
    for(s in 1:S){
    censored_z[[j]][s,,][outer(1:dim(SARI_z[s,(1+N_now-W):(N_now),])[1], 0:(dim(SARI_z[s,(1+N_now-W):(N_now),])[2]-1), FUN = "+") > N] <- NA
  }
  
  # Censor the totals.
  censored_y[[j]] <- apply(censored_z[[j]][,,1:D_max],c(1,2),sum,na.rm=TRUE)
  
  censored_mu_samples[[j]]<-censored_x_samples[[j]]/t(censored_y[[j]])
  censored_mu_samples[[j]]<-melt( censored_mu_samples[[j]], value.name = "cen_mu",variable.name = "s" )%>%mutate(t=rep(1:N,S),Nowcast_Date=paste(Nowcast_dates[j]))
  censored_x_melt[[j]]<-melt( censored_x_samples[[j]], value.name = "cen_x",variable.name = "s" )%>%mutate(t=rep(1:N,S),Nowcast_Date=paste(Nowcast_dates[j]))
  censored_y_melt[[j]]<-melt( t(censored_y[[j]]), value.name = "cen_y",variable.name = "s" )%>%mutate(t=rep(1:N,S),Nowcast_Date=paste(Nowcast_dates[j]),s=federal_units[Var2])
}


#### COMPARE MODEL OUTPUTS ####
# WARNING: ALL MODEL SCRIPTS MUST BE RUN FIRST AND OUTPUT SAVED.

# Load CLR model results.
load("~/media/alba/Disk 1/OneDrive/PhD/comp 2023 paper code/FInal_code/OMEGA CODE/cluster_clr_x.RData")
# Load CLR age model results.
load("~/media/alba/Disk 1/OneDrive/PhD/comp 2023 paper code/FInal_code/OMEGA CODE/cluster_clr_age_x.RData")
# Load Hazard results.
load("~/media/alba/Disk 1/OneDrive/PhD/comp 2023 paper code/FInal_code/OMEGA CODE/cluster_hazard_x.RData")
# Load Survivor results.
load("~/media/alba/Disk 1/OneDrive/PhD/comp 2023 paper code/FInal_code/OMEGA CODE/cluster_survivor_x.RData")
# Load Survivor age results.
load("~/media/alba/Disk 1/OneDrive/PhD/comp 2023 paper code/FInal_code/OMEGA CODE/cluster_survivor_age_x.RData")
# Load CLR no tensor model results.
load("~/media/alba/Disk 1/OneDrive/PhD/comp 2023 paper code/FInal_code/OMEGA CODE/cluster_notensor_x.RData")
# Load CLR no link model results.
load("~/media/alba/Disk 1/OneDrive/PhD/comp 2023 paper code/FInal_code/OMEGA CODE/cluster_clr_nolink_x.RData")

# Create lists so store MCMC samples.
SARI_clr<-list()
SARI_combined_clr<-list()
SARI_age<-list()
SARI_combined_age<-list()
SARI_hazard<-list()
SARI_combined_hazard<-list()
SARI_Survivor<-list()
SARI_combined_Survivor<-list()
SARI_nolink<-list()
SARI_combined_clr_nolink<-list()
SARI_Survivor_age<-list()
SARI_combined_Survivor_age<-list()
SARI_nointer<-list()
SARI_combined_clr_nointer<-list()

for(j in length(Nowcast_list):1){
  # CLR Model With Link
  SARI_clr[[j]]<-as.mcmc.list(SARI_clr_output2[[j]])
  SARI_combined_clr[[j]] <- as_tibble(do.call('rbind',SARI_clr[[j]]))
  # CLR age model
  SARI_age[[j]]<-as.mcmc.list(SARI_age_output[[j]])
  SARI_combined_age[[j]] <- as_tibble(do.call('rbind',SARI_age[[j]]))
  # Hazard model
  SARI_hazard[[j]]<-as.mcmc.list(SARI_clr_haz[[j]])
  SARI_combined_hazard[[j]] <- as_tibble(do.call('rbind',SARI_hazard[[j]]))
  # Survivor model
  SARI_Survivor[[j]]<-as.mcmc.list(SARI_surv[[j]])
  SARI_combined_Survivor[[j]] <- as_tibble(do.call('rbind',SARI_Survivor[[j]]))
  # Survivor age model
  SARI_Survivor_age[[j]]<-as.mcmc.list(SARI_age_output_surv_x[[j]])
  SARI_combined_Survivor_age[[j]] <- as_tibble(do.call('rbind',SARI_Survivor_age[[j]]))
  # CLR model no link model.
  SARI_nolink[[j]]<-as.mcmc.list(SARI_clr_nolink[[j]])
  SARI_combined_clr_nolink[[j]] <- as_tibble(do.call('rbind',SARI_nolink[[j]]))
  # CLR no tensormodel 
  SARI_nointer[[j]]<-as.mcmc.list(SARI_clr_notensor[[j]])
  SARI_combined_clr_nointer[[j]] <- as_tibble(do.call('rbind',SARI_nointer[[j]]))
  
}


# Compare mu for different regions and models 
n_sim <- dim(SARI_combined_Survivor_age[[1]])[1]

# Compare SARI & COVID-positive predictions:
# Create lists to save model samples. 
sari_clr_quant<-covid_clr_quant<-sari_age_quant<-covid_age_quant<-sari_nolink_quant<-covid_nolink_quant<-sari_nointer_quant<-covid_nointer_quant<-
  sari_hazard_quant<-covid_hazard_quant<-sari_survivor_quant<-covid_survivor_quant<-sari_surv_age_quant<-covid_surv_age_quant<-list()
x_over_y_age<-covid_age_raw<-sari_age_raw<-list()
x_over_y_clr<-covid_clr_raw<-sari_clr_raw<-list()
x_over_y_nolink<-covid_nolink_raw<-sari_nolink_raw<-list()
x_over_y_notensor<-covid_notensor_raw<-sari_notensor_raw<-list()
x_over_y_Surv_age<-covid_Surv_age_raw<-sari_Surv_age_raw<-list()
x_over_y_Surv<-covid_Surv_raw<-sari_Surv_raw<-list()
x_over_y_haz<-covid_haz_raw<-sari_haz_raw<-list()

# Compute COVID, SARI and COVID/SARI qunantiles. 
for(j in length(Nowcast_list):1){
  # CLR Model 
  sari_clr_quant[[j]]<-select(SARI_combined_clr[[j]],starts_with('y'))%>%as.matrix()%>%array(dim=c(n_sim,S,W))%>%
    apply(c(2,3),quantile,c(0.025,0.5,0.975))%>%melt(varnames=c('quantile','s','t'),value.name='y')%>%mutate(t=t+Nowcast_list[j]-W)%>%
    spread(quantile,y)%>%mutate(s=federal_units[s])
  
  covid_clr_quant[[j]]<-select(SARI_combined_clr[[j]],starts_with('x'))%>%as.matrix()%>%array(dim=c(n_sim,W,S))%>%
    apply(c(2,3),quantile,c(0.025,0.5,0.975))%>%melt(varnames=c('quantile','t','s'),value.name='y')%>%mutate(t=t+Nowcast_list[j]-W)%>%
    spread(quantile,y)%>%mutate(s=federal_units[s])
 
   sari_clr_raw[[j]]<-select(SARI_combined_clr[[j]],starts_with('y'))%>%unlist()%>%as.numeric()%>%as.matrix()%>%array(dim=c(n_sim,S,W))%>%aperm(c(1,3,2))
  covid_clr_raw[[j]]<-select(SARI_combined_clr[[j]],starts_with('x'))%>%unlist()%>%as.numeric()%>%as.matrix()%>%array(dim=c(n_sim,W,S))
  x_over_y_clr[[j]]<-covid_clr_raw[[j]]/(sari_clr_raw[[j]]+0.0000000001)
  x_over_y_clr[[j]]<-x_over_y_clr[[j]]%>%
    apply(c(2,3),quantile,c(0.025,0.5,0.975))%>%melt(varnames=c('quantile','t','s'),value.name='y')%>%mutate(t=t+Nowcast_list[j]-W)%>%
    spread(quantile,y)%>%mutate(s=federal_units[s])
  
  # CLR No Link
  sari_nolink_quant[[j]]<-select(SARI_combined_clr_nolink[[j]],starts_with('y'))%>%as.matrix()%>%array(dim=c(n_sim,S,W))%>%
   apply(c(2,3),quantile,c(0.025,0.5,0.975))%>%melt(varnames=c('quantile','s','t'),value.name='y')%>%mutate(t=t+Nowcast_list[j]-W)%>%
    spread(quantile,y)%>%mutate(s=federal_units[s])
  
    covid_nolink_quant[[j]]<-select(SARI_combined_clr_nolink[[j]],starts_with('x'))%>%as.matrix()%>%array(dim=c(n_sim,W,S))%>%
    apply(c(2,3),quantile,c(0.025,0.5,0.975))%>%melt(varnames=c('quantile','t','s'),value.name='y')%>%mutate(t=t+Nowcast_list[j]-W)%>%
    spread(quantile,y)%>%mutate(s=federal_units[s])
    
    sari_nolink_raw[[j]]<-select(SARI_combined_clr_nolink[[j]],starts_with('y'))%>%unlist()%>%as.numeric()%>%as.matrix()%>%array(dim=c(n_sim,S,W))%>%aperm(c(1,3,2))
    covid_nolink_raw[[j]]<-select(SARI_combined_clr_nolink[[j]],starts_with('x'))%>%unlist()%>%as.numeric()%>%as.matrix()%>%array(dim=c(n_sim,W,S))
    x_over_y_nolink[[j]]<-   covid_nolink_raw[[j]]/(sari_nolink_raw[[j]]+0.0000000001)
    x_over_y_nolink[[j]]<-x_over_y_nolink[[j]]%>%
      apply(c(2,3),quantile,c(0.025,0.5,0.975))%>%melt(varnames=c('quantile','t','s'),value.name='y')%>%mutate(t=t+Nowcast_list[j]-W)%>%
      spread(quantile,y)%>%mutate(s=federal_units[s])
    
  
  # CLR No Tensor
  sari_nointer_quant[[j]]<-select(SARI_combined_clr_nointer[[j]],starts_with('y'))%>%as.matrix()%>%array(dim=c(n_sim,S,W))%>%
   apply(c(2,3),quantile,c(0.025,0.5,0.975))%>%melt(varnames=c('quantile','s','t'),value.name='y')%>%mutate(t=t+Nowcast_list[j]-W)%>%
    spread(quantile,y)%>%mutate(s=federal_units[s])
  
  
    covid_nointer_quant[[j]]<-select(SARI_combined_clr_nointer[[j]],starts_with('x'))%>%as.matrix()%>%array(dim=c(n_sim,W,S))%>%
    apply(c(2,3),quantile,c(0.025,0.5,0.975))%>%melt(varnames=c('quantile','t','s'),value.name='y')%>%mutate(t=t+Nowcast_list[j]-W)%>%
    spread(quantile,y)%>%mutate(s=federal_units[s])
   
     sari_notensor_raw[[j]]<-select(SARI_combined_clr_nointer[[j]],starts_with('y'))%>%unlist()%>%as.numeric()%>%as.matrix()%>%array(dim=c(n_sim,S,W))%>%aperm(c(1,3,2))
    covid_notensor_raw[[j]]<-select(SARI_combined_clr_nointer[[j]],starts_with('x'))%>%unlist()%>%as.numeric()%>%as.matrix()%>%array(dim=c(n_sim,W,S))
    x_over_y_notensor[[j]]<-   covid_notensor_raw[[j]]/(sari_notensor_raw[[j]]+0.0000000001)
    x_over_y_notensor[[j]]<-x_over_y_notensor[[j]]%>%
      apply(c(2,3),quantile,c(0.025,0.5,0.975))%>%melt(varnames=c('quantile','t','s'),value.name='y')%>%mutate(t=t+Nowcast_list[j]-W)%>%
      spread(quantile,y)%>%mutate(s=federal_units[s])
    
  
  # CLR Age Model 
    sari_age_quant[[j]]<-select(SARI_combined_age[[j]],starts_with('y'))%>%as.matrix()%>%array(dim=c(n_sim,S,W))%>%
    apply(c(2,3),quantile,c(0.025,0.5,0.975))%>%melt(varnames=c('quantile','s','t'),value.name='y')%>%mutate(t=t+Nowcast_list[j]-W)%>%
    spread(quantile,y)%>%mutate(s=federal_units[s])

    covid_age_quant[[j]]<-select(SARI_combined_age[[j]],starts_with('x'))%>%as.matrix()%>%array(dim=c(n_sim,W,S))%>%
      apply(c(2,3),quantile,c(0.025,0.5,0.975))%>%melt(varnames=c('quantile','t','s'),value.name='y')%>%mutate(t=t+Nowcast_list[j]-W)%>%
      spread(quantile,y)%>%mutate(s=federal_units[s])
    
    sari_age_raw[[j]]<-select(SARI_combined_age[[j]],starts_with('y'))%>%unlist()%>%as.numeric()%>%as.matrix()%>%array(dim=c(n_sim,S,W))%>%aperm(c(1,3,2))
    covid_age_raw[[j]]<-select(SARI_combined_age[[j]],starts_with('x'))%>%unlist()%>%as.numeric()%>%as.matrix()%>%array(dim=c(n_sim,W,S))
    x_over_y_age[[j]]<-   covid_age_raw[[j]]/(sari_age_raw[[j]]+0.0000000001)
    x_over_y_age[[j]]<-x_over_y_age[[j]]%>%
      apply(c(2,3),quantile,c(0.025,0.5,0.975))%>%melt(varnames=c('quantile','t','s'),value.name='y')%>%mutate(t=t+Nowcast_list[j]-W)%>%
      spread(quantile,y)%>%mutate(s=federal_units[s])
    
  # Survivor Age Model 
    sari_surv_age_quant[[j]]<-select(SARI_combined_Survivor_age[[j]],starts_with('y'))%>%as.matrix()%>%array(dim=c(n_sim,S,W))%>%
      apply(c(2,3),quantile,c(0.025,0.5,0.975))%>%melt(varnames=c('quantile','s','t'),value.name='y')%>%mutate(t=t+Nowcast_list[j]-W)%>%
      spread(quantile,y)%>%mutate(s=federal_units[s])
    
    
    covid_surv_age_quant[[j]]<-select(SARI_combined_Survivor_age[[j]],starts_with('x'))%>%as.matrix()%>%array(dim=c(n_sim,W,S))%>%
      apply(c(2,3),quantile,c(0.025,0.5,0.975))%>%melt(varnames=c('quantile','t','s'),value.name='y')%>%mutate(t=t+Nowcast_list[j]-W)%>%
      spread(quantile,y)%>%mutate(s=federal_units[s])
    
    
    sari_Surv_age_raw[[j]]<-select(SARI_combined_Survivor_age[[j]],starts_with('y'))%>%unlist()%>%as.numeric()%>%as.matrix()%>%array(dim=c(n_sim,S,W))%>%aperm(c(1,3,2))
    covid_Surv_age_raw[[j]]<-select(SARI_combined_Survivor_age[[j]],starts_with('x'))%>%unlist()%>%as.numeric()%>%as.matrix()%>%array(dim=c(n_sim,W,S))
    x_over_y_Surv_age[[j]]<-   covid_Surv_age_raw[[j]]/(sari_Surv_age_raw[[j]]+0.0000000001)
    x_over_y_Surv_age[[j]]<-x_over_y_Surv_age[[j]]%>%
      apply(c(2,3),quantile,c(0.025,0.5,0.975))%>%melt(varnames=c('quantile','t','s'),value.name='y')%>%mutate(t=t+Nowcast_list[j]-W)%>%
      spread(quantile,y)%>%mutate(s=federal_units[s])
    
    # Hazard Model 
    sari_hazard_quant[[j]]<-select(SARI_combined_hazard[[j]],starts_with('y'))%>%as.matrix()%>%array(dim=c(n_sim,S,W))%>%
      apply(c(2,3),quantile,c(0.025,0.5,0.975))%>%melt(varnames=c('quantile','s','t'),value.name='y')%>%mutate(t=t+Nowcast_list[j]-W)%>%
      spread(quantile,y)%>%mutate(s=federal_units[s])
    
    covid_hazard_quant[[j]]<-select(SARI_combined_hazard[[j]],starts_with('x'))%>%as.matrix()%>%array(dim=c(n_sim,W,S))%>%
      apply(c(2,3),quantile,c(0.025,0.5,0.975))%>%melt(varnames=c('quantile','t','s'),value.name='y')%>%mutate(t=t+Nowcast_list[j]-W)%>%
      spread(quantile,y)%>%mutate(s=federal_units[s])
    
    sari_haz_raw[[j]]<-select(SARI_combined_hazard[[j]],starts_with('y'))%>%unlist()%>%as.numeric()%>%as.matrix()%>%array(dim=c(n_sim,S,W))%>%aperm(c(1,3,2))
    covid_haz_raw[[j]]<-select(SARI_combined_hazard[[j]],starts_with('x'))%>%unlist()%>%as.numeric()%>%as.matrix()%>%array(dim=c(n_sim,W,S))
    x_over_y_haz[[j]]<-   covid_haz_raw[[j]]/(sari_haz_raw[[j]]+0.0000000001)
    x_over_y_haz[[j]]<-x_over_y_haz[[j]]%>%
      apply(c(2,3),quantile,c(0.025,0.5,0.975))%>%melt(varnames=c('quantile','t','s'),value.name='y')%>%mutate(t=t+Nowcast_list[j]-W)%>%
      spread(quantile,y)%>%mutate(s=federal_units[s])
    
     # Survivor Model 
    sari_survivor_quant[[j]]<-select(SARI_combined_Survivor[[j]],starts_with('y'))%>%as.matrix()%>%array(dim=c(n_sim,S,W))%>%
      apply(c(2,3),quantile,c(0.025,0.5,0.975))%>%melt(varnames=c('quantile','s','t'),value.name='y')%>%mutate(t=t+Nowcast_list[j]-W)%>%
      spread(quantile,y)%>%mutate(s=federal_units[s])
    
    covid_survivor_quant[[j]]<-select(SARI_combined_Survivor[[j]],starts_with('x'))%>%as.matrix()%>%array(dim=c(n_sim,W,S))%>%
      apply(c(2,3),quantile,c(0.025,0.5,0.975))%>%melt(varnames=c('quantile','t','s'),value.name='y')%>%mutate(t=t+Nowcast_list[j]-W)%>%
      spread(quantile,y)%>%mutate(s=federal_units[s])
    
    sari_Surv_raw[[j]]<-select(SARI_combined_Survivor[[j]],starts_with('y'))%>%unlist()%>%as.numeric()%>%as.matrix()%>%array(dim=c(n_sim,S,W))%>%aperm(c(1,3,2))
    covid_Surv_raw[[j]]<-select(SARI_combined_Survivor[[j]],starts_with('x'))%>%unlist()%>%as.numeric()%>%as.matrix()%>%array(dim=c(n_sim,W,S))
    x_over_y_Surv[[j]]<-   covid_Surv_raw[[j]]/(sari_Surv_raw[[j]]+0.0000000001)
    x_over_y_Surv[[j]]<-x_over_y_Surv[[j]]%>%
      apply(c(2,3),quantile,c(0.025,0.5,0.975))%>%melt(varnames=c('quantile','t','s'),value.name='y')%>%mutate(t=t+Nowcast_list[j]-W)%>%
      spread(quantile,y)%>%mutate(s=federal_units[s])
    
    }

# Save quantiles of COVID and SARI samples. 
save(covid_survivor_quant,sari_survivor_quant,covid_surv_age_quant,sari_surv_age_quant,covid_age_quant,sari_age_quant,
     covid_clr_quant,sari_clr_quant, covid_nolink_quant,sari_nolink_quant, covid_nointer_quant,sari_nointer_quant, covid_hazard_quant,sari_hazard_quant,
     file="x_quantiles.RData")

# Plot proportion of SARI cases that are COVID-positive (x/y).
x_over_y_comp<-rbind( filter(xy_all,s%in%c("SP"),t<=Nowcast_list[length(Nowcast_list)],t>25)%>%
                        mutate(s=factor(s, levels=fu_abbr))%>%mutate(fu=fu_full[s], m=c("CLR Model No Link"),x_over_y=x/y),
                      filter(xy_all,s%in%c("SP"),t<=Nowcast_list[length(Nowcast_list)],t>25)%>%
                        mutate(s=factor(s, levels=fu_abbr))%>%mutate(fu=fu_full[s], m=c("CLR Model"),x_over_y=x/y))
x_over_y_comp$m<-factor(x_over_y_comp$m,levels=unique(x_over_y_comp$m))

compare_x_over_y<-ggplot()+
  geom_point(data=x_over_y_comp,aes(x=as.Date(first_date)+(t-1)*7,y=x_over_y))+ 
  facet_wrap(~m,nrow=1)+
  scale_color_viridis(option="D",discrete = TRUE)+
  scale_fill_viridis(option="D",discrete = TRUE)+
  theme_minimal()+labs(y=NULL, x=NULL, title="Comparison of regional COVID-positive SARI proportions", 
                       colour="Nowcast date", fill="Nowcast date")+
  theme(legend.position="bottom",axis.title.x = NULL)#+geom_vline(xintercept = 51)

quantiles_x_over_y<-list()
for(j in c(15,11,7,3)){
  quantiles_x_over_y[[j]]<-rbind( filter(x_over_y_nolink[[j]],s%in%c("SP"),t>25)%>%mutate(s=factor(s, levels=fu_abbr))%>%mutate(fu=fu_full[s])%>%
                                    mutate(Nowcast_Date=paste(Nowcast_dates[j]),m=c("CLR No Link Model")),
                                  filter(x_over_y_clr[[j]],s%in%c("SP"),t>25)%>%mutate(s=factor(s, levels=fu_abbr))%>%mutate(fu=fu_full[s])%>%
                                    mutate(Nowcast_Date=paste(Nowcast_dates[j]),m=c("CLR Model")))
  quantiles_x_over_y[[j]]$m<-factor(quantiles_x_over_y[[j]]$m,levels=unique(quantiles_x_over_y[[j]]$m))
  compare_x_over_y<-compare_x_over_y+
    geom_line(data=quantiles_x_over_y[[j]], aes(x=as.Date(first_date)+(t-1)*7, y=`50%`, colour=Nowcast_Date))+
    geom_line(data=filter(censored_mu_samples[[j]],s%in%c("SP"),t<=Nowcast_list[length(Nowcast_list)],t>25)%>%mutate(t=t+Nowcast_list[j]-W), aes(x=as.Date(first_date)+(t-1)*7, y=cen_mu, colour=Nowcast_Date),linetype="dashed")+
    geom_ribbon(data=quantiles_x_over_y[[j]], aes(x=as.Date(first_date)+(t-1)*7, ymin=`2.5%`, ymax=`97.5%`,fill=Nowcast_Date), alpha=0.2)
  
  
}
compare_x_over_y
# Save plot.
ggsave(compare_x_over_y, file='Plots/compare_x_over_y_clr.pdf',width=9,height=3)


# Plot nowcasts for a single nowcast date in São Paulo:
# COVID-positive nowcast predictions.
covid_comp<- filter(xy_all,s%in%fu_abbr,t<=Nowcast_list[length(Nowcast_list)],t>40)%>%
  mutate(s=factor(s, levels=fu_abbr))%>%mutate(fu=fu_full[s],hos="COVID-positive", m=c("Survivor Age Model"),y=x)%>%
  filter(s=='SP')
# SARI nowcast predictions.
sari_comp<- filter(xy_all,s%in%fu_abbr,t<=Nowcast_list[length(Nowcast_list)],t>40)%>%
  mutate(s=factor(s, levels=fu_abbr))%>%mutate(fu=fu_full[s],hos="All SARI", m=c("Survivor Age Model"),y=y)%>%filter(s=='SP')
# All nowcast predictions. 
all_comp<-rbind(covid_comp,sari_comp)

compare_SP<-ggplot()+
  geom_point(data=all_comp,aes(x=as.Date(first_date)+(t-1)*7,y=(y)))+ 
  facet_grid2(~hos, scales = 'free',independent = "all")+
  scale_color_viridis(option="D",discrete = TRUE)+
  scale_fill_viridis(option="D",discrete = TRUE)+
  theme_minimal()+labs(y=NULL, x=NULL, title="Nowcasts of SARI hospitalisations in São Paulo", 
                       colour="Nowcast Date", fill="Nowcast Date",subtitle = paste(Nowcast_dates[13]))+
  theme(legend.position="none")

quantiles_covid<-quantiles_all<-quantiles_sari<-list()
for(j in c(12)){
  quantiles_covid[[j]]<-filter(covid_surv_age_quant[[j]],s%in%fu_abbr, t>40)%>%mutate(s=factor(s, levels=fu_abbr))%>%mutate(fu=fu_full[s])%>%
    mutate(Nowcast_Date=paste(Nowcast_dates[j]),m=c("Survivor Age Model"),hos="COVID-positive")%>%
    filter(s=='SP')
  quantiles_sari[[j]]<-filter(sari_surv_age_quant[[j]],s%in%fu_abbr,t>40)%>%mutate(s=factor(s, levels=fu_abbr))%>%mutate(fu=fu_full[s])%>%
    mutate(Nowcast_Date=paste(Nowcast_dates[j]),m=c("Survivor Model"), hos="All SARI")%>%filter(s=='SP')
  quantiles_all[[j]]<-rbind(  quantiles_covid[[j]],quantiles_sari[[j]])
  
  # Plot of COVID and SARI nowcasts for São Paulo
  compare_SP<-compare_SP+
    geom_line(data=filter(censored_x_melt[[j]],s%in%c("SP"),t<=Nowcast_list[length(Nowcast_list)],t>30)%>%mutate(t=t+Nowcast_list[j]-W,hos="COVID-positive"), aes(x=as.Date(first_date)+(t-1)*7, y=cen_x, colour=Nowcast_Date),linetype="dashed")+
    geom_line(data=filter(censored_y_melt[[j]],s%in%c("SP"),t<=Nowcast_list[length(Nowcast_list)],t>30)%>%mutate(t=t+Nowcast_list[j]-W,hos="All SARI"), aes(x=as.Date(first_date)+(t-1)*7, y=cen_y, colour=Nowcast_Date),linetype="dashed")+
    geom_line(data=quantiles_all[[j]], aes(x=as.Date(first_date)+(t-1)*7, y=(`50%`), colour=Nowcast_Date))+
    geom_ribbon(data=quantiles_all[[j]], aes(x=as.Date(first_date)+(t-1)*7, ymin=(`2.5%`), ymax=(`97.5%`),fill=Nowcast_Date), alpha=0.15)
  
}

compare_SP
# Save plot.
ggsave(compare_SP, file='Plots/compare_sari_covid.pdf',width=9,height=3.5)


# Compare Nowcast performance of models for main article:

# SARI nowcasts.
nowcasts_compare<-list()
SARI_y_totals<-apply(SARI_y,2,sum)
SARI_z_totals<-apply(SARI_z,c(2,3),sum)
SARI_x_totals<-apply(SARI_x,1,sum)
for(j in 1:length(Nowcast_list)){
  nowcasts_compare[[j]] <- rbind(full_join(sari_clr_quant[[j]],filter(xy_all[,1:3], t%in%c((Nowcast_list[j]-W+1):Nowcast_list[j])),by=c('s','t'))%>%mutate(m=c('CLR'),d=t-Nowcast_list[j],now=j),
                                 full_join(sari_age_quant[[j]],filter(xy_all[,1:3], t%in%c((Nowcast_list[j]-W+1):Nowcast_list[j])),by=c('s','t'))%>%mutate(m=c('CLR Age'),d=t-Nowcast_list[j],now=j),
                                 full_join(sari_survivor_quant[[j]],filter(xy_all[,1:3], t%in%c((Nowcast_list[j]-W+1):Nowcast_list[j])),by=c('s','t'))%>%mutate(m=c('Survivor'),d=t-Nowcast_list[j],now=j),
                                 full_join(sari_surv_age_quant[[j]],filter(xy_all[,1:3], t%in%c((Nowcast_list[j]-W+1):Nowcast_list[j])),by=c('s','t'))%>%mutate(m=c('Survivor Age'),d=t-Nowcast_list[j],now=j)
                                 )%>%filter(d>-D_max)
}
nowcast_dmax<-array(unlist(nowcasts_compare), dim=c(D_max*4*S,ncol(nowcasts_compare[[1]]),length(Nowcast_list)))%>%
  apply(c(2),c)%>%as_tibble()
colnames(nowcast_dmax)<-colnames(nowcasts_compare[[1]])                
nowcast_dmax<-nowcast_dmax%>%mutate(nowcast_dmax,t=as.numeric(t),`2.5%`=as.numeric(`2.5%`),`50%`=as.numeric(`50%`),`97.5%`=as.numeric(`97.5%`),
                                    m=factor(m,levels=c('CLR','CLR Age','Survivor','Survivor Age' )),y=as.numeric(y),
                                    d=as.numeric(d),now=as.numeric(now))
nowcast_dmax_summary<-nowcast_dmax%>%group_by(m,d,s)%>%
  #group_by(m,d)%>%
  summarise(Coverage=mean(y>=`2.5%`&y<=`97.5%`),
            `Prediction interval width`=as.numeric(mean(`97.5%`-`2.5%`)),
            `Mean absolute error`=as.numeric(mean(abs(`50%`-y))))%>%
  pivot_longer(cols=c('Coverage','Mean absolute error','Prediction interval width'),names_to = 'metric')%>%
  mutate(m=factor(m,levels=c('CLR','CLR Age',  'Survivor','Survivor Age')),
         metric=factor(metric,levels=c('Mean absolute error','Prediction interval width','Coverage')))

nowcast_dmax_summary_average<-nowcast_dmax_summary%>%group_by(m,d,metric)%>%
  summarise(value.mean=mean(value), value.median=median(value))
nowcast_dmax_summary_average$value.mean[nowcast_dmax_summary_average$metric%in%c("Prediction Interval Width","Mean Absolute Error")]<-NA
nowcast_dmax_summary_average$value.median[nowcast_dmax_summary_average$metric%in%c("Coverage")]<-NA
nowcast_dmax_summary_average<-nowcast_dmax_summary_average%>%group_by(m,d,metric)%>%mutate(value=sum(value.mean,value.median,na.rm=TRUE),hos="All SARI")


# COVID-19 nowcasts.
covid_nowcasts_compare<-list()
for(j in 1:length(Nowcast_list)){
  covid_nowcasts_compare[[j]] <- rbind(
    full_join(covid_clr_quant[[j]],filter(xy_all[,c(1,2,4)], t%in%c((Nowcast_list[j]-W+1):Nowcast_list[j])),by=c('s','t'))%>%mutate(m=c('CLR'),d=t-Nowcast_list[j],now=j),
    full_join(covid_age_quant[[j]],filter(xy_all[,c(1,2,4)], t%in%c((Nowcast_list[j]-W+1):Nowcast_list[j])),by=c('s','t'))%>%mutate(m=c('CLR Age'),d=t-Nowcast_list[j],now=j),
    full_join(covid_survivor_quant[[j]],filter(xy_all[,c(1,2,4)], t%in%c((Nowcast_list[j]-W+1):Nowcast_list[j])),by=c('s','t'))%>%mutate(m=c('Survivor'),d=t-Nowcast_list[j],now=j),
    full_join(covid_surv_age_quant[[j]],filter(xy_all[,c(1,2,4)], t%in%c((Nowcast_list[j]-W+1):Nowcast_list[j])),by=c('s','t'))%>%mutate(m=c('Survivor Age'),d=t-Nowcast_list[j],now=j)
  )%>%filter(d>-R)
}
covid_nowcast_dmax<-array(unlist(covid_nowcasts_compare), dim=c(S*R*4,ncol(covid_nowcasts_compare[[1]]),length(Nowcast_list)))%>%
  apply(c(2),c)%>%as_tibble()
colnames(covid_nowcast_dmax)<-colnames(covid_nowcasts_compare[[1]])   
covid_nowcast_dmax<-covid_nowcast_dmax%>%mutate(t=as.numeric(t),`2.5%`=as.numeric(`2.5%`),`50%`=as.numeric(`50%`),`97.5%`=as.numeric(`97.5%`),
                                                m=factor(m,levels=c('CLR','CLR Age', 'Survivor','Survivor Age'
                                                )),x=as.numeric(x),
                                                d=as.numeric(d),now=as.numeric(now))
covid_nowcast_dmax_summary<-covid_nowcast_dmax%>%group_by(m,d,s)%>%
  summarise(Coverage=mean(x>=`2.5%`&x<=`97.5%`),
            `Prediction interval width`=mean(`97.5%`-`2.5%`),
            `Mean absolute error`=mean(abs(`50%`-x)))%>%
  pivot_longer(cols=c('Coverage','Mean absolute error','Prediction interval width'),names_to = 'metric')%>%
  mutate(m=factor(m,levels=c('CLR','CLR Age', 'Survivor','Survivor Age'
  )),
  metric=factor(metric,levels=c('Mean absolute error','Prediction interval width','Coverage')))


#plot_data_covid<-filter(covid_nowcast_dmax_summary, d==0)%>%pivot_wider(names_from=c('m'))%>%filter(metric=="Mean Absolute Error") 
#ggplot(plot_data_covid)+   geom_text(aes(x=`CLR`,y=`CLR Age`/`CLR`,label=s))+geom_hline(yintercept = 1)

covid_nowcast_dmax_summary_average<-covid_nowcast_dmax_summary%>%group_by(m,d,metric)%>%
  summarise(value.mean=mean(value), value.median=median(value))
covid_nowcast_dmax_summary_average$value.mean[covid_nowcast_dmax_summary_average$metric%in%c("Prediction Interval Width","Mean Absolute Error")]<-NA
covid_nowcast_dmax_summary_average$value.median[covid_nowcast_dmax_summary_average$metric%in%c("Coverage")]<-NA
covid_nowcast_dmax_summary_average<-covid_nowcast_dmax_summary_average%>%group_by(m,d,metric)%>%mutate(value=sum(value.mean,value.median,na.rm=TRUE), hos="COVID-positive SARI")

nowcast_dmax_summary<-rbind(nowcast_dmax_summary_average,covid_nowcast_dmax_summary_average)



# Plot of prediction performance metrics for main paper.
totals_compare_age <- ggplot(nowcast_dmax_summary)+
  geom_line(aes(x=d,y=value,colour=m),alpha=0.6)+
  geom_point(aes(x=d,y=value,colour=m,shape=m),alpha=0.6)+
  facet_grid2(hos~metric, scales = c("free"),independent = c("all"))+
  scale_colour_manual(name=NULL,values = c( "#0072B2","#009E73","#D55E00","#E69F00"))+
  scale_shape_manual(name=NULL,values=c(15,8,17,3))+
  scale_x_continuous()+
  labs(x='Prediction time difference (weeks)',y=NULL,
       title='Performance metrics from the rolling prediction experiment')+
  theme_minimal()+
  theme(legend.position = "bottom")
totals_compare_age
# Save plot. 
ggsave(totals_compare_age,file='Plots/totals_compare_paper.pdf',width=9,height=5)



# Compare nowcast performance of models for appendix:

# SARI nowcasts.
nowcasts_compare<-list()
SARI_y_totals<-apply(SARI_y,2,sum)
SARI_z_totals<-apply(SARI_z,c(2,3),sum)
SARI_x_totals<-apply(SARI_x,1,sum)
for(j in 1:length(Nowcast_list)){
  nowcasts_compare[[j]] <- rbind(full_join(sari_clr_quant[[j]],filter(xy_all[,1:3], t%in%c((Nowcast_list[j]-W+1):Nowcast_list[j])),by=c('s','t'))%>%mutate(m=c('CLR'),d=t-Nowcast_list[j],now=j),
                                 full_join(sari_nolink_quant[[j]],filter(xy_all[,1:3], t%in%c((Nowcast_list[j]-W+1):Nowcast_list[j])),by=c('s','t'))%>%mutate(m=c('CLR No Link'),d=t-Nowcast_list[j],now=j),
                                  full_join(sari_nointer_quant[[j]],filter(xy_all[,1:3], t%in%c((Nowcast_list[j]-W+1):Nowcast_list[j])),by=c('s','t'))%>%mutate(m=c('CLR No Tensor'),d=t-Nowcast_list[j],now=j),
                                 full_join(sari_hazard_quant[[j]],filter(xy_all[,1:3], t%in%c((Nowcast_list[j]-W+1):Nowcast_list[j])),by=c('s','t'))%>%mutate(m=c('Hazard'),d=t-Nowcast_list[j],now=j),
                                 full_join(sari_survivor_quant[[j]],filter(xy_all[,1:3], t%in%c((Nowcast_list[j]-W+1):Nowcast_list[j])),by=c('s','t'))%>%mutate(m=c('Survivor'),d=t-Nowcast_list[j],now=j)
  )%>%filter(d>-D_max)
}
nowcast_dmax<-array(unlist(nowcasts_compare), dim=c(D_max*5*S,ncol(nowcasts_compare[[1]]),length(Nowcast_list)))%>%
  apply(c(2),c)%>%as_tibble()
colnames(nowcast_dmax)<-colnames(nowcasts_compare[[1]])                
nowcast_dmax<-nowcast_dmax%>%mutate(nowcast_dmax,t=as.numeric(t),`2.5%`=as.numeric(`2.5%`),`50%`=as.numeric(`50%`),`97.5%`=as.numeric(`97.5%`),
                                    m=factor(m,levels=c('CLR','CLR No Link','CLR No Tensor','Hazard','Survivor')),
                                    y=as.numeric(y),
                                    d=as.numeric(d),now=as.numeric(now))
nowcast_dmax_summary<-nowcast_dmax%>%group_by(m,d,s)%>%
  #group_by(m,d)%>%
  summarise(Coverage=mean(y>=`2.5%`&y<=`97.5%`),
            `Prediction interval width`=as.numeric(mean(`97.5%`-`2.5%`)),
            `Mean absolute error`=as.numeric(mean(abs(`50%`-y))))%>%
  pivot_longer(cols=c('Coverage','Mean absolute error','Prediction interval width'),names_to = 'metric')%>%
  mutate(m=factor(m,levels=c('CLR', 'CLR No Link','CLR No Tensor','Hazard','Survivor')),
 metric=factor(metric,levels=c('Mean absolute error','Prediction interval width','Coverage')))

nowcast_dmax_summary_average<-nowcast_dmax_summary%>%group_by(m,d,metric)%>%
  summarise(value.mean=mean(value), value.median=median(value))
nowcast_dmax_summary_average$value.mean[nowcast_dmax_summary_average$metric%in%c("Prediction Interval Width","Mean Absolute Error")]<-NA
nowcast_dmax_summary_average$value.median[nowcast_dmax_summary_average$metric%in%c("Coverage")]<-NA
nowcast_dmax_summary_average<-nowcast_dmax_summary_average%>%group_by(m,d,metric)%>%mutate(value=sum(value.mean,value.median,na.rm=TRUE),hos="All SARI")


# COVID-19 nowcasts.
covid_nowcasts_compare<-list()
for(j in 1:length(Nowcast_list)){
  covid_nowcasts_compare[[j]] <- rbind(
    full_join(covid_clr_quant[[j]],filter(xy_all[,c(1,2,4)], t%in%c((Nowcast_list[j]-W+1):Nowcast_list[j])),by=c('s','t'))%>%mutate(m=c('CLR'),d=t-Nowcast_list[j],now=j),
    full_join(covid_nolink_quant[[j]],filter(xy_all[,c(1,2,4)], t%in%c((Nowcast_list[j]-W+1):Nowcast_list[j])),by=c('s','t'))%>%mutate(m=c('CLR No Link'),d=t-Nowcast_list[j],now=j),
     full_join(covid_nointer_quant[[j]],filter(xy_all[,c(1,2,4)], t%in%c((Nowcast_list[j]-W+1):Nowcast_list[j])),by=c('s','t'))%>%mutate(m=c('CLR No Tensor'),d=t-Nowcast_list[j],now=j),
     full_join(covid_hazard_quant[[j]],filter(xy_all[,c(1,2,4)], t%in%c((Nowcast_list[j]-W+1):Nowcast_list[j])),by=c('s','t'))%>%mutate(m=c('Hazard'),d=t-Nowcast_list[j],now=j),
    full_join(covid_survivor_quant[[j]],filter(xy_all[,c(1,2,4)], t%in%c((Nowcast_list[j]-W+1):Nowcast_list[j])),by=c('s','t'))%>%mutate(m=c('Survivor'),d=t-Nowcast_list[j],now=j)
  )%>%filter(d>-R)
}
covid_nowcast_dmax<-array(unlist(covid_nowcasts_compare), dim=c(S*R*5,ncol(covid_nowcasts_compare[[1]]),length(Nowcast_list)))%>%
  apply(c(2),c)%>%as_tibble()
colnames(covid_nowcast_dmax)<-colnames(covid_nowcasts_compare[[1]])   
covid_nowcast_dmax<-covid_nowcast_dmax%>%mutate(t=as.numeric(t),`2.5%`=as.numeric(`2.5%`),`50%`=as.numeric(`50%`),`97.5%`=as.numeric(`97.5%`),
                                                m=factor(m,levels=c('CLR', 'CLR No Link','CLR No Tensor','Hazard',
                                                                    'Survivor')),x=as.numeric(x),
                                                d=as.numeric(d),now=as.numeric(now))
covid_nowcast_dmax_summary<-covid_nowcast_dmax%>%group_by(m,d,s)%>%
  summarise(Coverage=mean(x>=`2.5%`&x<=`97.5%`),
            `Prediction interval width`=mean(`97.5%`-`2.5%`),
            `Mean absolute error`=mean(abs(`50%`-x)))%>%
  pivot_longer(cols=c('Coverage','Mean absolute error','Prediction interval width'),names_to = 'metric')%>%
  mutate(m=factor(m,levels=c('CLR','CLR No Link','CLR No Tensor','Hazard', 'Survivor')),
  metric=factor(metric,levels=c('Mean absolute error','Prediction interval width','Coverage')))



covid_nowcast_dmax_summary_average<-covid_nowcast_dmax_summary%>%group_by(m,d,metric)%>%
  summarise(value.mean=mean(value), value.median=median(value))
covid_nowcast_dmax_summary_average$value.mean[covid_nowcast_dmax_summary_average$metric%in%c("Prediction Interval Width","Mean Absolute Error")]<-NA
covid_nowcast_dmax_summary_average$value.median[covid_nowcast_dmax_summary_average$metric%in%c("Coverage")]<-NA
covid_nowcast_dmax_summary_average<-covid_nowcast_dmax_summary_average%>%group_by(m,d,metric)%>%mutate(value=sum(value.mean,value.median,na.rm=TRUE), hos="COVID-positive SARI")

nowcast_dmax_summary<-rbind(nowcast_dmax_summary_average,covid_nowcast_dmax_summary_average)



# Plot of prediction performance metrics for appendix.
totals_compare_appendix <- ggplot(nowcast_dmax_summary)+
  geom_line(aes(x=d,y=value,colour=m),alpha=0.6)+
  geom_point(aes(x=d,y=value,colour=m,shape=m),alpha=0.6)+
  facet_grid2(hos~metric, scales = c("free"),independent = c("all"))+
  scale_colour_manual(name=NULL,values = c( "#0072B2","#56B4E9","#9137E1","#CC79A7","#D55E00"))+
  scale_shape_manual(name=NULL,values=c(15,4,6,16,17))+
  scale_x_continuous()+
  labs(x='Prediction time difference (weeks)',y=NULL,
       title='Performance metrics from the rolling prediction experiment')+
  theme_minimal()+
  theme(legend.position = "bottom")
totals_compare_appendix
# Save plot:
ggsave(totals_compare_appendix,file='Plots/totals_compare_appendix.pdf',width=9,height=5)


##### AGE EFFECTS PLOTS #####
# 
# # Hospitalization data grouped by age:
# 
# #Read in data..
# INFLUD21 <- read_delim("INFLUD21-30-01-2023.csv", delim = ";", escape_double = FALSE, trim_ws = TRUE)
# 
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
# 
# SARI_2021<- INFLUD21[,c(1:16, col_hospital, col_collection, col_results, col_pcr_result, col_notif_date, col_pcr, col_an_result, col_ant)]%>%
#   mutate(DT_NOTIFIC=as.Date((DT_NOTIFIC),format="%d/%m/%Y"),DT_DIGITA=as.Date((DT_DIGITA),format="%d/%m/%Y"),DT_SIN_PRI=as.Date((DT_SIN_PRI),format="%d/%m/%Y"))
# dates_2021<-sort(as.Date(unique(SARI_2021$DT_SIN_PRI),format="%d/%m/%Y"))
# 
# SARI_2022<- INFLUD22[,c(1:16, col_hospital, col_collection, col_results, col_pcr_result, col_notif_date, col_pcr, col_an_result, col_ant)]%>%
#   mutate(DT_NOTIFIC=as.Date((DT_NOTIFIC),format="%d/%m/%Y"),DT_DIGITA=as.Date((DT_DIGITA),format="%d/%m/%Y"),DT_SIN_PRI=as.Date((DT_SIN_PRI),format="%d/%m/%Y"))%>%
#   filter(DT_SIN_PRI<as.Date("2023-01-01")) # cut off at end of 2022 to ensure fully reported data
# dates_2022<-sort(as.Date(unique(SARI_2022$DT_SIN_PRI),format="%d/%m/%Y"))
# 
# SARI_2122<-rbind(SARI_2021,mutate(SARI_2022, SEM_PRI=SEM_PRI+52,SEM_NOT=SEM_NOT+52))%>%
#   mutate(delay=difftime(DT_DIGITA,DT_SIN_PRI, unit=c("days")),delay_noti=difftime(DT_NOTIFIC,DT_SIN_PRI, unit=c("days")))%>%
#   mutate(DT_COLETA=as.Date((DT_COLETA),format="%d/%m/%Y"),DT_PCR=as.Date((DT_PCR),format="%d/%m/%Y"),DT_NASC=as.Date((DT_NASC),format="%d/%m/%Y"))%>%
#   mutate(pcr_delay=difftime(DT_PCR,DT_COLETA, unit=c("days")))%>%
#   mutate(delay_diffweeks=as.numeric(difftime(DT_DIGITA,DT_SIN_PRI, unit=c("weeks")))) #
# 
# 
# # Earliest date we have data for:
# first_date<-min(SARI_2122$DT_SIN_PRI)
# 
# # We shall use the difference between digitalization date (DT_DIGITA) and date of onset of symptoms (DT_SIN_PRI).
# # Determine delay in terms of weeks not days:
# SARI_2122$delay_week<-ceiling(SARI_2122$delay_diffweeks+1/7)
# 
# # Create age groups
# SARI_2122<-SARI_2122%>%mutate(age=round(as.numeric(time_length(difftime(DT_SIN_PRI,DT_NASC), "years"))))%>%
#   mutate(age_group=cut(age, c(0, 18, 60, Inf), c("0-18", "19-60", ">60"), include.lowest=TRUE))
# 
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

x_all<-melt(apply(SARI_age_x,c(1,2),sum))
colnames(x_all)<-c("s","t","x")
x_all$s<-federal_units[x_all$s]

y_all<-melt(apply(SARI_age_y,c(1,2),sum))
colnames(y_all)<-c("s","t","y")
y_all$s<-federal_units[y_all$s]


#### Plots #####

library(geobr)
states<-read_state(year=2013)


# Set up plots:
# Observed data for SARI and COVID cases by time and region.
xy_all<-inner_join(y_all,x_all,by=c('t','s') )%>%mutate(not_cov=y-x)%>%filter(s%in%federal_units)#%>%filter(s%in%fu_abbr)
# Maximum week with observed SARI hospitalizations.
N_max<-dim(SARI_age_y)[2]
# Number of age groups.
A<-3
# Range of each age group (years).
age_group<-c('0 to 18', '19 to 60','over 60')

# Set up ggplots:
pi_plot<-ggplot()+
  theme_minimal()+labs(y=NULL, x="Prediction time difference (weeks)", title="Expected proportion of COVID-19 reported") +
  facet_wrap(~Nowcast_date, nrow=1)+
  scale_color_viridis(option="D",discrete = TRUE)+
  scale_fill_viridis(option="D",discrete = TRUE)+theme(legend.position = "none")

beta_plot<-ggplot()+
 labs(y=NULL, x="Proportion of reported SARI cases in age group", title="Effect of age distribution on COVID-poisitve proportion",
                     color = "Nowcast date", fill = "Nowcast date")+
  scale_color_viridis(option="D",discrete = TRUE)+
  scale_fill_viridis(option="D",discrete = TRUE)+
  theme_minimal()

# Make lists to save samples from model. 
pi<-pi_mean<-list()
age_beta_1<-list()
age_beta_2<-list()
age_beta_cubic1<-list()
age_beta_cubic2<-list()
cubic_beta<-list()
cubic_quant1<-list()
cubic_quant2<-list()
simulate_age_prop<-(seq(from=0,to=1,length=100))
simulate_reported<-tibble(prop=c(1,seq(from=0.95, to=0.25, length=R)),d=-30:0)
for(j in c(15,11,7,3)){
  # Plot Model Variables.
  # Beta cubic polynomials for each age group.
  age_beta_1[[j]] <- select(SARI_combined_Survivor_age[[j]],starts_with('beta_cov_0to18'))%>%as.matrix()%>%array(dim=c(n_sim,A))
   age_beta_2[[j]] <- select(SARI_combined_Survivor_age[[j]],starts_with('beta_cov_19to60'))%>%as.matrix()%>%array(dim=c(n_sim,A))
   age_beta_cubic1[[j]]<-matrix(nrow=100,ncol=n_sim)
   for(i in 1:100){
     age_beta_cubic1[[j]][i,] <-age_beta_1[[j]][,1]*scale(simulate_age_prop)[i]+  
       age_beta_1[[j]][,2]*(scale(simulate_age_prop)^2)[i]+
       age_beta_1[[j]][,3]*(scale(simulate_age_prop)^3)[i]
   }
   age_beta_cubic2[[j]]<-matrix(nrow=100,ncol=n_sim)
   for(i in 1:100){
     age_beta_cubic2[[j]][i,] <-age_beta_2[[j]][,1]*scale(simulate_age_prop)[i]+  
       age_beta_2[[j]][,2]*(scale(simulate_age_prop)^2)[i]+
     age_beta_2[[j]][,3]*(scale(simulate_age_prop)^3)[i]
   }
   
   cubic_quant1[[j]]<-age_beta_cubic1[[j]]%>%apply(c(1),quantile,c(0.025,0.5,0.975))%>%melt(varnames=c('quantile','prop'),value.name='y')%>%mutate(age=age_group[1])
   cubic_quant2[[j]]<-age_beta_cubic2[[j]]%>%apply(c(1),quantile,c(0.025,0.5,0.975))%>%melt(varnames=c('quantile','prop'),value.name='y')%>%mutate(age=age_group[2])
   cubic_beta[[j]]<-cubic_quant1[[j]]%>%
    full_join(cubic_quant2[[j]], by=c('age','quantile','prop','y'))%>%
     spread(quantile,y)%>%mutate(Nowcast_date=paste(Nowcast_dates[j]))
   
   beta_plot<-beta_plot+geom_line(data=cubic_beta[[j]],
                                   aes(x=simulate_age_prop[prop], y=`50%`,colour=Nowcast_date))+
     geom_ribbon(data=cubic_beta[[j]],
                 aes(x=simulate_age_prop[prop] , ymin=`2.5%`, ymax=`97.5%`,fill=Nowcast_date),alpha=0.25)+
     facet_wrap(~age,scales="free")
  
   
  # Reporting rate for censored COVID-19 hospitalizations.
  pi[[j]] <- select(SARI_combined_Survivor_age[[j]],starts_with('pi'))%>%as.matrix()%>%array(dim=c(n_sim,W,S))%>%
    apply(c(2,3),quantile,c(0.025,0.5,0.975))%>%melt(varnames=c('quantile','t','s'),value.name='y')%>%
    spread(quantile,y)%>%mutate(t=t+Nowcast_list[j]-W)%>%
    mutate(d=t-Nowcast_list[j])%>%mutate(Nowcast_date=paste(Nowcast_dates[j]),s=federal_units[s])%>%
    filter(d>-31)
  pi_mean[[j]]<-  pi[[j]]%>%group_by(t,d,Nowcast_date)%>%summarise(mean=mean(`50%`))
  pi_plot<-pi_plot+
    geom_line(data=simulate_reported, aes(x=d, y=prop),colour='grey20',linetype='dashed',alpha=0.75)+
    geom_line(data=pi[[j]], aes(x=d, y=`50%`,colour=Nowcast_date,shape=s),linetype='dotted',alpha=0.85)+
    geom_line(data=pi_mean[[j]], aes(x=d, y=mean,colour=Nowcast_date))
  
  
}
pi_plot
beta_plot
# Save plots:
ggsave(pi_plot,file='Plots/pi_plot.pdf',width=9,height=3.5)
ggsave(beta_plot,file='Plots/cubic_beta_age_survivor.pdf',width=9,height=3.5)

# Map plots of delta for one nowcast date:
j=12
delta[[j]] <- select(SARI_combined_Survivor_age[[j]],starts_with('delta'))%>%as.matrix()%>%array(dim=c(n_sim,S))%>%
  apply(c(2),quantile,c(0.025,0.5,0.975))%>%melt(varnames=c('quantile','s'),value.name='y')%>%
  spread(quantile,y)%>%mutate(Nowcast_date=paste(Nowcast_dates[j]),s=federal_units[s])#%>%filter(s%in%fu_abbr)
delta_all<-rbind(delta[[j]])
states_geom<-states[,c(2,6)]
colnames(states_geom)[1]<-c("s")
delta_states<-full_join(delta_all,states_geom, by="s")

limit <- max(abs(delta_all$`50%`)) * c(-1, 1)
delta_map<-ggplot() +
  geom_sf(data=delta_states,aes(geometry=geom, fill=`50%`), size=.15) +
  labs(title="Effect of total SARI incidence on the COVID-positive proportion by Brazil federative regions", 
       subtitle = paste(Nowcast_dates[j]),legend=NULL) +
  scale_fill_fermenter(palette = "PRGn",limit=limit,name=NULL,direction = -1)+ 
    theme_minimal()# +theme(panel.grid.major = element_line(colour = "white"),  axis.text = element_blank(), )

delta_map
ggsave(delta_map,file='Plots/delta_map.png',width=9,height=8,dpi = 600)

# Map plots of delta for four nowcast date:
delta_all<-rbind(delta[[15]],delta[[11]],delta[[7]],delta[[3]])
limit <- max(abs(delta_all$`50%`)) * c(-1, 1)
states_geom<-states[,c(2,6)]
colnames(states_geom)[1]<-c("s")
delta_states<-full_join(delta_all,states_geom, by="s")
delta_map_all<-  ggplot() +      
  geom_sf(data=delta_states,aes(geometry=geom, fill=`50%`), size=.15) +
  labs(title="Effect of total SARI incidence on the COVID-positive proportion by Brazil federative regions", 
       subtitle = paste(Nowcast_dates[j]),legend=NULL) +
  scale_fill_fermenter(palette = "PRGn",limit=limit,name=NULL,direction = -1)+ 
    facet_wrap(~Nowcast_date)+
  theme_minimal()
ggsave(delta_map_all,file='Plots/delta_map_all.png',width=9,height=7,dpi = 600)


