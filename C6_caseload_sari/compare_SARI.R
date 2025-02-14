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

# List of weeks to perform the nowcasts for rolling predictions experiment.
Nowcast_list<-floor(seq(from=100, to=N-D_max, length=4))
# Dates for the nowcasts.
Nowcast_dates<- as.Date(first_date)+(Nowcast_list-1)*7

# load("~/media/alba/Disk 1/OneDrive/PhD/comp 2023 paper code/FInal_code/OMEGA CODE/SARI totals linear/hazard_linear27_delay1.RData")
# load("~/media/alba/Disk 1/OneDrive/PhD/comp 2023 paper code/FInal_code/OMEGA CODE/SARI totals linear/sari_spline27.RData")
# load("~/media/alba/Disk 1/OneDrive/PhD/comp 2023 paper code/FInal_code/OMEGA CODE/SARI totals linear/survivor_sea27.RData")
# load("~/media/alba/Disk 1/OneDrive/PhD/comp 2023 paper code/FInal_code/OMEGA CODE/SARI totals linear/hazard_sea27.RData")
# load("~/media/alba/Disk 1/OneDrive/PhD/comp 2023 paper code/FInal_code/OMEGA CODE/SARI totals linear/survivor_linear27.RData")

load("~/media/alba/Disk 1/OneDrive/PhD/Case load case studies/SARI totals link/survivor_link27.RData")
load("~/media/alba/Disk 1/OneDrive/PhD/Case load case studies/SARI totals link/sari_spline27.RData")
load("~/media/alba/Disk 1/OneDrive/PhD/Case load case studies/SARI totals link/survivor_sea27.RData")

# Explore output
SARI_logy<-list()
SARI_combined_logy<-list()
SARI_linear<-list()
SARI_combined_linear<-list()
SARI_haz<-list()
SARI_combined_haz<-list()
SARI_sea<-SARI_combined_sea<-list()
SARI_haz_linear<-list()
SARI_combined_haz_linear<-list()
SARI_haz_logspline<-list()
SARI_combined_haz_logspline<-list()

SARI_combined_sea[[1]] <- as_tibble(do.call('rbind',SARI_sea_output_delta[[1]]))
n_sim <- dim(SARI_combined_sea[[1]])[1]

SARI_linear_quantiles<-SARI_haz_linear_quantiles<-list()
SARI_sea_quantiles<-SARI_haz_sea_quantiles<-list()
SARI_logysmooth_quantiles<-SARI_haz_spline_quantiles<-list()


for(j in length(Nowcast_list):1){
  #survivor models 
  SARI_sea[[j]]<-as.mcmc.list(SARI_sea_output_delta[[j]])
  SARI_combined_sea[[j]] <- as_tibble(do.call('rbind',SARI_sea[[j]]))
  SARI_sea_quantiles[[j]]<- dplyr::select(SARI_combined_sea[[j]],starts_with('y'))%>%as.matrix()%>%array(dim=c(n_sim,S,Nowcast_list[j]))%>%
    apply(c(2,3),quantile,c(0.025,0.5,0.975))%>%melt(varnames=c('quantile','s','t'),value.name='y')%>%#mutate(t=t+Nowcast_list[j]-W)%>%
    spread(quantile,y)%>%mutate(Nowcast_date=paste(Nowcast_list[j]),s=federal_units[s])
  
  SARI_linear[[j]]<-as.mcmc.list(SARI_link_output_delta[[j]])
  SARI_combined_linear[[j]] <- as_tibble(do.call('rbind',SARI_linear[[j]]))
  SARI_linear_quantiles[[j]]<-dplyr::select(SARI_combined_linear[[j]],starts_with('y'))%>%as.matrix()%>%array(dim=c(n_sim,S,Nowcast_list[j]))%>%
    apply(c(2,3),quantile,c(0.025,0.5,0.975))%>%melt(varnames=c('quantile','s','t'),value.name='y')%>%#mutate(t=t+Nowcast_list[j]-W)%>%
    spread(quantile,y)%>%mutate(Nowcast_date=paste(Nowcast_list[j]),s=federal_units[s])
  
  SARI_logy[[j]]<-as.mcmc.list(SARI_spline_output_delta[[j]])
  SARI_combined_logy[[j]] <- as_tibble(do.call('rbind',SARI_logy[[j]]))
  SARI_logysmooth_quantiles[[j]]<-dplyr::select(SARI_combined_logy[[j]],starts_with('y'))%>%as.matrix()%>%array(dim=c(n_sim,S,Nowcast_list[j]))%>%
    apply(c(2,3),quantile,c(0.025,0.5,0.975))%>%melt(varnames=c('quantile','s','t'),value.name='y')%>%#mutate(t=t+Nowcast_list[j]-W)%>%
    spread(quantile,y)%>%mutate(Nowcast_date=paste(Nowcast_list[j]),s=federal_units[s])
  
#   # hazard models 
#   SARI_haz[[j]]<-as.mcmc.list(SARI_haz_sea[[j]])
#   SARI_combined_haz[[j]] <- as_tibble(do.call('rbind',SARI_haz[[j]]))
#   SARI_haz_sea_quantiles[[j]]<-dplyr::select(SARI_combined_haz[[j]],starts_with('y'))%>%as.matrix()%>%array(dim=c(n_sim,S,Nowcast_list[j]))%>%
#     apply(c(2,3),quantile,c(0.025,0.5,0.975))%>%melt(varnames=c('quantile','s','t'),value.name='y')%>%#mutate(t=t+Nowcast_list[j]-W)%>%
#     spread(quantile,y)%>%mutate(Nowcast_date=paste(Nowcast_list[j]),s=federal_units[s])
#   
#   SARI_haz_linear[[j]]<-as.mcmc.list(SARI_linear_haz[[j]])
#   SARI_combined_haz_linear[[j]] <- as_tibble(do.call('rbind',SARI_haz_linear[[j]]))
#   SARI_haz_linear_quantiles[[j]]<-dplyr::select(SARI_combined_haz_linear[[j]],starts_with('y'))%>%as.matrix()%>%array(dim=c(n_sim,S,Nowcast_list[j]))%>%
#     apply(c(2,3),quantile,c(0.025,0.5,0.975))%>%melt(varnames=c('quantile','s','t'),value.name='y')%>%#mutate(t=t+Nowcast_list[j]-W)%>%
#     spread(quantile,y)%>%mutate(Nowcast_date=paste(Nowcast_list[j]),s=federal_units[s])
# # 
# #   SARI_haz_logspline[[j]]<-as.mcmc.list(SARI_spline_haz[[j]])
# #   SARI_combined_haz_logspline[[j]] <- as_tibble(do.call('rbind',SARI_haz_logspline[[j]]))
# #   SARI_haz_spline_quantiles[[j]]<-dplyr::select(SARI_combined_haz_logspline[[j]],starts_with('y'))%>%as.matrix()%>%array(dim=c(n_sim,S,Nowcast_list[j]))%>%
# #     apply(c(2,3),quantile,c(0.025,0.5,0.975))%>%melt(varnames=c('quantile','s','t'),value.name='y')%>%#mutate(t=t+Nowcast_list[j]-W)%>%
# #     spread(quantile,y)%>%mutate(Nowcast_date=paste(Nowcast_list[j]),s=federal_units[s])

  
}


#compare models 
library(ggh4x)


# 
# # Survivor and hazard sea 
# SARI_nowcasts_compare<-list()
# for(j in 1:length(Nowcast_list)){
#   SARI_nowcasts_compare[[j]] <- rbind(
#     full_join(SARI_haz_sea_quantiles[[j]],filter(y_all, t%in%c((1):Nowcast_list[j])),by=c('s','t'))%>%mutate(m=c('SARI hazard'),d=t-Nowcast_list[j],now=j),
#     full_join(SARI_sea_quantiles[[j]],filter(y_all, t%in%c((1):Nowcast_list[j])),by=c('s','t'))%>%mutate(m=c('SARI survivor'),d=t-Nowcast_list[j],now=j)
#   )%>%filter(d>-D)
# }
# SARI_nowcast_dmax<-array(unlist(SARI_nowcasts_compare), dim=c(S*D*2,ncol(SARI_nowcasts_compare[[1]]),length(Nowcast_list)))%>%
#   apply(c(2),c)%>%as_tibble()
# colnames(SARI_nowcast_dmax)<-colnames(SARI_nowcasts_compare[[1]])   
# SARI_nowcast_dmax<-SARI_nowcast_dmax%>%mutate(t=as.numeric(t),`2.5%`=as.numeric(`2.5%`),`50%`=as.numeric(`50%`),`97.5%`=as.numeric(`97.5%`),
#                                                   m=factor(m,levels=c('SARI hazard','SARI survivor'
#                                                   )),y=as.numeric(y),
#                                                   d=as.numeric(d),now=as.numeric(now))
# SARI_nowcast_dmax_summary<-SARI_nowcast_dmax%>%group_by(m,d,s)%>%
#   summarise(Coverage=mean(y>=`2.5%`&y<=`97.5%`),
#             `Prediction interval width`=mean(`97.5%`-`2.5%`),
#             `Mean absolute error`=mean(abs(`50%`-y)))%>%
#   pivot_longer(cols=c('Coverage','Mean absolute error','Prediction interval width'),names_to = 'metric')%>%
#   mutate(m=factor(m,levels=c('SARI hazard','SARI survivor'
#   )),
#   metric=factor(metric,levels=c('Mean absolute error','Prediction interval width','Coverage')))
# 
# 
# #plot_data_SARI<-filter(SARI_nowcast_dmax_summary, d==0)%>%pivot_wider(names_from=c('m'))%>%filter(metric=="Mean Absolute Error") 
# #ggplot(plot_data_SARI)+   geom_text(aes(x=`CLR`,y=`CLR Age`/`CLR`,label=s))+geom_hline(yintercept = 1)
# SARI_nowcast_dmax_summary<-filter(SARI_nowcast_dmax_summary, s!="ES")
# SARI_nowcast_dmax_summary_average<-SARI_nowcast_dmax_summary%>%group_by(m,d,metric)%>%
#   summarise(value.mean=mean(value), value.median=median(value))
# SARI_nowcast_dmax_summary_average$value.mean[SARI_nowcast_dmax_summary_average$metric%in%c("Prediction Interval Width","Mean Absolute Error")]<-NA
# SARI_nowcast_dmax_summary_average$value.median[SARI_nowcast_dmax_summary_average$metric%in%c("Coverage")]<-NA
# SARI_nowcast_dmax_summary_average<-SARI_nowcast_dmax_summary_average%>%group_by(m,d,metric)%>%mutate(value=sum(value.mean,value.median,na.rm=TRUE))
# 
# # Plot of prediction performance metrics.
# sea_hazard_survivor <- ggplot(SARI_nowcast_dmax_summary_average)+
#   geom_line(aes(x=d,y=value,colour=m),alpha=0.6)+
#   geom_point(aes(x=d,y=value,colour=m,shape=m),alpha=0.6)+
#   facet_grid2(~metric, scales = c("free"),independent = c("all"))+
#   scale_colour_manual(name=NULL,values = c( "#0072B2", "#D55E00"))+
#   scale_shape_manual(name=NULL,values=c(15,3))+
#   scale_x_continuous()+
#   labs(x='Prediction time difference (weeks)',y=NULL,
#        title='Performance metrics from the rolling prediction experiment')+
#   theme_minimal()+
#   theme(legend.position = "bottom")
# sea_hazard_survivor

# Survivor sea and linear 
SARI_nowcasts_compare<-list()
for(j in 1:length(Nowcast_list)){
  SARI_nowcasts_compare[[j]] <- rbind(
    full_join(SARI_sea_quantiles[[j]],filter(y_all, t%in%c((1):Nowcast_list[j])),by=c('s','t'))%>%mutate(m=c('SARI survivor'),d=t-Nowcast_list[j],now=j),
    full_join(SARI_linear_quantiles[[j]],filter(y_all, t%in%c((1):Nowcast_list[j])),by=c('s','t'))%>%mutate(m=c('SARI survivor incidence-delay'),d=t-Nowcast_list[j],now=j)
  )%>%filter(d>-D)
}
SARI_nowcast_dmax<-array(unlist(SARI_nowcasts_compare), dim=c(S*D*2,ncol(SARI_nowcasts_compare[[1]]),length(Nowcast_list)))%>%
  apply(c(2),c)%>%as_tibble()
colnames(SARI_nowcast_dmax)<-colnames(SARI_nowcasts_compare[[1]])   
SARI_nowcast_dmax<-SARI_nowcast_dmax%>%mutate(t=as.numeric(t),`2.5%`=as.numeric(`2.5%`),`50%`=as.numeric(`50%`),`97.5%`=as.numeric(`97.5%`),
                                                  m=factor(m,levels=c('SARI survivor','SARI survivor incidence-delay'
                                                  )),y=as.numeric(y),
                                                  d=as.numeric(d),now=as.numeric(now))
SARI_nowcast_dmax_summary<-SARI_nowcast_dmax%>%group_by(m,d,s)%>%
  summarise(Coverage=mean(y>=`2.5%`&y<=`97.5%`),
            `Prediction interval width`=mean(`97.5%`-`2.5%`),
            `Mean absolute error`=mean(abs(`50%`-y)))%>%
  pivot_longer(cols=c('Coverage','Mean absolute error','Prediction interval width'),names_to = 'metric')%>%
  mutate(m=factor(m,levels=c('SARI survivor','SARI survivor incidence-delay'
  )),
  metric=factor(metric,levels=c('Mean absolute error','Prediction interval width','Coverage')))


#plot_data_SARI<-filter(SARI_nowcast_dmax_summary, d==0)%>%pivot_wider(names_from=c('m'))%>%filter(metric=="Mean Absolute Error") 
#ggplot(plot_data_SARI)+   geom_text(aes(x=`CLR`,y=`CLR Age`/`CLR`,label=s))+geom_hline(yintercept = 1)
SARI_nowcast_dmax_summary<-filter(SARI_nowcast_dmax_summary, s!="ES")
SARI_nowcast_dmax_summary_average<-SARI_nowcast_dmax_summary%>%group_by(m,d,metric)%>%
  summarise(value.mean=mean(value), value.median=median(value))
SARI_nowcast_dmax_summary_average$value.mean[SARI_nowcast_dmax_summary_average$metric%in%c("Prediction Interval Width","Mean Absolute Error")]<-NA
SARI_nowcast_dmax_summary_average$value.median[SARI_nowcast_dmax_summary_average$metric%in%c("Coverage")]<-NA
SARI_nowcast_dmax_summary_average<-SARI_nowcast_dmax_summary_average%>%group_by(m,d,metric)%>%mutate(value=sum(value.mean,value.median,na.rm=TRUE))

# Plot of prediction performance metrics.
SARI_nowcast_dmax_summary_average<-SARI_nowcast_dmax_summary_average%>%mutate(yintercept=case_when(metric=="Coverage"~0.95))
sea_link_survivor <- ggplot(SARI_nowcast_dmax_summary_average)+
  geom_hline(aes(yintercept = yintercept), colour="darkgrey")+
  geom_line(aes(x=d,y=value,colour=m),alpha=0.6)+
  geom_point(aes(x=d,y=value,colour=m,shape=m),alpha=0.6)+
  facet_grid2(~metric, scales = c("free"),independent = c("all"))+
  scale_colour_manual(name=NULL,values = c( "#0072B2", "#D55E00"))+
  scale_shape_manual(name=NULL,values=c(15,3))+
  scale_x_continuous()+
  labs(x='Prediction time difference (weeks)',y=NULL,
       title="Brazil SARI hospitalisations",
       subtitle='Performance metrics from the rolling prediction experiment')+
  theme_minimal()+
  theme(legend.position = "bottom",text = element_text(size=12), plot.title = element_text(size=16),
        plot.subtitle=element_text(size=14), legend.title = element_text(size=14),
        strip.text = element_text(size=14))
sea_link_survivor
ggsave(sea_link_survivor, file="plots/SARI_sea_link_survivor.pdf", width=9, height=4)
# 
# 
# 
# # hazard linear and sea
# SARI_nowcasts_compare<-list()
# for(j in 1:length(Nowcast_list)){
#   SARI_nowcasts_compare[[j]] <- rbind(
#     full_join(SARI_haz_sea_quantiles[[j]],filter(y_all, t%in%c((1):Nowcast_list[j])),by=c('s','t'))%>%mutate(m=c('SARI hazard'),d=t-Nowcast_list[j],now=j),
#     full_join(SARI_haz_linear_quantiles[[j]],filter(y_all, t%in%c((1):Nowcast_list[j])),by=c('s','t'))%>%mutate(m=c('SARI hazard linear'),d=t-Nowcast_list[j],now=j)
#   )%>%filter(d>-D)
# }
# SARI_nowcast_dmax<-array(unlist(SARI_nowcasts_compare), dim=c(S*D*2,ncol(SARI_nowcasts_compare[[1]]),length(Nowcast_list)))%>%
#   apply(c(2),c)%>%as_tibble()
# colnames(SARI_nowcast_dmax)<-colnames(SARI_nowcasts_compare[[1]])   
# SARI_nowcast_dmax<-SARI_nowcast_dmax%>%mutate(t=as.numeric(t),`2.5%`=as.numeric(`2.5%`),`50%`=as.numeric(`50%`),`97.5%`=as.numeric(`97.5%`),
#                                                   m=factor(m,levels=c('SARI hazard','SARI hazard linear'
#                                                   )),y=as.numeric(y),
#                                                   d=as.numeric(d),now=as.numeric(now))
# SARI_nowcast_dmax_summary<-SARI_nowcast_dmax%>%group_by(m,d,s)%>%
#   summarise(Coverage=mean(y>=`2.5%`&y<=`97.5%`),
#             `Prediction interval width`=mean(`97.5%`-`2.5%`),
#             `Mean absolute error`=mean(abs(`50%`-y)))%>%
#   pivot_longer(cols=c('Coverage','Mean absolute error','Prediction interval width'),names_to = 'metric')%>%
#   mutate(m=factor(m,levels=c('SARI hazard','SARI hazard linear'
#   )),
#   metric=factor(metric,levels=c('Mean absolute error','Prediction interval width','Coverage')))
# 
# 
# #plot_data_SARI<-filter(SARI_nowcast_dmax_summary, d==0)%>%pivot_wider(names_from=c('m'))%>%filter(metric=="Mean Absolute Error") 
# #ggplot(plot_data_SARI)+   geom_text(aes(x=`CLR`,y=`CLR Age`/`CLR`,label=s))+geom_hline(yintercept = 1)
# SARI_nowcast_dmax_summary<-filter(SARI_nowcast_dmax_summary, s!="ES")
# SARI_nowcast_dmax_summary_average<-SARI_nowcast_dmax_summary%>%group_by(m,d,metric)%>%
#   summarise(value.mean=mean(value), value.median=median(value))
# SARI_nowcast_dmax_summary_average$value.mean[SARI_nowcast_dmax_summary_average$metric%in%c("Prediction Interval Width","Mean Absolute Error")]<-NA
# SARI_nowcast_dmax_summary_average$value.median[SARI_nowcast_dmax_summary_average$metric%in%c("Coverage")]<-NA
# SARI_nowcast_dmax_summary_average<-SARI_nowcast_dmax_summary_average%>%group_by(m,d,metric)%>%mutate(value=sum(value.mean,value.median,na.rm=TRUE))
# 
# # Plot of prediction performance metrics.
# linear_sea_hazard <- ggplot(SARI_nowcast_dmax_summary_average)+
#   geom_line(aes(x=d,y=value,colour=m),alpha=0.6)+
#   geom_point(aes(x=d,y=value,colour=m,shape=m),alpha=0.6)+
#   facet_grid2(~metric, scales = c("free"),independent = c("all"))+
#   scale_colour_manual(name=NULL,values = c( "#0072B2", "#D55E00"))+
#   scale_shape_manual(name=NULL,values=c(15,3))+
#   scale_x_continuous()+
#   labs(x='Prediction time difference (weeks)',y=NULL,
#        title='Performance metrics from the rolling prediction experiment')+
#   theme_minimal()+
#   theme(legend.position = "bottom")
# linear_sea_hazard
# 
# 
# # Survivor linear and hazard linear 
# SARI_nowcasts_compare<-list()
# for(j in 1:length(Nowcast_list)){
#   SARI_nowcasts_compare[[j]] <- rbind(
#     full_join(SARI_linear_quantiles[[j]],filter(y_all, t%in%c((1):Nowcast_list[j])),by=c('s','t'))%>%mutate(m=c('SARI survivor linear'),d=t-Nowcast_list[j],now=j),
#     full_join(SARI_haz_linear_quantiles[[j]],filter(y_all, t%in%c((1):Nowcast_list[j])),by=c('s','t'))%>%mutate(m=c('SARI hazard linear'),d=t-Nowcast_list[j],now=j)
#   )%>%filter(d>-D)
# }
# SARI_nowcast_dmax<-array(unlist(SARI_nowcasts_compare), dim=c(S*D*2,ncol(SARI_nowcasts_compare[[1]]),length(Nowcast_list)))%>%
#   apply(c(2),c)%>%as_tibble()
# colnames(SARI_nowcast_dmax)<-colnames(SARI_nowcasts_compare[[1]])   
# SARI_nowcast_dmax<-SARI_nowcast_dmax%>%mutate(t=as.numeric(t),`2.5%`=as.numeric(`2.5%`),`50%`=as.numeric(`50%`),`97.5%`=as.numeric(`97.5%`),
#                                                   m=factor(m,levels=c('SARI survivor linear','SARI hazard linear'
#                                                   )),y=as.numeric(y),
#                                                   d=as.numeric(d),now=as.numeric(now))
# SARI_nowcast_dmax_summary<-SARI_nowcast_dmax%>%group_by(m,d,s)%>%
#   summarise(Coverage=mean(y>=`2.5%`&y<=`97.5%`),
#             `Prediction interval width`=mean(`97.5%`-`2.5%`),
#             `Mean absolute error`=mean(abs(`50%`-y)))%>%
#   pivot_longer(cols=c('Coverage','Mean absolute error','Prediction interval width'),names_to = 'metric')%>%
#   mutate(m=factor(m,levels=c('SARI survivor linear','SARI hazard linear'
#   )),
#   metric=factor(metric,levels=c('Mean absolute error','Prediction interval width','Coverage')))
# 
# 
# #plot_data_SARI<-filter(SARI_nowcast_dmax_summary, d==0)%>%pivot_wider(names_from=c('m'))%>%filter(metric=="Mean Absolute Error") 
# #ggplot(plot_data_SARI)+   geom_text(aes(x=`CLR`,y=`CLR Age`/`CLR`,label=s))+geom_hline(yintercept = 1)
# SARI_nowcast_dmax_summary<-filter(SARI_nowcast_dmax_summary, s!="ES")
# SARI_nowcast_dmax_summary_average<-SARI_nowcast_dmax_summary%>%group_by(m,d,metric)%>%
#   summarise(value.mean=mean(value), value.median=median(value))
# SARI_nowcast_dmax_summary_average$value.mean[SARI_nowcast_dmax_summary_average$metric%in%c("Prediction Interval Width","Mean Absolute Error")]<-NA
# SARI_nowcast_dmax_summary_average$value.median[SARI_nowcast_dmax_summary_average$metric%in%c("Coverage")]<-NA
# SARI_nowcast_dmax_summary_average<-SARI_nowcast_dmax_summary_average%>%group_by(m,d,metric)%>%mutate(value=sum(value.mean,value.median,na.rm=TRUE))
# 
# # Plot of prediction performance metrics.
# linear_survivor_hazard <- ggplot(SARI_nowcast_dmax_summary_average)+
#   geom_line(aes(x=d,y=value,colour=m),alpha=0.6)+
#   geom_point(aes(x=d,y=value,colour=m,shape=m),alpha=0.6)+
#   facet_grid2(~metric, scales = c("free"),independent = c("all"))+
#   scale_colour_manual(name=NULL,values = c( "#0072B2", "#D55E00"))+
#   scale_shape_manual(name=NULL,values=c(15,3))+
#   scale_x_continuous()+
#   labs(x='Prediction time difference (weeks)',y=NULL,
#        title='Performance metrics from the rolling prediction experiment')+
#   theme_minimal()+
#   theme(legend.position = "bottom")
# linear_survivor_hazard


# Survivor linear and survivor spline
SARI_nowcasts_compare<-list()
for(j in 1:length(Nowcast_list)){
  SARI_nowcasts_compare[[j]] <- rbind(
    full_join(SARI_linear_quantiles[[j]],filter(y_all, t%in%c((1):Nowcast_list[j])),by=c('s','t'))%>%mutate(m=c('SARI survivor linear'),d=t-Nowcast_list[j],now=j),
    full_join(SARI_logysmooth_quantiles[[j]],filter(y_all, t%in%c((1):Nowcast_list[j])),by=c('s','t'))%>%mutate(m=c('SARI survivor spline'),d=t-Nowcast_list[j],now=j)
  )%>%filter(d>-D)
}
SARI_nowcast_dmax<-array(unlist(SARI_nowcasts_compare), dim=c(S*D*2,ncol(SARI_nowcasts_compare[[1]]),length(Nowcast_list)))%>%
  apply(c(2),c)%>%as_tibble()
colnames(SARI_nowcast_dmax)<-colnames(SARI_nowcasts_compare[[1]])   
SARI_nowcast_dmax<-SARI_nowcast_dmax%>%mutate(t=as.numeric(t),`2.5%`=as.numeric(`2.5%`),`50%`=as.numeric(`50%`),`97.5%`=as.numeric(`97.5%`),
                                              m=factor(m,levels=c('SARI survivor linear','SARI survivor spline'
                                              )),y=as.numeric(y),
                                              d=as.numeric(d),now=as.numeric(now))
SARI_nowcast_dmax_summary<-SARI_nowcast_dmax%>%group_by(m,d,s)%>%
  summarise(Coverage=mean(y>=`2.5%`&y<=`97.5%`),
            `Prediction interval width`=mean(`97.5%`-`2.5%`),
            `Mean absolute error`=mean(abs(`50%`-y)))%>%
  pivot_longer(cols=c('Coverage','Mean absolute error','Prediction interval width'),names_to = 'metric')%>%
  mutate(m=factor(m,levels=c('SARI survivor linear','SARI survivor spline'
  )),
  metric=factor(metric,levels=c('Mean absolute error','Prediction interval width','Coverage')))


#plot_data_SARI<-filter(SARI_nowcast_dmax_summary, d==0)%>%pivot_wider(names_from=c('m'))%>%filter(metric=="Mean Absolute Error") 
#ggplot(plot_data_SARI)+   geom_text(aes(x=`CLR`,y=`CLR Age`/`CLR`,label=s))+geom_hline(yintercept = 1)
SARI_nowcast_dmax_summary<-filter(SARI_nowcast_dmax_summary, s!="ES")
SARI_nowcast_dmax_summary_average<-SARI_nowcast_dmax_summary%>%group_by(m,d,metric)%>%
  summarise(value.mean=mean(value), value.median=median(value))
SARI_nowcast_dmax_summary_average$value.mean[SARI_nowcast_dmax_summary_average$metric%in%c("Prediction Interval Width","Mean Absolute Error")]<-NA
SARI_nowcast_dmax_summary_average$value.median[SARI_nowcast_dmax_summary_average$metric%in%c("Coverage")]<-NA
SARI_nowcast_dmax_summary_average<-SARI_nowcast_dmax_summary_average%>%group_by(m,d,metric)%>%mutate(value=sum(value.mean,value.median,na.rm=TRUE))

# Plot of prediction performance metrics.
linear_spline_survivor <- ggplot(SARI_nowcast_dmax_summary_average)+
  geom_line(aes(x=d,y=value,colour=m),alpha=0.6)+
  geom_point(aes(x=d,y=value,colour=m,shape=m),alpha=0.6)+
  facet_grid2(~metric, scales = c("free"),independent = c("all"))+
  scale_colour_manual(name=NULL,values = c( "#0072B2", "#D55E00"))+
  scale_shape_manual(name=NULL,values=c(15,3))+
  scale_x_continuous()+
  labs(x='Prediction time difference (weeks)',y=NULL,
       title='Performance metrics from the rolling prediction experiment')+
  theme_minimal()+
  theme(legend.position = "bottom")
linear_spline_survivor

# Survivor sea and survivor spline
SARI_nowcasts_compare<-list()
for(j in 1:length(Nowcast_list)){
  SARI_nowcasts_compare[[j]] <- rbind(
    full_join(SARI_sea_quantiles[[j]],filter(y_all, t%in%c((1):Nowcast_list[j])),by=c('s','t'))%>%mutate(m=c('SARI survivor'),d=t-Nowcast_list[j],now=j),
    full_join(SARI_logysmooth_quantiles[[j]],filter(y_all, t%in%c((1):Nowcast_list[j])),by=c('s','t'))%>%mutate(m=c('SARI survivor spline'),d=t-Nowcast_list[j],now=j)
  )%>%filter(d>-D)
}
SARI_nowcast_dmax<-array(unlist(SARI_nowcasts_compare), dim=c(S*D*2,ncol(SARI_nowcasts_compare[[1]]),length(Nowcast_list)))%>%
  apply(c(2),c)%>%as_tibble()
colnames(SARI_nowcast_dmax)<-colnames(SARI_nowcasts_compare[[1]])   
SARI_nowcast_dmax<-SARI_nowcast_dmax%>%mutate(t=as.numeric(t),`2.5%`=as.numeric(`2.5%`),`50%`=as.numeric(`50%`),`97.5%`=as.numeric(`97.5%`),
                                              m=factor(m,levels=c('SARI survivor','SARI survivor spline'
                                              )),y=as.numeric(y),
                                              d=as.numeric(d),now=as.numeric(now))
SARI_nowcast_dmax_summary<-SARI_nowcast_dmax%>%group_by(m,d,s)%>%
  summarise(Coverage=mean(y>=`2.5%`&y<=`97.5%`),
            `Prediction interval width`=mean(`97.5%`-`2.5%`),
            `Mean absolute error`=mean(abs(`50%`-y)))%>%
  pivot_longer(cols=c('Coverage','Mean absolute error','Prediction interval width'),names_to = 'metric')%>%
  mutate(m=factor(m,levels=c('SARI survivor','SARI survivor spline'
  )),
  metric=factor(metric,levels=c('Mean absolute error','Prediction interval width','Coverage')))


#plot_data_SARI<-filter(SARI_nowcast_dmax_summary, d==0)%>%pivot_wider(names_from=c('m'))%>%filter(metric=="Mean Absolute Error") 
#ggplot(plot_data_SARI)+   geom_text(aes(x=`CLR`,y=`CLR Age`/`CLR`,label=s))+geom_hline(yintercept = 1)
#SARI_nowcast_dmax_summary<-filter(SARI_nowcast_dmax_summary, s!="ES")
SARI_nowcast_dmax_summary_average<-SARI_nowcast_dmax_summary%>%group_by(m,d,metric)%>%
  summarise(value.mean=mean(value), value.median=median(value))
SARI_nowcast_dmax_summary_average$value.mean[SARI_nowcast_dmax_summary_average$metric%in%c("Prediction Interval Width","Mean Absolute Error")]<-NA
SARI_nowcast_dmax_summary_average$value.median[SARI_nowcast_dmax_summary_average$metric%in%c("Coverage")]<-NA
SARI_nowcast_dmax_summary_average<-SARI_nowcast_dmax_summary_average%>%group_by(m,d,metric)%>%mutate(value=sum(value.mean,value.median,na.rm=TRUE))

# Plot of prediction performance metrics.
sea_spline_survivor <- ggplot(SARI_nowcast_dmax_summary_average)+
  geom_line(aes(x=d,y=value,colour=m),alpha=0.6)+
  geom_point(aes(x=d,y=value,colour=m,shape=m),alpha=0.6)+
  facet_grid2(~metric, scales = c("free"),independent = c("all"))+
  scale_colour_manual(name=NULL,values = c( "#0072B2", "#D55E00"))+
  scale_shape_manual(name=NULL,values=c(15,3))+
  scale_x_continuous()+
  labs(x='Prediction time difference (weeks)',y=NULL,
       title='Performance metrics from the rolling prediction experiment')+
  theme_minimal()+
  theme(legend.position = "bottom")
sea_spline_survivor

# COMPARE FOR EACH REGION 

# Survivor sea and survivor spline
SARI_nowcasts_compare<-list()
for(j in 1:length(Nowcast_list)){
  SARI_nowcasts_compare[[j]] <- rbind(
    full_join(SARI_sea_quantiles[[j]],filter(y_all, t%in%c((1):Nowcast_list[j])),by=c('s','t'))%>%mutate(m=c('SARI survivor'),d=t-Nowcast_list[j],now=j),
    full_join(SARI_logysmooth_quantiles[[j]],filter(y_all, t%in%c((1):Nowcast_list[j])),by=c('s','t'))%>%mutate(m=c('SARI survivor spline link'),d=t-Nowcast_list[j],now=j)
  )%>%filter(d>-D)
}
SARI_nowcast_dmax<-array(unlist(SARI_nowcasts_compare), dim=c(S*D*2,ncol(SARI_nowcasts_compare[[1]]),length(Nowcast_list)))%>%
  apply(c(2),c)%>%as_tibble()
colnames(SARI_nowcast_dmax)<-colnames(SARI_nowcasts_compare[[1]])   
SARI_nowcast_dmax<-SARI_nowcast_dmax%>%mutate(t=as.numeric(t),`2.5%`=as.numeric(`2.5%`),`50%`=as.numeric(`50%`),`97.5%`=as.numeric(`97.5%`),
                                              m=factor(m,levels=c('SARI survivor','SARI survivor spline link'
                                              )),y=as.numeric(y),
                                              d=as.numeric(d),now=as.numeric(now))
SARI_nowcast_dmax_summary<-SARI_nowcast_dmax%>%group_by(m,d,s)%>%
  summarise(Coverage=mean(y>=`2.5%`&y<=`97.5%`),
            `Prediction interval width`=mean(`97.5%`-`2.5%`),
            `Mean absolute error`=mean(abs(`50%`-y)))%>%
  pivot_longer(cols=c('Coverage','Mean absolute error','Prediction interval width'),names_to = 'metric')%>%
  mutate(m=factor(m,levels=c('SARI survivor','SARI survivor spline link'
  )),
  metric=factor(metric,levels=c('Mean absolute error','Prediction interval width','Coverage')))

SARI_nowcast_dmax_summary_average<-SARI_nowcast_dmax_summary%>%group_by(m,d,metric)%>%
  summarise(value.mean=mean(value), value.median=median(value))
SARI_nowcast_dmax_summary_average$value.mean[SARI_nowcast_dmax_summary_average$metric%in%c("Prediction Interval Width","Mean Absolute Error")]<-NA
SARI_nowcast_dmax_summary_average$value.median[SARI_nowcast_dmax_summary_average$metric%in%c("Coverage")]<-NA
SARI_nowcast_dmax_summary_average<-SARI_nowcast_dmax_summary_average%>%group_by(m,d,metric)%>%mutate(value=sum(value.mean,value.median,na.rm=TRUE))

#plot_data_SARI<-filter(SARI_nowcast_dmax_summary, d==0)%>%pivot_wider(names_from=c('m'))%>%filter(metric=="Mean Absolute Error") 
#ggplot(plot_data_SARI)+   geom_text(aes(x=`CLR`,y=`CLR Age`/`CLR`,label=s))+geom_hline(yintercept = 1)

# Plot of prediction performance metrics.
sea_spline_survivor <- ggplot(SARI_nowcast_dmax_summary_average)+
  geom_line(aes(x=d,y=value,colour=m),alpha=0.6)+
  geom_point(aes(x=d,y=value,colour=m,shape=m),alpha=0.6)+
  facet_grid2(~metric, scales = c("free"),independent = c("all"))+
  scale_colour_manual(name=NULL,values = c( "#0072B2", "#ed9add"))+
  scale_shape_manual(name=NULL,values=c(15,3))+
  scale_x_continuous()+
  labs(x='Prediction time difference (weeks)',y=NULL,
       title="Brazil SARI hospitalisations",
       subtitle='Performance metrics from the rolling prediction experiment')+
  theme_minimal()+
  theme(legend.position = "bottom",text = element_text(size=12), plot.title = element_text(size=16),
        plot.subtitle=element_text(size=14), legend.title = element_text(size=14),
        strip.text = element_text(size=14))

ggsave(sea_spline_survivor, file="SARI_sea_spline_survivor.pdf", width=9, height=4)

sea_spline_survivor_notitle <- ggplot(SARI_nowcast_dmax_summary_average)+
  geom_line(aes(x=d,y=value,colour=m),alpha=0.6)+
  geom_point(aes(x=d,y=value,colour=m,shape=m),alpha=0.6)+
  facet_grid2(~metric, scales = c("free"),independent = c("all"))+
  scale_colour_manual(name=NULL,values = c( "#0072B2", "#ed9add"))+
  scale_shape_manual(name=NULL,values=c(15,4))+
  scale_x_continuous()+
  labs(x='Prediction time difference (weeks)',y=NULL,
       title=NULL,
       subtitle=NULL)+
  theme_minimal()+
  theme(legend.position = "bottom",text = element_text(size=12), plot.title = element_text(size=16),
        plot.subtitle=element_text(size=14), legend.title = element_text(size=14),
        strip.text = element_text(size=14))
ggsave(sea_spline_survivor_notitle, file="SARI_sea_spline_survivor_notitle.pdf", width=9, height=3.5)


SARI_nowcast_dmax_summary<-filter(SARI_nowcast_dmax_summary, metric=='Mean absolute error')


# Survivor sea and survivor spline
SARI_nowcasts_compare<-list()
for(j in 1:length(Nowcast_list)){
  SARI_nowcasts_compare[[j]] <- rbind(
    full_join(SARI_sea_quantiles[[j]],filter(y_all, t%in%c((1):Nowcast_list[j])),by=c('s','t'))%>%mutate(m=c('SARI survivor'),d=t-Nowcast_list[j],now=j),
    full_join(SARI_linear_quantiles[[j]],filter(y_all, t%in%c((1):Nowcast_list[j])),by=c('s','t'))%>%mutate(m=c('SARI survivor incidence-delay'),d=t-Nowcast_list[j],now=j),
    full_join(SARI_logysmooth_quantiles[[j]],filter(y_all, t%in%c((1):Nowcast_list[j])),by=c('s','t'))%>%mutate(m=c('SARI survivor spline link'),d=t-Nowcast_list[j],now=j)
  )%>%filter(d>-D)
}
SARI_nowcast_dmax<-array(unlist(SARI_nowcasts_compare), dim=c(S*D*3,ncol(SARI_nowcasts_compare[[1]]),length(Nowcast_list)))%>%
  apply(c(2),c)%>%as_tibble()
colnames(SARI_nowcast_dmax)<-colnames(SARI_nowcasts_compare[[1]])   
SARI_nowcast_dmax<-SARI_nowcast_dmax%>%mutate(t=as.numeric(t),`2.5%`=as.numeric(`2.5%`),`50%`=as.numeric(`50%`),`97.5%`=as.numeric(`97.5%`),
                                              m=factor(m,levels=c('SARI survivor','SARI survivor incidence-delay','SARI survivor spline link'
                                              )),y=as.numeric(y),
                                              d=as.numeric(d),now=as.numeric(now))
SARI_nowcast_dmax_summary<-SARI_nowcast_dmax%>%group_by(m,d,s)%>%
  summarise(Coverage=mean(y>=`2.5%`&y<=`97.5%`),
            `Prediction interval width`=mean(`97.5%`-`2.F5%`),
            `Mean absolute error`=mean(abs(`50%`-y)))%>%
  pivot_longer(cols=c('Coverage','Mean absolute error','Prediction interval width'),names_to = 'metric')%>%
  mutate(m=factor(m,levels=c('SARI survivor','SARI survivor incidence-delay','SARI survivor spline link'
  )),
  metric=factor(metric,levels=c('Mean absolute error','Prediction interval width','Coverage')))

SARI_nowcast_dmax_summary_average<-SARI_nowcast_dmax_summary%>%group_by(m,d,metric)%>%
  summarise(value.mean=mean(value), value.median=median(value))
SARI_nowcast_dmax_summary_average$value.mean[SARI_nowcast_dmax_summary_average$metric%in%c("Prediction Interval Width","Mean Absolute Error")]<-NA
SARI_nowcast_dmax_summary_average$value.median[SARI_nowcast_dmax_summary_average$metric%in%c("Coverage")]<-NA
SARI_nowcast_dmax_summary_average<-SARI_nowcast_dmax_summary_average%>%group_by(m,d,metric)%>%mutate(value=sum(value.mean,value.median,na.rm=TRUE))

# Plot of prediction performance metrics.
sea_spline_linear <- ggplot(SARI_nowcast_dmax_summary_average)+
  geom_line(aes(x=d,y=value,colour=m),alpha=0.6)+
  geom_point(aes(x=d,y=value,colour=m,shape=m),alpha=0.6)+
  facet_grid2(~metric, scales = c("free"),independent = c("all"))+
  scale_colour_manual(name=NULL,values = c( "#0072B2", "#D55E00","#ed9add" ))+
  scale_shape_manual(name=NULL,values=c(15,3,4))+
  scale_x_continuous()+
  labs(x='Prediction time difference (weeks)',y=NULL,
       title="Brazil SARI hospitalisations",
       subtitle='Performance metrics from the rolling prediction experiment')+
  theme_minimal()+
  theme(legend.position = "bottom",text = element_text(size=12), plot.title = element_text(size=16),
        plot.subtitle=element_text(size=14), legend.title = element_text(size=14),
        strip.text = element_text(size=14))

ggsave(sea_spline_linear, file="SARI_all_survivor.pdf", width=9, height=4)


# Plot of prediction performance metrics.
mae_survivor <- ggplot(SARI_nowcast_dmax_summary)+
  geom_line(aes(x=d,y=value,colour=m),alpha=0.6)+
  geom_point(aes(x=d,y=value,colour=m,shape=m),alpha=0.6)+
  facet_wrap(s~metric, scales = 'free')+
  scale_colour_manual(name=NULL,values = c( "#0072B2", "#D55E00"))+
  scale_shape_manual(name=NULL,values=c(15,3))+
  scale_x_continuous()+
  labs(x='Prediction time difference (weeks)',y=NULL,
       title='Performance metrics from the rolling prediction experiment')+
  theme_minimal()+
  theme(legend.position = "bottom")
mae_survivor


# compare full delta to prediction experiment delta 
load("~/media/alba/Disk 1/OneDrive/PhD/Case load case studies/SARI totals link/survivor_FULL_linear.RData")
load("~/media/alba/Disk 1/OneDrive/PhD/Case load case studies/SARI totals link/survivor_link27.RData")
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
  delta_slope[[j]]<-dplyr::select(SARI_combined_Survivor_linear_x[[j]],starts_with('delta'))%>%as.matrix()%>%array(dim=c(n_sim))%>%
    quantile(c(0.025,0.5,0.975))%>%melt(varnames=c('quantile'),value.name='y')%>%mutate(s=j,quantile=c(0.025,0.5,0.975))%>%
    spread(quantile,y)
  
}
#delta matrix
delta_matrix<-matrix(unlist(delta_slope), byrow = TRUE, ncol = 4)%>%as_tibble()
colnames(delta_matrix)<-c('s','2.5%','50%','97.5%')
delta_matrix$s_abb<-factor(federal_units[delta_matrix$s],levels=federal_units[delta_matrix$s])
delta_matrix<-full_join(delta_matrix,state_guide[,2:3], by='s_abb')
delta_matrix$s_full<-factor(delta_matrix$s_full, levels=delta_matrix$s_full)
delta_matrix<-delta_matrix%>%mutate(model='Parameter inference experiment')

# prediction delta:
n_sim<-dim(SARI_link_output_delta[[1]]$chain1)[1]*2
# Survivor age model
SARI_Survivor_link_x<-as.mcmc.list(SARI_link_output_delta[[4]])
SARI_combined_Survivor_link_x <- as_tibble(do.call('rbind',SARI_Survivor_link_x))
delta_slope_pred<-dplyr::select(SARI_combined_Survivor_link_x,starts_with('delta'))%>%as.matrix()%>%array(dim=c(n_sim,S))%>%
  apply(c(2),quantile,c(0.025,0.5,0.975))%>%melt(varnames=c('quantile','s'),value.name='y')%>%
  spread(quantile,y)%>%mutate(Nowcast_date=as.numeric(paste(Nowcast_list[4])),s=federal_units[s])
delta_matrix_pred<-matrix(unlist(delta_slope_pred), byrow = FALSE, ncol = 5)%>%as_tibble()
colnames(delta_matrix_pred)<-c('s_abb','2.5%','50%','97.5%','Nowcast')
delta_matrix_pred$s_abb<-factor(delta_matrix_pred$s_abb,levels=delta_matrix_pred$s_abb)
delta_matrix_pred<-full_join(delta_matrix_pred,state_guide[,2:3], by='s_abb')
delta_matrix_pred$s_full<-factor(delta_matrix_pred$s_full, levels=delta_matrix_pred$s_full)
delta_matrix_pred<-delta_matrix_pred%>%mutate(model='Prediction accuracy experiment')
delta_matrix_pred[,2]<-as.numeric(unlist(delta_matrix_pred[,2]))
delta_matrix_pred[,3]<-as.numeric(unlist(delta_matrix_pred[,3]))
delta_matrix_pred[,4]<-as.numeric(unlist(delta_matrix_pred[,4]))

delta_matrix_all<-full_join(delta_matrix_pred,delta_matrix, by=c('s_abb','s_full','model','50%','97.5%','2.5%'))

delta_comp_plot<-ggplot()+geom_hline(yintercept = 0,colour='darkgrey',linetype='dashed')+
  geom_errorbar(data=delta_matrix_all,
                aes(x=s_abb,  ymin=`2.5%`,ymax=`97.5%`,colour=model))+
  geom_point(data=delta_matrix_all,aes(x=s_abb, y=`50%`,colour=model),shape=3,size=5)+
  scale_colour_viridis_d(name='Model',begin=0,end=0.7)+
  scale_fill_viridis_d(name='Model',begin=0,end=0.7)+
  labs(x="Brazilian federative unit (s)",y=expression(delta['s']),
       title='SARI hospitalisations',
       subtitle='Linear effect of case load on reporting delay', legend=NULL)+
  theme_minimal()+
  theme(text = element_text(size=12), plot.title = element_text(size=16),plot.subtitle=element_text(size=14),
        legend.title = element_text(size=14),strip.text = element_text(size=14), 
        legend.position = "bottom")
delta_comp_plot
ggsave(delta_comp_plot,file='full_delta_sari_compare.pdf',width=9,height=4)


