
# EXPLORATORY PLOTS FOR SARI HOSPITALISATIONS IN BRAZIL #

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
library(xml2)
library(rvest)

#setwd("~/OneDrive/PhD/comp 2023 paper code")
# Define beta binomial distribution._linear_log
rbetabin=nimbleFunction(run=function(n=integer(0),mu=double(0),phi=double(0),size=double(0)){
  pi=rbeta(n,mu*phi,(1-mu)*phi)
  returnType(double(1))
  return(rbinom(n,size,pi))
})

#Read in data..
INFLUD21 <- read_delim("INFLUD21-30-01-2023.csv", delim = ";", escape_double = FALSE, trim_ws = TRUE)

INFLUD22 <- read_delim("INFLUD22-30-01-2023.csv", 
                       delim = ";", escape_double = FALSE, 
                       col_types = cols(DT_NOTIFIC = col_date(format = "%d/%m/%Y"), 
                                        DT_SIN_PRI = col_date(format = "%d/%m/%Y")), 
                       trim_ws = TRUE)

col_ant<-which(colnames(INFLUD21)=="AN_SARS2")
col_pcr<-which(colnames(INFLUD21)=="PCR_SARS2")
col_results<-which(colnames(INFLUD21)=="PCR_RESUL")
col_an_results<-which(colnames(INFLUD21)=="RES_AN")
col_hospital<-which(colnames(INFLUD21)=="DT_INTERNA")
col_collection<-which(colnames(INFLUD21)=="DT_COLETA")
col_an_result<-which(colnames(INFLUD21)=="DT_RES_AN")
col_pcr_result<-which(colnames(INFLUD21)=="DT_PCR")
col_notif_date<-which(colnames(INFLUD21)=="DT_DIGITA")

SARI_2021<- INFLUD21[,c(1:16, col_hospital, col_collection, col_results,col_an_results, col_pcr_result, col_notif_date, col_pcr, col_an_result, col_ant)]%>%
  mutate(DT_NOTIFIC=as.Date((DT_NOTIFIC),format="%d/%m/%Y"),DT_DIGITA=as.Date((DT_DIGITA),format="%d/%m/%Y"),DT_SIN_PRI=as.Date((DT_SIN_PRI),format="%d/%m/%Y"))
dates_2021<-sort(as.Date(unique(SARI_2021$DT_SIN_PRI),format="%d/%m/%Y"))

SARI_2022<- INFLUD22[,c(1:16, col_hospital, col_collection, col_results, col_an_results,col_pcr_result, col_notif_date, col_pcr, col_an_result, col_ant)]%>%
  mutate(DT_NOTIFIC=as.Date((DT_NOTIFIC),format="%d/%m/%Y"),DT_DIGITA=as.Date((DT_DIGITA),format="%d/%m/%Y"),DT_SIN_PRI=as.Date((DT_SIN_PRI),format="%d/%m/%Y"))
dates_2022<-sort(as.Date(unique(SARI_2022$DT_SIN_PRI),format="%d/%m/%Y"))

SARI_2122<-rbind(SARI_2021,mutate(SARI_2022, SEM_PRI=SEM_PRI+52,SEM_NOT=SEM_NOT+52))%>%
  mutate(delay=difftime(DT_DIGITA,DT_SIN_PRI, unit=c("days")),delay_noti=difftime(DT_NOTIFIC,DT_SIN_PRI, unit=c("days")))%>%
  mutate(DT_COLETA=as.Date((DT_COLETA),format="%d/%m/%Y"),DT_PCR=as.Date((DT_PCR),format="%d/%m/%Y"))%>%
  mutate(pcr_delay=difftime(DT_PCR,DT_COLETA, unit=c("days")))%>%mutate(delay_diffweeks=as.numeric(difftime(DT_DIGITA,DT_SIN_PRI, unit=c("weeks"))))


#DT_DIGITA
sum(SARI_2122$delay<0)
sum(SARI_2122$delay_diffweeks<0)

#DT_NOTIFIC
sum(SARI_2122$delay_noti<0)

# delay of zero happened in first week and 1 happened in second week
SARI_2122$delay_week<-ceiling(SARI_2122$delay_diffweeks+1/7)

# Make matrix with onset week as rows and delay weeks as columns
SARI_delay_long<-filter(SARI_2122, delay_week>0)%>%group_by(SEM_PRI,delay_week,SG_UF_NOT)%>%summarise(cases=n())

SARI_delay<-SARI_delay_long%>%ungroup()%>%
  pivot_wider(id_cols=c(SEM_PRI,SG_UF_NOT),names_from = delay_week,values_from = cases)
SARI_delay[is.na(SARI_delay)]<-0

SARI_delay_totals<-SARI_delay
SARI_delay_totals$total<-rowSums(SARI_delay[,-(1:2)], na.rm=T)

N_raw<-max(SARI_delay$SEM_PRI)
# remove latest week as only 3 regions have observations 
SARI_delay<-filter(SARI_delay, SEM_PRI<N_raw)
N<-max(SARI_delay$SEM_PRI)
D_N<-dim(SARI_delay)[2]-2

SARI_delay_ordered<-as.matrix(SARI_delay)
SARI_delay_ordered<-SARI_delay_ordered[order(SARI_delay_ordered[,2]),]
SARI_delay_ordered<-SARI_delay_ordered[order(SARI_delay_ordered[,1]),]

# Save spatial constants
federal_units<-sort(unique(SARI_delay_ordered[,2]))
S<-length(federal_units)

# Make an array including delay:
# SARI_array[REGION, TIME, DELAY]
SARI_array<-SARI_delay_ordered[,-(1:2)]%>%as.matrix()%>%as.numeric()%>%array(dim=c(S,N,D_N))



# total SARI cases for each region
total_SARI<-apply(SARI_array,c(1,2),sum)
SARI_delay_totals<-melt(total_SARI)
colnames(SARI_delay_totals)<-c('SG_UF_NOT',"SARI","totals")


# set D_max
D_max<-20
D<-8

SARI_z<-SARI_array[,,1:(D_max)]
SARI_y<-apply(SARI_z,c(1,2),sum)


# DATA frame for COVID cases

SARI_2021<- INFLUD21[,c(1:16,col_results,col_an_results,col_pcr, col_ant)]
dates_2021<-sort(as.Date(unique(SARI_2021$DT_SIN_PRI),format="%d/%m/%Y"))

SARI_counts_21<-SARI_2021%>%group_by(SEM_PRI,SG_UF_NOT)%>%
  summarise(SARI=n(), PCR=sum(PCR_SARS2 == 1,na.rm = TRUE ),
            AN=sum(AN_SARS2 == 1,na.rm = TRUE ),
            SARS2=sum((PCR_SARS2==1)|(AN_SARS2==1),na.rm=TRUE),
            notTESTS=sum(((is.na(RES_AN))|(RES_AN==4))&((is.na(PCR_RESUL))|(PCR_RESUL==4))))%>%
  mutate(nu=SARS2/SARI)
SARI_counts_21<-SARI_counts_21%>%mutate( WEEK=SEM_PRI)

SARI_counts_total_21<-SARI_2021%>%group_by(SEM_PRI)%>%
  summarise(SARI=n(), PCR=sum(PCR_SARS2 == 1,na.rm = TRUE ), 
            AN=sum(AN_SARS2 == 1,na.rm = TRUE),
            SARS2=sum((PCR_SARS2==1)|(AN_SARS2==1),na.rm = TRUE),
            notTESTS=sum(((is.na(RES_AN))|(RES_AN==4))&((is.na(PCR_RESUL))|(PCR_RESUL==4))))%>%
  mutate(nu=SARS2/SARI)
SARI_counts_total_21<-SARI_counts_total_21%>%mutate( WEEK=SEM_PRI)

SARI_2022<- INFLUD22[,c(1:16,col_results,col_an_results,col_pcr, col_ant)]
dates_2022<-sort(as.Date(unique(SARI_2022$DT_SIN_PRI),format="%d/%m/%Y"))

SARI_counts_22<-SARI_2022%>%group_by(SEM_PRI,SG_UF_NOT)%>%summarise(SARI=n(), PCR=sum(PCR_SARS2 == 1,na.rm = TRUE ),
                                                                    AN=sum(AN_SARS2 == 1,na.rm = TRUE ),
                                                                    SARS2=sum((PCR_SARS2==1)|(AN_SARS2==1),na.rm = TRUE),
                                                                    notTESTS=sum(((is.na(RES_AN))|(RES_AN==4))&((is.na(PCR_RESUL))|(PCR_RESUL==4))))%>%
  mutate(nu=SARS2/SARI)
SARI_counts_22<-SARI_counts_22%>%mutate( WEEK=SEM_PRI+52)

SARI_counts<-rbind(SARI_counts_21,SARI_counts_22)
first_date<-min(dates_2021)
dates_all<-c(dates_2021,dates_2022)
ggplot(SARI_counts)+
  geom_line(aes(x=as.Date(first_date)+(WEEK-1)*7,y=nu))+
  labs(x='2021/2022', y='proportion of COVID-19 cases')+
  facet_wrap(~SG_UF_NOT)+geom_hline(yintercept = 1)

ggplot(SARI_counts)+
  geom_line(aes(x=as.Date(first_date)+(WEEK-1)*7,y=SARI),colour='red')+
  geom_line(aes(x=as.Date(first_date)+(WEEK-1)*7,y=SARS2),colour='blue')+
  labs(x='2021/2022', y='proportion of COVID-19 cases')+
  facet_wrap(~SG_UF_NOT,scale='free')


ggplot(SARI_counts)+
  geom_line(aes(x=SARI,y=SARS2/SARI),alpha=0.2,colour='red')+
  geom_point(aes(x=SARI,y=SARS2/SARI))+
  labs(x='SARI cases', y='Proportion of COVID-19 cases')+
  facet_wrap(~SG_UF_NOT,scale='free')

ggplot(SARI_counts)+
  geom_line(aes(x=SARI,y=notTESTS/SARI),alpha=0.2,colour='purple')+
  geom_point(aes(x=SARI,y=notTESTS/SARI))+
  labs(x='SARI cases', y='proportion NOT TESTED')+
  facet_wrap(~SG_UF_NOT,scale='free')


## EXPLORE DATA BY AGE GROUP
SARI_2021<- INFLUD21[,c(1:16, col_hospital, col_collection, col_results, col_pcr_result, col_notif_date, col_pcr, col_an_result, col_ant)]%>%
  mutate(DT_NOTIFIC=as.Date((DT_NOTIFIC),format="%d/%m/%Y"),DT_DIGITA=as.Date((DT_DIGITA),format="%d/%m/%Y"),DT_SIN_PRI=as.Date((DT_SIN_PRI),format="%d/%m/%Y"))
dates_2021<-sort(as.Date(unique(SARI_2021$DT_SIN_PRI),format="%d/%m/%Y"))

SARI_2022<- INFLUD22[,c(1:16, col_hospital, col_collection, col_results, col_pcr_result, col_notif_date, col_pcr, col_an_result, col_ant)]%>%
  mutate(DT_NOTIFIC=as.Date((DT_NOTIFIC),format="%d/%m/%Y"),DT_DIGITA=as.Date((DT_DIGITA),format="%d/%m/%Y"),DT_SIN_PRI=as.Date((DT_SIN_PRI),format="%d/%m/%Y"))%>%
  filter(DT_SIN_PRI<as.Date("2023-01-01")) # cut off at end of 2022 to enhaze fully reported data
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
SARI_link_z<-SARI_array[,,1:D_max,] #[REGION, TIME, DELAY, AGE GROUP]
# Matrix for total SARI cases.
SARI_link_y<-apply(SARI_link_z,c(1,2,4), sum) #[REGION, TIME, AGE GROUP]
# Dates of all weeks in data frame. 
SARI_dates<-as.Date(first_date)+(1:ncol(SARI_link_y)-1)*7

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
SARI_link_x<-COVID_wide # [REGION, TIME, AGE]

SARI_x<-t(apply(SARI_link_x,c(1,2),sum))
SARI_z<-apply(SARI_link_z,c(1,2,3),sum)
SARI_y<-t(apply(SARI_link_y,c(1,2),sum))

# Melt arrays/matrix to create data frames:
# Data frame for COVID cases.
x_all<-melt(apply(SARI_link_x,c(1,2),sum))
colnames(x_all)<-c("s","t","x")
x_all$s<-federal_units[x_all$s]
# Data frame for SARI cases.
y_all<-melt(apply(SARI_link_y,c(1,2),sum))
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

# Plot of total SARI & COVID-19 hospitalisations
# fit a gam over space and time:
all_data<-full_join(full_join(x_all,y_all,by=c('t','s')),population, by='s')
covid_lines<-numeric()
for(i in 1:S){
covid_gam<-gam(data = filter(all_data,s==federal_units[i]), x~s(t,k=20), family = poisson(link = 'log'))
covid_lines<-rbind(covid_lines,left_join(tibble(x=(covid_gam$fitted.values),t=1:103)%>%mutate(s=federal_units[i]),population, by='s'))
}
non_covid_lines<-numeric()
for(i in 1:S){
  non_covid_gam<-gam(data = filter(all_data,s==federal_units[i])%>%mutate(x2=y-x), x2~s(t,k=20), family = poisson(link = 'log'))
  non_covid_lines<-rbind(non_covid_lines,left_join(tibble(x=(non_covid_gam$fitted.values),t=1:103)%>%mutate(s=federal_units[i]),population, by='s'))
}
totals_plot<-ggplot(data=all_data)+
  geom_point(aes(x=as.Date(first_date)+t*7,y=(y-x),colour='non-COVID-19 SARI'),alpha=0.5)+ 
  geom_point(aes(x=as.Date(first_date)+t*7,y=x,colour='COVID-19-positive SARI'),alpha=0.5)+ 
 # geom_point(aes(x=as.Date(first_date)+t*7,y=y,colour='SARI'),alpha=0.5)+ 
  geom_line(data=non_covid_lines,aes(x=as.Date(first_date)+t*7,y=x,colour='non-COVID-19 SARI'))+ 
  geom_line(data=covid_lines, aes(x=as.Date(first_date)+t*7,y=x,colour='COVID-19-positive SARI'))+ 
  facet_wrap(~s_full, scales = 'free', nrow = 9)+
  scale_x_date(breaks = "4 month", date_labels =  "%b %y")+
  theme_minimal()+labs(y="Hospitalisations", x="Date", title="Brazil hospitalisations by federative unit")+#+geom_vline(xintercept = 51)
  scale_colour_discrete(name='Virus')+
    theme(text = element_text(size=14), axis.title = element_text(size=16), plot.title = element_text(size=20),plot.subtitle=element_text(size=18),
          legend.title = element_text(size=16),strip.text = element_text(size=14),
          axis.text = element_text(size=9),legend.position="bottom")
totals_plot
ggsave(totals_plot,file='SARI_totals_plot.pdf',width=10,height=12)



sari_wide<-melt(SARI_z)
colnames(sari_wide)<-c("s","t","d","z")
Data_wide<-sari_wide%>%ungroup()%>%
  pivot_wider(id_cols=c(t,s),names_from = d,values_from = z)%>%as.matrix()
Data_wide<-Data_wide[order(Data_wide[,2],Data_wide[,1]),]

Data_delay_gam<-as.tibble(Data_wide)
totals<-Data_delay_gam%>%mutate(total=as.numeric(unlist((apply(Data_wide[,3:22],1, sum, na.rm=TRUE)))))%>%mutate(s=paste(federal_units[s]))
Data_delay_gam$prop_0<-as.numeric(unlist((Data_wide[,3])))
Data_delay_gam$prop_1<-as.numeric(unlist((Data_wide[,3]+Data_wide[,4])))
Data_delay_gam$prop_2<-as.numeric(unlist((Data_wide[,3]+Data_wide[,4]+Data_wide[,5])))
Data_delay_gam$prop_3<-as.numeric(unlist((Data_wide[,3]+Data_wide[,4]+Data_wide[,5]+Data_wide[,6])))
Data_delay_gam$prop_4<-as.numeric(unlist((Data_wide[,3]+Data_wide[,4]+Data_wide[,5]+Data_wide[,6]+Data_wide[,7])))

Data_delay_long_gam<-as.tibble(Data_delay_gam%>%pivot_longer(!c(t,s), names_to = "delay", values_to = "count")%>%mutate(s=paste(federal_units[s])))
Data_delay_pop_gam<-full_join(Data_delay_long_gam,population, by="s")
Data_delay_prop_gam<-Data_delay_pop_gam%>%filter(delay%in%c("prop_0","prop_1","prop_2","prop_3","prop_4"))
Data_delay_prop_gam<-full_join(Data_delay_prop_gam,totals[,c(1,2,dim(totals)[2])],by = join_by(t, s))

Data_delay_prop_gam_brazil<-Data_delay_prop_gam%>%group_by(t,delay,cum_prop)%>%summarise(total=sum(total),population=sum(population), count=sum(count))

Data_delay_prop_gam<-Data_delay_prop_gam%>%mutate(total.per.capita=total/population, cum_prop=rep(0:4,N*S))

# scale<-filter(Data_delay_prop,delay=='prop_0')%>%group_by(s)%>%summarise(scale=scale(total))
# scale<-scale%>%mutate(t=rep(1:N))
# colnames(scale)[2]<-c("scale")
# Data_delay_prop2<-left_join(Data_delay_prop,scale, by=c("s", "t"))
# 
# logtotal_regions<-gam(data = Data_delay_prop2, (count)~log(total.per.capita)*s+delay, family = gaussian(link = log))
# scale_regions<-gam(data = Data_delay_prop2, (count)~scale*s+delay, family = gaussian(link = log))

# plot gam of relationship 
Data_delay_prop_gam<-mutate(Data_delay_prop_gam,delay=as.factor(delay),s=as.factor(s))
Data_delay_prop_gam_brazil<-mutate(Data_delay_prop_gam_brazil,delay=as.factor(delay))
cumu_gam_all<-gam(data = Data_delay_prop_gam_brazil,probit(count/total)~s(log(total),k=10,by=delay)+delay, family = gaussian)
cumu_gam_all_linear<-gam(data = Data_delay_prop_gam_brazil,probit(count/total)~log(total)*delay, family = gaussian)
#plot(cumu_gam_all)

overall_smooth_trend<-tibble(y=cumu_gam_all$fitted.values,totals=Data_delay_prop_gam_brazil$total,cum_prop=Data_delay_prop_gam_brazil$cum_prop)
overall_linear_trend<-tibble(y=cumu_gam_all_linear$fitted.values,totals=Data_delay_prop_gam_brazil$total,cum_prop=Data_delay_prop_gam_brazil$cum_prop)

pred <- predict(cumu_gam_all, se.fit = TRUE)
Data_delay_prop_gam_brazil$fit <- pred$fit
Data_delay_prop_gam_brazil$upper <- pred$fit + 1.96 * pred$se.fit
Data_delay_prop_gam_brazil$lower <- pred$fit - 1.96 * pred$se.fit

ggplot(Data_delay_prop_gam_brazil, aes(x=log(total), y=probit(count/total), colour=as.factor(cum_prop))) +
  geom_point(alpha=0.5) +
  geom_line(aes(y=fit_time), linetype='dashed') +
  geom_ribbon(aes(ymin=lower_time, ymax=upper_time, fill=as.factor(cum_prop)), alpha=0.2) +
  scale_colour_viridis_d(name="Cumulative proportions reported at delay week:", begin=0.1, end=0.85) +
  scale_fill_viridis_d(name="Cumulative proportions reported at delay week:", begin=0.1, end=0.85) +
  labs(x="Eventual total SARI hospitalisations reported (log transform)", 
       y="Cumulative proportions reported (probit transform)", 
       title="SARI hospitalisations for Brazil (Adjusted for Time)") +
  theme_minimal()


ggplot(Data_delay_prop_gam_brazil, aes(x=log(total), y=probit(count/total), colour=as.factor(cum_prop))) +
  geom_point(alpha=0.5) +
  geom_line(aes(y=fit), linetype='dashed') +
  geom_ribbon(aes(ymin=lower, ymax=upper, fill=as.factor(cum_prop)), alpha=0.2) +
  scale_colour_viridis_d(name="Cumulative proportions reported at delay week:", begin=0.1, end=0.85) +
  scale_fill_viridis_d(name="Cumulative proportions reported at delay week:", begin=0.1, end=0.85) +
  labs(x="Eventual total SARI hospitalisations reported (log transform)", 
       y="Cumulative proportions reported (probit transform)", 
       title="SARI hospitalisations for Brazil") +
  theme_minimal()

brazil_cumulative_totals<-ggplot()+
  geom_point(data=filter(Data_delay_prop_gam_brazil,cum_prop%in%c(0:4)),aes(x=log(total),y=probit(count/total), colour=as.factor(cum_prop)),alpha=0.5)+#facet_wrap(~s_full,scales="free", ncol=3)+
  geom_line(data=overall_smooth_trend,aes(x=log(totals),y=y, colour=as.factor(cum_prop)),linetype='dashed')+
  geom_smooth(data=filter(Data_delay_prop_gam_brazil, cum_prop %in% c(0:4)), 
              aes(x=log(total), y=probit(count/total), colour=as.factor(cum_prop)), 
              level=0.95, method='gam', formula=y~s(x, k=10))+
   #geom_smooth(data=filter(Data_delay_prop_gam_brazil,cum_prop%in%c(0:4)),aes(x=log(total),y=probit(count/total), colour=as.factor(cum_prop)),level=0.95,method='lm')+
  #geom_line(data=overall_linear_trend,aes(x=log(totals),y=y, colour=as.factor(cum_prop)))+
  scale_colour_viridis_d(name="Cumulative proportions reported at delay week:", begin=0.1, end=0.85)+
  labs(x="Eventual total SARI hospitaltions reported (log transform)",y="Cumulative proportions reported (probit transform)",
       title = 'SARI hospitalisations for Brazil')+theme_minimal()+
  theme(text = element_text(size=14), axis.title = element_text(size=16), plot.title = element_text(size=20),plot.subtitle=element_text(size=18),
        legend.title = element_text(size=16),strip.text = element_text(size=14),
        axis.text = element_text(size=9),legend.position="bottom")
brazil_cumulative_totals
ggsave(brazil_cumulative_totals,file='brazil_cumulative_totals.pdf',width=10,height=6)

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
Data_delay_prop<-Data_delay_prop%>%mutate(total.per.capita=total/population, cum_prop=rep(0:4,N*S))




# plot cumulative proportions against totals: smooths
sari_cumulative_gams<-ggplot(filter(Data_delay_prop,cum_prop%in%c(0,2,4)),aes(x=log(total),y=(count), colour=as.factor(cum_prop)))+
  geom_line(alpha=0.2,aes(group=s))+facet_wrap(~s_full,scales="free", ncol=3)+
  geom_smooth(alpha=0.1)+
  scale_colour_viridis_d(name="Cumulative proportions reported at delay week:", begin=0.1, end=0.85)+
  labs(x="Eventual total SARI hospitaltions reported (log scale)",y="Cumulative proportions of SARI hospitalisations reported",
       title = 'SARI hospitalisations cumulative proportions reported')+theme_minimal()+
  theme(legend.position = 'bottom')+ 
  theme(text = element_text(size=14), axis.title = element_text(size=16), plot.title = element_text(size=20),plot.subtitle=element_text(size=18),
        legend.title = element_text(size=16),strip.text = element_text(size=14),
        axis.text = element_text(size=9),legend.position="bottom")
sari_cumulative_gams
ggsave(sari_cumulative_gams,file='sari_cumulative_gams.pdf',width=9,height=12)



sari_cumulative_exploratory_all<-ggplot(Data_delay_prop,aes(x=total,y=(count), colour=as.factor(cum_prop)))+
  geom_point(alpha=0.25)+facet_wrap(~s_full,scales="free", ncol=3)+geom_smooth(alpha=0.1)+scale_x_log10()+
  scale_colour_viridis_d(name="Cumulative proportions reported at delay week:", begin=0.1, end=0.85)+
  labs(x="Eventual total SARI hospitaltions reported (log scale)",y="Cumulative proportions of SARI hospitalisations reported",
       title = 'SARI hospitalisations cumulative proportions reported')+theme_minimal()+
  theme(legend.position = 'bottom')+
  theme(text = element_text(size=14), axis.title = element_text(size=16), plot.title = element_text(size=20),plot.subtitle=element_text(size=18),
        legend.title = element_text(size=16),strip.text = element_text(size=14),
        axis.text = element_text(size=9),legend.position="bottom")
ggsave(sari_cumulative_exploratory_all,file='sari_cumulative_exploratory_all.pdf',width=10,height=16.5)


mean.sari=floor((apply(SARI_y[1:(N-D_max),],2,mean)))
sari_cumulative_exploratory_3<-ggplot(filter(Data_delay_prop,cum_prop%in%c(0,2,4)),aes(x=log(total.per.capita),y=(count), colour=as.factor(cum_prop)))+
  geom_point(alpha=0.25)+facet_wrap(~s_full,scales="free", ncol=3)+
  geom_smooth(alpha=0.1)+
  scale_colour_viridis_d(name="Cumulative proportions reported at delay week:", begin=0.1, end=0.85)+
  labs(x="Eventual total SARI hospitaltions reported per capita (log scale)",y="Cumulative proportions of SARI hospitalisations reported",
       title = 'SARI hospitalisations cumulative proportions reported')+theme_minimal()+
  theme(legend.position = 'bottom')+ 
  theme(text = element_text(size=14), axis.title = element_text(size=16), plot.title = element_text(size=20),plot.subtitle=element_text(size=18),
        legend.title = element_text(size=16),strip.text = element_text(size=14),
        axis.text = element_text(size=9),legend.position="bottom")
ggsave(sari_cumulative_exploratory_3,file='sari_cumulative_exploratory_3.pdf',width=9,height=12)

# plot cumulative proportions against time
cumu_lines<-numeric()
for(i in 1:S){
  for(d in 1:5){
  cumu_gam<-gam(data = filter(Data_delay_prop,s==federal_units[i],cum_prop==(d-1)), count~s(t,k=20), family = gaussian(link = log))
  cumu_lines<-rbind(cumu_lines,left_join(tibble(count=(cumu_gam$fitted.values),t=1:103)%>%mutate(s=federal_units[i], cum_prop=d-1),population, by='s'))
  }
  # cumu_lines<-rbind(cumu_lines,left_join(matrix((cumu_gam$fitted.values),ncol=103)%>%melt(varnames=c('cum_prop','t'),value.name="count")%>%mutate(s=federal_units[i],cum_prop=cum_prop-1),population, by='s'))
}

# Data_delay_prop<-mutate(Data_delay_prop,delay=as.factor(delay),s=as.factor(s))
# cumu_gam_all<-gam(data = Data_delay_prop, count~s(t,k=10,by=s), family = gaussian(link = log))


cumulative_plot_time<-ggplot(Data_delay_prop,aes(x=as.Date(first_date)+t*7,y=(count), colour=as.factor(cum_prop)))+
  geom_point(alpha=0.5)+facet_wrap(~s_full,scales="free", ncol=3)+
  geom_line(data=cumu_lines, aes(x=as.Date(first_date)+t*7, y=probit(count)))+
  scale_colour_viridis_d(name="Cumulative proportions reported at delay week:", begin=0.1, end=0.85)+
  labs(x="Date",y="Cumulative proportions of hospitalisations reported",
       title = 'Cumulative proportions of SARI hospitalisations reported by federative unit')+theme_minimal()+
  theme(legend.position = 'bottom')+
  scale_x_date(breaks = "4 month", date_labels =  "%b %y")+
  theme(text = element_text(size=14), axis.title = element_text(size=16), plot.title = element_text(size=20),plot.subtitle=element_text(size=18),
        legend.title = element_text(size=16),strip.text = element_text(size=14),
        axis.text = element_text(size=9),legend.position="bottom")
ggsave(cumulative_plot_time,file='cumulative_plot_time.pdf',width=10,height=12)

SARI_counts_full<-SARI_counts%>%mutate(s=SG_UF_NOT)%>%full_join(population, by='s')
testing_rates_plots<-ggplot(SARI_counts_full)+
  geom_smooth(aes(x=log(SARI),y=(SARI-notTESTS)/SARI),alpha=0.5,colour="#4363D8")+
  geom_point(aes(x=log(SARI),y=(SARI-notTESTS)/SARI),alpha=0.4)+
  labs(x='Total SARI hospitaltions reported (log scale)', y='Proportion of SARI hospitalisations tested for COVID-19',title = 'SARI hospitalisations testing rates for COVID-19')+
  facet_wrap(~s_full,scale='free',ncol=3)+theme_minimal()+
  theme(text = element_text(size=14), axis.title = element_text(size=16), plot.title = element_text(size=20),plot.subtitle=element_text(size=18),
        legend.title = element_text(size=16),strip.text = element_text(size=14),
        axis.text = element_text(size=9),legend.position="bottom")
ggsave(testing_rates_plots,file='testing_rates_plots.pdf',width=10,height=16)

Data_delay_prop_total<-Data_delay_prop%>%mutate(cum_counts=count*total)%>%group_by(t,cum_prop)%>%
  summarise(total=sum(total),cum_counts=sum(cum_counts) )%>%mutate(count=cum_counts/total)

plot1<-ggplot(filter(Data_delay_prop_total,cum_prop%in%c(0,2,4)),aes(x=as.Date(first_date)+t*7,y=total))+
  geom_line()+#facet_wrap(~s_full,scales="free", ncol=3)+
  #geom_smooth(alpha=0.1)+scale_y_log10()+
  labs(x="Date of hospitalisation",y="SARI hospitalisations",
       subtitle = 'Total hospitalisations ',title='SARI hospitalisations for Brazil')+theme_minimal()+
  theme(legend.position = 'bottom')+
  scale_x_date(date_breaks = "2 month", date_labels =  "%b %Y")+
  theme(text = element_text(size=14), axis.title = element_text(size=16), plot.title = element_text(size=20),plot.subtitle=element_text(size=18),
        legend.title = element_text(size=16),strip.text = element_text(size=14),
        axis.text = element_text(size=9),legend.position="bottom")

plot2<-ggplot(filter(Data_delay_prop_total,cum_prop%in%c(0,2,4)),aes(x=as.Date(first_date)+t*7,y=(count), colour=as.factor(cum_prop)))+
  geom_point(alpha=0.5)+#facet_wrap(~s_full,scales="free", ncol=3)+
  geom_smooth( method='lm',level=0.95)+#scale_x_log10()+
  scale_colour_viridis_d(name="Cumulative proportions reported at delay week:", begin=0.1, end=0.85)+
  labs(x="Date of hospitalisation",y="Cumulative proportions",
       subtitle = 'Cumulative proportions of SARI hospitalisations reported ')+theme_minimal()+
  theme(legend.position = 'bottom')+
  scale_x_date(date_breaks = "2 month", date_labels =  "%b %Y")+
  theme(text = element_text(size=14), axis.title = element_text(size=16), plot.title = element_text(size=20),plot.subtitle=element_text(size=18),
        legend.title = element_text(size=16),strip.text = element_text(size=14),
        axis.text = element_text(size=9),legend.position="bottom")

sari_time_plots<-grid.arrange(plot1,plot2,ncol=1)
ggsave(sari_time_plots,file='sari_time_plots.pdf',width=9,height=8)

# check to see if delay in cumulative proportions are related to totals:

stripe_indicator=function(x=double(1)){
  stripe<-x%in%c(1:20,40:60,80:100,120:140,160:180)
  return(stripe)
}

  

ggplot()+
  geom_point(data=filter(Data_delay_prop_total,cum_prop%in%c(0)),aes(x=as.Date(first_date)+t*7,y=log(total)/10),alpha=0.9)+#facet_wrap(~s_full,scales="free", ncol=3)+
  geom_line(data=filter(Data_delay_prop_total,cum_prop%in%c(0)),aes(x=as.Date(first_date)+t*7,y=probit(count)*(-1), colour=as.factor(cum_prop)),alpha=0.9)+
#  geom_line(data=filter(Data_delay_prop_total,cum_prop%in%c(0)),aes(x=as.Date(first_date)+t*7,y=probit(count)*(-1)-(log(total)/10), colour='difference'),alpha=0.9)+
  #facet_wrap(~s_full,scales="free", ncol=3)+
 # geom_smooth(data=filter(Data_delay_prop_total,cum_prop%in%c(0)),aes(x=as.Date(first_date)+t*7,y=probit(count)*(-1), colour=as.factor(cum_prop)),alpha=0.1)+#scale_y_log10()+
  scale_x_date(breaks = dates_all[seq(from=1, to=N,length=10)])+
  labs(x="Date of hospitalisation",y="Total SARI hospitalisations",
       subtitle = 'SARI hospitalisations overtime')+theme_bw()+
  theme(legend.position = 'bottom')+
  theme(text = element_text(size=12), axis.title = element_text(size=14), plot.title = element_text(size=16),plot.subtitle=element_text(size=14),
        legend.title = element_text(size=14),strip.text = element_text(size=12),legend.position="bottom")


# with month lag
ggplot()+
  geom_point(data=filter(Data_delay_prop_total,cum_prop%in%c(0)),aes(x=as.Date(first_date)+t*7-4*7,y=log(total)/10),alpha=0.9)+#facet_wrap(~s_full,scales="free", ncol=3)+
  geom_line(data=filter(Data_delay_prop_total,cum_prop%in%c(0)),aes(x=as.Date(first_date)+t*7,y=probit(count)*(-1), colour=as.factor(cum_prop)),alpha=0.9)+
  #facet_wrap(~s_full,scales="free", ncol=3)+
  #geom_line(data=filter(Data_delay_prop_total,cum_prop%in%c(0)),aes(x=as.Date(first_date)+t*7,y=probit(count)*(-1)-(log(total)/10), colour='difference'),alpha=0.9)+
 # geom_smooth(data=filter(Data_delay_prop_total,cum_prop%in%c(0)),aes(x=as.Date(first_date)+t*7,y=probit(count)*(-1), colour=as.factor(cum_prop)),alpha=0.1)+#scale_y_log10()+
  scale_x_date(breaks = dates_all[seq(from=1, to=N,length=10)])+
  labs(x="Date of hospitalisation",y="Total SARI hospitalisations",
       subtitle = 'SARI hospitalisations overtime lagged by a month')+theme_bw()+
  theme(legend.position = 'bottom')+
  theme(text = element_text(size=12), axis.title = element_text(size=14), plot.title = element_text(size=16),plot.subtitle=element_text(size=14),
        legend.title = element_text(size=14),strip.text = element_text(size=12),legend.position="bottom")

gradient<-numeric()
for(i in 2:N){
  gradient[i-1]<-filter(Data_delay_prop_total,cum_prop==0)$total[i]-filter(Data_delay_prop_total,cum_prop==0)$total[i-1]
}

#plot gradient of total counts:
ggplot()+
  geom_point(data=filter(Data_delay_prop_total,cum_prop%in%c(0)),aes(x=as.Date(first_date)+t*7,y=log(total)/10),alpha=0.9)+#facet_wrap(~s_full,scales="free", ncol=3)+
  geom_line(data=filter(Data_delay_prop_total,cum_prop%in%c(0)),aes(x=as.Date(first_date)+t*7,y=probit(count)*(-1), colour=as.factor(cum_prop)),alpha=0.9)+
  #facet_wrap(~s_full,scales="free", ncol=3)+
  geom_smooth(aes(x=dates_all[1:(N-1)],y=gradient/10000,colour='gradient'))+
  geom_line(aes(x=dates_all[1:(N-1)],y=gradient/10000,colour='gradient'))+
  geom_smooth(data=filter(Data_delay_prop_total,cum_prop%in%c(0)),aes(x=as.Date(first_date)+t*7,y=probit(count)*(-1), colour=as.factor(cum_prop)),alpha=0.1)+#scale_y_log10()+
  scale_x_date(breaks = dates_all[seq(from=1, to=N,length=10)])+
  labs(x="Date of hospitalisation",y="Total SARI hospitalisations",
       subtitle = 'SARI hospitalisations overtime')+theme_bw()+
  theme(text = element_text(size=14), axis.title = element_text(size=16), plot.title = element_text(size=20),plot.subtitle=element_text(size=18),
        legend.title = element_text(size=16),strip.text = element_text(size=14),
        axis.text = element_text(size=9),legend.position="bottom")
cum_plot_time<-ggplot(filter(Data_delay_prop_total),aes(x=as.Date(first_date)+t*7,y=probit(count), colour=as.factor(cum_prop)))+
  geom_point(alpha=0.5)+#facet_wrap(~s_full,scales="free", ncol=3)+
  geom_smooth(alpha=0.1)+#scale_x_log10()+
  scale_colour_viridis_d(name="Cumulative proportions reported at delay week:", begin=0.1, end=0.85)+
  labs(x="Date of hospitalisation",y="Cumulative proportions (probit scale)",
       title = 'Cumulative proportions of SARI hospitalisations reported for Brazil',subtitle='Reported ')+theme_minimal()+
  theme(legend.position = 'bottom')+
  scale_x_date(date_breaks = "2 month", date_labels =  "%b %Y")+
  theme(text = element_text(size=14), axis.title = element_text(size=16), plot.title = element_text(size=20),plot.subtitle=element_text(size=18),
        legend.title = element_text(size=16),strip.text = element_text(size=14),
        axis.text = element_text(size=9),legend.position="none")

sari_brazil_plot<-ggplot()+
  scale_colour_viridis_d(name="Cumulative proportions reported at delay week:", begin=0.1, end=0.85)+
  geom_point(data=Data_delay_prop_total,aes(x=log(total),y=probit(count), colour=as.factor(cum_prop)),alpha=0.5)+#facet_wrap(~s_full,scales="free", ncol=3)+
   geom_smooth(data=Data_delay_prop_total,alpha=0.1,aes(x=log(total),y=probit(count), colour=as.factor(cum_prop)))+
  labs(x="Total SARI hospitalisations (log scale)",y="Cumulative proportions (probit scale)",
       subtitle = 'Against total SARI hospitalisations')+theme_minimal()+
  theme(legend.position = 'bottom')+
  theme(text = element_text(size=14), axis.title = element_text(size=16), plot.title = element_text(size=20),plot.subtitle=element_text(size=18),
        legend.title = element_text(size=16),strip.text = element_text(size=14),
        axis.text = element_text(size=9),legend.position="bottom")
#ggsave(sari_brazil_plot,file='sari_brazil_plot.pdf',width=9,height=5)

cumulative_prop_plots<-grid.arrange(cum_plot_time,sari_brazil_plot,ncol=1)
ggsave(cumulative_prop_plots,file='cumulative_prop_plots.pdf',width=9,height=8)

Data_delay_sum<-as.tibble(Data_wide)
# Data_delay_sum$prop_0<-as.numeric(unlist((Data_wide[,3])/apply(Data_wide[,3:22],1, sum, na.rm=TRUE)))
# Data_delay_sum$prop_1<-as.numeric(unlist((Data_wide[,4])/(apply(Data_wide[,3:22],1, sum, na.rm=TRUE)-Data_wide[,3])))
# Data_delay_sum$prop_2<-as.numeric(unlist((Data_wide[,5])/(apply(Data_wide[,3:22],1, sum, na.rm=TRUE)-apply(Data_wide[,3:4],1, sum, na.rm=TRUE))))
# Data_delay_sum$prop_3<-as.numeric(unlist((Data_wide[,6])/(apply(Data_wide[,3:22],1, sum, na.rm=TRUE)-apply(Data_wide[,3:5],1, sum, na.rm=TRUE))))
# Data_delay_sum$prop_4<-as.numeric(unlist((Data_wide[,7])/(apply(Data_wide[,3:22],1, sum, na.rm=TRUE)-apply(Data_wide[,3:6],1, sum, na.rm=TRUE))))
Data_delay_sum$prop_0<-as.numeric(unlist((Data_wide[,3])/apply(Data_wide[,3:22],1, sum, na.rm=TRUE)))
Data_delay_sum$prop_1<-as.numeric(unlist((Data_wide[,4])/(apply(Data_wide[,3:22],1, sum, na.rm=TRUE))))
Data_delay_sum$prop_2<-as.numeric(unlist((Data_wide[,5])/(apply(Data_wide[,3:22],1, sum, na.rm=TRUE))))
Data_delay_sum$prop_3<-as.numeric(unlist((Data_wide[,6])/(apply(Data_wide[,3:22],1, sum, na.rm=TRUE))))
Data_delay_sum$prop_4<-as.numeric(unlist((Data_wide[,7])/(apply(Data_wide[,3:22],1, sum, na.rm=TRUE))))

Data_delay_long<-as.tibble(Data_delay_sum%>%pivot_longer(!c(t,s), names_to = "delay", values_to = "count")%>%mutate(s=paste(federal_units[s])))
Data_delay_pop<-full_join(Data_delay_long,population, by="s")
Data_delay_prop<-Data_delay_pop%>%filter(delay%in%c("prop_0","prop_1","prop_2","prop_3","prop_4"))
Data_delay_prop<-full_join(Data_delay_prop,totals[,c(1,2,dim(totals)[2])],by = join_by(t, s))
Data_delay_prop<-Data_delay_prop%>%mutate(total.per.capita=total/population, cum_prop=rep(0:4,N*S))

Data_delay_prop_total<-Data_delay_prop%>%mutate(cum_counts=count*total)%>%group_by(t,cum_prop)%>%
  summarise(total=sum(total),cum_counts=sum(cum_counts) )%>%mutate(count=cum_counts/total)

sari_brazil_plot<-ggplot()+
  scale_colour_viridis_d(name="Absolute proportions reported at delay week:", begin=0.1, end=0.85)+
  geom_point(data=Data_delay_prop_total,aes(x=log(total),y=probit(count), colour=as.factor(cum_prop)),alpha=0.5)+#facet_wrap(~s_full,scales="free", ncol=3)+
  geom_smooth(data=Data_delay_prop_total,alpha=0.1,aes(x=log(total),y=probit(count), colour=as.factor(cum_prop)))+
  labs(x="Total SARI hospitalisations",y="Absolute proportion reported",
       subtitle = 'SARI hospitalisations in Brazil')+theme_minimal()+
  theme(legend.position = 'bottom')+
  theme(text = element_text(size=14), axis.title = element_text(size=16), plot.title = element_text(size=20),plot.subtitle=element_text(size=18),
        legend.title = element_text(size=16),strip.text = element_text(size=14),
        axis.text = element_text(size=9),legend.position="bottom")
sari_brazil_plot
ggsave(sari_brazil_plot,file='sari_brazil_plot_absolute.pdf',width=9,height=5)

sari_cumulative_exploratory_3<-ggplot(filter(Data_delay_prop,cum_prop%in%c(0,2,4)),aes(x=total,y=probit(count), colour=as.factor(cum_prop)))+
  geom_point(alpha=0.25)+facet_wrap(~s_full,scales="free", ncol=3)+
  geom_smooth(alpha=0.1)+scale_x_log10()+
  scale_colour_viridis_d(name="Absolute proportions reported at delay week:", begin=0.1, end=0.85)+
  labs(x="Eventual total SARI hospitaltions reported (log scale)",y="Absolute proportions of SARI hospitalisations reported (probit scale)",
       title = 'SARI hospitalisations cumulative proportions reported')+theme_minimal()+
  theme(legend.position = 'bottom')+
  theme(text = element_text(size=14), axis.title = element_text(size=16), plot.title = element_text(size=20),plot.subtitle=element_text(size=18),
        legend.title = element_text(size=16),strip.text = element_text(size=14),
        axis.text = element_text(size=9),legend.position="bottom")
sari_cumulative_exploratory_3
ggsave(sari_cumulative_exploratory_3,file='sari_cumulative_exploratory_3_absolute.pdf',width=10,height=16.5)


abs_lines<-numeric()
for(i in 1:S){
  for(d in 1:5){
    abs_gam<-gam(data = filter(Data_delay_prop,s==federal_units[i],cum_prop==(d-1)), count~s(t,k=20), family = gaussian(link = 'log'))
    abs_lines<-rbind(abs_lines,left_join(tibble(count=(abs_gam$fitted.values),t=1:103)%>%mutate(s=federal_units[i], cum_prop=d-1),population, by='s'))
  }
  # abs_lines<-rbind(abs_lines,left_join(matrix((abs_gam$fitted.values),ncol=103)%>%melt(varnames=c('cum_prop','t'),value.name="count")%>%mutate(s=federal_units[i],cum_prop=cum_prop-1),population, by='s'))
}
abs_plot_time<-ggplot(Data_delay_prop,aes(x=as.Date(first_date)+t*7,y=(count), colour=as.factor(cum_prop)))+
  geom_point(alpha=0.5)+facet_wrap(~s_full,scales="free", ncol=3)+
  geom_line(data=abs_lines, aes(x=as.Date(first_date)+t*7, y=(count)))+
  scale_colour_viridis_d(name="Absolute proportions reported at delay week:", begin=0.1, end=0.85)+
  labs(x="Date",y="Absolute proportions of hospitalisations reported",
       title = 'Absolute proportions of SARI hospitalisations reported by federative unit')+theme_minimal()+
  theme(legend.position = 'bottom')+
  scale_x_date(breaks = "4 month", date_labels =  "%b %y")+
  theme(text = element_text(size=14), axis.title = element_text(size=16), plot.title = element_text(size=20),plot.subtitle=element_text(size=18),
        legend.title = element_text(size=16),strip.text = element_text(size=14),
        axis.text = element_text(size=9),legend.position="bottom")
ggsave(abs_plot_time,file='absolute_plot_time.pdf',width=10,height=12)

sari_brazil_plot<-ggplot()+
  scale_colour_viridis_d(name="Absolute proportions reported at delay week:", begin=0.1, end=0.85)+
  geom_point(data=Data_delay_prop_total,aes(x=log(total),y=(count), colour=as.factor(cum_prop)),alpha=0.5)+#facet_wrap(~s_full,scales="free", ncol=3)+
  geom_smooth(data=Data_delay_prop_total,alpha=0.1,aes(x=log(total),y=(count), colour=as.factor(cum_prop)))+
  labs(x="Total SARI hospitalisations",y="Absolute proportion reported",
       subtitle = 'Against total SARI hospitalisations')+theme_minimal()+
  theme(legend.position = 'bottom')+
  theme(text = element_text(size=14), axis.title = element_text(size=16), plot.title = element_text(size=20),plot.subtitle=element_text(size=18),
        legend.title = element_text(size=16),strip.text = element_text(size=14),
        axis.text = element_text(size=9),legend.position="bottom")
#ggsave(sari_brazil_plot,file='sari_brazil_plot.pdf',width=9,height=5)

absolute_prop_plots<-grid.arrange(abs_plot_time,sari_brazil_plot,ncol=1)
ggsave(absolute_prop_plots,file='absolute_prop_plots.pdf',width=9,height=8)

