## survivor MODEL ##
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
library(geobr)
library(rvest)

# Reading in the table from Wikipedia
page = read_html("https://en.wikipedia.org/wiki/Federative_units_of_Brazil")
# Obtain the piece of the web page that corresponds to the "wikitable" node
my.table = html_node(page, ".wikitable")
# Convert the html table element into a data frame
my.table = html_table(my.table, fill = TRUE)


#setwd("~/OneDrive/PhD/Dengue 2023")
#setwd("~/OneDrive - University of Glasgow/PhD/Dengue 2023")
#setwd("~/Documents/Dengue 2023")

states<-read_state(year=2013)

data_set <- read_csv("data set.csv")
data_set<-data_set[,-c(22:48)]


data<-data_set%>%mutate(start_of_symptoms=as.Date((DT_SIN_PRI),format="%Y-%m-%d"),
                        notification_date=as.Date((DT_NOTIFIC),format="%Y-%m-%d"))%>%
  mutate(delay=difftime(notification_date,start_of_symptoms, unit=c("days")))%>%
  select(start_of_symptoms,notification_date,delay,CLASSI_FIN,SG_UF_NOT)%>%
  drop_na(start_of_symptoms)%>%filter(start_of_symptoms>"2012-12-31",start_of_symptoms<"2020-12-27")

# Make matrix with onset day as rows and delay days as columns
first_day<-min(data$start_of_symptoms)
data_long<-filter(data, delay>=0)%>%
  group_by(start_of_symptoms,notification_date,delay,CLASSI_FIN,SG_UF_NOT)%>%
  summarise(cases=n())%>%
  mutate(delay=as.numeric(delay), 
         day=as.numeric(difftime(start_of_symptoms,first_day, unit=c("days"))))%>%
  as.data.frame()

# calculate weekly delays
data_long$delay_weeks<-floor(data_long$delay/7)
data_long$week<-floor(data_long$day/7)


# Order data frame 
data_long<-data_long[order(data_long[,5],data_long[,9]),]

data_all<-data_long%>%group_by(delay_weeks,week,CLASSI_FIN,SG_UF_NOT)%>%summarise(cases=sum(cases))
colnames(data_all)[4]<-"code_state"
data_regions<-full_join(states,data_all,by="code_state")

ggplot(filter(data_all,delay_weeks<10))+
  geom_point(aes(x=week,y=cases,colour=delay_weeks))+facet_wrap(~code_state)+
  scale_y_log10()




# group by delay week
data_dengue<-data_long%>%group_by(delay_weeks,week,SG_UF_NOT)%>%summarise(cases=sum(cases))
data_remainder<-data_dengue%>%filter(delay_weeks>30)%>%group_by(week,SG_UF_NOT)%>%summarise(cases=sum(cases))%>%mutate(delay_weeks=31)
Dengue_remainder<-full_join(filter(data_dengue,delay_weeks<=30),data_remainder,by=c('week','SG_UF_NOT','delay_weeks','cases'))
colnames(Dengue_remainder)<-c('d','t','s','z')


# Number of regions.
S<-length(unique(Dengue_remainder$s))
# Number of delays.
D_N<-length(unique(Dengue_remainder$d))
# Number of weeks.
N<-length(unique(Dengue_remainder$t))

# Create guide to convert between code_state and name_state.
federal_data<-tibble(s=c(1:S),code_state=unique(sort(unique(data_long$SG_UF_NOT))))
states_short<-as.matrix(states)[,1:3]%>%as.tibble()
states_short$code_state<-as.numeric(states_short$code_state)
states_short$abbrev_state<-as.character(states_short$abbrev_state)
federal_guide<-full_join((federal_data),states_short,by="code_state")

# convert s to numeric 1:27.
s_convert<-numeric(length(Dengue_remainder$s))
for(i in 1:length(Dengue_remainder$s)){
  s_convert[i]<-which(federal_guide$code_state==Dengue_remainder$s[i])
}

# Fill any un-reported weeks with zeros.
Dengue_left<-Dengue_remainder%>%mutate(t=t+1,d=d+1)
Dengue_left$s<-s_convert
all_tsd<-tibble(s=rep(1:S,N*D_N),
                t=rep(sort(rep(1:N,S)),D_N),
                d=sort(rep(1:D_N,N*S)),
                z=0)
Dengue_all<-full_join(Dengue_left,all_tsd,by=c('t','s','d'))%>%mutate(z=z.x+z.y)%>%mutate(z.x=NULL,z.y=NULL)


Data_wide<-Dengue_all%>%ungroup()%>%
  pivot_wider(id_cols=c(t,s),names_from = d,values_from = z)%>%as.matrix()
Data_wide<-Data_wide[order(Data_wide[,2],Data_wide[,1]),]
Dengue_z<-as.numeric(Data_wide[,-(1:2)])%>%array(dim=c(N,S,D_N))%>%aperm(c(2,1,3))
Dengue_z[is.na(Dengue_z)]<-0
# Matrix for total counts [REGION,TIME]
Dengue_y<-apply(Dengue_z,c(1,2),sum)
# 
# load("~/Library/CloudStorage/OneDrive-UniversityofGlasgow/PhD/Dengue 2023/Dengue_haz_spline27_delay1.RData")
# load("~/Library/CloudStorage/OneDrive-UniversityofGlasgow/PhD/Dengue 2023/Dengue_haz_bin27_1delay.RData")
# load("~/Library/CloudStorage/OneDrive-UniversityofGlasgow/PhD/Dengue 2023/Dengue_haz_sea27.RData")
# load("~/Library/CloudStorage/OneDrive-UniversityofGlasgow/PhD/Dengue 2023/Dengue_sea27.RData")
# load("~/Library/CloudStorage/OneDrive-UniversityofGlasgow/PhD/Dengue 2023/Dengue_surv_bin27.RData")
load("~/media/alba/Disk 1/OneDrive/PhD/Case load case studies/Dengue 2023/Dengue_surv_spline_D4.RData")
load("~/media/alba/Disk 1/OneDrive/PhD/Case load case studies/Dengue 2023/Dengue_surv_sea_D4.RData")
load("~/media/alba/Disk 1/OneDrive/PhD/Case load case studies/Dengue 2023/Dengue_surv_bin_D4.RData")
load("~/media/alba/Disk 1/OneDrive/PhD/Case load case studies/Dengue 2023/Dengue_surv_link_D4.RData")


# Explore output
Dengue_logy<-list()
Dengue_combined_logy<-list()
Dengue_bin<-list()
Dengue_combined_bin<-list()
Dengue_haz<-list()
Dengue_combined_haz<-list()
Dengue_sea<-Dengue_combined_sea<-list()
Dengue_bin_hazard<-list()
Dengue_combined_haz_bin<-list()
Dengue_haz_logspline<-list()
Dengue_combined_haz_logspline<-list()

Dengue_combined_sea[[1]] <- as_tibble(do.call('rbind',Dengue_surv_sea[[1]]))
n_sim <- dim(Dengue_combined_sea[[1]])[1]

Dengue_bin_quantiles<-Dengue_haz_bin_quantiles<-list()
Dengue_sea_quantiles<-Dengue_haz_sea_quantiles<-list()
Dengue_spline_quantiles<-Dengue_haz_spline_quantiles<-list()

D_all=32
D_max=30
n_iter=5000
n_burn=2500
n_thin=5
n_chains=2
n_knots=c(18,7,10) 
D=5
N_max<-ncol(Dengue_y)

Nowcast_list=floor(seq(from=2*52, to=375, length=4))


# seasonal model (no incidence-delay) quantiles

Dengue_sea<-list()
Dengue_combined_sea<-Dengue_sea_quantiles<-list()
for(j in length(Nowcast_list):1){
  Dengue_sea[[j]]<-as.mcmc.list(Dengue_surv_sea[[j]])
  # Combine all MCMC chains.
  Dengue_combined_sea[[j]] <- as_tibble(do.call('rbind',Dengue_sea[[j]]))
}

n_sim <- dim(Dengue_combined_sea[[1]])[1]

for(j in length(Nowcast_list):1){
  Dengue_sea_quantiles[[j]]<-select(Dengue_combined_sea[[j]],starts_with('y'))%>%as.matrix()%>%array(dim=c(n_sim,S,Nowcast_list[j]))%>%
    apply(c(2,3),quantile,c(0.025,0.5,0.975))%>%melt(varnames=c('quantile','s','t'),value.name='y')%>%#mutate(t=t+Nowcast_list[j]-W)%>%
    spread(quantile,y)%>%mutate(Nowcast_date=paste(Nowcast_list[j]),s=federal_guide$abbrev_state[s])
  
  
}
# 
# # incidence-delay model quantiles 
Dengue_link<-list()
Dengue_combined_link<-Dengue_link_quantiles<-list()
for(j in length(Nowcast_list):1){
  Dengue_link[[j]]<-as.mcmc.list(Dengue_surv_link[[j]])
  # Combine all MCMC chains.
  Dengue_combined_link[[j]] <- as_tibble(do.call('rbind',Dengue_link[[j]]))
}

n_sim <- dim(Dengue_combined_link[[1]])[1]

for(j in length(Nowcast_list):1){
  Dengue_link_quantiles[[j]]<-select(Dengue_combined_link[[j]],starts_with('y'))%>%as.matrix()%>%array(dim=c(n_sim,S,Nowcast_list[j]))%>%
    apply(c(2,3),quantile,c(0.025,0.5,0.975))%>%melt(varnames=c('quantile','s','t'),value.name='y')%>%#mutate(t=t+Nowcast_list[j]-W)%>%
    spread(quantile,y)%>%mutate(Nowcast_date=paste(Nowcast_list[j]),s=federal_guide$abbrev_state[s])


}



# time bin model quantiles

Dengue_bin<-list()
Dengue_combined_bin<-Dengue_bin_quantiles<-list()
for(j in length(Nowcast_list):1){
  Dengue_bin[[j]]<-as.mcmc.list(Dengue_surv_bin[[j]])
  # Combine all MCMC chains.
  Dengue_combined_bin[[j]] <- as_tibble(do.call('rbind',Dengue_bin[[j]]))
}

n_sim <- dim(Dengue_combined_bin[[1]])[1]

for(j in length(Nowcast_list):1){
  Dengue_bin_quantiles[[j]]<-select(Dengue_combined_bin[[j]],starts_with('y'))%>%as.matrix()%>%array(dim=c(n_sim,S,Nowcast_list[j]))%>%
    apply(c(2,3),quantile,c(0.025,0.5,0.975))%>%melt(varnames=c('quantile','s','t'),value.name='y')%>%#mutate(t=t+Nowcast_list[j]-W)%>%
    spread(quantile,y)%>%mutate(Nowcast_date=paste(Nowcast_list[j]),s=federal_guide$abbrev_state[s])
  
  
}

# spline model quantiles 
Dengue_spline<-list()
Dengue_combined_spline<-Dengue_spline_quantiles<-list()
for(j in length(Nowcast_list):1){
  Dengue_spline[[j]]<-as.mcmc.list(Dengue_surv_logysmooth[[j]])
  # Combine all MCMC chains.
  Dengue_combined_spline[[j]] <- as_tibble(do.call('rbind',Dengue_spline[[j]]))
}

n_sim <- dim(Dengue_combined_spline[[1]])[1]

for(j in length(Nowcast_list):1){
  Dengue_spline_quantiles[[j]]<-select(Dengue_combined_spline[[j]],starts_with('y'))%>%as.matrix()%>%array(dim=c(n_sim,S,Nowcast_list[j]))%>%
    apply(c(2,3),quantile,c(0.025,0.5,0.975))%>%melt(varnames=c('quantile','s','t'),value.name='y')%>%#mutate(t=t+Nowcast_list[j]-W)%>%
    spread(quantile,y)%>%mutate(Nowcast_date=paste(Nowcast_list[j]),s=federal_guide$abbrev_state[s])
  
  
}


#compare models 
library(ggh4x)
y_all<-melt(t(Dengue_y[1:S,]))
colnames(y_all)<-c("t","s","y")
y_all$s<-federal_guide$abbrev_state[y_all$s]

# Set up plots:
N_max<-dim(Dengue_y)[2]

# 
# # 
# # # Survivor and hazard sea 
# dengue_nowcasts_compare<-list()
# for(j in 1:length(Nowcast_list)){
#   dengue_nowcasts_compare[[j]] <- rbind(
#     full_join(Dengue_haz_sea_quantiles[[j]],filter(y_all, t%in%c((1):Nowcast_list[j])),by=c('s','t'))%>%mutate(m=c('Dengue hazard'),d=t-Nowcast_list[j],now=j),
#     full_join(Dengue_sea_quantiles[[j]],filter(y_all, t%in%c((1):Nowcast_list[j])),by=c('s','t'))%>%mutate(m=c('Dengue survivor'),d=t-Nowcast_list[j],now=j)
#   )%>%filter(d>-D)
# }
# dengue_nowcast_dmax<-array(unlist(dengue_nowcasts_compare), dim=c(S*D*2,ncol(dengue_nowcasts_compare[[1]]),length(Nowcast_list)))%>%
#   apply(c(2),c)%>%as_tibble()
# colnames(dengue_nowcast_dmax)<-colnames(dengue_nowcasts_compare[[1]])
# dengue_nowcast_dmax<-dengue_nowcast_dmax%>%mutate(t=as.numeric(t),`2.5%`=as.numeric(`2.5%`),`50%`=as.numeric(`50%`),`97.5%`=as.numeric(`97.5%`),
#                                                   m=factor(m,levels=c('Dengue hazard','Dengue survivor'
#                                                   )),y=as.numeric(y),
#                                                   d=as.numeric(d),now=as.numeric(now))
# dengue_nowcast_dmax_summary<-dengue_nowcast_dmax%>%group_by(m,d,s)%>%
#   summarise(Coverage=mean(y>=`2.5%`&y<=`97.5%`),
#             `Prediction interval width`=mean(`97.5%`-`2.5%`),
#             `Mean absolute error`=mean(abs(`50%`-y)))%>%
#   pivot_longer(cols=c('Coverage','Mean absolute error','Prediction interval width'),names_to = 'metric')%>%
#   mutate(m=factor(m,levels=c('Dengue hazard','Dengue survivor'
#   )),
#   metric=factor(metric,levels=c('Mean absolute error','Prediction interval width','Coverage')))
# 
# 
# #plot_data_dengue<-filter(dengue_nowcast_dmax_summary, d==0)%>%pivot_wider(names_from=c('m'))%>%filter(metric=="Mean Absolute Error")
# #ggplot(plot_data_dengue)+   geom_text(aes(x=`CLR`,y=`CLR Age`/`CLR`,label=s))+geom_hline(yintercept = 1)
# dengue_nowcast_dmax_summary<-filter(dengue_nowcast_dmax_summary, s!="ES")
# dengue_nowcast_dmax_summary_average<-dengue_nowcast_dmax_summary%>%group_by(m,d,metric)%>%
#   summarise(value.mean=mean(value), value.median=median(value))
# dengue_nowcast_dmax_summary_average$value.mean[dengue_nowcast_dmax_summary_average$metric%in%c("Prediction Interval Width","Mean Absolute Error")]<-NA
# dengue_nowcast_dmax_summary_average$value.median[dengue_nowcast_dmax_summary_average$metric%in%c("Coverage")]<-NA
# dengue_nowcast_dmax_summary_average<-dengue_nowcast_dmax_summary_average%>%group_by(m,d,metric)%>%mutate(value=sum(value.mean,value.median,na.rm=TRUE))
# 
# # Plot of prediction performance metrics.
# sea_hazard_survivor <- ggplot(dengue_nowcast_dmax_summary_average)+
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

# Survior sea and link 
dengue_nowcasts_compare<-list()
for(j in 1:length(Nowcast_list)){
  dengue_nowcasts_compare[[j]] <- rbind(
    full_join(Dengue_sea_quantiles[[j]],filter(y_all, t%in%c((1):Nowcast_list[j])),by=c('s','t'))%>%mutate(m=c('Dengue survivor'),d=t-Nowcast_list[j],now=j),
    full_join(Dengue_bin_quantiles[[j]],filter(y_all, t%in%c((1):Nowcast_list[j])),by=c('s','t'))%>%mutate(m=c('Dengue survivor incidence-delay (time category)'),d=t-Nowcast_list[j],now=j)
    )%>%filter(d>-D)
}
dengue_nowcast_dmax<-array(unlist(dengue_nowcasts_compare), dim=c(S*D*2,ncol(dengue_nowcasts_compare[[1]]),length(Nowcast_list)))%>%
  apply(c(2),c)%>%as_tibble()
colnames(dengue_nowcast_dmax)<-colnames(dengue_nowcasts_compare[[1]])   
dengue_nowcast_dmax2<-dengue_nowcast_dmax%>%mutate(t=as.numeric(t),`2.5%`=as.numeric(`2.5%`),`50%`=as.numeric(`50%`),`97.5%`=as.numeric(`97.5%`),
                                                  m=factor(m,levels=c('Dengue survivor','Dengue survivor incidence-delay (time category)'
                                                  )),y=as.numeric(y),
                                                  d=as.numeric(d),now=as.numeric(now))
dengue_nowcast_dmax_summary<-dengue_nowcast_dmax2%>%group_by(m,d,s)%>%
  summarise(Coverage=mean(y>=`2.5%`&y<=`97.5%`),
            `Prediction interval width`=mean(`97.5%`-`2.5%`),
            `Mean absolute error`=mean(abs(`50%`-y)))%>%
  pivot_longer(cols=c('Coverage','Mean absolute error','Prediction interval width'),names_to = 'metric')%>%
  mutate(m=factor(m,levels=c('Dengue survivor','Dengue survivor incidence-delay (time category)'
  )),
  metric=factor(metric,levels=c('Mean absolute error','Prediction interval width','Coverage')))


#plot_data_dengue<-filter(dengue_nowcast_dmax_summary, d==0)%>%pivot_wider(names_from=c('m'))%>%filter(metric=="Mean Absolute Error") 
#ggplot(plot_data_dengue)+   geom_text(aes(x=`CLR`,y=`CLR Age`/`CLR`,label=s))+geom_hline(yintercept = 1)
#dengue_nowcast_dmax_summary<-filter(dengue_nowcast_dmax_summary, s!="ES")
dengue_nowcast_dmax_summary_average<-dengue_nowcast_dmax_summary%>%group_by(m,d,metric)%>%
  summarise(value.mean=mean(value), value.median=median(value))
dengue_nowcast_dmax_summary_average$value.mean[dengue_nowcast_dmax_summary_average$metric%in%c("Prediction Interval Width","Mean Absolute Error")]<-NA
dengue_nowcast_dmax_summary_average$value.median[dengue_nowcast_dmax_summary_average$metric%in%c("Coverage")]<-NA
dengue_nowcast_dmax_summary_average<-dengue_nowcast_dmax_summary_average%>%group_by(m,d,metric)%>%mutate(value=sum(value.mean,value.median,na.rm=TRUE))

# Plot of prediction performance metrics.
sea_bin_survivor <- ggplot(dengue_nowcast_dmax_summary_average)+
  geom_line(aes(x=d,y=value,colour=m),alpha=0.6)+
  geom_point(aes(x=d,y=value,colour=m,shape=m),alpha=0.6)+
  facet_grid2(~metric, scales = c("free"),independent = c("all"))+
  scale_colour_manual(name=NULL,values = c( "#0072B2", "darkgreen"))+
  scale_shape_manual(name=NULL,values=c(15,3))+
  scale_x_continuous()+
  labs(x='Prediction time difference (weeks)',y=NULL,
       title="Brazil arbovirus cases",
       subtitle='Performance metrics from the rolling prediction experiment')+
  theme_minimal()+
  theme(legend.position = "bottom",text = element_text(size=12), plot.title = element_text(size=16),
        plot.subtitle=element_text(size=14), legend.title = element_text(size=14),
        strip.text = element_text(size=14))
sea_bin_survivor
ggsave(sea_bin_survivor, file="plots/sea_bin_survivor.pdf", width=9, height=4)


# 
# # # hazard link and sea
# dengue_nowcasts_compare<-list()
# for(j in 1:length(Nowcast_list)){
#   dengue_nowcasts_compare[[j]] <- rbind(
#     full_join(Dengue_haz_sea_quantiles[[j]],filter(y_all, t%in%c((1):Nowcast_list[j])),by=c('s','t'))%>%mutate(m=c('Dengue hazard'),d=t-Nowcast_list[j],now=j),
#     full_join(Dengue_haz_bin_quantiles[[j]],filter(y_all, t%in%c((1):Nowcast_list[j])),by=c('s','t'))%>%mutate(m=c('Dengue hazard link'),d=t-Nowcast_list[j],now=j)
#   )%>%filter(d>-D)
# }
# dengue_nowcast_dmax<-array(unlist(dengue_nowcasts_compare), dim=c(S*D*2,ncol(dengue_nowcasts_compare[[1]]),length(Nowcast_list)))%>%
#   apply(c(2),c)%>%as_tibble()
# colnames(dengue_nowcast_dmax)<-colnames(dengue_nowcasts_compare[[1]])
# dengue_nowcast_dmax<-dengue_nowcast_dmax%>%mutate(t=as.numeric(t),`2.5%`=as.numeric(`2.5%`),`50%`=as.numeric(`50%`),`97.5%`=as.numeric(`97.5%`),
#                                                   m=factor(m,levels=c('Dengue hazard','Dengue hazard link'
#                                                   )),y=as.numeric(y),
#                                                   d=as.numeric(d),now=as.numeric(now))
# dengue_nowcast_dmax_summary<-dengue_nowcast_dmax%>%group_by(m,d,s)%>%
#   summarise(Coverage=mean(y>=`2.5%`&y<=`97.5%`),
#             `Prediction interval width`=mean(`97.5%`-`2.5%`),
#             `Mean absolute error`=mean(abs(`50%`-y)))%>%
#   pivot_longer(cols=c('Coverage','Mean absolute error','Prediction interval width'),names_to = 'metric')%>%
#   mutate(m=factor(m,levels=c('Dengue hazard','Dengue hazard link'
#   )),
#   metric=factor(metric,levels=c('Mean absolute error','Prediction interval width','Coverage')))
# 
# 
# #plot_data_dengue<-filter(dengue_nowcast_dmax_summary, d==0)%>%pivot_wider(names_from=c('m'))%>%filter(metric=="Mean Absolute Error")
# #ggplot(plot_data_dengue)+   geom_text(aes(x=`CLR`,y=`CLR Age`/`CLR`,label=s))+geom_hline(yintercept = 1)
# dengue_nowcast_dmax_summary<-filter(dengue_nowcast_dmax_summary, s!="ES")
# dengue_nowcast_dmax_summary_average<-dengue_nowcast_dmax_summary%>%group_by(m,d,metric)%>%
#   summarise(value.mean=mean(value), value.median=median(value))
# dengue_nowcast_dmax_summary_average$value.mean[dengue_nowcast_dmax_summary_average$metric%in%c("Prediction Interval Width","Mean Absolute Error")]<-NA
# dengue_nowcast_dmax_summary_average$value.median[dengue_nowcast_dmax_summary_average$metric%in%c("Coverage")]<-NA
# dengue_nowcast_dmax_summary_average<-dengue_nowcast_dmax_summary_average%>%group_by(m,d,metric)%>%mutate(value=sum(value.mean,value.median,na.rm=TRUE))
# 
# # Plot of prediction performance metrics.
# link_sea_hazard <- ggplot(dengue_nowcast_dmax_summary_average)+
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
# link_sea_hazard

# 
# # Survivor link and hazard link 
# dengue_nowcasts_compare<-list()
# for(j in 1:length(Nowcast_list)){
#   dengue_nowcasts_compare[[j]] <- rbind(
#     full_join(Dengue_bin_quantiles[[j]],filter(y_all, t%in%c((1):Nowcast_list[j])),by=c('s','t'))%>%mutate(m=c('Dengue survivor link'),d=t-Nowcast_list[j],now=j),
#     full_join(Dengue_haz_bin_quantiles[[j]],filter(y_all, t%in%c((1):Nowcast_list[j])),by=c('s','t'))%>%mutate(m=c('Dengue hazard link'),d=t-Nowcast_list[j],now=j)
#   )%>%filter(d>-D)
# }
# dengue_nowcast_dmax<-array(unlist(dengue_nowcasts_compare), dim=c(S*D*2,ncol(dengue_nowcasts_compare[[1]]),length(Nowcast_list)))%>%
#   apply(c(2),c)%>%as_tibble()
# colnames(dengue_nowcast_dmax)<-colnames(dengue_nowcasts_compare[[1]])   
# dengue_nowcast_dmax<-dengue_nowcast_dmax%>%mutate(t=as.numeric(t),`2.5%`=as.numeric(`2.5%`),`50%`=as.numeric(`50%`),`97.5%`=as.numeric(`97.5%`),
#                                                   m=factor(m,levels=c('Dengue survivor link','Dengue hazard link'
#                                                   )),y=as.numeric(y),
#                                                   d=as.numeric(d),now=as.numeric(now))
# dengue_nowcast_dmax_summary<-dengue_nowcast_dmax%>%group_by(m,d,s)%>%
#   summarise(Coverage=mean(y>=`2.5%`&y<=`97.5%`),
#             `Prediction interval width`=mean(`97.5%`-`2.5%`),
#             `Mean absolute error`=mean(abs(`50%`-y)))%>%
#   pivot_longer(cols=c('Coverage','Mean absolute error','Prediction interval width'),names_to = 'metric')%>%
#   mutate(m=factor(m,levels=c('Dengue survivor link','Dengue hazard link'
#   )),
#   metric=factor(metric,levels=c('Mean absolute error','Prediction interval width','Coverage')))
# 
# 
# #plot_data_dengue<-filter(dengue_nowcast_dmax_summary, d==0)%>%pivot_wider(names_from=c('m'))%>%filter(metric=="Mean Absolute Error") 
# #ggplot(plot_data_dengue)+   geom_text(aes(x=`CLR`,y=`CLR Age`/`CLR`,label=s))+geom_hline(yintercept = 1)
# dengue_nowcast_dmax_summary<-filter(dengue_nowcast_dmax_summary, s!="ES")
# dengue_nowcast_dmax_summary_average<-dengue_nowcast_dmax_summary%>%group_by(m,d,metric)%>%
#   summarise(value.mean=mean(value), value.median=median(value))
# dengue_nowcast_dmax_summary_average$value.mean[dengue_nowcast_dmax_summary_average$metric%in%c("Prediction Interval Width","Mean Absolute Error")]<-NA
# dengue_nowcast_dmax_summary_average$value.median[dengue_nowcast_dmax_summary_average$metric%in%c("Coverage")]<-NA
# dengue_nowcast_dmax_summary_average<-dengue_nowcast_dmax_summary_average%>%group_by(m,d,metric)%>%mutate(value=sum(value.mean,value.median,na.rm=TRUE))
# 
# # Plot of prediction performance metrics.
# link_survivor_hazard <- ggplot(dengue_nowcast_dmax_summary_average)+
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
# link_survivor_hazard

# Survivor sea and spline
dengue_nowcasts_compare<-list()
for(j in 1:length(Nowcast_list)){
  dengue_nowcasts_compare[[j]] <- rbind(
    full_join(Dengue_sea_quantiles[[j]],filter(y_all, t%in%c((1):Nowcast_list[j])),by=c('s','t'))%>%mutate(m=c('Dengue survivor'),d=t-Nowcast_list[j],now=j),
    full_join(Dengue_spline_quantiles[[j]],filter(y_all, t%in%c((1):Nowcast_list[j])),by=c('s','t'))%>%mutate(m=c('Dengue survivor spline link'),d=t-Nowcast_list[j],now=j)
  )%>%filter(d>-D)
}
dengue_nowcast_dmax<-array(unlist(dengue_nowcasts_compare), dim=c(S*D*2,ncol(dengue_nowcasts_compare[[1]]),length(Nowcast_list)))%>%
  apply(c(2),c)%>%as_tibble()
colnames(dengue_nowcast_dmax)<-colnames(dengue_nowcasts_compare[[1]])   
dengue_nowcast_dmax<-dengue_nowcast_dmax%>%mutate(t=as.numeric(t),`2.5%`=as.numeric(`2.5%`),`50%`=as.numeric(`50%`),`97.5%`=as.numeric(`97.5%`),
                                                  m=factor(m,levels=c('Dengue survivor','Dengue survivor spline link'
                                                  )),y=as.numeric(y),
                                                  d=as.numeric(d),now=as.numeric(now))
dengue_nowcast_dmax_summary<-dengue_nowcast_dmax%>%group_by(m,d,s)%>%
  summarise(Coverage=mean(y>=`2.5%`&y<=`97.5%`),
            `Prediction interval width`=mean(`97.5%`-`2.5%`),
            `Mean absolute error`=mean(abs(`50%`-y)))%>%
  pivot_longer(cols=c('Coverage','Mean absolute error','Prediction interval width'),names_to = 'metric')%>%
  mutate(m=factor(m,levels=c('Dengue survivor','Dengue survivor spline link'
  )),
  metric=factor(metric,levels=c('Mean absolute error','Prediction interval width','Coverage')))


#plot_data_dengue<-filter(dengue_nowcast_dmax_summary, d==0)%>%pivot_wider(names_from=c('m'))%>%filter(metric=="Mean Absolute Error") 
#ggplot(plot_data_dengue)+   geom_text(aes(x=`CLR`,y=`CLR Age`/`CLR`,label=s))+geom_hline(yintercept = 1)
#dengue_nowcast_dmax_summary<-filter(dengue_nowcast_dmax_summary, s!="ES")
dengue_nowcast_dmax_summary_average<-dengue_nowcast_dmax_summary%>%group_by(m,d,metric)%>%
  summarise(value.mean=mean(value), value.median=median(value))
dengue_nowcast_dmax_summary_average$value.mean[dengue_nowcast_dmax_summary_average$metric%in%c("Prediction Interval Width","Mean Absolute Error")]<-NA
dengue_nowcast_dmax_summary_average$value.median[dengue_nowcast_dmax_summary_average$metric%in%c("Coverage")]<-NA
dengue_nowcast_dmax_summary_average<-dengue_nowcast_dmax_summary_average%>%group_by(m,d,metric)%>%mutate(value=sum(value.mean,value.median,na.rm=TRUE))

# Plot of prediction performance metrics.
sea_spline_survivor <- ggplot(dengue_nowcast_dmax_summary_average)+
  geom_line(aes(x=d,y=value,colour=m),alpha=0.6)+
  geom_point(aes(x=d,y=value,colour=m,shape=m),alpha=0.6)+
  facet_grid2(~metric, scales = c("free"),independent = c("all"))+
  scale_colour_manual(name=NULL,values = c( "#0072B2", "purple"))+
  scale_shape_manual(name=NULL,values=c(15,3))+
  scale_x_continuous()+
  labs(x='Prediction time difference (weeks)',y=NULL,
       title="Brazil arbovirus cases",
       subtitle='Performance metrics from the rolling prediction experiment')+
  theme_minimal()+
  theme(legend.position = "bottom",text = element_text(size=12), plot.title = element_text(size=16),
        plot.subtitle=element_text(size=14), legend.title = element_text(size=14),
        strip.text = element_text(size=14))
sea_spline_survivor
ggsave(sea_spline_survivor, file="sea_spline_survivor.pdf", width=9, height=4)

sea_spline_survivor_notitle <- ggplot(dengue_nowcast_dmax_summary_average)+
  geom_line(aes(x=d,y=value,colour=m),alpha=0.6)+
  geom_point(aes(x=d,y=value,colour=m,shape=m),alpha=0.6)+
  facet_grid2(~metric, scales = c("free"),independent = c("all"))+
  scale_colour_manual(name=NULL,values = c( "#0072B2", "purple"))+
  scale_shape_manual(name=NULL,values=c(15,3))+
  scale_x_continuous()+
  labs(x='Prediction time difference (weeks)',y=NULL,
       title=NULL,
       subtitle=NULL)+
  theme_minimal()+
  theme(legend.position = "bottom",text = element_text(size=12), plot.title = element_text(size=16),
        plot.subtitle=element_text(size=14), legend.title = element_text(size=14),
        strip.text = element_text(size=14))
ggsave(sea_spline_survivor_notitle, file="sea_spline_survivor_notitle.pdf", width=9, height=3.5)


# Survivor sea and survivor spline
# Survivor sea and spline
dengue_nowcasts_compare<-list()
for(j in 1:length(Nowcast_list)){
  dengue_nowcasts_compare[[j]] <- rbind(
    full_join(Dengue_sea_quantiles[[j]],filter(y_all, t%in%c((1):Nowcast_list[j])),by=c('s','t'))%>%mutate(m=c('Arbovirus survivor'),d=t-Nowcast_list[j],now=j),
    full_join(Dengue_bin_quantiles[[j]],filter(y_all, t%in%c((1):Nowcast_list[j])),by=c('s','t'))%>%mutate(m=c('Arbovirus survivor incidence-delay'),d=t-Nowcast_list[j],now=j),
    full_join(Dengue_spline_quantiles[[j]],filter(y_all, t%in%c((1):Nowcast_list[j])),by=c('s','t'))%>%mutate(m=c('Arbovirus survivor spline link'),d=t-Nowcast_list[j],now=j)
  )%>%filter(d>-D)
}
dengue_nowcast_dmax<-array(unlist(dengue_nowcasts_compare), dim=c(S*D*3,ncol(dengue_nowcasts_compare[[1]]),length(Nowcast_list)))%>%
  apply(c(2),c)%>%as_tibble()
colnames(dengue_nowcast_dmax)<-colnames(dengue_nowcasts_compare[[1]])   
dengue_nowcast_dmax<-dengue_nowcast_dmax%>%mutate(t=as.numeric(t),`2.5%`=as.numeric(`2.5%`),`50%`=as.numeric(`50%`),`97.5%`=as.numeric(`97.5%`),
                                                  m=factor(m,levels=c('Arbovirus survivor','Arbovirus survivor incidence-delay','Arbovirus survivor spline link'
                                                  )),y=as.numeric(y),
                                                  d=as.numeric(d),now=as.numeric(now))
dengue_nowcast_dmax_summary<-dengue_nowcast_dmax%>%group_by(m,d,s)%>%
  summarise(Coverage=mean(y>=`2.5%`&y<=`97.5%`),
            `Prediction interval width`=mean(`97.5%`-`2.5%`),
            `Mean absolute error`=mean(abs(`50%`-y)))%>%
  pivot_longer(cols=c('Coverage','Mean absolute error','Prediction interval width'),names_to = 'metric')%>%
  mutate(m=factor(m,levels=c('Arbovirus survivor','Arbovirus survivor incidence-delay','Arbovirus survivor spline link'
  )),
  metric=factor(metric,levels=c('Mean absolute error','Prediction interval width','Coverage')))


#plot_data_dengue<-filter(dengue_nowcast_dmax_summary, d==0)%>%pivot_wider(names_from=c('m'))%>%filter(metric=="Mean Absolute Error") 
#ggplot(plot_data_dengue)+   geom_text(aes(x=`CLR`,y=`CLR Age`/`CLR`,label=s))+geom_hline(yintercept = 1)
#dengue_nowcast_dmax_summary<-filter(dengue_nowcast_dmax_summary, s!="ES")
dengue_nowcast_dmax_summary_average<-dengue_nowcast_dmax_summary%>%group_by(m,d,metric)%>%
  summarise(value.mean=mean(value), value.median=median(value))
dengue_nowcast_dmax_summary_average$value.mean[dengue_nowcast_dmax_summary_average$metric%in%c("Prediction Interval Width","Mean Absolute Error")]<-NA
dengue_nowcast_dmax_summary_average$value.median[dengue_nowcast_dmax_summary_average$metric%in%c("Coverage")]<-NA
dengue_nowcast_dmax_summary_average<-dengue_nowcast_dmax_summary_average%>%group_by(m,d,metric)%>%mutate(value=sum(value.mean,value.median,na.rm=TRUE))

# Plot of prediction performance metrics.
sea_spline_linear <- ggplot(dengue_nowcast_dmax_summary_average)+
  geom_line(aes(x=d,y=value,colour=m),alpha=0.6)+
  geom_point(aes(x=d,y=value,colour=m,shape=m),alpha=0.6)+
  facet_grid2(~metric, scales = c("free"),independent = c("all"))+
  scale_colour_manual(name=NULL,values = c( "#0072B2", "darkgreen","purple" ))+
  scale_shape_manual(name=NULL,values=c(15,3,4))+
  scale_x_continuous()+
  labs(x='Prediction time difference (weeks)',y=NULL,
       title="Brazil arbovirus cases",
       subtitle='Performance metrics from the rolling prediction experiment')+
  theme_minimal()+
  theme(legend.position = "bottom",text = element_text(size=12), plot.title = element_text(size=16),
        plot.subtitle=element_text(size=14), legend.title = element_text(size=14),
        strip.text = element_text(size=14))
sea_spline_linear
ggsave(sea_spline_linear, file="plots/dengue_all_survivor.pdf", width=9, height=4)

# Survivor link and spline
dengue_nowcasts_compare<-list()
for(j in 1:length(Nowcast_list)){
  dengue_nowcasts_compare[[j]] <- rbind(
    full_join(Dengue_bin_quantiles[[j]],filter(y_all, t%in%c((1):Nowcast_list[j])),by=c('s','t'))%>%mutate(m=c('Arbovirus survivor link (time category)'),d=t-Nowcast_list[j],now=j),
    full_join(Dengue_link_quantiles[[j]],filter(y_all, t%in%c((1):Nowcast_list[j])),by=c('s','t'))%>%mutate(m=c('Arbovirus survivor link'),d=t-Nowcast_list[j],now=j)
  )%>%filter(d>-D)
}
dengue_nowcast_dmax<-array(unlist(dengue_nowcasts_compare), dim=c(S*D*2,ncol(dengue_nowcasts_compare[[1]]),length(Nowcast_list)))%>%
  apply(c(2),c)%>%as_tibble()
colnames(dengue_nowcast_dmax)<-colnames(dengue_nowcasts_compare[[1]])   
dengue_nowcast_dmax<-dengue_nowcast_dmax%>%mutate(t=as.numeric(t),`2.5%`=as.numeric(`2.5%`),`50%`=as.numeric(`50%`),`97.5%`=as.numeric(`97.5%`),
                                                  m=factor(m,levels=c('Arbovirus survivor link (time category)','Arbovirus survivor link'
                                                  )),y=as.numeric(y),
                                                  d=as.numeric(d),now=as.numeric(now))
dengue_nowcast_dmax_summary<-dengue_nowcast_dmax%>%group_by(m,d,s)%>%
  summarise(Coverage=mean(y>=`2.5%`&y<=`97.5%`),
            `Prediction interval width`=mean(`97.5%`-`2.5%`),
            `Mean absolute error`=mean(abs(`50%`-y)))%>%
  pivot_longer(cols=c('Coverage','Mean absolute error','Prediction interval width'),names_to = 'metric')%>%
  mutate(m=factor(m,levels=c('Arbovirus survivor link (time category)','Arbovirus survivor link'
  )),
  metric=factor(metric,levels=c('Mean absolute error','Prediction interval width','Coverage')))


#plot_data_dengue<-filter(dengue_nowcast_dmax_summary, d==0)%>%pivot_wider(names_from=c('m'))%>%filter(metric=="Mean Absolute Error") 
#ggplot(plot_data_dengue)+   geom_text(aes(x=`CLR`,y=`CLR Age`/`CLR`,label=s))+geom_hline(yintercept = 1)
dengue_nowcast_dmax_summary<-filter(dengue_nowcast_dmax_summary, s!="ES")
dengue_nowcast_dmax_summary_average<-dengue_nowcast_dmax_summary%>%group_by(m,d,metric)%>%
  summarise(value.mean=mean(value), value.median=median(value))
dengue_nowcast_dmax_summary_average$value.mean[dengue_nowcast_dmax_summary_average$metric%in%c("Prediction Interval Width","Mean Absolute Error")]<-NA
dengue_nowcast_dmax_summary_average$value.median[dengue_nowcast_dmax_summary_average$metric%in%c("Coverage")]<-NA
dengue_nowcast_dmax_summary_average<-dengue_nowcast_dmax_summary_average%>%group_by(m,d,metric)%>%mutate(value=sum(value.mean,value.median,na.rm=TRUE))

# Plot of prediction performance metrics.
link_bin_survivor <- ggplot(dengue_nowcast_dmax_summary_average)+
  geom_line(aes(x=d,y=value,colour=m),alpha=0.6)+
  geom_point(aes(x=d,y=value,colour=m,shape=m),alpha=0.6)+
  facet_grid2(~metric, scales = c("free"),independent = c("all"))+
  scale_colour_manual(name=NULL,values = c( "darkgreen", "#D55E00"))+
  scale_shape_manual(name=NULL,values=c(15,3))+
  scale_x_continuous()+
  labs(x='Prediction time difference (weeks)',y=NULL,
       title='Performance metrics from the rolling prediction experiment')+
  theme_minimal()+
  theme(legend.position = "bottom")
link_bin_survivor



# Survivor sea and survivor lin and survuivor bin
dengue_nowcasts_compare<-list()
for(j in 1:length(Nowcast_list)){
  dengue_nowcasts_compare[[j]] <- rbind(
    full_join(Dengue_sea_quantiles[[j]],filter(y_all, t%in%c((1):Nowcast_list[j])),by=c('s','t'))%>%mutate(m=c('Arbovirus survivor'),d=t-Nowcast_list[j],now=j),
    full_join(Dengue_link_quantiles[[j]],filter(y_all, t%in%c((1):Nowcast_list[j])),by=c('s','t'))%>%mutate(m=c('Arbovirus survivor incidence-delay'),d=t-Nowcast_list[j],now=j),
    full_join(Dengue_bin_quantiles[[j]],filter(y_all, t%in%c((1):Nowcast_list[j])),by=c('s','t'))%>%mutate(m=c('Arbovirus survivor incidence-delay (time category)'),d=t-Nowcast_list[j],now=j)
  )%>%filter(d>-D)
}
dengue_nowcast_dmax<-array(unlist(dengue_nowcasts_compare), dim=c(S*D*3,ncol(dengue_nowcasts_compare[[1]]),length(Nowcast_list)))%>%
  apply(c(2),c)%>%as_tibble()
colnames(dengue_nowcast_dmax)<-colnames(dengue_nowcasts_compare[[1]])   
dengue_nowcast_dmax<-dengue_nowcast_dmax%>%mutate(t=as.numeric(t),`2.5%`=as.numeric(`2.5%`),`50%`=as.numeric(`50%`),`97.5%`=as.numeric(`97.5%`),
                                                  m=factor(m,levels=c('Arbovirus survivor','Arbovirus survivor incidence-delay','Arbovirus survivor incidence-delay (time category)'
                                                  )),y=as.numeric(y),
                                                  d=as.numeric(d),now=as.numeric(now))
dengue_nowcast_dmax_summary<-dengue_nowcast_dmax%>%group_by(m,d,s)%>%
  summarise(Coverage=mean(y>=`2.5%`&y<=`97.5%`),
            `Prediction interval width`=mean(`97.5%`-`2.5%`),
            `Mean absolute error`=mean(abs(`50%`-y)))%>%
  pivot_longer(cols=c('Coverage','Mean absolute error','Prediction interval width'),names_to = 'metric')%>%
  mutate(m=factor(m,levels=c('Arbovirus survivor','Arbovirus survivor incidence-delay','Arbovirus survivor incidence-delay (time category)'
  )),
  metric=factor(metric,levels=c('Mean absolute error','Prediction interval width','Coverage')))


#plot_data_dengue<-filter(dengue_nowcast_dmax_summary, d==0)%>%pivot_wider(names_from=c('m'))%>%filter(metric=="Mean Absolute Error") 
#ggplot(plot_data_dengue)+   geom_text(aes(x=`CLR`,y=`CLR Age`/`CLR`,label=s))+geom_hline(yintercept = 1)
#dengue_nowcast_dmax_summary<-filter(dengue_nowcast_dmax_summary, s!="ES")
dengue_nowcast_dmax_summary_average<-dengue_nowcast_dmax_summary%>%group_by(m,d,metric)%>%
  summarise(value.mean=mean(value), value.median=median(value))
dengue_nowcast_dmax_summary_average$value.mean[dengue_nowcast_dmax_summary_average$metric%in%c("Prediction Interval Width","Mean Absolute Error")]<-NA
dengue_nowcast_dmax_summary_average$value.median[dengue_nowcast_dmax_summary_average$metric%in%c("Coverage")]<-NA
dengue_nowcast_dmax_summary_average<-dengue_nowcast_dmax_summary_average%>%group_by(m,d,metric)%>%mutate(value=sum(value.mean,value.median,na.rm=TRUE))

# Plot of prediction performance metrics.
sea_spline_linear <- ggplot(dengue_nowcast_dmax_summary_average)+
  geom_line(aes(x=d,y=value,colour=m),alpha=0.6)+
  geom_point(aes(x=d,y=value,colour=m,shape=m),alpha=0.6)+
  facet_grid2(~metric, scales = c("free"),independent = c("all"))+
  scale_colour_manual(name=NULL,values = c( "#0072B2", "#D55E00", "darkgreen" ))+
  scale_shape_manual(name=NULL,values=c(15,3,1))+
  scale_x_continuous()+
  labs(x='Prediction time difference (weeks)',y=NULL,
       title="Brazil arbovirus cases",
       subtitle='Performance metrics from the rolling prediction experiment')+
  theme_minimal()+
  theme(legend.position = "bottom",text = element_text(size=12), plot.title = element_text(size=16),
        plot.subtitle=element_text(size=14), legend.title = element_text(size=14),
        strip.text = element_text(size=14))
sea_spline_linear
ggsave(sea_spline_linear, file="plots/dengue_all_survivor.pdf", width=9, height=4)

# 
# 
# 
# # hazard link and spline
# dengue_nowcasts_compare<-list()
# for(j in 1:length(Nowcast_list)){
#   dengue_nowcasts_compare[[j]] <- rbind(
#     full_join(Dengue_haz_bin_quantiles[[j]],filter(y_all, t%in%c((1):Nowcast_list[j])),by=c('s','t'))%>%mutate(m=c('Dengue hazard link'),d=t-Nowcast_list[j],now=j),
#     full_join(Dengue_haz_spline_quantiles[[j]],filter(y_all, t%in%c((1):Nowcast_list[j])),by=c('s','t'))%>%mutate(m=c('Dengue hazard spline'),d=t-Nowcast_list[j],now=j)
#   )%>%filter(d>-D)
# }
# dengue_nowcast_dmax<-array(unlist(dengue_nowcasts_compare), dim=c(S*D*2,ncol(dengue_nowcasts_compare[[1]]),length(Nowcast_list)))%>%
#   apply(c(2),c)%>%as_tibble()
# colnames(dengue_nowcast_dmax)<-colnames(dengue_nowcasts_compare[[1]])   
# dengue_nowcast_dmax<-dengue_nowcast_dmax%>%mutate(t=as.numeric(t),`2.5%`=as.numeric(`2.5%`),`50%`=as.numeric(`50%`),`97.5%`=as.numeric(`97.5%`),
#                                                   m=factor(m,levels=c('Dengue hazard link','Dengue hazard spline'
#                                                   )),y=as.numeric(y),
#                                                   d=as.numeric(d),now=as.numeric(now))
# dengue_nowcast_dmax_summary<-dengue_nowcast_dmax%>%group_by(m,d,s)%>%
#   summarise(Coverage=mean(y>=`2.5%`&y<=`97.5%`),
#             `Prediction interval width`=mean(`97.5%`-`2.5%`),
#             `Mean absolute error`=mean(abs(`50%`-y)))%>%
#   pivot_longer(cols=c('Coverage','Mean absolute error','Prediction interval width'),names_to = 'metric')%>%
#   mutate(m=factor(m,levels=c('Dengue hazard link','Dengue hazard spline'
#   )),
#   metric=factor(metric,levels=c('Mean absolute error','Prediction interval width','Coverage')))
# 
# 
# #plot_data_dengue<-filter(dengue_nowcast_dmax_summary, d==0)%>%pivot_wider(names_from=c('m'))%>%filter(metric=="Mean Absolute Error") 
# #ggplot(plot_data_dengue)+   geom_text(aes(x=`CLR`,y=`CLR Age`/`CLR`,label=s))+geom_hline(yintercept = 1)
# dengue_nowcast_dmax_summary<-filter(dengue_nowcast_dmax_summary, s!="ES")
# dengue_nowcast_dmax_summary_average<-dengue_nowcast_dmax_summary%>%group_by(m,d,metric)%>%
#   summarise(value.mean=mean(value), value.median=median(value))
# dengue_nowcast_dmax_summary_average$value.mean[dengue_nowcast_dmax_summary_average$metric%in%c("Prediction Interval Width","Mean Absolute Error")]<-NA
# dengue_nowcast_dmax_summary_average$value.median[dengue_nowcast_dmax_summary_average$metric%in%c("Coverage")]<-NA
# dengue_nowcast_dmax_summary_average<-dengue_nowcast_dmax_summary_average%>%group_by(m,d,metric)%>%mutate(value=sum(value.mean,value.median,na.rm=TRUE))
# 
# # Plot of prediction performance metrics.
# link_spline_hazard <- ggplot(dengue_nowcast_dmax_summary_average)+
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
# link_spline_hazard



# COMPARE FOR EACH REGION 

# Survivor sea and survivor spline
SARI_nowcasts_compare<-list()
for(j in 1:length(Nowcast_list)){
  SARI_nowcasts_compare[[j]] <- rbind(
    full_join(Dengue_sea_quantiles[[j]],filter(y_all, t%in%c((1):Nowcast_list[j])),by=c('s','t'))%>%mutate(m=c('Arbovirus survivor'),d=t-Nowcast_list[j],now=j),
    full_join(Dengue_spline_quantiles[[j]],filter(y_all, t%in%c((1):Nowcast_list[j])),by=c('s','t'))%>%mutate(m=c('Arbovirus survivor spline'),d=t-Nowcast_list[j],now=j)
  )%>%filter(d>-D)
}
SARI_nowcast_dmax<-array(unlist(SARI_nowcasts_compare), dim=c(S*D*2,ncol(SARI_nowcasts_compare[[1]]),length(Nowcast_list)))%>%
  apply(c(2),c)%>%as_tibble()
colnames(SARI_nowcast_dmax)<-colnames(SARI_nowcasts_compare[[1]])   
SARI_nowcast_dmax<-SARI_nowcast_dmax%>%mutate(t=as.numeric(t),`2.5%`=as.numeric(`2.5%`),`50%`=as.numeric(`50%`),`97.5%`=as.numeric(`97.5%`),
                                              m=factor(m,levels=c('Arbovirus survivor','Arbovirus survivor spline'
                                              )),y=as.numeric(y),
                                              d=as.numeric(d),now=as.numeric(now))
SARI_nowcast_dmax_summary<-SARI_nowcast_dmax%>%group_by(m,d,s)%>%
  summarise(Coverage=mean(y>=`2.5%`&y<=`97.5%`),
            `Prediction interval width`=mean(`97.5%`-`2.5%`),
            `Mean absolute error`=mean(abs(`50%`-y)))%>%
  pivot_longer(cols=c('Coverage','Mean absolute error','Prediction interval width'),names_to = 'metric')%>%
  mutate(m=factor(m,levels=c('Arbovirus survivor','Arbovirus survivor spline'
  )),
  metric=factor(metric,levels=c('Mean absolute error','Prediction interval width','Coverage')))


#plot_data_SARI<-filter(SARI_nowcast_dmax_summary, d==0)%>%pivot_wider(names_from=c('m'))%>%filter(metric=="Mean Absolute Error") 
#ggplot(plot_data_SARI)+   geom_text(aes(x=`CLR`,y=`CLR Age`/`CLR`,label=s))+geom_hline(yintercept = 1)
SARI_nowcast_dmax_summary<-filter(SARI_nowcast_dmax_summary, metric=='Mean absolute error')

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



# Survivor sea and survivor spline
SARI_nowcasts_compare<-list()
for(j in 1:length(Nowcast_list)){
  SARI_nowcasts_compare[[j]] <- rbind(
    full_join(Dengue_sea_quantiles[[j]],filter(y_all, t%in%c((1):Nowcast_list[j])),by=c('s','t'))%>%mutate(m=c('Arbovirus survivor'),d=t-Nowcast_list[j],now=j),
    full_join(Dengue_bin_quantiles[[j]],filter(y_all, t%in%c((1):Nowcast_list[j])),by=c('s','t'))%>%mutate(m=c('Arbovirus survivor linear'),d=t-Nowcast_list[j],now=j)
  )%>%filter(d>-D)
}
SARI_nowcast_dmax<-array(unlist(SARI_nowcasts_compare), dim=c(S*D*2,ncol(SARI_nowcasts_compare[[1]]),length(Nowcast_list)))%>%
  apply(c(2),c)%>%as_tibble()
colnames(SARI_nowcast_dmax)<-colnames(SARI_nowcasts_compare[[1]])   
SARI_nowcast_dmax<-SARI_nowcast_dmax%>%mutate(t=as.numeric(t),`2.5%`=as.numeric(`2.5%`),`50%`=as.numeric(`50%`),`97.5%`=as.numeric(`97.5%`),
                                              m=factor(m,levels=c('Arbovirus survivor','Arbovirus survivor linear'
                                              )),y=as.numeric(y),
                                              d=as.numeric(d),now=as.numeric(now))
SARI_nowcast_dmax_summary<-SARI_nowcast_dmax%>%group_by(m,d,s)%>%
  summarise(Coverage=mean(y>=`2.5%`&y<=`97.5%`),
            `Prediction interval width`=mean(`97.5%`-`2.5%`),
            `Mean absolute error`=mean(abs(`50%`-y)))%>%
  pivot_longer(cols=c('Coverage','Mean absolute error','Prediction interval width'),names_to = 'metric')%>%
  mutate(m=factor(m,levels=c('Arbovirus survivor','Arbovirus survivor linear'
  )),
  metric=factor(metric,levels=c('Mean absolute error','Prediction interval width','Coverage')))


#plot_data_SARI<-filter(SARI_nowcast_dmax_summary, d==0)%>%pivot_wider(names_from=c('m'))%>%filter(metric=="Mean Absolute Error") 
#ggplot(plot_data_SARI)+   geom_text(aes(x=`CLR`,y=`CLR Age`/`CLR`,label=s))+geom_hline(yintercept = 1)
SARI_nowcast_dmax_summary<-filter(SARI_nowcast_dmax_summary, metric=='Mean absolute error')

# Plot of prediction performance metrics.
mae_survivor_linear <- ggplot(SARI_nowcast_dmax_summary)+
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
mae_survivor_linear



# Compare parameter inference delta to prediction experiment delta.
# Plots:
library(geobr)
states<-read_state(year=2013)
state_guide<-tibble(s=states$code_state, s_abb=states$abbrev_state,s_full=states$name_state)
geom_guide<-tibble(s=states$abbrev_state,s_full=states$name_state,geom=states$geom)
n_sim<-dim(Dengue_surv_full[[1]]$chain1)[1]*2
Dengue_Survivor_linear_x<-list()
Dengue_combined_Survivor_linear_x<-list()
surv_linear_mu_x<-covid_surv_linear_quant_x<-Dengue_surv_linear_quant_x<-list()
Dengue_pi_samples_surv_age<-Dengue_mu_samples_surv_age<-Dengue_y_samples_surv_age<-Dengue_chi_surv_age<-Dengue_x_samples_surv_age<-list()
for(j in 1:S){
  # Survivor age model
  Dengue_Survivor_linear_x[[j]]<-as.mcmc.list(Dengue_surv_full[[j]])
  Dengue_combined_Survivor_linear_x[[j]] <- as_tibble(do.call('rbind',Dengue_Survivor_linear_x[[j]]))
}

# Model coefficient output
delta_slope<-zeta_spline<-eta_spline<-list()
for(j in 1:S){
  delta_slope[[j]]<-select(Dengue_combined_Survivor_linear_x[[j]],starts_with('delta'))%>%as.matrix()%>%array(dim=c(n_sim))%>%
    quantile(c(0.025,0.5,0.975))%>%melt(varnames=c('quantile'),value.name='y')%>%mutate(s=j,quantile=c(0.025,0.5,0.975))%>%
    spread(quantile,y)%>%mutate(s_abb=federal_guide$abbrev_state[s])
  
}
#delta matrix
delta_matrix<-matrix(unlist(delta_slope), byrow = TRUE, ncol = 5)%>%as_tibble()
colnames(delta_matrix)<-c('s','2.5%','50%','97.5%','s_abb')
delta_matrix$s_full<-factor(federal_guide$name_state[as.numeric(delta_matrix$s)],levels=(federal_guide$name_state[as.numeric(delta_matrix$s)]))
delta_matrix<-full_join(delta_matrix,state_guide[,2:3], by=c('s_abb','s_full'))
delta_matrix$s_full<-factor(delta_matrix$s_full, levels=delta_matrix$s_full)
delta_matrix<-delta_matrix%>%mutate(model='Parameter inference experiment')
delta_matrix[,2]<-as.numeric(unlist(delta_matrix[,2]))
delta_matrix[,3]<-as.numeric(unlist(delta_matrix[,3]))
delta_matrix[,4]<-as.numeric(unlist(delta_matrix[,4]))

# prediction delta:
n_sim<-dim(Dengue_surv_bin[[1]]$chain1)[1]*2
# Survivor age model
Dengue_Survivor_bin_x<-as.mcmc.list(Dengue_surv_bin[[4]])
Dengue_combined_Survivor_bin_x <- as_tibble(do.call('rbind',Dengue_Survivor_bin_x))
delta_slope_pred<-select(Dengue_combined_Survivor_bin_x,starts_with('delta'))%>%as.matrix()%>%array(dim=c(n_sim,S))%>%
  apply(c(2),quantile,c(0.025,0.5,0.975))%>%melt(varnames=c('quantile','s'),value.name='y')%>%
  spread(quantile,y)%>%mutate(Nowcast_date=as.numeric(paste(Nowcast_list[4])),s=federal_guide$abbrev_state[s])
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
       title='Arbovirus cases',
       subtitle='Linear effect of case load on reporting delay', legend=NULL)+
  theme_minimal()+
  theme(text = element_text(size=12), plot.title = element_text(size=16),plot.subtitle=element_text(size=14),
        legend.title = element_text(size=14),strip.text = element_text(size=14), 
        legend.position = "bottom")
delta_comp_plot
ggsave(delta_comp_plot,file='full_delta_Dengue_compare.pdf',width=9,height=4)



