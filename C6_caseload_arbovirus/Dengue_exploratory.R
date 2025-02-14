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
library(gridExtra)

# Reading in the table from Wikipedia
page = read_html("https://en.wikipedia.org/wiki/Federative_units_of_Brazil")
# Obtain the piece of the web page that corresponds to the "wikitable" node
my.table = html_node(page, ".wikitable")
# Convert the html table eleument into a data frame
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


Data_delay_gam<-as.tibble(Data_wide)
Data_wide[is.na(Data_wide)]=0
totals<-Data_delay_gam%>%mutate(total=as.numeric(unlist((apply(Data_wide[,3:34],1, sum, na.rm=TRUE)))))%>%mutate(s=paste(federal_guide$abbrev_state[s]))
Data_delay_gam$prop_0<-as.numeric(unlist((Data_wide[,3])))
Data_delay_gam$prop_1<-as.numeric(unlist((Data_wide[,3]+Data_wide[,4])))
Data_delay_gam$prop_2<-as.numeric(unlist((Data_wide[,3]+Data_wide[,4]+Data_wide[,5])))
Data_delay_gam$prop_3<-as.numeric(unlist((Data_wide[,3]+Data_wide[,4]+Data_wide[,5]+Data_wide[,6])))
Data_delay_gam$prop_4<-as.numeric(unlist((Data_wide[,3]+Data_wide[,4]+Data_wide[,5]+Data_wide[,6]+Data_wide[,7])))

Data_delay_long_gam<-as.tibble(Data_delay_gam%>%pivot_longer(!c(t,s), names_to = "delay", values_to = "count")%>%mutate(s=paste(federal_guide$abbrev_state[s])))
population<-my.table[-1,c(1,2,6)]
colnames(population)<-c("s_full","s","population")
population[,3]<-as.numeric(str_replace_all(as.character(unlist(population[,3])),",",""))
population[,3]<-as.numeric(str_replace_all(as.character(unlist(population[,3])),",",""))

population_order<-full_join(federal_guide[,-1],population%>%mutate(abbrev_state=s), by="abbrev_state")

Data_delay_pop_gam<-full_join(Data_delay_long_gam,population_order, by="s")
Data_delay_prop_gam<-Data_delay_pop_gam%>%filter(delay%in%c("prop_0","prop_1","prop_2","prop_3","prop_4"))
Data_delay_prop_gam<-full_join(Data_delay_prop_gam,totals[,c(1,2,dim(totals)[2])],by = join_by(t, s))
Data_delay_prop_gam_brazil<-Data_delay_prop_gam%>%group_by(t,delay)%>%summarise(total=sum(total,na.rm=TRUE),population=sum(population), count=sum(count,na.rm =TRUE))
Data_delay_prop_gam<-Data_delay_prop_gam%>%mutate(total.per.capita=total/population, cum_prop=rep(0:4,N*S))
Data_delay_prop_gam_brazil<-ungroup(Data_delay_prop_gam_brazil)%>%mutate(cum_prop=rep(0:4,N))

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

cumu_gam_all<-gam(data = Data_delay_prop_gam_brazil,probit((count-1)/total)~s(log(total+1),k=20,by=delay)+delay, family = gaussian)
cumu_gam_all_linear<-gam(data = Data_delay_prop_gam_brazil,probit((count-1)/total)~log(total+1)+delay, family = gaussian)
#plot(cumu_gam_all)
overall_smooth_trend<-tibble(y=cumu_gam_all$fitted.values,totals=Data_delay_prop_gam_brazil$total,delay=as.factor(Data_delay_prop_gam_brazil$cum_prop))
overall_linear_trend<-tibble(y=cumu_gam_all_linear$fitted.values,totals=Data_delay_prop_gam_brazil$total,delay=as.factor(Data_delay_prop_gam_brazil$cum_prop))

brazil_cumulative_totals<-ggplot()+
  geom_point(data=filter(Data_delay_prop_gam_brazil, cum_prop%in%c(0,2,4)),aes(x=log(total+1),y=probit((count-1)/total), colour=as.factor(cum_prop)),alpha=0.5)+#facet_wrap(~s_full,scales="free", ncol=3)+
  geom_line(data=filter(overall_smooth_trend, delay%in%c(0,2,4)),aes(x=log(totals+1),y=y, colour=delay),linetype='dashed')+
  geom_smooth(data=filter(Data_delay_prop_gam_brazil, cum_prop%in%c(0,2,4), total>exp(7)),aes(x=log(total+1),y=probit((count-1)/total), colour=as.factor(cum_prop)),level=0.95,method='lm')+
 # geom_line(data=overall_linear_trend,aes(x=log(totals),y=y, colour=delay))+
  scale_colour_viridis_d(name="Cumulative proportions reported at delay week:", begin=0.1, end=0.85)+
  labs(x="Total arbovirus cases reported (log transform)",y="Cumulative proportions reported (probit transform)",
       title = 'Arbovirus cases for Brazil')+theme_minimal()+
  theme(text = element_text(size=14), axis.title = element_text(size=16), plot.title = element_text(size=20),plot.subtitle=element_text(size=18),
        legend.title = element_text(size=16),strip.text = element_text(size=14),
        axis.text = element_text(size=9),legend.position="bottom")
brazil_cumulative_totals
ggsave(brazil_cumulative_totals,file='Plots/brazil_cumulative_totals_dengue_7.pdf',width=10,height=6)

Data_delay_prop_gam_brazil<-Data_delay_prop_gam_brazil%>%mutate(time_bin=paste(as.character((floor(t/100))*100+1),"weeks","≤ t <",as.character((floor(t/100)+1)*100+1),"weeks"))
brazil_cumulative_timebin<-ggplot()+
  geom_point(data=filter(Data_delay_prop_gam_brazil,cum_prop==0),aes(x=log(total+1),y=probit((count-1)/total), colour=time_bin),alpha=0.5)+#facet_wrap(~s_full,scales="free", ncol=3)+
  #geom_line(data=overall_smooth_trend,aes(x=log(totals+1),y=y, colour=delay),linetype='dashed')+
  geom_smooth(data=filter(Data_delay_prop_gam_brazil,cum_prop==0),aes(x=log(total+1),y=probit((count-1)/total), colour=time_bin),level=0.95,method='lm')+
  # geom_line(data=overall_linear_trend,aes(x=log(totals),y=y, colour=delay))+
  scale_colour_viridis_d(name="Data in time category:", begin=0.1, end=0.85, option = 1)+
  labs(x="Total arbovirus cases reported (log transform)",y="Cumulative proportions reported (probit transform)",
       title = 'Arbovirus cases for Brazil',
       subtitle = 'Proportions of cases reported the same week as occurrence')+theme_minimal()+
  theme(text = element_text(size=14), axis.title = element_text(size=16), plot.title = element_text(size=20),plot.subtitle=element_text(size=18),
        legend.title = element_text(size=14),strip.text = element_text(size=14),
        axis.text = element_text(size=9),legend.position="right")
brazil_cumulative_timebin
ggsave(brazil_cumulative_timebin,file='Plots/brazil_cumulative_timebin_dengue.png',width=12,height=6)


Data_delay_sum<-as.tibble(Data_wide)
totals<-Data_delay_sum%>%mutate(total=as.numeric(unlist((apply(Data_wide[,3:34],1, sum, na.rm=TRUE)))))%>%mutate(s=paste(federal_guide$abbrev_state[s]))
Data_delay_sum$prop_0<-as.numeric(unlist((Data_wide[,3])/apply(Data_wide[,3:34],1, sum, na.rm=TRUE)))
Data_delay_sum$prop_1<-as.numeric(unlist((Data_wide[,3]+Data_wide[,4])/apply(Data_wide[,3:34],1, sum, na.rm=TRUE)))
Data_delay_sum$prop_2<-as.numeric(unlist((Data_wide[,3]+Data_wide[,4]+Data_wide[,5])/apply(Data_wide[,3:34],1, sum, na.rm=TRUE)))
Data_delay_sum$prop_3<-as.numeric(unlist((Data_wide[,3]+Data_wide[,4]+Data_wide[,5]+Data_wide[,6])/apply(Data_wide[,3:34],1, sum, na.rm=TRUE)))
Data_delay_sum$prop_4<-as.numeric(unlist((Data_wide[,3]+Data_wide[,4]+Data_wide[,5]+Data_wide[,6]+Data_wide[,7])/apply(Data_wide[,3:34],1, sum, na.rm=TRUE)))

Data_delay_long<-as.tibble(Data_delay_sum%>%pivot_longer(!c(t,s), names_to = "delay", values_to = "count")%>%mutate(s=paste(federal_guide$abbrev_state[s])))
population<-my.table[,c(1,2,6)]
colnames(population)<-c("s_full","s","population")
population[,3]<-as.numeric(str_replace_all(as.character(unlist(my.table[,6])),",",""))
population_order<-full_join(federal_guide,population%>%mutate(abbrev_state=s), by="abbrev_state")

Data_delay_pop<-full_join(Data_delay_long,population, by="s")
Data_delay_prop<-Data_delay_pop%>%filter(delay%in%c("prop_0","prop_1","prop_2","prop_3","prop_4"))
Data_delay_prop<-full_join(Data_delay_prop,totals[,c(1,2,dim(totals)[2])],by = join_by(t, s))
Data_delay_prop<-Data_delay_prop%>%mutate(total.per.capita=total/population, cum_prop=rep(0:4,N*S))
Data_delay_prop$count[is.na(Data_delay_prop$count)]<-1
# plot cumulative proportions against totals
ggplot(Data_delay_prop,aes(x=total,y=probit(count), colour=delay))+
  geom_point(alpha=0.1)+facet_wrap(~s,scales="free")+geom_smooth(alpha=0.1)+scale_x_log10()+
  labs(x="Total counts",y="Cumulative proportions")


# plot cumulative proportions against totals
dengue_cumulative_exploratory_all<-ggplot(Data_delay_prop,aes(x=total,y=probit(count), colour=as.factor(cum_prop)))+
  geom_point(alpha=0.25)+facet_wrap(~s_full,scales="free", ncol=3)+geom_smooth(alpha=0.1)+scale_x_log10()+
  scale_colour_viridis_d(name="Cumulative proportions reported at delay week:", begin=0.1, end=0.85)+
  labs(x="Eventual total dengue hospitaltions reported (log scale)",y="Cumulative proportions of Arbovirus cases reported (probit scale)",
       title = 'Arbovirus cases cumulative proportions reported')+theme_minimal()+
  theme(legend.position = 'bottom')+
  theme(text = element_text(size=12), axis.title = element_text(size=14), plot.title = element_text(size=16),plot.subtitle=element_text(size=14),
        legend.title = element_text(size=14),strip.text = element_text(size=12),legend.position="bottom")
ggsave(dengue_cumulative_exploratory_all,file='dengue_cumulative_exploratory_all.pdf',width=10,height=16.5)


dengue_cumulative_exploratory_3<-ggplot(filter(Data_delay_prop,cum_prop%in%c(0,2,4)),aes(x=total+1,y=probit(count), colour=as.factor(cum_prop)))+
  geom_point(alpha=0.25)+facet_wrap(~s_full,scales="free", ncol=3)+
  geom_smooth(alpha=0.1)+scale_x_log10()+
  scale_colour_viridis_d(name="Cumulative proportions reported at delay week:", begin=0.1, end=0.85)+
  labs(x="Eventual total dengue hospitaltions reported (log scale)",y="Cumulative proportions of Arbovirus cases reported (probit scale)",
       title = 'Arbovirus cases cumulative proportions reported')+theme_minimal()+
  theme(legend.position = 'bottom')+
  theme(text = element_text(size=12), axis.title = element_text(size=14), plot.title = element_text(size=16),plot.subtitle=element_text(size=14),
        legend.title = element_text(size=14),strip.text = element_text(size=12),legend.position="bottom")
ggsave(dengue_cumulative_exploratory_3,file='dengue_cumulative_exploratory_3.pdf',width=10,height=16.5)

dengue_cumulative_per_capita<-ggplot(filter(Data_delay_prop,cum_prop%in%c(0,2,4)),aes(x=total.per.capita,y=probit(count), colour=as.factor(cum_prop)))+
  geom_point(alpha=0.25)+facet_wrap(~s_full,scales="free", ncol=3)+
  geom_smooth(alpha=0.1)+scale_x_log10()+
  scale_colour_viridis_d(name="Cumulative proportions reported at delay week:", begin=0.1, end=0.85)+
  labs(x="Eventual total dengue hospitaltions reported per capita (log scale)",y="Cumulative proportions of Arbovirus cases reported (probit scale)",
       title = 'Arbovirus cases cumulative proportions reported')+theme_minimal()+
  theme(legend.position = 'bottom')+
  theme(text = element_text(size=12), axis.title = element_text(size=14), plot.title = element_text(size=16),plot.subtitle=element_text(size=14),
        legend.title = element_text(size=14),strip.text = element_text(size=12),legend.position="bottom")

ggsave(dengue_cumulative_per_capita,file='dengue_cumulative_per_capita.pdf',width=10,height=16.5)

Data_delay_prop_total<-Data_delay_prop%>%mutate(cum_counts=count*total)%>%group_by(t,cum_prop)%>%
  summarise(total=sum(total),cum_counts=sum(cum_counts))%>%mutate(count=cum_counts/total)

plot1<-ggplot(filter(Data_delay_prop_total,cum_prop%in%c(0,2,4)),aes(x=dates_all[t],y=total))+
  geom_line()+#facet_wrap(~s_full,scales="free", ncol=3)+
  #geom_smooth(alpha=0.1)+scale_y_log10()+
  labs(x="Date of hospitalisation",y="Arbovirus cases",
       subtitle = 'Total hospitalisations',title='Arbovirus cases for Brazil')+theme_minimal()+
  theme(legend.position = 'bottom')+
  scale_x_date(date_breaks = "2 month", date_labels =  "%b %Y")+
  theme(text = element_text(size=12), axis.title = element_text(size=14), plot.title = element_text(size=16),plot.subtitle=element_text(size=14),
        legend.title = element_text(size=14),strip.text = element_text(size=12),legend.position="bottom")

plot2<-ggplot(filter(Data_delay_prop_total,cum_prop%in%c(0,2,4)),aes(x=dates_all[t],y=probit(count), colour=as.factor(cum_prop)))+
  geom_point(alpha=0.5)+#facet_wrap(~s_full,scales="free", ncol=3)+
  geom_smooth(alpha=0.1)+#scale_x_log10()+
  scale_colour_discrete(name="Cumulative proportions reported at delay week:")+
  labs(x="Date of hospitalisation",y="Cumulative proportions (probit scale)",
       subtitle = 'Cumulative proportions of Arbovirus cases reported')+theme_minimal()+
  theme(legend.position = 'bottom')+
  scale_x_date(date_breaks = "2 month", date_labels =  "%b %Y")+
  theme(text = element_text(size=12), axis.title = element_text(size=14), plot.title = element_text(size=16),plot.subtitle=element_text(size=14),
        legend.title = element_text(size=14),strip.text = element_text(size=12),legend.position="bottom")

Dengue_time_Plots<-grid.arrange(plot1,plot2,ncol=1)
ggsave(Dengue_time_Plots,file='Dengue_time_plots.pdf',width=9,height=8)


dengue_brazil_plot<-ggplot()+
  scale_colour_viridis_d(name="Cumulative proportions reported at delay week:", begin=0.1, end=0.85)+
  geom_point(data=Data_delay_prop_total,aes(x=log(total),y=(count), colour=as.factor(cum_prop)),alpha=0.5)+#facet_wrap(~s_full,scales="free", ncol=3)+
  geom_smooth(data=Data_delay_prop_total,alpha=0.1,aes(x=log(total),y=(count), colour=as.factor(cum_prop)))+
  labs(x="Total Arbovirus cases",y="Cumulative proportion reported",
       subtitle = 'Arbovirus cases in Brazil')+theme_minimal()+
  theme(legend.position = 'bottom')+
  theme(text = element_text(size=12), axis.title = element_text(size=14), plot.title = element_text(size=16),plot.subtitle=element_text(size=14),
        legend.title = element_text(size=14),strip.text = element_text(size=12),legend.position="bottom")
ggsave(dengue_brazil_plot,file='dengue_brazil_plot.pdf',width=9,height=5)



