# plot((seq(from=-5,to=5,length=100)),iprobit(seq(from=-5,to=5,length=100)))
# 
# plot((seq(from=-5,to=5,length=100)),icloglog(seq(from=-5,to=5,length=100)))
# 
# plot((seq(from=-5,to=5,length=100)),expit(seq(from=-5,to=5,length=100)))



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
censored_z_melt_cum<-censored_z_brazil%>%melt(value.name = 'z', varnames=c('t','d'))%>%mutate(d=as.factor(d-1))%>% 
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
         
z_plot<-ggplot()+
  geom_point(data=filter(y_all_brazil,t<101),aes(x=as.Date(first_date)+(t)*7,y=y))+ 
  theme_minimal()+labs(y="Cumulative count reported", x="Date of hospitalisation", title="Severe acute respiratory illness (SARI) hospitalisations",subtitle = "Brazil")+#+geom_vline(xintercept = 51)
  geom_line(data=censored_z_melt_cum, aes(x=as.Date(first_date)+(t)*7, y=CumulativeSum, colour=d))+#, linetype = 'dashed')+
  scale_colour_discrete(name='Delay (weeks)')+
  scale_x_date(limits = c(as.Date("2022-04-15"),as.Date("2022-12-25")),breaks = "1 month", date_labels =  "%b %Y")+
  scale_y_continuous(limits = c(0,15000))+ 
  theme(text = element_text(size=12), axis.title = element_text(size=14), plot.title = element_text(size=16),plot.subtitle=element_text(size=14),
                                                 legend.title = element_text(size=14),strip.text = element_text(size=12),legend.position="bottom")+
  guides(color = guide_legend(nrow = 2, byrow=TRUE))
z_plot
ggsave(z_plot,file='SARI_z_plot.pdf',width=9,height=5)
  
              
y_plot<-ggplot()+
  guides(fill="none", colour="none")+
  scale_linetype(name="Counts")+
  geom_point(data=filter(y_all_brazil,t<101),aes(x=as.Date(first_date)+(t)*7,y=y))+ 
  theme_minimal()+labs(y="Hospitalisations", x="Date of hospitalisation", title="Severe acute respiratory illness (SARI) hospitalisations",subtitle = "Brazil")+#+geom_vline(xintercept = 51)
   geom_line(data=censored_y_brazil, aes(x=as.Date(first_date)+(t)*7, y=y, linetype='Reported'))+#, linetype = 'dashed')+
  scale_colour_discrete(name='Nowcast Date')+scale_fill_discrete(name='Nowcast Date')+
  geom_line(data=y_quantiles, aes(x=as.Date(first_date)+(t)*7, y=`50%`, colour=Nowcast_date, linetype='Predicted'))+
   scale_x_date(limits = c(as.Date("2022-04-15"),as.Date("2022-12-25")),breaks = "1 month", date_labels =  "%b %Y")+
  scale_y_continuous(limits = c(0,15000))+ theme(text = element_text(size=12), axis.title = element_text(size=14), plot.title = element_text(size=16),plot.subtitle=element_text(size=14),
                                                 legend.title = element_text(size=14),strip.text = element_text(size=12),legend.position="bottom")+
  geom_ribbon(data=y_quantiles, aes(x=as.Date(first_date)+(t)*7, ymin=`2.5%`, ymax=`97.5%`,fill=Nowcast_date),alpha=0.5)
y_plot
ggsave(y_plot,file='SARI_y_plot_seasonal.pdf',width=9,height=4.5)

