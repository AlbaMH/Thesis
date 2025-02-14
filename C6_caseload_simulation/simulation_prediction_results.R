
# SIMULATION PERFORMANCE METRICS #
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
library(compositions)
library(MASS)
library(dplyr)
library(ggh4x)
library(grid)
library(gridExtra)
# Simulated lambda plots

load("~/media/alba/Disk 1/OneDrive/PhD/Case load simulations/Simulation studies/new data/Simulation_data_y.RData")

sims=100
# Number of weeks in the time series.
N <- 100
# Number of delays.
D <- 7
# Number of scenarios
S <- 27
# Simulate GDM lambda.
# # Create temporal trend.
# t_poly <- poly(1:N,8)
# sim_zeta <- matrix(nrow=N, ncol=S)
# temp<-matrix(nrow=N, ncol=S)
# sim_zeta_coefficients <- matrix(c(c(4,0,-2,0,0,0,2,0.5),c(0,0,-2,0,0,0,2,0.5),c(-4,0,-2,0,0,0,2,0.5)),ncol=S)#,)
# sim_zeta <- t_poly%*%sim_zeta_coefficients
# sim_zeta_data <- melt(sim_zeta,varnames = c('t','s'),value.name='zeta')
# # Simulate mean temporal trend.
# sim_lambda <- matrix(nrow=N,ncol=S)
# sim_iota<- rep(5,9) # Spatial intercept of mean reported cases - iota 
# for(s in 1:S){
#   sim_lambda[,s]<-exp(sim_iota[s]+sim_zeta[,s])
# }
# 
# sim_lambda_data <- melt(sim_lambda,varnames=c('t','s'),value.name='lambda')%>%
#   mutate(s=as.factor(s))
# zeta_ind<-c("-0.2","0","0.2")
# 
# # Plot of mean reported cases:
# lambda_all <- ggplot(sim_lambda_data)+geom_line(aes(x=t,y=lambda, colour=zeta_ind[s],linetype=zeta_ind[s]))+
#   theme_minimal()+
#   labs(title = "Simulated mean of the total cases",
#        subtitle="For three scenarios of linear temporal trends in the total cases ",x="Time",y=expression(lambda['t']))+
#   scale_colour_viridis_d(name=expression(zeta),begin=0.1,end=0.9)+
#   scale_linetype(name=expression(zeta))+
#   theme_minimal()+
#   theme(text = element_text(size=12), plot.title = element_text(size=16),plot.subtitle=element_text(size=14),
#         legend.title = element_text(size=14),strip.text = element_text(size=14), legend.position = "right")
# lambda_all
# ggsave(lambda_all,file='Plots/simulated_lambda.pdf',width=9,height=4)

# 
# blank_data=tibble(y=rnorm(N,0,1),t=1:N)
# blank_jagam=jagam(y~s(t,k=floor(N/5), bs='tp'),
#                   data=blank_data,file='blank.jags',
#                   knots=list(t=seq(1,N,length=floor(N/5))))
# 
# 
# K_t=dim(blank_jagam$jags.data$S1)[1]
# S_t=blank_jagam$jags.data$S1
# X_t=blank_jagam$jags.data$X[,2:(K_t+1)]
# 
# # simulate means of total cases
# set.seed(24)
# zeta_value <- c(rep(-0.1,3),rep(0,3),rep(0.1,3)) # Chosen temporal trend in totals. 
# sigma_alpha_sim <- c(5,5)
# Omega_alpha_sim<-array(dim=c(K_t,K_t,S))
# kappa_alpha_sim<-matrix(nrow=K_t,ncol=S)
# sim_alpha_raw<-sim_alpha<-sim_lambda<-matrix(nrow=N,ncol=S)
# alpha_linear_log_trend<-numeric(S)
# # Spatial intercept of mean reported cases - iota 
# sim_iota<- rnorm(S,5,sd=0.25)
# for(s in 1:S){
#   # Generate thin plate spline:
#   Omega_alpha_sim[,,s] <- S_t[1:(K_t),1:(K_t)]/sigma_alpha_sim[1]^2 + S_t[1:(K_t),(K_t+1):(2*K_t)]/sigma_alpha_sim[2]^2
#   kappa_alpha_sim[,s] <- rmnorm_chol(n=1,mean=rep(0,K_t), chol(Omega_alpha_sim[1:(K_t),1:(K_t),s]))
#   sim_alpha_raw[,s]<- X_t[1:N,1:(K_t)]%*%kappa_alpha_sim[1:(K_t),s]
#   # Scale spline to have mean zero and no linear trend.
#   alpha_linear_log_trend[s]<-round(as.numeric(lm(data=tibble(y=sim_alpha_raw[,s],x=(1:N)-mean(1:N)))$coefficients[2]),3)
#   sim_alpha[,s]<-sim_alpha_raw[,s] - ((1:N)-mean(1:N))*alpha_linear_log_trend[s]
#   # Add in chosen linear trend in the totals.
#   for(t in 1:N){
#     sim_lambda[t,s]<-exp(sim_iota[s]+sim_alpha[t,s]+zeta_value[s]*(t-mean(1:N))/sd(1:N))
#   }
# }
# 
# plot(sim_lambda[,1])
# 
# sim_lambda_data <- melt(sim_lambda,varnames=c('t','s'),value.name='lambda')%>%
#   mutate(s=as.factor(s),zeta=as.factor(zeta_value[s]),line=as.factor(rep(c(rep(1,N),rep(2,N),rep(3,N)),3)))
# # Plot of mean reported cases:
# lambda_sim_plot <- ggplot(sim_lambda_data)+geom_line(aes(x=t,y=lambda, colour=zeta, linetype=line))+
#   theme_minimal()+
#   facet_wrap(~zeta, scales='free')+
#   labs(title = "Simulation examples of the mean total cases temporal trend",
#        subtitle="For three different overall linear temporal trends in the total cases", x="Time (t)",y=expression(lambda['t,s']))+
#   scale_colour_viridis_d(name='Linear temporal trend in total cases (zeta)',begin=0.1,end=0.8)+
#   scale_linetype(name='Simulation example')+
#   theme_minimal()+
#   theme(text = element_text(size=12), plot.title = element_text(size=16),plot.subtitle=element_text(size=14),
#         legend.title = element_text(size=12),strip.text = element_text(size=14), legend.position = "bottom")
# lambda_sim_plot
# ggsave(lambda_sim_plot,file='Plots/lambda_sim_plot.pdf',width=9,height=4)       

# Simulated data from simulation_data.R script:
sim_y<-list()
sim_lambda<-list()
sim_p<-list()
for(i in 1:sims){
sim_y[[i]]<-matrix(simulation_data_list[[i]]$y,nrow=N,ncol=S)%>%melt(varnames=c('t','s'),value.name='y')%>%
  mutate(sim=i)
sim_lambda[[i]]<-matrix(simulation_data_list[[i]]$lambda,nrow=N,ncol=S)%>%melt(varnames=c('t','s'),value.name='lambda')%>%
  mutate(sim=i)
sim_p[[i]]<-array(simulation_data_list[[i]]$p,dim=dim(simulation_data_list[[i]]$p))%>%melt(varnames=c('t','s','d'),value.name='p')%>%
  mutate(sim=i)
}
# Table of simulated total counts (y)
sim_y_long<-unlist(sim_y)%>%array(dim=c(dim(sim_y[[1]]),sims))%>%apply(2,c)
colnames(sim_y_long)<- colnames(sim_y[[1]])
sim_y_long<-as_tibble(sim_y_long)%>%
  mutate(eta=simulation_data_list[[1]]$eta[s],zeta=simulation_data_list[[1]]$zeta[s],delta=simulation_data_list[[1]]$delta[s])%>%
  mutate(sim=as.factor(sim),zeta=factor(zeta,levels=unique(zeta)),
         eta=factor(eta,levels=unique(eta)),delta=factor(delta,levels=unique(delta)))
# Table of simulated mean total counts (lambda).
sim_lambda_long<-unlist(sim_lambda)%>%array(dim=c(dim(sim_lambda[[1]]),sims))%>%apply(2,c)
colnames(sim_lambda_long)<- colnames(sim_lambda[[1]])
sim_lambda_long<-as_tibble(sim_lambda_long)%>%
  mutate(eta=simulation_data_list[[1]]$eta[s],zeta=simulation_data_list[[1]]$zeta[s],delta=simulation_data_list[[1]]$delta[s])%>%
  mutate(s=as.factor(s),sim=as.factor(sim),zeta=factor(zeta,levels=unique(zeta)),
          eta=factor(eta,levels=unique(eta)),delta=factor(delta,levels=unique(delta)))
# Table of simulated cumulative proportions (p).
sim_p_long<-unlist(sim_p)%>%array(dim=c(dim(sim_p[[1]]),sims))%>%apply(2,c)
colnames(sim_p_long)<- colnames(sim_p[[1]])
sim_p_long<-as_tibble(sim_p_long)%>%mutate(eta=simulation_data_list[[1]]$eta[s],zeta=simulation_data_list[[1]]$zeta[s],delta=simulation_data_list[[1]]$delta[s])%>%
  mutate(s=as.factor(s),sim=as.factor(sim),zeta=factor(zeta,levels=unique(zeta)),
         eta=factor(eta,levels=unique(eta)),delta=factor(delta,levels=unique(delta)))

# Plot of mean reported cases:
lambda_sim_plot <- ggplot(filter(sim_lambda_long, sim%in%c(1:11),s%in%c(1,10,19)))+geom_line(aes(x=t,y=lambda, colour=zeta, linetype=(sim)))+
  theme_minimal()+
  facet_wrap(~zeta, scales='free')+
  labs(title = "Simulation examples of the mean total cases temporal trend",
       subtitle="For three different overall linear temporal trends in the total cases", x="Time (t)",y=expression(lambda['t,s']))+
  scale_colour_viridis_d(name='Linear temporal trend in total cases (zeta)',begin=0.1,end=0.8)+
  scale_linetype(name='Simulation')+
  theme_minimal()+
  theme(text = element_text(size=12), plot.title = element_text(size=16),plot.subtitle=element_text(size=14),
        legend.title = element_text(size=12),strip.text = element_text(size=14), legend.position = "bottom")
lambda_sim_plot
#ggsave(lambda_sim_plot,file='Plots/lambda_sim_plot.pdf',width=9,height=4)       

# Plot of cumulative proportions
p_sim_plot <- ggplot(filter(sim_p_long, sim%in%c(1),s%in%c(10:18)))+geom_line(aes(x=t,y=p, colour=as.factor(d)))+
  theme_minimal()+
  facet_grid(paste('Case load effect:',delta)~paste('Delay temporal trend:',eta))+
  labs(title = "Simulated cumulative proportions reported",
       subtitle="For three values of case load effect and three values of temporal trend in the delay" ,x="Time (t)",
       y=expression(S['t,d']))+
  scale_colour_viridis_d(name='Delay',begin=0.1,end=0.8)+
  theme_minimal()+
  scale_y_continuous(labels = scales::label_percent(), limits = c(0,1))+
  theme(text = element_text(size=12), plot.title = element_text(size=16),plot.subtitle=element_text(size=14),
        legend.title = element_text(size=12),strip.text = element_text(size=14))
p_sim_plot
#ggsave(p_sim_plot,file='Plots/p_sim_plot.pdf',width=9,height=7)       


load("~/media/alba/Disk 1/OneDrive/PhD/Case load simulations/Simulation studies/new data/matrix_simulations_prediction_nolink.RData") 
load("~/media/alba/Disk 1/OneDrive/PhD/Case load simulations/Simulation studies/new data/matrix_simulations_prediction_linear_log.RData")
load("~/media/alba/Disk 1/OneDrive/PhD/Thesis code/Chapter 6 Case Load Effect/Section 6_3 Simulation Experiments/Prediction Experiment/matrix_simulations_prediction_fixed_delta.RData")

#delta_slope_combined<-rbind(delta_slope_all, delta_slope_all_linear_log)
eta_slope_combined<-rbind(mutate(eta_slope_all,model='No incidence-delay'), mutate(eta_slope_all_linear_log,model='Incidence-delay'))
zeta_slope_combined<-rbind(mutate(zeta_slope_all,model='No incidence-delay'), mutate(zeta_slope_all_linear_log,model='Incidence-delay'))

# Simulation performance metrics 
# delta_metric<-delta_slope_combined%>%group_by(s,zeta,eta,delta)%>%
#   summarise(Coverage=mean(actual>=`2.5%`&actual<=`97.5%`),
#             `Prediction interval width`=mean(`97.5%`-`2.5%`),
#             `Mean absolute percentage error`=mean(abs((`50%`-actual)/actual)),
#             `Bias`=mean(`50%`-actual))%>%
#   pivot_longer(cols=c('Coverage','Mean absolute percentage error','Prediction interval width','Bias'),names_to = 'metric')

eta_metric<-eta_slope_combined%>%group_by(s,zeta,eta,delta,model)%>%
  summarise(Coverage=mean(actual>=`2.5%`&actual<=`97.5%`),
            `Prediction interval width`=mean(`97.5%`-`2.5%`),
            `Mean absolute percentage error`=mean(abs((`50%`-actual)/actual)),
            `Bias`=mean(`50%`-actual))%>%
  pivot_longer(cols=c('Coverage','Mean absolute percentage error','Prediction interval width','Bias'),names_to = 'metric')

zeta_metric<-zeta_slope_combined%>%group_by(s,zeta,eta,delta,model)%>%
  summarise(Coverage=mean(actual>=`2.5%`&actual<=`97.5%`),
            `Prediction interval width`=mean(`97.5%`-`2.5%`),
            `Mean absolute percentage error`=mean(abs((`50%`-actual)/actual)),
            `Bias`=mean(`50%`-actual))%>%
  pivot_longer(cols=c('Coverage','Mean absolute percentage error','Prediction interval width','Bias'),names_to = 'metric')


# Custom colours for tables.
lightGreen = "#9BE2B7"
lightRed = "#F19999"
lightBlue="#9CD4EF"
darkGreen = "#32B767"
darkRed = "#DD5A5A"
darkBlue="#1491CD"
# width=930, height=920
red <- "#E6194B"
green <- "#2AA198"
blue <- "#4363D8"
yellow <- "#FFE119"
gray <- "#A9A9A9"

# Nowcasts from Survivor model.
n_sim<-dim(y_quant[[1]])[1]/S
y_nolink<-y_quant%>%unlist()%>%array(dim=c(dim(y_quant[[1]]),sims))
colnames(y_nolink)<- colnames(y_quant[[1]])
y_nolink<-apply(y_nolink,2,c)%>%as_tibble()%>%mutate(`2.5%`=as.numeric(`2.5%`),
                                                     `50%`=as.numeric(`50%`), 
                                                     `97.5%`=as.numeric(`97.5%`),
                                                     t=as.numeric(t),
                                                     s=as.numeric(s))%>%mutate(model='Survivor')
# Nowcast from Survivor linear model.
y_linear_log<-y_quant_linear_log%>%unlist()%>%array(dim=c(dim(y_quant_linear_log[[1]]),sims))
colnames(y_linear_log)<- colnames(y_quant_linear_log[[1]])
y_linear_log<-apply(y_linear_log,2,c)%>%as_tibble()%>%mutate(`2.5%`=as.numeric(`2.5%`),
                                                          `50%`=as.numeric(`50%`), 
                                                          `97.5%`=as.numeric(`97.5%`),t=as.numeric(t),
                                                          s=as.numeric(s))%>%mutate(model='Survivor incidence-delay')

# Nowcasts from both models:
y_all<-rbind(y_nolink,y_linear_log)
# Add true y values:
y_all_sim<-right_join(sim_y_long,y_all, by=c('t','s','eta','zeta','delta','sim'))

ggplot()+geom_line(data=filter(y_all_sim,sim==5,delta==c(-0.2)),aes(x=t,y=y))+
  geom_line(data=filter(y_all_sim,sim==5,delta==c(-0.2)),aes(x=t,y=`50%`,colour=model))+
  facet_wrap(zeta~eta)+
  geom_ribbon(data=filter(y_all_sim,sim==5,delta==c(-0.2)),aes(x=t,ymin=`2.5%`,ymax=`97.5%`,fill=model),alpha=0.5)

# Filter for dates that are censored:
y_all_d<-y_all_sim%>%mutate(d=t-N)%>%filter(d>-D)

y_all_summary<-y_all_d%>%mutate(model=factor(model,levels=c('Survivor','Survivor incidence-delay'
                                             )),y=as.numeric(y),
                                             d=as.numeric(d),sim=as.numeric(sim))
y_all_summary_metrics<-y_all_summary%>%group_by(model,d,s,delta)%>%
  summarise(Coverage=mean(y>=`2.5%`&y<=`97.5%`),
            `Prediction interval width`=mean(`97.5%`-`2.5%`),
            `Mean absolute error`=mean(abs(`50%`-y)))%>%
  pivot_longer(cols=c('Coverage','Mean absolute error','Prediction interval width'),names_to = 'metric')%>%
  mutate(model=factor(model,levels=c('Survivor','Survivor incidence-delay'
  )),
  metric=factor(metric,levels=c('Mean absolute error','Prediction interval width','Coverage')))


#plot_data_SARI<-filter(SARI_nowcast_dmax_summary, d==0)%>%pivot_wider(names_from=c('m'))%>%filter(metric=="Mean Absolute Error") 
#ggplot(plot_data_SARI)+   geom_text(aes(x=`CLR`,y=`CLR Age`/`CLR`,label=s))+geom_hline(yintercept = 1)
y_all_summary_metrics<-filter(y_all_summary_metrics,delta!='0')
y_all_summary_average<-y_all_summary_metrics%>%group_by(model,d,metric)%>%
  summarise(value.mean=mean(as.numeric(value)), value.median=median(as.numeric(value)))

y_all_summary_average$value.mean[y_all_summary_average$metric%in%c("Prediction Interval Width","Mean Absolute Error")]<-NA
y_all_summary_average$value.median[y_all_summary_average$metric%in%c("Coverage")]<-NA
y_all_summary_average<-y_all_summary_average%>%group_by(model,d,metric)%>%mutate(value=sum(value.mean,value.median,na.rm=TRUE))
# Plot of prediction performance metrics.
compare_survivor_non_zero <- ggplot(y_all_summary_average)+
  geom_line(aes(x=d,y=value,colour=model),alpha=0.6)+
  geom_point(aes(x=d,y=value,colour=model,shape=model),alpha=0.6)+
  facet_grid2(~metric, scales = c("free"),independent = c("all"))+
  scale_colour_manual(name=NULL,values = c( "#4363D8", "#E6194B"))+
  scale_shape_manual(name=NULL,values=c(15,3))+
  scale_x_continuous()+
  labs(x='Prediction time difference',y=NULL,
       title="Simulation experiment",
       subtitle='Performance metrics from the prediction experiment')+
  theme_minimal()+
  theme(legend.position = "bottom",text = element_text(size=12), plot.title = element_text(size=16),
        plot.subtitle=element_text(size=14), legend.title = element_text(size=14),
        strip.text = element_text(size=14))
compare_survivor_non_zero
ggsave(compare_survivor_non_zero, file="Plots/Survivor_prediction_compare.pdf", width=9, height=4)
# 

# compare bias for positive and negative delta
y_all_summary_metrics_delta<-y_all_summary%>%group_by(model,d,s,delta)%>%
  summarise(Coverage=mean(y>=`2.5%`&y<=`97.5%`),
            `Prediction interval width`=mean(`97.5%`-`2.5%`),
            `Mean absolute error`=mean(abs(`50%`-y)),
          #  `Bias`=mean(`50%`-y)
            )%>%
  pivot_longer(cols=c('Coverage','Prediction interval width','Mean absolute error'),names_to = 'metric')%>%
  mutate(model=factor(model,levels=c('Survivor','Survivor incidence-delay'
  )),
  metric=factor(metric,levels=c('Mean absolute error','Prediction interval width','Coverage')))

y_all_summary_average_delta<-y_all_summary_metrics_delta%>%group_by(model,d,metric,delta)%>%
  summarise(value.mean=mean(as.numeric(value)), value.median=median(as.numeric(value)))

y_all_summary_average_delta$value.mean[y_all_summary_average_delta$metric%in%c("Prediction interval width","Mean absolute error")]<-NA
y_all_summary_average_delta$value.median[y_all_summary_average_delta$metric%in%c("Coverage")]<-NA
y_all_summary_average_delta<-y_all_summary_average_delta%>%group_by(model,d,metric,delta)%>%mutate(value=sum(value.mean,value.median,na.rm=TRUE))


# Plot of prediction performance metrics.
compare_delta_mae <- ggplot(filter(y_all_summary_average_delta,metric=="Mean absolute error"))+
  geom_line(aes(x=d,y=value,colour=model),alpha=0.6)+
  geom_point(aes(x=d,y=value,colour=model,shape=model),alpha=0.6)+
  facet_grid2(~paste('Case load effect:',delta),switch=c('y'))+
  scale_colour_manual(name=NULL,values = c( "#4363D8", "#E6194B"))+
  scale_shape_manual(name=NULL,values=c(15,3))+
  scale_x_continuous()+
  labs(x=NULL,y="Mean absolute error",
       title="Simulation experiment",
       subtitle='Performance metrics from the prediction experiment')+
  theme_minimal()+guides(colour="none", shape='none')+
  theme(legend.position = "bottom",text = element_text(size=14), axis.title = element_text(size=16), plot.title = element_text(size=20),plot.subtitle=element_text(size=18),
        legend.title = element_text(size=16),strip.text = element_text(size=14))

compare_delta_piw <- ggplot(filter(y_all_summary_average_delta,metric=="Prediction interval width"))+
  geom_line(aes(x=d,y=value,colour=model),alpha=0.6)+
  geom_point(aes(x=d,y=value,colour=model,shape=model),alpha=0.6)+
  facet_grid2(~paste('Case load effect:',delta))+
  scale_colour_manual(name=NULL,values = c( "#4363D8", "#E6194B"))+
  scale_shape_manual(name=NULL,values=c(15,3))+
  scale_x_continuous()+
  labs(x=NULL,y="Prediction interval width",
       title=NULL,
       subtitle=NULL)+
  theme_minimal()+guides(colour="none", shape='none')+
  theme(legend.position = "bottom",text = element_text(size=14), axis.title = element_text(size=16), plot.title = element_text(size=20),plot.subtitle=element_text(size=18),
        legend.title = element_text(size=16),strip.text = element_text(size=14))

compare_delta_cov <- ggplot(filter(y_all_summary_average_delta,metric=="Coverage"))+
  geom_line(aes(x=d,y=value,colour=model),alpha=0.6)+
  geom_point(aes(x=d,y=value,colour=model,shape=model),alpha=0.6)+
  facet_grid2(~paste('Case load effect:',delta))+
  scale_colour_manual(name=NULL,values = c( "#4363D8", "#E6194B"))+
  scale_shape_manual(name=NULL,values=c(15,3))+
  scale_x_continuous()+
  labs(x='Prediction time difference',y="Coverage",
       title=NULL,
       subtitle=NULL)+
  theme_minimal()+ 
  theme(legend.position = "bottom",text = element_text(size=14), axis.title = element_text(size=16), plot.title = element_text(size=20),plot.subtitle=element_text(size=18),
        legend.title = element_text(size=16),strip.text = element_text(size=14))

compare_delta_sep<-grid.arrange(compare_delta_mae,compare_delta_piw,compare_delta_cov,nrow=3,heights=c(4.1,3.2,4.1))


ggsave(compare_delta_sep, file="Plots/Survivor_delta_pred_sep.pdf", width=9, height=11)

# separate by eta and zeta 
# negative delta
y_all_summary_metrics_neg<-filter(y_all_summary,delta==(-0.2))%>%group_by(model,d,s,eta,zeta)%>%
  summarise(Coverage=mean(y>=`2.5%`&y<=`97.5%`),
            `Prediction interval width`=mean(`97.5%`-`2.5%`),
            `Mean absolute error`=mean(abs(`50%`-y)),
            `Bias`=mean(`50%`-y))%>%
  pivot_longer(cols=c('Coverage','Prediction interval width','Mean absolute error','Bias'),names_to = 'metric')%>%
  mutate(model=factor(model,levels=c('Survivor','Survivor incidence-delay'
  )),
  metric=factor(metric,levels=c('Mean absolute error','Prediction interval width','Coverage','Bias')))

y_all_summary_average_neg<-y_all_summary_metrics_neg%>%group_by(model,d,metric,eta,zeta)%>%
  summarise(value.mean=mean(as.numeric(value)), value.median=median(as.numeric(value)))

y_all_summary_average_neg$value.mean[y_all_summary_average_neg$metric%in%c("Prediction interval width","Mean absolute error","Bias")]<-NA
y_all_summary_average_neg$value.median[y_all_summary_average_neg$metric%in%c("Coverage")]<-NA
y_all_summary_average_neg<-y_all_summary_average_neg%>%group_by(model,d,metric,eta,zeta)%>%mutate(value=sum(value.mean,value.median,na.rm=TRUE))
compare_delta_mae_neg <- ggplot(filter(y_all_summary_average_neg,metric=="Mean absolute error"))+
  geom_line(aes(x=d,y=value,colour=model),alpha=0.6)+
  geom_point(aes(x=d,y=value,colour=model,shape=model),alpha=0.6)+
  facet_grid2(paste('Delay temporal trend:',eta)~paste('Totals temporal trend:',zeta))+
  scale_colour_manual(name=NULL,values = c( "#4363D8", "#E6194B"))+
  scale_shape_manual(name=NULL,values=c(15,3))+
  scale_x_continuous()+
  labs(x='Prediction time difference',y="Mean absolute error",
       title="Simulation experiment",
       subtitle='Mean absolute error for negative case load relationships (delta=-0.2)')+
  theme_minimal()+
  theme(legend.position = "bottom",text = element_text(size=14), axis.title = element_text(size=16), plot.title = element_text(size=20),plot.subtitle=element_text(size=18),
        legend.title = element_text(size=16),strip.text = element_text(size=14))
compare_delta_mae_neg
ggsave(compare_delta_mae_neg, file="Plots/neg_delta_mae.pdf", width=9, height=10)

compare_delta_piw_neg <- ggplot(filter(y_all_summary_average_neg,metric=="Prediction interval width"))+
  geom_line(aes(x=d,y=value,colour=model),alpha=0.6)+
  geom_point(aes(x=d,y=value,colour=model,shape=model),alpha=0.6)+
  facet_grid2(paste('Delay temporal trend:',eta)~paste('Totals temporal trend:',zeta))+
  scale_colour_manual(name=NULL,values = c( "#4363D8", "#E6194B"))+
  scale_shape_manual(name=NULL,values=c(15,3))+
  scale_x_continuous()+
  labs(x='Prediction time difference',y="Prediction interval width",
       title="Simulation experiment",
       subtitle='Prediction interval width for negative case load relationships (delta=-0.2)')+
  theme_minimal()+
  theme(legend.position = "bottom",text = element_text(size=14), axis.title = element_text(size=16), plot.title = element_text(size=20),plot.subtitle=element_text(size=18),
        legend.title = element_text(size=16),strip.text = element_text(size=14))
compare_delta_piw_neg
ggsave(compare_delta_piw_neg, file="Plots/neg_delta_piw.pdf", width=9, height=10)


compare_delta_cov_neg <- ggplot(filter(y_all_summary_average_neg,metric=="Coverage"))+
  geom_line(aes(x=d,y=value,colour=model),alpha=0.6)+
  geom_point(aes(x=d,y=value,colour=model,shape=model),alpha=0.6)+
  facet_grid2(paste('Delay temporal trend:',eta)~paste('Totals temporal trend:',zeta))+
  scale_colour_manual(name=NULL,values = c( "#4363D8", "#E6194B"))+
  scale_shape_manual(name=NULL,values=c(15,3))+
  scale_x_continuous()+
  labs(x='Prediction time difference',y="Coverage",
       title="Simulation experiment",
       subtitle='Coverage for negative case load relationships (delta=-0.2)')+
  theme_minimal()+
  theme(legend.position = "bottom",text = element_text(size=14), axis.title = element_text(size=16), plot.title = element_text(size=20),plot.subtitle=element_text(size=18),
        legend.title = element_text(size=16),strip.text = element_text(size=14))
ggsave(compare_delta_cov_neg, file="Plots/neg_delta_cov.pdf", width=9, height=10)



# separate by eta and zeta 
#positive delta
y_all_summary_metrics_pos<-filter(y_all_summary,delta==(0.2))%>%group_by(model,d,s,eta,zeta)%>%
  summarise(Coverage=mean(y>=`2.5%`&y<=`97.5%`),
            `Prediction interval width`=mean(`97.5%`-`2.5%`),
            `Mean absolute error`=mean(abs(`50%`-y)),
            `Bias`=mean(`50%`-y))%>%
  pivot_longer(cols=c('Coverage','Prediction interval width','Mean absolute error','Bias'),names_to = 'metric')%>%
  mutate(model=factor(model,levels=c('Survivor','Survivor incidence-delay'
  )),
  metric=factor(metric,levels=c('Mean absolute error','Prediction interval width','Coverage','Bias')))

y_all_summary_average_pos<-y_all_summary_metrics_pos%>%group_by(model,d,metric,eta,zeta)%>%
  summarise(value.mean=mean(as.numeric(value)), value.median=median(as.numeric(value)))

y_all_summary_average_pos$value.mean[y_all_summary_average_pos$metric%in%c("Prediction interval width","Mean absolute error","Bias")]<-NA
y_all_summary_average_pos$value.median[y_all_summary_average_pos$metric%in%c("Coverage")]<-NA
y_all_summary_average_pos<-y_all_summary_average_pos%>%group_by(model,d,metric,eta,zeta)%>%mutate(value=sum(value.mean,value.median,na.rm=TRUE))
compare_delta_mae_pos <- ggplot(filter(y_all_summary_average_pos,metric=="Mean absolute error"))+
  geom_line(aes(x=d,y=value,colour=model),alpha=0.6)+
  geom_point(aes(x=d,y=value,colour=model,shape=model),alpha=0.6)+
  facet_grid2(paste('Delay temporal trend:',eta)~paste('Totals temporal trend:',zeta))+
  scale_colour_manual(name=NULL,values = c( "#4363D8", "#E6194B"))+
  scale_shape_manual(name=NULL,values=c(15,3))+
  scale_x_continuous()+
  labs(x='Prediction time difference',y="Mean absolute error",
       title="Simulation experiment",
       subtitle='Mean absolute error for positive case load relationships (delta=0.2)')+
  theme_minimal()+
  theme(legend.position = "bottom",text = element_text(size=14), axis.title = element_text(size=16), plot.title = element_text(size=20),plot.subtitle=element_text(size=18),
        legend.title = element_text(size=16),strip.text = element_text(size=14))
compare_delta_mae_pos
ggsave(compare_delta_mae_pos, file="Plots/pos_delta_mae.pdf", width=9, height=10)

compare_delta_piw_pos <- ggplot(filter(y_all_summary_average_pos,metric=="Prediction interval width"))+
  geom_line(aes(x=d,y=value,colour=model),alpha=0.6)+
  geom_point(aes(x=d,y=value,colour=model,shape=model),alpha=0.6)+
  facet_grid2(paste('Delay temporal trend:',eta)~paste('Totals temporal trend:',zeta))+
  scale_colour_manual(name=NULL,values = c( "#4363D8", "#E6194B"))+
  scale_shape_manual(name=NULL,values=c(15,3))+
  scale_x_continuous()+
  labs(x='Prediction time difference',y="Prediction interval width",
       title="Simulation experiment",
       subtitle='Prediction interval width for positive case load relationships (delta=0.2)')+
  theme_minimal()+
  theme(legend.position = "bottom",text = element_text(size=14), axis.title = element_text(size=16), plot.title = element_text(size=20),plot.subtitle=element_text(size=18),
        legend.title = element_text(size=16),strip.text = element_text(size=14))
compare_delta_piw_pos
ggsave(compare_delta_piw_pos, file="Plots/pos_delta_piw.pdf", width=9, height=10)


compare_delta_cov_pos <- ggplot(filter(y_all_summary_average_pos,metric=="Coverage"))+
  geom_line(aes(x=d,y=value,colour=model),alpha=0.6)+
  geom_point(aes(x=d,y=value,colour=model,shape=model),alpha=0.6)+
  facet_grid2(paste('Delay temporal trend:',eta)~paste('Totals temporal trend:',zeta))+
  scale_colour_manual(name=NULL,values = c( "#4363D8", "#E6194B"))+
  scale_shape_manual(name=NULL,values=c(15,3))+
  scale_x_continuous()+
  labs(x='Prediction time difference',y="Coverage",
       title="Simulation experiment",
       subtitle='Coverage for positive case load relationships (delta=0.2)')+
  theme_minimal()+
  theme(legend.position = "bottom",text = element_text(size=14), axis.title = element_text(size=16), plot.title = element_text(size=20),plot.subtitle=element_text(size=18),
        legend.title = element_text(size=16),strip.text = element_text(size=14))
ggsave(compare_delta_cov_pos, file="Plots/pos_delta_cov.pdf", width=9, height=10)





# Compare metrics for all models:
n_sim<-dim(y_quant_fixed[[1]])[1]/S
y_fixed<-y_quant_fixed%>%unlist()%>%array(dim=c(dim(y_quant_fixed[[1]]),sims))
colnames(y_fixed)<- colnames(y_quant_fixed[[1]])
y_fixed<-apply(y_fixed,2,c)%>%as_tibble()%>%mutate(`2.5%`=as.numeric(`2.5%`),
                                                           `50%`=as.numeric(`50%`),
                                                           `97.5%`=as.numeric(`97.5%`),
                                                           t=as.numeric(t),
                                                           s=as.numeric(s))%>%mutate(model='Survivor incidence-delay (delta fixed)')

# Nowcasts from all models:
y_all<-rbind(y_nolink,y_linear_log,y_fixed)
# Add true y values:
y_all_sim<-right_join(sim_y_long,y_all, by=c('t','s','eta','zeta','delta','sim'))

ggplot()+geom_line(data=filter(y_all_sim,sim==5,delta==c(-0.2)),aes(x=t,y=y))+
  geom_line(data=filter(y_all_sim,sim==5,delta==c(-0.2)),aes(x=t,y=`50%`,colour=model))+
  facet_wrap(zeta~eta)+
  geom_ribbon(data=filter(y_all_sim,sim==5,delta==c(-0.2)),aes(x=t,ymin=`2.5%`,ymax=`97.5%`,fill=model),alpha=0.5)

# Filter for dates that are censored:
y_all_d<-y_all_sim%>%mutate(d=t-N)%>%filter(d>-D)

y_all_summary<-y_all_d%>%mutate(model=factor(model,levels=c('Survivor','Survivor incidence-delay','Survivor incidence-delay (delta fixed)'
)),y=as.numeric(y),
d=as.numeric(d),sim=as.numeric(sim))
y_all_summary_metrics<-y_all_summary%>%group_by(model,d,s,delta)%>%
  summarise(Coverage=mean(y>=`2.5%`&y<=`97.5%`),
            `Prediction interval width`=mean(`97.5%`-`2.5%`),
            `Mean absolute error`=mean(abs(`50%`-y)))%>%
  pivot_longer(cols=c('Coverage','Mean absolute error','Prediction interval width'),names_to = 'metric')%>%
  mutate(model=factor(model,levels=c('Survivor','Survivor incidence-delay','Survivor incidence-delay (delta fixed)'
  )),
  metric=factor(metric,levels=c('Mean absolute error','Prediction interval width','Coverage')))



y_all_summary_metrics_all<-y_all_summary%>%group_by(model,d,s,delta)%>%
  summarise(Coverage=mean(y>=`2.5%`&y<=`97.5%`),
            `Prediction interval width`=mean(`97.5%`-`2.5%`),
            `Mean absolute error`=mean(abs(`50%`-y)),
            #  `Bias`=mean(`50%`-y)
  )%>%
  pivot_longer(cols=c('Coverage','Prediction interval width','Mean absolute error'),names_to = 'metric')%>%
  mutate(model=factor(model,levels=c('Survivor','Survivor incidence-delay','Survivor incidence-delay (delta fixed)'
  )),
  metric=factor(metric,levels=c('Mean absolute error','Prediction interval width','Coverage')))

y_all_summary_average_all<-y_all_summary_metrics_all%>%group_by(model,d,metric,delta)%>%
  summarise(value.mean=mean(as.numeric(value)), value.median=median(as.numeric(value)))

y_all_summary_average_all$value.mean[y_all_summary_average_all$metric%in%c("Prediction interval width","Mean absolute error")]<-NA
y_all_summary_average_all$value.median[y_all_summary_average_all$metric%in%c("Coverage")]<-NA
y_all_summary_average_all<-y_all_summary_average_all%>%group_by(model,d,metric,delta)%>%mutate(value=sum(value.mean,value.median,na.rm=TRUE))
# Plot of prediction performance metrics.

compare_delta_mae <- ggplot(filter(y_all_summary_average_all,metric=="Mean absolute error"))+
  geom_line(aes(x=d,y=value,colour=model),alpha=0.6)+
  geom_point(aes(x=d,y=value,colour=model,shape=model),alpha=0.6)+
  facet_grid2(~paste('Case load effect:',delta),switch=c('y'))+
  scale_colour_manual(name=NULL,values = c( "#4363D8", "#E6194B","#2AA198"))+
  scale_shape_manual(name=NULL,values=c(15,3,2))+
  scale_x_continuous()+
  labs(x=NULL,y="Mean absolute error",
       title="Simulation experiment",
       subtitle='Performance metrics from the prediction experiment')+
  theme_light()+ 
  guides(colour="none", shape='none')+
  theme(legend.position = "bottom",text = element_text(size=14), axis.title = element_text(size=16), plot.title = element_text(size=20),plot.subtitle=element_text(size=18),
        legend.title = element_text(size=16),strip.text = element_text(size=14))

compare_delta_piw <- ggplot(filter(y_all_summary_average_all,metric=="Prediction interval width"))+
  geom_line(aes(x=d,y=value,colour=model),alpha=0.6)+
  geom_point(aes(x=d,y=value,colour=model,shape=model),alpha=0.6)+
  facet_grid2(~paste('Case load effect:',delta))+
  scale_colour_manual(name=NULL,values = c( "#4363D8", "#E6194B","#2AA198"))+
  scale_shape_manual(name=NULL,values=c(15,3,2))+
  scale_x_continuous()+
  labs(x=NULL,y="Prediction interval width",
       title=NULL,
       subtitle=NULL)+
  theme_light()+ 
  guides(colour="none", shape='none')+
  theme(legend.position = "bottom",text = element_text(size=14), axis.title = element_text(size=16), plot.title = element_text(size=20),plot.subtitle=element_text(size=18),
        legend.title = element_text(size=16),strip.text = element_text(size=14))

compare_delta_cov <- ggplot(filter(y_all_summary_average_all,metric=="Coverage"))+
  geom_hline(yintercept=0.95, colour="darkgrey")+
  geom_line(aes(x=d,y=value,colour=model),alpha=0.6)+
  geom_point(aes(x=d,y=value,colour=model,shape=model),alpha=0.6)+
  facet_grid2(~paste('Case load effect:',delta))+
  scale_colour_manual(name=NULL,values = c( "#4363D8", "#E6194B","#2AA198"))+
  scale_shape_manual(name=NULL,values=c(15,3,2))+
  scale_x_continuous()+
  labs(x='Prediction time difference',y="Coverage",
       title=NULL,
       subtitle=NULL)+
  theme_light()+ 
  theme(legend.position = "bottom",text = element_text(size=14), axis.title = element_text(size=16), plot.title = element_text(size=20),plot.subtitle=element_text(size=18),
        legend.title = element_text(size=16),strip.text = element_text(size=14))

compare_delta_sep<-grid.arrange(compare_delta_mae,compare_delta_piw,compare_delta_cov,nrow=3,heights=c(4.1,3.2,4.1))


ggsave(compare_delta_sep, file="Plots/ALL_pred_sep.pdf", width=9, height=11)


