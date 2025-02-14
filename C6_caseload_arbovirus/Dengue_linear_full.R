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
# Convert the html table eleument into a data frame
my.table = html_table(my.table, fill = TRUE)

my.table<-my.table[-1,]
#setwd("~/OneDrive/PhD/Dengue 2023")
#setwd("~/OneDrive - University of Glasgow/PhD/Dengue 2023")
#setwd("~/Documents/Dengue 2023")

states<-read_state(year=2013)


data_set <- read_csv("data set.csv")
data_set<-data_set[,-c(22:48)]


data<-data_set%>%mutate(start_of_symptoms=as.Date((DT_SIN_PRI),format="%Y-%m-%d"),
                        notification_date=as.Date((DT_NOTIFIC),format="%Y-%m-%d"))%>%
  mutate(delay=difftime(notification_date,start_of_symptoms, unit=c("days")))%>%
 # select(start_of_symptoms,notification_date,delay,CLASSI_FIN,SG_UF_NOT)%>%
  drop_na(start_of_symptoms)%>%filter(start_of_symptoms>"2012-12-31",start_of_symptoms<"2020-12-27")
data<-data[,which(colnames(data)%in%c("start_of_symptoms","notification_date","delay","CLASSI_FIN","SG_UF_NOT"))]
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
states_short$name_state<-as.character(states_short$name_state)
# create guide between different federative unit labels
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

Data_delay_sum<-as.tibble(Data_wide)
totals<-Data_delay_sum%>%mutate(total=as.numeric(unlist((apply(Data_wide[,3:34],1, sum, na.rm=TRUE)))))%>%mutate(s=paste(federal_guide$abbrev_state[s]))
Data_delay_sum$prop_0<-as.numeric(unlist((Data_wide[,3])/apply(Data_wide[,3:34],1, sum, na.rm=TRUE)))
Data_delay_sum$prop_1<-as.numeric(unlist((Data_wide[,3]+Data_wide[,4])/apply(Data_wide[,3:34],1, sum, na.rm=TRUE)))
Data_delay_sum$prop_2<-as.numeric(unlist((Data_wide[,3]+Data_wide[,4]+Data_wide[,5])/apply(Data_wide[,3:34],1, sum, na.rm=TRUE)))
Data_delay_sum$prop_3<-as.numeric(unlist((Data_wide[,3]+Data_wide[,4]+Data_wide[,5]+Data_wide[,6])/apply(Data_wide[,3:34],1, sum, na.rm=TRUE)))
Data_delay_sum$prop_4<-as.numeric(unlist((Data_wide[,3]+Data_wide[,4]+Data_wide[,5]+Data_wide[,6]+Data_wide[,7])/apply(Data_wide[,3:34],1, sum, na.rm=TRUE)))

# Data_delay_sum$prop_0<-as.numeric(unlist((Data_wide[,3])/apply(Data_wide[,3:34],1, sum, na.rm=TRUE)))
# Data_delay_sum$prop_1<-as.numeric(unlist((Data_wide[,4])/(apply(Data_wide[,3:34],1, sum, na.rm=TRUE)-Data_wide[,3])))
# Data_delay_sum$prop_2<-as.numeric(unlist((Data_wide[,5])/(apply(Data_wide[,3:34],1, sum, na.rm=TRUE)-apply(Data_wide[,3:4],1, sum, na.rm=TRUE))))
# Data_delay_sum$prop_3<-as.numeric(unlist((Data_wide[,6])/(apply(Data_wide[,3:34],1, sum, na.rm=TRUE)-apply(Data_wide[,3:5],1, sum, na.rm=TRUE))))
# Data_delay_sum$prop_4<-as.numeric(unlist((Data_wide[,7])/(apply(Data_wide[,3:34],1, sum, na.rm=TRUE)-apply(Data_wide[,3:6],1, sum, na.rm=TRUE))))

Data_delay_long<-as.tibble(Data_delay_sum%>%pivot_longer(!c(t,s), names_to = "delay", values_to = "count")%>%mutate(s=paste(federal_guide$abbrev_state[s])))
population<-my.table[,c(1,2,6)]
colnames(population)<-c("s_full","s","population")
population[,3]<-as.numeric(str_replace_all(as.character(unlist(my.table[,6])),",",""))
population_order<-full_join(federal_guide,population%>%mutate(abbrev_state=s), by="abbrev_state")

Data_delay_pop<-full_join(Data_delay_long,population, by="s")
Data_delay_prop<-Data_delay_pop%>%filter(delay%in%c("prop_0","prop_1","prop_2","prop_3","prop_4"))
Data_delay_prop<-full_join(Data_delay_prop,totals[,c(1,2,dim(totals)[2])],by = join_by(t, s))
Data_delay_prop<-Data_delay_prop%>%mutate(total.per.capita=total/population)
# plot cumulative proportions against totals
ggplot(Data_delay_prop,aes(x=total,y=probit(count), colour=delay))+
  geom_point(alpha=0.1)+facet_wrap(~s,scales="free")+geom_smooth(alpha=0.1)+scale_x_log10()+
  labs(x="Total counts",y="Cumulative proportions")
ggplot(Data_delay_prop,aes(x=total,y=probit(count), colour=s))+
  geom_point(alpha=0.1)+facet_wrap(~delay,scales="free")+geom_smooth(alpha=0.1)+scale_x_log10()+
  labs(x="Total counts",y="Cumulative proportions")


# plot cumulative proportions against totals per capita
ggplot(Data_delay_prop,aes(x=total.per.capita,y=count, colour=delay))+
  geom_point()+facet_wrap(~s,scales="free")+geom_smooth()+
  labs(x="Total counts per capita",y="Cumulative proportions")

test_model <- gam(count~population+s(total)+s,data=filter(Data_delay_prop,delay=="prop_1"))
summary(test_model)
plot(test_model)

test_model <- gam(count~population+s(total)+s,data=filter(Data_delay_prop,delay=="prop_3",t>200))
summary(test_model)
plot(test_model)
# Matrix for partial counts [REGION,TIME,DELAY]
Dengue_z<-as.numeric(Data_wide[,-(1:2)])%>%array(dim=c(N,S,D_N))%>%aperm(c(2,1,3))
Dengue_z[is.na(Dengue_z)]<-0
# Matrix for total counts [REGION,TIME]
Dengue_y<-apply(Dengue_z,c(1,2),sum)

y_long<-melt(Dengue_y, value.name = 'y', varnames = c("s","t"))

# Define beta binomial distribution.
rbetabin=nimbleFunction(run=function(n=integer(0),mu=double(0),phi=double(0),size=double(0)){
  pi=rbeta(n,mu*phi,(1-mu)*phi)
  returnType(double(1))
  return(rbinom(n,size,pi))
})




Dengue_surv_full_region <- function(seed, Dengue_y, Dengue_z, population_order, Nowcast_list, S, N, D, D_max, n_knots, n_iter, n_burn, n_thin, n_chains){
  # Load libraries.
  library(nimble)
  library(tidyverse)
  library(mgcv)
  
  # Define the Beta-Binomial as a distribution for NIMBLE.
  dbetabin=nimbleFunction(run=function(x=double(0),mu=double(0),phi=double(0),size=double(0),log=integer(0)){
    phi <- min(phi,1e+04) # Hard upper limit on phi for computational stability.
    returnType(double(0))
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
  
  # Register the Beta-Binomial as a distribution.
  registerDistributions(list(dbetabin=list(
    BUGSdist='dbetabin(mu,phi,size)',discrete=TRUE)))
  
  # NIMBLE model code. 
  Dengue_code <- nimbleCode({
      for(t in 1:N){
        # Negative Binomial Model for total Dengue cases
        y[t] ~ dnegbin(prob=theta/(theta+lambda[t]),size=theta)  
        
        # Model for delayed Dengue counts
        z[t,1] ~ dbetabin(nu[t,1],phi[1],y[t])
      }
      for(d in 2:D){
        for(t in 1:N){
          z[t,d] ~ dbetabin(nu[t,d],phi[d],y[t]-sum(z[t,1:(d-1)]))
        }
      }
      
      for(t in 1:N){
        # Linear predictor (Beta-Binomial means).
        # index_delta[t]<-delta_converter(x=(log((y[t]+1)/pop)),grid =delta_grid[1:length_delta], grid_length = length_delta)[1]
        for(d in 1:D){
          # Expected cumulative proportions.
          probit(p[t,d]) <- psi[d] + eta[t] + delta*log((y[t]+1)/pop)
        }
        # Relative proportions (Beta-Binomial means).
        nu[t,1] <- p[t,1]
        for(d in 2:D){
          nu[t,d] <- (p[t,d]-p[t,d-1])/(1-p[t,d-1])
        }
        # Mean for total Dengue cases 
        log(lambda[t]) <-  log(pop) + int_zeta + zeta[t] + xi[weeks[t]] 
        
      }
      
      
      # logsplines of time in cumulative proportion reported.
      Omega_eta[1:K_t,1:K_t] <- S_t[1:K_t,1:K_t]/sigma_eta^2
      kappa_eta[1:K_t] ~ dmnorm(zeros[1:K_t],Omega_eta[1:K_t,1:K_t])
      eta[1:N] <- X_t[1:N,1:K_t]%*%kappa_eta[1:K_t]
      
      # Seasonal logspline 
      Omega_xi[1:K_w,1:K_w]<-S_w[1:K_w,1:K_w]/sigma_xi^2
      kappa_xi[1:K_w]~dmnorm(zeros[1:K_w],Omega_xi[1:K_w,1:K_w])
      xi[1:52]<-X_w[1:52,1:K_w]%*%kappa_xi[1:K_w]
      
      # logsplines of time for Dengue cases.
      Omega_zeta[1:K_t,1:K_t] <- S_t[1:K_t,1:K_t]/sigma_zeta^2
      kappa_zeta[1:K_t] ~ dmnorm(zeros[1:K_t],Omega_zeta[1:K_t,1:K_t])
      zeta[1:N] <- X_t[1:N,1:K_t]%*%kappa_zeta[1:K_t]
      
      ##for(d in 1:D){
      #  Omega_delta[1:K_c,1:K_c] <- S_c[1:K_c,1:K_c]/sigma_delta[1]^2 + S_c[1:K_c,(K_c+1):(2*K_c)]/sigma_delta[2]^2 #non-linear plus linear spline 
      # kappa_delta[1:K_c] ~ dmnorm(zeros[1:K_c],Omega_delta[1:K_c,1:K_c])
      # delta[1:length_delta] <- X_c[1:length_delta,1:K_c]%*%kappa_delta[1:K_c]
      
      
      #for(k in 1:2){
      # sigma_delta[k] ~ T(dnorm(0,1),0,) 
      ## }
      # }
      # Smoothing parameter priors.
      sigma_zeta ~  T(dnorm(0,sd=10),0,)
      sigma_eta ~  T(dnorm(0,sd=10),0,)
      sigma_xi ~ T(dnorm(0,sd=10),0,) 
      
      # Cumulative delay curves.
      psi[1] ~ dnorm(0,sd=5)
      for(d in 2:D){
        psi[d] ~ T(dnorm(psi[d-1],sd=5),psi[d-1],)
      }
      
      # Beta-Binomial dispersion parameters.
      for(d in 1:D){
        phi[d] ~ dgamma(2,0.02) 
      }
      
      # Intercepts
      int_zeta  ~ dnorm(-10,sd=10)
      theta ~ dgamma(2,0.02) # Negative-Binomial dispersion parameters.
      
      # Linear link between case load and delay. 
      delta ~ dnorm(0,5)
      
    
    
    
  })
  
  # Data for current nowcast date 
 
  # Censor the totals.
  week_index<-rep(1:52,9) #first date is 01-01-2013
  # Set up logspline bases.
  blank_data=tibble(y=rnorm(N,0,1),t=1:N,w=week_index[1:N]) 
  blank_jagam=jagam(y~s(t,k=floor(N/20), bs='cs')+
                      s(w,bs='cc',k=n_knots[2]),data=blank_data,file='blank.jags',
                    knots=list(t=seq(1,N,length=floor(N/20)),w=seq(0,52,length=n_knots[2]))) # #
  
  
  
  
  # Constants
  Dengue_constants <- list(N=N,
                           D=D,
                           weeks=week_index[1:N],
                           pop=as.numeric(population_order$population)[seed],
                           K_w=dim(blank_jagam$jags.data$S2)[1],
                           S_w=blank_jagam$jags.data$S2,
                           K_t=dim(blank_jagam$jags.data$S1)[1],
                           S_t=blank_jagam$jags.data$S1)
  
  
  Dengue_constants$X_t <- blank_jagam$jags.data$X[,2:(Dengue_constants$K_t+1)]
  Dengue_constants$X_w <- blank_jagam$jags.data$X[,(Dengue_constants$K_t+2):(Dengue_constants$K_t+Dengue_constants$K_w+1)]
  Dengue_constants$zeros <- rep(0,max(Dengue_constants$K_t,Dengue_constants$K_w)) 
  
  # Data 
  Dengue_data <- list(y=Dengue_y[seed,1:N], z=Dengue_z[seed,1:N,1:D]) 
  
  # Generate random initial values.
  Dengue_inits <- list(
    kappa_zeta=rnorm((Dengue_constants$K_t),0,0.1),
    kappa_eta=rnorm(S*(Dengue_constants$K_t),0,0.1),
    kappa_xi=rnorm(S*(Dengue_constants$K_w),0,0.1),
    # kappa_delta=rnorm((Dengue_constants$K_c),0,0.1),
    delta=rnorm(1,0,0.1),#matrix(rnorm(D*S,0,1),nrow=D),
    #   kappa_delta=matrix(rnorm(S*(Dengue_constants$K_c),0,0.1),ncol=S),
    phi=abs(rnorm(D,30,10)),
    psi=sort(rnorm(D,0,0.1)),
    sigma_eta=runif(1,0,1),
    sigma_zeta=runif(1,0,1),
    sigma_xi=runif(1,0,1),
    #    sigma_delta=runif(2,0,1),#matrix(runif(2*S,0,1),nrow=2),
    theta=abs(rnorm(1,50,10)),
    int_zeta=rnorm(1,-10,1))

  
  # Build the model.
  Dengue_model <- nimbleModel(Dengue_code,Dengue_constants,Dengue_data,Dengue_inits)
  # Compile the model.
  Dengue_compiled_model <- compileNimble(Dengue_model)
  
  # Set up the MCMC.
  Dengue_mcmc_config <- configureMCMC(Dengue_model,monitors=c("lambda","int_zeta","xi",
                                                              "psi","theta","phi", 
                                                              "sigma_eta","sigma_zeta",
                                                              "eta","zeta","delta",
                                                              "kappa_zeta","kappa_eta"),
                                      useConjugacy = FALSE)
  
  Dengue_mcmc_config$removeSamplers('kappa_zeta','sigma_zeta','int_zeta',
                                    'kappa_eta','sigma_eta',"psi",
                                    'delta',
                                    'kappa_xi','sigma_xi')

    Dengue_mcmc_config$addSampler(target=c(paste('kappa_zeta[1:',(Dengue_constants$K_t),']', sep=''),paste('sigma_zeta',sep=''),
                                           paste('kappa_xi[1:',(Dengue_constants$K_w),']', sep=''),paste('sigma_xi',sep=''),
                                           paste('int_zeta',sep='')),type='AF_slice') 
    Dengue_mcmc_config$addSampler(target=c(paste('delta', sep=''),
                                           paste('kappa_eta[1:',(Dengue_constants$K_t),']', sep=''),paste('sigma_eta',sep=''),
                                           paste('psi[1:',D,']', sep='')),type='AF_slice')

  
  Dengue_mcmc<- buildMCMC(Dengue_mcmc_config)
  # Compile the MCMC.
  Dengue_compiled_mcmc <- compileNimble(Dengue_mcmc,project=Dengue_model)
  
  Dengue_output_cluster<- runMCMC(Dengue_compiled_mcmc,niter=n_iter,nburnin=n_burn,thin=n_thin,inits=Dengue_inits,nchains=n_chains,samplesAsCodaMCMC = TRUE)
  
  return(Dengue_output_cluster)
}

D_all=32
D_max=30
n_iter=10000
n_burn=5000
n_thin=5
n_chains=2
n_knots=c(18,7,10) 
D=4
N_max<-ncol(Dengue_y)
N<-375
y_long<-filter(y_long,t<=N)
# lower_delta= min(log((1/population_order$population)))
# upper_delta=  max(log((Dengue_y+1)/population_order$population))+4
# by_delta=0.1
# length_delta=1000
# delta_grid<-seq(from=lower_delta, to=upper_delta, length=length_delta)
# grid_length<-length(delta_grid)
#Nowcast_list=floor(seq(from=2*52, to=375, length=4))

population_order$population_scale<-scale(population_order$population)
# Number of available computer cores 
n_cores<-min(detectCores(),S,5)
# Make Cluster for MCMC. 
this_cluster_full<- makeCluster(n_cores)



# Note: takes approximately 
time_full<-system.time({
  Dengue_surv_full <- parLapply(cl = this_cluster_full, X = 1:S, 
                                fun = Dengue_surv_full_region,
                                Dengue_y=as.array(Dengue_y)[1:S,1:N],
                                Dengue_z=as.array(Dengue_z)[1:S,1:N,],
                                population_order=population_order,
                                N=N,#
                                S=S,
                                D=D,
                                D_max=D_max,
                                n_iter=n_iter,
                                n_burn=n_burn,
                                n_thin=n_thin,
                                n_chains=n_chains,
                                n_knots=n_knots)
})


stopCluster(this_cluster_full)

#SAVE Samples
save(Dengue_surv_full,time_full,file='Plots/Dengue_full_linear.RData')
#


# Check PSRF of model parameters.
Dengue_parameter_group=c("int_zeta","xi","lambda",
                         "psi",#"delta",
                         "sigma_eta","sigma_zeta", "zeta",
                         "kappa_zeta","kappa_eta")
Dengue_psrf<-list()
Dengue_mean_psrf<-Dengue_max_psrf<-Dengue_leq_105 <-Dengue_leq_110<-array(NA,dim=c(length(Dengue_parameter_group),length(Nowcast_list)))
for(k in 1:length(Dengue_parameter_group)){
  Dengue_psrf[[k]]<-list()  
  for(j in 1:S){
    Dengue_parameter_names <- colnames((as.matrix(Dengue_surv_full[[j]])))
    Dengue_index <- which(startsWith(Dengue_parameter_names,Dengue_parameter_group[k]))
    Dengue_psrf[[k]][[j]]<- sapply(Dengue_index,function(v)gelman.diag(Dengue_surv_full[[j]][,v],autoburnin = FALSE,transform = TRUE,multivariate = FALSE)$psrf[,1])
    Dengue_mean_psrf[k,j] <- mean(Dengue_psrf[[k]][[j]])
    Dengue_max_psrf[k,j] <- max(Dengue_psrf[[k]][[j]])
    Dengue_leq_105[k,j] <- mean(Dengue_psrf[[k]][[j]]<=1.05)
    Dengue_leq_110[k,j] <- mean(Dengue_psrf[[k]][[j]][j]<=1.1)
  }
}
Dengue_MAX_psrf<-cbind(Dengue_parameter_group,apply(Dengue_max_psrf,1,max))
Dengue_MEAN_psrf<-cbind(Dengue_parameter_group,apply(Dengue_mean_psrf,1,mean))



plot(Dengue_surv_full[[1]][,"delta"],main="delta")
plot(Dengue_surv_full[[2]][,"delta"],main="delta")

#### Plots #####
library(geobr)
states<-read_state(year=2013)
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
delta_slope<-zeta_spline<-eta_spline<-psi_d<-list()
for(j in 1:S){
  delta_slope[[j]]<-select(Dengue_combined_Survivor_linear_x[[j]],starts_with('delta'))%>%as.matrix()%>%array(dim=c(n_sim))%>%
    quantile(c(0.025,0.5,0.975))%>%melt(varnames=c('quantile'),value.name='y')%>%mutate(s=j,quantile=c(0.025,0.5,0.975))%>%
    spread(quantile,y)

  zeta_spline[[j]]<-select(Dengue_combined_Survivor_linear_x[[j]],starts_with('zeta['))%>%as.matrix()%>%array(dim=c(n_sim,N))%>%
    apply(c(2),quantile,c(0.025,0.5,0.975))%>%melt(varnames=c('quantile','t'),value.name='y')%>%mutate(s=j)%>%
    spread(quantile,y)
  
  eta_spline[[j]]<-select(Dengue_combined_Survivor_linear_x[[j]],starts_with('eta'))%>%as.matrix()%>%array(dim=c(n_sim,N))%>%
    apply(c(2),quantile,c(0.025,0.5,0.975))%>%melt(varnames=c('quantile','t'),value.name='y')%>%mutate(s=j)%>%
    spread(quantile,y)
  
  psi_d[[j]]<-select(Dengue_combined_Survivor_linear_x[[j]],starts_with('psi'))%>%as.matrix()%>%array(dim=c(n_sim,D))%>%
    apply(c(2),quantile,c(0.025,0.5,0.975))%>%melt(varnames=c('quantile','d'),value.name='y')%>%mutate(s=j)%>%
    spread(quantile,y)
  
  
  
}

# delta matrix:
delta_matrix<-matrix(unlist(delta_slope), byrow = TRUE, ncol = 4)%>%as_tibble()
colnames(delta_matrix)<-c('s','2.5%','50%','97.5%')
delta_matrix$s_abb<-factor(federal_guide$abbrev_state[delta_matrix$s],levels = sort(federal_guide$abbrev_state[delta_matrix$s]))
delta_matrix<-delta_matrix[order(delta_matrix$s_abb),]
delta_matrix$s_full<-factor(federal_guide$name_state[delta_matrix$s], levels = (federal_guide$name_state[delta_matrix$s]))
# eta matrix:
eta_matrix<-as_tibble(eta_spline[[1]],ncol = 5)
for(s in 2:S){
  eta_matrix<-rbind(eta_matrix, eta_spline[[s]])
  
}
colnames(eta_matrix)<-c('t','s','2.5%','50%','97.5%')

# zeta matrix:
zeta_matrix<-as_tibble(zeta_spline[[1]],ncol = 5)
for(s in 2:S){
  zeta_matrix<-rbind(zeta_matrix, zeta_spline[[s]])
  
}
colnames(zeta_matrix)<-c('t','s','2.5%','50%','97.5%')

# 

delta_comp_plot<-ggplot()+geom_hline(yintercept = 0,colour='darkgrey',linetype='dashed')+
  geom_errorbar(data=delta_matrix,
                aes(x=s_abb,  ymin=`2.5%`,ymax=`97.5%`,colour=s_full))+
  geom_point(data=delta_matrix,aes(x=s_abb, y=`50%`,colour=s_full),shape=3,size=5)+
  scale_colour_viridis_d(name='Region',begin=0.1,end=0.9)+
  scale_fill_viridis_d(name='Region',begin=0.1,end=0.9)+
  labs(x="Brazilian federative unit (s)",y=expression(delta['s']),
       title='Arbovirus cases',
       subtitle='Linear effect of case load on reporting delay', legend=NULL)+
  theme_minimal()+
  theme(text = element_text(size=12), plot.title = element_text(size=16),plot.subtitle=element_text(size=14),
        legend.title = element_text(size=14),strip.text = element_text(size=14), 
        legend.position = "bottom",legend.box="vertical")+guides(colour=guide_legend(nrow=6,byrow=TRUE))
delta_comp_plot
ggsave(delta_comp_plot,file='Plots/full_delta_dengue.pdf',width=9,height=6)

geom_guide<-tibble(s_abb=states$abbrev_state,s_full=states$name_state,geom=states$geom)
delta_matrix$s_abb<-as.character(delta_matrix$s_abb)
delta_matrix_geom<-full_join(delta_matrix,geom_guide, by='s_abb')

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
 # scale_fill_fermenter(palette = "Spectral",limit=limit,name=expression(delta['s']),direction = 1)+ 
  scale_fill_gradient2(low = red,
                       mid = yellow,
                       high = green,name=expression(delta['s']))+  labs(title="Arbovirus cases",
                                              subtitle = "Linear case load effect of reporting delay by federative unit", size=8) +
  theme_minimal()
delta_link_map
ggsave(delta_link_map,file='Plots/delta_link_map_dengue_cb.pdf',width=9,height=7,dpi = 600)





# 
my.table<-my.table[,c(1,2,3,4,5,7,8,10,11)]
my.table[,5]<-as.numeric(str_replace_all(as.character(unlist(my.table[,5])),",",""))
my.table[,6]<-as.numeric(str_replace_all(as.character(unlist(my.table[,6])),",",""))
my.table[,7]<-as.numeric(unlist(my.table[,7]))
my.table[,8]<-as.numeric(str_replace_all(as.character(unlist(my.table[,8])),",",""))
my.table[,9]<-as.numeric(unlist(my.table[,9]))
colnames(my.table)<-c('name.full','s_abb','capital','largest.city','area.km2','population','density.perkm2','GDP','HDI')
#load("~/media/alba/Disk 1/OneDrive/PhD/comp 2023 paper code/FInal_code/OMEGA CODE/SARI totals link/health_expenditure.RData")

delta_explore<-full_join(delta_matrix,my.table,by='s_abb')
delta_explore<-full_join(delta_explore,health_expenditure, by='name.full')
delta_explore<-full_join(delta_explore,geom_guide, by=c('s_abb','s_full'))

delta_explore<-delta_explore%>%mutate(mean.cases=apply(Dengue_y,1,mean))

delta_cases<-ggplot()+geom_hline(yintercept = 0,colour='darkgrey',linetype='dashed')+
  geom_errorbar(data=delta_explore,
                aes(x=(mean.cases),  ymin=`2.5%`,ymax=`97.5%`,colour=s_abb))+
  geom_point(data=delta_explore,aes(x=(mean.cases), y=`50%`,colour=s_abb),shape=3,size=4)+
  scale_colour_viridis_d(name='Region',begin=0.1,end=0.9)+
  scale_fill_viridis_d(name='Region',begin=0.1,end=0.9)+
  scale_x_log10()+
  labs(x="Mean total cases",y=expression(delta['s']),title='Linear effect of case load on delay', legend=NULL)+
  theme(text = element_text(size=12), plot.title = element_text(size=16),plot.subtitle=element_text(size=14),
        legend.title = element_text(size=14),strip.text = element_text(size=14))
delta_cases
ggsave(delta_cases,file='Plots/delta_cases.pdf',width=9,height=3)

delta_pop<-ggplot()+geom_hline(yintercept = 0,colour='darkgrey',linetype='dashed')+
  geom_errorbar(data=delta_explore,
                aes(x=population,  ymin=`2.5%`,ymax=`97.5%`,colour=s_abb))+
  geom_point(data=delta_explore,aes(x=population, y=`50%`,colour=s_abb),shape=3,size=4)+
  scale_colour_viridis_d(name='Region',begin=0.1,end=0.9)+
  scale_fill_viridis_d(name='Region',begin=0.1,end=0.9)+
  labs(x="Population (log)",y=expression(delta['s']),title='Linear effect of case load on delay', legend=NULL)+
  theme(text = element_text(size=12), plot.title = element_text(size=16),plot.subtitle=element_text(size=14),
        legend.title = element_text(size=14),strip.text = element_text(size=14))+scale_x_log10()
delta_pop
ggsave(delta_pop,file='Plots/delta_pop.pdf',width=9,height=3)

delta_area<-ggplot()+geom_hline(yintercept = 0,colour='darkgrey',linetype='dashed')+
  geom_errorbar(data=delta_explore,
                aes(x=area.km2,  ymin=`2.5%`,ymax=`97.5%`,colour=s_abb))+
  geom_point(data=delta_explore,aes(x=area.km2, y=`50%`,colour=s_abb),shape=3,size=4)+
  scale_colour_viridis_d(name='Region',begin=0.1,end=0.9)+
  scale_fill_viridis_d(name='Region',begin=0.1,end=0.9)+
  labs(x="Area (kM squared log scale)",y=expression(delta['s']),title='Linear effect of case load on delay', legend=NULL)+
  theme(text = element_text(size=12), plot.title = element_text(size=16),plot.subtitle=element_text(size=14),
        legend.title = element_text(size=14),strip.text = element_text(size=14))+scale_x_log10()
delta_area
ggsave(delta_area,file='Plots/delta_area.pdf',width=9,height=3)

delta_density<-ggplot()+geom_hline(yintercept = 0,colour='darkgrey',linetype='dashed')+
  geom_errorbar(data=delta_explore,
                aes(x=density.perkm2,  ymin=`2.5%`,ymax=`97.5%`,colour=s_abb))+
  geom_point(data=delta_explore,aes(x=density.perkm2, y=`50%`,colour=s_abb),shape=3,size=4)+
  scale_colour_viridis_d(name='Region',begin=0.1,end=0.9)+
  scale_fill_viridis_d(name='Region',begin=0.1,end=0.9)+
  labs(x="Density (per Km squared log scale)",y=expression(delta['s']),title='Linear effect of case load on delay', legend=NULL)+
  theme(text = element_text(size=12), plot.title = element_text(size=16),plot.subtitle=element_text(size=14),
        legend.title = element_text(size=14),strip.text = element_text(size=14))+scale_x_log10()
delta_density
ggsave(delta_density,file='Plots/delta_density.pdf',width=9,height=3)

delta_GDP<-ggplot()+geom_hline(yintercept = 0,colour='darkgrey',linetype='dashed')+
  geom_errorbar(data=delta_explore,
                aes(x=GDP,  ymin=`2.5%`,ymax=`97.5%`,colour=s_abb))+
  geom_point(data=delta_explore,aes(x=GDP, y=`50%`,colour=s_abb),shape=3,size=4)+
  scale_colour_viridis_d(name='Region',begin=0.1,end=0.9)+
  scale_fill_viridis_d(name='Region',begin=0.1,end=0.9)+
  labs(x="GDP (log)",y=expression(delta['s']),title='Linear effect of case load on delay', legend=NULL)+
  theme(text = element_text(size=12), plot.title = element_text(size=16),plot.subtitle=element_text(size=14),
        legend.title = element_text(size=14),strip.text = element_text(size=14))+scale_x_log10()
ggsave(delta_GDP,file='Plots/delta_GDP.pdf',width=9,height=3)

delta_HDI<-ggplot()+geom_hline(yintercept = 0,colour='darkgrey',linetype='dashed')+
  geom_errorbar(data=delta_explore,
                aes(x=HDI,  ymin=`2.5%`,ymax=`97.5%`,colour=s_abb))+
  geom_point(data=delta_explore,aes(x=HDI, y=`50%`,colour=s_abb),shape=3,size=4)+
  scale_colour_viridis_d(name='Region',begin=0.1,end=0.9)+
  scale_fill_viridis_d(name='Region',begin=0.1,end=0.9)+
  labs(x="HDI",y=expression(delta['s']),title='Linear effect of case load on delay', legend=NULL)+
  theme(text = element_text(size=12), plot.title = element_text(size=16),plot.subtitle=element_text(size=14),
        legend.title = element_text(size=14),strip.text = element_text(size=14))
ggsave(delta_HDI,file='Plots/delta_GDP.pdf',width=9,height=3)


delta_expend<-ggplot()+geom_hline(yintercept = 0,colour='darkgrey',linetype='dashed')+
  geom_errorbar(data=delta_explore,
                aes(x=expenditure,  ymin=`2.5%`,ymax=`97.5%`,colour=s_abb))+
  geom_point(data=delta_explore,aes(x=expenditure, y=`50%`,colour=s_abb),shape=3,size=4)+
  scale_colour_viridis_d(name='Region',begin=0.1,end=0.9)+
  scale_fill_viridis_d(name='Region',begin=0.1,end=0.9)+
  labs(x="Goverment health expenditure (billion Brazillian reals log scale)",y=expression(delta['s']),title='Linear effect of case load on delay', legend=NULL)+
  theme(text = element_text(size=12), plot.title = element_text(size=16),plot.subtitle=element_text(size=14),
        legend.title = element_text(size=14),strip.text = element_text(size=14))#+scale_x_log10()
ggsave(delta_expend,file='Plots/delta_expend.pdf',width=9,height=3)


delta_explore_all<-delta_explore[,c(3,5,7,12,13,15,17)]%>%pivot_longer(cols = 4:7)

ALL_dengue<-ggplot()+
  geom_smooth(data=delta_explore_all,aes(x=value, y=`50%`), method='lm', alpha=0.2, colour="grey")+
  geom_hline(yintercept = 0,colour='darkgrey',linetype='dashed')+
  facet_wrap(~name, ncol=1, scales = "free_x",  strip.position="top",
             labeller = as_labeller(c(density.perkm2='Population density (per Km squared on log scale):',
                                                                        expenditure='Goverment health expenditure (billion Brazillian reals on log scale):',
                                                                        GDP='Gross domestic product (on log scale):',
                                                                        mean.cases='Mean arbovirus cases over all years (on log scale):')))+
  geom_label(data=delta_explore_all,aes(x=value, y=`50%`,label=s_abb,colour=name.full),alpha=0.5)+
  scale_colour_viridis_d(name='Region',begin=0.1,end=0.9)+
  scale_fill_viridis_d(name='Region',begin=0.1,end=0.9)+
  labs(x=NULL,y=expression(delta['s']),title='Arbovirus cases',subtitle='Linear effect of case load on cumulative proportions reported', legend=NULL)+
  scale_x_log10()+theme_minimal()+ 
  theme(text = element_text(size=12), plot.title = element_text(size=16),plot.subtitle=element_text(size=14),
        legend.title = element_text(size=12),strip.text = element_text(size=12), legend.position = "right")+guides(colour=guide_legend(ncol=1,byrow=TRUE))
ggsave(ALL_dengue,file='Plots/delta_dengue_all.pdf',width=9,height=11)



mean_cases_dengue<-ggplot()+geom_hline(yintercept = 0,colour='darkgrey',linetype='dashed')+
  geom_label(data=delta_explore,aes(x=mean.cases, y=`50%`,label=s_abb,colour=name.full))+
  scale_colour_viridis_d(name='Region',begin=0.1,end=0.9)+
  scale_fill_viridis_d(name='Region',begin=0.1,end=0.9)+
  labs(x="Mean number of arbovirus cases over all years (log scale)",y=expression(delta['s']),title='Linear effect of case load on delay', legend=NULL)+
  theme(text = element_text(size=12), plot.title = element_text(size=16),plot.subtitle=element_text(size=14),
        legend.title = element_text(size=14),strip.text = element_text(size=14))+scale_x_log10()+theme_minimal()+ 
  theme(text = element_text(size=12), plot.title = element_text(size=16),plot.subtitle=element_text(size=14),
        legend.title = element_text(size=12),strip.text = element_text(size=14), legend.position = "bottom")

ggplot()+geom_hline(yintercept = 0,colour='darkgrey',linetype='dashed')+
  geom_label(data=delta_explore,aes(x=expenditure/mean.cases, y=`50%`,label=s_abb,colour=name.full))+
  scale_colour_viridis_d(name='Region',begin=0.1,end=0.9)+
  scale_fill_viridis_d(name='Region',begin=0.1,end=0.9)+
  labs(x="Goverment health expenditure/mean cases (billion Brazillian reals log scale)",y=expression(delta['s']),title='Linear effect of case load on delay', legend=NULL)+
  theme(text = element_text(size=12), plot.title = element_text(size=16),plot.subtitle=element_text(size=14),
        legend.title = element_text(size=14),strip.text = element_text(size=14))+scale_x_log10()+theme_minimal()+ 
  theme(text = element_text(size=12), plot.title = element_text(size=16),plot.subtitle=element_text(size=14),
        legend.title = element_text(size=12),strip.text = element_text(size=14), legend.position = "bottom")

ggplot()+
  geom_label(data=delta_explore,aes(y=expenditure, x=mean.cases,label=s_abb,colour=name.full))+
  scale_colour_viridis_d(name='Region',begin=0.1,end=0.9)+
  scale_fill_viridis_d(name='Region',begin=0.1,end=0.9)+
  labs(y="Goverment health expenditure (billion Brazillian reals log scale)",x="Mean number of arbovirus cases over all years",
       title='Linear effect of case load on delay', legend=NULL)+
  theme(text = element_text(size=12), plot.title = element_text(size=16),plot.subtitle=element_text(size=14),
        legend.title = element_text(size=14),strip.text = element_text(size=14))+scale_x_log10()+scale_y_log10()+theme_minimal()+ 
  theme(text = element_text(size=12), plot.title = element_text(size=16),plot.subtitle=element_text(size=14),
        legend.title = element_text(size=12),strip.text = element_text(size=14), legend.position = "bottom")


ggplot() +
  geom_sf(data=states,colour="grey20", size=.05) +
  geom_sf(data=delta_explore,aes(geometry=geom, fill=s_abb, colour=s_abb),colour=NA, size=.15) +
  # scale_fill_fermenter(palette = "Spectral",limit=limit,name=expression(delta['s']),direction = 1)+ 
#  scale_fill_gradient2(low = lightRed,
      #                 mid = lightBlue,
           #            high = darkGreen,name=expression(delta['s']))+  
  labs(title="Arbovirus cases",
                                                                            subtitle = "Linear case load effect of reporting delay by federative unit", size=8) +
  theme_minimal()




eta_comp_plot<-ggplot()+
  geom_ribbon(data=eta_matrix, 
              aes(x=t,  ymin=`2.5%`,ymax=`97.5%`,fill=as.factor(federal_guide$abbrev_state[s])),alpha=0.2)+
  geom_line(data=eta_matrix,aes(x=t, y=`50%`,colour=as.factor(federal_guide$abbrev_state[s])))+
  facet_wrap(~federal_guide$abbrev_state[s], scales = "free")+
  scale_colour_viridis_d(name='Region',begin=0.1,end=0.9)+
  scale_fill_viridis_d(name='Region',begin=0.1,end=0.9)+
  labs(x="Time",y=expression(eta['s']),title='Temporal trend in delay', legend=NULL)+
  theme(text = element_text(size=12), plot.title = element_text(size=16),plot.subtitle=element_text(size=14),
        legend.title = element_text(size=14),strip.text = element_text(size=14))
ggsave(eta_comp_plot,file='Plots/full_eta.pdf',width=9,height=3)


zeta_comp_plot<-ggplot()+
  geom_ribbon(data=zeta_matrix, 
              aes(x=t,  ymin=`2.5%`,ymax=`97.5%`,fill=federal_guide$abbrev_state[s]),alpha=0.2)+
  geom_line(data=zeta_matrix,aes(x=t, y=`50%`,colour=federal_guide$abbrev_state[s]))+
  facet_wrap(~federal_guide$abbrev_state[s], scales = "free")+
  scale_colour_viridis_d(name='Region',begin=0.1,end=0.9)+
  scale_fill_viridis_d(name='Region',begin=0.1,end=0.9)+
  labs(x="Time",y=expression(zeta['s']),title='Temporal trend in total counts', legend=NULL)+
  theme(text = element_text(size=12), plot.title = element_text(size=16),plot.subtitle=element_text(size=14),
        legend.title = element_text(size=14),strip.text = element_text(size=14))
ggsave(zeta_comp_plot,file='Plots/full_zeta.pdf',width=9,height=3)


stripe_indicator=function(x=double(1)){
  stripe<-x%in%c(1:50,100:150,200:250,300:350,400:450)
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
         scale_colour_discrete(name='Region') +
  scale_fill_manual(values = c("white", "grey50"), name=NULL) +
  geom_point( aes(x=t,y=y))+facet_wrap(m~federal_guide$abbrev_state[s],scales="free",nrow=2)+
  geom_ribbon( aes(x=t,  ymin=`2.5%`,ymax=`97.5%`),alpha=0.2)+
  geom_line(aes(x=t, y=`50%`,colour=as.factor(federal_guide$abbrev_state[s]))))

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

