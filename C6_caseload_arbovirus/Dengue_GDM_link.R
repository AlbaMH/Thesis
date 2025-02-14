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


#setwd("~/OneDrive/PhD/Dengue 2023")
#setwd("~/OneDrive - University of Glasgow/PhD/Dengue 2023")
#setwd("~/Documents/Dengue 2023")

states<-read_state(year=2013)


data_set <- read_csv("data set.csv")
data_set<-data_set[,-c(22:48)]


data<-data_set%>%mutate(start_of_symptoms=as.Date((DT_SIN_PRI),format="%Y-%m-%d"),
                        notification_date=as.Date((DT_NOTIFIC),format="%Y-%m-%d"))%>%
  mutate(delay=difftime(notification_date,start_of_symptoms, unit=c("days")))%>%
  dplyr::select(start_of_symptoms,notification_date,delay,CLASSI_FIN,SG_UF_NOT)%>%
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

median(Dengue_z[,,1]/Dengue_y, na.rm=TRUE)
max(Dengue_z[,,2]/Dengue_y, na.rm=TRUE)
max(Dengue_z[,,3]/Dengue_y, na.rm=TRUE)
median(Dengue_z[,,4]/Dengue_y, na.rm=TRUE)

# Define beta binomial distribution.
rbetabin=nimbleFunction(run=function(n=integer(0),mu=double(0),phi=double(0),size=double(0)){
  pi=rbeta(n,mu*phi,(1-mu)*phi)
  returnType(double(1))
  return(rbinom(n,size,pi))
})




Dengue_surv_link_region <- function(seed, Dengue_z, population_order, Nowcast_list, S,  D, D_max, n_knots, n_iter, n_burn, n_thin, n_chains){
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
    for(s in 1:S){
      for(t in 1:N){
        # Negative Binomial Model for total Dengue cases
        y[s,t] ~ dnegbin(prob=theta[s]/(theta[s]+lambda[t,s]),size=theta[s])  
        
        # Model for delayed Dengue counts
        z[s,t,1] ~ dbetabin(nu[t,1,s],phi[s,1],y[s,t])
      }
      for(d in 2:D){
        for(t in 1:obs_index[s,d]){
          z[s,t,d] ~ dbetabin(nu[t,d,s],phi[s,d],y[s,t]-sum(z[s,t,1:(d-1)]))
        }
      }
      
      for(t in 1:N){
        # Linear predictor (Beta-Binomial means).
       # index_delta[t,s]<-delta_converter(x=(log((y[s,t]+1)/pop[s])),grid =delta_grid[1:length_delta], grid_length = length_delta)[1]
        for(d in 1:D){
          # Expected cumulative proportions.
          probit(p[t,d,s]) <- psi[d,s] + eta[t,s] + delta[s]*log((y[s,t]+1)/pop[s])
        }
        # Relative proportions (Beta-Binomial means).
        nu[t,1,s] <- p[t,1,s]
        for(d in 2:D){
          nu[t,d,s] <- (p[t,d,s]-p[t,d-1,s])/(1-p[t,d-1,s])
        }
        # Mean for total Dengue cases 
        log(lambda[t,s]) <-  log(pop[s]) + int_zeta[s] + zeta[t,s] + xi[weeks[t],s] 
        
      }
      
      
      # logsplines of time in cumulative proportion reported.
      Omega_eta[1:K_t,1:K_t,s] <- S_t[1:K_t,1:K_t]/sigma_eta[s]^2
      kappa_eta[1:K_t,s] ~ dmnorm(zeros[1:K_t],Omega_eta[1:K_t,1:K_t,s])
      eta[1:N,s] <- X_t[1:N,1:K_t]%*%kappa_eta[1:K_t,s]
      
      # Seasonal logspline 
      Omega_xi[1:K_w,1:K_w,s]<-S_w[1:K_w,1:K_w]/sigma_xi[s]^2
      kappa_xi[1:K_w,s]~dmnorm(zeros[1:K_w],Omega_xi[1:K_w,1:K_w,s])
      xi[1:52,s]<-X_w[1:52,1:K_w]%*%kappa_xi[1:K_w,s]
      
      # logsplines of time for Dengue cases.
      Omega_zeta[1:K_t,1:K_t,s] <- S_t[1:K_t,1:K_t]/sigma_zeta[s]^2
      kappa_zeta[1:K_t,s] ~ dmnorm(zeros[1:K_t],Omega_zeta[1:K_t,1:K_t,s])
      zeta[1:N,s] <- X_t[1:N,1:K_t]%*%kappa_zeta[1:K_t,s]
      
      ##for(d in 1:D){
    #  Omega_delta[1:K_c,1:K_c,s] <- S_c[1:K_c,1:K_c,s]/sigma_delta[1]^2 + S_c[1:K_c,(K_c+1):(2*K_c),s]/sigma_delta[2]^2 #non-linear plus linear spline 
     # kappa_delta[1:K_c,s] ~ dmnorm(zeros[1:K_c],Omega_delta[1:K_c,1:K_c,s])
     # delta[1:length_delta,s] <- X_c[1:length_delta,1:K_c,s]%*%kappa_delta[1:K_c,s]

      
      #for(k in 1:2){
      # sigma_delta[k,s] ~ T(dnorm(0,1),0,) 
      ## }
      # }
      
      # Smoothing parameter priors.
      sigma_zeta[s] ~  T(dnorm(0,sd=10),0,)
      sigma_eta[s] ~  T(dnorm(0,sd=10),0,)
      sigma_xi[s] ~  T(dnorm(0,sd=10),0,)
      
      # Cumulative delay curves.
      psi[1,s] ~ dnorm(0,sd=5)
      for(d in 2:D){
        psi[d,s] ~ T(dnorm(psi[d-1,s],sd=5),psi[d-1,s],)
      }
      
      # Beta-Binomial dispersion parameters.
      for(d in 1:D){
        phi[s,d] ~ dgamma(2,0.02) 
      }
      
      # Intercepts
      int_zeta[s]  ~ dnorm(-10,sd=10)
      theta[s] ~ dgamma(2,0.02) # Negative-Binomial dispersion parameters.
      
      # Linear link between case load and delay. 
      delta[s] ~ dnorm(0,5)

    }
      
    
  })
  
  # Data for current nowcast date 
  N_now<-Nowcast_list[seed]
  N<-N_now
  C<-N-D_max
  z<- Dengue_z[1:S,(1):(N_now),]
  
  # Reported partial counts:
  censored_z <- z
  
  for(s in 1:S){
    censored_z[s,,][outer(1:dim(z[s,,])[1], 0:(dim(z[s,,])[2]-1), FUN = "+") > N] <- NA
  }
  
  # Maximum observed index for each region at each delay.
  obs_index<-matrix(NA, nrow=S, ncol=D)
  for(s in 1:S){
    for(d in 1:D){
      obs_index[s,d]<-which(is.na(censored_z[s,,d])==TRUE)[1]-1
    }
  }
  
  
  # Censor the totals.
  censored_y <- apply(censored_z[,,1:D_max],c(1,2),sum)
  week_index<-rep(1:52,9) #first date is 01-01-2013
  # Set up logspline bases.
  y_spline<-c(apply(censored_y[,1:(N-D_max)],2,mean),apply(censored_z[,(N-D_max+1):N,],2,sum,na.rm=T))
  blank_data=tibble(y=rnorm(N,0,1),t=1:N,w=week_index[1:N]) 
  blank_jagam=jagam(y~s(t,k=floor(N/20), bs='cs',pc=round(C/2))+
                      s(w,bs='cc',k=n_knots[2]),data=blank_data,file='blank.jags',
                    knots=list(t=seq(1,N,length=floor(N/20)),w=seq(0,52,length=n_knots[2]))) # #
  

  
  
  # Constants
  Dengue_constants <- list(N=N,
                           S=S,
                           D=D,
                           weeks=week_index[1:N],
                           pop=as.numeric(population_order$population),
                           obs_index=(obs_index),
                           K_w=dim(blank_jagam$jags.data$S2)[1],
                           S_w=blank_jagam$jags.data$S2,
                           K_t=dim(blank_jagam$jags.data$S1)[1],
                           S_t=blank_jagam$jags.data$S1)
  
 
  Dengue_constants$X_t <- blank_jagam$jags.data$X[,2:(Dengue_constants$K_t+1)]
  Dengue_constants$X_w <- blank_jagam$jags.data$X[,(Dengue_constants$K_t+2):(Dengue_constants$K_t+Dengue_constants$K_w+1)]
  Dengue_constants$zeros <- rep(0,max(Dengue_constants$K_t,Dengue_constants$K_w)) 
  
  # Data 
  Dengue_data <- list(y=censored_y[,1:N], z=censored_z[,,1:D]) 
  
  # Generate random initial values.
  Dengue_inits <- list(
    kappa_zeta=matrix(rnorm(S*(Dengue_constants$K_t),0,0.1),ncol=S),
    kappa_eta=matrix(rnorm(S*(Dengue_constants$K_t),0,0.1),ncol=S),
    kappa_xi=matrix(rnorm(S*(Dengue_constants$K_w),0,0.1),ncol=S),
    # kappa_delta=rnorm((Dengue_constants$K_c),0,0.1),
    delta=rnorm(S,0,0.1),#matrix(rnorm(D*S,0,1),nrow=D),
 #   kappa_delta=matrix(rnorm(S*(Dengue_constants$K_c),0,0.1),ncol=S),
    phi=matrix(abs(rnorm(S*D,30,10)),nrow=S,ncol=D),
    psi=matrix(sort(rnorm(D*S,0,0.1)),ncol=S, byrow = TRUE),
    sigma_eta=runif(S,0,1),
    sigma_zeta=runif(S,0,1),
    sigma_xi=runif(S,0,1),
#    sigma_delta=runif(2,0,1),#matrix(runif(2*S,0,1),nrow=2),
    theta=abs(rnorm(S,50,10)),
    int_zeta=rnorm(S,-10,1),
    y=matrix(NA,nrow=S,ncol=N) )
  
  for(t in 1:N){
    for(s in 1:S){
      if(is.na(censored_y[s,t])) Dengue_inits$y[s,t]=sum(Dengue_data$z[s,t,],rpois(1,median(Dengue_data$y[s,]-rowSums(Dengue_data$z[s,,]),na.rm=T)),na.rm=TRUE)
    }
  }
  
  # Build the model.
  Dengue_model <- nimbleModel(Dengue_code,Dengue_constants,Dengue_data,Dengue_inits)
  # Compile the model.
  Dengue_compiled_model <- compileNimble(Dengue_model)
  
  # Set up the MCMC.
  Dengue_mcmc_config <- configureMCMC(Dengue_model,monitors=c("lambda","int_zeta","y","xi",
                                                              "psi","theta","phi", 
                                                              "sigma_eta","sigma_zeta",
                                                              "eta","zeta","delta",
                                                              "kappa_zeta","kappa_eta"),
                                      useConjugacy = FALSE)
  
  Dengue_mcmc_config$removeSamplers('kappa_zeta','sigma_zeta','int_zeta',
                                    'kappa_eta','sigma_eta',"psi",
                                     'delta',
                                    'kappa_xi','sigma_xi')
  for(s in 1:S){
    Dengue_mcmc_config$addSampler(target=c(paste('kappa_zeta[1:',(Dengue_constants$K_t),', ',s,']', sep=''),paste('sigma_zeta[',s,']',sep=''),
                                           paste('kappa_xi[1:',(Dengue_constants$K_w),', ',s,']', sep=''),paste('sigma_xi[',s,']',sep=''),
                                           paste('int_zeta[',s,']',sep='')),type='AF_slice') 
    Dengue_mcmc_config$addSampler(target=c(paste('delta[',s,']', sep=''),
                                           paste('kappa_eta[1:',(Dengue_constants$K_t),', ',s,']', sep=''),paste('sigma_eta[',s,']',sep=''),
                                           paste('psi[1:',D,', ',s,']', sep='')),type='AF_slice')
  }
  
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

# lower_delta= min(log((1/population_order$population)))
# upper_delta=  max(log((Dengue_y+1)/population_order$population))+4
# by_delta=0.1
# length_delta=1000
# delta_grid<-seq(from=lower_delta, to=upper_delta, length=length_delta)
# grid_length<-length(delta_grid)
Nowcast_list=floor(seq(from=2*52, to=375, length=4))

population_order$population_scale<-scale(population_order$population)
# Number of available computer cores 
n_cores<-min(detectCores(),length(Nowcast_list))
# Make Cluster for MCMC. 
this_cluster_link<- makeCluster(n_cores)



# Note: takes approximately 
time_link<-system.time({
  Dengue_surv_link <- parLapply(cl = this_cluster_link, X = 1:length(Nowcast_list), 
                                      fun = Dengue_surv_link_region,
                                      Dengue_z=as.array(Dengue_z)[1:S,,],
                                      population_order=population_order,
                                      Nowcast_list=Nowcast_list,#
                                      S=S,
                                      D=D,
                                      D_max=D_max,
                                      n_iter=n_iter,
                                      n_burn=n_burn,
                                      n_thin=n_thin,
                                      n_chains=n_chains,
                                      n_knots=n_knots)
})


stopCluster(this_cluster_link)

#SAVE Samples
save(Dengue_surv_link,time_link,file='Dengue_surv_link_D4.RData')
#


plot(Dengue_surv_link[[2]][,"zeta[3, 1]"],main="zeta[3, 1]")
plot(Dengue_surv_link[[2]][,"int_zeta[1]"],main="int_zeta[1]")
plot(Dengue_surv_link[[2]][,"zeta[100, 1]"],main="zeta[100, 1]")
plot(Dengue_surv_link[[2]][,"xi[10, 1]"],main="xi[10, 1]")
plot(Dengue_surv_link[[2]][,"xi[1, 1]"],main="xi[1, 1]")
plot(Dengue_surv_link[[2]][,"delta[5]"],main="delta[5]")
plot(Dengue_surv_link[[2]][,"delta[3]"],main="delta[3]")
plot(Dengue_surv_link[[1]][,"delta[1]"],main="delta[1]")
plot(Dengue_surv_link[[1]][,"delta[2]"],main="delta[2]") 
plot(Dengue_surv_link[[1]][,"delta[5]"],main="delta[5]")
plot(Dengue_surv_link[[1]][,"delta[6]"],main="delta[6]")
plot(Dengue_surv_link[[1]][,"psi[1, 1]"],main="psi[1, 1]")
plot(Dengue_surv_link[[1]][,"eta[3, 1]"],main="eta[1, 1]")

plot(Dengue_surv_link[[2]][,"psi[2, 1]"],main="psi[2, 1]")
plot(Dengue_surv_link[[2]][,"psi[3, 1]"],main="psi[3, 1]")
plot(Dengue_surv_link[[2]][,"psi[4, 1]"],main="psi[4, 1]")
plot(Dengue_surv_link[[2]][,"psi[4, 1]"],main="psi[4, 1]")
plot(Dengue_surv_link[[2]][,"psi[2, 2]"],main="psi[2, 2]")


#min(dplyr::select(as_tibble(do.call('rbind',as.mcmc.list(Dengue_output_cluster[[1]]))),starts_with('zeta')))
##plot(Dengue_output_cluster[[2]][,"int_zeta[3]"],main="int_zeta[1]")
#plot(Dengue_output_cluster[[2]][,"zeta[104, 1]"],main="zeta[1]")


Dengue_parameter_group=c("lambda","int_zeta","delta",
                         "psi","theta","phi",
                         "sigma_eta","sigma_zeta", 
                         "zeta", "kappa_zeta","kappa_eta")
Dengue_psrf<-list()
Dengue_mean_psrf<-Dengue_max_psrf<-Dengue_leq_105 <-Dengue_leq_110<-array(NA,dim=c(length(Dengue_parameter_group),length(Nowcast_list)))
for(k in 1:length(Dengue_parameter_group)){
  Dengue_psrf[[k]]<-list()
  for(j in 1:length(Nowcast_list)){
    Dengue_parameter_names <- colnames((as.matrix(Dengue_surv_link[[j]])))
    Dengue_index <- which(startsWith(Dengue_parameter_names,Dengue_parameter_group[k]))
    Dengue_psrf[[k]][[j]]<- sapply(Dengue_index,function(v)gelman.diag(Dengue_surv_link[[j]][,v],autoburnin = FALSE,transform = TRUE,multivariate = FALSE)$psrf[,1])
    Dengue_mean_psrf[k,j] <- mean(Dengue_psrf[[k]][[j]])
    Dengue_max_psrf[k,j] <- max(Dengue_psrf[[k]][[j]])
    Dengue_leq_105[k,j] <- mean(Dengue_psrf[[k]][[j]]<=1.05)
    Dengue_leq_110[k,j] <- mean(Dengue_psrf[[k]][[j]][j]<=1.1)
  }
}
Dengue_MAX_psrf<-cbind(Dengue_parameter_group,apply(Dengue_max_psrf,1,max))
Dengue_MEAN_psrf<-cbind(Dengue_parameter_group,apply(Dengue_mean_psrf,1,mean))





#S=3
# Observed Dengue hospitalisations.
y_all<-melt(t(Dengue_y[1:S,]))
colnames(y_all)<-c("t","s","y")
y_all$s<-federal_guide$abbrev_state[y_all$s]

# Set up plots:
N_max<-dim(Dengue_y)[2]

lambda_plot<-ggplot()+
  geom_point(data=y_all,aes(x=t,y=y))+
  facet_wrap(~s, scales = 'free')+
  theme_minimal()+labs(y="Mean Dengue Hospitalisations", x="Time", title="CLR Regional Model")+
  scale_colour_discrete(name='Nowcast Date')+scale_fill_discrete(name='Nowcast Date')

y_plot<-ggplot()+
  geom_point(data=y_all,aes(x=t,y=y))+ 
  facet_wrap(~s, scales = 'free')+
  theme_minimal()+labs(y="Dengue Hospitalisations", x="Time", title="CLR Regional Model")+#+geom_vline(xintercept = 51)
  scale_colour_discrete(name='Nowcast Date')+scale_fill_discrete(name='Nowcast Date')

delta_plot<-ggplot()+  theme_minimal()+
  scale_colour_viridis_c(name='Nowcast Date',begin=0.1,end=0.9, option="magma")+
  labs(x="Brazilian federative unit (s)",y=expression(delta['s']),
       title='Brazil dengue hospitalisations',
       subtitle='Linear effect of case load on reporting delay from rolling prediction experiment', legend=NULL)+
  theme(text = element_text(size=12), plot.title = element_text(size=16),plot.subtitle=element_text(size=14),
        legend.title = element_text(size=14),strip.text = element_text(size=14), 
        legend.position = "bottom",legend.box="vertical")

delta_plot2<-ggplot()+
  theme_minimal()+labs(y="Effect size (probit scale)", x="lambda (expected mean dengue counts)", title="Effect of dengue incidence on delay",
                       color = "Nowcast date") + theme(legend.position = "bottom")+
  facet_wrap(~s,scales='free')+
  scale_colour_discrete(name='Nowcast Date')+scale_fill_discrete(name='Nowcast Date')

omega_plot<-ggplot()+
  facet_wrap(~s, scales = 'free')+
  theme_minimal()+labs(y="zeta", x="Time", title="Temporal logspline for mean dengue cases")+
  scale_colour_discrete(name='Nowcast Date')+scale_fill_discrete(name='Nowcast Date')

zeta_plot<-ggplot()+
  facet_wrap(~s, scales = 'free')+
  theme_minimal()+labs(y="zeta", x="Time", title="Temporal logspline for mean dengue cases")
xi_plot<-ggplot()+
  facet_wrap(~s, scales = 'free')+
  theme_minimal()+labs(y="xi", x="Weeks", title="Seasonal logspline for mean dengue cases")
psi_plot<-ggplot()+
  facet_wrap(~s, scales = 'free')+
  theme_minimal()+labs(y="psi", x="delay", title="Delay trend on cumulative proportions")

Dengue_y_samples<-list()
lambda<-lambda_raw<-list()
Dengue_pred<-list()
Dengue_link_quantiles<-list()
theta<-delta<-zeta<-xi<-psi <- list()
delta<-delta_cubic<-delta_quantiles<-omega<-list()
# Explore output
Dengue_logy<-list()
Dengue_combined_logy<-list()
for(j in length(Nowcast_list):1){
  Dengue_logy[[j]]<-as.mcmc.list(Dengue_surv_link[[j]])
  # Combine all MCMC chains.
  Dengue_combined_logy[[j]] <- as_tibble(do.call('rbind',Dengue_logy[[j]]))
}

n_sim <- dim(Dengue_combined_logy[[1]])[1]
rbetabin=nimbleFunction(run=function(n=integer(0),mu=double(0),phi=double(0),size=double(0)){
  pi=rbeta(n,mu*phi,(1-mu)*phi)
  returnType(double(0))
  return(rbinom(n,size,pi))
})
week_index<-rep(1:52,9)
# plot output
Dengue_logy_raw<-Dengue_logy_rug<-list()
for(j in length(Nowcast_list):1){
  # Plot Model Variables 
  lambda[[j]] <- dplyr::select(Dengue_combined_logy[[j]],starts_with('lambda'))%>%as.matrix()%>%array(dim=c(n_sim,Nowcast_list[j],S))%>%
    apply(c(2,3),quantile,c(0.025,0.5,0.975))%>%melt(varnames=c('quantile','t','s'),value.name='y')%>%#mutate(t=t+Nowcast_list[j]-W)%>%
    spread(quantile,y)%>%mutate(Nowcast_date=paste(Nowcast_list[j]),s=federal_guide$abbrev_state[s])
  
  lambda_plot<-lambda_plot+
    geom_line(data=lambda[[j]], aes(x=t, y=`50%`, colour=Nowcast_date))+
    geom_ribbon(data=lambda[[j]], aes(x=t, ymin=`2.5%`, ymax=`97.5%`,fill=Nowcast_date), alpha=0.2)
  
  
  Dengue_link_quantiles[[j]]<-dplyr::select(Dengue_combined_logy[[j]],starts_with('y'))%>%as.matrix()%>%array(dim=c(n_sim,S,Nowcast_list[j]))%>%
    apply(c(2,3),quantile,c(0.025,0.5,0.975))%>%melt(varnames=c('quantile','s','t'),value.name='y')%>%#mutate(t=t+Nowcast_list[j]-W)%>%
    spread(quantile,y)%>%mutate(Nowcast_date=paste(Nowcast_list[j]),s=federal_guide$abbrev_state[s])
  
  y_plot<-y_plot+
    geom_line(data=Dengue_link_quantiles[[j]], aes(x=t, y=`50%`, colour=Nowcast_date))+
    geom_ribbon(data=Dengue_link_quantiles[[j]], aes(x=t, ymin=`2.5%`, ymax=`97.5%`,fill=Nowcast_date))
  
  
  # Plot Zeta Variables 
  zeta[[j]] <- dplyr::select(Dengue_combined_logy[[j]],starts_with('zeta'))%>%as.matrix()%>%array(dim=c(n_sim,Nowcast_list[j],S))%>%
    apply(c(2,3),quantile,c(0.025,0.5,0.975))%>%melt(varnames=c('quantile','t','s'),value.name='y')%>%#mutate(t=t+Nowcast_list[j]-W)%>%
    spread(quantile,y)%>%mutate(Nowcast_date=paste(Nowcast_list[j]),s=federal_guide$abbrev_state[s])
  
  zeta_plot<-zeta_plot+
    geom_line(data=zeta[[j]], aes(x=t, y=`50%`, colour=Nowcast_date))+
    geom_ribbon(data=zeta[[j]], aes(x=t, ymin=`2.5%`, ymax=`97.5%`,fill=Nowcast_date), alpha=0.2)
  
  xi[[j]] <- dplyr::select(Dengue_combined_logy[[j]],starts_with('xi'))%>%as.matrix()%>%array(dim=c(n_sim,52,S))%>%
    apply(c(2,3),quantile,c(0.025,0.5,0.975))%>%melt(varnames=c('quantile','w','s'),value.name='y')%>%#mutate(t=t+Nowcast_list[j]-W)%>%
    spread(quantile,y)%>%mutate(Nowcast_date=paste(Nowcast_list[j]),s=federal_guide$abbrev_state[s])
  
  xi_plot<-xi_plot+
    geom_line(data=xi[[j]], aes(x=w, y=`50%`, colour=Nowcast_date))+
    geom_ribbon(data=xi[[j]], aes(x=w, ymin=`2.5%`, ymax=`97.5%`,fill=Nowcast_date), alpha=0.2)
  

  
  delta[[j]] <- dplyr::select(Dengue_combined_logy[[j]],starts_with('delta'))%>%as.matrix()%>%array(dim=c(n_sim,S))%>%
    apply(c(2),quantile,c(0.025,0.5,0.975))%>%melt(varnames=c('quantile','s'),value.name='y')%>%
    spread(quantile,y)%>%mutate(Nowcast_date=as.numeric(paste(Nowcast_list[j])),s=federal_guide$abbrev_state[s])#%>%filter(s%in%fu_abbr)
  delta_plot<-delta_plot+geom_point(data=delta[[j]], aes(x=s, y=`50%`,colour=Nowcast_date))+
    geom_errorbar(data=delta[[j]], aes(x=s, ymin=`2.5%`, ymax=`97.5%`,colour=Nowcast_date),alpha=0.75)
  
  #omega[[j]] <- dplyr::select(Dengue_combined_logy[[j]],starts_with('omega'))%>%as.matrix()%>%array(dim=c(n_sim,W,S))%>%
  #  apply(c(2,3),quantile,c(0.025,0.5,0.975))%>%melt(varnames=c('quantile','t','s'),value.name='y')%>%mutate(t=t+Nowcast_list[j]-W)%>%
  # spread(quantile,y)%>%mutate(Nowcast_date=paste(Nowcast_list[j]),s=federal_guide$abbrev_state[s])
  #  omega_plot<-omega_plot+
  #    geom_line(data=omega[[j]], aes(x=t, y=`50%`, colour=Nowcast_date))+
  #   geom_ribbon(data=omega[[j]], aes(x=t, ymin=`2.5%`, ymax=`97.5%`,fill=Nowcast_date), alpha=0.2)
  
  
  psi[[j]] <- dplyr::select(Dengue_combined_logy[[j]],starts_with('psi'))%>%as.matrix()%>%array(dim=c(n_sim,D,S))%>%
    apply(c(2,3),quantile,c(0.025,0.5,0.975))%>%melt(varnames=c('quantile','d','s'),value.name='y')%>%#mutate(t=t+Nowcast_list[j]-W)%>%
    spread(quantile,y)%>%mutate(Nowcast_date=paste(Nowcast_list[j]),s=federal_guide$abbrev_state[s])
  
  psi_plot<-psi_plot+
    geom_line(data=psi[[j]], aes(x=d, y=`50%`, colour=Nowcast_date))+
    geom_ribbon(data=psi[[j]], aes(x=d, ymin=`2.5%`, ymax=`97.5%`,fill=Nowcast_date), alpha=0.2)
  
}

# SAVE PLOTS
lambda_plot
y_plot
delta_plot
ggsave(delta_plot,file='predictions_delta_dengue.pdf',width=9,height=6)
zeta_plot
xi_plot
psi_plot
#ggsave(y_plot,file='Plots/pi_surv_Dengue.pdf',width=10,height=7)

check_ae_plot<-ggplot()+
  geom_point(data=filter(y_all,t%in%c(371:376)),aes(x=t,y=y, colour=s))+ 
  theme_minimal()+labs(y="Dengue Hospitalisations", x="Time", title="CLR Regional Model")+#+geom_vline(xintercept = 51)
  scale_colour_discrete(name='Nowcast Date')+scale_fill_discrete(name='Nowcast Date')+
  geom_line(data=filter(Dengue_link_quantiles[[4]],t%in%c(371:376)), aes(x=t, y=`50%`, colour=s))+
  facet_wrap(~s)+scale_x_continuous(limits=c(371,376))+
  geom_ribbon(data=filter(Dengue_link_quantiles[[4]],t%in%c(371:376)), aes(x=t, ymin=`2.5%`, ymax=`97.5%`,fill=s), alpha=0.2)

check_ae<-full_join(filter(Dengue_link_quantiles[[4]],t%in%c(373)),filter(y_all,t%in%c(373)))%>%mutate(ae=`50%`-y)

# compare to previous version 
load("~/media/alba/Disk 1/OneDrive/PhD/Dengue 2023/Dengue_surv_sea_D4.RData")
Dengue_sea<-list()
Dengue_combined_sea<-Dengue_sea_quantiles<-list()
for(j in length(Nowcast_list):1){
  Dengue_sea[[j]]<-as.mcmc.list(Dengue_surv_sea[[j]])
  # Combine all MCMC chains.
  Dengue_combined_sea[[j]] <- as_tibble(do.call('rbind',Dengue_sea[[j]]))
  Dengue_sea_quantiles[[j]]<-dplyr::select(Dengue_combined_sea[[j]],starts_with('y'))%>%as.matrix()%>%array(dim=c(n_sim,S,Nowcast_list[j]))%>%
    apply(c(2,3),quantile,c(0.025,0.5,0.975))%>%melt(varnames=c('quantile','s','t'),value.name='y')%>%#mutate(t=t+Nowcast_list[j]-W)%>%
    spread(quantile,y)%>%mutate(Nowcast_date=paste(Nowcast_list[j]),s=federal_guide$abbrev_state[s])
  
  
}

library(ggh4x)
# dengue nowcasts.
dengue_nowcasts_compare<-list()
for(j in 1:length(Nowcast_list)){
  dengue_nowcasts_compare[[j]] <- rbind(
    full_join(Dengue_link_quantiles[[j]],filter(y_all, t%in%c((1):Nowcast_list[j])),by=c('s','t'))%>%mutate(m=c('Dengue delta link'),d=t-Nowcast_list[j],now=j)
   # full_join(Dengue_sea_quantiles[[j]],filter(y_all, t%in%c((1):Nowcast_list[j])),by=c('s','t'))%>%mutate(m=c('Dengue no delta'),d=t-Nowcast_list[j],now=j)
  )%>%filter(d>-D)
}
dengue_nowcast_dmax<-array(unlist(dengue_nowcasts_compare), dim=c(S*D*1,ncol(dengue_nowcasts_compare[[1]]),length(Nowcast_list)))%>%
  apply(c(2),c)%>%as_tibble()
colnames(dengue_nowcast_dmax)<-colnames(dengue_nowcasts_compare[[1]])   
dengue_nowcast_dmax<-dengue_nowcast_dmax%>%mutate(t=as.numeric(t),`2.5%`=as.numeric(`2.5%`),`50%`=as.numeric(`50%`),`97.5%`=as.numeric(`97.5%`),
                                                  m=factor(m,levels=c('Dengue delta link'#,'Dengue no delta'
                                                  )),y=as.numeric(y),
                                                  d=as.numeric(d),now=as.numeric(now))
dengue_nowcast_dmax_summary<-dengue_nowcast_dmax%>%group_by(m,d,s)%>%
  summarise(Coverage=mean(y>=`2.5%`&y<=`97.5%`),
            `Prediction interval width`=mean(`97.5%`-`2.5%`),
            `Mean absolute error`=mean(abs(`50%`-y)))%>%
  pivot_longer(cols=c('Coverage','Mean absolute error','Prediction interval width'),names_to = 'metric')%>%
  mutate(m=factor(m,levels=c('Dengue delta link'#,'Dengue no delta'
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
totals_dengue <- ggplot(dengue_nowcast_dmax_summary_average)+
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
totals_dengue

