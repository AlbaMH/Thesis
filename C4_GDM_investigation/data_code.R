# Downloading and formatting NHS England covid deaths data

#DOWNLOAD DATA
library(readxl)
library(openxlsx)
library(lubridate)
library(tidyverse)
library(mgcv)
library(coda)
library(stringr)
library(abind)


# function to get all url link from a webpage 
all_urls <- function(url){
  # Create an html document from the url
  webpage <- xml2::read_html(url)
  # Extract the URLs
  url_ <- webpage %>%
    rvest::html_nodes("a") %>%
    rvest::html_attr("href") 
  return(tibble(url = url_))
}

# get urls from archive data page
archive_urls<-all_urls("https://www.england.nhs.uk/statistics/statistical-work-areas/covid-19-daily-deaths/covid-19-daily-announced-deaths-archive/")
# get url from latest data page
latest_urls<-all_urls("https://www.england.nhs.uk/statistics/statistical-work-areas/covid-19-daily-deaths/")
# combine url's from both pages 
all_pages_urls<-rbind(latest_urls,archive_urls)

# Specify URL where file is stored
nchar_start<-nchar('https://www.england.nhs.uk/statistics/wp-content/uploads/sites/2/2020/04/')
nchar_end<-nchar("https://www.england.nhs.uk/statistics/wp-content/uploads/sites/2/2021/09/COVID-19-daily-announced-deaths-")
is_DAD<-sapply(all_pages_urls,function(x)substr(x,nchar_start,nchar_end))=='/COVID-19-daily-announced-deaths-'
all_urls_data<-all_pages_urls[is_DAD]

# get dates of each data 
urls_nchar<-nchar(all_urls_data)
urls_dates<-as.character(substr(all_urls_data,nchar_end+1,urls_nchar-5))

# clean dates
is_version<-substr(urls_dates,nchar(urls_dates)-1,nchar(urls_dates)-1)=='-'
urls_dates[is_version]<-substr(urls_dates[is_version],1,nchar(urls_dates)[is_version]-2)
is_revised<-substr(urls_dates,nchar(urls_dates),nchar(urls_dates))=='d'
urls_dates[is_revised]<-substr(urls_dates[is_revised],1,nchar(urls_dates)[is_revised]-8)

# format dates for file name
urls_dates_file<-format(as.Date(urls_dates, format = "%d-%B-%Y"),"%d_%b_%y")
n_files<-length(urls_dates)
destfile<-numeric(n_files)
# Specify destination where file should be saved
for(i in 1:n_files){
  destfile[i] <- paste0("~/Dropbox/Covid delay Biometrics revisions/data download/DailyAnnouncedDeaths_",urls_dates_file[i],".xlsx")
}
# Apply download.file function in R
#for(i in 1:n_files){
#  download.file(all_urls_data[i], destfile[i])
#}

# Read file page required
for(i in 1:n_files){
  assign(paste0("DAD_",urls_dates_file[i]),tibble(openxlsx::read.xlsx(paste(all_urls_data[i]),sheet='Tab1 Deaths by region',detectDates=TRUE))[10:18,])
  if(get(paste0("DAD_",urls_dates_file[i]))[9,1]=="South East"){
    assign(paste0("DAD_",urls_dates_file[i]),tibble(openxlsx::read.xlsx(paste(all_urls_data[i]),sheet='Tab1 Deaths by region',detectDates=TRUE))[c(10:11,13:19),])
  }
  m=i
}

for(i in n_files:1){
assign(paste0("DAD_",urls_dates_file[i]),tibble(openxlsx::read.xlsx(paste(all_urls_data[i]),sheet='COVID19 daily deaths by region',detectDates=TRUE,colNames = FALSE))[c(11,12,14:20),])
n=i
}

for(i in (m+1):(n-1)){
  assign(paste0("DAD_",urls_dates_file[i]),tibble(openxlsx::read.xlsx(paste(all_urls_data[i]),sheet='Tab1 Deaths by region',detectDates=TRUE,colNames = FALSE))[c(11,12,14:20),])
}

# FORMAT DATA
#remove NA's
min_date<-numeric(0)
col.na<-list()
for(i in  1:n_files){
 #columns with no values  
 col.na[[i]]<-which(is.na(get(paste0("DAD_",urls_dates_file[i]))[1,]))
 #remove columns
 assign(paste0("DAD2_",urls_dates_file[i]),get(paste0("DAD_",urls_dates_file[i]))[,-col.na[[i]]])
 #first death in each report 
 min_date=c(min_date,as.numeric(get(paste0("DAD2_",urls_dates_file[i]))[1,2]))
 }

# calculate first and last death
first_death_days<-min(na.omit(min_date))
latest_death_days<-as.numeric(get(paste0("DAD2_",urls_dates_file[1]))[1,ncol(get(paste0("DAD2_",urls_dates_file[1])))-2])
origin<-as.Date(urls_dates[1], format="%d-%B-%Y")-1-latest_death_days
latest_death<-origin+latest_death_days
first_death<-origin+first_death_days

DoD_n<-as.numeric(latest_death-first_death)+1
DA_n<-as.numeric(as.Date(urls_dates[1], format="%d-%B-%Y")-as.Date(urls_dates[n_files], format="%d-%B-%Y"))+1
S<-nrow(get(paste0("DAD2_",urls_dates_file[1])))-1
DA_dates<-seq(from=as.Date(urls_dates[n_files], format="%d-%B-%Y")-1, to=as.Date(urls_dates[1], format="%d-%B-%Y")-1, by='days')
DoD_dates<-seq(from=first_death, to=latest_death, by='days')  

data_matrix<-matrix(0, nrow=DoD_n,ncol=DA_n)
data_frame<-as_tibble(data_matrix)
colnames(data_frame)<-DA_dates
rownames(data_frame)<-DoD_dates

n_regions<-nrow(get(paste0("DAD2_",urls_dates_file[1])))-1
regions_compact<-numeric(n_regions)
for(i in 1:n_regions){
  regions_compact[i]<-gsub(" ", "",get(paste0("DAD2_",urls_dates_file[1]))[i+1,1])
  assign(paste0(regions_compact[i]),data_frame)
}

# pub data in table
for(i in 1:n_files){
  table<-get(paste0("DAD2_",urls_dates_file[i]))
  table_ncol<-ncol(table)
  if(table_ncol>3){
    for(c in 2:(table_ncol-2)){
      DOD<-origin+as.numeric(table[1,c])
      for(r in 2:(n_regions+1)){
        new_frame<-get(paste0(regions_compact[r-1]))
        new_frame[which(DoD_dates==DOD),which(DA_dates==(as.Date(urls_dates[i], format="%d-%B-%Y")-1))]<-as.numeric(table[r,c])
        rownames(new_frame)<-DoD_dates
        assign(paste0(regions_compact[r-1]), new_frame)
      }
    }
  }
}

#comibine all regions into array
all_regions<-get(regions_compact[1])
for (r in 2:n_regions) {
  all_regions<-abind(all_regions,get(regions_compact[r]),along=3)
}

#regions_raw for COVID_Master
regions_raw<-all_regions[-(1:32),-1,]
save(regions_raw, file = "regions_raw.Rdata")


