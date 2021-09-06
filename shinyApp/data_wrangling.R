
#rm(list=ls())

pacman::p_load(dplyr, tidyverse, zoo, shiny, leaflet, RColorBrewer, xts, rgdal)

# owid_Data<-read.csv(
#     "https://raw.githubusercontent.com/owid/covid-19-data/master/public/data/owid-covid-data.csv",#    fileEncoding="UTF-8-BOM", stringsAsFactors = FALSE
#     ) %>%
#     select(
#       iso_code, location, date, total_cases_per_million, new_cases_per_million, total_vaccinations_per_hundred,
#       stringency_index, population
#     )

owid_Data$date<-as.Date(owid_Data$date)
owid_Data$date7<-as.Date(cut(owid_Data$date,breaks = "1 week",start.on.monday = FALSE))


owid_Data[c(
  "total_cases_per_million", "new_cases_per_million", "total_vaccinations_per_hundred"
)][is.na(owid_Data[c(
  "total_cases_per_million", "new_cases_per_million", "total_vaccinations_per_hundred")])] <- 0

#drop anything that starts with OWID_
removeList<-unique(owid_Data$iso_code[startsWith(owid_Data$iso_code, "OWID_")]) 

owid_Data1 <- owid_Data %>%
  subset(., ! iso_code %in% removeList)%>% 
  dplyr::group_by(iso_code) %>% 
  dplyr::arrange(date) %>% 
  dplyr::mutate(
    smoothed_new_cases_per_million = zoo::rollmean(new_cases_per_million, k = 14, fill = NA)) %>%   #biweekly rolling mean - for smoothing
  dplyr::mutate(
    smoothed_total_cases_per_million = zoo::rollmean(total_cases_per_million, k = 14, fill = NA)) %>%   #biweekly rolling mean - for smoothing
  dplyr::ungroup() 
# %>%
#   select(-c(total_cases_per_million, new_cases_per_million, population))

head(owid_Data)
