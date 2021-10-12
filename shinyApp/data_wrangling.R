

#rm(list=ls())
#packages####
pacman::p_load(dplyr, tidyverse, zoo, shiny, leaflet, RColorBrewer, xts, rgdal, here)

#OWID data#####

# owid_Data<-read.csv(
#     "https://raw.githubusercontent.com/owid/covid-19-data/master/public/data/owid-covid-data.csv",#    fileEncoding="UTF-8-BOM", stringsAsFactors = FALSE
#     ) %>%
#     select(
#       iso_code, location, date, total_cases_per_million, new_cases_per_million, total_vaccinations_per_hundred,
#       stringency_index, population
#     )

##clean up OWID data####
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

#GISAID data#####

here::here()
gisaid <- read.csv(here::here("../data/gisaid_metadata_290621.csv"), header = F) %>%
  dplyr::select(V2, V10, V5, V8, V11) 

column_headings <- c("gisaid_ID", "country", "lineage", "date", "continent")
names(gisaid)<-column_headings

##clean up GISAID data####
gisaid$date<-as.Date(gisaid$date)

gisaid <- gisaid %>% 
  mutate(
    #replace underscore with space & change continent to sentence case
    continent = continent %>% gsub("_", " ", .) %>% 
      str_to_title(gisaid$continent), 
    country = country %>% gsub("_", " ", .) %>% 
      str_to_title(gisaid$continent),
    #options in time intervals
    date2 = as.Date(cut(gisaid$date,
                        breaks = "2 weeks",start.on.monday = FALSE)),
    date3 = as.Date(cut(gisaid$date,
                        breaks = "1 month",start.on.monday = FALSE))
  )



