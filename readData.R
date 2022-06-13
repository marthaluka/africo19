
# This script reads data and should be run before:
    # lineage_infectionRate_plot.R
    # regression_continent.R
    # regression_country.R


#Read data

#rm(list = ls())

## packages
require("pacman")
pacman::p_load(data.table, tidyverse, zoo, RColorBrewer, directlabels,corrplot, data.table,
               stringr, spData, ggrepel, ggcorrplot, patchwork, cowplot, egg, here, visreg,
               vroom)

# OWID data ######### 
# read data
owidData0<- vroom::vroom(
  "https://raw.githubusercontent.com/owid/covid-19-data/master/public/data/owid-covid-data.csv",
  #stringsAsFactors = FALSE   #for fread
  ) %>%                             
  dplyr::select(continent, location, date, new_cases_per_million,total_vaccinations_per_hundred, 
                stringency_index, new_tests_per_thousand, tests_per_case)

names(owidData0)


head(owidData0)
#owidData0[owidData0=='']=NA #if you used fread above 
owidData0 <- drop_na(owidData0,continent)
head(owidData0)

#dim(owidData0) # Note: use vroom/read.csv not fread for summary stats
length(owidData0$new_tests_per_thousand)
sum(is.na(owidData0$new_tests_per_thousand))  #drop number of tests per capita  due to missing data (>50% entries are NA)

#date formats
owidData0$date<-as.Date(owidData0$date)
owidData0$date7<-as.Date(cut(owidData0$date,breaks = "1 week",start.on.monday = TRUE))


#Pre-processing
#biweekly rolling mean - for smoothing

owidData1 <- owidData0 %>%
  dplyr::mutate(across(where(is.character),str_trim)) %>% 
  dplyr::filter(new_cases_per_million>=0) %>%    #drop any negative values 
  dplyr::group_by(location) %>% 
  dplyr::arrange(date) %>% 
  dplyr::mutate(smoothed_cases_per_million = zoo::rollmean(new_cases_per_million, k = 14, fill = 0)) %>%  
  dplyr::mutate(smoothed_vaccines_per_hundred = zoo::rollmean(total_vaccinations_per_hundred, k = 14, fill = 0))  %>%
  dplyr::ungroup() 
  # dplyr::select(-c(new_cases_per_million, total_vaccinations_per_hundred))

head(owidData1)
owidData1 <- setDT(data)                           # Convert data.frame to data.table(owidData1)

#clean up data. We have both continent-level and country-level entries. 
removeList<-c("World","Africa","Asia", "Oceania", "Africa", "North America", "South America", 
              "Europe", "European Union", "International")  

continentList<-c("Africa","Asia","Oceania","Africa","North America","South America","Europe")

##drop/select for continents#######

continentData1 <- owidData1 %>% 
  group_by(continent, date) %>%
  mutate_at(c("new_cases_per_million","new_tests_per_thousand",
              "total_vaccinations_per_hundred", "stringency_index"), ~replace(., is.na(.), 0)) %>% 
  dplyr::summarise(
    cases = sum(new_cases_per_million),
    tests = sum(new_tests_per_thousand),
    vaccine_doses = sum(total_vaccinations_per_hundred),
    govt_index = mean(stringency_index)
  )
continentData1$date7<-as.Date(cut(continentData1$date,breaks = "1 week",
                                start.on.monday = TRUE))


countryData <- owidData1 %>% 
  subset(., ! location %in% removeList) %>% #ensure only country-level entries
  group_by(location, date7) %>%
  mutate_at(c("new_cases_per_million","new_tests_per_thousand",
              "total_vaccinations_per_hundred", "stringency_index"), ~replace(., is.na(.), 0)) %>% 
  dplyr::summarise(
    cases = mean(new_cases_per_million),
    tests = mean(new_tests_per_thousand),
    vaccine_doses = mean(total_vaccinations_per_hundred),
    govt_index = mean(stringency_index)
  ) 

table(countryData$location)
head(continentData1)

ggplot(continentData1)+
  geom_line(aes(x=date, y=cases))+
  facet_wrap(vars(continent))


#GISAID data##########
#Write a file containing the latest data using the following bash commands
#cat ../../../home5/nCov/Richard/gisaid/20210907/7_aligned.fasta | grep ">" > ./fileName
#sed 's/.//;s/|$//;s/|/,/g' 10052021 | cut -d, -f12 --complement | awk 'BEGIN{FS=OFS=","} \n
#{for (i=11;i<=NF;i++) sub(/\//,",",$i)} 1' > metadata_date.csv


#read data 
setwd(here::here())

mutationData <- vroom::vroom("../data/metadata_20220128.csv", col_names = F, 
                             col_select = c('X2', 'X10', 'X5', 'X8', 'X11'))

column_headings <- c("gisaid_ID", "country", "lineage", "date", "continent")
names(mutationData)<-column_headings

mutationData$date<-as.Date(mutationData$date)
mutationData$date7<-as.Date(cut(mutationData$date,breaks = "1 week",
                                 start.on.monday = TRUE))
mutationData$date14<-as.Date(cut(mutationData$date,breaks = "2 weeks",
                                 start.on.monday = TRUE))

mutationData <- drop_na(mutationData,lineage) #drop un-assigned lineages
dim(mutationData)

mutationData <- setDT(mutationData)         # Convert data.frame to data.table
#rename countries and continents in gisaid data to match owid data
renameFunction<-function(locationName, na.rm=T){
  gsub("_", " ", locationName) %>%    #replace underscore with space
    str_to_title(locationName) %>%    #change country to sentence case
    gsub("And ", "and ", .)  %>%      #clean up the tough names 
    gsub("Of ", "of ", .)  %>%
    gsub(" The", " the", .) %>%
    gsub("Usa", "United States", .)%>%
    gsub("Gambia", "The Gambia", .)
}


mutationData<- mutationData %>% 
  dplyr::mutate_at(
    c("continent", "country"), renameFunction
  ) %>%
  mutate(across(where(is.character),str_trim)) %>%
  drop_na()


#clean up UK naming 
mutationData$country[startsWith(mutationData$country, "Uk-")] <- "United Kingdom"

#clean up the (now confusing) lineage naming system
mutationData$lineage[startsWith(mutationData$lineage, "AY.")] <- "Delta"        #alias
mutationData$lineage[startsWith(mutationData$lineage, "B.1.617.2")] <- "Delta"

mutationData$lineage[startsWith(mutationData$lineage, "B.1.1.7")] <- "Alpha"
mutationData$lineage[startsWith(mutationData$lineage, "Q.")] <- "Alpha"         #alias

mutationData$lineage[startsWith(mutationData$lineage, "B.1.351.")] <- "Beta"

mutationData$lineage[startsWith(mutationData$lineage, "P.1")] <- "Gamma"

mutationData$lineage[startsWith(mutationData$lineage, "B.1.1.529")] <- "Omicron"
mutationData$lineage[startsWith(mutationData$lineage, "BA.")] <- "Omicron"      #alias


dict <- c("B.1.1.7" ="Alpha",
          "B.1.351" = "Beta",
          "P.1" = "Gamma",
          "B.1.617.2" = "Delta",
          "B.1.1.529" = "Omicron")  #the dict is now redundant as this is now taken care above, but hold on to it

mutationData$lineage2 <- as.character(dict[mutationData$lineage])

mutationData<-mutationData %>%
  mutate(lineage = coalesce(lineage2,lineage)) %>%
  dplyr::select(-c(lineage2))

table(mutationData$lineage)


#unique PANGO lineages in all data
unique(mutationData$lineage) %>% length()
