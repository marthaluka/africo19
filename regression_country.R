#Regression Analysis at country-level
# This script does data pre-processing and model selection of predictors of infection rates

#Read data

#rm(list = ls())

## packages
require("pacman")
pacman::p_load(data.table, tidyverse, zoo, RColorBrewer, directlabels,corrplot, broom,
               stringr, spData, ggrepel, corrr, patchwork, cowplot, egg, here,
               imputeTS) 


#OWID data######### 
#read data
owidData0<- read.csv(
  "https://raw.githubusercontent.com/owid/covid-19-data/master/public/data/owid-covid-data.csv",
  fileEncoding="UTF-8-BOM", stringsAsFactors = FALSE) %>%                                   #select relevant fields
  dplyr::select(continent, location, date, new_cases_per_million,total_vaccinations_per_hundred, 
                stringency_index, new_tests_per_thousand) 

head(owidData0)

length(owidData0$new_tests_per_thousand)
sum(is.na(owidData0$new_tests_per_thousand))  #drop number of tests per capita  due to missing data (>50% entries are NA)

#date formats
owidData0$date<-as.Date(owidData0$date)

#Pre-processing
#biweekly rolling mean - for smoothing

owidData1 <- owidData0 %>%
  mutate(across(where(is.character),str_trim)) %>% 
  dplyr::filter(new_cases_per_million>=0) %>%    #drop any negative values 
  dplyr::group_by(location) %>% 
  dplyr::arrange(date) %>% 
  dplyr::mutate(smoothed_cases_per_million = zoo::rollmean(new_cases_per_million, k = 14, fill = 0)) %>%  
  dplyr::mutate(smoothed_vaccines_per_hundred = zoo::rollmean(total_vaccinations_per_hundred, k = 14, fill = 0)) %>% 
  dplyr::ungroup()  %>%
  dplyr::select(-c(new_cases_per_million, total_vaccinations_per_hundred))

head(owidData1)

#clean up data. We have both continent-level and country-level entries. 
removeList<-c("World","Africa","Asia", "Oceania", "North America", "South America", 
              "Europe", "European Union", "International")  

 
###country level data#######
countryData<-owidData1 %>% 
  subset(., ! location %in% removeList) #ensure only country-level entries


#visualize 
ggplot(data = countryData)+
  geom_line(aes(x=date, y=smoothed_cases_per_million, colour = "Infection rate"))+
  geom_line(aes(x=date, y=stringency_index, colour="stringency"))+
  facet_wrap(vars(location))


#GISAID data##########
#Write a file containing the latest data using the following bash commands
#cat ../../../home5/nCov/Richard/gisaid/20210907/7_aligned.fasta | grep ">" > ./fileName
#sed 's/.//;s/|$//;s/|/,/g' 10052021 | cut -d, -f12 --complement | awk 'BEGIN{FS=OFS=","} \n
#{for (i=11;i<=NF;i++) sub(/\//,",",$i)} 1' > metadata_date.csv


#read data 
here::here()
mutationData<-read.csv("../data/metadata_20210907.csv", header = F) %>%
  dplyr::select(V2, V10, V5, V8, V11) 

column_headings <- c("gisaid_ID", "country", "lineage", "date", "continent")
names(mutationData)<-column_headings

mutationData$date<-as.Date(mutationData$date)

head(mutationData)

#rename countries and continents in gisaid data to match owid data
renameFunction<-function(locationName, na.rm=T){
  gsub("_", " ", locationName) %>%    #replace underscore with space
    str_to_title(locationName) %>%    #change country to sentence case
    gsub("And ", "and ", .)  %>%      #clean up the tough names 
    gsub("Of ", "of ", .)  %>%
    gsub(" The", " the", .) %>%
    gsub("Usa", "United States", .)#%>%
    #gsub("Gambia", "The Gambia", .) #dont need this
  
  # "OWID, GISAID"
  # "Czechia", "Czech Republic" Czech Republic
  # "Guyana", "French Guiana" French Guiana
  # "Democratic Republic of Congo", "Democratic Republic of the Congo"
  # "Gambia", "The Gambia" Gambia
  # "Sint Maarten (Dutch part)","Sint Maarten" Saint Martin
  
}


mutationData<- mutationData %>% 
  dplyr::mutate_at(
    c("continent", "country"), renameFunction
  ) %>%
  mutate(across(where(is.character),str_trim)) %>%
  drop_na()


# #drop countries with less than 500 sequences (makes it easier for more stringent filtering downstream)
table1<-table(mutationData$country)
mutationData<- mutationData %>% 
subset(.,country %in% names(table1[table1>500])) #drop

#clean up country names and the ever increasing in complexity lineage naming system
mutationData$country[startsWith(mutationData$country, "Uk-")] <- "United Kingdom"
mutationData$lineage[startsWith(mutationData$lineage, "AY.")] <- "Delta"
mutationData$lineage[startsWith(mutationData$lineage, "Q.")] <- "Alpha"

dict<-c("B.1.1.7" ="Alpha",
        "B.1.351" = "Beta",
        "P.1" = "Gamma",
        "B.1.617.2" = "Delta")

mutationData$lineage2 <- as.character(dict[mutationData$lineage])

mutationData<-mutationData %>%
  mutate(lineage = coalesce(lineage2,lineage)) %>%
  dplyr::select(-c(lineage2))

table(mutationData$country)
table(owidData0$location)
# "OWID, GISAID"
# "Czechia", "Czech Republic"
# "Guyana", "French Guiana"
# "Democratic Republic of Congo", "Democratic Republic of the Congo"
# "Gambia", "The Gambia"
# "Sint Maarten (Dutch part)","Sint Maarten"

head(mutationData)

##1. Virus Diversity ie Distinct lineages per time period#################################
#diversity
#because some regions sequence more, distinct lineages may biased to more sequences. 
#We divide distinct lineages by total number of sequences to achieve a #diversity_score
# We will smooth the numbers afterwards

#total number of seqs per contiinent per day
seq_count<-dplyr::count(mutationData, country, date)


pacman::p_load(plyr) 

country_mutationData1<- mutationData%>%
  group_by(country, date) %>% 
  ddply(~date+country,summarise,distinct_lineages_per_country_per_day=length(unique(lineage))) %>% 
  ungroup() %>%
  left_join (mutationData, 
             by= c("country", "date"))  

lineage_diversity_df<-country_mutationData1 %>%
  left_join (seq_count, 
             by= c("country", "date")) %>%
  mutate(
    diversity_score=distinct_lineages_per_country_per_day/n,  #diversity taking into account total number of sequences. 
    #this doesn't work too well at the beginning (pre-May 2020) when we have very few sequences. Diversity appears high while instaed the sequences were few
    # drop data before 1st April 2020 before the linear regression
  )  %>%
  dplyr::select(country, date, diversity_score, distinct_lineages_per_country_per_day) %>%
  distinct() %>%    #drop duplicated rows
  group_by(country) %>%
  dplyr::arrange(date) %>%
  mutate(
    imputed_diversity_score = na_ma(diversity_score, k = 40, weighting="simple", maxgap=Inf),#impute missing data using rolling average
    smoothed_diveristy = zoo::rollmean(imputed_diversity_score, k = 14, na.pad = T, align = "right"),  #biweekly rolling mean - for smoothing,
    smoothed_distinct_lineages = zoo::rollmean(distinct_lineages_per_country_per_day, k = 14, na.pad = T, align = "right")
  ) 

pacman::p_unload(plyr) #will conflict with dplyr downstream



##2. VOC proportions per time period########

#define VOCs 
vocs<-c("Alpha", "Beta", "Gamma", "Delta")

### Summary of lineages per country per day
dplyr::count(mutationData, country, date, lineage) -> country_daily_lineage_summary
# #count vocs per country per day
vocs_count <- country_daily_lineage_summary %>% 
  dplyr::filter(lineage %in% vocs) %>% 
  group_by(country, date) %>% 
  summarise(
    voc_total = sum(n, na.rm = T)
  )
# #total sequences per country per day
data_set_country <- country_daily_lineage_summary %>% 
  # left_join(vocs_count,
  #           by = c("country", "date")) %>% 
  group_by(country, date) %>% 
  summarise(
    tot_sequences = sum(n)
  )

# #proportion of vocs per country per day
voc_proportions_country <- data_set_country %>% 
  left_join(vocs_count,
            by = c("country", "date")) %>% 
  mutate(
    voc_total = replace(voc_total, which(is.na(voc_total)), NA),   
    voc_prop = (voc_total/tot_sequences)
  ) %>%
  #dplyr::select(-c(tot_sequences, voc_total)) %>%  
  dplyr::group_by(country) %>% 
  dplyr::arrange(date) %>%
  mutate(
    imputed_voc_prop = na_ma(voc_prop, k = 14, weighting="simple", maxgap=50), #impute missing values with rolling average
    smoothed_voc_prop = zoo::rollmean(imputed_voc_prop, k = 14, na.pad = T, align = "right")  #biweekly rolling mean - for smoothing
  ) 


####Comprehensive mutation data - for GLM purposes only (no OWID just yet) ####
head(voc_proportions_country)
head(lineage_diversity_df)

comprehensive_mutation_df<-lineage_diversity_df %>% 
  left_join(voc_proportions_country, by=c("country", "date"))


##Create df #################################
## ie Merge gisaid with owid data
head(countryData)
head(comprehensive_mutation_df)

df_country<-countryData %>% 
  dplyr::rename(., country = location) %>% 
  left_join(comprehensive_mutation_df, by=c("country", "date"))

#filter out countries before regression
  #criteria is data completeness

df_data_completeness <- df_country%>%
  dplyr::filter(
    ! country == "Palau"
    ) %>%
  dplyr::group_by(continent,country) %>%
  dplyr::summarise(missingSequenceData = sum(is.na(tot_sequences)))#use this

df_data_completeness %>%
  dplyr::filter(
    missingSequenceData<100,
    )%>%
  ggplot()+
  geom_col(aes(y=country, x=missingSequenceData))+
  facet_wrap(vars(continent))+
  theme_bw()


table2<- df_data_completeness%>%
  dplyr::filter(
    missingSequenceData<150
  )

df_country_filtered <- df_country %>% 
  subset(.,country %in% table2$country) #drop



df_country_filtered<-mutate_at(df_country_filtered, c("smoothed_cases_per_million", "smoothed_vaccines_per_hundred",
                                        "stringency_index", "smoothed_diveristy", "smoothed_voc_prop",
                                        "smoothed_distinct_lineages"),
                        ~replace(., is.na(.), 0))


head(df_country_filtered)

#visualize 
ggplot(data = df_country_filtered)+
  geom_line(aes(x=date, y=smoothed_cases_per_million, colour = "Infection rate"))+
  geom_line(aes(x=date, y=stringency_index, colour="stringency"))+
  facet_wrap(vars(country))

df_country_filtered$smoothed_voc_prop

ggplot(data= df_country_filtered)+
  geom_line(aes(x=date, y=smoothed_voc_prop, color="voc"))+
  geom_line(aes(x=date, y=smoothed_diveristy))+
  facet_wrap(vars(country))+
  theme_bw()
#drop palau


df_country_filtered<-df_country_filtered[df_country_filtered$smoothed_diveristy !=0,]

regressionData<- df_country_filtered %>%
  dplyr::select(country, date, smoothed_cases_per_million, smoothed_voc_prop, smoothed_diveristy, stringency_index,
                smoothed_vaccines_per_hundred)%>%
  dplyr::filter(date >"2020-04-01" & date <"2021-09-01") #data seems a bit unstable before this date




#all looks good! Time to fit a linear model :-)


# Multiple Linear Regression
#lm([target variable] ~ [predictor variables], data = [data source])

# Use the plyr package to run regressions per continent instead of filtering each time
library(plyr)
#defn function for linear model
lm_with_plyr<-function(df){
  lm(smoothed_cases_per_million ~ stringency_index +
       smoothed_vaccines_per_hundred +
       smoothed_voc_prop +
       smoothed_diveristy,
     data = df)
}
#divide df into smaller dfs (creates a list of dfs) based on country and run lm for each of the dfs
models<-dlply(regressionData, "country", lm_with_plyr)
# Apply coef to each model and return a data frame
ldply(models, coef)
# Print the summary of each model
l_ply(models, summary, .print = TRUE)
# R Squared
laply(models, function(mod) summary(mod)$r.squared)

pacman::p_unload(plyr)
# Find fitted values
summary(models$France)



fitAndPlotFunction<-function(selected_country){
  regressionData %>%
    dplyr::filter(country==selected_country) %>%
    #get fitted infections per capita and smooth using rolling mean
    dplyr::mutate(fitted_cases = fitted(models[[selected_country]]),
           smoothed_fitted = zoo::rollmean(fitted_cases, k = 14, fill = NA)) %>%
    #dplyr::arrange(date) %>%
    #mutate(smoothed_fitted = zoo::rollmean(fitted_cases, k = 14, fill = 0))
    ggplot()+
    geom_line(aes(x=date, y=smoothed_cases_per_million, color="Data"))+
    geom_line(aes(x=date, y=smoothed_fitted, color="Fitted"))+
    theme_bw()+
    ylab('Cases') +
    scale_x_date(date_labels = "%b \n %Y",date_breaks = "6 months", limits = as.Date(c('2020-01-01','2021-10-01')))+
    theme(#panel.background = element_rect(fill = "transparent"), # bg of the panel
          #plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
          axis.text.x = element_text(color="black", size=4),
          axis.title.x = element_blank(),
          axis.title.y = element_text(color="black", size=4),
          axis.text.y = element_text(color="black", size=4),
          plot.title = element_text(size=5),
          legend.position = "top",
          legend.text = element_text(color="black", size = 5),
          legend.justification = "left",
          legend.margin = margin(0,0,0,0),
          legend.box.margin = margin(-10,-10,-10,-10),
          legend.key.width =  unit(2, "mm"),
          legend.background = element_rect(fill = "transparent"), # get rid of legend bg
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())+ 
    guides(color=guide_legend(title=NULL))+
    ggtitle(selected_country)+
    scale_color_manual(values=c("black", "violetred1"))
}















