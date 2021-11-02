#Regression Analysis
# This script does data pre-processing and model selection of predictors of infection rates

#Read data

rm(list = ls())

## packages
require("pacman")
pacman::p_load(data.table, tidyverse, zoo, RColorBrewer, directlabels,corrplot, broom,
               stringr, spData, ggrepel, ggcorrplot, patchwork, cowplot, egg, here, visreg,
               glmmTMB, lme4, DHARMa, performance, MASS, mgcv) 


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
removeList<-c("World","Africa","Asia", "Oceania", "Africa", "North America", "South America", 
              "Europe", "European Union", "International")  

continentList<-c("Africa","Asia","Oceania","Africa","North America","South America","Europe")

##drop/select for continents#######
continentData0<-subset(owidData1, location %in% continentList) %>%   #continent level data
  dplyr::select(-c(stringency_index)) #drop the country-specific stringency field. Recalculate for continent-level below
continentData0$continent<-continentData0$location

###calculate average stringency index per continent and add to continent data#######
countryData<-owidData1 %>% 
  subset(., ! location %in% removeList) #ensure only country-level entries


contStringency<-countryData %>%
  dplyr::select(continent, date, stringency_index)%>%
  mutate_at("stringency_index", ~replace(., is.na(.), 0)) %>%  #replace NA with zero
  group_by(continent, date) %>% 
  summarise_each(list(stringency_index = mean))%>%         #continental average daily stringency index
  subset(continent %in% continentList)  # drop unknown/misspelt contient names

continentData1<-merge(continentData0, contStringency, by=c("continent", "date"), all.x = T)

head(continentData1)

#smooth govt stringency at the continent level too. Note, this isn't smoothed at country level
continentData2 <- continentData1 %>%
  dplyr::group_by(continent) %>% 
  dplyr::arrange(date) %>% 
  dplyr::mutate(smoothed_stringency_index = zoo::rollmean(stringency_index, k = 7, fill = 0)) %>%   #weekly rolling mean - for smoothing
  dplyr::ungroup() %>%
  dplyr::select(-c(stringency_index,location)) %>%
  distinct()

#visualize data to check if all is ok
ggplot(data=continentData2)+
  #geom_line(aes(x=date, y=new_tests_per_thousand, colour = "Tests per 1M")) +
  geom_line(aes(x=date, y=smoothed_stringency_index, colour = "Smoothed stringency Index"))+
  geom_line(aes(x=date, y=smoothed_cases_per_million, colour = "Smoothed Cases per 1M"))+
  #scale_x_date(date_labels = "%b-%y",date_breaks = "1 month", limits = as.Date(c('2019-12-01','2021-10-01')))+
  facet_wrap(vars(continent))+ ylab(NULL)+
  theme_bw()

#plot country data
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

#unique lineages in all data
unique(mutationData$lineage) %>% length()

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


# #drop countries with less than 100 sequences
# table1<-table(mutationData$country)
# mutationData<- mutationData %>% 
#   subset(.,country %in% names(table1[table1>100])) #drop

#clean up country names and the ever increasing in complexity lineage naming system
mutationData$country[startsWith(mutationData$country, "Uk-")] <- "United Kingdom"
mutationData$lineage[startsWith(mutationData$lineage, "AY")] <- "Delta"
mutationData$lineage[startsWith(mutationData$lineage, "Q")] <- "Alpha"

dict<-c("B.1.1.7" ="Alpha",
        "B.1.351" = "Beta",
        "P.1" = "Gamma",
        "B.1.617.2" = "Delta")

mutationData$lineage2 <- as.character(dict[mutationData$lineage])

mutationData<-mutationData %>%
  mutate(lineage = coalesce(lineage2,lineage)) %>%
  dplyr::select(-c(lineage2))

table(mutationData$lineage)

head(mutationData)

##1. Virus Diversity ie Distinct lineages per time period#################################
#diversity
#because some regions sequence more, distinct lineages may biased to more sequences. 
#We divide distinct lineages by total number of sequences to achieve a #diversity_score
# We will smooth the numbers afterwards

#total number of seqs per contiinent per day
seq_count<-dplyr::count(mutationData, continent, date)


pacman::p_load(plyr) 

cont_mutationData1<- mutationData%>%
  group_by(continent, date) %>% 
  ddply(~date+continent,summarise,distinct_lineages_per_continent_per_day=length(unique(lineage))) %>% 
  ungroup() %>%
  left_join (mutationData, 
             by= c("continent", "date"))  

lineage_diversity_df<-cont_mutationData1 %>%
  left_join (seq_count, 
             by= c("continent", "date")) %>%
  mutate(
    diversity_score=distinct_lineages_per_continent_per_day/n,  #diversity taking into account total number of sequences. 
    #this doesn't work too well at the beginning (pre-May 2020) when we have very few sequences. Diversity appears high while instaed the sequences were few
    #should we stick with distinct_lineages_per_continent_per_day?
    #New plan (19-10-2021), we drop data before 1st April 2020 before the linear regression
  )  %>%
  dplyr::select(continent, date, diversity_score, distinct_lineages_per_continent_per_day) %>%
  distinct() %>%    #drop duplicated rows
  group_by(continent) %>%
  dplyr::arrange(date) %>%
  mutate(
    smoothed_diveristy = zoo::rollmean(diversity_score, k = 14, fill = 0),  #biweekly rolling mean - for smoothing,
    smoothed_distinct_lineages = zoo::rollmean(distinct_lineages_per_continent_per_day, k = 14, fill = 0)
  ) 

pacman::p_unload(plyr) #will conflict with dplyr downstream



##2. VOC proportions per time period########

#define VOCs 
vocs<-c("Alpha", "Beta", "Gamma", "Delta")

### Summary of lineages per continent per day
dplyr::count(mutationData, continent, date, lineage) -> continent_daily_lineage_summary
# #count vocs per continent per day
vocs_count <- continent_daily_lineage_summary %>% 
  dplyr::filter(lineage %in% vocs) %>% 
  group_by(continent, date) %>% 
  summarise(
    voc_total = sum(n, na.rm = T)
  )
# #total sequences per continent per day
data_set_cont <- continent_daily_lineage_summary %>% 
  # left_join(vocs_count,
  #           by = c("continent", "date")) %>% 
  group_by(continent, date) %>% 
  summarise(
    tot_sequences = sum(n)
  )

# #proportion of vocs per continent per day
voc_proportions_continent <- data_set_cont %>% 
  left_join(vocs_count,
            by = c("continent", "date")) %>% 
  mutate(
    voc_total = replace(voc_total, which(is.na(voc_total)), 0),   
    voc_prop_cont = (voc_total/tot_sequences)
  ) %>%
  #dplyr::select(-c(tot_sequences, voc_total)) %>%  
  dplyr::group_by(continent) %>% 
  dplyr::arrange(date) %>%
  mutate(
    smoothed_voc_prop = zoo::rollmean(voc_prop_cont, k = 14, fill = 0),  #biweekly rolling mean - for smoothing
  ) 


####Comprehensive mutation data - for GLM purposes only (no OWID just yet) ####
head(voc_proportions_continent)
head(lineage_diversity_df)

comprehensive_mutation_df<-lineage_diversity_df %>% 
  left_join(voc_proportions_continent, by=c("continent", "date"))


##Create df #################################
## ie Merge gisaid with owid data
head(continentData2)
head(comprehensive_mutation_df)

df_continent<-continentData2 %>% 
  left_join(comprehensive_mutation_df, by=c("continent", "date"))

df_continent<-mutate_at(df_continent, c("smoothed_cases_per_million", "smoothed_vaccines_per_hundred", 
                                        "smoothed_stringency_index", "smoothed_diveristy", "smoothed_voc_prop",
                                        "smoothed_distinct_lineages"),
                        ~replace(., is.na(.), 0))
                        
                
head(df_continent)

regressionData<- df_continent %>%
  dplyr::group_by(continent) %>% 
  dplyr::arrange(date) %>%
  mutate(
    smoothed_diversity = zoo::rollmean(smoothed_diveristy, k = 7, fill = 0),  #we repeat because data is still noisy
  ) %>%
  dplyr::select(continent, date, smoothed_cases_per_million, smoothed_voc_prop, smoothed_diversity, smoothed_stringency_index,
                smoothed_vaccines_per_hundred, smoothed_distinct_lineages)%>%
  dplyr::filter(date >"2020-04-01" & date <"2021-09-01") #data seems a bit unstable before this date

  

ggplot(data= regressionData)+
  geom_line(aes(x=date, y=smoothed_voc_prop, color="voc"))+
  geom_line(aes(x=date, y=smoothed_diversity, color="diversity"))+
  facet_wrap(vars(continent))+
  theme_bw()


#all looks good! Time to fit a linear model 

# Wait, correlation test first 

regressionData %>% 
  nest(data = -continent) %>% 
  mutate(
    voc = map(data, ~ cor.test(.x$smoothed_cases_per_million, .x$smoothed_voc_prop)), # S3 list-col
    tidied = map(voc, tidy)
  ) %>% 
  unnest(tidied)  #repeat for other variables?. Actually just do this globally once and visualize



# corr data
corr_data<-regressionData %>% 
  ungroup () %>% 
  dplyr::select(smoothed_cases_per_million, smoothed_voc_prop, smoothed_diversity, 
                smoothed_stringency_index, smoothed_vaccines_per_hundred) %>% 
  dplyr::rename(VOCs = smoothed_voc_prop,
                `Infection rate` = smoothed_cases_per_million,
                Vaccination =smoothed_vaccines_per_hundred,
                `Virus diversity` =smoothed_diversity,
                `Stringency index` =smoothed_stringency_index)

# corr values
corr <- round(cor(corr_data),2)
# p values  
p.mat <- cor_pmat(corr_data)
# plot

corr_plot<-ggcorrplot(corr, hc.order = F, type = "lower",lab = T,
           outline.col = "white",
           #insig = "blank", 
           colors = c("royalblue4", "white", "#E46726"))+
  theme(panel.grid = element_blank())

pdf(paste0("./figures/correlation.pdf"), width=8, height = 8)
corr_plot
dev.off()



# Time to fit a linear model, for real this time :-)

# Multiple Linear Regression
    #lm([target variable] ~ [predictor variables], data = [data source])

# Use the plyr package to run regressions per continent instead of filtering each time
library(plyr)
#defn function for linear model
lm_with_plyr<-function(df){
  lm(smoothed_cases_per_million ~ smoothed_stringency_index +
       smoothed_vaccines_per_hundred +
       smoothed_voc_prop +
       smoothed_diversity,
     data = df)
}
#divide df into smaller dfs (creates a list of dfs) based on continent and run lm for each of the dfs
models<-dlply(regressionData, "continent", lm_with_plyr)
# Apply coef to each model and return a data frame
ldply(models, coef)
# Print the summary of each model
l_ply(models, summary, .print = TRUE)
# R Squared
laply(models, function(mod) summary(mod)$r.squared)

summary(models$Africa)

# Find fitted values

fitAndPlotFunction<-function(selected_continent){
  regressionData %>%
    dplyr::filter(continent==selected_continent) %>%
    mutate(fitted_cases = fitted(models[[selected_continent]]),
           smoothed_fitted = zoo::rollmean(fitted_cases, k = 14, fill = NA)) %>%
    #dplyr::arrange(date) %>%
    #mutate(smoothed_fitted = zoo::rollmean(fitted_cases, k = 14, fill = 0))
    ggplot()+
    geom_line(aes(x=date, y=smoothed_cases_per_million, color="Data"))+
    geom_line(aes(x=date, y=smoothed_fitted, color="Fitted"))+
    theme_bw()+
    scale_x_date(date_labels = "%b",date_breaks = "6 months", limits = as.Date(c('2020-01-01','2021-10-01')))+
    theme(panel.background = element_rect(fill = "transparent"), # bg of the panel
          plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
          axis.text.x = element_text(color="black", size=9.5),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.y = element_text(color="black", size=10),
          legend.position = "top",
          legend.justification = "left",
          legend.margin = margin(0,0,0,0),
          legend.box.margin = margin(-10,-10,-10,-10),
          legend.background = element_rect(fill = "transparent"), # get rid of legend bg
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())+ 
    guides(color=guide_legend(title=NULL))+
    scale_color_manual(values=c("black", "violetred1"))
  }


Africa<-fitAndPlotFunction("Africa")
Asia<-fitAndPlotFunction("Asia")
Europe<-fitAndPlotFunction("Europe")
North_America<-fitAndPlotFunction("North America")
Oceania<-fitAndPlotFunction("Oceania")
South_America<-fitAndPlotFunction("South America")















