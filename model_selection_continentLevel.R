#Regression Analysis
# This script does data pre-processing and model selection of predictors of infection rates

#Read data

rm(list = ls())

## packages
require("pacman")
pacman::p_load(data.table, tidyverse, zoo, RColorBrewer, directlabels,corrplot,
               stringr, spData, ggrepel, corrr, patchwork, cowplot, egg, here, visreg,
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
  geom_line(aes(x=date, y=new_tests_per_thousand, colour = "Tests per 1M")) +
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

mutationData$country[startsWith(mutationData$country, "Uk-")] <- "United Kingdom"

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
    #we stick with distinct_lineages_per_continent_per_day
  )  %>%
  dplyr::select(continent, date, diversity_score, distinct_lineages_per_continent_per_day) %>%
  distinct() %>%    #drop duplicated rows
  group_by(continent, date) %>%
  mutate(
    smoothed_diveristy = zoo::rollmean(diversity_score, k = 14, fill = 0),  #biweekly rolling mean - for smoothing,
    smoothed_distinct_lineages = zoo::rollmean(distinct_lineages_per_continent_per_day, k = 14, fill = 0)
  )
  
  

pacman::p_unload(plyr) #will conflict with dplyr downstream



##2. VOC proportions per time period########

#define VOCs 
vocs<-c("B.1.1.7", "P.1", "B.1.351", "B.1.617.2")

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
    smoothed_voc_prop = zoo::rollmean(voc_prop_cont, k = 14, fill = NA),  #biweekly rolling mean - for smoothing
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

# Model 1 
# Multiple Linear Regression
    #lm([target variable] ~ [predictor variables], data = [data source])

# df_continent$log_smoothed_cases_per_million<-log(df_continent$smoothed_cases_per_million)
# df_continent$log_smoothed_cases_per_million[which(is.infinite(df_continent$log_smoothed_cases_per_million))] <- 0. #replace inf with zero
# df_continent<-df_continent %>%
#   filter_all(all_vars(!is.infinite(.))). #drop infinite values

      #using log infection rates reduces variation but residuals are NOT normally distributed

###GLM additive
glm_additive1<- lm(smoothed_cases_per_million ~ smoothed_stringency_index +
                     smoothed_voc_prop +
                     smoothed_distinct_lineages *
                    continent, 
                  data = df_continent)
summary(glm_additive1)
#hist(residuals(glm_additive1))
# par(mfrow=c(2,2))
# visreg(glm_additive1)

glm_additive2<- lm(smoothed_cases_per_million ~ smoothed_stringency_index +
                     smoothed_vaccines_per_hundred +
                     smoothed_voc_prop +
                     smoothed_distinct_lineages *
                     continent, 
                   data = df_continent)
summary(glm_additive2)


##Predict ####
#get predicted values#####################
#use fitted, not predict
head(df_continent)
df_continent$allNoVaccine <- fitted(glm_additive1)

df_continent$all <- fitted(glm_additive2)

head(df_continent)

#smooth the predictions above
df_continent2 <- df_continent %>%
  group_by(continent) %>%
  dplyr::arrange(date) %>% 
  dplyr::mutate(all_smoothed = zoo::rollmean(all, k = 14, fill = NA)) %>%   #biweekly rolling mean - for smoothing
  dplyr::mutate(allNoVaccine_smoothed = zoo::rollmean(allNoVaccine, k = 14, fill = NA)) 


#Plot 1#################################
#select per continent to plot

#selected<-"South America"
comprehensive_df3<-df_continent2 %>%
  filter(continent==selected)

#pdf(paste0("../figures/Fitted", selected, ".pdf"), width=8, height = 4.5)
plot1<-ggplot(data=df_continent2)+
  geom_line(aes(x=date, y=smoothed_cases_per_million, color="Data"))+
  geom_line(aes(x=date, y=allNoVaccine_smoothed, color="Fitted 1"))+
  geom_line(aes(x=date, y=all_smoothed, color="Fitted 2"))+
  theme_bw()+
  theme(legend.position = "right") +
  scale_x_date(date_labels = "%b\n%Y",date_breaks = "3 months", limits = as.Date(c('2019-12-01','2021-08-01')))+
  #theme(axis.title.x = element_text(color="black", size=15, face="bold"))+
  theme(axis.text.x = element_text(color="black", size=9.5),
        axis.title.x = element_blank(),
        axis.title.y = element_text(color="black", size=11),
        axis.text.y = element_text(color="black", size=10),
        legend.position = "top",
        #panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+ 
  ylab("Cases per million")+
  scale_fill_brewer(palette = "Set3")+
  #ggtitle(selected)+
  guides(color=guide_legend(title=NULL))+
  facet_wrap(vars(continent))
#dev.off()

plot1










##########################################################

        #Data vs model look better when using log infections, but normality of residuals is now outside the window.

# Model 2
###GLM multiplicative
#Similar estimates when continent is additive  

# Model 3
### GMEM

gmem1 <- lmer(smoothed_cases_per_million ~ smoothed_stringency_index + 
                (smoothed_diveristy) + 
                (smoothed_voc_prop)+
                continent+   
                (1| smoothed_stringency_index), 
              data = df_continent)

summary(gmem1)
hist(residuals(gmem1))
visreg(gmem1)



gmem2 <- lmer(smoothed_cases_per_million ~ smoothed_stringency_index + 
                smoothed_diveristy + 
                smoothed_voc_prop+
                continent+
                (1 | smoothed_voc_prop),#+
                #(1 | continent), 
              data = df_continent)

summary(gmem2)
hist(residuals(gmem2))
visreg(gmem2)


gmem3 <- lmer(log_smoothed_cases_per_million ~ smoothed_stringency_index + 
                smoothed_diveristy + 
                smoothed_voc_prop+
                continent+
                (1 |continent)+
                (1 |smoothed_diveristy), 
              #(1 | continent), 
              data = df_continent)

summary(gmem3)
hist(residuals(gmem3)) #residuals slightly skewed to left
visreg(gmem3)
pp_check(gmem3)

simulationOutput<-simulateResiduals(gmem3, plot = T)
testDispersion(simulationOutput)  #there is a clear overdispersion

###


compare_performance(glm_additive, gmem1, gmem2, gmem3)














