#Regression Analysis
# This script does data pre-processing and model selection of predictors of infection rates

#Read data

rm(list = ls())

## packages
require("pacman")
pacman::p_load(data.table, tidyverse, zoo, RColorBrewer, directlabels, 
               stringr, spData, ggrepel, corrr, patchwork, cowplot, egg, here, visreg,
               lme4, DHARMa, performance, MASS, mgcv) 


#OWID data######### 
#read data
owidData0<- read.csv(
  "https://raw.githubusercontent.com/owid/covid-19-data/master/public/data/owid-covid-data.csv",
  fileEncoding="UTF-8-BOM", stringsAsFactors = FALSE) %>%                                   #select relevant fields
  dplyr::select(continent, location, date, new_cases_per_million,total_vaccinations_per_hundred, 
                stringency_index) 

head(owidData0)

#date formats
owidData0$date<-as.Date(owidData0$date)
owidData0$date7<-as.Date(cut(owidData0$date,breaks = "1 weeks",start.on.monday = FALSE))
owidData0$date14<-as.Date(cut(owidData0$date,breaks = "2 weeks",start.on.monday = FALSE))

#Pre-processing
#biweekly rolling mean - for smoothing

owidData1 <- owidData0 %>%
  mutate(across(where(is.character),str_trim)) %>% 
  filter(new_cases_per_million>0) %>%    #drop any negative/zero values ie treat zero as missing data
  dplyr::group_by(location) %>% 
  dplyr::arrange(date) %>% 
  dplyr::mutate(smoothed_cases_per_million = zoo::rollmean(new_cases_per_million, k = 14, fill = NA)) %>%  
  dplyr::mutate(smoothed_vaccines_per_hundred = zoo::rollmean(total_vaccinations_per_hundred, k = 14, fill = NA)) %>%   
  dplyr::ungroup()  %>%
  dplyr::select(-c(new_cases_per_million, total_vaccinations_per_hundred))

head(owidData1)

#clean up data. We have both continent-level and country-level entries. 
removeList<-c("World","Africa","Asia", "Oceania", "Africa", "North America", "South America", 
              "Europe", "European Union", "International")  

continentList<-c("Africa","Asia","Oceania","Africa","North America","South America","Europe")

##drop/select for continents#######
freqTable<-table(owidData1$location) #frequency of country-specific entries

owidData2<-owidData1 %>% 
  subset(., ! location %in% removeList) %>% #ensure only country-level entries
  subset(.,location %in% names(freqTable[freqTable>200])) #drop

head(owidData2)

countryLevelData_owid<-owidData2

continentData0<-subset(owidData1, location %in% continentList) %>%   #continent level data
  dplyr::select(-c(stringency_index)) #drop the country-specific stringency field. Recalculate for continent-level below
continentData0$continent<-continentData0$location

###calculate average stringency index per continent and add to continent data#######
contStringency<-owidData2 %>%
  dplyr::select(continent, date, stringency_index)%>%
  mutate_at("stringency_index", ~replace(., is.na(.), 0)) %>%  #replace NA with zero
  group_by(continent, date) %>% 
  summarise_each(funs(mean))%>%  #continental average of stringency index
  subset(continent %in% continentList)  # drop unknown/misspelt contient names

continentData1<-merge(continentData0, contStringency, by=c("continent", "date"), all.x = T)

head(continentData1)

#smooth govt stringency at the continent level too. Note, this isn't smoothed at country level
continentData2 <- continentData1 %>%
  dplyr::group_by(continent) %>% 
  dplyr::arrange(date) %>% 
  dplyr::mutate(smoothed_stringency_index = zoo::rollmean(stringency_index, k = 7, fill = NA)) %>%   #weekly rolling mean - for smoothing
  dplyr::ungroup() %>%
  dplyr::select(-c(stringency_index))

#visualize data to check if all is ok
ggplot(data=continentData2)+
  #geom_line(aes(x=date, y=new_cases_per_million, colour = "Cases per 1M"))+
  geom_line(aes(x=date, y=smoothed_stringency_index, colour = "Smoothed stringency Index"))+
  geom_line(aes(x=date, y=smoothed_cases_per_million, colour = "Smoothed Cases per 1M"))+
  #scale_x_date(date_labels = "%b-%y",date_breaks = "1 month", limits = as.Date(c('2019-12-01','2021-10-01')))+
  facet_wrap(vars(continent))+ ylab(NULL)+
  theme_bw()

#plot country data
ggplot(data = countryLevelData_owid)+
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
mutationData$date7<-as.Date(cut(mutationData$date,breaks = "1 weeks",start.on.monday = FALSE))
mutationData$date14<-as.Date(cut(mutationData$date,breaks = "2 weeks",start.on.monday = FALSE))

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
  mutate(across(where(is.character),str_trim)) 

        # #drop countries with less than 100 sequences
        # table1<-table(mutationData$country)
        # mutationData<- mutationData %>% 
        #   subset(.,country %in% names(table1[table1>100])) #drop

mutationData$country[startsWith(mutationData$country, "Uk-")] <- "United Kingdom"

table(mutationData$country)
#table(owidData0$location)

#From here, we branch into two analyses: country-specific and continent-specific
  #(where we lump every country in the continent together)

##1. Distinct lineages per time period#################################
  #At the country level, we will use an interval of 2 weeks since some countries dont sequence so much
  #At continent level, we will smooth the numbers afterwards

pacman::p_load(plyr) 

#country
tmp<- mutationData %>% group_by(country, date14) %>% 
  ddply(~date14+country,summarise,distinct_lineages_per_country_biweekly=length(unique(lineage))) %>% 
  ungroup()  %>% 
  left_join (mutationData, 
             by= c("country", "date14"))
   
#continent

cont_mutationData1<- mutationData%>%
  group_by(continent, date) %>% 
  ddply(~date+continent,summarise,distinct_lineages_per_continent_per_day=length(unique(lineage))) %>% 
  ungroup() %>%
  left_join (mutationData, 
             by= c("continent", "date")) 

pacman::p_unload(plyr) #will conflict with dplyr downstream

##2. VOC proportions per time period########

#define VOCs 
vocs<-c("B.1.1.7", "P.1", "B.1.351", "B.1.617.2")

### Summary of lineages per location per day

#country

dplyr::count(mutationData1, country, date14, lineage) -> country_daily_lineage_summary

# count vocs per country every 2 weeks
vocs_count1 <- country_daily_lineage_summary %>% 
  filter(lineage %in% vocs) %>% 
  group_by(country, date14) %>% 
  summarise(
    voc_total = sum(n, na.rm = T)
  )
# #total sequences per country every 2 weeks
data_set3 <- country_daily_lineage_summary %>% 
  left_join(vocs_count1,
            by = c("country", "date14")) %>% 
  group_by(country, date14) %>% 
  summarise(
    tot_sequences = sum(n)
  )
# #percent of vocs per country every 2 weeks
voc_proportions_country <- data_set3 %>% 
  left_join(vocs_count1,
            by = c("country", "date14")) %>% 
  mutate(
    voc_total = replace(voc_total, which(is.na(voc_total)), 0),
    voc_prop_country = (voc_total/tot_sequences)*100
  )%>%
  dplyr::select(-c(tot_sequences, voc_total))


#continent
dplyr::count(cont_mutationData1, continent, date, lineage) -> continent_daily_lineage_summary
# #count vocs per continent per day
vocs_count2 <- continent_daily_lineage_summary %>% 
  filter(lineage %in% vocs) %>% 
  group_by(continent, date) %>% 
  summarise(
    voc_total = sum(n, na.rm = T)
  )
# #total sequences per continent per day
data_set_cont <- continent_daily_lineage_summary %>% 
  left_join(vocs_count2,
            by = c("continent", "date")) %>% 
  group_by(continent, date) %>% 
  summarise(
    tot_sequences = sum(n)
  )
# #proportion of vocs per continent per day
voc_proportions_continent <- data_set_cont %>% 
  left_join(vocs_count2,
            by = c("continent", "date")) %>% 
  mutate(
    voc_total = replace(voc_total, which(is.na(voc_total)), 0),
    voc_prop_cont = (voc_total/tot_sequences)*100
  ) %>%
  dplyr::select(-c(tot_sequences, voc_total))


####Comprehensive mutation data (no OWID just yet) ####
head(mutationData1)
head(voc_proportions_country)
head(voc_proportions_continent)

comprehensive_mutation_df<-mutationData1 %>% 
  #left_join(voc_proportions_country, by=c("country", "date14")) %>% 
  left_join(voc_proportions_continent, by=c("continent", "date"))


##Create df #################################
## ie Merge gisaid with owid data
head(countryLevelData_owid)
head(continentData2)
head(comprehensive_mutation_df)


#country level df
country_level<- comprehensive_mutation_df  %>%  
  dplyr::select(-c(date14, distinct_lineages_per_continent_per_day, gisaid_ID, 
                   lineage, date7, voc_prop_cont))%>% 
  inner_join(
    (rename(countryLevelData_owid, country = location)), by=c("country", "date")
  ) %>% unique()

#continent df
continent_level<-comprehensive_mutation_df  %>% 
  dplyr::select(-c(date14, distinct_lineages_per_country_biweekly, country, gisaid_ID, 
                   lineage, date7, voc_prop_country))  %>% 
  inner_join(
    continentData2, by=c("continent", "date")
  ) %>% unique() %>%
           #smooth voc proportions and distinct lineages  #for regression
  dplyr::mutate(
    smoothed_distinct_lineages_per_day = zoo::rollmean(distinct_lineages_per_continent_per_day, k = 7, fill = NA),
    smoothed_voc_prop = zoo::rollmean(voc_prop_cont, k = 7, fill = NA)
  ) %>%
  dplyr::select(-c(voc_prop_cont, distinct_lineages_per_continent_per_day)
  )
  
head(country_level)
head(continent_level)

rm(list=c("comprehensive_mutation_df", "continent_daily_lineage_summary", continent))

                      #And now we have data to fit a linear, or mixed effects model 


#Model selection###########
##GLMs###########
    #all combined minus vaccination
    #lm([target variable] ~ [predictor variables], data = [data source])

# selected<-"South America"
# comprehensive_df2<-comprehensive_df %>%
#   filter(continent==selected)

###GLM additive
glm_additive<- lm(new_smoothed_cases_per_million ~ smoothed_voc_prop +
               smoothed_distinct_lineages_per_day +
               smoothed_stringency_index +
               continent, 
             data = comprehensive_df2)
summary(glm_additive)
hist(residuals(glm_additive))

par(mfrow=c(2,2))
visreg(glm_additive) #glms dont suit this data very well
simulateResiduals(glm_additive, plot = TRUE)
performance::pp_check(glm_additive)

###GLM interactive
glm_interactive <- lm(new_smoothed_cases_per_million ~ smoothed_voc_prop +
               smoothed_distinct_lineages_per_day +
               smoothed_stringency_index *
               continent, 
             data = comprehensive_df2)
summary(glm_interactive)
hist(residuals(glm_interactive))
visreg(glm_interactive) #glms dont suit this data very well
simulateResiduals(glm_interactive, plot = TRUE)

performance::pp_check(glm_interactive)

##GMEM###########
#unconditional

country_level[is.na(country_level)]<-0

gmem1 <- lmer(smoothed_cases_per_million ~ voc_prop_country + 
                distinct_lineages_per_country_biweekly + 
                stringency_index+
                (1 + stringency_index | country), 
              data = country_level_model_df)



###With continent as a random effect
gmem1 <- lmer(new_smoothed_cases_per_million ~ smoothed_voc_prop + 
                smoothed_stringency_index + 
                smoothed_distinct_lineages_per_day+
                (1 + smoothed_distinct_lineages_per_day | continent), 
              data = comprehensive_df2)


summary(gmem1)
visreg(gmem1, xvar ="new_smoothed_cases_per_million", by="continent")
simulateResiduals(gmem1, plot = TRUE). #still not the best

#another glm
glm.poi <- glm(new_smoothed_cases_per_million ~ smoothed_voc_prop+
                   smoothed_stringency_index + 
                   smoothed_distinct_lineages_per_day*
                   continent,
                 data = comprehensive_df2,
                 family = poisson)
summary(glm.poi) #highly overdispersed. Let's try quasipoisson ie nb

#quasipoisson ie nb 

glm.quasipoi <- glm.nb(new_smoothed_cases_per_million ~ smoothed_voc_prop+
                         smoothed_stringency_index + 
                         smoothed_distinct_lineages_per_day+
                         continent, 
                       data = comprehensive_df2)

summary(glm.quasipoi)
visreg(glm.quasipoi). #pthoo non-linear

###With continent + vaccination + diversity as a random effects


#Three level GMEM
head(country_level)
head(continent_level)

#filter out countries with inadequate data entries
table1<-table(country_level$country)
country_level_model_df<-subset(country_level, country %in% names(table1[table1>200]))
table(country_level_model_df$country)

lmer(smoothed_cases_per_million ~ voc_prop_country+
       (distinct_lineages_per_country_biweekly|continent.x) +
      (stringency_index|country),
     data=country_level_model_df,
     REML = FALSE, 
     control = lmerControl(optimizer ="Nelder_Mead"))




stringr::str_squish("   Enda n    kalale ")



##GAMs###########
  #change date to numeric
library(lubridate)
comprehensive_df2$numeric.date<-decimal_date(comprehensive_df2$date)

gam_m1 <- gam(new_smoothed_cases_per_million ~ continent + #continent still not random
                s(numeric.date, 
                                                 smoothed_stringency_index,
                                                 smoothed_voc_prop,
                                                 smoothed_distinct_lineages_per_day,
                                                 k = 50), 
              data = comprehensive_df2, method = "REML")

simulateResiduals(gam_m1, plot = TRUE)
summary(gam_m1)
visreg(gam_m1)

###GAM model 2
comprehensive_df2$continent<-as.factor(comprehensive_df2$continent)
    #Need to specify this because R reads categorical variables in as characters, not factors

gam_m2 <- gam(new_smoothed_cases_per_million ~ s(numeric.date, 
                  smoothed_stringency_index,
                  smoothed_voc_prop,
                  smoothed_distinct_lineages_per_day,
                  k = 50,
                  by = continent), 
              data = comprehensive_df2, method = "REML")

summary(gam_m2)
visreg(gam_m2)
visreg(gam_m2, xvar = "new_smoothed_cases_per_million", by = "continent")
simulateResiduals(gam_m2, plot = TRUE)
performance::pp_check(gam_m2)

#compare

compare_performance(gmem1, glm.quasipoi, gam_m1, gam_m2)


#all combined including vaccination
lmInfectionsXX <- lm(new_smoothed_cases_per_million ~ voc_prop+distinct_lineages_per_day+
                       stringency_index+total_vaccinations_per_hundred, 
                     data = comprehensive_df2)
sm_XX<-summary(lmInfectionsXX)
sm_XX
mean(sm_XX$residuals^2)  #root mean square error

#get predicted values#####################
#use fitted, not predict
head(comprehensive_df2)
comprehensive_df2$allNoVaccine <- predict(object = lmInfectionsX, newdata = comprehensive_df2)
comprehensive_df2$all <- predict(object = lmInfectionsXX, newdata = comprehensive_df2)

head(comprehensive_df2)

#smooth the predictions above
comprehensive_df3 <- comprehensive_df2 %>%
  dplyr::arrange(date) %>% 
  dplyr::mutate(all_smoothed = zoo::rollmean(all, k = 14, fill = NA)) %>%   #biweekly rolling mean - for smoothing
  dplyr::mutate(allNoVaccine_smoothed = zoo::rollmean(allNoVaccine, k = 14, fill = NA)) 



#Plot 1#################################
#select per continent to plot

selected<-"South America"
comprehensive_df3<-comprehensive_df %>%
  filter(continent==selected)

#pdf(paste0("../figures/Fitted", selected, ".pdf"), width=8, height = 4.5)
plot1<-ggplot(data=comprehensive_df3)+
  geom_line(aes(x=date, y=new_smoothed_cases_per_million, color="Data"))+
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
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+ 
  ylab("Cases per million")+
  scale_fill_brewer(palette = "Set3")+
  #ggtitle(selected)+
  guides(color=guide_legend(title=NULL))
#dev.off()

plot1


###Plot topLineage dynamics and infection rates################################################
head(continentData1)
head(mutationData1)
mutationData1$date2<-as.Date(cut(mutationData1$date,breaks = "2 weeks",start.on.monday = FALSE))

head(comprehensive_mutation_df)


lineages<-c("B.1.351","B.1.617.2","B.1","AY.4","B.1.1.7","B.1.351.2","B.1.1",
            "B.1.1.214","B.1.1.284","R.1","B.1.177","B.1.160","B.1.2","B.1.429",
            "B.1.526","P.1","D.2","A.2.2","C.37","P.2","B.1.1.28","B.1.1.33")
  
typeof(lineages) 
 
#Top X lineages per continent
summ_lineages<-comprehensive_mutation_df %>%
  na.omit() %>%
  count(continent, lineage) %>%
  group_by(continent)  %>%
  arrange(desc(n)) %>%
  slice(1:7) 

#find proprtions using this technique and facet wrap

ourTopLineages<-as.list(unique(summ_lineages$lineage))

#create a custom pallete to fix lineage colours across continents  
colour_palette<-c("azure1", "darksalmon", "dodgerblue", "gold1", "hotpink", "ivory4", 
                  "lavenderblush", "lightpink4", "orange1", "palegoldenrod", "plum4",
                  "purple3", "turquoise")





selection = selected            #select for continent #defined above
selected_mutationData<-mutationData1[mutationData1$continent==selection,]
selected_continentData1<-continentData1[continentData1$continent==selection,]

#######Count of top lineages
lineageCount<-as.data.frame(table(selected_mutationData$lineage))
lineageCount<-lineageCount[order(lineageCount$Freq),]
lineageCountTop5<-tail(lineageCount,n=5)
lineageCountTop5

#create a tmp df to add VOCs in case they aren't among top 5 lineages
Var1<-c("B.1.1.7", "P.1", "B.1.351", "B.1.617.2")
Freq<-c("","","","")
df<-data.frame(Var1, Freq)

#merge the vocs to top 5
lineageCountTopX<-rbind(lineageCountTop5, df)

lineageCountTopX<-lineageCountTopX %>% 
  dplyr::rename(
    lineage = Var1,
  )

lineageCountTopX$topLineages<-lineageCountTopX$lineage
dayta2topLineagesOnly <- merge(selected_mutationData,lineageCountTopX,by="lineage",all.x = TRUE)


highestCasesInADay<- selected_continentData1$new_smoothed_cases_per_million[!is.na(selected_continentData1$new_smoothed_cases_per_million)] %>% 
  max()

#####
##Sum of topLineages per 2 weeks#######
#data, group_by, count
biweekly_lineage_summary<-dplyr::count(dayta2topLineagesOnly, date2, topLineages) #count toplineages ever 2 weeks
head(biweekly_lineage_summary)

# #total sequences per continent per 2 weeks
total_biweekly_seqs<-biweekly_lineage_summary %>% 
  group_by(date2) %>%
  summarise(
    tot_sequences = sum(n)                   #total sequences every 2 weeks
  ) %>% 
  left_join (biweekly_lineage_summary,        #merge
             by = "date2") %>% 
  mutate(
    topLineages_prop = (n/tot_sequences)      #find proportions of top lineages
  ) %>% 
  drop_na(topLineages)                         #drop "other" lineages so they don't appear on graph


# Plot 2 ##########
#pdf(paste0("../figures/", selection, "3.pdf"), width=8, height = 4.5)
plot2<-ggplot() + 
  geom_bar(data=total_biweekly_seqs, aes(x = date2, y=topLineages_prop,fill=topLineages), #fill=forcats::fct_rev(topLineages) #to reverse order
           position = "stack", stat="identity", width=14,color='black')+
  geom_line(data=selected_continentData1, aes(x=date,y= new_smoothed_cases_per_million/highestCasesInADay,color='data'), size=1,color="black") +
  #geom_line(data=comprehensive_df3, aes(x=date, y=allNoVaccine_smoothed/highestCasesInADay, color="Fitted - no vaccines"), size=1, color="grey")+
  #geom_line(data=comprehensive_df3, aes(x=date, y=all_smoothed/highestCasesInADay, color="Fitted all"), size=1, color="brown")+
  theme_classic() + theme_bw() +
  ylab('Proportion') +
  scale_y_continuous(breaks = c(0,0.25, 0.50, 0.75, 1),
                     labels = scales::comma,
                     sec.axis = sec_axis(~.*highestCasesInADay, name = "Cases per million", labels = scales::number))+
  theme(legend.position = "right") +
  scale_x_date(date_labels = "%b-%y",date_breaks = "1 month", limits = as.Date(c('2019-12-01','2021-08-01')))+
  #theme(axis.title.x = element_text(color="black", size=15, face="bold"))+
  theme(axis.text.x = element_text(angle=90, color="black", size=13),
        axis.title.x = element_blank(),
        axis.title.y = element_text(color="black", size=15, face="bold"),
        axis.text.y = element_text(color="black", size=14),
        plot.title = element_text(size = 15, face = "bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+ 
  scale_fill_brewer(palette = "Set3")+
  ggtitle(selection)+
  guides(fill=guide_legend(title="Lineage"))

#dev.off()

#Arrange plots in set##########

#pdf(paste0("../figures/", selection, "3.pdf"), width=8, height = 4.5)
plot2 + inset_element(plot1, left = 0.01,right = 0.5, bottom = 0.6, top = 0.99)
#dev.off()

#clear memory
rm(list = ls())
dev.off()
