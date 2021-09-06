#Regression Analysis#######################################################
#OWID#########

#setwd
setwd("~/Documents/CVR/data")
rm(list = ls())

##required packages
require("pacman")
pacman::p_load(dplyr, ggplot2, data.table, tidyverse, zoo, RColorBrewer, directlabels, 
               stringr, spData, ggrepel, corrr, patchwork, cowplot, egg, here) 

#Non-mutation data##### 
owidData0<- read.csv(
  "https://raw.githubusercontent.com/owid/covid-19-data/master/public/data/owid-covid-data.csv",
  fileEncoding="UTF-8-BOM", stringsAsFactors = FALSE
) %>% #select relevant fields 
  dplyr::select(continent, location, date, total_cases, new_cases, total_cases_per_million, 
                new_cases_per_million, total_vaccinations, total_vaccinations_per_hundred, 
                population, stringency_index)
head(owidData0)

##biweekly rolling mean - for smoothing####### Also lag stringency/ lead cases
owidData0$date<-as.Date(owidData0$date)
owidData1 <- owidData0 %>%
  filter(new_cases>0) %>%    #drop any negative/zero values ie treat zero as missing data
  dplyr::group_by(location) %>% 
  dplyr::arrange(date) %>% 
  dplyr::mutate(cases_14da = zoo::rollmean(new_cases, k = 14, fill = NA)) %>%   #biweekly rolling mean - for smoothing
  dplyr::mutate(new_smoothed_cases_per_million= (cases_14da/population)*1000000) %>% # smoothed infection rates per 1M
  #dplyr::mutate(cases_lead_3days = lag (new_smoothed_cases_per_million, n=3)) %>%  #"lead" cases by 3 days of the stringency index
  #dplyr::mutate(cases_lead_7days = lag (new_smoothed_cases_per_million, n=7)) %>%  #"lead" cases by 7 days of the stringency index
  #dplyr::mutate(cases_lead_14days = lag (new_smoothed_cases_per_million, n=14)) %>%  #"lead" cases by 14 days of the stringency index
  dplyr::ungroup()

head(owidData1)
#write_csv(owidData2, "shinyData.csv")

owidData1$date2<-as.Date(cut(owidData1$date,breaks = "2 weeks",start.on.monday = FALSE))
#clean up data. We have both continent-level and country-level entries. 
removeList<-c("World","Africa","Asia", "Oceania", "Africa", "North America", "South America", 
              "Europe", "European Union", "International")  

continentList<-c("Africa","Asia","Oceania","Africa","North America","South America","Europe")

##drop/select for continents#######
owidData2<-subset(owidData1, ! location %in% removeList) #ensure only country-level entries
head(owidData2)

continentData0<-subset(owidData1, location %in% continentList) %>%   #continent level data
  select(-c(stringency_index)) #drop the country specific stringency field. Recalculate for continent-level below
continentData0$continent<-continentData0$location

###calculate average stringency index per continent and add to continent data#######
contStringency<-owidData2 %>%
  dplyr::select(continent, date, stringency_index)%>%
  mutate_at("stringency_index", ~replace(., is.na(.), 0)) %>% 
  group_by(continent, date) %>% 
  summarise_each(funs(mean))%>% 
  subset(continent %in% continentList)  #drop NAs rows within the continent field

continentData1<-merge(continentData0, contStringency, by=c("continent", "date"), all.x = T)
head(continentData1)

#visualize data to check if all is ok
ggplot(data=continentData1)+
  #geom_line(aes(x=date, y=new_cases_per_million, colour = "Cases per 1M"))+
  geom_line(aes(x=date, y=stringency_index, colour = "Stringency Index"))+
  #geom_line(aes(x=date, y=cases_lead_3days, colour= "cases_lead_3days"))+
  #geom_line(aes(x=date, y=cases_lead_7days, colour= "cases_lead_7days"))+
  #geom_line(aes(x=date, y=cases_lead_14days, colour= "cases_lead_14days"))+
  geom_line(aes(x=date, y=new_smoothed_cases_per_million, colour = "Smoothed Cases per 1M"))+
  #scale_x_date(date_labels = "%b-%y",date_breaks = "1 month", limits = as.Date(c('2019-12-01','2021-07-01')))+
  facet_wrap(vars(continent))+ ylab(NULL)+
  theme_bw()


#GISAID data##########
#read data 
mutationData<-read.csv("./alpha_data/metadata29062021.csv", header = F) %>%
  dplyr::select(V2, V10, V5, V8, V11) 

column_headings <- c("gisaid_ID", "country", "lineage", "date", "continent")
names(mutationData)<-column_headings
mutationData$date<-as.Date(mutationData$date)
head(mutationData)
unique(mutationData$lineage) %>% length()

mutationData$continent <- mutationData$continent %>% gsub("_", " ", .) %>% 
  str_to_title(mutationData$continent) #replace underscore with space & change continent to sentence case

pacman::p_load(plyr) 

##Distinct lineages #################################
mutationData1<-mutationData %>% group_by(continent, date) %>% 
  ddply(~date+continent,summarise,distinct_lineages_per_day=length(unique(lineage))) %>% 
  ungroup()  %>% 
  left_join (mutationData, 
             by= c("continent", "date"))

pacman::p_unload(plyr) #will conflict with dplyr

#define VOCs 
vocs<-c("B.1.1.7", "P.1", "B.1.351", "B.1.617.2")

##Summary of lineages per country/location per day#######
dplyr::count(mutationData1, continent, date, lineage) -> daily_lineage_summary

# #count vocs per continent per day
vocs_spec <- daily_lineage_summary %>% 
  filter(lineage %in% vocs) %>% 
  group_by(continent, date) %>% 
  summarise(
    voc_total = sum(n, na.rm = T)
  )

# #total sequences per continent per day
data_set3 <- daily_lineage_summary %>% 
  left_join(vocs_spec,
            by = c("continent", "date")) %>% 
  group_by(continent, date) %>% 
  summarise(
    tot_sequences = sum(n)
  )

#View(data_set3)
# #proportion of vocs per continent per day
voc_proportions <- data_set3 %>% 
  left_join(vocs_spec,
            by = c("continent", "date")) %>% 
  mutate(
    voc_total = replace(voc_total, which(is.na(voc_total)), 0),
    voc_prop = (voc_total/tot_sequences)*100
  )

#merge voc proportions with OWID data
head(continentData1)
head(voc_proportions)
mergedData<-merge(continentData1, voc_proportions, by=c("continent", "date")) %>%
  select(continent, date, new_smoothed_cases_per_million, stringency_index, voc_prop, 
         total_vaccinations_per_hundred) #cases_lead_3days, cases_lead_7days, cases_lead_14days)
head(mergedData)


##Create df #################################
#merge voc proportions with OWID data
head(continentData1)
head(mutationData1)
head(mergedData)

table(mutationData1$continent)

comprehensive_df<-left_join(mergedData, mutationData1, by= c("continent", "date")) %>%
  select(-c(gisaid_ID, country, lineage))  %>% unique()
head(comprehensive_df)  #for regression

#Feature scaling prior to regression analysis########### Not done as differences not too big?
### Predictor variables at this point are (i) VOC proportion (range 0-100), 
# (ii) distinct lineages per day (range 1- ~250), (iii) stringency index (range 0-100) and 
# (iv) total vaccination per 100 people (range 0 - ~100)


#Regression###########
###Regression I - Stringency index #################################
##lm([target variable] ~ [predictor variables], data = [data source])
head(comprehensive_df)
table(comprehensive_df$continent)

selected<-"South America"
comprehensive_df2<-comprehensive_df %>%
  filter(continent==selected)

clearList<- ls() %>% str_subset("^lmInfections") 
rm(list = clearList)

#all combined minus vaccination
lmInfectionsX <- lm(new_smoothed_cases_per_million ~ voc_prop+distinct_lineages_per_day+
                      stringency_index, 
                    data = comprehensive_df2)
sm_X<-summary(lmInfectionsX)
sm_X
mean(sm_X$residuals^2)  #root mean square error

#all combined including vaccination
lmInfectionsXX <- lm(new_smoothed_cases_per_million ~ voc_prop+distinct_lineages_per_day+
                       stringency_index+total_vaccinations_per_hundred, 
                     data = comprehensive_df2)
sm_XX<-summary(lmInfectionsXX)
sm_XX
mean(sm_XX$residuals^2)  #root mean square error

#get predicted values#####################
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

table(mutationData1$continent)
table(continentData1$continent)

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
