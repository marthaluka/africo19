#Regression Analysis
# This script does data pre-processing and model selection of predictors of infection rates

#Read data

rm(list = ls())

## packages
require("pacman")
pacman::p_load(data.table, tidyverse, zoo, RColorBrewer, directlabels, 
               stringr, spData, ggrepel, corrr, patchwork, cowplot, egg, here) 


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

#Pre-processing
#biweekly rolling mean - for smoothing

owidData1 <- owidData0 %>%
  mutate(across(where(is.character),str_trim)) %>% 
  dplyr::filter(new_cases_per_million >= 0) %>%    #drop any negative/zero values ie treat zero as missing data
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

countryLevelData_owid<-owidData1 %>% 
  subset(., ! location %in% removeList) #ensure only country-level entries


continentData0<-subset(owidData1, location %in% continentList) %>%   #continent level data
  dplyr::select(-c(stringency_index)) #drop the country-specific stringency field. Recalculate for continent-level below

continentData0$continent<-continentData0$location

###calculate average stringency index per continent and add to continent data#######
contStringency<-countryLevelData_owid %>%
  dplyr::select(continent, date, stringency_index)%>%
  mutate_at("stringency_index", ~replace(., is.na(.), 0)) %>%  #replace NA with zero
  group_by(continent, date) %>% 
  #summarise_each(funs(mean)) %>%  #continental average of stringency index. #using next command as this is deprecated
  summarise_each(list(stringency_index = mean)) %>%   #continental average of stringency index
  subset(continent %in% continentList)  # drop unknown/misspelt contient names

continentData1<-merge(continentData0, contStringency, by=c("continent", "date"), all.x = T)


continentData1<-continentData1%>%
  group_by(continent) %>%
  mutate(yaxis = smoothed_cases_per_million/(max(smoothed_cases_per_million)))


head(continentData1)

#visualize data to check if all is ok
ggplot(data=continentData1)+
  #geom_line(aes(x=date, y=new_cases_per_million, colour = "Cases per 1M"))+
  geom_line(aes(x=date, y=stringency_index, colour = "Stringency Index"))+
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
  mutate(across(where(is.character),str_trim)) %>%
  dplyr::filter(continent != is.na(.))

        #drop countries with less than 100 sequences
        # table1<-table(mutationData$country)
        # mutationData<- mutationData %>%
        #   subset(.,country %in% names(table1[table1>100])) #drop

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

table(mutationData$lineage)


###Plot topLineage dynamics and infection rates################################################
head(continentData1)
head(mutationData)



#Top X (5 in this case) lineages per continent
summ_lineages<- mutationData %>%
  na.omit() %>%
  group_by(continent, lineage) %>%
  dplyr::summarise(n=n()) %>%
  arrange(desc(n)) %>%
  slice(1:5) 

topLineages<-unique(summ_lineages$lineage)


#create a custom pallete to fix lineage colours across continents  #############
my_colour_palette<-c("#8DD3C7","#FFFFB3","royalblue","#FB8072","#80B1D3",
                  "#B3DE69","#FCCDE5","deeppink3","#BC80BD","#CCEBC5","darkorange", 
                  "goldenrod1","#BEBADA",
                  "peachpuff","gray40", "cyan", "deepskyblue3", "tan4","gray95", "firebrick1")

lineage_levels<-c("Alpha", "Beta", "Gamma", "Delta", "A.2.2", "B.1", "B.1.1", "B.1.1.214", "B.1.1.284",
          "B.1.177", "B.1.2", "B.1.351.2", "B.1.429", "C.37", "D.2", "P.2", "R.1")

df_topLineages<- mutationData %>%
  dplyr::group_by(continent, lineage, date14)  %>%
  dplyr::summarise(topTotal=n()) %>%
  subset(lineage %in% topLineages)  #filter based on top
  

df_biweeklyTotal<-mutationData %>%
  drop_na(.)  %>%
  dplyr::group_by(continent, date14)  %>%
  dplyr::summarise(allSeqs=n()) %>%
  left_join(df_topLineages,
            by=c("date14", "continent"))

df_biweeklyTotal$topLineage_prop <- df_biweeklyTotal$topTotal/df_biweeklyTotal$allSeqs  

df_biweeklyTotal$lineage<-factor(df_biweeklyTotal$lineage, levels = lineage_levels)

head(continentData1)

highestCasesInADay<-max(continentData1$smoothed_cases_per_million)

# Plot ##########
#pdf(paste0("../figures/", selection, "3.pdf"), width=8, height = 4.5)
plot2<-ggplot() + 
  geom_bar(data=df_biweeklyTotal, aes(x = date14, y=topLineage_prop,fill=lineage), #fill=forcats::fct_rev(topLineages) #to reverse order
           position = "stack", stat="identity", width=14,color='black')+
  geom_line(data=continentData1, aes(x=date,y= smoothed_cases_per_million/highestCasesInADay,color='data'), size=1,color="black") +
  theme_classic() + theme_bw() +
  ylab('Proportion') +
  scale_y_continuous(breaks = c(0,0.25, 0.50, 0.75, 1),
                     labels = scales::comma,
                     sec.axis = sec_axis(~.*highestCasesInADay, name = "Cases per million", labels = scales::number))+
  theme(legend.position = "right") +
  scale_x_date(date_labels = "%b \n %Y",date_breaks = "3 months", limits = as.Date(c('2019-12-01','2021-10-01')))+
  #theme(axis.title.x = element_text(color="black", size=15, face="bold"))+
  theme(axis.text.x = element_text(angle=0, color="black", size=13),
        axis.title.x = element_blank(),
        axis.title.y = element_text(color="black", size=15, face="bold"),
        axis.text.y = element_text(color="black", size=14),
        plot.title = element_text(size = 10, face = "bold"),
        legend.text=element_text(size=13),
        legend.title=element_text(size=15),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = 15))+ 
  scale_fill_manual(values = my_colour_palette, name ="Lineage")+
  facet_wrap(vars(continent), ncol = 3)+
  guides(fill=guide_legend(title="Lineage"))

plot2
#dev.off()


#run regression_continent.R script at this stage to end up with ggplot objects in the environment

#Arrange plots in set##########

africa <- plot2 + inset_element(Africa, left = 0.018,right = 0.25, bottom = 0.88, top = 0.993)
africa_asia <- africa + inset_element(Asia, left = 0.5,right = 0.75, bottom = 0.88, top = 0.993)
africa_asia_europe <- africa_asia + inset_element(Europe, left =  -0.01,right = 0.24, bottom = 0.54, top = 0.642)
africa_asia_europe_namerica <- africa_asia_europe + inset_element(North_America, left = 0.5,right = 0.75, bottom = 0.54, top = 0.642)
africa_asia_europe_namerica_oceania<-africa_asia_europe_namerica + inset_element(Oceania, left = -0.01,right = 0.175, bottom = 0.183, top = 0.297)
all <-africa_asia_europe_namerica_oceania + inset_element(South_America, left = 0.5,right = 0.75, bottom = 0.18, top = 0.297)
all


pdf(paste0("./figures/Figure_1.pdf"), width=14, height = 10)
all
dev.off()

# #clear memory
# rm(list = ls())
# dev.off()
