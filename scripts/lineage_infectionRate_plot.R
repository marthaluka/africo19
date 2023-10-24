#Regression Analysis
# This script does data pre-processing and model selection of predictors of infection rates

# Read data
# source("readData.R")

require(pacman)
pacman::p_load(flextable)

# Continent plot ################
  ###Plot topLineage dynamics and infection rates
head(continentData1)
head(lineageData)

#rolling average to smooth
continentData2 <- continentData1 %>% 
  group_by(continent) %>% 
  arrange(date) %>% 
  mutate(
    smoothed_cases = zoo::rollmean(cases, k = 14, fill = "extend"), # rolling mean - for smoothing
    smoothed_tests = zoo::rollmean(tests, k = 14, fill = "extend"),
    smoothed_vaccine_doses = zoo::rollmean(vaccine_doses, k = 14, fill = "extend"),
    stringency_index = zoo::rollmean(govt_index, k = 4, fill = "extend") 
    ) %>% 
  dplyr::select(-c(cases, tests, vaccine_doses, govt_index))


  #Top X (5 in this case) lineages per continent
summ_lineages<- lineageData %>%
  na.omit() %>%
  group_by(continent, lineage) %>%
  dplyr::summarise(n=n()) %>%
  arrange(desc(n)) %>%
  slice(1:5) 

topLineages<-unique(summ_lineages$lineage) 
topLineages %>% sort() #use this to order the legend, with vocs coming first

  #Create table of top lineages
table_df <- lineageData %>%
  rename(Lineage=lineage) %>%
  dplyr::filter(Lineage %in% topLineages) %>%
  group_by(Lineage) %>%
  dplyr::summarise(Sum=n()) %>%
  ungroup() %>%
  arrange(desc(Sum))

table_df$Proportion <- table_df$Sum/nrow(lineageData) %>% round(.,digits=3)
table_df$Proportion <- table_df$Proportion  %>% round(.,digits=3)
sum(table_df$Proportion)

tot_sequences<-sum(table_df$Sum)
tot_proportion <- sum(table_df$Proportion)
table_df[nrow(table_df) + 1,] = list("Total", as.integer(tot_sequences), as.numeric(tot_proportion))

table_df %>%
  flextable::flextable() %>%
  flextable::autofit() %>%
  flextable::save_as_docx(path = "figures/top_lineages.docx")
  
#create a df  - frequency of top lineages 
df_topLineages<- lineageData %>%
  dplyr::group_by(continent, lineage, date14)  %>%  #date21
  dplyr::summarise(topTotal=n()) %>%
  subset(lineage %in% topLineages)  #filter based on top
  

df_biweeklyTotal<-lineageData %>%
  drop_na(.)  %>%
  dplyr::group_by(continent, date14)  %>% #date21
  dplyr::summarise(allSeqs=n()) %>%
  left_join(df_topLineages,
            by=c("date14", "continent")) %>% #date21
  drop_na(.)

df_biweeklyTotal$topLineage_prop <- df_biweeklyTotal$topTotal/df_biweeklyTotal$allSeqs  

# Create a custom pallete to fix lineage colours across continents  #############
my_colour_palette<-c("#8DD3C7","#FFFFB3","royalblue","#FB8072","#80B1D3",
                     "#B3DE69","#FCCDE5","deeppink3","#BC80BD","#CCEBC5","darkorange", 
                     "goldenrod1","#BEBADA",
                     "peachpuff","gray40", "cyan", "deepskyblue3", "tan4","gray95", "firebrick1")
#order legend
lineage_levels<-c("Alpha", "Beta", "Gamma", "Delta", "Omicron", "B.1", "B.1.1.214", "B.1.1.284",
                  "B.1.177", "B.1.2", "B.1.621", "C.37", "D.2")


highestCasesInADay<-max(continentData2$smoothed_cases)

df_biweeklyTotal$lineage<-factor(df_biweeklyTotal$lineage, levels = lineage_levels)

table(df_biweeklyTotal$lineage)

head(continentData2)






#######################
#summary #not part of workflow

#1. length of data
tmp<-lineageData %>%
  dplyr::filter(date >= "2020-01-01") 

tmp%>%
  dplyr::summarise(n())

#2. unique lineages in all data
unique(tmp$lineage) %>% length()

#3. Proportion of sequences (read infections) of topLineages
(tmp %>%
  subset(lineage %in% topLineages) %>%  dplyr::summarise(n())) / (tmp%>%
                                                                    dplyr::summarise(n()))

#######################


# Plot ##########
plot1 <-ggplot() +  #date21
  geom_bar(data=df_biweeklyTotal, aes(x = date14, y=topLineage_prop,fill=lineage), #fill=forcats::fct_rev(topLineages) #to reverse order
           position = "stack", stat="identity", width=14,color='black') + #date21
  geom_line(data=continentData2, aes(x=date,y= smoothed_cases/highestCasesInADay,color='data'), size=1,color="black") +
  theme_bw() +
  ylab('Proportion') +
  scale_y_continuous(breaks = c(0,0.25, 0.50, 0.75, 1),
                     labels = scales::comma,
                     sec.axis = sec_axis(~.*highestCasesInADay, name = "Cases per million", 
                                         labels = scales::unit_format(unit = "k", scale = 1e-3)))+
  scale_x_date(date_labels = "%b \n %Y",date_breaks = "4 months", limits = as.Date(c('2019-12-01','2022-02-02')))+
  #theme(axis.title.x = element_text(color="black", size=15, face="bold"))+
  theme(axis.text.x = element_text(angle=0, color="black", size=13),
        axis.title.x = element_blank(),
        axis.title.y = element_text(color="black", size=15, face="bold"),
        axis.text.y = element_text(color="black", size=14),
        plot.title = element_text(size = 10, face = "bold"),
        legend.position = "right",
        legend.text=element_text(size=13),
        legend.title=element_text(size=15),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = 15))+ 
  scale_fill_manual(values = my_colour_palette, name ="Lineage")+
  facet_wrap(vars(continent), ncol = 2)+
  guides(fill=guide_legend(title="Lineage"))

plot1

# pdf(paste0("./figures/Figure1.pdf"), width=9, height = 10)
# plot1
# dev.off()



rm(list = c("summ_lineages", "table_df", "df_biweeklyTotal", "highestCasesInADay"))

# Country plot #######

head(owidData1)
head(lineageData)

# summarise data
#Top X (5 in this case) lineages per continent
summ_lineages<- lineageData %>%
  dplyr::filter(country %in% countries_with_data) %>%
  na.omit() %>%
  group_by(country, lineage) %>%
  dplyr::summarise(n=n()) %>%
  arrange(desc(n)) %>%
  slice(1:5) 

# use same top lineages identified at cont level
topLineages %>% sort() #use this to order the legend, with vocs coming first

#Create table of top lineages

table_df <- lineageData %>%
  dplyr::filter(country %in% countries_with_data) %>%
  rename(Lineage=lineage) %>%
  dplyr::filter(Lineage %in% topLineages) %>%
  group_by(Lineage) %>%
  dplyr::summarise(Sum=n()) %>%
  ungroup() %>%
  arrange(desc(Sum))

table_df$Proportion <- table_df$Sum/nrow(lineageData %>%
                                           dplyr::filter(country %in% countries_with_data)) %>% 
  round(.,digits=3)
table_df$Proportion <- table_df$Proportion  %>% round(.,digits=3)
sum(table_df$Proportion)

tot_sequences<-sum(table_df$Sum)
tot_proportion <- sum(table_df$Proportion)
table_df[nrow(table_df) + 1,] = list("Total", as.integer(tot_sequences), as.numeric(tot_proportion))

#create a df  - frequency of top lineages 
df_topLineages<- lineageData %>%
  dplyr::filter(country %in% countries_with_data) %>%
  dplyr::group_by(country, lineage, date14)  %>%  #date21
  dplyr::summarise(topTotal=n()) %>%
  subset(lineage %in% topLineages)  #filter based on top


df_biweeklyTotal<-lineageData %>%
  dplyr::filter(country %in% countries_with_data) %>%
  drop_na(.)  %>%
  dplyr::group_by(country, date14)  %>% #date21
  dplyr::summarise(allSeqs=n()) %>%
  left_join(df_topLineages,
            by=c("date14", "country")) %>% #date21
  drop_na(.)  %>%
  ungroup()  %>%
  dplyr::mutate(
    topLineage_prop = topTotal/allSeqs,
    lineage = factor(lineage, levels = lineage_levels)
  )

# df_biweeklyTotal$topLineage_prop <- df_biweeklyTotal$topTotal/df_biweeklyTotal$allSeqs  
# 
# df_biweeklyTotal$lineage<-factor(df_biweeklyTotal$lineage, levels = lineage_levels)

table(df_biweeklyTotal$lineage)

countryData <- owidData1 %>%
  dplyr::filter(location%in% countries_with_data)%>%
  dplyr::rename(country=location)

highestCasesInADay<-max(countryData$smoothed_cases_per_million)



# plot
plot2<-ggplot() + 
  geom_bar(data=df_biweeklyTotal, aes(x = date14, y=topLineage_prop,fill=lineage), #fill=forcats::fct_rev(topLineages) #to reverse order
           position = "stack", stat="identity", width=14,color='black') + #date21
  geom_line(data=countryData, aes(x=date,y= smoothed_cases_per_million/highestCasesInADay,color='data'), size=1,color="black") +
  theme_bw() +
  ylab('Proportion') +
  scale_y_continuous(breaks = c(0,0.25, 0.50, 0.75, 1),
                     labels = scales::comma,
                     sec.axis = sec_axis(~.*highestCasesInADay, name = "Cases per million", 
                                         labels = scales::unit_format(unit = "k", scale = 1e-3)))+
  scale_x_date(date_labels = "%b \n %Y",date_breaks = "6 months", limits = as.Date(c('2019-12-01','2022-02-02')))+
  #theme(axis.title.x = element_text(color="black", size=15, face="bold"))+
  theme(axis.text.x = element_text(angle=0, color="black", size=12),
        axis.title.x = element_blank(),
        axis.title.y = element_text(color="black", size=15, face="bold"),
        axis.text.y = element_text(color="black", size=14),
        plot.title = element_text(size = 10, face = "bold"),
        legend.position = "right",
        legend.text=element_text(size=13),
        legend.title=element_text(size=15),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = 15))+ 
  scale_fill_manual(values = my_colour_palette, name ="Lineage")+
  facet_wrap(vars(country))+
  guides(fill=guide_legend(title="Lineage"))

plot2


pdf(paste0("./figures/FigureA1.pdf"), width=14, height = 10)
plot2
dev.off()



