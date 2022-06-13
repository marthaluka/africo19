#Regression Analysis
# This script does data pre-processing and model selection of predictors of infection rates

# Read data
source("readData.R")

# rm(list=setdiff(ls(), c("mutationData", "continentData1", "regressionData", "models")))
head(continentData1)

###Plot topLineage dynamics and infection rates################################################
head(continentData1)
head(mutationData)

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


ggplot(continentData2)+
  geom_line(aes(x=date, y=stringency_index))+
  facet_wrap(vars(continent))



#Top X (5 in this case) lineages per continent
summ_lineages<- mutationData %>%
  na.omit() %>%
  group_by(continent, lineage) %>%
  dplyr::summarise(n=n()) %>%
  arrange(desc(n)) %>%
  slice(1:5) 

head(summ_lineages)
topLineages<-unique(summ_lineages$lineage) 
topLineages %>% sort() #use this to order the legend, with vocs coming first

#Create table of top lineages
library(flextable)      #to export tables

table_df <- mutationData %>%
  rename(Lineage=lineage) %>%
  dplyr::filter(Lineage %in% topLineages) %>%
  group_by(Lineage) %>%
  dplyr::summarise(Sum=n()) %>%
  ungroup() %>%
  arrange(desc(Sum))

table_df$Proportion <- table_df$Sum/nrow(mutationData) %>% round(.,digits=3)
table_df$Proportion <- table_df$Proportion  %>% round(.,digits=3)
sum(table_df$Proportion)

tot_sequences<-sum(table_df$Sum)
tot_proportion <- sum(table_df$Proportion)
table_df[nrow(table_df) + 1,] = list("Total", as.integer(tot_sequences), as.numeric(tot_proportion))

table_df %>%
  flextable::flextable() %>%
  flextable::autofit() %>%
  flextable::save_as_docx(path = "figures/top_lineages.docx")
  

# Create a custom pallete to fix lineage colours across continents  #############
my_colour_palette<-c("#8DD3C7","#FFFFB3","royalblue","#FB8072","#80B1D3",
                  "#B3DE69","#FCCDE5","deeppink3","#BC80BD","#CCEBC5","darkorange", 
                  "goldenrod1","#BEBADA",
                  "peachpuff","gray40", "cyan", "deepskyblue3", "tan4","gray95", "firebrick1")
#order legend
lineage_levels<-c("Alpha", "Beta", "Gamma", "Delta", "Omicron", "B.1", "B.1.1.214", "B.1.1.284",
          "B.1.177", "B.1.2", "B.1.621", "C.37", "D.2")

#create a df  - frequency of top lineages 
df_topLineages<- mutationData %>%
  dplyr::group_by(continent, lineage, date14)  %>%  #date21
  dplyr::summarise(topTotal=n()) %>%
  subset(lineage %in% topLineages)  #filter based on top
  

df_biweeklyTotal<-mutationData %>%
  drop_na(.)  %>%
  dplyr::group_by(continent, date14)  %>% #date21
  dplyr::summarise(allSeqs=n()) %>%
  left_join(df_topLineages,
            by=c("date14", "continent")) %>% #date21
  drop_na(.)

df_biweeklyTotal$topLineage_prop <- df_biweeklyTotal$topTotal/df_biweeklyTotal$allSeqs  

df_biweeklyTotal$lineage<-factor(df_biweeklyTotal$lineage, levels = lineage_levels)

table(df_biweeklyTotal$lineage)

head(continentData2)
highestCasesInADay<-max(continentData2$smoothed_cases)

#######################
#summary #not part of workflow

#1. length of data
tmp<-mutationData %>%
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
#pdf(paste0("../figures/", selection, "3.pdf"), width=8, height = 4.5)
plot2<-ggplot() +  #date21
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

plot2
#dev.off()


#run regression_continent.R script at this stage to end up with ggplot objects in the environment

# source("regression_continent.R")

#Arrange plots in set##########

africa <- plot2 + inset_element(Africa, left = 0.045,right = 0.25, bottom = 0.88, top = 0.993)
africa_asia <- africa + inset_element(Asia, left = 0.5,right = 0.75, bottom = 0.88, top = 0.993)
africa_asia_europe <- africa_asia + inset_element(Europe, left =  -0.005,right = 0.24, bottom = 0.52, top = 0.642)
africa_asia_europe_namerica <- africa_asia_europe + inset_element(North_America, left = 0.5,right = 0.75, bottom = 0.52, top = 0.642)
africa_asia_europe_namerica_oceania<-africa_asia_europe_namerica + inset_element(Oceania, left = -0.005,right = 0.172, bottom = 0.183, top = 0.297)
all <-africa_asia_europe_namerica_oceania + inset_element(South_America, left = 0.5,right = 0.75, bottom = 0.18, top = 0.297)
all


pdf(paste0("./figures/Figure_1.pdf"), width=12, height = 12)
all
dev.off()

# #clear memory
# rm(list = ls())
# pacman::p_unload(negate=T)
# dev.off()
