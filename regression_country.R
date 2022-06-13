#Regression Analysis at country-level
# This script does data pre-processing and model selection of predictors of infection rates



## packages
require("pacman")
pacman::p_load(data.table, tidyverse, zoo, RColorBrewer, directlabels,corrplot, broom,
               stringr, spData, ggrepel, corrr, patchwork, cowplot, egg, here,
               vroom) 

#Read data
source("readData.R")

table(mutationData$country)
table(owidData0$location)
# "OWID, GISAID"
# "Czechia", "Czech Republic"
# "Guyana", "French Guiana"
# "Democratic Republic of Congo", "Democratic Republic of the Congo"
# "Gambia", "The Gambia"
# "Sint Maarten (Dutch part)","Sint Maarten"

head(mutationData)

mutationData <- mutationData[, if(.N >= 1000) .SD, by=country] #drop countries with less than 1k sequences

##1. Virus Diversity ie Distinct lineages per time period#################################
#diversity
#because some regions sequence more, distinct lineages may biased to more sequences. 
#We divide distinct lineages by total number of sequences to achieve a #diversity_score
# We will smooth the numbers afterwards

#total number of seqs per country per week
seq_count <- dplyr::count(mutationData, country, date7)


pacman::p_load(plyr) 

country_mutationData1<- mutationData%>%
  group_by(country, date7) %>% 
  ddply(~date7+country,summarise,distinct_lineages_per_country_per_week=length(unique(lineage))) %>% 
  ungroup() 

lineage_diversity_df<- seq_count  %>%
  dplyr::rename(tot_genomes_per_week = n) %>%
  left_join (country_mutationData1, 
             by= c("country", "date7")) %>%
  mutate(
    diversity_score=distinct_lineages_per_country_per_week/tot_genomes_per_week,  #diversity taking into account total number of sequences. 
    # drop data before 1st April 2020 before the linear regression
  )  %>%
  dplyr::select(country, date7, diversity_score, distinct_lineages_per_country_per_week, tot_genomes_per_week) %>%
  distinct() %>%    #drop duplicated rows
  group_by(country) %>%
  dplyr::arrange(date7) %>%
  mutate(
    diversity_score = zoo::rollmean(diversity_score, k = 3, fill = "extend")  
  ) %>%  
    dplyr::ungroup()

pacman::p_unload(plyr) #will conflict with dplyr downstream



##2. VOC proportions per time period########

#define VOCs 
vocs<-c("Alpha", "Beta", "Gamma", "Delta")

### Summary of lineages per country per day
dplyr::count(mutationData, country, date7, lineage) -> country_weekly_lineage_summary

              # # #count vocs per country per day
              # vocs_count <- country_weekly_lineage_summary %>% 
              #   dplyr::filter(lineage %in% vocs) %>% 
              #   group_by(country, date7) %>% 
              #   summarise(
              #     voc_total = sum(n, na.rm = T)
              #   )


# #total sequences per country per day
data_set_country <- country_weekly_lineage_summary %>% 
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





##2. Fitness defined using mutations #####

pacman::p_load(data.table, tidyverse, vroom, parallel)

clean_mutations<-vroom::vroom("../data/mutations_clean_file2.csv") %>%
  select(Gisaid, Date, Country, NewMutation)

# tibble to data.table for ? faster wrangling
clean_mutations<-data.table(clean_mutations)

# Date
clean_mutations$Date <- as.Date(clean_mutations$Date)
clean_mutations$date7<-as.Date(cut(clean_mutations$Date,breaks = "1 week",
                                   start.on.monday = TRUE))

#total genomes per country per week
tmp <- clean_mutations %>% 
  dplyr::select(Gisaid, Country, date7) %>% 
  distinct()

dplyr::count(tmp, Country, date7) -> genome_count_per_cont_per_week

rm(list=c('tmp', 'continentList', 'column_headings', 'dict', 'removeList', 'seq_count', 'continent_weekly_mutations'))

#total mutations per continent week
dplyr::count(clean_mutations, Country, date7) -> country_weekly_mutations

#specific mutations count per cont per week
dplyr::count(clean_mutations, Country, date7, NewMutation) -> country_weekly_specific_mutation_summary


# Mutations from www.medrxiv.org/content/10.1101/2021.09.07.21263228v2.full
mutations_increasing_fitness <- str_to_upper(c("S:H655Y", "S:T95I", "ORF1a:P3395H", "S:N764K", "ORF1a:K856R", "S:S371L",
                                               "E:T9I", "S:Q954H", "ORF9b:P10S", "S:L981F", "N:P13L", "S:G339D", "S:S375F",
                                               "S:S477N", "S:N679L", "S:S373P", "M:Q19E", "S:D796Y", "S:N969K", "S:T547K",
                                               "ORF1b:I1566V", "M:D3G", "S:G446S", "S:N440K", "M:A63T", "S:N856K", "ORF1a:A2710T",
                                               "ORF1a:I3758V", "S:E484A", "S:A67V", "S:K417N", "S:Q493R", "S:N501Y","S:Y505H", 
                                               "S:L452R", "S:P681H", "S:Q498R","S:G496S", "ORF1a:T3255I", "ORF14:G50W", "S:P681R",
                                               "N:R203M", "ORF1b:P1000L", "ORF1a:P2287S", "M:I82T", "ORF3a:S26L", "N:D63G", "N:G215C",
                                               "ORF1a:V3718A", "ORF9b:T60A"))

mutations_of_interest <- mutations_increasing_fitness


#fit_mutations_count
fit_mutations_count <- country_weekly_specific_mutation_summary %>% 
  dplyr::filter(NewMutation %in% mutations_of_interest) %>% 
  group_by(Country, date7) %>% 
  dplyr::summarise(
    weekly_sum_of_fit_mutations = sum(n, na.rm = T)
  )


fitness_df <- fit_mutations_count %>% 
  left_join(country_weekly_mutations, by = c("Country", "date7")) %>% 
  dplyr::rename(total_mutations_per_week = n) %>% 
  left_join(genome_count_per_cont_per_week, by = c("Country", "date7")) %>% 
  dplyr::rename(total_seqs_per_week = n) %>% 
  dplyr::mutate(#fitness_score1 = weekly_sum_of_fit_mutations/total_seqs_per_week,
    #fitness_score2 = weekly_sum_of_fit_mutations/total_mutations_per_week,
    fitness_score3 = weekly_sum_of_fit_mutations/(total_seqs_per_week + log(total_mutations_per_week)),
    #fitness_score4 = (weekly_sum_of_fit_mutations/total_seqs_per_week)* (1/log(total_mutations_per_week))
  )

# match country names

renameFunction<-function(locationName, na.rm=T){
  gsub("_", " ", locationName) %>%    #replace underscore with space
    str_to_title(locationName) %>%    #change country to sentence case
    gsub("And ", "and ", .)  %>%      #clean up the tough names 
    gsub("Of ", "of ", .)  %>%
    gsub(" The", " the", .) %>%
    gsub("Usa", "United States", .)%>%
    gsub("Gambia", "The Gambia", .)
}


fitness_df <- fitness_df %>% 
  dplyr::rename(country = Country)%>% 
  dplyr::mutate_at(
    c("country"), renameFunction
  ) %>%
  mutate(across(where(is.character),str_trim)) %>%
  drop_na()

####Comprehensive mutation data  ####
head(fitness_df)
head(lineage_diversity_df)

comprehensive_mutation_df<-lineage_diversity_df %>% 
  left_join(fitness_df, by=c("country", "date7")) %>% 
  dplyr::select(country, date7, diversity_score, fitness_score3, weekly_sum_of_fit_mutations, 
                total_mutations_per_week, distinct_lineages_per_country_per_week, total_seqs_per_week, )


##Create df #################################
## ie Merge gisaid with owid data
head(countryData)
head(comprehensive_mutation_df)

df_country<-countryData %>% 
  dplyr::rename(., country = location) %>% 
  left_join(comprehensive_mutation_df, by=c("country", "date7"))

rm(list=setdiff(ls(), c("mutationData", "owidData0", "df_country", "fitness_df", "countryData")))

#filter out countries before regression
  #criteria is data completeness

df_data_completeness <- df_country%>%
  dplyr::group_by(country) %>%
  dplyr::summarise(weeks_with_sequenceData = sum(!is.na(total_seqs_per_week)))#use this

df_data_completeness %>%
  dplyr::filter(
    weeks_with_sequenceData>70,
    )%>%
  ggplot()+
  geom_col(aes(y=country, x=weeks_with_sequenceData))+
  theme_bw()


table2<- df_data_completeness%>%
  dplyr::filter(
    weeks_with_sequenceData>70
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
  #geom_line(aes(x=date7, y=cases, colour = "Infection rate"))+
  geom_line(aes(x=date7, y=govt_index, colour="govt_index"))+
  facet_wrap(vars(country))


ggplot(data= df_country_filtered)+
  #geom_line(aes(x=date7, y=fitness_score3, color="fitness_score3"))+
  geom_line(aes(x=date7, y=diversity_score))+
  facet_wrap(vars(country))+
  theme_bw()


df_country_filtered <- df_country_filtered %>%
  dplyr::group_by(country) %>%
  dplyr::arrange(date7) %>%
  dplyr::mutate(
    diversity_score = zoo::rollmean(diversity_score, k = 4, fill = "extend"),  # smoothing
    fitness_score3 = zoo::rollmean(fitness_score3, k = 3, fill = "extend")  # smoothing
    ) %>%
  dplyr::ungroup()


regressionData<- df_country_filtered %>%
  dplyr::select(country, date7, cases, fitness_score3, diversity_score, govt_index,
                vaccine_doses)%>%
  dplyr::filter(between (date7, as.Date("2020-04-01"), as.Date("2022-01-10"))) %>% #data seems a bit unstable outside these dates 
  drop_na()

#all looks good! Time to fit a linear model :-)

# Correlation ####

p<-regressionData %>% 
  group_by (country) %>% 
  dplyr::select(-c(date7)) %>% 
  dplyr::rename(`Virus fitness` = fitness_score3,
                `Infection rate` = cases,
                Vaccination =vaccine_doses,
                `Virus diversity` =diversity_score,
                `Stringency index` =govt_index) %>%
  do(data.frame(Cor=t(cor(.[,2:6], .[,2], method = "pearson"))))

p %>%
  flextable::flextable() %>%
  flextable::autofit() %>%
  flextable::save_as_html(path = "figures/correlation_countries.html")

# Multiple Linear Regression
#lm([target variable] ~ [predictor variables], data = [data source])

# Use the plyr package to run regressions per country instead of filtering each time
library(plyr)
#defn function for linear model
lm_with_plyr<-function(df){
  lm(cases ~ govt_index +
       vaccine_doses +
       fitness_score3 +
       diversity_score,
     data = df)
}
#divide df into smaller dfs (creates a list of dfs) based on country and run lm for each of the dfs
models<-dlply(regressionData, "country", lm_with_plyr)
# Apply coef to each model and return a data frame
ldply(models, coef)
# Print the summary of each model
l_ply(models, summary, .print = TRUE)
# R Squared
p<-laply(models, function(mod) summary(mod)$r.squared)
p
pacman::p_unload(plyr)
# Find fitted values
summary(models$France)



fitModelFunction <- function(selected_country){
  regressionData %>%
    dplyr::filter(country==selected_country) %>%
    #get fitted infections per capita and smooth using rolling mean
    dplyr::mutate(fitted_cases = fitted(models[[selected_country]]),
           smoothed_fitted = zoo::rollmean(fitted_cases, k = 3, fill = 'extend')
    )
  }

countryList <- names(p)

fittedCountries <- lapply(countryList, fitModelFunction)
fittedCountries <- bind_rows(fittedCountries)

library(trelliscopejs)

p<-ggplot(data = fittedCountries)+
  geom_line(aes(x=date7, y=cases, color="Data"))+
  geom_line(aes(x=date7, y=smoothed_fitted, color="Fitted"))+
  theme_bw()+
  ylab('Cases') +
  scale_x_date(date_labels = "%b \n %Y",date_breaks = "6 months", limits = as.Date(c('2020-01-01','2022-01-10')))+
  theme(
    panel.background = element_rect(fill = "transparent"), # bg of the panel
    plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
    axis.text.x = element_text(color="black", size=10),
    axis.title.x = element_blank(),
    axis.title.y = element_text(color="black", size=12),
    axis.text.y = element_text(color="black", size=10),
    plot.title = element_text(size=12),
    legend.position = "top",
    legend.text = element_text(color="black", size = 10),
    legend.justification = "left",
    strip.background = element_blank(),
    # legend.margin = margin(0,0,0,0),
    # legend.box.margin = margin(-10,-10,-10,-10),
    # legend.key.width =  unit(2, "mm"),
    # legend.background = element_rect(fill = "transparent"), # get rid of legend bg
    panel.grid.major = element_blank(),
    #panel.grid.minor = element_blank()
    )+ 
  guides(color=guide_legend(title=NULL))+
  scale_color_manual(values=c("black", "violetred1")) +
  # facet_trelliscope(
  #   ~country,
  #   self_contained = TRUE
  # )

  facet_wrap(vars(country))


pdf(paste0("./figures/Supp_Fig2.pdf"), width=12, height = 12)
p
dev.off()


#########
head(mutationData)
select_mutationData<- mutationData %>%
  dplyr::filter(country %in% countryList) %>%
  dplyr::select(-c(gisaid_ID, date, date7))


top_ten_lineages <- dplyr::count(select_mutationData, lineage) %>% 
  arrange(desc(n)) %>%
  slice(1:10) 

genomes_per_fortnight_per_country <- dplyr::count(select_mutationData,country, date14) %>% 
  dplyr::rename(total_seqs = n)

lineage_specific_count_per_fortnight_per_country <- dplyr::count(select_mutationData,country, date14, lineage) %>% 
  dplyr::filter(lineage %in% top_ten_lineages$lineage) %>% 
  left_join(genomes_per_fortnight_per_country, by = c('country', 'date14')) %>% 
  dplyr::mutate(lineage_prop = n/total_seqs)


# Create a custom pallete to fix lineage colours across continents  
my_colour_palette<-c("#8DD3C7","#FFFFB3","royalblue","#FB8072","#80B1D3",
                     "#B3DE69","#FCCDE5","deeppink3","#BC80BD","#CCEBC5","darkorange", 
                     "goldenrod1","#BEBADA",
                     "peachpuff","gray40", "cyan", "deepskyblue3", "tan4","gray95", "firebrick1")

#order legend
lineage_levels<-c("Alpha", "Beta", "Gamma", "Delta", "Omicron", "B.1", "B.1.1.214", "B.1.1.284",
                  "B.1.177", "B.1.2", "B.1.621", "C.37", "D.2")

countryData_filtered<- countryData %>%
  dplyr::rename(country=location) %>%
  dplyr::filter(country %in% countryList)

highestCasesInADay<- max(countryData_filtered$cases)

another_plot<-ggplot() +  
  geom_bar(data=lineage_specific_count_per_fortnight_per_country, aes(x = date14, y=lineage_prop,fill=lineage), #fill=forcats::fct_rev(topLineages) #to reverse order
           position = "stack", stat="identity", width=14,color='black') + #date21
  geom_line(data=countryData_filtered, aes(x=date7,y= cases/highestCasesInADay, color='data'), size=1,color="black") +
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
        axis.text.y = element_text(color="black", size=12),
        plot.title = element_text(size = 10, face = "bold"),
        legend.position = "right",
        legend.text=element_text(size=13),
        legend.title=element_text(size=15),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = 15))+ 
  scale_fill_manual(values = my_colour_palette, name ="Lineage")+
  facet_wrap(vars(country), ncol = 5)+
  guides(fill=guide_legend(title="Lineage"))


pdf(paste0("./figures/Supp_Fig3.pdf"), width=12, height = 12)
another_plot
dev.off()

