#Regression Analysis
# This script does data pre-processing and regression of predictors of infection rates

# Read data
source("readData.R")
source("mutations.R")

require("pacman")
pacman::p_load(scales)

# rm(list=setdiff(ls(), c("mutationData", "continentData1", "fitness_df")))

head(continentData1)

#modify data to read as weekly value (average of weekly)
## 1. OWID data #########

owid_data <- continentData1 %>% 
  group_by(continent, date7) %>% 
  dplyr::summarise(
    smoothed_cases = mean(cases),
    smoothed_tests = mean(tests),
    smoothed_vaccine_doses = mean(vaccine_doses),
    stringency_index = mean(govt_index)
    )


## 2a. Virus Diversity ie Distinct lineages per time period (7 days)##############
#diversity
  #because some regions sequence more, distinct lineages may biased to more sequences. 
  #We divide distinct lineages by total number of sequences to achieve a #diversity_score
  # We will smooth the numbers afterwards

#total number of seqs per contiinent per day
pacman::p_load(plyr) 

cont_df1 <- mutationData %>%
  ddply(~date7+continent,summarise,
        distinct_lineages_per_continent_per_week=length(unique(lineage)))


seq_count <- dplyr::count(mutationData, continent, date7) %>%
  dplyr::rename(seq_count = n)


lineage_diversity_df<- cont_df1 %>%
  left_join (seq_count, 
             by= c("continent", "date7")) %>%
  mutate(
    viral_diversity=distinct_lineages_per_continent_per_week/seq_count,  
    #diversity taking into account total number of sequences. 
    #this doesn't look good at the beginning, diversity 'appears' high as sequences were few
          # drop data before 1st April 2020 before the linear regression
  )  %>%
  dplyr::select(continent, date7, viral_diversity) %>%
  #distinct() %>%    #drop duplicated rows
  group_by(continent) %>%
  dplyr::arrange(date7) %>%
  mutate(
    diversity_score = zoo::rollmean(viral_diversity, k = 4, fill = "extend"),  #biweekly rolling mean - for smoothing,
  ) %>%
  select(-c(viral_diversity))

pacman::p_unload(plyr) # avoid conflict with dplyr 

## 2b. Define fitness #######

### Approach 1 - use vocs (they are 'fit lineages) ####
fitLineages<-c("Alpha", "Beta", "Gamma", "Delta", "Omicron")

### Approach 2 ######
    # Get proportion of lineages per time period
    # Highest proportion (of total sequences) a lineage has ever achieved 
    # Lineages above a certain threshold of percent infections would be considered to have fitness advantage

dplyr::count(mutationData, continent, date7, lineage) -> continent_weekly_lineage_summary

# #total sequences per continent per week
data_set_cont <- continent_weekly_lineage_summary %>% 
  group_by(continent, date7) %>% 
  dplyr::summarise(
    tot_sequences = sum(n)) %>%
  left_join(continent_weekly_lineage_summary, by = c("continent", "date7")) %>%
  mutate (lineage_specific_proportion = n/tot_sequences)  %>%
  dplyr::filter(date7 >= "2020-04-01")

#top lineages - if a lineage has achieved a lineage of x in a given week, it'd be consisered fit
df_range <- 1:length(data_set_cont$lineage_specific_proportion)

fitLineages <-list() #initialize list
counter <- 1

for (i in df_range) {
  if (data_set_cont$lineage_specific_proportion[i] > 0.65){
    fitLineages[[counter]] <- data_set_cont$lineage[i]
    counter <- counter +1
    #print(data_set_cont$proportion[i])
  }
}
    
    #unique lineages
fitLineages <- unique(fitLineages) 

### Approach 1 - use vocs (they are 'fit lineages)
fitLineages<-c("Alpha", "Beta", "Gamma", "Delta", "Omicron")

fit_lineage_count <- continent_weekly_lineage_summary %>% 
  dplyr::filter(lineage %in% fitLineages) %>% 
  group_by(continent, date7) %>% 
  dplyr::summarise(
  total_fit_lineages = sum(n, na.rm = T)
  )

# #proportion of fit lineages per continent per week
fit_lineages_proportions_continent <- data_set_cont %>% 
  left_join(fit_lineage_count,
            by = c("continent", "date7")) %>% 
  dplyr::select(-c(lineage, n, lineage_specific_proportion)) %>% 
  distinct() %>% 
  mutate(
    total_fit_lineages = replace(total_fit_lineages, which(is.na(total_fit_lineages)), 0),   
    virus_fitness = (total_fit_lineages/tot_sequences)
  ) %>%
  #dplyr::select(-c(tot_sequences, total_fit_lineages)) %>%  
  dplyr::group_by(continent) %>% 
  dplyr::arrange(date7) %>%
  mutate(
    virus_fitness = zoo::rollmean(virus_fitness, k = 4, fill = "extend"),  # rolling mean - for smoothing
  ) 

### Approach 3 ######
# count 'fit' mutations and find their proportion
# This is done by mutations.R

  #outcome is fitness_df

#rename countries and continents in gisaid data to match owid data
renameFunction<-function(locationName, na.rm=T){
  gsub("_", " ", locationName) %>%    #replace underscore with space
    str_to_title(locationName)   #change country to sentence case
}


fitness_df <- fitness_df %>% 
  dplyr::mutate_at(
    c("Continent"), renameFunction
  ) %>%
  mutate(across(where(is.character),str_trim)) %>%
  drop_na()


head(lineage_diversity_df)
head(fitness_df)
head(fit_lineages_proportions_continent) # approaches 1 and 2

## 2c.  Merge mutation data - diversity and fitness ####
comprehensive_mutation_df <- fitness_df %>% 
  dplyr::select(Continent, date7, fitness_score3) %>% 
  dplyr::rename(virus_fitness = fitness_score3,
                continent = Continent) %>% 
  left_join(lineage_diversity_df, by=c("continent", "date7")) %>% 
  dplyr::filter(date7 >= "2020-03-01")


head(comprehensive_mutation_df)

##3. Regression df #################################
## ie Merge gisaid with owid data
head(owid_data)
head(comprehensive_mutation_df)

df_continent<-owid_data %>% 
  left_join(comprehensive_mutation_df, by=c("continent", "date7"))

names(df_continent)

df_continent<-mutate_at(df_continent, c("smoothed_cases", "smoothed_vaccine_doses",
                                        "stringency_index", "diversity_score", "virus_fitness"),
                        ~replace(., is.na(.), 0))
                        
rm(list=setdiff(ls(), c("df_continent","mutationData", "continentData1")))   

head(df_continent)

regressionData<- df_continent %>% 
  dplyr::group_by(continent) %>% 
  dplyr::arrange(date7) %>%
  mutate(
    smoothed_cases = zoo::rollmean(smoothed_cases, k = 2, fill = "extend"),  # rolling mean - for smoothing
    diversity_score = zoo::rollmean(diversity_score, k = 2, fill = "extend")
  ) %>%
  dplyr::ungroup() %>%
  dplyr::filter(date7 >"2020-04-01" & date7 <"2022-01-10") #data seems a bit unstable outside these dates

ggplot(data= regressionData)+
  geom_line(aes(x=date7, y=virus_fitness, color="fitness"))+
  geom_line(aes(x=date7, y=diversity_score*10, color="diversity"))+
  facet_wrap(vars(continent))+
  theme_bw()

ggplot(data= regressionData)+
  #geom_line(aes(x=date7, y=smoothed_cases, color="cases"))+
  geom_line(aes(x=date7, y=stringency_index, color="govt_index"))+
  #geom_line(aes(x=date7, y=smoothed_vaccine_doses, color="vaccines"))+
  facet_wrap(vars(continent))+
  theme_bw()

#all looks good! Time to fit a linear model 

##  Correlation test first 
## correlation #######

 # per continent
regressionData %>% 
  group_by (continent) %>% 
  dplyr::select(smoothed_cases, virus_fitness, diversity_score, 
                stringency_index, smoothed_vaccine_doses) %>% 
  dplyr::rename(`Virus fitness` = virus_fitness,
                `Infection rate` = smoothed_cases,
                Vaccination =smoothed_vaccine_doses,
                `Virus diversity` =diversity_score,
                `Stringency index` =stringency_index) %>%
  do(data.frame(Cor=t(cor(.[,2:6], .[,2], method = "pearson"))))


# corr data
corr_data<-regressionData %>% 
  ungroup () %>% 
  dplyr::select(smoothed_cases, virus_fitness, diversity_score, 
                stringency_index, smoothed_vaccine_doses) %>% 
  dplyr::rename(`Virus fitness` = virus_fitness,
                `Infection rate` = smoothed_cases,
                Vaccination = smoothed_vaccine_doses,
                `Virus diversity` =diversity_score,
                `Stringency index` =stringency_index)

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

pdf(paste0("./figures/Supplementary_Figure_1.pdf"), width=6, height = 6)
corr_plot
dev.off()



# Time to fit a linear model, for real this time :-)

# Multiple Linear Regression
    #lm([target variable] ~ [predictor variables], data = [data source])

# Use the plyr package to run regression per continent easily across groups
library(plyr)
#defn function for linear model
lm_with_plyr<-function(df){
  lm(smoothed_cases ~ stringency_index +
       smoothed_vaccine_doses +
       virus_fitness +
       diversity_score,
     data = df)
}
#divide df into smaller dfs (creates a list/panel of dfs) based on continent and run lm for each of the dfs
models<-dlply(regressionData, "continent", lm_with_plyr)

# Apply coef to each model and return a data frame
ldply(models, coef)

# Print the summary of each model
l_ply(models, summary, .print = TRUE)

# R Squared
laply(models, function(mod) summary(mod)$r.squared)

pacman::p_unload(plyr)

summary(models$Africa)

# Find fitted values

fitAndPlotFunction <- function(selected_continent){
  regressionData %>%
    dplyr::filter(continent == selected_continent) %>%
    mutate(fitted_cases = fitted(models[[selected_continent]]),
           smoothed_fitted = zoo::rollmean(fitted_cases, k = 3, fill = "extend")) %>%
    #dplyr::arrange(date) %>%
    #mutate(smoothed_fitted = zoo::rollmean(fitted_cases, k = 3, fill = "extend"))
    ggplot()+
    geom_line(aes(x=date7, y=smoothed_cases, color="Data"))+
    geom_line(aes(x=date7, y=smoothed_fitted, color="Fitted"))+
    theme_bw()+
    scale_x_date(date_labels = "%b",date_breaks = "6 months", limits = as.Date(c('2020-01-01','2022-02-10')))+ #date filter
    scale_y_continuous(labels = scales::unit_format(unit = "k", scale = 1e-3)) +
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






obermeyers<- models
sergeis<-models







