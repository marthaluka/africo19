
# This script reads data and should be run before:
    # lineage_infectionRate_plot.R
    # regression_continent.R
    # regression_country.R

#Read data

#rm(list = ls())

## packages
require("pacman")
pacman::p_load(data.table, tidyverse, zoo, RColorBrewer, 
               directlabels, corrplot, data.table, stringr, 
               spData, ggrepel, ggcorrplot, patchwork, corrr,
               cowplot, egg, here, visreg, broom,
               vroom, stringi, parallel, lubridate,
               slider
               )

# OWID data ######### 

## read data
owidData1<- vroom::vroom(
  "./data/owid-covid-data.csv"
  #"https://raw.githubusercontent.com/owid/covid-19-data/master/public/data/owid-covid-data.csv",
  #stringsAsFactors = FALSE   #for fread
  ) %>%                             
  dplyr::select(continent, location, date, new_cases_per_million,total_vaccinations_per_hundred,
                stringency_index, new_tests_per_thousand, tests_per_case, 
                total_cases_per_million
                ) %>%
  drop_na(. ,c(location, continent) # drop with misiing loaction data
          ) %>%
  dplyr::mutate(date=as.Date(date),  # date formats
                across(where(is.character),str_trim) # strip extra spaces
  )      
        # People vaccinated (to handle boosters), hospital/ icu admissions etc are important variables 
               # too but we simply lack suffient data in lots of countries 


## find the sum of NAs per column
find_NAs <- function(df) {
  apply(df, 2, function(x) sum(is.na(x)))
}

find_NAs(owidData1)

#Pre-processing
    #biweekly rolling mean - for smoothing
owidData1 <- owidData1 %>%
  dplyr::filter(new_cases_per_million>=0) %>%    #drop any negative values 
  mutate_at(c("new_cases_per_million","new_tests_per_thousand",
              "total_vaccinations_per_hundred", "stringency_index"), ~replace(., is.na(.), 0)) %>%
  dplyr::group_by(location) %>% 
  dplyr::arrange(date) %>% 
  dplyr::mutate(smoothed_cases_per_million = zoo::rollmean(new_cases_per_million, k = 14, fill = 'extend'),
                smoothed_vaccines_per_hundred = zoo::rollmean(total_vaccinations_per_hundred, k = 14, fill = 'extend')) %>%  
  dplyr::ungroup() 
  # dplyr::select(-c(new_cases_per_million, total_vaccinations_per_hundred))


## hdi
hdi <- read.csv("./data/hdi_stats_2021.csv")%>%
  dplyr::select(-c(`HDI.rank`)) %>%
  dplyr::rename(hdi = X2021_HDI)

## continents#######
# Join the datasets to get HDI values

continentData1 <- owidData1 %>% 
  inner_join(., hdi, by = c("location" = "Country"), relationship = "many-to-many") %>%
  dplyr::mutate(
    # adjust cases with hdi (lower reported cases in countries with low hdi)
    reported_cases = new_cases_per_million,
    new_cases_per_million = exp(log(new_cases_per_million + 1e-10) + (0.64 * (1-hdi))) # (cases + 1e-10) just avoids NaNs when finding log
  )  %>%
  dplyr::group_by(continent, date) %>% 
  dplyr::summarise(
    cases = sum(new_cases_per_million),# / length(unique(location)),
    #tests = sum(new_tests_per_thousand), #/ length(unique(location)),
    vaccine_doses = sum(total_vaccinations_per_hundred),
    reported_cases = sum(reported_cases),
    govt_index = ifelse(cases != 0, sum((new_cases_per_million/ sum(new_cases_per_million))* stringency_index), 0) # adjust govt_index to no of cases a country contributes
  ) %>% 
  dplyr::mutate(govt_index = zoo::rollmean(govt_index, k = 14, fill = 'extend'),
                vaccine_doses = zoo::rollmean(vaccine_doses, k = 14, fill = 'extend'),
                reported_cases = zoo::rollmean(reported_cases, k = 14, fill = 'extend'),
                cases = zoo::rollmean(cases, k = 14, fill = 'extend')
                
  ) %>%  
  dplyr::ungroup() %>% 
  # find previous burden by summing up previous cases in last 97 days but excluding the most recent 7 days
  dplyr::group_by(continent) %>% 
  dplyr::arrange(date) %>% 
  dplyr::mutate(previous_burden = slider::slide_dbl(cases, sum, .before = 90, .after = -1, .complete = FALSE),
                # positivity starts to decline by week 3 and subsequently becomes undetectable (https://jamanetwork.com/journals/jama/fullarticle/2765837)
                vaccine_doses = slider::slide_dbl(vaccine_doses, sum, .before = 101, .after = -12, .complete = FALSE) 
                #vaccine (Pfizer) started to offer partial protection ~12 days after the first dose (https://www.nejm.org/doi/full/10.1056/NEJMoa2034577)
  )

head(continentData1)



ggplot(data=continentData1)+
  geom_line(aes(x=date, y=reported_cases, color='reported_cases'))+
  geom_line(aes(x=date, y=cases, color='adj_cases'))+
  facet_wrap(~continent, scales = "free_y")

# continentData1 <- owidData1 %>% 
#   dplyr::group_by(continent, date) %>% 
#   dplyr::summarise(
#     cases = sum(new_cases_per_million),
#     #tests = sum(new_tests_per_thousand)/ length(unique(location)),
#     tests = sum(new_tests_per_thousand),
#     vaccine_doses = sum(total_vaccinations_per_hundred),
#     #each country's contribution of govt index is relative to no of cases
#     govt_index = ifelse(cases != 0, sum((new_cases_per_million/ sum(new_cases_per_million))* stringency_index), 0)
#   ) %>% 
#   dplyr::mutate(govt_index = zoo::rollmean(govt_index, k = 14, fill = 'extend'),
#                 vaccine_doses = zoo::rollmean(vaccine_doses, k = 14, fill = 'extend'),
#                 cases = zoo::rollmean(cases, k = 14, fill = 'extend')
#                 ) %>%  
#   dplyr::ungroup() %>% 
#   # find previous burden by summing up previous cases in last 97 days but excluding the most recent 7 days
#   dplyr::group_by(continent) %>% 
#   dplyr::arrange(date) %>% 
#   dplyr::mutate(previous_burden = slider::slide_dbl(cases, sum, .before = 90, .after = -1, .complete = FALSE),
#                     # positivity starts to decline by week 3 and subsequently becomes undetectable (https://jamanetwork.com/journals/jama/fullarticle/2765837)
#                 vaccine_doses = slider::slide_dbl(vaccine_doses, sum, .before = 101, .after = -12, .complete = FALSE) 
#                     #vaccine (Pfizer) started to offer partial protection ~12 days after the first dose (https://www.nejm.org/doi/full/10.1056/NEJMoa2034577)
#   )


## Country level ######
# countries with less than 180 days of missing sequence data (before trimming for extreme dates)
countries_with_data <- c("Argentina", "Australia", "Austria", "Belgium", "Brazil", 
                         "Canada", "Chile", "France", "Germany", "India", "Indonesia", "Ireland", "Italy", 
                         "Japan", "Kenya", "Luxembourg", "Mexico", "Netherlands", "Norway", "Russia", 
                         "Singapore", "Slovenia", "South Africa", "South Korea", "Spain", "Sweden", 
                         "Switzerland", "United Kingdom", "United States")

# Argentina, Chile, Ireland, United States
#UK fitness score
# check gaps in sequences-- rollmean. nas to zero first


countryData <- owidData1 %>% 
  dplyr::filter(location %in% countries_with_data)  %>% 
  inner_join(., hdi, by = c("location" = "Country"), relationship = "many-to-many") %>%
  dplyr::mutate(
    # adjust cases with hdi (lower reported cases in countries with low hdi)
    adj_smoothed_cases_per_million = exp(log(smoothed_cases_per_million + 1e-10) + (0.64 * (1-hdi))) # (cases + 1e-10) just avoids NaNs when finding log
  ) %>%
  # find previous burden by summing up previous cases in last 90 days
  dplyr::group_by(location) %>% 
  dplyr::arrange(date) %>% 
  dplyr::mutate(#previous_burden = (slider::slide_dbl(smoothed_cases_per_million, sum, .before = 89, .complete = FALSE))
    previous_burden = slider::slide_dbl(adj_smoothed_cases_per_million, sum, .before = 90, .after = -1, .complete = FALSE),
    # positivity starts to decline by week 3 and subsequently becomes undetectable (https://jamanetwork.com/journals/jama/fullarticle/2765837)
    vaccine_doses = slider::slide_dbl(smoothed_vaccines_per_hundred, sum, .before = 101, .after = -12, .complete = FALSE), 
    #vaccine (Pfizer) started to offer partial protection ~12 days after the first dose (https://www.nejm.org/doi/full/10.1056/NEJMoa2034577)
    cases = smoothed_cases_per_million,
    adj_cases = adj_smoothed_cases_per_million) %>%
  dplyr::select(-c(continent, new_tests_per_thousand, total_cases_per_million, 
                   new_cases_per_million, total_vaccinations_per_hundred,
                   tests_per_case, smoothed_vaccines_per_hundred, smoothed_cases_per_million, 
                   adj_smoothed_cases_per_million, hdi)
  ) %>% 
  dplyr::mutate(
    date7=as.Date(cut(date, breaks = "1 weeks", start.on.monday = TRUE)),
  ) 


# Convert tibble to data.table
countryData<- setDT(countryData)

# Compute desired statistics (weekly)
countryData <- countryData[, .(
  stringency_index = mean(stringency_index, na.rm = TRUE),
  previous_burden = mean(previous_burden, na.rm = TRUE),
  vaccine_doses = mean(vaccine_doses, na.rm = TRUE),
  adj_cases = mean(adj_cases, na.rm = TRUE)
), by = .(location, date7)]

head(countryData)



# GISAID / Lineages data##########

# Read the data
lineageData <- vroom::vroom("./data/metadata_20220128.csv", col_names = F, 
                             col_select = c('X2', 'X10', 'X5', 'X8', 'X11')) %>%
  setNames(c("gisaid_ID", "country", "lineage", "date", "continent")) %>%
  dplyr::mutate(
    date = as.Date(date),
    date14 = as.Date(cut(date, breaks = "2 weeks", start.on.monday = TRUE))
  ) %>%
  drop_na(lineage)  #drop un-assigned lineages


dim(lineageData)

lineageData <- setDT(lineageData)         # Convert data.frame to data.table

# clean up country names
renameFunction <- function(locationName) {
    replacements <- c(
      "_" = " ",
      "And " = "and ",
      "Of " = "of ",
      " The" = " the",
      "Usa" = "United States",
      "Gambia" = "The Gambia"
    )
    
    for (pattern in names(replacements)) {
      locationName <- stringi::stri_replace_all_regex(locationName, pattern, replacements[pattern], 
                                                      opts_regex = stringi::stri_opts_regex(case_insensitive = TRUE))
    }
    
    stringi::stri_trans_totitle(locationName)
  }

# run function
lineageData <- lineageData %>%
  mutate(across(c(continent, country), renameFunction)) %>%
  mutate(across(where(is.character), stringi::stri_trim_both)) %>%
  drop_na()

# Use WHO names for lineage
lineageData <- lineageData %>%
  mutate(
    country = ifelse(startsWith(country, "Uk-"), "United Kingdom", country),
    lineage = case_when(
      startsWith(lineage, "AY.") | startsWith(lineage, "B.1.617.2") ~ "Delta",
      startsWith(lineage, "B.1.1.7") | startsWith(lineage, "Q.") ~ "Alpha",
      startsWith(lineage, "B.1.351") ~ "Beta",
      startsWith(lineage, "P.1") ~ "Gamma",
      startsWith(lineage, "B.1.1.529") | startsWith(lineage, "BA.") ~ "Omicron",
      TRUE ~ lineage
    )
  )

## continent level (regression data) #########
cont_lineageData <- lineageData%>%
  # Count sequences per continent per date
  count(continent, date) %>%
  rename(tot_genomes_per_day = n) %>%
  
  # Calculate distinct lineages per continent per day
  left_join(
    lineageData %>%
      group_by(continent, date) %>%
      summarise(distinct_lineages_per_continent_per_day = n_distinct(lineage)) %>%
      ungroup(),
    by= c("continent", "date")
  ) %>%
  
  # Compute the diversity score
  mutate(
    diversity_score = distinct_lineages_per_continent_per_day / tot_genomes_per_day
  ) %>%
  
  # Select necessary columns and order by date
  select(continent, date, diversity_score) %>%
  distinct() %>%
  group_by(continent) %>%
  arrange(date) %>%
  
  # Smooth the diversity score using rolling mean
  mutate(
    diversity_score = zoo::rollmean(diversity_score, k = 14, fill = "extend")
  ) %>%
  ungroup()




## country level (regression data) ################

country_lineageData <- lineageData%>%
  dplyr::filter(country %in% countries_with_data) %>%
  # Count sequences per country per date
  count(country, date) %>%
  rename(tot_genomes_per_day = n) %>%
  
  # Calculate distinct lineages per country per day
  left_join(
    lineageData %>%
      group_by(country, date) %>%
      summarise(distinct_lineages_per_country_per_day = n_distinct(lineage)) %>%
      ungroup(),
    by= c("country", "date")
  ) %>%
  
  # Compute the diversity score
  mutate(
    diversity_score = distinct_lineages_per_country_per_day / tot_genomes_per_day
  ) %>%
  
  # Select necessary columns and order by date
  select(country, date, diversity_score) %>%
  distinct() %>%
  group_by(country) %>%
  arrange(date) %>%
  
  # Smooth the diversity score using rolling mean
  dplyr::mutate(
    diversity_score = zoo::rollmean(diversity_score, k = 14, fill = "extend"),
    date7=as.Date(cut(date, breaks = "1 weeks", start.on.monday = TRUE))
  ) %>%
  ungroup()


country_lineageData


# Convert tibble to data.table
country_lineageData<- setDT(country_lineageData)

# Compute desired statistics (weekly)
country_lineageData <- country_lineageData[, .(
  diversity_score = mean(diversity_score, na.rm = TRUE)
), by = .(country, date7)]

head(country_lineageData)


# MUTATIONS data ########

  # Mutations from www.medrxiv.org/content/10.1101/2021.09.07.21263228v2.full
mutations_increasing_fitness <- str_to_upper(c("S:H655Y", "S:T95I", "ORF1a:P3395H", "S:N764K", "ORF1a:K856R", "S:S371L",
                                               "E:T9I", "S:Q954H", "ORF9b:P10S", "S:L981F", "N:P13L", "S:G339D", "S:S375F",
                                               "S:S477N", "S:N679L", "S:S373P", "M:Q19E", "S:D796Y", "S:N969K", "S:T547K",
                                               "ORF1b:I1566V", "M:D3G", "S:G446S", "S:N440K", "M:A63T", "S:N856K", "ORF1a:A2710T",
                                               "ORF1a:I3758V", "S:E484A", "S:A67V", "S:K417N", "S:Q493R", "S:N501Y","S:Y505H", 
                                               "S:L452R", "S:P681H", "S:Q498R","S:G496S", "ORF1a:T3255I", "ORF14:G50W", "S:P681R",
                                               "N:R203M", "ORF1b:P1000L", "ORF1a:P2287S", "M:I82T", "ORF3a:S26L", "N:D63G", "N:G215C",
                                               "ORF1a:V3718A", "ORF9b:T60A"))

# Read data 
# Count Mutations by continent and date
clean_mutations <- fread("./data/mutations_clean_file2.csv")[, .(Gisaid, Date, Continent, Country, NewMutation)]
clean_mutations <- clean_mutations[!is.na(Date)]

## Continent fitness score #####

# 1. Calculate Genomes and NewMutation totals for each Continent and day grouping.
totals <- clean_mutations[
  , date := as.Date(Date)][
    , .(Total_Genomes = uniqueN(Gisaid), Total_Mutations = .N), by = .(Continent, date)]

# 2. Summarize the filtered data for mutations_increasing_fitness
summary_data <- clean_mutations[
  NewMutation %in% mutations_increasing_fitness
][, .(Mutation_Count = .N, fitGenomes = uniqueN(Gisaid)), by = .(Continent, date)]

# 3. Calculate fitness_score using a join.
continent_fit_mutation_summary <- summary_data[
  totals, on = .(Continent, date)
][, fitness_score := Mutation_Count/(Total_Genomes + log(Total_Mutations))
][order(Continent, date)]

# Replace NAs in fitness_score with zero
continent_fit_mutation_summary <- continent_fit_mutation_summary[, fitness_score := replace(fitness_score, is.na(fitness_score), 0)
][, fitness_score := rollmean(fitness_score, k = 14, fill = 'extend', align = "right"), by = .(Continent)]
                              

# rename
continent_fit_mutation_summary<- continent_fit_mutation_summary %>%
  dplyr::mutate(across(c(Continent), renameFunction)) %>%
  mutate(across(where(is.character), stringi::stri_trim_both)) %>%
  drop_na()

                                                                                                                                                                

rm(list = c("totals", "summary_data"))


## Country fitness score #####



# # 1. Calculate Genomes and NewMutation totals for each Country and date grouping.
# totals <- clean_mutations[
#   , date := as.Date(Date)][
#     , .(Total_Genomes = uniqueN(Gisaid), Total_Mutations = .N), by = .(Country, date)]
# 
# # 2. Summarize the filtered data for mutations_increasing_fitness
# summary_data <- clean_mutations[
#   NewMutation %in% mutations_increasing_fitness
# ][, .(Mutation_Count = .N, fitGenomes = uniqueN(Gisaid)), by = .(Country, date)]
# 
# 
# # 3. Calculate fitness_score using a join.
# country_fit_mutation_summary <- summary_data[
#   totals, on = .(Country, date)
# ][, fitness_score := Mutation_Count/(Total_Genomes + log(Total_Mutations))
# ][order(Country, date)]
# 
# # Replace NAs in fitness_score with zero
# country_fit_mutation_summary <- country_fit_mutation_summary[, fitness_score := replace(fitness_score, is.na(fitness_score), 0)]
# 
# 
# country_fit_mutation_summary <- country_fit_mutation_summary %>%
#   mutate(across(c(Country), renameFunction)) %>%
#   mutate(across(where(is.character), stringi::stri_trim_both)) %>%
#   drop_na() %>% 
#   dplyr::filter(Country %in% countries_with_data)
# 
# # Compute the 14-day rolling average for each Continent group
# country_fit_mutation_summary <-
#   country_fit_mutation_summary[, fitness_score := rollmean(fitness_score, k = 14, 
#                                                                     fill = 'extend', align = "right"), by = .(Country)]
# 
# 
# rm(list = c("totals", "summary_data"))


country_fit_mutation_summary <- read.csv("./data/country_fit_mutation_summary.csv", header=T) %>%
  dplyr::mutate(date7= as.Date(date7))

# check
ggplot(data = country_fit_mutation_summary) +
  geom_line(aes(x=date7, y=fitness_score))+
  facet_wrap(vars(Country))



