#Regression Analysis at country-level
# This script does data pre-processing and model selection of predictors of infection rates


#Read data
#source("readData.R")

### Comprehensive mutation data  ####
# data #####
head(countryData)   # owid
head(country_lineageData)   
head(country_fit_mutation_summary)

# Merging datasets
merged_data <- countryData %>%
  dplyr::rename(country = location) %>%
  left_join(country_lineageData, by = c("country", "date7")) %>%
  left_join(country_fit_mutation_summary, by = c("country" = "Country", "date7")) %>%
  dplyr::mutate(
    diversity_score = zoo::rollmean(diversity_score, k = 4, fill = 'extend'),
    fitness_score = zoo::rollmean(fitness_score, k = 4, fill = 'extend')
  ) %>%
  dplyr::filter(country %in% countries_with_data, 
                date7 >"2020-04-01" & date7 <"2022-01-10") 



ggplot(data=merged_data)+
  geom_line(aes(x=date7, y=previous_burden))+
  facet_wrap(vars(country))

# rm(list=setdiff(ls(), c("lineageData", "owidData0", "df_country", "fitness_df", "fit_mutations_count", "countryData")))



# Select data####### 
#
regressionData<- merged_data %>%
  dplyr::select(country, date7, adj_cases, previous_burden, fitness_score, diversity_score, 
                stringency_index, vaccine_doses 
                ) %>%
  drop_na() # we shouldn't be having any NAs at this point given the data processing

#all looks good! Time to fit a linear model :-)

# Correlation ####

p <- regressionData %>% 
  group_by (country) %>% 
  dplyr::select(-c(date7)) %>% 
  dplyr::rename(`Virus fitness` = fitness_score,
                `Infection rate` = adj_cases,
                Vaccination = vaccine_doses,
                `Virus diversity` = diversity_score,
                `Stringency index` = stringency_index,
                `Previous burden` = previous_burden
                ) %>%
  do(data.frame(Cor=t(cor(.[,2:6], .[,2], method = "pearson"))))

# Corr coefficients
p %>%
  flextable::flextable() %>%
  flextable::autofit() %>%
  flextable::save_as_html(path = "figures/correlation_countries.html")


# corr data (all countries)
corr_data<-regressionData %>% 
  ungroup () %>% 
  dplyr::select(adj_cases, fitness_score, diversity_score, 
                stringency_index, vaccine_doses, previous_burden) %>% 
  dplyr::rename(`Virus fitness` = fitness_score,
                `Infection rate` = adj_cases,
                Vaccination =vaccine_doses,
                `Virus diversity` =diversity_score,
                `Stringency index` =stringency_index,
                `Previous burden` =previous_burden)

# p-values
p.mat <- cor_pmat(corr_data)

# Multiple Linear Regression
#lm([target variable] ~ [predictor variables], data = [data source])

# Use the plyr package to run regressions per country instead of filtering each time
library(plyr)
#defn function for linear model
lm_with_plyr<-function(df){
  lm(adj_cases ~ stringency_index +
       vaccine_doses +
       fitness_score +
       diversity_score+
       previous_burden,
     data = df)
}
#divide df into smaller dfs (creates a list of dfs) based on country and run lm for each of the dfs
models<-dlply(regressionData, "country", lm_with_plyr)

# Apply coef to each model and return a data frame
coefficients_df <- ldply(models, coef)
coefficients_df$R_squared <- laply(models, function(mod) summary(mod)$r.squared)

format(coefficients_df, digits = 2) %>%
  flextable::flextable() %>%
  flextable::autofit() %>%
  flextable::save_as_html(path = "figures/regression_coefficients.html")

# Print the summary of each model
l_ply(models, summary, .print = TRUE)
# R Squared
p<-laply(models, function(mod) summary(mod)$r.squared)
p
data.table(country = names(p), value = p)


pacman::p_unload(plyr)
# Find fitted values
summary(models$France)


# Plot Q-Q plots
# Step 1: Extract residuals
resid_data <- lapply(names(models), function(country) {
  data.frame(
    Country = country,
    Residuals = resid(models[[country]])
  )
})

# Step 2: Combine into single dataframe
resid_data_combined <- do.call(rbind, resid_data)

# Step 3: Create QQ plots using ggplot2
qq_plot<- ggplot(resid_data_combined, aes(sample = Residuals)) +
  geom_qq() +
  facet_wrap(~Country, scales = "free") +
  labs(
    x = "Theoretical Quantiles",
    y = "Sample Quantiles")+
  theme_bw()

pdf(paste0("./figures/FigureA3.pdf"), width=12, height = 10)
qq_plot
dev.off()




fitModelFunction <- function(selected_country){
  data_for_country <- regressionData %>%
    dplyr::filter(country == selected_country)
  
  # Get predictions and confidence intervals
  preds <- predict(models[[selected_country]], newdata = data_for_country, 
                   interval = "confidence", level = 0.95)
  
  # Merge predictions and confidence intervals with original data
  data_for_country %>%
    dplyr::mutate(
      fitted_cases = preds[, "fit"],
      lower_conf_int = preds[, "lwr"],
      upper_conf_int = preds[, "upr"]
    )
}

countryList <- names(p)

# run function & combine
fittedCountries <- lapply(countryList, fitModelFunction)
fittedCountries <- bind_rows(fittedCountries)

# plot
plot2<-ggplot(data = fittedCountries)+
  geom_ribbon(aes(x = date7, ymin = lower_conf_int, ymax = upper_conf_int), fill = "grey70", alpha = 0.5) +
  geom_line(aes(x=date7, y=adj_cases, color="Data"))+
  geom_line(aes(x=date7, y=fitted_cases, color="Fitted"))+
  theme_bw()+
  ylab('Adjusted daily new cases per million') +
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
  facet_wrap(~country, scales = "free_y")

plot2

# pdf(paste0("./figures/FigureA2.pdf"), width=12, height = 10)
# plot2
# dev.off()


