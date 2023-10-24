#Regression Analysis
# This script does data pre-processing and regression of predictors of infection rates

# Read data
#source("./scripts/readData.R")

require("pacman")
pacman::p_load(scales)

# rm(list=setdiff(ls(), c("lineageData", "continentData1", "fitness_df")))


## 1. Datasets ##########
head(continentData1)        # owid
head(continent_fit_mutation_summary) # fit mutations
head(cont_lineageData)        # lineage diversity



# check
ggplot(data = cont_lineageData) +
  geom_line(aes(x=date, y=diversity_score))+
  facet_wrap(vars(continent))


##3. Regression df #################################
# Combine data

combined_data <- continent_fit_mutation_summary %>%
  dplyr::rename(continent = Continent)  %>%
  left_join(continentData1, by = c("continent", "date")) %>%
  left_join(cont_lineageData, by = c("continent", "date")) %>%
  dplyr::select(-c(Mutation_Count, fitGenomes, Total_Genomes, Total_Mutations)) %>%
  dplyr::filter(date >"2020-04-01" & date <"2022-01-10") #data seems a bit unstable outside these dates


head(combined_data)

# plot to check
ggplot(data= combined_data)+
  geom_line(aes(x=date, y=cases, color="cases"))+
  # geom_line(aes(x=date, y=govt_index, color="govt_index"))+
  # geom_line(aes(x=date, y=vaccine_doses, color="vaccines"))+
  # geom_line(aes(x=date, y=previous_burden, color="previous_burden"))+
  facet_wrap(vars(continent))+
  theme_bw()

                      

#all looks good! Time to fit a linear model 

##  Correlation test first 
## correlation #######

 # per continent
combined_data %>% 
  group_by (continent) %>% 
  dplyr::select(cases, fitness_score, diversity_score, 
                govt_index, vaccine_doses, previous_burden) %>% 
  dplyr::rename(`Virus fitness` = fitness_score,
                `Infection rate` = cases,
                Vaccination =vaccine_doses,
                `Virus diversity` =diversity_score,
                `Stringency index` =govt_index,
                `Previous burden` =previous_burden) %>%
  do(data.frame(Cor=t(cor(.[,2:7], .[,2], method = "pearson"))))


# corr data
corr_data<-combined_data %>% 
  ungroup () %>% 
  dplyr::select(cases, fitness_score, diversity_score, 
                govt_index, vaccine_doses, previous_burden) %>% 
  dplyr::rename(`Virus fitness` = fitness_score,
                 `Infection rate` = cases,
                 Vaccination =vaccine_doses,
                 `Virus diversity` =diversity_score,
                 `Stringency index` =govt_index,
                 `Previous burden` =previous_burden)

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

corr_plot


# Time to fit a linear model

# Use the plyr package to run regression per continent easily across groups
library(plyr)
#defn function for linear model
lm_with_plyr<-function(df){
  lm(cases ~ govt_index +
       vaccine_doses +
       previous_burden +
       fitness_score +
       diversity_score,
     data = df)
}
#divide df into smaller dfs (creates a list/panel of dfs) based on continent and run lm for each of the dfs
models_cont<-dlply(combined_data, "continent", lm_with_plyr)

# Apply coef to each model and return a data frame
ldply(models_cont, coef)

# Print the summary of each model
l_ply(models_cont, summary, .print = TRUE)

# R Squared
laply(models_cont, function(mod) summary(mod)$r.squared)

# data to plot residuals
forest_plot_data <- ldply(models_cont, function(model) {
  coef_summary <- summary(model)$coefficients
  
  data.frame(
    variable = rownames(coef_summary),
    coef = coef_summary[, "Estimate"],
    se = coef_summary[, "Std. Error"],
    pvalue = coef_summary[, "Pr(>|t|)"] # or "Pr(>|z|)" depending on model type
  )
}) %>%
  dplyr::filter(!variable =="(Intercept)")  %>%
  dplyr::rename(trait=continent,
                name=variable,
                beta=coef)

# R-squared
laply(models_cont, function(mod) summary(mod)$r.squared)
# residual std error
laply(models_cont, function(mod) summary(mod)$sigma)

laply(models_cont, function(mod) summary(mod)$df[2])

# avoid conflict with dplyr
pacman::p_unload(plyr)


# Forest plot########
forest_plot_data <- forest_plot_data %>%
  mutate(name = recode(name, 
                       govt_index = "Government Index",
                       vaccine_doses = "Vaccine Doses",
                       previous_burden = "Previous Burden",
                       fitness_score = "Virus fitness",
                       diversity_score = "Virus diversity"))

ggforestplot::forestplot(
  df = forest_plot_data,
  estimate = beta,
  pvalue = pvalue,
  psignif = 0.005,
  xlab = "Regression coefficients",
  colour = trait
)



summary(models_cont$Africa)


library(purrr)

# Find fitted values

fitAndPlotFunction <- function(models_to_fit=models_cont){
  
  # fit
  all_continents_combined <- map_dfr(names(models_cont), ~ {
    data_for_continent <- combined_data %>% dplyr::filter(continent == .x)

    # Get fitted values and standard errors for predictions
    predictions <- predict(models_cont[[.x]], newdata = data_for_continent, interval = "confidence", level = 0.95)

    # # Calculate 95% confidence intervals
    # ci_lower <- predictions$fit - 1.96 * predictions$se.fit
    # ci_upper <- predictions$fit + 1.96 * predictions$se.fit

    data_for_continent %>%
      mutate(fitted_cases = predictions[,"fit"],
             ci_lower = predictions[, "lwr"],
             ci_upper = predictions[, "upr"],
             #fitted_cases = fitted_cases,
             smoothed_fitted = zoo::rollmean(fitted_cases, k = 14, fill = "extend")
             ) %>%
      dplyr::select(continent, date, cases, smoothed_fitted, ci_lower, ci_upper)
  })
  
  
  # plot
  ggplot(all_continents_combined) +
    # Ribbon for confidence intervals
    geom_ribbon(aes(x = date, ymin = ci_lower, ymax = ci_upper), fill = "grey70", alpha = 0.4) +
    
    # Original lines for data and fitted values
    geom_line(aes(x = date, y = cases, color = "Data")) +
    geom_line(aes(x = date, y = smoothed_fitted, color = "Fitted")) +
    
    theme_bw() +
    scale_x_date(date_labels = "%b \n %Y", date_breaks = "6 months", limits = as.Date(c('2020-01-01','2022-02-10'))) +
    scale_y_continuous(labels = scales::unit_format(unit = "k", scale = 1e-3)) +
    theme(axis.text.x = element_text(angle=0, color="black", size=13),
          axis.title.x = element_blank(),
          axis.title.y = element_text(color="black", size=15, face="bold"),
          axis.text.y = element_text(color="black", size=14),
          legend.position = "top",
          legend.text = element_text(size = 13),
          strip.background = element_rect(fill = "white"),
          strip.text = element_text(size = 15)
          ) + 
    guides(color = guide_legend(title = NULL)) +
    scale_color_manual(values = c("black", "violetred1")) +
    labs(y = "Adjusted daily new cases per million") +
    #theme_bw() +
    facet_wrap(~continent, scales = "free_y",ncol = 2)
  }



fig2A<- fitAndPlotFunction(models_cont)
fig2A


fig2A <- fig2A + theme(plot.tag = element_text(size = 20))
corr_plot <- corr_plot + theme(plot.tag = element_text(size = 20))


fig2<- corr_plot + fig2A + plot_layout(widths = c(3/10, 7/10)) + plot_annotation(tag_levels = 'a', tag_prefix = "(", tag_suffix = ")")
fig2

# pdf(paste0("./figures/Figure2.pdf"), width=11, height = 8)
# fig2
# dev.off()


# Coefficients table ###########
# Function to extract desired statistics for each model
extract_stats <- function(model, continent_name){
  s <- summary(model)
  data.frame(
    Estimate = coef(model),
    Std.Error = s$coefficients[, "Std. Error"],
    t.value = s$coefficients[, "t value"],
    Pr = s$coefficients[, "Pr(>|t|)"],
    continent = continent_name
  )
}

# Extract coefficients, standard errors, t-values, and p-values and transform to a data frame
coefficients_df <- imap_dfr(models_cont, ~ extract_stats(.x, .y), .id = "variable") %>%
  dplyr::mutate(coef = rownames(.)) 

rownames(coefficients_df) <- NULL
coefficients_df$coef <- gsub("\\(|\\)|\\.\\.\\.\\d+", "", coefficients_df$coef)



coefficients_df2<-coefficients_df %>%
  gather(key = "metric", value = "value", -continent, -coef) %>%
  unite("key", coef, metric, sep = " ") %>%
  spread(key = "key", value = "value") %>%
  as.data.frame()

coefficients_df2<-t(coefficients_df2)

###write.csv(coefficients_df2, "./figures/regression_coefs.csv")





########

# confidence intervals 
coefficients_df

# Function to compute confidence intervals
compute_CI <- function(model, data, conf_level = 0.95) {
  alpha <- 1 - conf_level
  
  # Degrees of freedom
  df <- df.residual(model)
  
  # t-value for desired confidence level
  t_value <- qt(1 - alpha/2, df)
  
  # Compute predicted values
  predicted <- predict(model, data, se.fit = TRUE)
  
  # Confidence intervals
  CI_lower <- predicted$fit - t_value * predicted$se.fit
  CI_upper <- predicted$fit + t_value * predicted$se.fit
  
  data.frame(Fitted = predicted$fit, CI_lower = CI_lower, CI_upper = CI_upper)
} 
results_list <- lapply(models_cont, compute_CI, data = coefficients_df2)
names(results_list) <- names(models_cont)

# Check results for a continent
results_list$Africa











