
# Part 1: Prepare the file #####
# Load libraries
require("pacman")
pacman::p_load(data.table, tidyverse, vroom, parallel)

# create textConnection (file doesn't meet requirements for vroom or fread)
input <- file("../data/pseudoref_aa_jan2022_cleaned.csv", "r", blocking = FALSE) 

# Read data
m_data <- readLines(input)

# delete header row
m_data = m_data[-(1)]

# Create a clean data table with just what we need
clean_df = data.table(
  Gisaid=character(), 
  Date=character(), 
  Country=character(),
  Continent=character(),  
  Mutation=character()
)

# Function to extract data
counter = 0
process_line <- function(i) {
  split_row = strsplit(i, ",") 
  if (! length(split_row[[1]]) == 17) {
    counter=counter+1
  }
  tmp = paste0(sapply(split_row, tail,1), sep=":",split_row[[1]][[1]], split_row[[1]][[4]], split_row[[1]][[2]])
  row_w_detailed_mutation = c(split_row, tmp)
  
  extract_continent = str_split(row_w_detailed_mutation[[1]][[16]],"/" , simplify = T)
  
  of_importance = data.table(Gisaid = row_w_detailed_mutation[[1]][[7]], Date = row_w_detailed_mutation[[1]][[13]], 
                             Country = row_w_detailed_mutation[[1]][[15]],
                             Continent= extract_continent[[1]][[1]], Mutation = row_w_detailed_mutation[[2]][[1]])
}


# parallelize function using mclapply
clean_df2 <- mclapply(m_data, process_line, mc.cores = 14) # Running on 14 cores

# write this to file to release some memory + be safe should R crash (append is faster than binding rows first)
fwrite(transpose(clean_df2), "../data/mutations_clean_file.csv", col.names=F, sep=",", append=TRUE)

rm(list=ls())

# Read file to memory
clean_mutations<-vroom::vroom("../data/mutations_clean_file.csv")
names(clean_mutations) <- c("Gisaid", "Date", "Country", "Continent", "Mutation")
clean_mutations$Date <- as.Date(clean_mutations$Date)
clean_mutations$Month<-as.Date(cut(clean_mutations$Date,breaks = "1 month"))
clean_mutations$Mutation <- str_to_upper(clean_mutations$Mutation)
clean_mutations<- clean_mutations %>%
  dplyr::filter(Continent %in% c("AFRICA", "ASIA", "EUROPE", "SOUTH_AMERICA", "OCEANIA", "NORTH_AMERICA"))

# Correct residue numbering which is currently continuous as opposed to per ORF 

library(splitstackshape)
#clean_mutations<-data.table(clean_mutations)
clean_mutations<-cSplit(clean_mutations, splitCols = "Mutation", sep = ":", direction = "wide", drop = FALSE) %>% dplyr::select(-c(Mutation))

clean_mutations <- clean_mutations %>% 
  separate(Mutation_2, into = c('Pos_tmp', 'ALT'), sep = -1, convert = TRUE) %>% 
  separate(Pos_tmp, into = c('REF', 'Pos'), sep = 1, convert = TRUE)

clean_mutations <- clean_mutations %>% mutate(New_position = case_when(
  Mutation_1 == "ORF1AB" ~ Pos - 0,
  Mutation_1 == "S" ~ Pos - 7096,
  Mutation_1 == "ORF3A" ~ Pos - 8370,
  Mutation_1 == "E" ~ Pos - 8646,
  Mutation_1 == "M" ~ Pos - 8722,
  Mutation_1 == "ORF6" ~ Pos - 8945,
  Mutation_1 == "ORF7A" ~ Pos - 9007,
  Mutation_1 == "ORF8" ~ Pos - 9129,
  Mutation_1 == "N" ~ Pos - 9251,
  Mutation_1 == "ORF10" ~ Pos - 9671
))

clean_mutations$NewMutation <- str_c(clean_mutations$Mutation_1,":", clean_mutations$REF, 
                                     clean_mutations$New_position,
                                     clean_mutations$ALT)
head(clean_mutations)

#Split orf1ab to distinct a and b
FILTERED <- clean_mutations %>%dplyr::filter(Mutation_1 == "ORF1AB")

FILTERED <- FILTERED %>% mutate(
  New_position = case_when(
    Pos <= 4401 ~ Pos - 0,
    Pos > 4401 ~ Pos - 4401),
  Mutation_1 = case_when(
    Pos <= 4401 ~ "ORF1A",
    Pos > 4401 ~ "ORF1B")
)

FILTERED$NewMutation <- str_c(FILTERED$Mutation_1,":", FILTERED$REF, FILTERED$New_position,
                              FILTERED$ALT)
head(FILTERED)


clean_mutations <- clean_mutations %>% 
  dplyr::filter(! Mutation_1 == "ORF1AB")

fwrite(clean_mutations, "../data/mutations_no_orf1.csv")
fwrite(FILTERED, "../data/mutations_just_orf1.csv")

       # merge the two using cat into one file (../data/mutations_clean_file2.csv)



# We start here #####

require("pacman")
pacman::p_load(data.table, tidyverse, vroom, parallel)

clean_mutations<-vroom::vroom("../data/mutations_clean_file2.csv") %>%
  select(Gisaid, Date, Continent, Country, NewMutation)

# tibble to data.table for ? faster wrangling
# clean_mutations<-data.table(clean_mutations)

# Date
clean_mutations$Date <- as.Date(clean_mutations$Date)
clean_mutations$date7<- as.Date(cut(clean_mutations$Date,breaks = "1 week",
                                 start.on.monday = TRUE))

#total genomes per continent per week
tmp <- clean_mutations %>% 
  dplyr::select(Gisaid, Continent, date7) %>% 
  distinct()

dplyr::count(tmp, Continent, date7) -> genome_count_per_cont_per_week

#total mutations per continent week
dplyr::count(clean_mutations, Continent, date7) -> continent_weekly_mutations

#specific mutations count per cont per week
dplyr::count(clean_mutations, Continent, date7, NewMutation) -> continent_weekly_specific_mutation_summary


# Mutations from www.medrxiv.org/content/10.1101/2021.09.07.21263228v2.full
mutations_increasing_fitness <- str_to_upper(c("S:H655Y", "S:T95I", "ORF1a:P3395H", "S:N764K", "ORF1a:K856R", "S:S371L",
                                               "E:T9I", "S:Q954H", "ORF9b:P10S", "S:L981F", "N:P13L", "S:G339D", "S:S375F",
                                               "S:S477N", "S:N679L", "S:S373P", "M:Q19E", "S:D796Y", "S:N969K", "S:T547K",
                                               "ORF1b:I1566V", "M:D3G", "S:G446S", "S:N440K", "M:A63T", "S:N856K", "ORF1a:A2710T",
                                               "ORF1a:I3758V", "S:E484A", "S:A67V", "S:K417N", "S:Q493R", "S:N501Y","S:Y505H", 
                                               "S:L452R", "S:P681H", "S:Q498R","S:G496S", "ORF1a:T3255I", "ORF14:G50W", "S:P681R",
                                               "N:R203M", "ORF1b:P1000L", "ORF1a:P2287S", "M:I82T", "ORF3a:S26L", "N:D63G", "N:G215C",
                                               "ORF1a:V3718A", "ORF9b:T60A"))

# mutations from https://observablehq.com/@spond/n501y-metasignature
epistatic_mutations <- c("ORF1A:.264I", "ORF1A:.680F", "ORF1A:.1654N", "ORF1A:.1000I", "ORF1A:.1187L", "ORF1A:.1707D",
                         "ORF1A:.3675-", "ORF1A:.3828F", "ORF1B:.1521I", "ORF1B:.969S", "ORF1B:.969L",
                         "ORF1B:.2020V", "S:.1176F","S:.1264L","S:.681H","S:.681R","S:.655Y","S:.18F","S:.69-","S98F","S138Y",
                         "S:.144-", "S:.501Y","S:.716I","S:.1118H","S:.243-","S:.138H","S:.701V","S:.1027I","S:.484K","S26L",
                         "S:.215Y","S:.26S",
                         "S:.215H","S:.215V","S:.681L","S:.26R","S:.417N","S:.417T","S:.215G","ORF3A:.57H","ORF3A:.171L",
                         "ORF3A:.57Y","ORF3A:.57L","E:.71L","E:.71R","E:.71T","N:.205I","N:.235F","N:.205-","N:.235L")
###
# mutations_of_interest <- epistatic_mutations
# exploration (looking for specific mutations)
indices <- with(continent_weekly_specific_mutation_summary, 
                grepl("ORF1A:.264I|ORF1A:.680F|ORF1A:.1654N|ORF1A:.1000I|ORF1A:.1187L|ORF1A:.1707D|ORF1A:.3675-|ORF1A:.3828F|ORF1B:.1521I|ORF1B:.969S|ORF1B:.969L|ORF1B:.2020V|S:.1176F|S:.1264L|S:.681H|S:.681R|S:.655Y|S:.18F|S:.69-|S98F|S138Y|S:.144-|S:.501Y|S:.716I|S:.1118H|S:.243-|S:.138H|S:.701V|S:.1027I|S:.484K|S26L|S:.215Y|S:.26S|S:.215H|S:.215V|S:.681L|S:.26R|S:.417N|S:.417T|S:.215G|ORF3A:.57H|ORF3A:.171L|ORF3A:.57Y|ORF3A:.57L|E:.71L|E:.71R|E:.71T|N:.205I|N:.235F|N:.205-|N:.235L", 
                      NewMutation)
                )

continent_weekly_specific_mutation_summary[indices, ]

#fit_mutations_count
fit_mutations_count <-  continent_weekly_specific_mutation_summary[indices, ] %>% 
  group_by(Continent, date7) %>% 
  dplyr::summarise(
    weekly_sum_of_fit_mutations = sum(n, na.rm = T)
  )


###
mutations_of_interest <- mutations_increasing_fitness

#fit_mutations_count
fit_mutations_count <- continent_weekly_specific_mutation_summary %>% 
  dplyr::filter(NewMutation %in% mutations_of_interest) %>% 
  group_by(Continent, date7) %>% 
  dplyr::summarise(
    weekly_sum_of_fit_mutations = sum(n, na.rm = T)
  )


ggplot(data = fit_mutations_count) +
  geom_line(aes(x=date7, y=weekly_sum_of_fit_mutations))+
  facet_wrap(vars(Continent))

# exploration (looking for specific mutations)
index1 <- with(continent_weekly_specific_mutation_summary, grepl("S:.18F", NewMutation))
continent_weekly_specific_mutation_summary[index1, ]

length(unique(continent_weekly_specific_mutation_summary$NewMutation))


# fitness score per week
head(fit_mutations_count)
head(continent_weekly_mutations)
# head(continent_weekly_specific_mutation_summary)
head(genome_count_per_cont_per_week)


fitness_df <- fit_mutations_count %>% 
  left_join(continent_weekly_mutations, by = c("Continent", "date7")) %>% 
  dplyr::rename(total_mutations_per_week = n) %>% 
  left_join(genome_count_per_cont_per_week, by = c("Continent", "date7")) %>% 
  dplyr::rename(total_seqs_per_week = n) %>% 
  dplyr::mutate(#fitness_score1 = weekly_sum_of_fit_mutations/total_seqs_per_week,
                #fitness_score2 = weekly_sum_of_fit_mutations/total_mutations_per_week,
                fitness_score3 = weekly_sum_of_fit_mutations/(total_seqs_per_week + log(total_mutations_per_week)),
                #fitness_score4 = (weekly_sum_of_fit_mutations/total_seqs_per_week)* (1/log(total_mutations_per_week))
                )

#rm(list=c('tmp','fit_mutations_count', 'continent_weekly_mutations', 'genome_count_per_cont_per_week'))
fitness_df %>%
  dplyr::filter(date7 <"2022-01-18")%>%
  ggplot()+
    # geom_line(aes(x=date7, y = fitness_score1, color = 'fitness_score1'))+
    # geom_line(aes(x=date7, y = fitness_score2*5, color = 'fitness_score2'))+
     geom_line(aes(x=date7, y = fitness_score3, color = 'fitness_score3')) +
    # geom_line(aes(x=date7, y = fitness_score4*5, color = 'fitness_score4')) +
    facet_wrap(vars(Continent))


#sed 's/\t/,/g;s/|/,/g' 7_seqmuts.txt | cut -d, -f12 --complement | awk 'BEGIN{FS=OFS=","} {for (i=11;i<=NF;i++) sub(/\//,",",$i)} 1' >  ~/data/seqmuts_20220128.csv


