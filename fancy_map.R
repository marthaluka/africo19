
#To create a fancy map to visualize global lineage dynamics

rm(list = ls())

## packages
require(pacman)
pacman::p_load(tidyverse, rgdal, sf, rmapshaper, ggforce, cowplot, patchwork, rnaturalearth,
               ggtext, showtext, rnaturalearthdata)


font_add_google("Pacifico", "pacifico")
font_add_google("Source Sans Pro", "source")
showtext_opts(dpi = 320)
showtext_auto(enable = TRUE)


## load data ##################################

mutationData<-read.csv("../data/metadata_20210907.csv", header = F) %>%
  dplyr::select(V2, V10, V5, V8, V11) 

column_headings <- c("gisaid_ID", "country", "lineage", "date", "continent")
names(mutationData)<-column_headings

mutationData$date<-as.Date(mutationData$date)
mutationData$date14<-as.Date(cut(mutationData$date,breaks = "2 weeks",start.on.monday = FALSE))

#clean and process data ###########
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

#
mutationData<- mutationData %>% 
  dplyr::mutate_at(
    c("continent", "country"), renameFunction
  ) %>%
  mutate(across(where(is.character),str_trim)) %>%
  dplyr::filter(continent != is.na(.))


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

table(mutationData$country)

#our geographical selction
#africo19<- c("Kenya", "Uganda", "The Gambia")
africo19 <-c("South Africa", "India", "Japan", "South Korea", "Belgium", "France",        
             "Germany", "Italy",  "Netherlands", "Norway", "Spain",  "Sweden",        
             "Switzerland", "United Kingdom", "Canada","Mexico", "United States", "Australia",     
              "Brazil")


selectData<- mutationData %>%
  dplyr::filter(
    country %in% africo19
  )

table(selectData$country)

#Top X (5 in this case) lineages per continent
vocs<-c("Alpha", "Beta", "Gamma", "Delta")

summ_lineages<- selectData %>%
  na.omit() %>%
  group_by(country, lineage) %>%
  dplyr::summarise(n=n()) %>%
  arrange(desc(n)) %>%
  slice(1:5) 
                        
                          #######add vocs if not in top 5 per country?


topLineages<-unique(summ_lineages$lineage)

#create a custom pallete to fix lineage colours across continents  #############
my_colour_palette<-c("#8DD3C7","#FFFFB3","#BEBADA","#FB8072","#80B1D3",
                     "#B3DE69","#FCCDE5","deeppink3","#BC80BD","#CCEBC5","darkorange", 
                     "goldenrod1","royalblue",
                     "cyan", "gray30", "firebrick1","gray95", "peachpuff",  "magenta")

number_of_colors_needed <- length(africo19)
  #extract qualitative colour palettes from RColorBrewer
qualitative_color_palettes<- brewer.pal.info[brewer.pal.info$category == "qual",]
  #create mash-up of colours
colour_palette <- unlist(mapply(brewer.pal, qualitative_color_palettes$maxcolors, 
                                   rownames(qualitative_color_palettes)))
  #subset to what we need
my_colour_palette<- colour_palette[c(1:number_of_colors_needed)] 
  #or randomly (will change everytime)
my_colour_palette<-sample(colour_palette,18)


df_topLineages<- selectData %>%
  dplyr::group_by(country, lineage, date14)  %>%
  dplyr::summarise(topTotal=n()) %>%
  subset(lineage %in% topLineages)  #filter based on top


df_biweeklyTotal<-selectData %>%
  drop_na(.)  %>%
  dplyr::group_by(country, date14)  %>%
  dplyr::summarise(allSeqs=n()) %>%
  left_join(df_topLineages,
            by=c("date14", "country"))

df_biweeklyTotal$topLineage_prop <- df_biweeklyTotal$topTotal/df_biweeklyTotal$allSeqs 


## map ##########
#download.file("http://thematicmapping.org/downloads/TM_WORLD_BORDERS_SIMPL-0.3.zip" , destfile="world_shape_file.zip")
#system("unzip world_shape_file.zip")

globe_shp <- st_read(dsn="./mapShapefiles/", layer="TM_WORLD_BORDERS_SIMPL-0.3")
globe <- st_transform(globe_shp, "+init=epsg:4326") 
globe <- ms_simplify(globe) # Simplify as otherwise is massive and laptop gets sad

globe$NAME

selectRegions <- globe %>% 
  dplyr::filter(NAME %in% africo19)


no_borders<- globe %>% 
  dplyr::filter(! NAME %in% africo19)
  
ken <- ne_countries(scale = "medium", returnclass = "sf", country = "Kenya")


(map <- ggplot() +
    geom_sf(data = globe, fill = NA) +
    geom_sf(data = no_borders, color = "#FAFAFA", alpha= 1) +
   #geom_sf(data=ibra_regs, fill = "#E1E3D4") +
   geom_sf(data=selectRegions, fill = my_colour_palette)+
   xlab("Longitude") +
   ylab("Latitude") + 
    #xlim(-135, 135)+
    ylim(-55, 80) +
   theme_bw() +
   theme(plot.background = element_blank(),
         panel.background = element_rect(fill = "lightblue1", colour = "powderblue"),
         axis.title = element_blank(),
         axis.text = element_blank(),
         panel.grid = element_blank(),
         axis.ticks = element_blank(),
         panel.border = element_blank())
)



## plots############
plotFunction<-function(big_data, selection, na.rm=T) {
  big_data %>% 
    dplyr::filter(country == selection) %>%
    ggplot() + 
    geom_bar(aes(x = date14, y=topLineage_prop,fill=lineage), #fill=forcats::fct_rev(topLineages) #to reverse order
             position = "stack", stat="identity", width=14,color='black')+
    theme_classic() + theme_bw() +
    ylab('Proportion') +
    scale_y_continuous(breaks = c(0,0.25, 0.50, 0.75, 1))+
    theme(legend.position = "right") +
    scale_x_date(date_labels = "%b \n \n %Y",date_breaks = "6 months", limits = as.Date(c('2019-12-01','2021-10-01')))+
    theme(axis.title = element_blank(),
          axis.text = element_text(color="black", size=2.5),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.text = element_text(size = 2.5),
          legend.title = element_text(size = 2.5))+ 
    scale_fill_brewer(palette ="Set3")+
    #scale_fill_manual(values = my_colour_palette,name ="Lineage")+
    guides(fill=guide_legend(title="Lineage"))
}


australia <- plotFunction(df_biweeklyTotal, "Australia")
  
uganda <- plotFunction(df_biweeklyTotal, "Uganda")

gambia <- plotFunction(df_biweeklyTotal, "The Gambia")




# bring map together
maps_2 <- ggdraw(xlim = c(0, 35), ylim = c(0, 30)) +
  draw_plot(map, x = 0, y = 0, width = 30, height = 30) +
  geom_polygon(aes(x = c(17,13,17.5), y = c(15.2, 16.2, 16.2)),
               fill = "lightpink2", alpha = 0.6) +
  geom_polygon(aes(x = c(18.5,20,20.5), y = c(15, 17, 12)),
               fill = "lightblue", alpha = 0.6) +
  draw_plot(kenya, x = 20, y = 11.8, width = 5.5, height = 5.5) +
  draw_plot(uganda, x = 12.5, y = 16.2, width = 5.5, height = 5) 




