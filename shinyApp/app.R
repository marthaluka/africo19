


#helper function for choropleth animation
setShapeStyle <- function( map, data = getMapData(map), layerId,
                           stroke = NULL, color = NULL,
                           weight = NULL, opacity = NULL,
                           fill = NULL, fillColor = NULL,
                           fillOpacity = NULL, dashArray = NULL,
                           smoothFactor = NULL, noClip = NULL, label = NULL,
                           options = NULL){
    
    options <- c(list(layerId = layerId),
                 options,
                 filterNULL(list(stroke = stroke, color = color,
                                 weight = weight, opacity = opacity,
                                 fill = fill, fillColor = fillColor,
                                 fillOpacity = fillOpacity, dashArray = dashArray,
                                 smoothFactor = smoothFactor, noClip = noClip, label = label
                 )))
    
    options <- evalFormula(options, data = data)
    options <- do.call(data.frame, c(options, list(stringsAsFactors=FALSE)))
    
    layerId <- options[[1]]
    style <- options[-1]
    if("label" %in% colnames(style)){
        labelData = style[,"label", FALSE]
        style = style[,-which(colnames(style)=="label"), FALSE]
        leaflet::invokeMethod(map, data, "setLabel", "shape", layerId, label)
    }
    leaflet::invokeMethod(map, data, "setStyle", "shape", layerId, style);
}

#helper function in JS for choropleth animation
leafletjs <-  tags$head(
    tags$script(HTML('
  
window.LeafletWidget.methods.setStyle = function(category, layerId, style){
  var map = this;
  if (!layerId){
    return;
  } else if (!(typeof(layerId) === "object" && layerId.length)){
    layerId = [layerId];
  }
  style = HTMLWidgets.dataframeToD3(style);
  layerId.forEach(function(d,i){
    var layer = map.layerManager.getLayer(category, d);
    if (layer){
      layer.setStyle(style[i]);
    }
  });
};
window.LeafletWidget.methods.setLabel = function(category, layerId, label){
  var map = this;
  if (!layerId){
    return;
  } else if (!(typeof(layerId) === "object" && layerId.length)){
    layerId = [layerId];
  }
  layerId.forEach(function(d,i){
    var layer = map.layerManager.getLayer(category, d);
    if (layer){
      layer.unbindTooltip();
      layer.bindTooltip(label[i])
    }
  });
};
'
    ))
)


#you only have to do this once!
#download.file("http://thematicmapping.org/downloads/TM_WORLD_BORDERS_SIMPL-0.3.zip" , destfile="world_shape_file.zip")
#system("unzip world_shape_file.zip")

#load spatial data
world_spdf <- readOGR( 
    dsn = getwd() , 
    layer = "TM_WORLD_BORDERS_SIMPL-0.3",
    verbose = FALSE
)


#load covid data set
source("./data_wrangling.R")

#select a certain date
selectedData <- owid_Data1[owid_Data1$date == "2020-07-15", ]

#match cases and spatial data via ISO3/iso_code
world_spdf$Total_Cases <- selectedData$total_cases_per_million[match(world_spdf$ISO3, selectedData$iso_code)]

#create label texts
world_spdf@data$LabelText <- paste0(
    "<b>Country:</b> ", world_spdf@data$NAME,"<br>", 
    "<b>Cases:</b> ", format(world_spdf@data$Total_Cases, nsmall=0, big.mark=","))

#define colorpalette for chart legend
#paletteBins <- c(0, 50000, 100000, 250000, 500000, 1000000, 2500000, 5000000, 10000000)
colorPalette <- colorBin(palette = "YlOrBr", 
                         domain = owid_Data1$total_cases_per_million,
                         na.color = "transparent") 
#na.color = "transparent",
#bins = paletteBins)

#shiny UI
ui <- fluidPage(
    leafletjs,
    titlePanel("Total Cases"),
    
    sidebarPanel(width = 2,
                 
                 radioButtons(inputId = "mapType",
                              label = "Select Map Type",
                              choices = c("Markers", "Choropleth"),
                              selected = "Markers",
                              inline = TRUE),
                 
                 radioButtons(inputId = "frequency",
                              label = "Select Data Frequency",
                              choices = c("days", "weeks"),
                              selected = "weeks",
                              inline = TRUE
                 ),
                 
                 uiOutput("dateUI")
                 
    ),
    
    mainPanel(width = 15,
              
              leafletOutput("map", width = "70%", height = "750px")
              
    )
)


#shiny server
server <- function(input, output, session) {
    
    #create slider input depending on data frequency
    observe({
        
        allDates <- unique(owid_Data1$date)
        eligibleDates <- allDates[xts::endpoints(allDates, on = input$frequency)]
        
        if(input$frequency == "weeks"){
            stepSize = 7
        }else{
            stepSize = 1
        }
        
        output$dateUI <- renderUI({
            sliderInput("dateSel", "Date",
                        min = min(eligibleDates),
                        max = max(eligibleDates),
                        value = min(eligibleDates),
                        step = stepSize,
                        timeFormat = "%b %y",
                        animate = animationOptions(interval = 500, loop = FALSE)
            )
        })
    })
    
    #filter data depending on selected date
    filteredData <- reactive({
        req(input$dateSel)
        owid_Data1[owid_Data1$date == input$dateSel, ]
    })
    
    #create the base leaflet map
    output$map <- renderLeaflet({
        
        leaflet(world_spdf) %>% 
            addTiles()  %>% 
            setView(lat = 0, lng = 0, zoom = 2) %>%
            
            addPolygons( 
                layerId = ~ISO2,
                fillColor = "lightgray", 
                stroke = TRUE, 
                fillOpacity = 1, 
                color = "white", 
                weight = 1
            ) %>%
            
            #need to specify the leaflet::addLegend function here to avoid ambiguity with the xts::addLegend function
            leaflet::addLegend(pal = colorPalette, 
                               values = owid_Data1$total_cases_per_million, 
                               opacity = 0.9, title = "Total Cases per million", 
                               position = "bottomleft",
                               na.label = "No data")
        
    })
    
    
    #prepare data depending on selected date and draw either markers or update polygons depending on the selected map type
    observe({
        
        world_spdf$Total_Cases <- filteredData()$total_cases_per_million[match(world_spdf$ISO3, filteredData()$iso_code)]
        
        world_spdf@data$LabelText <- paste0(
            "<b>Country:</b> ", world_spdf@data$NAME,"<br>", 
            "<b>Cases:</b> ", format(world_spdf@data$Total_Cases, nsmall=0, big.mark=","))
        
        if(input$mapType == "Markers"){
            
            leafletProxy("map", data = world_spdf) %>%
                clearMarkers() %>%
                setShapeStyle(layerId = ~ISO2, fillColor = "lightgray") %>%
                addCircleMarkers(lng = ~LON,
                                 lat = ~LAT,
                                 radius = ~log(Total_Cases) * 2,
                                 weight = 1,
                                 opacity = 1,
                                 color = ~ifelse(Total_Cases > 0, "black", "transparent"),
                                 fillColor = ~ifelse(Total_Cases > 0, colorPalette(Total_Cases), "transparent"),
                                 fillOpacity = 0.8,
                                 label = ~lapply(LabelText, htmltools::HTML))
            
        }else if(input$mapType == "Choropleth"){
            
            leafletProxy("map", data = world_spdf) %>%
                clearMarkers() %>%
                setShapeStyle(layerId = ~ISO2, fillColor = ~ifelse(Total_Cases > 0, colorPalette(Total_Cases), "lightgray"), label = world_spdf$LabelText)
            
        }
    })
}

shinyApp(ui, server) 