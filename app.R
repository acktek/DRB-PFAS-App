################################################################################
#### DELAWARE RIVER BASIN PFAS APP 
################################################################################

#Suppress all warning messages
options(warn = -1, shiny.sanitize.errors = TRUE)

# Load dependencies
library(shiny)
library(leaflet)
library(dplyr)
library(sf)

# Load and preprocess data
load_data <- function() {
  simplify <- function(df, tissue = F) {
    
    if(tissue){
      
      df <- data.frame(
        agency = ifelse(
          grepl("DRBC|Delaware River Basin Commission", df$OrganizationFormalName),
          "DRBC",
          ifelse(
            grepl("USGS|U.S. Geological Survey", df$OrganizationFormalName),
            "USGS",
            ifelse(
              grepl("USEPA", df$OrganizationFormalName),
              "USEPA",
              ifelse(
                grepl("NJDEP|New Jersey Department of Environmental Protection", df$OrganizationFormalName),
                "NJDEP","Other")
            )
          )
        ),
        loc = df$Location.Name,
        lat = round(df$Latitude,3),
        lon = round(df$Longitude,3),
        yr = df$Year,
        conc = df$Result.Measure.Value..ppt.,
        chem = df$PFAS.Chemical.Name,
        abbrev = df$PFAS.abbrev,
        group = df$PFAS.group,
        species = df$Fish.Species.Common #ifelse(is.na(df$Fish.Species.Common),"No ID",df$Fish.Species.Common)
      )
      
      df <- df[!duplicated(df[, c(
        "agency",
        "loc",
        "yr",
        "conc",
        "chem",
        "abbrev",
        "group",
        "species"
      )]), ]
      
      df <- aggregate(conc ~ loc + lat + lon + yr + abbrev + group + agency + species, df, mean, na.rm = TRUE)
      df.frac <- aggregate(conc ~ loc + lat + lon + yr + abbrev + group + agency + species, df, function (x) x/sum(x, na.rm = TRUE))
      df.group <- aggregate(conc ~ loc + lat + lon + yr + group + agency + species, df, sum, na.rm = TRUE)
      df.group.frac <- aggregate(conc ~ loc + lat + lon + yr + group + agency + species, df, function (x) x/sum(x, na.rm = TRUE))
      df.all <- aggregate(conc ~ loc + lat + lon + yr + agency + species, df, sum, na.rm = TRUE)
      list(raw = df, grouped = df.group, all = df.all)
      
    } else {
      
      df <- data.frame(
        agency = ifelse(
          grepl("DRBC|Delaware River Basin Commission", df$OrganizationFormalName),
          "DRBC",
          ifelse(
            grepl("USGS|U.S. Geological Survey", df$OrganizationFormalName),
            "USGS",
            ifelse(
              grepl("USEPA", df$OrganizationFormalName),
              "USEPA",
              ifelse(
                grepl("NJDEP|New Jersey Department of Environmental Protection", df$OrganizationFormalName),
                "NJDEP","Other")
            )
          )
        ),
        loc = df$Location.Name,
        lat = round(df$Latitude,3),
        lon = round(df$Longitude,3),
        yr = df$Year,
        conc = df$Result.Measure.Value..ppt.,
        chem = df$PFAS.Chemical.Name,
        abbrev = df$PFAS.abbrev,
        group = df$PFAS.group
      )
      
      df <- df[!duplicated(df[, c(
        "agency",
        "loc",
        "yr",
        "conc",
        "chem",
        "abbrev",
        "group"
      )]), ]
      
      df <- aggregate(conc ~ loc + lat + lon + yr + abbrev + group + agency, df, mean, na.rm = TRUE)
      df.frac <- aggregate(conc ~ loc + lat + lon + yr + abbrev + group + agency, df, function (x) x/sum(x, na.rm = TRUE))
      df.group <- aggregate(conc ~ loc + lat + lon + yr + group + agency, df, sum, na.rm = TRUE)
      df.group.frac <- aggregate(conc ~ loc + lat + lon + yr + group + agency, df, function (x) x/sum(x, na.rm = TRUE))
      df.all <- aggregate(conc ~ loc + lat + lon + yr + agency, df, sum, na.rm = TRUE)
      list(raw = df, grouped = df.group, all = df.all)
      
    }
    
  }
  
  # Load raw data
  water_df    <- read.csv(list.files(pattern = "^PFAS_water_data_DRB"))
  gw_df    <- read.csv(list.files(pattern = "^PFAS_ground_water_data_DRB"))
  sed_df      <- read.csv(list.files(pattern = "^PFAS_sediment_data_DRB"))
  tissue_df   <- read.csv(list.files(pattern = "^PFAS_tissue_data_DRB"))
  
  # Convert ppt → ppb (divide by 1000)
  sed_df$Result.Measure.Value..ppt.    <- sed_df$Result.Measure.Value..ppt. / 1000
  tissue_df$Result.Measure.Value..ppt. <- tissue_df$Result.Measure.Value..ppt. / 1000
  
  # Return processed
  list(
    Water    = simplify(water_df),
    GW    = simplify(gw_df),
    Sediment = simplify(sed_df),
    Tissue   = simplify(tissue_df, tissue = T)
  )
  
}

# Legend breaks and labels per dataset
legend_breaks <- list(
  Water   = c(0, 1, 2, 5, 10, 25, 50, 100, 1000, Inf),
  GW   = c(0, 1, 2, 5, 10, 25, 50, 100, 1000, Inf),
  Sediment = c(0, 10, 50, 100, 250, 500, 1000, 2500, 5000, Inf)/1000,
  Tissue   = c(0, 100, 500, 1000, 5000, 10000, 25000, 50000, 100000, Inf)/1000
)

legend_labels <- list(
  Water   = c("B.D. - 1", ">1 - 2", ">2 - 5", ">5 - 10", ">10 - 25", ">25 - 50", ">50 - 100", ">100 - 1000", ">1000"),
  GW   = c("B.D. - 1", ">1 - 2", ">2 - 5", ">5 - 10", ">10 - 25", ">25 - 50", ">50 - 100", ">100 - 1000", ">1000"),
  Sediment = c("B.D. - 0.01", ">0.01 - 0.05", ">0.05 - 0.10", ">0.10 - 0.25", ">0.25 - 0.50", ">0.50 - 1.00", ">1.00 - 2.50", ">2.50 - 5.00", ">5.00"),
  Tissue   = c("B.D. - 0.10", ">0.10 - 0.50", ">0.50 - 1.00", ">1.00 - 5.00", ">5.00 - 10.0", ">10.0 - 25.0", ">25.0 - 50.0", ">50.0 - 100.0", ">100.0")
)

all_data <- load_data()

# List of agencies
agencies <- sort(unique(unlist(
  lapply(all_data, function(media) media$raw$agency)
)))


drbbnd <- read_sf("drb_bnd_arc.shp")
drbbnd <- st_transform(drbbnd, 4326)
huc <- read_sf("drbhuc12.shp")
huc <- st_transform(huc, 4326)
rm <- read_sf("delrivRM10th.shp") 
rm <- st_transform(rm, 4326)
delrivbay <- read_sf("delrivbay.shp")
delrivbay <- st_transform(delrivbay, 4326) 
#pacz <- st_read("pacz_shapefile.shp")
#pacz <- st_transform(pacz, 4326)
#pacz <- pacz[st_coordinates(st_centroid(pacz))[, "X"] >= st_bbox(drbbnd)["xmin"], ]

# Function to plot PFAS compound by River Mile
pfas_by_rm = function(df, chem, data_type, dataset){
  
  # # data frame
  # #df <- df[df$agency == "DRBC", ]
  # if (data_type == "raw") {
  #   df <- df[df$abbrev == chem, ]
  # } 
  # if (data_type == "grouped") {
  #   df <- df[df$group == chem, ]
  # } 
  # 
  # req(nrow(df) > 0)
  
  # convert water to sf (X = lon, Y = lat)
  df_sf <- suppressMessages(st_as_sf(
    df,
    coords = c("lon", "lat"),
    crs = 4326,
    remove = FALSE
  ))
  
  # determine if points are on the main stem
  inside_idx <- st_within(df_sf, st_make_valid(delrivbay), sparse = FALSE)[, 1]
  df_sf <- df_sf[inside_idx, ]
  
  if(nrow(df_sf) == 0){
    # y axis label
    ylab_text <- if (data_type == "all") {
      paste0("∑PFAS")
    } else {
      paste0(chem)
    }
    
    if(dataset == "Water" | dataset == "GW"){
      ylab_text <- paste0(ylab_text, " (ng/L)")
    } else {
      ylab_text <- paste0(ylab_text, " (ng/g)")
    }
    
    par(mar = c(8, 5, 4, 2) + 0.1,
        cex.axis = 1.4,
        cex.lab  = 1.6) 
    plot(
      0,0, 
      xaxs = "i", yaxs = "i",
      xlim = c(140, 0), 
      ylim = c(0, 1),
      xlab = "River Mile", 
      ylab = ylab_text,
      col = rgb(1,1,1,0)
    )
    text(70,0.5, "Currently, there is no data for this selection!")
  } else {
    # assign nearest river mile to each water point
    nearest_idx <- st_nearest_feature(df_sf, rm)
    df_sf$RM <- rm$RM[nearest_idx]
    
    # return to regular dataframe if desired
    df <- st_drop_geometry(df_sf)
    # 
    # # years present
    # yrs <- sort(unique(df$yr))
    # yrs <- as.character(yrs)
    # 
    # # Spectral palette (n colors)
    # cols <- colorRampPalette(c(
    #   "#9E0142", "#F46D43", "#FDAE61",
    #   "#FEE08B", "#ABDDA4", "#66C2A5",
    #   "#3288BD", "#5E4FA2"))(length(yrs))
    # 
    # # named color vector
    # yr_cols <- setNames(cols, yrs)
    
    # define full year range
    all_yrs <- as.character(seq(2004, as.numeric(substr(Sys.Date(),1,4)), 1))
    
    # fixed palette with one color per year
    all_cols <- colorRampPalette(c(
      "#9E0142", "#F46D43", "#FDAE61",
      "#FEE08B", "#ABDDA4", "#66C2A5",
      "#3288BD", "#5E4FA2"
    ))(length(all_yrs))
    
    # named vector: names are years
    yr_cols_all <- setNames(all_cols, all_yrs)
    
    # years present in data
    yrs <- as.character(sort(unique(df$yr)))
    
    # subset palette to years actually present
    yr_cols <- yr_cols_all[yrs]
    
    
    # y axis label
    ylab_text <- if (data_type == "all") {
      paste0("∑PFAS")
    } else {
      paste0(chem)
    }
    
    if(dataset == "Water" | dataset == "GW"){
      ylab_text <- paste0(ylab_text, " (ng/L)")
    } else {
      ylab_text <- paste0(ylab_text, " (ng/g)")
    }
    
    # Set arbitrary below detection value
    bd.val <- ifelse(sum(df$conc) == 0,
                     (((0.01/0.25) + 0.01) * 0.5) * 0.5,
                     (min(df$conc[df$conc!=0], na.rm = T) * 0.25 + min(df$conc[df$conc!=0], na.rm = T) * 0.5)/2)
    
    # Set plot limits
    ymin <- ifelse(sum(df$conc) == 0,
                   0.01,
                   min(df$conc[df$conc!=0], na.rm = T) * 0.25)
    
    if(dataset == "Water" | dataset == "GW"){
      ymax <- ifelse(sum(df$conc) == 0,
                     0.1,
                     round((max(df$conc[df$conc!=0], na.rm = T) + 1500) / 1000, 1) * 1000)
    } else if (dataset == "Sediment"){
      ymax <- ifelse(sum(df$conc) == 0,
                     0.1,
                     round((max(df$conc[df$conc!=0], na.rm = T) + 25) / 10, 1) * 10)
    } else {
      ymax <- ifelse(sum(df$conc) == 0,
                     0.1,
                     round((max(df$conc[df$conc!=0], na.rm = T) + 100) / 100, 1) * 100)
    }
    
    
    # River Mile plot
    par(mar = c(8, 5, 4, 2) + 0.1,
        cex.axis = 1.4,
        cex.lab  = 1.6) 
    plot(
      df$RM, ifelse(df$conc == 0, bd.val, df$conc), 
      log = 'y', yaxt = 'n',
      xaxs = "i", yaxs = "i",
      xlim = c(140, 0), 
      ylim = c(ymin, ymax),
      xlab = "River Mile", 
      ylab = ylab_text,
      pch = 1, cex = 2, lwd = 2,
      col = yr_cols[as.character(df$yr)]
    )
    
    polygon(c(140,0,0,140),c(ymin,ymin,((ymin/0.25) + ymin) * 0.5, ((ymin/0.25) + ymin) * 0.5), 
            border = NULL, col = rgb(0.1,0.1,0.1,0.1))
    abline(h = ((ymin/0.25) + ymin) * 0.5, lty = 2,lwd = 1.5)
    text(5, ((ymin/0.25) + ymin) * 0.6, "B.D.", cex = 1.5)
    
    
    mtext("Delaware River Estuary", 3, adj = 0.05, cex = 2)
    #mtext("NOTE: Below-detection samples are displayed at an arbitrary position\nbeneath detected values and are not plotted at their actual concentrations.", 3, adj = 0.95)
    mtext(
      "NOTE: In order to be transparent, non-detected samples are shown at the\nbottom of the graph. However, detection limits vary across the dataset;\ntherefore, B.D. samples are shown for informational purposes only.",
      side = 3,
      adj = 0.95,
      cex = 0.9
    )
    
    maj.ticks <- c(0.01, 0.1, 1, 10, 100, 1000, 10000)
    maj.ticks <- maj.ticks[maj.ticks > ((ymin/0.25) + ymin) * 0.5]
    
    min.ticks <- c(
      seq(0.01,0.09,0.01),
      seq(0.1,0.9,0.1),
      1:9,
      10 * (1:9),
      100 * (1:9),
      1000 * (1:9),
      10000 * (1:9)
    )
    
    min.ticks <- min.ticks[min.ticks > ((ymin/0.25) + ymin) * 0.5]
    
    axis(
      side = 2,
      at = maj.ticks,
      labels = maj.ticks,
      las = 1,      # horizontal labels
      cex.axis = 1.4
    )
    
    axis(
      side = 4,
      at = maj.ticks,
      labels = FALSE,
      las = 1,      # horizontal labels
      cex.axis = 1.4
    )
    
    axis(
      side = 2,
      at = min.ticks,
      labels = FALSE,
      tcl = -0.3     # shorter tick marks
    )
    
    axis(
      side = 4,
      at = min.ticks,
      labels = FALSE,
      tcl = -0.3     # shorter tick marks
    )
    
    legend(
      "bottom",
      legend = names(yr_cols),
      col = yr_cols,
      pch = 1, pt.cex = 2, pt.lwd = 2, cex = 1.25,
      title = "",
      bty = "n", horiz = T, x.intersp = 0.8,
      inset = -0.25, xpd = NA
    )
    
    # --------------------------
    # Head of Tide (top-left, arrow points right)
    # --------------------------
    arrows(
      x0 = 120, x1 = 135,  # arrow goes to smaller RM (right visually)
      y0 = ifelse(dataset != "Sediment", 10^(log10(ymax) * 0.9), 10^(log10(ymax) * 0.825)), 
      y1 = ifelse(dataset != "Sediment", 10^(log10(ymax) * 0.9), 10^(log10(ymax) * 0.825)),
      length = 0.08, lwd = 2
    )
    
    text(
      x = 135,
      y = ifelse(dataset != "Sediment", 10^(log10(ymax) * 0.95), 10^(log10(ymax) * 0.9)),
      labels = "Delaware River @ Trenton",
      adj = c(0, 0),
      cex = 1.2
    )
    
    # --------------------------
    # Mouth (top-right, arrow points left)
    # --------------------------
    arrows(
      x0 = 20, x1 = 5,   # arrow goes to larger RM (left visually)
      y0 = ifelse(dataset != "Sediment", 10^(log10(ymax) * 0.9), 10^(log10(ymax) * 0.825)), 
      y1 = ifelse(dataset != "Sediment", 10^(log10(ymax) * 0.9), 10^(log10(ymax) * 0.825)),
      length = 0.08, lwd = 2
    )
    
    text(
      x = 5,
      y = ifelse(dataset != "Sediment", 10^(log10(ymax) * 0.95), 10^(log10(ymax) * 0.9)),
      labels = "Atlantic Ocean",
      adj = c(1, 0),
      cex = 1.2
    )
  }
}

##### Shiny App User Interface
ui <- fluidPage(
  #titlePanel("PFAS in the Delaware River Basin"),
  fluidRow(
    column(
      width = 9,
      h1("PFAS in the Delaware River Basin")
    ),
    column(
      width = 3,
      div(
        img(src = "DRBC_logo.png", height = "90px", width = "171px"),
        style = "text-align: right; margin-top: 10px;"  # aligns top-right
      )
    )
  ),
  sidebarLayout(
    sidebarPanel(
      selectInput("dataset", "Media:", choices = c("Surface Water" = "Water",
                                                   "Groundwater" = "GW",
                                                   "Sediment" = "Sediment",
                                                   "Tissue" = "Tissue")),#names(all_data)),
      radioButtons("data_type", "PFAS Data Type:", choices = c("ΣPFAS" = "all",
                                                               "Compounds" = "raw", 
                                                               "Groups" = "grouped")),
      selectInput("agency_filter", "Agency:", choices = c("All", agencies), selected = "All"),
      uiOutput("filter_select"),
      uiOutput("year_slider"),
      selectInput("mapped_sample", "Value displayed at sampling locations with multiple sampling years:", choices = c("Most Recent Sample","Maximum Observed","Average Across Years","Minimum Observed"), selected = "Most Recent Sample"),
      checkboxInput("show_huc", "Show HUC12 Averages (based on value displayed at each sampling location)", value = FALSE),
      #checkboxInput("show_trend", "Show Trends Over Time (using most recent samples for each location)", value = FALSE),
      #checkboxInput("show_pacz", "Show PA Coastal Zone GIS layer", value = FALSE),
      checkboxInput("show_rm", "Show River Miles", value = FALSE),
      #downloadButton("downloadData", "Download Filtered Data"),
      #p(),br(),br(),br(),
      #div(
      #  img(src = "DRBC_logo.png", height = "135px", width = "257px"),
      #  style = "text-align: center;"
      #),
      p(),p(),
      div(
        style = "background-color: #f8f9f9;  /* light neutral background */
           padding: 12px; 
           border-radius: 6px; 
           border: 1px solid #d1d5d8;  /* subtle gray border */
           box-shadow: 0 1px 3px rgba(0,0,0,0.05); 
           margin-top: 20px;",
        p("Note: Data were averaged in instances where multiple samples were collected at the same location within the same calendar year.",
          style = "font-size: 13px; line-height: 1.4;"),
        p("ΣPFAS = Sum of detected PFAS compounds in a sample; not all samples were analyzed for the same compounds.",
          style = "font-size: 13px; line-height: 1.4;"),
        p("B.D. = Below analytical detection limits; varies by sample.",
          style = "font-size: 13px; line-height: 1.4;")
      ),
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("Map",
                 column(width = 6, leafletOutput("map", width = "100%", height = "700px")),
                 column(width = 6, plotOutput("pfas_hist", height = "350px"), plotOutput("yr_hist", height = "350px"))
        ),
        tabPanel("Estuary Analysis",
                 plotOutput("pfas_rm_plot", width = "100%", height = "600px")),
        tabPanel(
          "Criteria",
          tags$div(
            style = "
      background: linear-gradient(180deg, #f7f9fb 0%, #eef2f5 100%);
      padding: 30px;
      border-radius: 8px;
      box-shadow: 0 2px 8px rgba(0,0,0,0.08);
      max-width: 1100px;
      margin: 20px auto;
    ",
            
            h2("EPA Drinking Water Maximum Contaminant Level (MCL) for PFAS"),
            
            div(
              style = "width: 500px;",
              tags$table(
                border = "1", cellpadding = "6",
                style = "border-collapse: collapse; width: 100%; text-align: center;",
                
                tags$tr(
                  tags$th("Compound",style = "text-align: center;width: 40%;"),
                  tags$th("MCL Goal (ng/L)",style = "text-align: center;width: 30%;"),
                  tags$th("Final MCL (ng/L)",style = "text-align: center;width: 30%;")
                ),
                
                tags$tr(tags$td("PFOA"),tags$td("0"),tags$td("4.0")),
                tags$tr(tags$td("PFOS"),tags$td("0"),tags$td("4.0")),
                tags$tr(tags$td("PFHxS"),tags$td("10"),tags$td("10")),
                tags$tr(tags$td("PFNA"),tags$td("10"),tags$td("10")),
                tags$tr(tags$td("HFPO-DA (GenX)"),tags$td("10"),tags$td("10")),
                tags$tr(
                  tags$td("Mixture of PFHxS, PFNA, HFPO-DA, PFBS"),
                  tags$td("1 (Hazard Index)"),
                  tags$td("1 (Hazard Index)")
                )
              )
            ),
            
            br(),
            
            h2("EPA Draft Human Health Water Quality Criteria (HHC) for PFAS"),
            
            tags$table(
              border = "1", cellpadding = "6",
              style = "width: 450px; border-collapse: collapse; text-align: center;",
              tags$tr(
                tags$th("Compound",  style = "text-align: center;min-width: 150px;"),
                tags$th("Water + Organism HHC (ng/L)", style = "text-align: center;min-width: 150px;"),
                tags$th("Organism Only HHC (ng/L)", style = "text-align: center;min-width: 150px;")
              ),
              tags$tr(tags$td("PFOA"), tags$td("0.0009"), tags$td("0.0036")),
              tags$tr(tags$td("PFOS"), tags$td("0.06"), tags$td("0.07")),
              tags$tr(tags$td("PFBS"), tags$td("400"), tags$td("500"))
            ),
            
            br(),
            
            p("Note: These HHC values are draft, non-regulatory criteria. States may adopt their own values."),
            
            br(),
            
            p("The following resources provide more information about these criteria and what the numbers mean:"),
            
            tags$ul(
              tags$li(
                tags$a(
                  href = "https://www.epa.gov/ground-water-and-drinking-water/national-primary-drinking-water-regulations",
                  "National Primary Drinking Water Regulations | US EPA",
                  target = "_blank"
                )
              ),
              tags$li(
                tags$a(
                  href = "https://www.epa.gov/system/files/documents/2024-12/draft-hhc-pfas-tech-fact-sheet.pdf",
                  "Technical Fact Sheet: Draft National Recommended Human Health Ambient Water Quality Criteria for PFOA, PFOS, and PFBS",
                  target = "_blank"
                )
              ),
              tags$li(
                tags$a(
                  href = "https://www.epa.gov/wqc/human-health-water-quality-criteria-pfas",
                  "Human Health Water Quality Criteria for PFAS | US EPA",
                  target = "_blank"
                )
              )
            )
          )
        ),
        tabPanel(
          "History",
          tags$div(
            style = "
      background: linear-gradient(180deg, #f7f9fb 0%, #eef2f5 100%);
      padding: 30px;
      border-radius: 8px;
      box-shadow: 0 2px 8px rgba(0,0,0,0.08);
      max-width: 1100px;
      margin: 20px auto;
      position: relative;
      overflow: hidden;
    ",
            
            tags$div(
              style = "
        float: right;
        width: 480px;
        margin: 10px 0 10px 20px;
        text-align: center;
      ",
              
              tags$img(
                src = "pfas_archive_image.png",
                style = "width: 100%; height: auto; border: 1px solid #ccc;"
              ),
              
              tags$div(
                style = "
          font-size: 11px;
          font-style: italic;
          margin-top: 4px;
          color: #555;
        ",
                "Source: Hagley Museum and Library Digital Archive"
              )
            ),
            
            h2("History of PFAS in the Basin"),
            
            p(
              "Per- and polyfluoroalkyl substances (PFAS) are a group of more than 13,000 synthetic organic compounds characterized by exceptionally strong carbon–fluorine bonds. These bonds make PFAS highly resistant to heat, water, and oil, and protect them from natural biological and chemical degradation. As a result, PFAS persist for long periods in the environment and are often referred to as 'forever chemicals.' This persistence, together with the wide range of physical and chemical properties across PFAS compounds, allows some to accumulate in sediments and organisms while others are readily transported by air or water, ensuring they can be found almost everywhere in the environment."
            ),
            
            p(
              "The Delaware River Basin holds a unique place in the history of PFAS. It was at DuPont’s Chambers Works facility in Deepwater, New Jersey, that the first commercially produced PFAS compound—polytetrafluoroethylene (PTFE)—was accidentally discovered in 1938. Later marketed as Teflon™, this discovery helped establish the region as a global center for fluoropolymer research and production. For decades, the basin served as both a laboratory and manufacturing hub for the expansion of PFAS into thousands of consumer and industrial applications."
            ),
            
            p(
              "This long history of industrial and commerical activity across the Delaware River Basin has left a legacy of PFAS pollution that is heavily concentrated in the Lehigh Valley and extends south through the estuary. According to the EPA’s ECHO database, nearly 3,000 industrial sites in the basin have been identified as potentially involved in the production, manufacture, or use of PFAS, with more than 1,500 facilities still active today. Together, these historical and ongoing sources present a significant challenge to protecting water resources across the basin, which is shared by more than 13 million people across four states."
            ),
            
            
            br(),
            
            p(
              tags$strong("Historical reference: "),
              tags$a(
                href = "https://digital.hagley.org/PAM_08024596",
                "Polytetrafluoroethylene (Teflon), 1938 – Hagley Museum and Library Digital Archive",
                target = "_blank"
              )
            )
          )
        )
        ,
        tabPanel(
          "More Information",
          tags$div(
            style = "
      background: linear-gradient(180deg, #f7f9fb 0%, #eef2f5 100%);
      padding: 30px;
      border-radius: 8px;
      box-shadow: 0 2px 8px rgba(0,0,0,0.08);
      max-width: 1100px;
      margin: 20px auto;
    ",
            
            tags$h2(
              "About the Delaware River Basin PFAS Data App",
              style = "color: #003366; margin-top: 0;"
            ),
            
            tags$p(
              style = "font-size: 15px; line-height: 1.6;",
              "This interactive application provides basin-wide access to publicly available monitoring data for ",
              tags$strong("per- and polyfluoroalkyl substances (PFAS)"),
              " in surface water, groundwater, sediment, and biological tissue across the Delaware River Basin."
            ),
            
            tags$p(
              style = "font-size: 15px; line-height: 1.6;",
              "The app was originally developed by scientists at the ",
              tags$strong("Delaware River Basin Commission (DRBC)"),
              " as an internal tool to better understand the scope, distribution, and variability of PFAS contamination throughout the basin. 
      Recognizing the importance of transparency and shared science, DRBC has made this tool publicly available to support state agencies, researchers, utilities, and other stakeholders."
            ),
            
            tags$hr(style = "margin: 25px 0;"),
            
            tags$h3(
              "Data Sources and Scope",
              style = "color: #003366;"
            ),
            
            tags$p(
              style = "font-size: 15px; line-height: 1.6;",
              "All data displayed in this application are ",
              tags$strong("publicly available"),
              " and were retrieved from the U.S. Environmental Protection Agency’s ",
              tags$em("Water Quality Portal"),
              " or the U.S. Geological Survey's ",tags$em("Water Data API."),"The dataset currently includes ",
              tags$strong("over 900 ΣPFAS samples"),
              ", most of which are surface and ground water samples (>500), in addition to almost 300 tissue samples from 17 different aquatic species."
            ),
            
            tags$div(
              style = "
        background-color: #ffffff;
        border-left: 5px solid #005ea2;
        padding: 15px 20px;
        margin: 20px 0;
      ",
              tags$p(
                style = "font-size: 14px; margin: 0; line-height: 1.6;",
                tags$strong("Important note on agency attribution: "),
                "The",
                tags$em("Agency"),
                "filter in this app reflects the organization that",
                tags$strong("handled the sampling and/or uploaded the data to WQP."),
                "In many cases this will appear as DRBC or USGS, even when sampling programs were coordinated with or funded by other federal, state, or local agencies. 
        The listed agency does not necessarily represent sole program ownership."
              )
            ),
            
            tags$hr(style = "margin: 25px 0;"),
            
            tags$h3(
              "How This App Supports DRBC’s PFAS Roadmap",
              style = "color: #003366;"
            ),
            
            tags$p(
              style = "font-size: 15px; line-height: 1.6;",
              "PFAS are highly persistent, mobile, and bioaccumulative chemicals that present long-term challenges for water quality management. 
      DRBC’s PFAS program emphasizes a ",
              tags$strong("phased, science-based approach"),
              " that includes expanding monitoring, improving data integration, identifying potential sources, and supporting regulatory and management decisions across the basin."
            ),
            
            tags$p(
              style = "font-size: 15px; line-height: 1.6;",
              "This application directly supports those goals by:",
              tags$ul(
                style = "margin-top: 10px;",
                tags$li("Centralizing PFAS data across various agencies and jurisdictions"),
                tags$li("Allowing users to explore spatial and temporal patterns"),
                tags$li("Improving transparency and accessibility of public PFAS data"),
                tags$li("Providing a foundation for future analyses, prioritization, and decision-making")
              )
            ),
            
            tags$p(
              style = "font-size: 15px; line-height: 1.6;",
              "As monitoring continues and analytical methods evolve, this tool is expected to grow alongside DRBC’s long-term PFAS strategy."
            ),
            
            tags$hr(style = "margin: 5px 0;"),
            
            tags$p(
              style = "font-size: 15px; line-height: 1.6;",
              "For additional background on PFAS in the Delaware River Basin, including DRBC’s long-term roadmap and ongoing initiatives, please visit the official webpage: ",
              tags$a(
                href = "https://www.nj.gov/drbc/programs/quality/pfas.html",
                target = "_blank",
                style = "font-weight: bold; color: #005ea2; text-decoration: none;",
                "https://www.nj.gov/drbc/programs/quality/pfas.html"
              )
            )
          )
        )
        
      )
    )
  ),
  
  # Text describing PFAS app
  #hr(), # Optional horizontal rule for separation
  div(
    p("WARNING: Development of this app is ongoing; information is subject to change."),
    p("Development of this app was funded by EPA's Section 106 Water Pollution Control Grant Program (I-98339317-3)"),
    br(),
    p("For more information, contact:"),
    p("Jeremy Conkle, Sr. Chemist/Toxicologist (jeremy.conkle@drbc.gov)"),
    p("Matthew Amato, Water Resource Scientist (matthew.amato@drbc.gov)"),
    div(
      p("App Last Modified on January 12, 2026"),
      p("Dataset Last Modified on January 12, 2026"),
      style = "text-align: right; font-style: italic; width: 100%;"
    )
  )
)


server <- function(input, output, session) {
  
  map_bounds <- reactive({
    input$map_bounds  # This will be a named list: list(north, south, east, west)
  })
  
  
  get_full_data <- reactive({
    all_data[[input$dataset]]
  })
  
  # =========================
  # Filter Select UI
  # =========================
  output$filter_select <- renderUI({
    df <- get_full_data()[[input$data_type]]
    if (nrow(df) == 0) return(NULL)
    
    tagList(
      
      # Compound filter (raw)
      if (input$data_type == "raw") {
        selectInput(
          "chem_filter",
          "Compound Name:",
          choices = sort(unique(df$abbrev)),
          selected = "PFOA"
        )
      },
      
      # Group filter (grouped)
      if (input$data_type == "grouped") {
        selectInput(
          "group_filter",
          "Group Name:",
          choices = sort(unique(df$group)),
          selected = "PFCA"
        )
      },
      
      # Species filter (Tissue only)
      if (input$dataset == "Tissue" && "species" %in% names(df)) {
        selectInput(
          "species_filter",
          "Species:",
          choices = c("All", sort(unique(df$species))),
          selected = "All"
        )
      }
    )
  })
  
  
  # =========================
  # Year Slider UI
  # =========================
  output$year_slider <- renderUI({
    df <- get_full_data()[[input$data_type]]
    
    if (nrow(df) == 0 || is.null(df$yr)) return(NULL)
    
    years <- range(df$yr, na.rm = TRUE)
    
    sliderInput(
      "year_range",
      "Years:",
      min = years[1],
      max = years[2],
      value = years,
      sep = "",
      step = 1
    )
  })
  
  # =========================
  # Reactive filtered data
  # =========================
  filtered_data <- reactive({
    data_list <- get_full_data()
    df_raw <- data_list$raw
    df_grouped <- data_list$grouped
    df_all <- data_list$all
    
    # Filter by species (Tissue only)
    if (input$dataset == "Tissue" &&
        !is.null(input$species_filter) &&
        input$species_filter != "All") {
      
      df_raw     <- df_raw[df_raw$species == input$species_filter, ]
      df_grouped <- df_grouped[df_grouped$species == input$species_filter, ]
      df_all     <- df_all[df_all$species == input$species_filter, ]
    }
    
    # Filter by year if slider exists
    if (!is.null(input$year_range)) {
      df_raw <- df_raw[df_raw$yr >= input$year_range[1] & df_raw$yr <= input$year_range[2], ]
      df_grouped <- df_grouped[df_grouped$yr >= input$year_range[1] & df_grouped$yr <= input$year_range[2], ]
      df_all <- df_all[df_all$yr >= input$year_range[1] & df_all$yr <= input$year_range[2], ]
    }
    
    # Filter by agency
    if (!is.null(input$agency_filter) && input$agency_filter != "All") {
      df_raw <- df_raw[df_raw$agency == input$agency_filter, ]
      df_grouped <- df_grouped[df_grouped$agency == input$agency_filter, ]
      df_all <- df_all[df_all$agency == input$agency_filter, ]
    }
    
    # Filter by chemical/group only if input exists
    if (input$data_type == "raw" && !is.null(input$chem_filter) && input$chem_filter != "All") {
      df_raw <- df_raw[df_raw$abbrev == input$chem_filter, ]
    }
    if (input$data_type == "grouped" && !is.null(input$group_filter) && input$group_filter != "All") {
      df_grouped <- df_grouped[df_grouped$group == input$group_filter, ]
    }
    
    list(raw = df_raw, grouped = df_grouped, all = df_all)
  })
  
  
  combined_coords <- reactive({
    df_list <- filtered_data()
    df <- df_list[[input$data_type]]
    
    if(input$mapped_sample == "Most Recent Sample"){
      locs <- df %>% select(lat, lon) %>% distinct()
      df_latest <- suppressWarnings(df %>%
                                      group_by(lat, lon) %>%
                                      filter(!is.na(yr)) %>%
                                      filter(yr == max(yr, na.rm = TRUE)) %>%
                                      summarise(conc = mean(conc, na.rm = TRUE), .groups = "drop"))
      
      left_join(locs, df_latest, by = c("lat", "lon"))
    } else if (input$mapped_sample == "Maximum Observed"){
      locs <- df %>% select(lat, lon) %>% distinct()
      df_latest <- suppressWarnings(df %>%
                                      group_by(lat, lon) %>%
                                      summarise(conc = max(conc, na.rm = TRUE), .groups = "drop"))
      
      left_join(locs, df_latest, by = c("lat", "lon"))
    } else if (input$mapped_sample == "Average Across Years"){
      locs <- df %>% select(lat, lon) %>% distinct()
      df_latest <- suppressWarnings(df %>%
                                      group_by(lat, lon) %>%
                                      summarise(conc = mean(conc, na.rm = TRUE), .groups = "drop"))
      
      left_join(locs, df_latest, by = c("lat", "lon"))
    }else if (input$mapped_sample == "Minimum Observed"){
      locs <- df %>% select(lat, lon) %>% distinct()
      df_latest <- suppressWarnings(df %>%
                                      group_by(lat, lon) %>%
                                      summarise(conc = min(conc, na.rm = TRUE), .groups = "drop"))
      
      left_join(locs, df_latest, by = c("lat", "lon"))
    }
    
    # locs <- df %>% select(lat, lon) %>% distinct()
    # df_latest <- suppressWarnings(df %>%
    #   group_by(lat, lon) %>%
    #   filter(!is.na(yr)) %>%
    #   filter(yr == max(yr, na.rm = TRUE)) %>%
    #   summarise(conc = mean(conc, na.rm = TRUE), .groups = "drop"))
    # 
    # left_join(locs, df_latest, by = c("lat", "lon"))
  })
  
  # calculate_trend_direction <- function(df) {
  #   df %>%
  #     filter(!is.na(yr) & !is.na(conc)) %>%
  #     arrange(lat, lon, yr) %>%
  #     group_by(lat, lon) %>%
  #     summarise(
  #       first_yr = first(yr),
  #       last_yr = last(yr),
  #       first_conc = first(conc),
  #       last_conc = last(conc),
  #       trend = ifelse(last_yr != first_yr, (last_conc - first_conc) / (last_yr - first_yr), NA_real_),
  #       .groups = "drop"
  #     )
  # }
  
  # popup_summary <- function(lat, lon, df) {
  #   df_point <- df[df$lat == lat & df$lon == lon, ]
  #   if (nrow(df_point) == 0) return("No data")
  #   
  #   df_point <- df_point[order(df_point$yr), ]
  #   rows <- apply(df_point, 1, function(row) {
  #     if ("abbrev" %in% names(df_point)) {
  #       paste0("Agency : ", row["agency"], ", Year: ", row["yr"], ", Conc: ", round(as.numeric(row["conc"]), 2), " ppt")
  #     } else if ("group" %in% names(df_point)) {
  #       paste0("Agency: ", row["agency"], ", Year: ", row["yr"], ", Conc: ", round(as.numeric(row["conc"]), 2), " ppt")
  #     } else {
  #       paste0("Agency: ", row["agency"], ", Year: ", row["yr"], ", Conc: ", round(as.numeric(row["conc"]), 2), " ppt")
  #     }
  #   })
  #   
  #   paste(rows, collapse = "<br>")
  # }
  
  ##### River Mile Plot
  output$pfas_rm_plot <- renderPlot({
    
    if(nrow(filtered_data()[[input$data_type]]) > 0){
      
      df <- filtered_data()[[input$data_type]]
      
      chem_arg <- switch(
        input$data_type,
        raw     = input$chem_filter,
        grouped = input$group_filter,
        all     = NULL
      )
      
      pfas_by_rm(
        df        = df,
        chem      = chem_arg,
        data_type = input$data_type,
        dataset = input$dataset
      )
    } else {
      chem_arg <- switch(
        input$data_type,
        raw     = input$chem_filter,
        grouped = input$group_filter,
        all     = NULL
      )
      
      # y axis label
      ylab_text <- if (input$data_type == "all") {
        paste0("∑PFAS")
      } else {
        paste0(chem_arg)
      }
      
      if(input$dataset == "Water"){
        ylab_text <- paste0(ylab_text, " (ng/L)")
      } else {
        ylab_text <- paste0(ylab_text, " (ng/g)")
      }
      
      par(mar = c(8, 5, 4, 2) + 0.1,
          cex.axis = 1.4,
          cex.lab  = 1.6) 
      plot(
        0,0, 
        xaxs = "i", yaxs = "i",
        xlim = c(140, 0), 
        ylim = c(0, 1),
        xlab = "River Mile", 
        ylab = ylab_text,
        col = rgb(1,1,1,0)
      )
      text(70,0.5, "Currently, there is no data for this selection!")
    }
  })
  
  #Find data points currently in mapview
  map_filtered_data <- reactive({
    df <- filtered_data()[[input$data_type]]
    
    # If no bounds, return full filtered data
    if (is.null(input$map_bounds) || nrow(df) == 0) return(df)
    
    bounds <- input$map_bounds
    df[df$lat >= bounds$south & df$lat <= bounds$north &
         df$lon >= bounds$west  & df$lon <= bounds$east, ]
  })
  
  
  output$pfas_hist <- renderPlot({
    
    df <- map_filtered_data()
    if (nrow(df) == 0 || all(is.na(df$conc))) return(NULL)
    
    vals <- df$conc
    vals <- vals[!is.na(vals)]
    
    unit <- ifelse(input$dataset %in% c("Sediment", "Tissue"), "ng/g", "ng/L")
    
    # ---- Summary statistics (always safe) ----
    n_tot <- length(vals)
    n_bd  <- sum(vals == 0)
    pct_bd <- if (n_tot > 0) round(100 * n_bd / n_tot, 1) else 0
    
    # ---- Expand right margin to fit legend ----
    par(mar = c(5, 4, 4, 8.5) + 0.1)
    
    # ============================================================
    # CASE 1: No detected values (all B.D. or zeros)
    # ============================================================
    if (all(vals == 0)) {
      
      stats_txt <- c(
        paste0("n = ", n_tot),
        paste0("B.D. = ", n_bd, " (", pct_bd, "%)")
      )
      
      plot(
        NA, NA,
        xlim = c(0, 1),
        ylim = c(0, 1),
        main = "Concentration Distribution",
        xlab = paste0(
          ifelse(input$data_type == "raw",
                 input$chem_filter,
                 ifelse(input$data_type == "grouped",
                        input$group_filter, "ΣPFAS")),
          " (", unit, ")"
        ),
        ylab = "", yaxt = 'n', xaxt = 'n', frame.plot = FALSE,
        xaxs = "i", yaxs = "i"
      )
      
      text(
        0.5, 0.5,
        "There are no samples with detections\nbased on the current selection.",
        cex = 1
      )
      
      legend(
        "topright",
        legend = stats_txt,
        title = "Summary Statistics",
        bty = "o",
        inset = c(-0.3, 0),
        xpd = TRUE,
        bg = "grey90"
      )
      
      return(invisible())
    }
    
    # ============================================================
    # CASE 2: At least one detected value → histogram
    # ============================================================
    
    vals_pos <- vals[vals > 0]
    
    mean_v <- mean(vals_pos)
    med_v  <- median(vals_pos)
    min_v  <- min(vals_pos)
    max_v  <- max(vals_pos)
    
    stats_txt <- c(
      paste0("n = ", n_tot),
      paste0("B.D. = ", n_bd, " (", pct_bd, "%)"),
      paste0("Mean = ", round(mean_v, 2)),
      paste0("Median = ", round(med_v, 2)),
      paste0("Min = ", round(min_v, 2)),
      paste0("Max = ", round(max_v, 2))
    )
    
    h <- hist(
      vals_pos,
      breaks = max(25, ceiling(max(vals_pos) / 10)),
      plot = FALSE
    )
    
    hist(
      vals_pos,
      breaks = h$breaks,
      col = "#5D6B8D",
      border = "white",
      main = "Concentration Distribution",
      xlab = paste0(
        ifelse(input$data_type == "raw",
               input$chem_filter,
               ifelse(input$data_type == "grouped",
                      input$group_filter, "ΣPFAS")),
        " (", unit, ")"
      ),
      xaxs = "i", yaxs = "i",
      xlim = c(0, max(vals_pos) * 1.25),
      ylim = c(0, max(h$counts) * 1.2)
    )
    
    legend(
      "topright",
      legend = stats_txt,
      title = "Summary Statistics",
      bty = "o",
      inset = c(-0.3, 0),
      xpd = TRUE,
      bg = "grey90"
    )
  })
  
  
  output$yr_hist <- renderPlot({
    
    df <- map_filtered_data()
    if (nrow(df) == 0 || all(is.na(df$yr))) return(NULL)
    
    yrs <- df$yr[!is.na(df$yr)]
    
    # ---- Year counts ----
    yr_tab <- table(yrs)
    tab_txt <- paste0(names(yr_tab), ": ", as.integer(yr_tab))
    
    # ---- Expand right margin ----
    par(mar = c(5, 4, 4, 8.5) + 0.1)
    
    # ---- Histogram ----
    hist(
      yrs,
      breaks = seq(min(yrs) - 0.5,
                   max(yrs) + 0.5,
                   by = 1),
      col = "#3A6B35",
      border = "white",
      main = "Sampling Years",
      xlab = "Year",
      xaxs = "i", yaxs = "i",
      xlim = c(min(yrs) - 1, max(yrs) + 1),
      ylim = c(0, max(table(yrs)) * 1.2),
      xaxt = "n"
    )
    axis(
      side = 1,
      at = seq(min(yrs) - 1, max(yrs) + 1, by = 1),
      labels = seq(min(yrs) - 1, max(yrs) + 1, by = 1)
    )
    
    # ---- Add vertical legend on right ----
    legend(
      "topright",
      legend = tab_txt,
      bty = "o",
      inset = c(-0.3, 0),
      xpd = TRUE,
      title = "Samples per Year",
      bg = "grey90"
    )
    
  })
  
  
  popup_summary <- function(lat, lon, df) {
    df_point <- df[df$lat == lat & df$lon == lon, ]
    if (nrow(df_point) == 0) return("No data")
    
    df_point <- df_point[order(df_point$yr), ]
    
    # This generates the textual summary for display at the top
    if (input$data_type == "raw") {
      data_type_name <- input$chem_filter
    } else if (input$data_type == "grouped") {
      data_type_name <- input$group_filter
    } else {
      data_type_name <- "ΣPFAS"
    }
    
    unit <- ifelse(input$dataset %in% c("Sediment","Tissue"), "ng/g", "ng/L")
    
    if(input$dataset != "Tissue") {
      # Header row
      header <- paste0(
        #"<b>Location</b>&nbsp;&nbsp;|&nbsp;&nbsp;",
        "<b>Agency</b>&nbsp;&nbsp;|&nbsp;&nbsp;",
        "<b>Year</b>&nbsp;&nbsp;|&nbsp;&nbsp;",
        "<b>", data_type_name, " (", unit, ")</b>"
      )
      
      # Data rows with zero values displayed as "B.D."
      summary_rows <- apply(df_point, 1, function(row) {
        conc_val <- as.numeric(row["conc"])
        conc_text <- ifelse(is.na(conc_val), "NA",
                            ifelse(conc_val == 0, "B.D.", round(conc_val, 2)))
        
        paste0(
          #row["loc"], "&nbsp;&nbsp;|&nbsp;&nbsp;",
          row["agency"], "&nbsp;&nbsp;|&nbsp;&nbsp;",
          row["yr"], "&nbsp;&nbsp;|&nbsp;&nbsp;",
          conc_text
        )
      })} else{
        header <- paste0(
          #"<b>Location</b>&nbsp;&nbsp;|&nbsp;&nbsp;",
          "<b>Agency</b>&nbsp;&nbsp;|&nbsp;&nbsp;",
          "<b>Year</b>&nbsp;&nbsp;|&nbsp;&nbsp;",
          "<b>Species</b>&nbsp;&nbsp;|&nbsp;&nbsp;",
          "<b>", data_type_name, " (", unit, ")</b>"
        )
        
        # Data rows with zero values displayed as "B.D."
        summary_rows <- apply(df_point, 1, function(row) {
          conc_val <- as.numeric(row["conc"])
          conc_text <- ifelse(is.na(conc_val), "NA",
                              ifelse(conc_val == 0, "B.D.", round(conc_val, 2)))
          
          paste0(
            #row["loc"], "&nbsp;&nbsp;|&nbsp;&nbsp;",
            row["agency"], "&nbsp;&nbsp;|&nbsp;&nbsp;",
            row["yr"], "&nbsp;&nbsp;|&nbsp;&nbsp;",
            row["species"], "&nbsp;&nbsp;|&nbsp;&nbsp;",
            conc_text
          )
        }) 
      }
    
    summary_html <- paste(c(header, summary_rows), collapse = "<br>")
    
    
    if (input$data_type == "all") {
      df_all <- get_full_data()
      df_frac <- df_all$raw
      df_group_frac <- df_all$grouped
      
      frac_data <- df_frac[df_frac$lat == lat & df_frac$lon == lon, ]
      group_frac_data <- df_group_frac[df_group_frac$lat == lat & df_group_frac$lon == lon, ]
      
      if (nrow(frac_data) == 0 || nrow(group_frac_data) == 0) return("No fraction data available")
      
      abbrev_data <- aggregate(conc ~ abbrev, frac_data, mean)
      group_data <- aggregate(conc ~ group, group_frac_data, mean)
      
      abbrev_data <- abbrev_data[order(-abbrev_data$conc), ]
      group_data <- group_data[order(-group_data$conc), ]
      
      abbrev_data$frac <- abbrev_data$conc / sum(abbrev_data$conc)
      group_data$frac <- group_data$conc / sum(group_data$conc)
      
      my_palette <- colorRampPalette(c(
        "#B0413E",  # Red clay
        "#D36C2D",  # Terracotta orange
        "#D9B44A",  # Goldenrod yellow
        "#3A6B35",  # Forest green
        "#4A6D8C",  # Steel blue
        "#5D6B8D",  # Slate blue
        "#6E5773"   # Muted violet
      ))
      
      abbrev_colors <- adjustcolor(my_palette(nrow(abbrev_data)), alpha = 0.6)
      group_colors <- adjustcolor(my_palette(nrow(group_data)), alpha = 0.6)
      
      create_pie_svg_with_labels <- function(data, label_col, colors) {
        radius <- 50
        cx <- 60
        cy <- 60
        angles <- cumsum(c(0, data$frac))
        paths <- c()
        labels <- c()
        
        for (i in seq_len(nrow(data))) {
          start <- angles[i]
          end <- angles[i + 1]
          mid <- (start + end) / 2
          large_arc <- ifelse(end - start > 0.5, 1, 0)
          
          x1 <- cx + radius * cos(2 * pi * start)
          y1 <- cy + radius * sin(2 * pi * start)
          x2 <- cx + radius * cos(2 * pi * end)
          y2 <- cy + radius * sin(2 * pi * end)
          
          path <- sprintf(
            "<path d='M%f,%f L%f,%f A%f,%f 0 %d,1 %f,%f Z' fill='%s' stroke='white' stroke-width='0.5'/>",
            cx, cy, x1, y1, radius, radius, large_arc, x2, y2, colors[i]
          )
          paths <- c(paths, path)
          
          # Label
          label_angle <- 2 * pi * mid
          label_x <- cx + 35 * cos(label_angle)
          label_y <- cy + 35 * sin(label_angle)
          if (!is.na(data$frac[i]) && data$frac[i] * 100 >= 1) {
            label_text <- sprintf("%s (%.1f%%)", data[[label_col]][i], data$frac[i] * 100)
            label <- sprintf(
              "<text x='%f' y='%f' font-size='4' text-anchor='middle' fill='black'>%s</text>",
              label_x, label_y, label_text
            )
            labels <- c(labels, label)
          }
        }
        
        paste0(
          "<svg width='280' height='280' viewBox='0 0 120 120'>",
          paste0(paths, collapse = ""),
          paste0(labels, collapse = ""),
          "</svg>"
        )
      }
      
      pie_abbrev <- create_pie_svg_with_labels(abbrev_data, "abbrev", abbrev_colors)
      pie_group <- create_pie_svg_with_labels(group_data, "group", group_colors)
      
      # Combine summary + pie charts
      html <- paste0(
        "<div style='
        width: 320px; 
        max-height: 600px; 
        font-size: 12px; 
        display: flex; 
        flex-direction: column; 
        align-items: center; 
        overflow-y: auto; 
        padding: 10px;
        box-sizing: border-box;
      '>",
        
        "<div style='margin-bottom: 10px; text-align: left; width: 100%;'>",
        summary_html,
        "</div>",
        
        "<div style='margin-bottom: 5px; text-align: center;'>
        <b>Avg. Composition by Compound</b><br>",
        pie_abbrev,
        "</div>",
        
        "<div style='margin-top: 5px; text-align: center;'>
        <b>Avg. Composition by Group</b><br>",
        pie_group,
        "</div>",
        
        "</div>"
      )
      
      return(html)
    }
    
    # Fallback for other data types
    return(summary_html)
  }
  
  # Plot PFAS data on Interactive Map
  output$map <- renderLeaflet({
    leaflet(options = leafletOptions(maxZoom = 16)) %>%
      setView(lng = -75.22281, lat = 40.61678, zoom = 7) %>%
      addPolylines(data = drbbnd, color = "black", weight = 2) %>%
      #addProviderTiles("CartoDB.Voyager") %>%
      #addProviderTiles("OpenStreetMap")
      addProviderTiles("USGS.USTopo", options = providerTileOptions(opacity = 0.6, maxZoom = 16))
  })
  
  observe({
    df_list <- filtered_data()
    df <- df_list[[input$data_type]]
    df_coords <- combined_coords()
    
    proxy <- leafletProxy("map")
    
    proxy %>% clearMarkers() %>% clearShapes() %>% clearControls()
    
    # add DRB boundary back after clearing shapes
    proxy %>% addPolylines(data = drbbnd, color = "black", weight = 2)
    
    if (nrow(df_coords) == 0) {
      # removed: addProviderTiles("CartoDB.Positron")
      return()
    }
    
    breaks <- legend_breaks[[input$dataset]]
    labs <- legend_labels[[input$dataset]]
    
    # if (input$show_pacz) {
    #   proxy %>%
    #     addPolygons(data = pacz, color = "black", weight = 0, 
    #                 fillColor = rgb(0.1,0.1,0.1,0.75))
    # }
    
    if (input$show_rm) {
      proxy %>%
        addLabelOnlyMarkers(
          data = rm[rm$RM == trunc(rm$RM), ],
          label = ~as.character(RM),
          labelOptions = labelOptions(
            noHide = TRUE,
            direction = "center",
            textOnly = TRUE,
            style = list(
              "font-size" = "10px",
              "font-weight" = "bold",
              "color" = "black",
              "text-shadow" = "1px 1px 2px white"
            )
          )
        )
    }
    
    if (input$show_huc) {
      # Convert points to sf
      pts_sf <- st_as_sf(df_coords, coords = c("lon", "lat"), crs = 4326, remove = FALSE)
      
      # Use st_intersects to get indices of points in polygons
      intersection_indices <- st_intersects(huc, pts_sf)
      
      # For each polygon, calculate mean concentration of intersecting points
      huc$conc <- sapply(intersection_indices, function(idxs) {
        if (length(idxs) == 0) {
          NA  # No points inside
        } else {
          mean(pts_sf$conc[idxs], na.rm = TRUE)
        }
      })
      
      # For each polygon, calculate the number of sampling locations
      huc$nsamps <- sapply(intersection_indices, function(idxs) {
        length(idxs)
      })
      
      # Create color palette from non-NA values
      values_vec <- na.omit(huc$conc)
      pal <- colorBin("inferno", bins = breaks, domain = values_vec, na.color = "#CCCCCC")
      
      # Add polygons with fill based on mean concentration
      proxy %>%
        addPolygons(
          data = huc, 
          fillColor = ~pal(conc),
          weight = 1,
          color = "grey45",
          fillOpacity = 0.5,
          popup = ~paste0("<b>HUC12 Name:</b> ", HU_12_NAME,
                          "<br><b>Avg Conc:</b> ", round(conc, 2), ifelse(input$dataset %in% c("Sediment","Tissue"), " ng/g", " ng/L"),
                          "<br><b># of Sampling Locations:<b> ", nsamps,
                          "<br><b>Special Protection Waters:</b> ", SPW,
                          "<br><b>Tidal:</b> ", Tidal,
                          "<br><b>Square Miles:</b> ", SQMI),
          label = ~{
            conc_label <- ifelse(is.na(conc), 
                                 "NA",
                                 ifelse(conc == 0, 
                                        "B.D.",
                                        paste0(round(conc, 2), 
                                               ifelse(input$dataset %in% c("Sediment","Tissue"), " ng/g", " ng/L"))
                                 )
            )
            conc_label
          }
        ) %>%
        addLegend(
          "bottomright",
          pal = pal,
          values = values_vec, 
          labFormat = function(type, cuts, p) labs,
          title = ifelse(input$dataset %in% c("Sediment","Tissue"), "Concentration ng/g", "Concentration ng/L")
        )
      
    } else {
      # if (input$show_trend) {
      #   # Add trend info
      #   df_trend <- calculate_trend_direction(df)
      #   df_coords <- left_join(df_coords, df_trend, by = c("lat", "lon"))
      #   
      #   # Define symmetric color scale around 0
      #   trend_breaks <- c(-Inf, -25, -20, -15, -10, -5, -1, 0, 1, 5, 10, 15, 20, 25, Inf)
      #   
      #   # Diverging palette: blue → white → red
      #   trend_palette_colors <- c(
      #     "#053061", "#2166AC", "#4393C3", "#92C5DE", "#D1E5F0", "#F7FBFF",
      #     "#FFFFFF",  # neutral 0
      #     "#FFF5F0", "#FEE0D2", "#FCBBA1", "#FC9272", "#FB6A4A", "#DE2D26", "#A50F15"
      #   )
      #   
      #   trend_pal <- colorBin(
      #     palette = trend_palette_colors,
      #     bins = trend_breaks,
      #     na.color = "#CCCCCC",
      #     right = FALSE  # include upper bound (so 0 maps to 0–1 bin)
      #   )
      #   
      #   proxy %>%
      #     addCircleMarkers(
      #       data = df_coords,
      #       lng = ~lon, lat = ~lat,
      #       radius = 7,
      #       color = ~trend_pal(trend),
      #       stroke = FALSE, fillOpacity = 0.8,
      #       popup = ~mapply(function(lat, lon) popup_summary(lat, lon, df), lat, lon),
      #       label = ~paste0("Trend: ", ifelse(is.na(trend), "NA", sprintf("%.3f ppt/yr", trend))),
      #       labelOptions = labelOptions(direction = "auto")
      #     ) %>%
      #     addLegend(
      #       position = "bottomright",
      #       pal = trend_pal,
      #       values = df_coords$trend,
      #       title = "Trend (ppt/year)"
      #     )
      #   
      # } else {
      # Default color by concentration
      pal <- colorBin("inferno", bins = breaks, domain = df_coords$conc, na.color = "#CCCCCC")
      proxy %>%
        addCircleMarkers(
          data = df_coords,
          lng = ~lon, lat = ~lat,
          radius = ~case_when(
            conc <= breaks[2] ~ 3,
            conc <= breaks[3] ~ 4,
            conc <= breaks[4] ~ 5,
            conc <= breaks[5] ~ 6,
            conc <= breaks[6] ~ 7,
            conc <= breaks[7] ~ 8,
            conc <= breaks[8] ~ 9,
            conc <= breaks[9] ~ 11,
            TRUE ~ 12
          ),
          color = ~pal(conc),
          stroke = FALSE, fillOpacity = 0.7,
          popup = ~mapply(function(lat, lon) popup_summary(lat, lon, df), lat, lon),
          label = ~{
            conc_label <- ifelse(is.na(conc), 
                                 "NA",
                                 ifelse(conc == 0, 
                                        "B.D.",
                                        paste0(round(conc, 2), 
                                               ifelse(input$dataset %in% c("Sediment","Tissue"), " ng/g", " ng/L"))
                                 )
            )
            conc_label
          },
          labelOptions = labelOptions(direction = "auto")
        ) %>%
        addLegend("bottomright", pal = pal, values = df_coords$conc, 
                  labFormat = function(type, cuts, p) labs,
                  title = ifelse(input$dataset %in% c("Sediment","Tissue"), "Concentration (ng/g)", "Concentration (ng/L)"))
      #}
    }
  })
  
  
  output$downloadData <- downloadHandler(
    filename = function() {
      paste("PFAS_", input$dataset, "_", input$data_type, "_filtered.csv", sep = "")
    },
    content = function(file) {
      df_list <- filtered_data()
      df <- df_list[[input$data_type]]
      write.csv(df, file, row.names = FALSE)
    }
  )
}

shinyApp(ui, server)
