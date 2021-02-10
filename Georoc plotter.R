# moorelr/georoc-parser is licensed under The MIT License
# Copyright 2019 Lowell R. Moore

# TO-DO
#   - Default to no beep
#   - Better commenting and text output
#   - Option to load from local georoc file
#   - All plots on one figure
#   - Setting for dot size, color, shape
#   - Default to more sensible elements
#   - Plotting in separate script

# Read info from settings file ####

# To import the settigs, I do a brute-force count of the rows in the text file. Sorry!
skip_TAS <- 66
nrows_settings <- 10

Settings <- read.table("Settings.txt", header = FALSE, sep = ":", skip = 2, nrows = nrows_settings, stringsAsFactors = FALSE)
rownames(Settings) <- Settings[,1]

path_in <- Settings["Import path",2]
location_filt <- Settings["Location filter", 2]
column_filt <- Settings["Column filter",2]
std_fe <- as.logical(Settings["Standardize Fe",2])
std_tr <- as.logical(Settings["Standardize minor elements",2])
figs_out <- Settings["Figures",2]
year_cutoff <- as.numeric(Settings["Reference year cutoff",2])
make_noises <- as.logical(Settings["Make noises",2])
plot_cex <- as.numeric(Settings["Plot point size",2])
plot_alpha <- as.numeric(Settings["Plot point opacity",2])

print("Settings file loaded.")

# Turn on/off noise notificaitons
if("beepr" %in% rownames(installed.packages())
   & make_noises
){
  library("beepr")
  make_noises <- TRUE
} else {make_noises <- FALSE}

# Load parsed GeoRoc file
Import <- read.csv(file = "Georoc_parsed.csv", stringsAsFactors = FALSE)

# Parse and plot figures in list
if(nchar(figs_out) > 2){

  # Parse list of desired figures
  figs_list <- strsplit(figs_out, split = ", ")[[1]]
  if("TAS" %in% figs_list){
    plot_TAS <- TRUE
    figs_list <- figs_list[-which(figs_list == "TAS")]
  } else {plot_TAS <- FALSE}
  
  # Plot XY figures
  for(i in 1:length(figs_list)){
    # Save .pdf file
    file_name_i <- paste(gsub(x = figs_list[i], pattern = "\\.", replacement = "")
                         , ".pdf", sep = "")
    plot_cols_i <- strsplit(figs_list[i], split = " vs ")[[1]]
    if(length(Import[,plot_cols_i[1]]) > 0
       & length(Import[,plot_cols_i[2]]) > 0
       & any(!is.na(Import[,plot_cols_i[1]] + Import[,plot_cols_i[2]]))
       ){
      
      print(paste("Drawing plot of ", figs_list[i], "...", sep = ""))
      pdf(file = file_name_i, width = 6, height = 6, useDingbats = FALSE)
      plot(Import[,plot_cols_i[1]], Import[,plot_cols_i[2]]
          , xlab = plot_cols_i[1], ylab = plot_cols_i[2]
          , pch = 19, col = rgb(0.8, 0, 0, plot_alpha), cex = plot_cex
           )
      
      dev.off()
    }else{print(paste("Omitting plot of ", figs_list[i], "-- no data available!"))}
  }

  # Draw TAS plot and add literature data
  if(plot_TAS
     & length(Import$`SIO2.WT..`) > 0
     & length(Import$`NA2O.WT..`) > 0
     & length(Import$`K2O.WT..`) > 0
     & any(!is.na(Import$`SIO2.WT..` + Import$`K2O.WT..` + Import$`K2O.WT..`))){
    print("Drawing TAS plot...")
    pdf(file = "TAS.pdf", width = 6, height = 6, useDingbats = FALSE)
    
    TAS <- read.table("Settings.txt", header = TRUE, sep = "\t", skip = skip_TAS, nrows = 51, stringsAsFactors = FALSE)
    plot(c(35, 77), c(0, 16), cex = 0, xlab = "SiO2, wt%", ylab = "Na2O + K2O, wt%"
         , xlim = c(35, 60), ylim = c(0, 12)
    )
    tas_lines <- which(TAS$Name == "line")
    tas_text <- which(TAS$Name != "line")
    
    # Draw lines for TAS plot
    for(i in tas_lines[tas_lines %% 2 == 0]){
      rows <- i:(i-1)
      lines(round(TAS$SiO2[rows], 0), round(TAS$NaK[rows], 0))
    }
    
    # Add labels to TAS plot
    # Note: I screwed up and switched the columns for this in the
    #   SETTINGS.txt file. Don't worry about it -- it's all fine!
    text(TAS$NaK[tas_text], TAS$SiO2[tas_text]
         , labels = TAS$Name[tas_text], adj = c(0.5, 0.5), cex = 0.5
         )
    
    # Add points to TAS plot
    xs <- Import$`SIO2.WT..`
    ys <- Import$`NA2O.WT..` + Import$`K2O.WT..`
    points(xs, ys, pch = 19, col = rgb(0.8, 0, 0, plot_alpha), cex = plot_cex)
    
    dev.off()
  }else{print(paste("Omitting TAS plot based on data availability or settings."))}
}

print("Done!")
if(make_noises){beep(1)}
Sys.sleep(2)