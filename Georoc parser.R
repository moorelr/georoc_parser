# moorelr/georoc-parser is licensed under The MIT License
# Copyright 2019 Lowell R. Moore

# Read info from settings file
Settings <- read.table("Settings.txt", header = FALSE, sep = ":", skip = 2, nrows = 8, stringsAsFactors = FALSE)
rownames(Settings) <- Settings[,1]

path_in <- Settings["Import path",2]
location_filt <- Settings["Location filter", 2]
column_filt <- Settings["Column filter",2]
std_fe <- as.logical(Settings["Standardize Fe",2])
std_tr <- as.logical(Settings["Standardize minor elements",2])
figs_out <- Settings["Figures",2]
year_cutoff <- as.numeric(Settings["Reference year cutoff",2])
make_noises <- as.logical(Settings["Make noises",2])

print("Settings file loaded.")

# Turn on/off noise notificaitons
if("beepr" %in% rownames(installed.packages())
   & make_noises
){
  library("beepr")
  make_noises <- TRUE
} else {make_noises <- FALSE}

# Unused code to prompt user for input path
if(FALSE){
  if(interactive()){
    georoc_url <- readline("Input URL for precompilled GeoRoc data:")
  } else{
    cat("Download GEOROC file from URL? (y/n):")
    georoc_url <- readLines("stdin", n = 1)
  }
}

# Unused code to run as demo
if(FALSE){#!grepl(pattern = "http://georoc.mpch-mainz.gwdg.de", path_in)){
  print("No GeoRoc URL was entered."); Sys.sleep(0.5)
  print("Running demo with Cascades precompiled file..."); Sys.sleep(1)
  georoc_url <- "http://georoc.mpch-mainz.gwdg.de/georoc/Csv_Downloads/Convergent_Margins_comp/CASCADES.csv"
}

# Import georoc data from URL
if(grepl(pattern = "http://georoc.mpch-mainz.gwdg.de", path_in)){
  print("Downloading GeoRoc data... (this could take a minute or two)")
  Import <- read.csv(url(path_in), stringsAsFactors = FALSE, na.strings = c(""))
} else {
  print(paste("Loading GeoRoc data :", path_in))
  Import <- read.csv(path_in, stringsAsFactors = FALSE, na.strings = c(""))
}

print(paste("Data successfully loaded with", nrow(Import), "rows."))

# Remove any commas from the parsed file
flag <- which(sapply(X = Import, FUN = grepl, pattern = ","), arr.ind = TRUE)
if(nrow(flag) > 0){
  Import[flag] <- gsub(Import[flag], pattern = ",", replacement = "")
  print(paste("Removed errant commas from", nrow(flag), "cells in dataset."))
}

print(paste("File loaded with", nrow(Import), "rows."))# and columns:", paste(colnames(Import), collapse = " ")))
print("Thanks Bärbel et al!")
if(make_noises){beep(1)}
Sys.sleep(2)

# Separate references and data into separate variables
print("Parsing references...")
# Make sure the references column is labeled consistently
if("CITATION" %in% colnames(Import)){
  colnames(Import)[colnames(Import) == "CITATION"] <- "CITATIONS"
}

# Remove the dummy row that marks the beginning of the references list
flag_str <- "(References)|(Abbreviations)"
if(any(grepl(x = Import$CITATIONS, pattern = flag_str))){
  Import <- Import[-which(grepl(x = Import$CITATIONS, pattern = flag_str)),]
}
# Flag the beginning of the references based on character length, and separate from main dataframe
ref_start <- which(nchar(Import$CITATIONS) > 50)[1]
Import_refs <- Import$CITATIONS[(ref_start):length(Import$CITATIONS)]
Import <- Import[1:ref_start-1,]

# Remove duplicate references (just in case)
Import_refs <- unique(Import_refs)

# Apply location filter filters
if(nchar(location_filt) > 0){
  print(paste("Filtering data for location :", location_filt))
  flag <- which(!grepl(x = tolower(Import$LOCATION), pattern = tolower(location_filt)))
  Import <- Import[-flag,]
  print(paste("Removed", length(flag), "rows based on location filter."))
}

# Make sure all elemental values are coerced to numeric
ele_start <- which(colnames(Import) == "SIO2.WT..")
for(i in ele_start:ncol(Import)){
  #print(i)
  Import[,i] <- as.numeric(Import[,i])
}

# Parse reference list to extract year and reference code
ref_year <- numeric(0)
ref_number <- numeric(0)
for(i in 1:length(Import_refs)){
  date_start <- gregexpr(pattern = " \\[", text = Import_refs[i])[[1]][1]
  ref_year[i] <- substr(Import_refs[i], start = date_start+2, stop = date_start + 5)
  
  num_start <- gregexpr(pattern = "\\[", text = Import_refs[i])[[1]][1]
  num_stop <- gregexpr(pattern = "\\]", text = Import_refs[i])[[1]][1]
  ref_number[i] <- substr(Import_refs[i], start = num_start, stop = num_stop)
}

# flag data collected pre-1990
# Debug: i <- 11; j <- 2
old_papers <- which(as.numeric(ref_year) < year_cutoff)
is_old <- rep(FALSE, nrow(Import))
if(!is.na(year_cutoff)){
  for(i in 1:length(old_papers)){
    print(paste("Flagging data from old paper ", i, " of ", length(old_papers), " pre-", year_cutoff, " papers.", sep = ""))
    paper_i <- old_papers[i]
    num_i <- ref_number[paper_i]
    for(j in 1:length(Import$CITATIONS)){
      if(grepl(num_i, Import$CITATIONS[j], fixed = TRUE)){
        is_old[j] <- TRUE
      }
    }
  }
}
Import <- cbind(Import, is_old)
print("Done parsing references"); Sys.sleep(0.5)

# Calculate FeO_total
if(std_fe){
  print("Standardizing Fe...")
  #colnames(Import)
  FeO_total <- numeric(0)
  # Debug : i <- 11
  for(i in 1:nrow(Import)){
    if(is.na(Import$FEOT.WT..[i]) & is.na(Import$FEO.WT..[i]) & is.na(Import$FE2O3.WT..[i])){
      FeO_total[i] <- NA
      next
    }
    else if(!is.na(Import$FEOT.WT..[i])){
      FeO_total[i] <- Import$FEOT.WT..[i]#1
      next
    }
    else if(is.na(Import$FEOT.WT..[i]) & !is.na(Import$FEO.WT..[i]) & !is.na(Import$FE2O3.WT..[i])){
      FeO_total[i] <- Import$FEO.WT..[i] + (0.8998*Import$FE2O3.WT..[i])#2
      next
    }
    else if(is.na(Import$FEOT.WT..[i]) & !is.na(Import$FEO.WT..[i]) & is.na(Import$FE2O3.WT..[i])){
      FeO_total[i] <- Import$FEO.WT..[i]#3
      next
    }
    else if(is.na(Import$FEOT.WT..[i]) & is.na(Import$FEO.WT..[i]) & !is.na(Import$FE2O3.WT..[i])){
      FeO_total[i] <- 0.8998*Import$FE2O3.WT..[i]#4
      next
    }
  }
Import <- cbind(Import, FeO_total)
}else{print("Not standardizing Fe."); Sys.sleep(0.5)}

# Debug: i <- round(runif(1, 1, nrow(Import)), 0)
if(std_tr){
  print("Standardizing other trace elements...")
  for(i in 1:nrow(Import)){
    if(is.na(Import$CR2O3.WT..[i]) & !is.na(Import$CR.PPM.[i])){
      Import$CR2O3.WT..[i] <- 1.462*Import$CR.PPM.[i]*10^-4
    }
    if(is.na(Import$TIO2.WT..[i]) & !is.na(Import$TI.PPM.[i])){
      Import$TIO2.WT..[i] <- 1.668*Import$TI.PPM.[i]*10^-4
    }
    if(is.na(Import$MNO.WT..[i]) & !is.na(Import$MN.PPM.[i])){
      Import$MNO.WT..[i] <- 1.291*Import$MN.PPM.[i]*10^-4
    }
    if(is.na(Import$NIO.WT..[i]) & !is.na(Import$NI.PPM.[i])){
      Import$NIO.WT..[i] <- 1.273*Import$NI.PPM.[i]*10^-4
    }
    if(is.na(Import$K2O.WT..[i]) & !is.na(Import$K.PPM.[i])){
      Import$K2O.WT..[i] <- 1.205*Import$K.PPM.[i]*10^-4
    }
    if(is.na(Import$P2O5.WT..[i]) & !is.na(Import$P.PPM.[i])){
      Import$P2O5.WT..[i] <- 2.292*Import$P.PPM.[i]*10^-4
    }
    
    # Also standardize to PPM for Ti and K
    if(!is.na(Import$K2O.WT..[i]) & is.na(Import$K.PPM.[i])){
      Import$K.PPM.[i] <- Import$K2O.WT..[i]/(1.205*10^-4)
    }
    if(!is.na(Import$TIO2.WT..[i]) & is.na(Import$TI.PPM.[i])){
      Import$TI.PPM.[i] <- Import$TIO2.WT..[i]/(1.668*10^-4)
    }
    if(!is.na(Import$P2O5.WT..[i]) & is.na(Import$P.PPM.[i])){
      Import$P.PPM.[i] <- Import$P2O5.WT..[i]/(2.292*10^-4)
    }
  }
} else {print("Not standardizing trace elements"); Sys.sleep(0.5)}

# Calculate totals
print("Calculating totals...")
#colnames(Import)
Totals <- numeric(0)
for(i in 1:nrow(Import)){
  majors_i <- c(Import$SIO2.WT..[i]
                , Import$TIO2.WT..[i]
                , Import$AL2O3.WT..[i]
                , Import$CR2O3.WT..[i]
                , FeO_total[i]
                , Import$CAO.WT..[i]
                , Import$MGO.WT..[i]
                , Import$MNO.WT..[i]
                , Import$NIO.WT..[i]
                , Import$K2O.WT..[i]
                , Import$NA2O.WT..[i]
                , Import$P2O5.WT..[i]
                #, Import$H2O.WT..[i]
                #, Import$CO2.WT..[i]
                #, Import$F.WT..[i]
                #, Import$CL.WT..[i]
                #, Import$SO2.WT..[i]
  )
  majors_i <- majors_i[!is.na(majors_i)]
  Totals[i] <- sum(as.numeric(majors_i))
} # hist(Totals)
print("Done!"); Sys.sleep(0.5)

# Filter for missing data in selected columns
if(nchar(column_filt) > 0){
  print(paste("Filtering data for missing values :", column_filt))
  filter_cols <- strsplit(column_filt, ", ")[[1]]
  for(j in filter_cols){
    flag <- which(is.na(Import[,j]))
    if(length(flag) > 0){Import <- Import[-flag,]}
    print(paste("Removed", length(flag), "rows for missing", j))
    print(paste(nrow(Import), "rows remaining in filtered dataset."))
  }
}

# Save parsed GeoRoc files
print(paste("Saving parsed georoc file (", nrow(Import), " rows) and References...", sep = ""))
write.csv(x = Import, file = "Georoc_parsed.csv"
          , quote = FALSE, row.names = FALSE)
write.csv(x = Import_refs, file = "Georoc_refs_parsed.csv"
          , quote = FALSE, row.names = FALSE)
Sys.sleep(2)

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
          , pch = 21, bg = rgb(0.8, 0, 0, 0.8), cex = 1.2
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
    
    TAS <- read.table("Settings.txt", header = TRUE, sep = "\t", skip = 64, nrows = 51, stringsAsFactors = FALSE)
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
    points(xs, ys, pch = 21, bg = rgb(0.8, 0, 0, 0.8), cex = 1.2)
    
    dev.off()
  }else{print(paste("Omitting TAS plot based on data availability or settings."))}
}

print("Done!")
if(make_noises){beep(1)}
Sys.sleep(2)