# Instructions ----

# This script is designed to combine a set of several data tables as .csv files,
#   calculate normalized compositions, and plot the data on a multielement ("spider")
#   diagram.
# The workflow is as follows:
#   - Configure settings file "Element norm.txt" for a list of elements to use.
#   - Import data and manually format as necessary (i.e. it is the user's
#     responsibility to make sure that K2O is converted to K ppm... etc).
#       > It may be necessary to use the "merge_dfs" and/or "get_chunk" functions
#           to facilitate combining existing data tables.
#       > Note: Each dataset to be processed must have a "Sample" column, which may or
#           may not contain a useful sample name (this will be passed onto the final
#           normalized dataset) and a "flag" column -- the contents of which will not
#           be passed on or affect anything (it's just a placeholder)
#   - Run "init_norm" function to create a dataframe to hold normalized data
#   - Run "calc_norm" function to append imported data to the above dataframe
#       > Note: it is necessary to identify sub-samples of data at this step based on
#           the "flag" and "flag_rows" arguments of the "calc_norm" function.
#           For example, a data frame with mixed basaltic and andesitic samples may
#           require a "Basalt" and "Andesite" flag corresponding to each.
#   - Run "spider_draw" function to draw a spider diagram
#   - Run one or several of the following functions to populate the spider diagram:
#       > "spider_lines" to draw lines for a specified flag
#       > "spider_poly" to draw a polygonal field for a specified flag (and optional
#           quantile interval)
#       > The plotting parameters (line width, color, quantiles, etc.) can be found in
#         the settings file "Element norm.txt"


# To do:
# - Configure paths so everything will run from a common directory

setwd("G:/My Drive/VT Volcanoes/Hawaii/Hawaii Rejuvenated/Hawaii paper/Hawaii paper 6 -- Reviewer comments/200423 light element plots/Georoc parser")

# Functions ----

# Load settings from file
load_settings <- function(settings_path = getwd()){
  # My intent was to put some kind of useful plotting parameters in the same
  #   settings file as the norm composition. I have not done that yet, so this
  #   function doesn't do anything.
  if(FALSE){
    read.table(settings_path, stringsAsFactors = FALSE, fill = TRUE)
    settings_start <- 1
  }
}

# Initialize dataframe to hold data to be normalized
init_norm <- function(){
  # This function takes no input, and generates a detaframe to hold normalized
  #   concentrations based on those specified in the settings file "Element norm.txt"
  
  # Open the settings file
  Ele_norm <- read.table("Element norm.txt", stringsAsFactors = FALSE, skip = 1, header = FALSE, fill = TRUE)
  Ele_norm <- Ele_norm[1:(grep("\\*", Ele_norm$V1)[1]-1),] # Stop reading at first row contianing asterisks
  
  # Create an empty dataframe to hold normalized compositions.
  # The column names need to match the elements from the settings file
  # and the first row will hold the normal composition.
  ree_norm <- as.data.frame(matrix(NA, nrow = 1, ncol = nrow(Ele_norm)))
  colnames(ree_norm) <- Ele_norm$V1
  ree_norm[1,] <- Ele_norm$V2
  # Create columns for flag and sample name which are needed to plot subsets of data
  ree_norm <- cbind(data.frame(flag = "Norm comp", Sample = "Norm comp", stringsAsFactors = FALSE), ree_norm)
  
  # Currently unused: create a directory for outputs
  if(FALSE){
    directory <- "C:/Users/"
    if(!dir.exists(directory)){
      dir.create(directory)
    }
  }
  
  # Make sure concentrations are formatted as numeric
  ree_norm[,3:ncol(ree_norm)] <- as.numeric(ree_norm[,3:ncol(ree_norm)])

  return(ree_norm)
}

# Append data to a normalized dataframe
calc_norm <- function(df_norm, df_new, flag = "df_new", norm_cols = numeric(0)){
  # This function requires an existing dataframe of normalized concentrations
  # (greated by the "init_norm" function) and an existing dataframe of elemental
  # concentrations that has been properly formatted to include:
  #   - a labeled column for each element being added to the spider diagram
  #   - a column labeled "Sample"
  #   - a column labeled "flag"
  # Note: It doesn't matter if there is data in any of these columns -- they just
  #   need to exist for the function to run!
  
  # Create a dataframe that holds data being added to the spider diagram, and
  #   label it with matching columns
  df_add <- as.data.frame(matrix(data = NA, ncol = ncol(df_norm), nrow = 1))
  colnames(df_add) <- colnames(df_norm)
  
  # If matched columns have not been supplied, identify matching columns here.
  if(length(norm_cols) == 0){
    norm_cols <- get_cols(df_new, colnames(df_add))
  }
  colnames(df_new)[norm_cols] <- colnames(df_add)
  
  # set missing and zero values to NA
  # Note: this is going to eat up time fast if df_new has more than a few 10s of
  #   rows, and should be updated to use mapply()
  for(i in 1:nrow(df_new)){# Debug : i <- 1; j <- 24
    for(j in 1:ncol(df_new)){
      if(is.na(df_new[i,j])){next}
      if(is.null(df_new[i,j])){
        df_new[i,j] <- NA}
      if(df_new[i,j] == 0){
        df_new[i,j] <- NA}
    }
  }
  
  # Append data to the dataframe to be appended so that the columns are in the
  #   proper order, and then calculate a normalized composition for each row
  #   (excluding the first two columns)
  df_add <- rbind(df_add, df_new[,norm_cols])[-1,]
  for(i in 1:nrow(df_add)){
    df_add[i,3:ncol(df_add)] <- as.numeric(df_add[i,3:ncol(df_add)])/as.numeric(df_norm[1,3:ncol(df_norm)])
  }
  
  # Label the "flag" column according to the supplied label.
  df_add$flag <- flag
  
  return(df_add)
}

# Draw spider diagram
spider_draw <- function(df_norm, ylim = c(2, 130)){
  # Note: This function assumes that the first two columns in the "df_norm" dataframe
  #   are "Sample" and "flag" columns which are not used.
  
  # Debug : df_norm <- norm_data
  plot(1, 1, type = "n"
       , xlab = "Element", ylab = "Concentration/Primitive mantle"
       , xlim = c(3, ncol(df_norm)), ylim = ylim
       , log = "y", xaxt = "n"
  )
  
  # offset element names for easier readability
  even_xs <- which(1:ncol(df_norm) %% 2 == 0)[-1]
  odd_xs <- which(1:ncol(df_norm) %% 2 != 0)[-1]
  axis(side = 1, at = odd_xs, labels = colnames(df_norm)[odd_xs], padj = -0.5)
  axis(side = 1, at = even_xs, labels = colnames(df_norm)[even_xs], padj = 1)
}

# Add polygon of multiple compositions
spider_poly <- function(df_norm, flag, border = "black", col = rgb(0, 0, 0, 0.2), lwd = 2
                        , filt_quantile = c(0, 1)){
  # Note: the "polygon" function needs an ordered set of points which go around the
  #   perimeter of the desired polygon. Here, these are calculated along the "top" and
  #   "bottom" of the polygon, and plotted in a "counterclockwise" fashion.
  
  
  # Debug : df_norm <- norm_data; flag <- "HAVF"
  
  # Empty vectors for top and bottom of polygon
  poly_bottom <- rep(NA, ncol(df_norm))
  poly_top <- rep(NA, ncol(df_norm))
  
  # Flag subset of normalized data to be plotted
  flag_rows <- grep(flag, df_norm$flag)
  
  # Calculate the top and bottom of the polygon as quantiles
  for(i in 3:ncol(df_norm)){
    poly_bottom[i] <- quantile(df_norm[flag_rows, i], probs = filt_quantile[1], na.rm = TRUE)
    poly_top[i] <- quantile(df_norm[flag_rows, i], probs = filt_quantile[2], na.rm = TRUE)
  }
  
  # Add polygon to existing plot
  polygon(c(3:ncol(df_norm), ncol(df_norm):3)
          , c(poly_bottom[-c(1, 2)], rev(poly_top[-c(1, 2)]))
          , col = col, border = border)
}

# Add lines
spider_lines <- function(df_norm, flag, col = "black", lwd = 1){
  # Debug : df_norm <- norm_data; flag <- "HA-16-6"
  flag_rows <- which(df_norm$flag == flag)
  for(j in flag_rows){
    # Suppress warnings due to missing values -- not a problem
    suppressWarnings(
      lines(1:ncol(df_norm), df_norm[j,], col = col, lwd = lwd)
      )
  }
}

# Function to merge dataframes together (with some shared columns)
merge_dfs <- function(df1, df2){
  # Note: This function was written to be used on GeoRoc data, which come with formatted
  #   columns. Therefore, this function assumes that there are no variations in the
  #   spelling of column names.
  
  # New data (rows) plus NAs to compensate for missing columns
  have_cols <- intersect(colnames(df1), colnames(df2))
  missing_cols <- setdiff(colnames(df1), colnames(df2))
  missing_block <- as.data.frame(matrix(data = NA, nrow = nrow(df2), ncol = length(missing_cols)))
  colnames(missing_block) <- missing_cols
  add_rows <- cbind(df2[,have_cols], missing_block)
  
  # New data (columns) plus NAs to compensate for missing columns in the original data
  new_cols <- setdiff(colnames(df2), colnames(df1))
  missing_block <- as.data.frame(matrix(data = NA, nrow = nrow(df1), ncol = length(new_cols)))
  colnames(missing_block) <- new_cols
  add_cols <- rbind(missing_block, df2[,new_cols])
  
  # Bind everyting together
  if(nrow(add_rows) > 0){
    df3 <- rbind(df1, add_rows)
  }
  if(ncol(add_cols) > 0){
    df3 <- cbind(df3, add_cols)
  }
  
  return(df3)
}

# Another function to merge dataframes together (with some shared rows)
get_chunk <- function(df1, df2){
  # Note: This function was designed for a bunch of melt inclusion datasets, which
  #   unlike the GeoRoc data, are more likely to have rows in common, but no shared
  #   columns.
  # Note: This function assumes that the first column contains some kind of
  #   unique sample label that can be used to identify matching rows.
  #   df1 is the dataframe to which new data (df2) are being added.
  add_chunk <- as.data.frame(matrix(data = NA, nrow = nrow(df1), ncol = ncol(df2)))
  colnames(add_chunk) <- colnames(df2)
  match_rows <- match(df2$Sample, df1$Sample)
  add_chunk[match_rows,] <- df2
  add_chunk <- add_chunk[,-1]
  return(add_chunk)
}

# Function to identify columns in a dataframe for desired elements
get_cols <- function(df_in, ele_names){
  # Note: This function involves a crazy workaround: I couldn't figure out how to handle
  #   the problem of "do I use P from P2O5 vs P from P_ppm vs P from H2OP," so rather
  #   than an elegant regular expression, I decided to do some rough brute force matching
  #   and then sort out any issues with a custom menu. ("Which one do you want to use?")
  
  # Debug : df_in <- df_new; ele_names <- colnames(df_add)
  
  # Create vector to hold matched column positions
  match_cols <- numeric(0)
  for(i in 1:length(ele_names)){
    # Match using regular expressions
    match_i <- which(tolower(ele_names[i]) == tolower(colnames(df_in))
                        | grepl(x = tolower(colnames(df_in)), pattern = paste(tolower(ele_names[i]), "\\d*o*\\d*\\.*ppm\\.*", sep = ""))
                        | grepl(x = tolower(colnames(df_in)), pattern = paste(tolower(ele_names[i]), "\\d*o*\\d*\\.*wt\\.*", sep = ""))
    )
    
    # Solve multiple matches with user input
    if(length(match_i) > 1){
      match_i <- match_i[menu(choices = colnames(df_in)[match_i], graphics = TRUE, title = "Use which column?")]
    }
    # Warn for failed match attempts
    if(length(match_i) < 1){
      warning(paste("Could not find matching column for", ele_names[i]))
    }
    # Output showing which matches were used
    print(paste(ele_names[i], paste(colnames(df_in)[match_i], collapse = " "), sep = " - "))
    
    # Append each match to the output vector
    match_cols <- c(match_cols, match_i)
  }
  return(match_cols)
}

# Import and format data ----

# Import everything
if(TRUE){
  # Import GeoRoc data
  
  # Note: I had each GeoRoc data file saved to a separate data file. This is probably
  #   because it was a more convenient way to use the "GeoRoc Parser" script.
  #   Here, I append a "folder" column so I can recall which precompiled file the
  #   data come from.
  folders <- list.files()[grep("\\.", list.files(), invert = TRUE)]
  
  # Import data from the first folder
  path_in <- paste(folders[1], "Georoc_parsed.csv", sep = "/")
  Import <- read.csv(path_in, stringsAsFactors = FALSE)
  Folder <- rep(folders[1], nrow(Import))
  Import <- cbind(Folder, Import)
  
  # Loop over each of the remaining folders, open the georoc parsed file from each,
  #   and append each next file onto the previous one
  for(i in 2:length(folders)){
    path_in <- paste(folders[i], "Georoc_parsed.csv", sep = "/")
    Import_add <- read.csv(path_in, stringsAsFactors = FALSE)
    Folder <- rep(folders[i], nrow(Import_add))
    Import_add <- cbind(Folder, Import_add)
    Import <- merge_dfs(Import, Import_add)
  }
  
  # Import my melt inclusion data
  MI_vol <- read.csv("G:/My Drive/VT Volcanoes/Hawaii/Hawaii Rejuvenated/Hawaii paper/Hawaii paper 6 -- Reviewer comments/200423 light element plots/Table S5 -- MI volatile elements simplified.csv", stringsAsFactors = FALSE)
  MI_maj <- read.csv("G:/My Drive/VT Volcanoes/Hawaii/Hawaii Rejuvenated/Hawaii paper/Hawaii paper 6 -- Reviewer comments/200423 light element plots/Table S3 -- MI major elements simplified.csv", stringsAsFactors = FALSE)
  MI_trace <- read.csv("G:/My Drive/VT Volcanoes/Hawaii/Hawaii Rejuvenated/Hawaii paper/Hawaii paper 6 -- Reviewer comments/200423 light element plots/Table S4 -- MI trace elements simplified.csv", stringsAsFactors = FALSE, na.strings = c("b.q.l.", "b.d.l."))
  
  # Add major element data to volatile element data
  My_MIs <- cbind(MI_vol, get_chunk(MI_vol, MI_maj))
  
  # Add trace element data to other data
  MI_trace <- MI_trace[-(22:29),]
  My_MIs <- cbind(My_MIs, get_chunk(My_MIs, MI_trace))
  
  # Import whole rock data
  WR_maj <- read.csv("G:/My Drive/VT Volcanoes/Hawaii/Hawaii Rejuvenated/Hawaii paper/Hawaii paper 6 -- Reviewer comments/200423 light element plots/Table S1 -- WR Major elements simplified.csv", stringsAsFactors = FALSE)
  WR_trace <- read.csv("G:/My Drive/VT Volcanoes/Hawaii/Hawaii Rejuvenated/Hawaii paper/Hawaii paper 6 -- Reviewer comments/200423 light element plots/Table S2 -- WR Trace elements simplified.csv", stringsAsFactors = FALSE)
  WR_data <- cbind(WR_maj, get_chunk(WR_maj, WR_trace))
}

# Format Sample, K, Ti, and P in Import file (GeoRoc data)
for(i in 1:nrow(Import)){
  if(is.na(Import$K.PPM.[i]) & !is.na(Import$K2O.WT..[i])){
    Import$K.PPM.[i] <- Import$K2O.WT..[i] / (1.205*(10^-4))
  }
  if(is.na(Import$P.PPM.[i]) & !is.na(Import$P2O5.WT..[i])){
    Import$P.PPM.[i] <- Import$P2O5.WT..[i] / (2.292*(10^-4))
  }
  if(is.na(Import$TI.PPM.[i]) & !is.na(Import$TIO2.WT..[i])){
    Import$TI.PPM.[i] <- Import$TIO2.WT..[i] / (1.668*(10^-4))
  }
}
# Add a "Sample column using the location as an identifier
sample_col <- data.frame(Sample = Import$LOCATION, stringsAsFactors = FALSE)
Import <- cbind(Import, sample_col)
# Add a "flag" column
flag_col <- data.frame(flag = rep(NA, nrow(Import)))
Import <- cbind(flag_col, Import)

# Format K, Ti, and P in MI file
for(i in 1:nrow(My_MIs)){
  My_MIs$K2O..wt.[i] <- My_MIs$K2O..wt.[i] / (1.205*(10^-4))
  My_MIs$P2O5..wt.[i] <- My_MIs$P2O5..wt.[i] / (2.292*(10^-4))
  My_MIs$TiO2..wt.[i] <- My_MIs$TiO2..wt.[i] / (1.668*(10^-4))
}
# Add a "flag" column
flag_col <- data.frame(flag = rep(NA, nrow(My_MIs)))
My_MIs <- cbind(flag_col, My_MIs)

# Format K, Ti, and P in WR data
for(i in 1:nrow(WR_data)){
  WR_data$X.K2O..wt.[i] <- WR_data$X.K2O..wt.[i] / (1.205*(10^-4))
  WR_data$X.P2O5..wt.[i] <- WR_data$X.P2O5..wt.[i] / (2.292*(10^-4))
  WR_data$X.TiO2..wt.[i] <- WR_data$X.TiO2..wt.[i] / (1.668*(10^-4))
}
# Add a "flag" column
flag_col <- data.frame(flag = rep(NA, nrow(WR_data)))
WR_data <- cbind(flag_col, WR_data)

# Generate norm df to use for flagging columns
norm_data <- init_norm()

# Pre-emptively identify norm columns for each dataset
norm_cols_mi <- get_cols(df_in = My_MIs, ele_names = colnames(norm_data))
norm_col_wr <- get_cols(df_in = WR_data, ele_names = colnames(norm_data))
norm_cols_import <- get_cols(df_in = Import, ele_names = colnames(norm_data))


# Flag data and calculate norms ----

load_settings("G:/My Drive/VT Volcanoes/Hawaii/Hawaii Rejuvenated/Hawaii paper/Hawaii paper 6 -- Reviewer comments/Figures in prep/Element norm.txt")
norm_data <- init_norm()

# Add MI data 
flag_sub <- grep("HA-16-6", My_MIs$Sample)
norm_data <- rbind(norm_data, calc_norm(norm_data, My_MIs[flag_sub,], flag = "HA-16-6"
                  , norm_cols = norm_cols_mi))
flag_sub <- grep("HA-16-7", My_MIs$Sample)
norm_data <- rbind(norm_data, calc_norm(norm_data, My_MIs[flag_sub,], flag = "HA-16-7"
                  , norm_cols = norm_cols_mi))

# Add GeoRoc data
flag_sub <- intersect(grep("HALEAKALA", Import$LOCATION)
                  , grep("/ SHIELD STAGE", Import$LOCATION)
                  )
norm_data <- rbind(norm_data, calc_norm(norm_data, Import[flag_sub,], flag = "HA-shield"
                  , norm_cols = norm_cols_import))

flag_sub <- intersect(grep("HALEAKALA", Import$LOCATION)
                      , union(grep("/ POST-SHIELD STAGE", Import$LOCATION)
                              , grep("/ REJUVENATED STAGE", Import$LOCATION)
                                  )
                      )
norm_data <- rbind(norm_data, calc_norm(norm_data, Import[flag_sub,], flag = "HA-post-shield"
                                        , norm_cols = norm_cols_import))

flag_sub <- grep("HAVF", Import$Folder)
norm_data <- rbind(norm_data, calc_norm(norm_data, Import[flag_sub,], flag = "HAVF"
                                        , norm_cols = norm_cols_import))

# Add WR data
flag_sub <- grep("HA", WR_data$Sample)
norm_data <- rbind(norm_data, calc_norm(norm_data, WR_data[flag_sub,], flag = "WR_HA"
                  , norm_cols = norm_col_wr))
flag_sub <- grep("MA", WR_data$Sample)
norm_data <- rbind(norm_data, calc_norm(norm_data, WR_data[flag_sub,], flag = "WR_MA"
                                        , norm_cols = norm_col_wr))

# Plot figures ----

spider_draw(norm_data, ylim = c(4, 130))

spider_poly(norm_data, "WR_MA",          col = rgb(1, 0.7, 0.3, 0.3))
spider_poly(norm_data, "WR_HA",          col = rgb(1, 1, 0.1, 0.3))
spider_poly(norm_data, "HA-shield",      col = rgb(0.2, 1, 0.3, 1)  , filt_quantile = c(0.25, 0.75))
spider_poly(norm_data, "HA-post-shield", col = rgb(0.2, 0.3, 1, 1)  , filt_quantile = c(0.25, 0.75))
spider_poly(norm_data, "HAVF",           col = rgb(0.8, 0.2, 0.8, 1), filt_quantile = c(0.25, 0.75))

spider_lines(norm_data, "HA-16-6", col = rgb(0.8, 0, 0), lwd = 2)
spider_lines(norm_data, "HA-16-7", col = rgb(0, 0.8, 0), lwd = 2)

