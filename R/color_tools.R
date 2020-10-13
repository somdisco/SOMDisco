#' Tableau 20 Color Palette 
#' @return a named vector containing Tableau 20 colors 
tableau20 = function(col = NULL) {
  mycolors = c("blue" = "#1F77B4", "lightblue" = "#73C2FB", 
               "orange" = "#ff7f0e", "lightorange" = "#ffd27f", 
               "green" = "#2ca02c", "lightgreen" = "#98fb98", 
               "red" = "#d62728", "lightred" = "#ffb09c", 
               "purple" = "#b660cd", "lightpurple" = "#e4a0f7", 
               "yellow" = "#ffdb58", "lightyellow" = "#fdfd96", 
               "teal" = "#17becf", "lightteal" = "#c8ffff", 
               "gray" = "#88807b", "lightgray" = "#c7c6c1", 
               "brown" = "#8c564b", "lightbrown" = "#ceb180", 
               "pink" = "#ff6fff", "lightpink" = "#fde6fa")
  
  if(is.null(col)) {
    return(mycolors)
  } else {
    return(mycolors[col])
  }
}


# Determine whether a color name is valid 
# 
# @param color_name (possibly a vector of) color names to test 
# 
# @return T if color_name can be converted to RGB codes, F otherwise
.is_valid_color = function(color_list) {
  
  
  # out = sapply(color_list, function(X) {
  #   tryCatch(is.matrix(col2rgb(X)), 
  #            error = function(e) FALSE)
  # })
  # 
  # out[is.na(color_names)] = NA
  # return(out)
  
  
  ## First, check if color_list contains HEX codes 
  if(is.character(color_list)) {
    
    # Determine whether each character string is a valid HEX color 
    valid_color = sapply(color_list, function(X) {tryCatch(is.matrix(col2rgb(X)), error = function(e) FALSE)})
    
    valid_color[is.na(color_list)] = NA 
    
    return(valid_color)
  } 
  
  
  ## Next, check if color_list contains RGB triplets 
  if(is.matrix(color_list)) {
    
    # RGB triplets must be in rows of matrix, with 3 columns 
    if(ncol(color_list) != 3) return(FALSE)
    
    return(TRUE)
  } 
}



.distinguishable_colors = function(ncolors, palette_seed = NULL, Lightness_range = c(0.3, 0.8), return_as = "HEX") {
  ## Returns a set of maximally distinguishable colors in L*ab colorspace 
  ## Function translated to R from Matlab. 
  ## Source:  https://www.mathworks.com/matlabcentral/fileexchange/29702-generate-maximally-perceptually-distinct-colors
  
  ## Inputs: 
  ## ncolors = # of colors requested 
  ## palette_seed = a list of colors used to initialize the distinguishable set. 
  ##   If NULL the tableau20 palette is used. 
  ##   Can be either a vector of HEX codes, or a matrix whose rows contain RGB codes. 
  ##   If given as RGB, assumed to be in [0-255]. 
  ## Lightness_range = range of Ligthness values (after converting to L*ab colorspace) for resulting 
  ##   colors to occupy. Candidate colors are select across a grid of RGB space and converted to L*ab. 
  ##   Those whose Lightness value is outside this range are removed prior to the selection process 
  ## return_as = one of {"HEX","RGB"}
  
  
  ## *** Process the starting palette 
  ## This produces a palette_seed matrix in RGB format 
  if(is.null(palette_seed)) {
    palette_seed = c("#1F77B4", # blue 
                     "#ff7f0e", # orange 
                     "#2ca02c", # green 
                     "#d62728", # red 
                     "#b660cd", # purple 
                     "#ffdb58", # yellow 
                     "#17becf", # teal 
                     "#88807b", # gray 
                     "#8c564b", # brown 
                     "#ff6fff", # pink 
                     "#73C2FB", # light blue 
                     "#ffd27f", # light orange
                     "#98fb98", # light green 
                     "#ffb09c", # light red 
                     "#e4a0f7", # light purple 
                     "#fdfd96", # light yellow 
                     "#c8ffff", # light teal
                     "#c7c6c1", # light gray 
                     "#ceb180", # light brown 
                     "#fde6fa") # light pink 
    # Convert to RGB 
    palette_seed = unname(t(col2rgb(palette_seed)))
    
  } else {
    
    # Check for HEX codes 
    if(is.character(palette_seed)) {
      
      # Make sure colors are valid 
      if(!all(na.omit(.is_valid_color(color_list = palette_seed)))) stop("palette_seed contains invalid color names.")
      
      # Convert to RGB 
      palette_seed = unname(t(col2rgb(palette_seed)))
      
    # Check for RGB triplets 
    } else {
      if(!(is.matrix(palette_seed) && ncol(palette_seed)==3)) {
        stop("Input palette_seed in unknown format.")
      }
    }
  }
  
  
  ## *** Return palette_seed, if ncolors < nrow(palette_seed)
  if(ncolors <= nrow(palette_seed)) {
    out = rgb(palette_seed[1:ncolors,,drop=F], maxColorValue = 255)
    return(out)
  }
  
  
  ## *** Otherwise, create a color dictionary in RGB, 
  ## Convert to L*ab space & filter by given lightness values 
  grid_size = 50
  grid_points = seq(0, 255, length.out = grid_size)
  color.dict = unname(as.matrix(expand.grid(grid_points, grid_points, grid_points)))
  color.dict.lab = grDevices::convertColor(color.dict/255, from = "sRGB", to = "Lab")
  remove.these = which(color.dict.lab[,1] < Lightness_range[1]*100 | color.dict.lab[,1] > Lightness_range[2]*100)
  if(length(remove.these) > 0) {
    color.dict = color.dict[-remove.these,,drop=F]
    color.dict.lab = color.dict.lab[-remove.these,,drop=F]
  }
  
  ## Convert the seed palette to L*ab and append it to the diciionaries
  ## Do this after filtering by Lightness above to avoid removing any user-specified colors 
  color.dict = rbind(palette_seed, color.dict)
  palette_seed.lab = grDevices::convertColor(palette_seed/255, from = "sRGB", to = "Lab")
  color.dict.lab = rbind(palette_seed.lab, color.dict.lab)
  
  
  ## *** Compute distances from the dictionary colors to the seed colors
  mindist2 = rep(Inf, nrow(color.dict))
  for(i in 1:nrow(palette_seed)) {
    # Dist of all colors from this seed color 
    dist2 = rowSums(sweep(color.dict.lab, MARGIN = 2, palette_seed.lab[i,])^2)
    # Dist to closest seed color 
    mindist2 = pmin(dist2, mindist2)
  } 
  
  
  ## *** Initialize the returned colors to the seeds
  out.color.idx = 1:nrow(palette_seed)
  ncolors.to.choose = ncolors - nrow(palette_seed)
  
  
  ## *** Iteratively pick the next color from the dictionary 
  ## that maximizes the distance to the nearest already-picked color
  prev.lab = palette_seed.lab[nrow(palette_seed.lab),] # initialize by making the "previous" color equal to last seed color 
  for(i in 1:ncolors.to.choose) {
    # Dist of all colors from the previous color 
    dist2 = rowSums(sweep(color.dict.lab, MARGIN = 2, STATS = prev.lab)^2) 
    # Dist to closest oreviously chosen color 
    mindist2 = pmin(dist2, mindist2)
    # Index in color.dict of farthest color from all previously-chosen colors
    index = which.max(mindist2)
    # Store 
    out.color.idx = c(out.color.idx, index)
    # Reset previously chosen color 
    prev.lab = color.dict.lab[index,]
  }
  
  
  ## *** Return colors, in their selection order 
  out.colors = color.dict[out.color.idx,,drop=F]
  if(return_as == "RGB") return(out.colors)
  
  return(rgb(out.colors, maxColorValue = 255))
}



#' @name build_ctab 
#' @title Build a color table
#' @description Returns a default \code{ctab} given a set of (possibly non-unique) labels. 
#' The returned object is a data frame with columns \code{label} and \code{color} defining the color mapping. 
#' @param labels a set of data labels to whose unique values colors will be assigned
#' @export 
build_ctab = function(labels) {
  
  ## Get & sort unique set of input labels 
  unq.labels = gtools::mixedsort(unique(na.omit(labels)))
  
  ctab = data.frame(label = unq.labels, 
                    color = .distinguishable_colors(ncolors = length(unq.labels), Lightness_range = c(0.3, 0.9), return_as = 'HEX'), 
                    stringsAsFactors = F)
  
  return(ctab)
}


# Check a color table for validity 
# Requirements for valid ctab: 
#   1. should be a data frame with columns 'label' and 'color'
#   2. 'color' should contain valid HEX color codes 
#   3. If a set of query labels is given, the ctab will be checked to ensure it contains every label in the list 
#      NA labels are ignored in the check (i.e., they do not trigger an error)
# If any of the above conditions are not met, this function will error, otherwise it returns TRUE 
.ctab_check = function(ctab, query_labels = NULL) {
  
  # Check that ctab is a data frame with columns 'label' and 'color' 
  if(class(ctab)!="data.frame") {
    stop("ctab must be a data frame")
    return(FALSE)
  } 
  
  # Check that ctab has correct columns 
  if(!all(c('label','color') %in% names(ctab))) {
    stop("ctab must contain columns 'label' and 'color'")  
    return(FALSE)
  } 
  
  # Check that all the colors are valid 
  color_check = SOMDisco:::.is_valid_color(ctab$color)
  if(any(!color_check)) {
    bad_colors = names(which(!color_check))
    bad_colors = paste(bad_colors, collapse = ",")
    stop(paste("ctab$color contains invalid color names",bad_colors))
    return(FALSE)
  }
  
  ## Match labels to color table, if labels were given 
  if(!is.null(query_labels)) {
    
    # Remove NAs 
    query_labels = gtools::mixedsort(unique(na.omit(query_labels)))
    
    # Match query labels to ctab$label
    match_idx = match(query_labels, ctab$label)
    
    # Return error if any query_labels are unmatched in ctab 
    missing_query_labels = query_labels[which(is.na(match_idx))]
    if(length(missing_query_labels) > 0) {
      missing_query_labels = paste(missing_query_labels, collapse = ",")
      stop(paste("ctab is missing colors for labels",missing_query_labels))
      return(FALSE)
    }
  }
  
  return(TRUE)
}

# ctab_check = function(ctab, query_labels = NULL) {
#   
#   # Check that ctab is a data frame with columns 'label' and 'color' 
#   if(class(ctab)!="data.frame") {
#     stop("ctab must be a data frame")
#     return(FALSE)
#   } 
#   
#   # Check that ctab has correct columns 
#   if(!all(c('label','color') %in% names(ctab))) {
#     stop("ctab must contain columns 'label' and 'color'")  
#     return(FALSE)
#   } 
#   
#   # Check that all the colors are valid 
#   color_check = SOMDisco:::.is_valid_color(ctab$color)
#   if(any(!color_check)) {
#     bad_colors = names(which(!color_check))
#     bad_colors = paste(bad_colors, collapse = ",")
#     stop(paste("ctab$color contains invalid color names",bad_colors))
#     return(FALSE)
#   }
#   
#   ## Match labels to color table, if labels were given 
#   if(!is.null(query_labels)) {
#     
#     # Remove NAs 
#     query_labels = gtools::mixedsort(unique(na.omit(query_labels)))
#     
#     # Match query labels to ctab$label
#     match_idx = match(query_labels, ctab$label)
#     
#     # Return error if any query_labels are unmatched in ctab 
#     missing_query_labels = query_labels[which(is.na(match_idx))]
#     if(length(missing_query_labels) > 0) {
#       missing_query_labels = paste(missing_query_labels, collapse = ",")
#       stop(paste("ctab is missing colors for labels",missing_query_labels))
#       return(FALSE)
#     }
#   }
#   
#   return(TRUE)
# }


# Query a color table given a list of labels
# Return the list of colors mapping the query labels 
# This function assumes that the ctab has already been checked for validity via .ctab_check
.ctab_query = function(ctab, query_labels) {
  return(ctab$color[match(query_labels, ctab$label)])
}
