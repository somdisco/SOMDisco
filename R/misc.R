

# Decode fille colors
.decode_som_fill = function(SOM, fill, ctab = NULL) {
  
  ## Case 1:  No fill color given 
  ## Return a vector NA (no filling)
  if(is.null(fill)) {
    fill = rep(NA, SOM$nW)
    return(fill)
  }
  
  ## *** Fill needs to be either length 1 or length = nW
  if(!(length(fill)==1) && !(length(fill)==SOM$nW)) stop("length(fill) must = 1 or nW")
  
  ## Case 2: If a ctab is given, we assume that fill contains corresponding labels 
  ## whose colors should be retrieving by matching to the color table 
  if(!is.null(ctab)) {
    # Check the ctab 
    SOMDisco:::.ctab_check(ctab = ctab, query_labels = fill)
    # Assign the colors 
    fill = SOMDisco:::.ctab_query(ctab = ctab, query_labels = fill)
    # If fill is length=1, assume that the color should be recycled 
    if(length(fill)==1) fill = rep(fill, SOM$nW)
    return(fill)
  }
  
  ## *** From here on, we can assume fill is intended to contain color names 
  
  ## Case 3:  fill contains color names
  fill_is_color_name = SOMDisco:::.is_valid_color(na.omit(fill))
  if(!all(fill_is_color_name)) stop("If fill contains color names they must all be valid color names.")
  if(length(fill)==1) fill = rep(fill, SOM$nW)
  return(fill)
}


.vis_pie_slice = function(start.degree = 0, end.degree = 360, radius = 1, center = c(0, 0), clock.wise = TRUE, 
                         col = NA,border = "black",lwd = par("lwd"),lty = par("lty")) {
  
  ## *** Internal Functions 
  polar2Cartesian = function(d) {
    theta = as.radian(d[, 1])
    rou = d[, 2]
    x = rou * cos(theta)
    y = rou * sin(theta)
    return(cbind(x, y))
  }
  
  as.radian = function(degree) {
    return(degree/180*pi)
  }
  
  is.circular = function(start.degree, end.degree) {
    (end.degree - start.degree) %% 360 == 0 && (end.degree - start.degree) != 0
  }
  
  degree_diff = function(start, end, clock.wise = TRUE) {
    if(is.circular(start, end)) {
      360
    } else {
      start = start %% 360
      end = end %% 360
      if(clock.wise) (start - end) %% 360
      else (end - start) %% 360
    }
  }
  
  # from start to end
  degree_seq = function(start, end, clock.wise = TRUE, ...) {
    if(is.circular(start, end)) {
      seq(0, 360, ...)
    } else {
      start = start %% 360
      end = end %% 360
      if(clock.wise) {
        # make start is larger than end, but the difference is less than 360
        if(start < end) start = start + 360
        seq(start, end, ...)
      } else {
        if(start > end) start = start - 360
        seq(start, end, ...)
      }
    }
  }
  
  d1 = NULL
  
  # calculate the number of segments of the up arc
  unit.circle.segments = 500
  l1 = as.radian(degree_diff(start.degree, end.degree, clock.wise)) * radius
  ncut1 = l1/ (2*pi/unit.circle.segments)
  ncut1 = floor(ncut1)
  ncut1 = ifelse(ncut1 < 2, 2, ncut1)
  
  # d1 is from the start.degree to end.degree
  d1 = rbind(d1, cbind(degree_seq(start.degree, end.degree, clock.wise, length.out = ncut1), rep(radius, ncut1)))
  
  m1 = polar2Cartesian(d1)
  if(is.circular(start.degree, end.degree)) {  # it is a circle
    m = m1
  } else {
    m = rbind(m1, c(0, 0))
  }
  
  # and shift to the center
  m[, 1] = m[, 1] + center[1]
  m[, 2] = m[, 2] + center[2]
  polygon(m, col = col, border = border, lwd = lwd, lty = lty)
  
  return(invisible(NULL))
}


## Lightness value from an RGB code 
.HSL_Lightness = function(R=NULL,G=NULL,B=NULL,hex = NULL) {
  if(!is.null(hex)) {
    RGB = col2rgb(hex)
    R = RGB[1,]
    G = RGB[2,]
    B = RGB[3,]
    hex_is_na = is.na(hex)
    R[hex_is_na] = G[hex_is_na] = B[hex_is_na] = NA 
  } else {
    if(is.null(R) && is.null(G) && is.null(B)) stop("either R+G+B, or hex, must be given")
  }
  return((pmin(R,G,B) + pmax(R,G,B))/(255*2))
}




# Map a vector of values to [0,1] via linear scaling 
# 
# @param x the values to map 
# @param clamp the type of clamping applied to x before mapping. 
# Options are "none" (default), meaning no clamping is used; 
# "identity", meaning the clamp.range is used directly; 
# "quantile", meaning the clamp.range quantiles are used.  
# @param clamp.range sets the clamping range. 
# If clamp = 'none', the min/max of values is used. 
# If clamp = 'identity', clamp.range should be given as an effective [lo,hi] range to clamp values to. 
# If clamp = 'quantile', clamp.range should contain the effective [lo,hi] probabilities whose quantiles define the clamp range. 
# Default = c(0,1), which should be changed if clamp = 'identity'
# 
# @return the mapped value 
.unitmap_linscale = function(x, clamp = "none", clamp.range = c(0.0, 1.0)) {
  # clamp can = "none", "identity", or "quantile" 
  
  # No clamping means use the min/max of x as scaling range
  if(clamp == "none") {
    clamp_lo = min(x, na.rm = T); 
    clamp_hi = max(x, na.rm = T); 
    
    # Identity clamping means use the given lo / hi values as scaling range  
  } else if(clamp == "identity") {
    
    # Quantile clamping means use the quantiles at the given lo / hi values as scaling range 
  } else if(clamp == "quantile") {
    clamp_lo = quantile(x, probs = clamp.range[1], na.rm = T); 
    clamp_hi = quantile(x, probs = clamp.range[2], na.rm = T); 
  } else {
    stop("clamp must be one of 'none', 'identity', or 'quantile'");
  }
  
  # Initialize the output vector as the x values, and clamp them to the given range 
  u = x; 
  if(clamp != "none") {
    u[u < clamp_lo] = clamp_lo; 
    u[u > clamp_hi] = clamp_hi; 
  }
  
  # Linearly scale to [0,1] using the given ranges 
  u = (u - clamp_lo) / (clamp_hi - clamp_lo);
  
  return(u)
}


# Map a vector of values to [0,1] via CDF of its standardized values 
# 
# @param x the values to map 
# @param clamp the type of clamping applied to x before mapping. 
# Options are "none" (default), meaning no clamping is used; 
# "identity", meaning the clamp.range is used directly; 
# "quantile", meaning the clamp.range quantiles are used.  
# @param clamp.range sets the clamping range. 
# If clamp = 'none', the min/max of values is used. 
# If clamp = 'identity', clamp.range should be given as an effective [lo,hi] range to clamp values to. 
# If clamp = 'quantile', clamp.range should contain the effective [lo,hi] probabilities whose quantiles define the clamp range. 
# Default = c(0,1), which should be changed if clamp = 'identity'
# 
# @return the mapped value 
.unitmap_zcdf = function(x, clamp = "none", clamp.range = c(0.0, 1.0)) {
  # clamp can = "none", "identity", or "quantile" 
  
  # No clamping means use the min/max of x as scaling range
  if(clamp == "none") {
    clamp_lo = min(x, na.rm = T); 
    clamp_hi = max(x, na.rm = T); 
    
    # Identity clamping means use the given lo / hi values as scaling range  
  } else if(clamp == "identity") {
    
    # Quantile clamping means use the quantiles at the given lo / hi values as scaling range 
  } else if(clamp == "quantile") {
    clamp_lo = quantile(x, probs = clamp.range[1], na.rm = T); 
    clamp_hi = quantile(x, probs = clamp.range[2], na.rm = T); 
  } else {
    stop("clamp must be one of 'none', 'identity', or 'quantile'");
  }
  
  # Initialize the output vector as the x values, and clamp them to the given range 
  u = x; 
  if(clamp != "none") {
    u[u < clamp_lo] = clamp_lo; 
    u[u > clamp_hi] = clamp_hi; 
  }
  
  u = scale(x, center = T, scale = T)
  
  u = pnorm(u)
  
  return(u)
}


# process_color.palette = function(color.palette) {
#   # input color.palette is either a string of form <palette>.<name>, or a vector of color names 
#   # output is a vector of color names which comprise the requested palette 
#   
#   ## *** Decode color palette 
#   color.palette = tolower(color.palette)
#   
#   ## Check for Viridis 
#   if(any(stringr::str_detect(color.palette, "viridis"))) {
#     color.palette = stringr::str_replace(color.palette, "viridis.", "")
#     
#     if(color.palette == "plasma") color.palette = viridis::plasma(n = 50)
#     else if(color.palette == "magma") color.palette = viridis::magma(n = 50)
#     else if(color.palette == "inferno") color.palette = viridis::inferno(n = 50)
#     else if(color.palette == "viridis") color.palette = viridis::viridis(n = 50) 
#     else if(color.palette == "cividis") color.palette = viridis::cividis(n = 50) 
#     else stop("Unknown viridis palette name. Check ?viridis")
#     
#     
#     ## Check for Brewer 
#   } else if(any(stringr::str_detect(color.palette, "brewer"))) {
#     color.palette = stringr::str_replace(color.palette, "brewer.", "")
#     
#     which_brewer_pal = match(color.palette, tolower(rownames(RColorBrewer::brewer.pal.info)))    
#     if(is.na(which_brewer_pal)) stop("Unknown brewer palette name. Check ?RColorBrewer")
#     
#     color.palette = RColorBrewer::brewer.pal(n = RColorBrewer::brewer.pal.info$maxcolors[which_brewer_pal], 
#                                              name = rownames(RColorBrewer::brewer.pal.info)[which_brewer_pal])
#     
#     ## Check for given color names   
#   } else if(!all(SOMDisco:::.is_valid_color(color.palette))) {
#     stop("Unknown color names in palette. Check.")
#   } 
#   
#   return(color.palette)
# }


#' View a color table
#'
#' @description This is a helper function to provide a quick visualization of the labels and their associated colors that are defined
#' in an input color table.
#'
#' @param ctab the color table dataframe to view. Must have fields  \code{$class}, \code{$R}, \code{$G}, \code{$B}.
#' @param label.cex the size of labels, as plotted on the visual color grid
#' @param label.font the font size of the label, 1=regular, 2 = bold
#' @param nrows_in_display sets the number of rows in the visualized color table. Default = NULL invokes automatic generation of this value,
#' with an attempt to make the resulting color grid as square as possible.
#' @return nothing, only used for visualization
#' @export
vis_ctab = function(ctab, label.cex = 0.9, label.font = 2, nrows_in_display = NULL) {
  
  ## Check the ctab for validity 
  SOMDisco:::.ctab_check(ctab = ctab)
  
  ### Dimensions of display table
  n = nrow(ctab)
  
  if(is.null(nrows_in_display)) {
    xdim = floor(sqrt(n))
    ydim = ceiling(n / xdim)
  } else {
    xdim = ceiling(n / nrows_in_display)
    ydim = ceiling(n / xdim)
  }
  
  ### x & y coordinates of center squares
  xcoord = 1:xdim - 0.5
  ycoord = ydim:1 - 0.5
  xy = as.matrix(expand.grid(xcoord, ycoord))
  xy = xy[1:n,]
  
  ### Colors of blocks from color table
  #block.colors = grDevices::rgb(red = ctab$R, green = ctab$G, blue = ctab$B, maxColorValue = 255)
  block.colors = ctab$color
  
  ### The text color of the labels is dependent on the corresponding block color
  ### Default text color is black, but if the block is too dark (toward black)
  ### the text should be white to be readable.
  text.colors = rep("black", times = n)
  Lightness = .HSL_Lightness(hex = ctab$color)
  text.colors[Lightness < 0.4] = "white" ## use default 0.4 as a cutoff if not given
  
  graphics::par(mar = rep(0,4), oma = rep(0,4), bty="n", xaxs = "i", yaxs = "i")
  plot(x = c(0, xdim), y = c(0, ydim), type="n", asp = 1, axes = F)
  graphics::rect(xleft = xy[,1]-0.5, ybottom = xy[,2] - 0.5, xright = xy[,1] + 0.5, ytop = xy[,2] + 0.5,
                 col = block.colors, border = NA)
  graphics::text(x = xy[,1], y = xy[,2], labels = ctab$label, cex = label.cex , offset = 0, col = text.colors, font = label.font)
  
}

