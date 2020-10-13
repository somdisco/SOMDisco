

#' Setup the SOM lattice for visualizations
#' 
#' @param SOM a SOM object
#' @param mar Optional, the plot margin around the lattice (will be recycled for all sides). 
#' Default = 0.1
#' @param lattice_coords whether to add the lattice (i,j) coordinates to the left & bottom of the plot 
#' @param coords.cex size of lattice coordinate text, if requested
#' @param active Optional, if the SOM object has been recalled,
#' restricts plotting only to active neurons (those whose RF_size > 0).
#' Default = TRUE.
#' @param subset Optional, a vector of neuron indices to restrict the plotting to. 
#' Default = NULL imposes no restriction (plots whole lattice)
#' @param change.par whether to allow the vis function to optimally change the \code{par} of the plot. Default = TRUE. 
#' If \code{par} is allowed to changed internally, it is always reset upon function exit. 
#' @export
vis_som_setup = function(SOM, mar = 0.1,  lattice_coords = F, coords.cex = 0.75, active = T, subset = NULL, change.par = TRUE) {
  if(!class(SOM)=="Rcpp_SOM") stop("Input SOM must be an object of class SOM")
  
  
  plot_these = rep(T, SOM$nW)
  if(active && SOM$is_recalled) {
    plot_these[SOM$RF_size==0] = FALSE
  }
  if(!is.null(subset)) {
    plot_these[setdiff(1:SOM$nW, subset)] = F
  }
  
  ## Compute and store the plotting limits 
  vert_xrng = range(apply(SOM$nu_verts[,,plot_these,drop=F],3,function(z) range(z[,1])))
  vert_yrng = range(apply(SOM$nu_verts[,,plot_these,drop=F],3,function(z) range(z[,2])))
  xlim = vert_xrng
  ylim = vert_yrng 
  if(lattice_coords) {
    xlim[1] = xlim[1] - 0.75
    ylim[1] = ylim[1] - 0.75
  }
  
  # Save old par(), define the som_vis par, and save it
  if(change.par) {
    opar = par(no.readonly=TRUE)
    on.exit(par(opar))
    par(mar = rep(mar, 4), xpd = NA, pty = 'm', xaxs = 'i', yaxs='i', xaxt = 'n', yaxt = 'n')  
  }
  
  plot(0, cex = 0, frame.plot = F, axes = F, xlab=NA, ylab=NA, asp = 1, xlim=xlim, ylim=ylim)
 
  if(lattice_coords) {
    coords.i = sort(unique(SOM$nu_ij[plot_these,1]))
    coords.j = sort(unique(SOM$nu_ij[plot_these,2]))
    coords.x = coords.y = NULL 
    for(i in coords.i) {
      coords.y = c(coords.y, median(SOM$nu_xy[SOM$nu_ij[,1]==i,2]))
    }
    for(j in coords.j) {
      coords.x = c(coords.x, median(SOM$nu_xy[SOM$nu_ij[,2]==j,1]))
    }
    
    text(x = coords.x, y = rep(ylim[1]+0.25, length(coords.x)), labels = coords.j, cex = coords.cex)
    text(x = rep(xlim[1]+0.25, length(coords.y)), y = coords.y, labels = coords.i, cex = coords.cex)
  }
  
  SOM$set_vis_par(par(no.readonly = T))
  SOM$set_vis_tile_bg(rep("white", SOM$nW))
  SOM$set_vis_xlim(xlim)
  SOM$set_vis_ylim(ylim)
  
  #.SOMDisco_vis_som_par <<- list()
  
  #.SOMDisco_vis_som_par$par <<- par(no.readonly = T)
  #.SOMDisco_vis_som_par$tile_fill <<- rep('white', SOM$nW)
  #.SOMDisco_vis_som_par$xlim <<- xlim
  #.SOMDisco_vis_som_par$ylim <<- ylim
}



#' Visualize SOM neurons on the lattice
#' 
#' @param SOM a SOM object
#' @param add whether to create a new SOM visualization panel (=FALSE, default), or add to an existing one (=TRUE)
#' @param nu.pch neuron marker type, default = 16
#' @param nu.cex neuron marker size, default = 1,
#' @param nu.col neuron marker color, default = "black"
#' @param active Optional, if the SOM object has been recalled,
#' restricts plotting only to active neurons (those whose RF_size > 0).
#' Default = TRUE.
#' @param subset Optional, a vector of neuron indices to restrict the plotting to. 
#' Default = NULL imposes no restriction (plots whole lattice)
#' @param change.par whether to allow the vis function to optimally change the \code{par} of the plot. Default = TRUE. 
#' If \code{par} is allowed to changed internally, it is always reset upon function exit. 
#' 
#' @details See ?points or ?text for help with plotting parameters cex, col, font, etc 
#' @export
vis_som_neurons = function(SOM, add = FALSE, nu.pch = 16, nu.cex = 1, nu.col = "black", 
                           active = TRUE, subset = NULL, change.par = TRUE) {
  if(!class(SOM)=="Rcpp_SOM") stop("Input SOM must be an object of class SOM")
  
  ## Determine which neurons to plot 
  plot_these = rep(TRUE,SOM$nW)
  if(active && SOM$is_recalled) {
    plot_these[SOM$RF_size==0] = FALSE
  }
  if(!is.null(subset)) {
    plot_these[setdiff(1:SOM$nW, subset)] = FALSE
  }
  
  ## Decode their color 
  if(length(nu.col) == 1) {
    nu.col = rep(nu.col, SOM$nW)
  } else if(length(nu.col) != SOM$nW) {
    stop("If nu.col is given as a vector it must have length = SOM$nW")
  }
  
  ## Setup plot, if not adding 
  if(!add) vis_som_setup(SOM = SOM, subset = subset, change.par = change.par)
  
  ## Plot neuron markers  
  if(change.par) {
    opar = par(no.readonly = T)
    on.exit(par(opar))
    par(SOM$vis_par)
  }
  points(x = SOM$nu_xy[plot_these,1], y = SOM$nu_xy[plot_these,2], pch = nu.pch, cex = nu.cex, col = nu.col[plot_these])  
}



#' Visualize SOM lattice tiles 
#' 
#' @param SOM a SOM object
#' @param add whether to create a new SOM visualization panel (=FALSE, default), or add to an existing one (=TRUE)
#' @param border.lwd tile border line width, default = 1. Set = 0 to suppress plotting borders. 
#' @param border.lty tile border line type, default = 1 
#' @param border.col tile border color, default = "black". 
#' If border.col = 'fill' and input fill is given, borders will inherit fill colors. 
#' @param fill Optional, controls fill color of plotted tiles. 
#' If given as a single color name (e.g., "black"), all tiles will be colored with this color. 
#' If given as a vector of color names with length = SOM$nW, each tile will be colored separately. 
#' Default = NULL means no fill 
#' @param active Optional, if the SOM object has been recalled,
#' restricts plotting only to active neurons (those whose RF_size > 0).
#' Default = TRUE.
#' @param subset Optional, a vector of neuron indices to restrict the plotting to. 
#' Default = NULL imposes no restriction (plots whole lattice)
#' @param change.par whether to allow the vis function to optimally change the \code{par} of the plot. Default = TRUE. 
#' If \code{par} is allowed to changed internally, it is always reset upon function exit. 
#' 
#' @export
vis_som_tiles = function(SOM, add = FALSE, border.lwd = 1, border.lty = 1, border.col = "black", 
                         fill = NULL, active = TRUE, subset = NULL, change.par = TRUE) {
  if(!class(SOM)=="Rcpp_SOM") stop("Input SOM must be an object of class SOM")
  
  ## Determine which neurons to plot 
  plot_these = rep(TRUE,SOM$nW)
  if(active && SOM$is_recalled) {
    plot_these[SOM$RF_size==0] = FALSE
  }
  if(!is.null(subset)) {
    plot_these[setdiff(1:SOM$nW, subset)] = FALSE
  }
  
  ## Setup plot, if not adding 
  if(!add) vis_som_setup(SOM, subset = subset, change.par = change.par)
  
  ## Decode fill colors 
  fill = SOMDisco:::.decode_som_fill(SOM, fill)
  tmpfill = fill; tmpfill[is.na(fill)] = "white"
  SOM$set_vis_tile_bg(tmpfill)

  
  ## Check if borders should inherit fill colors 
  if(border.col == "fill") {
    border.col = fill 
  } else {
    border.col = rep(border.col, SOM$nW)
  }
  
  ## Decode line width, in case no borders were requested 
  if(border.lwd == 0) {
    border.col = rep(NA, SOM$nW) # polygon border must be NA to suppress plotting
    border.lwd = 1 # polygon line width has to be > 0, otherwise it errors out 
  } 
  
    

  ## Plot each tile separately 
  if(change.par) {
    opar = par(no.readonly = T)
    on.exit(par(opar))
    par(SOM$vis_par)
  }
  
  for(i in 1:SOM$nW) {
    if(!plot_these[i]) next() 
    polygon(x = SOM$nu_verts[,1,i], y = SOM$nu_verts[,2,i], density = NULL, lwd = border.lwd, lty = border.lty, 
            border = border.col[i], col = fill[i])
  }
  

}



#' Visualize SOM fences
#'
#' @param SOM a SOM object
#' @param add whether to create a new SOM visualization panel (=FALSE, default), or add to an existing one (=TRUE)
#' @param fence.lwd fence line width
#' @param clamp controls how values are clamped before unit mapping.  Options are: 
#' "none" (default), meaning no clamping is used; 
#' "identity", meaning the clamp.range is used directly; 
#' "quantile", meaning the clamp.range quantiles are used.  
#' @param clamp.range sets the clamping range. 
#' If clamp = 'none', the min/max of values is used. 
#' If clamp = 'identity', clamp.range should be given as an effective [lo,hi] range to clamp values to. 
#' If clamp = 'quantile', clamp.range should contain the effective [lo,hi] probabilities whose quantiles define the clamp range. 
#' Default = c(0,1), which should be changed if clamp = 'identity'. 
#' @param color.mapping the mapping function used to map values to the range [0,1], which is needed for plotting. 
#' Default ='linear' for a linear mapping. 
#' Can also be 'zcdf', which first computes the z-score of the values, then maps the z-score to [0,1] via the Std. Normal CDF. 
#' @param color.palette the color palette used to represent the fence values. 
#' This should be given as a character vector of color names which will set the colors which represent the 
#' the min (first element) and max (last element) of the fence values. Colors for intermediate values between min and max 
#' will be interpolated via \link[grDevices]{colorRamp}.  
#' Should have at least two elements, but can have more. In this case, the colors will be interpolated throughout the 
#' range of those found in \code{color.palette}.  Example:  \code{c('black','yellow','white')} will map the min values to black, 
#' the mid values to yellow and the max values to white.  
#' Default = c('black','white'). 
#' @param color.range a range in [0,1] defining what portion of the color scale is visualized. 
#' This is useful if you want to restrict the plotted color range to the darker or lighter end of the spectrum. 
#' Default = c(0,1)
#' @param color.nbins allows binning of the (unit-mapped) values prior to color assignment. Options are: 
#' 'none' (default), which just maps the values directly to the color range without binning; 
#' 'auto', which computes a histogram of the (unit-mapped) values and uses the resulting bin counts; 
#' or an integer giving the desired number of bins. 
#' @param active Optional, if the SOM object has been recalled,
#' restricts plotting only to active neurons (those whose RF_size > 0).
#' Default = TRUE.
#' @param subset Optional, a vector of neuron indices to restrict the plotting to. 
#' Default = NULL imposes no restriction (plots whole lattice)
#' @param change.par whether to allow the vis function to optimally change the \code{par} of the plot. Default = TRUE. 
#' If \code{par} is allowed to changed internally, it is always reset upon function exit. 
#' 
#' @export
vis_som_fences = function(SOM, add = FALSE, fence.lwd = 2, clamp = "none", clamp.range = c(0,1), 
                          color.mapping = "linear", color.palette = c('black','white'), color.range = c(0,1), color.nbins = "none", 
                          active = T, subset = NULL, change.par = TRUE)  {
  if(!class(SOM)=="Rcpp_SOM") stop("Input SOM must be an object of class SOM")
  
  ## *** Perform input checks
  stopifnot(color.mapping == "linear" || color.mapping == "zcdf")
  stopifnot(length(color.range) == 2) 
  stopifnot(color.range[1] <= color.range[2])
  stopifnot(min(color.range) >= 0)
  stopifnot(max(color.range) <= 1)
  stopifnot(color.nbins == "none" || color.nbins == "auto" || is.numeric(color.nbins))
  stopifnot(clamp == "none" || clamp == "identity" || clamp == "quantile")
  stopifnot(length(clamp.range) == 2)
  stopifnot(clamp.range[1] <= clamp.range[2])
  
  
  ## *** Decode color palette and create color Ramp function based on this palette 
  #color.palette = process_color.palette(color.palette)
  #color.fxn = colorRamp(colors = color.palette, interpolate = "spline")
  if(any(!SOMDisco:::.is_valid_color(color.palette))) stop("color.palette contains invalid color names")
  color.fxn = colorRamp(colors = color.palette, interpolate = "spline")
  
  
  
  ## *** Extract the fence data frame, and discard any fences involving non-active neurons 
  fencedf = SOM$fences
  plot_these = rep(TRUE,nrow(SOM$fences))
  if(active && SOM$is_recalled) {
    keep_these = which(SOM$RF_size > 0)
    #fencedf = base::subset(fencedf, j %in% keep_these | k %in% keep_these)
    fencedf = base::subset(fencedf, i %in% keep_these & j %in% keep_these)
  }
  
  ## *** Map the fence values to unit range based on requesting mapping function, then subset 
  if(color.mapping == "linear") {
    fencedf$value = SOMDisco:::.unitmap_linscale(x = fencedf$value, clamp = clamp, clamp.range = clamp.range)
  } else {
    fencedf$value = SOMDisco:::.unitmap_zcdf(x = fencedf$value, clamp = clamp, clamp.range = clamp.range)
  }
  if(!is.null(subset)) {
    fencedf = base::subset(fencedf, i %in% subset | j %in% subset)
  }
  
  
  ## *** Map the unit fence values to the requested color range 
  fencedf$value = (fencedf$value - 0) / 1 * diff(range(color.range)) + min(color.range)
  
  ## *** Bin the fence values, if requested
  if(is.numeric(color.nbins)) {
    h = seq(0, 1, length.out = color.nbins+1)
    fencedf$value = h[cut(x = fencedf$value, breaks = h, include.lowest = T, labels = F) + 1]
  } else if(color.nbins == "auto") {
    h = hist(fencedf$value, breaks='scott', plot = F)$breaks
    fencedf$value = h[cut(x = fencedf$value, breaks = h, include.lowest = T, labels = F) + 1]
  }

  ## *** Assign colors based on value, and overwrite any NA values with the requested color 
  fencedf$color = character(nrow(fencedf))
  fencedf$color = rgb(color.fxn(fencedf$value), maxColorValue = 255)

  ## Setup plot, if not adding 
  if(!add) vis_som_setup(SOM, subset = subset, change.par = change.par)
  
  ## Add fence segments 
  if(change.par) {
    opar = par(no.readonly = T)
    on.exit(par(opar))
    #par(.SOMDisco_vis_som_par$par)
    par(SOM$vis_par)
  }
  
  segments(x0 = fencedf$x0, y0 = fencedf$y0, x1 = fencedf$x1, y1 = fencedf$y1, col = fencedf$color, lty = 1, lwd = fence.lwd)
}



#' Visualize SOM prototypes at their lattice locations
#' 
#' @param SOM a SOM object
#' @param add whether to create a new SOM visualization panel (=FALSE, default), or add to an existing one (=TRUE)
#' @param netrng whether to plot the prototypes in internal ("int") or external ("ext") network range. Default = "int". 
#' @param gutter the proportion of the tile width / height used for plotting axes in each tile. 
#' The actual space occupied by the axes = gutter * min(tile_width, tile_height)
#' @param wgt.lwd the line width used for plotting the weight vectors 
#' @param wgt.col the color used for plotting the weight vectors 
#' @param axes boolean, whether to plot individual axes inside each lattice tile 
#' @param axes.lwd the line width used for plotting axes 
#' @param axes.col the color used for plotting axes 
#' @param axes.nticksx the number of ticks on x-axis inside each lattice tile
#' @param axes.nticksy the number of ticks on y-axis inside each lattice tile 
#' @param active Optional, if the SOM object has been recalled,
#' restricts plotting only to active neurons (those whose RF_size > 0).
#' Default = TRUE.
#' @param subset Optional, a vector of neuron indices to restrict the plotting to. 
#' Default = NULL imposes no restriction (plots whole lattice)
#' @param change.par whether to allow the vis function to optimally change the \code{par} of the plot. Default = TRUE. 
#' If \code{par} is allowed to changed internally, it is always reset upon function exit. 
#' 
#' @export
vis_som_prototypes = function(SOM, add = FALSE, netrng = "int", gutter = 0.05, 
                               wgt.lwd = 1, wgt.col = tableau20("blue"), 
                               axes = T, axes.lwd = 0.5, axes.col = "black", axes.nticksx = 4, axes.nticksy = 4, 
                               active = T, subset = NULL, change.par = TRUE) {
  if(!class(SOM)=="Rcpp_SOM") stop("Input SOM must be an object of class SOM")
  if(nrow(SOM$W) == 0) stop("SOM weights do not exist; initialize or train them")
  netrng = tolower(netrng)
  if(!(netrng %in% c('int','ext'))) stop("netrng must be 'int' or 'ext'")
  
  ## Determine which neurons to plot 
  plot_these = rep(TRUE,SOM$nW)
  if(active && SOM$is_recalled) {
    plot_these[SOM$RF_size==0] = FALSE
  }
  if(!is.null(subset)) {
    plot_these[setdiff(1:SOM$nW, subset)] = FALSE
  }
  
  
  ## *** Each weight vector will be scaled to be plotted inside a "normalized" box
  ## whose dimensions depend on the type of lattice (roughly, [0,1] x [0,1]). 
  ## These boxes will then be shifted around to their proper lattice locations.  
  ## The 0.9 term represents a fudge factor to leave room for drawing axes in each box later 
  if(SOM$lattice_type == "grid") {
    tile_width = 1
    tile_height = 1
  } else if(SOM$lattice_type == "hex") {
    tile_width = 1
    tile_height = 1 / sqrt(3)
  }
  offset = min(gutter * tile_width, gutter * tile_height) ## leave space between tile perimeter and plot axes 
  box_width = tile_width - offset - 1/2*offset # width - left gutter - right gutter 
  box_height = tile_height - offset - 1/2*offset # height - bottom gutter - top gutter 
  
  ## *** Scale prototypes to the box 
  ## x-axis values 
  box_xvals = (1:SOM$d - 1) / (SOM$d-1) * box_width
  ## y-axis values, matrix with 1 row for every prototype 
  if(netrng == "int") {
    W = SOM$W
  } else {
    W = SOM$map_from_netrng(SOM$W)
  }
  box_yvals = (W - min(W)) / diff(range(W)) * box_height
  box_yvals = box_yvals[plot_these,,drop=F]
  
  
  ## *** Find coords of the bottom-left corner of each tile 
  boxes_x0y0 = SOM$nu_xy
  boxes_x0y0[,1] = boxes_x0y0[,1] - 1/2 * tile_width 
  boxes_x0y0[,2] = boxes_x0y0[,2] - 1/2 * tile_height
  boxes_x0y0 = boxes_x0y0[plot_these,,drop=F]
  
  ## *** Shift each box to its lattice location 
  boxes_xvals = matrix(box_xvals, byrow = T, nrow = sum(plot_these), ncol = SOM$d)
  boxes_xvals = boxes_xvals + boxes_x0y0[,1] + offset
  boxes_yvals = box_yvals + boxes_x0y0[,2] + offset 
  
  
  ## *** Plot the SOM weights
  if(!add) {
    vis_som_setup(SOM = SOM, subset = subset, change.par = change.par)
    vis_som_tiles(SOM = SOM, add = T, active = active, border.lwd = 0.25, border.col = tableau20("gray"), fill = NULL, subset = subset, change.par = change.par)
  }
  
  if(change.par) {
    opar = par(no.readonly = T)
    on.exit(par(opar))
    #par(.SOMDisco_vis_som_par$par)
    par(SOM$vis_par)
  }
  
  matplot(x = t(boxes_xvals), y = t(boxes_yvals), type = "l", lty = 1, lwd = wgt.lwd, col = wgt.col, add = T)
  
  
  if(axes) {
    ## Add axes
    ## x-axis
    segments(x0 = boxes_x0y0[,1] + offset,
             y0 = boxes_x0y0[,2] + offset,
             x1 = boxes_x0y0[,1] + tile_width - 1/2 * offset,
             y1 = boxes_x0y0[,2] + offset, lty = 1, lwd = axes.lwd, col = axes.col)
    ## y-axis
    segments(x0 = boxes_x0y0[,1] + offset,
             y0 = boxes_x0y0[,2] + offset,
             x1 = boxes_x0y0[,1] + offset,
             y1 = boxes_x0y0[,2] + tile_height - 1/2 * offset, lty = 1, lwd = axes.lwd, col = axes.col)
    
    ## Add x-tick marks 
    if(axes.nticksx > 0) {
      box_xticks = seq(offset, tile_width, length.out = axes.nticksx + 2) ## don't put a tick at x=0 or x=end, they're hard to see 
      box_xticks = box_xticks[2:(axes.nticksx+1)]
      ## x-axis 
      for(i in 1:axes.nticksx) {
        segments(x0 = boxes_x0y0[,1] + box_xticks[i], 
                 y0 = boxes_x0y0[,2] + 1/2 * offset, 
                 x1 = boxes_x0y0[,1] + box_xticks[i], 
                 y1 = boxes_x0y0[,2] + offset, 
                 lty = 1, lwd = axes.lwd, col = axes.col)
      }
    }
    
    ## Add y-tick marks
    if(axes.nticksy > 0) {
      box_yticks = seq(offset, tile_height, length.out = axes.nticksy + 2)
      box_yticks = box_yticks[2:(axes.nticksy+1)]
      for(i in 1:axes.nticksy) {
        segments(x0 = boxes_x0y0[,1] + 1/2 * offset,
                 y0 = boxes_x0y0[,2] + box_yticks[i],
                 x1 = boxes_x0y0[,1] + offset,
                 y1 = boxes_x0y0[,2] + box_yticks[i],
                 lty = 1, lwd = axes.lwd, col = axes.col)
      }
    }
    
  } # close axes 
  
  
}



#' Annotate each tile of the SOM lattice 
#' 
#' @param SOM a SOM object
#' @param add whether to create a new SOM visualization panel (=FALSE, default), or add to an existing one (=TRUE)
#' @param text vector of strings used to label each lattice tile 
#' @param text.cex the text size used for plotting the labels 
#' @param text.col the color used for plotting the labels. Can be one of: 
#' \itemize{
#' \item a single color name, which will be recycled across all lattice tiles
#' \item a vector of color names to be applied to each tile
#' \item the reserved keyword "auto", which will set a text color of either white or black, 
#' depending on the Lightness (from the HSL colorspace) values of any existing tile fill colors. 
#' This is useful if, for example, some tiles are dark (where black text would not show up well) and others are light (where light text would not show up well). 
#' }
#' @param text.font the font used for plotting the labels 
#' @param text.lightness the Lightness threshold used to determine whether the text.col in each tile is white or black. 
#' Only valid if text.col = "auto", default = 0.4. Higher thresholds create more white labels. 
#' @param theta the angle (in degrees) which specifies the direction of offset of each label, 
#' relative to each tile center. 
#' @param rprop the length of the offset of each label, relative to each tile center. 
#' This should be given as a proportion of the total distance from the center to the boundary of each tile, 
#' in the direction theta (i.e., \code{rprop} <= 1)
#' @param active Optional, if the SOM object has been recalled,
#' restricts plotting only to active neurons (those whose RF_size > 0).
#' Default = TRUE.
#' @param subset Optional, a vector of neuron indices to restrict the plotting to. 
#' Default = NULL imposes no restriction (plots whole lattice)
#' @param change.par whether to allow the vis function to optimally change the \code{par} of the plot. Default = TRUE. 
#' If \code{par} is allowed to changed internally, it is always reset upon function exit. 
#' @export
vis_som_annotate = function(SOM, add = FALSE, text, text.cex = 1, text.col = "auto", text.font = 1, text.lightness = 0.4, 
                            theta = 135, rprop = 0, active = T, subset = NULL, change.par = TRUE) {
  if(!class(SOM)=="Rcpp_SOM") stop("Input SOM must be an object of class SOM")
  if(length(text) != SOM$nW) stop("text must be a vector of length = SOM$nW")
  
  ## Determine which labels to plot 
  plot_these = rep(TRUE,SOM$nW)
  if(active && SOM$is_recalled) {
    plot_these[SOM$RF_size==0] = FALSE
  }
  if(!is.null(subset)) {
    plot_these[setdiff(1:SOM$nW, subset)] = FALSE
  }
  plot_these[is.na(text)] = F
  
  ## Setup plot, if requested 
  if(!add) {
    vis_som_setup(SOM = SOM, subset = subset, change.par = change.par)
    vis_som_tiles(SOM = SOM, add = T, active = active, border.lwd = 0.25, border.col = tableau20("gray"), fill = NULL, subset = subset, change.par = change.par)
  }
  
  ## Get (x,y) locations for labels 
  labels_xy = SOM$tile_interior_point(theta, rprop)
  
  ## Decode color 
  if(length(text.col) == 1) {
    if(text.col == "auto") {
      text.col = rep("black", SOM$nW)
      #Lightness = SOMDisco:::.HSL_Lightness(hex = .SOMDisco_vis_som_par$tile_fill)
      Lightness = SOMDisco:::.HSL_Lightness(hex = SOM$vis_tile_bg)
      text.col[Lightness < text.lightness] = "white"
    } else {
      if(!SOMDisco:::.is_valid_color(text.col)) stop("text.col is not a valid color name.")
      text.col = rep(text.col, SOM$nW)
    }
  } else {
    # If a vector, text.col must contain either NA or valid color names
    color_check = SOMDisco:::.is_valid_color(text.col)
    if(any(na.omit(color_check)==F)) stop("text.col contains invalid color names")
  }
  
  # Cannot plot any text whose color = NA. Remove these, if they exist
  #plot_these[is.na(text.col)] = F 
  
  ## Add labels 
  if(change.par) {
    opar = par(no.readonly = T)
    on.exit(par(opar))
    #par(.SOMDisco_vis_som_par$par)
    par(SOM$vis_par)
  }
  
  text(labels_xy[plot_these,,drop=F], labels = text[plot_these], cex = text.cex, col = text.col[plot_these], font = text.font)
  
}



#' Visualize a color gradient on SOM lattice tiles 
#' 
#' @param SOM a SOM object
#' @param add whether to create a new SOM visualization panel (=FALSE, default), or add to an existing one (=TRUE)
#' @param values the numeric values used for heatmap coloring of tiles  
#' @param clamp controls how values are clamped before unit mapping.  Options are: 
#' "none" (default), meaning no clamping is used; 
#' "identity", meaning the clamp.range is used directly; 
#' "quantile", meaning the clamp.range quantiles are used.  
#' @param clamp.range sets the clamping range. 
#' If clamp = 'none', the min/max of values is used. 
#' If clamp = 'identity', clamp.range should be given as an effective [lo,hi] range to clamp values to. 
#' If clamp = 'quantile', clamp.range should contain the effective [lo,hi] probabilities whose quantiles define the clamp range. 
#' Default = c(0,1), which should be changed if clamp = 'identity'. 
#' @param color.mapping the mapping function used to map values to the range [0,1], which is needed for plotting. 
#' Default ='linear' for a linear mapping. 
#' Can also be 'zcdf', which first computes the z-score of the values, then maps the z-score to [0,1] via the Std. Normal CDF. 
#' @param color.palette the color palette used to represent values. 
#' This should be given as a character vector of color names which will set the colors which represent the 
#' the min (first element) and max (last element) of the \code{values}. Colors for intermediate values between min and max 
#' will be interpolated via \link[grDevices]{colorRamp}.  
#' Should have at least two elements, but can have more. In this case, the colors will be interpolated throughout the 
#' range of those found in \code{color.palette}.  Example:  \code{c('black','lightblue','blue')} will map the min values to black, 
#' the mid values to lightblue and the max values to blue.  
#' Default = c('black',tableau20('lightblue')). 
#' @param color.range a range in [0,1] defining what portion of the color scale is visualized. 
#' This is useful if you want to restrict the plotted color range to the darker or lighter end of the spectrum. 
#' Default = c(0,1)
#' @param color.nbins allows binning of the (unit-mapped) values prior to color assignment. Options are: 
#' \itemize{
#' \item 'none' (default), which maps the values directly to the color range without binning; 
#' \item 'auto', which computes a histogram of the (unit-mapped) values and uses the resulting bins and counts; 
#' \item an integer giving the desired number of bins. 
#' }
#' @param color.NA specifies the color used for any NA values 
#' @param active Optional, if the SOM object has been recalled,
#' restricts plotting only to active neurons (those whose RF_size > 0).
#' Default = TRUE.
#' @param subset Optional, a vector of neuron indices to restrict the plotting to. 
#' Default = NULL imposes no restriction (plots whole lattice)
#' @param change.par whether to allow the vis function to optimally change the \code{par} of the plot. Default = TRUE. 
#' If \code{par} is allowed to changed internally, it is always reset upon function exit. 
#' 
#' @details This is a wrapper function for \code{vis_som_tiles}. The tile colors are deduced from the given values 
#' and mapping arguments, then passed along internally. 
#' 
#' @export
vis_som_gradient = function(SOM, add = FALSE, values, clamp = "none", clamp.range = c(0,1), 
                            color.mapping = "linear", color.palette = c('black',SOMDisco::tableau20("lightblue")), 
                            color.range = c(0,1), color.nbins = "none", color.NA = "black", 
                            active = T, subset = NULL, change.par = TRUE) {
  if(!class(SOM)=="Rcpp_SOM") stop("Input SOM must be an object of class SOM")
  
  ## *** Perform input checks
  stopifnot(length(values) == SOM$nW)
  stopifnot(color.mapping == "linear" || color.mapping == "zcdf")
  stopifnot(clamp == "none" || clamp == "identity" || clamp == "quantile")
  stopifnot(length(clamp.range) == 2)
  stopifnot(clamp.range[1] <= clamp.range[2])
  stopifnot(length(color.range) == 2) 
  stopifnot(color.range[1] <= color.range[2])
  stopifnot(min(color.range) >= 0)
  stopifnot(max(color.range) <= 1)
  stopifnot(color.nbins == "none" || color.nbins == "auto" || is.numeric(color.nbins))
  
  ## *** Decode color palette, create color interpolation function 
  #color.palette = process_color.palette(color.palette)
  #color.fxn = colorRamp(colors = color.palette, interpolate = "spline")
  if(any(!SOMDisco:::.is_valid_color(color.palette))) stop("color.palette contains invalid color names")
  color.fxn = colorRamp(colors = color.palette, interpolate = "spline")
  
  ## *** Overwrite any non-active neurons with NA to not influence the scaling statistics 
  if(active) {
    values[SOM$RF_size==0] = NA   
  }
  
  
  ## *** Map the values to unit range based on requesting mapping function 
  if(color.mapping == "linear") {
    values = SOMDisco:::.unitmap_linscale(x = values, clamp = clamp, clamp.range = clamp.range)
  } else {
    values = SOMDisco:::.unitmap_zcdf(x = values, clamp = clamp, clamp.range = clamp.range)
  }
  
  ## *** Map the unit values to the requested color range 
  values = (values - 0) / 1 * diff(range(color.range)) + min(color.range)

  
  ## *** Bin the values, if requested
  if(is.numeric(color.nbins)) {
    h = seq(0, 1, length.out = color.nbins+1)
    values = h[cut(x = values, breaks = h, include.lowest = T, labels = F) + 1]
  } else if(color.nbins == "auto") {
    h = hist(values, breaks='scott', plot = F)$breaks
    values = h[cut(x = values, breaks = h, include.lowest = T, labels = F) + 1]
  }
  
  
  ## *** Assign colors based on value, and overwrite any NA values with the requested color 
  colors = character(length(values))
  colors[!is.na(values)] = rgb(color.fxn(values[!is.na(values)]), maxColorValue = 255)
  colors[is.na(values)] = color.NA
  SOM$set_vis_tile_bg(colors)

  
  ## *** Plot the heatmap on the tiles using these colors 
  SOMDisco::vis_som_tiles(SOM = SOM, add = add, border.lwd = 0.25, border.col = "gray", fill = colors, active = active, subset = subset, change.par = change.par)
  
  
}

  
#' Visualize monitoring measure of SOM training
#' 
#' @param SOM a SOM object
#' @param vis.SOM boolean, whether to show the SOM with RF_size gradient, default = T
#' @param vis.mtr boolean, whether to show the monitoring history, if available. Default = T
#' @export
vis_som_training = function(SOM, vis.SOM = T, vis.mtr = T) {
  if(!class(SOM)=="Rcpp_SOM") stop("Input SOM must be an object of class SOM")
  if(vis.mtr && length(SOM$mtr_age)==0) {
    warning("No monitoring history found in SOM, resetting vis.mtr = F")
    vis.mtr = F
  }
  
  opar = par(no.readonly = T)
  on.exit(par(opar))
  
  #par(mar = c(4,4,1,1))
  
  ## Setup layout 
  if(vis.SOM && vis.mtr) {
    layoutmat = matrix(c(1,1,1,2,3,4), byrow=T, ncol=3)
    layout(mat = layoutmat, heights = c(2, 1))
  } else if(vis.SOM && !vis.mtr) {
    layoutmat = matrix(c(1), byrow=T, nrow=1)
    layout(mat = layoutmat, heights = c(1))  
  } else if(!vis.SOM && vis.mtr) {
    layoutmat = matrix(c(1,2,3), byrow=T, nrow=1)
    layout(mat = layoutmat, heights = c(1))  
  } else {
    stop("One of vis.SOM  or vis.mtr must be TRUE")
  }
  
  
  if(vis.SOM) {
    par(mar = rep(0.1, 4), xpd = NA, pty = 'm', xaxs = 'i', yaxs='i', xaxt = 'n', yaxt = 'n')  
    SOMDisco::vis_som_mUMatrix(SOM = SOM, add = F, change.par = F)
  }
  
  if(vis.mtr) {
    par(mar = c(4,2,2,1), xaxs = 'r', yaxs = 'r', xaxt = 's', yaxt = 's')
    plot(x = SOM$mtr_age/10000, y = SOM$mtr_RMSQE, lwd = 2, col = SOMDisco::tableau20("orange"), type='l', xlab = 'Age/10k', ylab = NA, main = 'RMSQE')
    plot(x = SOM$mtr_age/10000, y = SOM$mtr_QSI, type='l', col = SOMDisco::tableau20("orange"), xlab = 'Age/10k', ylab = NA, main = 'QSI', lwd = 2)
    plot(x = SOM$mtr_age/10000, y = SOM$mtr_Entropy, type='l', col = SOMDisco::tableau20("orange"), xlab = 'Age/10k', ylab = NA, main = 'Entropy', lwd = 2)
    
    #par(mfrow = opar$mfrow) 
  }
}


#' Visualize the Learning Rate Annealing Schedule
#' 
#' @param SOM a SOM object
#' @export
vis_LRAS = function(SOM) {
  if(!class(SOM)=="Rcpp_SOM") stop("Input SOM must be an object of class SOM")
  
  opar = par(no.readonly = T)
  on.exit(par(opar))
  #oldmar = par()$mar
  
  plotdf = SOM$get_LRAS()
  plotdf = rbind(plotdf[1,], plotdf, plotdf[nrow(plotdf),])
  plotdf$t[1] = 1; plotdf$t[nrow(plotdf)] = ceiling(plotdf$t[nrow(plotdf)-1]*1.05)
  
  par(mfrow=c(2,2), mar = c(2,2,3,1))
  plot(x = plotdf$t/10000, y = plotdf$alpha, type='s', xlab = NA, ylab = NA, main = expression(alpha), col = SOMDisco::tableau20("blue"), lwd = 2)
  plot(x = plotdf$t/10000, y = plotdf$beta, type='s', xlab = NA, ylab = NA, main = expression(beta), col = SOMDisco::tableau20("blue"), lwd = 2)
  
  par(mar = c(5,2,2,1))
  plot(x = plotdf$t/10000, y = plotdf$gamma, type='s', xlab = 'Age/10k', ylab = NA, main = expression(gamma), col = SOMDisco::tableau20("blue"), lwd = 2)
  
  nnhb = max(plotdf$sigma)
  
  
  eta = matrix(0, nrow = nrow(plotdf), ncol = nnhb+1)
  for(i in 1:nrow(plotdf)) {
    this_eta = SOM$calc_eta(plotdf$sigma[i])
    eta[i,1:length(this_eta)] = this_eta
  }
  
  matplot(x = plotdf$t/10000, y = eta, type='s', lty = 1, xlab = 'Age/10k', ylab = NA, main = expression(eta), col = SOMDisco::tableau20()[1:(nnhb+1)], lwd = 2)
  legend(x="bottomright", legend = 0:nnhb, col = SOMDisco::tableau20()[1:(nnhb+1)], lty = 1)
  
  #par(opar)
  #par(mfrow=c(1,1), mar = oldmar)
}


#' Visualize the distribution of labels in each RF 
#' 
#' @description The distribution of data labels mapped to each Receptive Field is visualized as a 
#' proportional pie chart residing in each lattice tile. This function is helpful for determining the 
#' organization of SOM neurons: assuming the data labels accurately reflect differences within the dataset, 
#' prototypes whose receptive fields contain data of the same label are considered 
#' to be better representative than those whose data member labels are mixed. 
#' 
#' @param SOM a SOM object
#' @param add whether to create a new SOM visualization panel (=FALSE, default), or add to an existing one (=TRUE)
#' @param subset Optional, a vector of neuron indices to restrict the plotting to. 
#' Default = NULL imposes no restriction (plots whole lattice)
#' @param change.par whether to allow the vis function to optimally change the \code{par} of the plot. Default = TRUE. 
#' If \code{par} is allowed to changed internally, it is always reset upon function exit. 
#' 
#' @details The internal color table \code{$ctab} is used to map labels to colors. Call \code{$set_ctab} to change this mapping. 
#' 
#' @export
vis_som_labeldist = function(SOM, add = FALSE, subset = NULL, change.par = TRUE) {
  ## CHECKS:
  if(!class(SOM)=="Rcpp_SOM") stop("Input SOM must be an object of class SOM")
  
  if(length(SOM$RF_label_dist) == 0) stop("length(SOM$RF_label_dist)==0. Call SOM$set_RF_label() first.")
  
  ## Determine which neurons to plot. 
  ## Empty RFs cannot have a label distribution (they have no data members), so skip them automatically 
  plot_these = rep(TRUE,SOM$nW)
  plot_these[SOM$RF_size==0] = FALSE
  if(!is.null(subset)) {
    plot_these[setdiff(1:SOM$nW, subset)] = FALSE
  }
  
  ## Setup plot, if not adding 
  if(!add) vis_som_setup(SOM, subset = subset, change.par = change.par)
  
  ## Determine all the unique labels found in the RF_label_dist list 
  ## and make sure their corresponding colors exist in the ctab 
  unq_labels = sort(unique(names(unlist(SOM$RF_label_dist))))
  SOMDisco:::.ctab_check(ctab = SOM$ctab, query_labels = unq_labels)
  
  ## Loop over each RF and plot its distribution separately 
  if(change.par) {
    opar = par(no.readonly = T)
    on.exit(par(opar))
    #par(.SOMDisco_vis_som_par$par)
    par(SOM$vis_par)
  }
  
  for(i in 1:SOM$nW) {
    
    # If this RF is empty, there is no label dist. Skip. 
    if(!plot_these[i]) next
    
    # Otherwise strip out the distribution and remove any counts whose labels are NA from it 
    # This would occur if the labels given to set_RF_label contained any NA values 
    this_dist = SOM$RF_label_dist[[i]]
    these_labels = names(this_dist)
    these_labels_na = is.na(these_labels)
    if(any(these_labels_na)) {
      this_dist = this_dist[!these_labels_na]
      these_labels = these_labels[!these_labels_na]
    }
    
    # Convert the distribution counts to cumulative percents. 
    # Also need to pad the beginning of this vector with 0 for plotting 
    this_dist = this_dist / sum(this_dist)
    this_dist = cumsum(this_dist)
    this_dist = c(0, this_dist)
    
    # Match the colors of these labels to the ctab 
    #colors = ctab$color[match(these_labels,ctab$label)]
    colors = SOMDisco:::.ctab_query(ctab = SOM$ctab, query_labels = these_labels)
    
    # Loop over each slide of the pie chart, plot it separately
    for(j in 1:(length(this_dist)-1)) {
      SOMDisco:::.vis_pie_slice(start.degree = this_dist[j]*360, end.degree = this_dist[j+1]*360, col = colors[j], clock.wise = F, radius = 0.4, center = SOM$nu_xy[i,], border = NA)
    }
  }
  
}


#' Visualize the labels of each RF 
#' 
#' @description If the training data is labeled then the SOM prototypes (and their neurons and RFs) will inherit a label 
#' from the learned mapping. Each prototype is labeled by a plurality vote of the labels of data mapped to their RFs.  
#' 
#' @param SOM a SOM object
#' @param add whether to create a new SOM visualization panel (=FALSE, default), or add to an existing one (=TRUE)
#' @param subset Optional, a vector of neuron indices to restrict the plotting to. 
#' Default = NULL imposes no restriction (plots whole lattice)
#' @param change.par whether to allow the vis function to optimally change the \code{par} of the plot. Default = TRUE. 
#' If \code{par} is allowed to changed internally, it is always reset upon function exit. 
#' 
#' @details This is a wrapped function which calls \code{vis_som_tiles} internally using the internal color table \code{$ctab} 
#' to map the \code{$RF_labels} to colors. Call \code{$set_ctab} to change this mapping. 
#' 
#' @export
vis_som_label = function(SOM, add = FALSE, text.cex = 1, text.font = 1, theta = 90, rprop = 0.75, subset = NULL, change.par = TRUE) {
  ## CHECKS:
  if(!class(SOM)=="Rcpp_SOM") stop("Input SOM must be an object of class SOM")
  #if(length(SOM$RF_label) == 0) stop("length(SOM$RF_label)==0. Call SOM$set_RF_label() first.")
  
  #vis_som_tiles(SOM, add = add, border.lwd = 1, border.lty = 1, border.col = "black", 
  #              fill = SOM$RF_label, ctab = ctab, active = TRUE, subset = subset, change.par = change.par)
  vis_som_tiles(SOM, add = add, border.lwd = 1, border.lty = 1, border.col = "black", 
                fill = SOMDisco:::.ctab_query(ctab = SOM$ctab, query_labels = SOM$RF_label), active = TRUE, subset = subset, change.par = change.par)
  if(text.cex > 0) {
    vis_som_annotate(SOM=SOM, add = T, text = SOM$RF_label, 
                     text.cex = text.cex, text.col = "auto", text.font = text.font, text.lightness = 0.4, 
                     theta = theta, rprop = rprop, active = T, subset = subset, change.par = change.par)  
  }
  
  
  
}



#' Visualize a color gradient on SOM lattice tiles 
#' 
#' @param SOM a SOM object
#' @param add whether to create a new SOM visualization panel (=FALSE, default), or add to an existing one (=TRUE)
#' @param grad.clamp see \code{vis_som_gradient}
#' @param grad.clamp.range see \code{vis_som_gradient}
#' @param grad.color.mapping see \code{vis_som_gradient}
#' @param grad.color.palette see \code{vis_som_gradient}
#' @param grad.color.range see \code{vis_som_gradient}
#' @param grad.color.nbins see \code{vis_som_gradient}
#' @param fence.lwd see \code{vis_som_fences}
#' @param fence.clamp see \code{vis_som_fences}
#' @param fence.clamp.range see \code{vis_som_fences}
#' @param fence.color.mapping see \code{vis_som_fences}
#' @param fence.color.palette see \code{vis_som_fences}
#' @param fence.color.range see \code{vis_som_fences}
#' @param fence.color.nbins see \code{vis_som_fences}
#' @param subset Optional, a vector of neuron indices to restrict the plotting to. 
#' Default = NULL imposes no restriction (plots whole lattice)
#' @param change.par whether to allow the vis function to optimally change the \code{par} of the plot. Default = TRUE. 
#' If \code{par} is allowed to changed internally, it is always reset upon function exit. 
#' 
#' @details This is a wrapper function which layers \code{vis_som_fences} (called with the given \code{fence.*} parameters)  atop 
#' \code{vis_som_gradient} (called with the given \code{grad.*} parameters). 
#' 
#' @references 
#' \insertRef{UMatrix}{SOMDisco}
#' \insertRef{MerenyiJainVillmann}{SOMDisco}
#' 
#' @export
vis_som_mUMatrix = function(SOM, add = FALSE, 
                            grad.clamp = "quantile", grad.clamp.range = c(0.05,0.95), 
                            grad.color.mapping = "linear", grad.color.palette = c('black',SOMDisco::tableau20("lightblue")), 
                            grad.color.range = c(0,1), grad.color.nbins = "auto", 
                            fence.lwd = 2, 
                            fence.clamp = "quantile", fence.clamp.range = c(0.05,0.95), 
                            fence.color.mapping = "linear", fence.color.palette = c('black','white'), 
                            fence.color.range = c(0,1), fence.color.nbins = "auto", 
                            subset = NULL, change.par = TRUE) {
  
  if(!class(SOM)=="Rcpp_SOM") stop("Input SOM must be an object of class SOM")
  #if(!SOM$is_recalled) stop("Must call $recall_SOM before calling vis_som_mUMatrix")
  
  SOMDisco::vis_som_gradient(SOM=SOM, add = add, values = SOM$RF_size, 
                             clamp = grad.clamp, clamp.range = grad.clamp.range, 
                             color.mapping = grad.color.mapping, color.palette = grad.color.palette, 
                             color.range = grad.color.range, color.nbins = grad.color.nbins, color.NA = "white", 
                             active = T, subset = subset, change.par = change.par)
  
  SOMDisco::vis_som_fences(SOM = SOM, add = T, fence.lwd = fence.lwd, clamp = fence.clamp, clamp.range = fence.clamp.range, 
           color.mapping = fence.color.mapping, color.palette = fence.color.palette, color.range = fence.color.range, color.nbins = fence.color.nbins, 
           active = T, subset = subset, change.par = change.par)
  
}




#' CONNvis Visualization 
#' 
#' @param SOM a SOM object
#' @param TRN a TRN object from package \code{TopoRNet}
#' @param add whether to create a new plotting device (=FALSE, default), or add to an existing one (=TRUE)
#' @param nu.pch see \code{vis_som_neurons}
#' @param nu.cex see \code{vis_som_neurons}
#' @param nu.col see \code{vis_som_neurons}
#' @param edge.lwd_range the min/max range of the plotted CONN edges, Default = c(1, 5). 
#' @param active Optional, if the SOM object has been recalled,
#' restricts plotting (of vertices and edges) only to active neurons (those whose RF_size > 0).
#' Default = TRUE.
#' @param subset Optional, a vector of neuron indices to restrict the plotting to. 
#' Default = NULL imposes no restriction (plots whole lattice)
#' 
#' @details The CONNvis visualization is documented in \code{TopoRNet::vis_CONNvis}.  
#' 
#' @references 
#' \insertRef{TasdemirMerenyi2009}{SOMDisco}
#'
#' @export
vis_som_CONNvis = function(SOM, TRN, add = F, nu.pch = 16, nu.cex = 1, nu.col = "black", edge.lwd_range = c(1,5), active = T, subset = NULL) {
  
  # Checks 
  if(!class(SOM)=="Rcpp_SOM") stop("Input SOM must be an object of class SOM")
  if(!SOM$is_recalled) stop("Must call $recall_SOM before calling vis_som_CONNvis")
  if(!any(installed.packages()[,1] == "TopoRNet")) stop("Package TopoRNet not installed, which is required for CONNvis visualization")
  
  if(!class(TRN)=="Rcpp_TRN") stop("Input TRN must be an object of class TRN from package TopoRNet")
  
  
  # Change par to SOM settings 
  if(!add) vis_som_setup(SOM=SOM, active = active, subset = subset)
  
  opar = par(no.readonly = T)
  on.exit(par(opar))
  
  par(SOM$vis_par)
  TopoRNet::vis_CONNvis(TRN = TRN, add = T, vertex.xy = SOM$nu_xy, 
                        vertex.pch = nu.pch, vertex.cex = nu.cex, vertex.col = nu.col, 
                        edge.lwd_range = edge.lwd_range, 
                        vertex.active = active, vertex.subset = subset)
}


#' CADJvis Visualization 
#' 
#' @param SOM a SOM object
#' @param TRN a TRN object from package \code{TopoRNet}
#' @param add whether to create a new plotting device (=FALSE, default), or add to an existing one (=TRUE)
#' @param nu.pch see \code{vis_som_neurons}
#' @param nu.cex see \code{vis_som_neurons}
#' @param nu.col see \code{vis_som_neurons}
#' @param edge.lwd_range the min/max range of the plotted CONN edges, Default = c(1, 5). 
#' @param active Optional, if the SOM object has been recalled,
#' restricts plotting (of vertices and edges) only to active neurons (those whose RF_size > 0).
#' Default = TRUE.
#' @param subset Optional, a vector of neuron indices to restrict the plotting to. 
#' Default = NULL imposes no restriction (plots whole lattice)
#' 
#' @details The CADJvis visualization is documented in \code{TopoRNet::vis_CADJvis}.  
#' 
#' @references 
#' \insertRef{TasdemirMerenyi2009}{SOMDisco}
#'
#' @export
vis_som_CADJvis = function(SOM, TRN, add = F, nu.pch = 16, nu.cex = 1, nu.col = "black", edge.lwd_range = c(1,5), active = T, subset = NULL) {
  
  # Checks 
  if(!class(SOM)=="Rcpp_SOM") stop("Input SOM must be an object of class SOM")
  if(!SOM$is_recalled) stop("Must call $recall_SOM before calling vis_som_CONNvis")
  if(!any(installed.packages()[,1] == "TopoRNet")) stop("Package TopoRNet not installed, which is required for CONNvis visualization")
  
  if(!class(TRN)=="Rcpp_TRN") stop("Input TRN must be an object of class TRN from package TopoRNet")
  
  
  # Change par to SOM settings 
  if(!add) vis_som_setup(SOM=SOM, active = active, subset = subset)
  
  opar = par(no.readonly = T)
  on.exit(par(opar))
  
  par(SOM$vis_par)
  TopoRNet::vis_CADJvis(TRN = TRN, add = T, vertex.xy = SOM$nu_xy, 
                        vertex.pch = nu.pch, vertex.cex = nu.cex, vertex.col = nu.col, 
                        edge.lwd_range = edge.lwd_range, 
                        vertex.active = active, vertex.subset = subset)
}


#' TopoView Visualization 
#' 
#' @param SOM a SOM object
#' @param ADJ an adjacency matrix defining edges connecting SOM neurons (or prototypes)
#' @param add whether to create a new plotting device (=FALSE, default), or add to an existing one (=TRUE)
#' @param nu.pch see \code{vis_som_neurons}
#' @param nu.cex see \code{vis_som_neurons}
#' @param nu.col see \code{vis_som_neurons}
#' @param edge.color line color of plotted edges, default = "darkorange".
#' @param edge.lwd_range the min/max range of the plotted CONN edges, Default = c(1, 5). 
#' @param active Optional, if the SOM object has been recalled,
#' restricts plotting (of vertices and edges) only to active neurons (those whose RF_size > 0).
#' Default = TRUE.
#' @param subset Optional, a vector of neuron indices to restrict the plotting to. 
#' Default = NULL imposes no restriction (plots whole lattice)
#' 
#' @details The TopoView visualization is documented in \code{TopoRNet::vis_TopoView}.  If the input \code{ADJ} is weighted and 
#' \code{edge.lwd_range} spans a non-empty set, the visualized edge widths will represent the edge weights (larger weights = thicker edges). 
#' 
#' @references 
#' \insertRef{Merenyietal2009}{SOMDisco}
#'
#' @export
vis_som_TopoView = function(SOM, ADJ, add = F, nu.pch = 16, nu.cex = 1, nu.col = "black", edge.col = "darkorange", edge.lwd_range = c(1,5), active = T, subset = NULL) {
  # Checks 
  if(!class(SOM)=="Rcpp_SOM") stop("Input SOM must be an object of class SOM")
  if(!SOM$is_recalled) stop("Must call $recall_SOM before calling vis_som_mUMatrix")
  if(!any(installed.packages()[,1] == "TopoRNet")) stop("Package TopoRNet not installed, which is required for TopoView visualization")
  #if(!class(TRN)=="Rcpp_TRN") stop("Input TRN must be an object of class TRN from package TopoRNet")
  if(!(nrow(ADJ) == ncol(ADJ))) stop("ADJ must be square")
  if(!(nrow(ADJ) == SOM$nW)) stop("ADJ must have nrows = ncols = SOM$nW")

  # Change par to SOM settings 
  if(!add) vis_som_setup(SOM=SOM, active = active, subset = subset)
  
  opar = par(no.readonly = T)
  on.exit(par(opar))
  
  par(SOM$vis_par)
  TopoRNet::vis_TopoView(ADJ = ADJ, add = T, 
                         vertex.xy = SOM$nu_xy, vertex.pch = nu.pch, vertex.cex = nu.cex, vertex.col = nu.col, 
                         edge.col = edge.col, edge.lwd_range = edge.lwd_range, 
                         vertex.active = active, vertex.subset = subset)
}
