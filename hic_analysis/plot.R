plot_hic <- function(mat, res, chr, label = '', plot_col_bar = TRUE){
  old_par <- par()

  #get the din, get min value of width and height

  x_mar_1 <- c(0.05, 0.13) # margin of x and y
  x_mar_2 <- c(0.90, 0.08)
  y_mar_1 <- c(0.09, 0.09)
  y_mar_2 <- c(0.2, 0.2)



  w_h_ratio <- old_par$pin[1]/old_par$pin[2]

  if(w_h_ratio < 1){ # adjust height
    y_mar_1 <- y_mar_1/w_h_ratio
    y_mar_2 <- y_mar_2/w_h_ratio
  } else { # adjust width
    x_mar_1 <- x_mar_1 * w_h_ratio
    x_mar_2 <- x_mar_2 * w_h_ratio
  }

  # adjust coordinate
  x_coor_1 <- c(x_mar_1[1], 1 - x_mar_1[2])
  x_coor_2 <- c(x_mar_2[1], 1 - x_mar_2[2])
  y_coor_1 <- c(y_mar_1[1], 1 - y_mar_1[2])
  y_coor_2 <- c(y_mar_2[1], 1 - y_mar_2[2])

  # plot main heatmap
  # set params
  par(plt=c(x_coor_1,y_coor_1))
  q95 <- quantile(mat,.95)[1]
  breaks <- c(seq(0,q95,length.out = 100),max(mat))
  image_col <- colorRampPalette(c('white','red'))(100)
  # plot
  image(rot90c(mat), breaks = breaks, col = image_col, axes = FALSE, useRaster = T, main = label)
  # add axes
  chr_len <- refine_chrlen(dim(mat)[1] * res)
  axis(side = 1, at = c(1), labels = c(chr_len), tick = FALSE, mgp = c(3,0.5,0))
  box()

  # plot colbar
  if(plot_col_bar){
  par(plt=c(x_coor_2, y_coor_2),new=TRUE)
  col_bar_matrix <- matrix(1:100,nrow = 1)

    vmin <- 0
    vmax <- q95
    y    <- seq(vmin,vmax,len=100)

    image(1,y,col_bar_matrix, col = image_col, axes = FALSE, useRaster = TRUE, xlab = "", ylab = "")
    box()
  #add value

    axis(side = 4, at = c(vmin,vmax), labels = round(c(vmin,vmax),0), tick = F, las = 1, mgp = c(3,0.2,0))
    invisible()
  }
}

simple_hic_plot <- function(mat, mapCol = 'red', vmax=.95){
  max_v <- quantile(mat, vmax)[1]
  breaks <- c(seq(0, max_v, length.out = 100), max(mat))
  image_col <- colorRampPalette(c('white', mapCol))(100)
  # plot
  image(rot90c(mat), breaks = breaks, col = image_col, axes = FALSE, useRaster = T)
  box()
  invisible()
}

add_loops <- function(loop_file, res, chr = NULL, region = NULL, pos='both',col = 'yellow', lty = 1, enlarge = T, blow = 5,...){
# read in file
  loops <- read.table(loop_file,TRUE)
  sta <- region[1]
  end <- region[2]
  region_len <- end - sta
  loops <- subset(loops,chr1 == chr & (((x1+x2)/2 > sta & (x1+x2)/2 < end) | ((y1+y2)/2 > sta & (y1+y2)/2 < end)))
  loops <- within(loops,{x1_x <- (x1 - sta)/region_len
                         x2_x <- (x2 - sta)/region_len
                         x1_y <- (region_len-(x1 - sta))/region_len
                         x2_y <- (region_len-(x2 - sta))/region_len
                         y1_y <- (region_len-(y1 - sta))/region_len
                         y2_y <- (region_len-(y2 - sta))/region_len
                         y1_x <- (y1 - sta)/region_len
                         y2_x <- (y2 - sta)/region_len})
  if(enlarge){#enlarge 5bin
    blow  <- blow*1/(region_len/res)
    loops <- within(loops,{x1_x <- x1_x - blow
                           x2_x <- x2_x + blow
                           y1_x <- y1_x - blow
                           y2_x <- y2_x + blow
                           x1_y <- x1_y + blow
                           x2_y <- x2_y - blow
                           y1_y <- y1_y + blow
                           y2_y <- y2_y - blow})
  }

  if(pos == 'both'){
    rect_fun <- function(x1_x,x2_x,y1_x,y2_x,x1_y,x2_y,y1_y,y2_y){
      rect(x1_x,y1_y,x2_x,y2_y,lty = lty, border = col, ...)
      rect(y1_x,x1_y,y2_x,x2_y,lty = lty, border = col, ...)
    }

  }
  if(pos == 'left'){
    rect_fun <- function(x1_x,x2_x,y1_x,y2_x,x1_y,x2_y,y1_y,y2_y){
      rect(x1_x,y1_y,x2_x,y2_y,lty = lty, border = col, ...)
    }
  }
  if(pos == 'right'){
    rect_fun <- function(x1_x,x2_x,y1_x,y2_x,x1_y,x2_y,y1_y,y2_y){
      rect(y1_x,x1_y,y2_x,x2_y,lty = lty, border = col, ...)
    }
  }
  mapply(rect_fun,x1_x=loops$x1_x,x2_x=loops$x2_x,y1_x=loops$y1_x,y2_x=loops$y2_x,
         x1_y=loops$x1_y,x2_y=loops$x2_y,y1_y=loops$y1_y,y2_y=loops$y2_y)
  invisible()
}

simple_pearson_plot <- function(mat, vmax=.95){
  q95_pos <- quantile(mat[mat > 0], .95)[1]
  q95_neg <- quantile(mat[mat < 0], .5)[1]
  q95     <- min(-q95_neg,q95_pos) #red and blue color range should be same
  minV    <- min(mat)
  maxV    <- max(mat)

  image_col <- colorRampPalette(c('blue','black','red'))(100)

  breaks    <- c(minV, seq(-q95, q95, length.out = 99), maxV)


  # plot
  image(rot90c(mat), breaks = breaks, col = image_col, axes = FALSE, useRaster = T)
  box()
}

plot_hic_from_file <- function(hic_file, chr, res, norm='KR', label=chr, plot_col_bar = TRUE){
  plot_hic(get_mat_from_hic(hic_file, norm, chr, res), res, chr, label, plot_col_bar)
  invisible()
}

# a modified version of rect function
rotate_square <- function (xleft, ybottom, xright, ytop, density = NULL, angle = 45,
                           col = NA, border = NULL, lty = par("lty"), lwd = par("lwd"),
                           ...)
{

    n <- range(length(xleft), length(xright), length(ybottom),
               length(ytop))
    if (n[1L] == 0)
      stop("invalid rectangle specification")
    n <- n[2L]
    x <- rbind(rep.int(NA, n), xleft, (xleft + xright)/2, xright, (xleft + xright)/2)[-1L]
    y <- rbind(rep.int(NA, n), (ybottom + ytop)/2, ybottom, (ybottom + ytop)/2, ytop)[-1L]
    polygon(x, y, col = col, border = border, lty = lty,
            lwd = lwd, density = density, angle = angle, ...)

  invisible()
}


# rotate 45 degree
# TODO: 1. need set plt to adjust device size
#       2. need set max distance

plot_rotate_mat <- function(mat, chr, res, start_bin = NA, end_bin = NA, max_dis = NA, label = '',
                            tad_file = NA,
                            plot_axis = TRUE){

  max_len <- max(c(mat[,1],mat[,2]))
  # 2. get coordinate of each point
  i       <- mat[,1]
  j       <- mat[,2]
  x_coor  <- (i + j)/2 - 0.5
  y_coor  <- abs(j - i)/2

  xleft   <- x_coor - 0.5
  xright  <- x_coor + 0.5
  ybottom <- y_coor - 0.5
  ytop    <- y_coor + 0.5

  # 3. get color of each point
  value   <- mat[,3]
  col_fun <- colorRampPalette(c('white','red'))
  breaks  <- c(seq(0,quantile(value,.95),length.out = 100),max(value))
  col_map <- col_fun(100)[as.numeric(cut(value,breaks = breaks))]

  # 4. plot
  ## set ymax to adjust max_dis
  if(is.na(max_dis)){ymax = max_len/2}
  else{ymax = max_dis/2}

  pos_pair   <- refine_pair_chrlen((start_bin - 1) * res, end_bin * res)
  start_pos  <- pos_pair[1]
  end_pos    <- pos_pair[2]
  plot_title <- paste(chr,': from', start_pos, 'to', end_pos, '  resolution:', refine_chrlen(res))

  plot(c(0,max_len),c(0, ymax), xaxs = 'i', yaxs = 'i',type = 'n', axes = F, xlab = '', ylab = label, mgp = c(1,1,0))

  if(plot_axis){
    axis(side = 1, at = c(0,max_len), tick = FALSE, labels = pos_pair, mgp = c(3,0.5,0))
  }

  rotate_square(xleft, ybottom, xright, ytop, col = col_map, border = NA)

  # TODO : add tad plot
  if(!is.na(tad_file)){
    tad_data <- read.table(tad_file)
    tad_data <- tad_data[tad_data[,1] == chr,] # extract chromosome
    tad_data[,2] <- tad_data[,2]/res # get relative coordinate
    tad_data[,3] <- tad_data[,3]/res # get relative coordinate
    tad_data <- tad_data[(tad_data[,2] >= start_bin & tad_data[,2] < end_bin),]
    tad_data <- cbind(tad_data[,2] - start_bin, tad_data[,3] - start_bin)

    # now col1 and col2 is start and end bin of tad
    x_coor <- rbind(tad_data[,1], (tad_data[,1] + tad_data[,2])/2, tad_data[,2], rep(NA, nrow(tad_data))) %>% as.vector()
    y_coor <- rbind(rep(0,nrow(tad_data)), (tad_data[,2] - tad_data[,1])/2, rep(0, nrow(tad_data)),rep(NA, nrow(tad_data))) %>% as.vector()
    lines(x_coor, y_coor, col = 'yellow', lty = 2)
  }
  # 5. add a color legend


}

plot_rotate_with_border <- function(chr, res, start_bin = NA, end_bin = NA, border_file){
  #plot_rotate_mat(mat, chr, res, start_bin, end_bin, max_dis, label, plot_axis = plot_axis)
  border <- read.table(border_file)
  border <- border[border$V1 == chr,]
  border_pos <- (border$V2 + border$V3)/2/res + 0.5
  border_pos <- border_pos[border_pos >= start_bin & border_pos <= end_bin]
  border_pos <- border_pos - start_bin
  abline(v = border_pos, lty = 3)

  invisible()
}

plot_rotate_with_loop <- function(chr, res, start_bin = NA, end_bin = NA, loop_file){

  loop_data <- read.table(loop_file, header = TRUE)
  loop_data <- loop_data[loop_data$chr1 == chr,]
  loop_x_leftpos   <- with(loop_data, (y1 + x1)/2/res - start_bin + 1)
  loop_x_rightpos  <- with(loop_data, (y2 + x2)/2/res - start_bin + 1)
  loop_y_bottompos <- with(loop_data, (y1 - x2)/2/res)
  loop_y_toppos    <- with(loop_data, (y2 - x1)/2/res)
  rect(loop_x_leftpos,loop_y_bottompos,loop_x_rightpos,loop_y_toppos, border = 'blue')

  invisible()
}

plot_rotate_with_bed <- function(chr, res, start_bin = NA, end_bin = NA, max_dis = NA, bed_file){
  bed_data <- read.table(bed_file)
  bed_data <- bed_data[bed_data$V1 == chr,]
  bed_data$V2 <- bed_data$V2/res
  bed_data$V3 <- bed_data$V3/res
  bed_data <- bed_data[bed_data$V2 >= start_bin & bed_data$V2 <= end_bin,]
  x_pos <- with(bed_data,rbind(V2 - start_bin, V3 - start_bin, V3 - start_bin, V2 - start_bin, rep(NA, nrow(bed_data)))) %>% as.vector()
  y_pos <- rep(c(0, 0, max_dis, max_dis, NA), nrow(bed_data))
  color <- rgb(0,0,1, alpha = 0.2)
  polygon(x_pos, y_pos, col = color, border = NA)
}

plot_rotate_from_hic <- function(hic_file, chr, res,
                                 start_bin = NA, end_bin = NA, max_dis = NA,
                                 label = '', plot_axis = TRUE, tad_file = NA, border_file = NA,
                                 loop_file = NA, bed_file = NA,norm= 'KR'){

  mat <- get_mat_from_hic(hic_file, norm, chr, res)

  if(!is.na(start_bin) & !is.na(end_bin)){ # extract expect mat
    mat <- mat[start_bin:end_bin, start_bin:end_bin] %>% summary()
  }
  else{
    start_bin <- 1
    end_bin   <- dim(mat)[1]
  }
  plot_rotate_mat(mat, chr, res, start_bin, end_bin, max_dis, label, tad_file, plot_axis)

  if(!is.na(border_file)){
    plot_rotate_with_border(chr,res,start_bin,end_bin,border_file)
  }

  if(!is.na(loop_file)){
    plot_rotate_with_loop(chr,res,start_bin,end_bin,loop_file)
  }

  if(!is.na(bed_file)){
    plot_rotate_with_bed(chr,res,start_bin,end_bin,max_dis,bed_file)
  }

  invisible()
}


plot_multi_rotate <- function(hic_file_list, chr, res, start_bin = NA, end_bin = NA, max_dis = NA,
                              label_list, norm = 'KR', pdf_file = 'plot_multi_rotate.pdf'){
  # adjust asp
  w_h <- (end_bin - start_bin + 1) / ((max_dis / 2) * length(hic_file_list))

  # only for pdf plot
  if(w_h > 1){
    pdf_w <- 7 * w_h
    pdf_h <- 7 * 3
  }
  else{
    pdf_h <- (7 / w_h) * 3
    pdf_w <- 7
  }


  pdf(pdf_file, pdf_w, pdf_h)

  par(mfrow = c(length(hic_file_list),1),
      mar = c(1,3,1,2) + 0.1,
      oma = c(1,0,0,0))


  i <- 1
  for(hic_file in hic_file_list){
    if(i < length(hic_file_list)){
      plot_rotate_from_hic(hic_file, chr, res, start_bin, end_bin, max_dis, label_list[i], FALSE, norm)
    }
    else{
      plot_rotate_from_hic(hic_file, chr, res, start_bin, end_bin, max_dis, label_list[i], TRUE, norm)
    }
    # last plot
    i <- i+1
  }

  dev.off()
  invisible()
}



# plot heatmap with border strengh or border index
plot_rotate_with_index <- function(hic_file, chr, res, start_bin = NA, end_bin = NA, max_dis = NA,
                                   label, norm = 'KR',
                                   index_file, index_res, index_lab = '', index_type = 'line',
                                   pdf_file = 'rotate_with_index.pdf'){
  # adjust asp
  w_h <- (end_bin - start_bin + 1) / ((max_dis / 2) * 2)

  # only for pdf plot
  if(w_h > 1){
    pdf_w <- 7 * w_h
    pdf_h <- 7
  }
  else{
    pdf_h <- 7 / w_h
    pdf_w <- 7
  }


  pdf(pdf_file, pdf_w, pdf_h)

  par(mfrow = c(2,1),
      mar = c(1,3,1,2) + 0.1,
      oma = c(1,0,0,0))
  plot_rotate_from_hic(hic_file, chr, res, start_bin, end_bin, max_dis, label, TRUE, norm)

  # adjust index resolution with heatmap resolution
  res_ratio       <- index_res/res
  index_start_bin <- start_bin/res_ratio
  index_end_bin   <- end_bin/res_ratio
  index_data      <- read.table(index_file)
  index_x         <- index_start_bin:index_end_bin
  index_y         <- index_data[index_x,1]

  if(index_type == 'line'){
    plot(index_x,index_y,type = 'l',xaxs = 'i', yaxs = 'i', xlab = '', ylab = index_lab, axes = FALSE)
    axis(side = 2)
  }
  else if(index_type == 'bar'){
    barplot(index_y,space = 0,xaxs = 'i', col = ifelse(index_y < 0,'blue','red'),
            border = NA,ylab = 'DI index',
            axes = FALSE)
    axis(side = 2)
  }
  dev.off()
}



# 1. define a curve function from two locate
get_arc_function <- function(loc1,loc2,height){

  # check distance, if larger than 2*height, use circle func, else use curve func

  dist <- loc2 - loc1

  if(dist > 2*height){
    # get centre point and radius

    circle_x <- (loc2+loc1)/2
    circle_y <- (height^2 - (dist/2)^2)/(2*height)
    radius   <- height - circle_y

    function(x){ (radius^2 - (x-circle_x)^2)^0.5 + circle_y }
  }
  else{
    a <- loc1
    b <- loc2
    k <- height/(a*b - (a+b)^2/4)

    function(x){ k*(x-a)*(x-b) }
  }
}

# 2. get arc plot function
# TODO: set threshold of lowest value ?
arc_plot <- function(mat, height, col = 'black'){
  # 1. set ylim and xlim

  xmin <- 0
  xmax <- max(mat[,1:2])
  ymin <- 0
  ymax <- 1.5 * height

  plot(c(xmin,xmax),c(ymin,ymax), type='n', xaxs='i', yaxs='i', xlab='', ylab='', axes=FALSE)

  # 2. plot arc

  col_function <- colorRampPalette(c('white',col))

  breaks       <- c(seq(0,quantile(mat[,3],.70),length.out = 10),max(mat[,3]))

  mat[,3]      <- col_function(10)[as.numeric(cut(mat[,3], breaks = breaks))]

  # 3. define plot function for apply

  curve_plot <- function(loc1,loc2,col,height){

    loc1       <- as.numeric(loc1) - .5
    loc2       <- as.numeric(loc2) - .5
    curve_func <- get_arc_function(loc1,loc2,height)

    curve(curve_func,loc1,loc2, add = TRUE, col = col)
    invisible()
  }

  ploting <- apply(mat, 1, function(x) curve_plot(x[1], x[2], x[3], height))
  invisible()
}







