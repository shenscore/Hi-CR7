# detect sparse region
# input : hic_file
# output: sparse region
# outputformat:
# chr start end
detect_sparse_region <- function(hic_file, chrlist, res = 5000, window_size = 10, step = 4){

  entry_num_list <- list()
  window_list    <- data.frame()

  for (chr in chrlist){
    hic_mat      <- get_mat_from_hic(hic_file, 'NONE', chr, res)
    mat_size     <- dim(hic_mat)[1]
    window_start <- seq(1,mat_size - window_size + 1, by = step)
    window_end   <- seq(window_size, mat_size, by = step)
    chromosome   <- rep(chr, length(window_start))

    window       <- data.frame(chromosome, window_start, window_end)
    window_list  <- rbind(window_list, window)
    # check each window if it's sparse
    get_entry_num <- function(start,end){
      nnzero(hic_mat[start:end,start:end])
    }
    entry_num <- mapply(get_entry_num, start = window_start, end = window_end, SIMPLIFY = TRUE)
    entry_num_list[[chr]] <- entry_num
  }

  all_entry <- do.call(c, entry_num_list)
  q10       <- quantile(all_entry, .1)
  orig_sparse_region <- window_list[which(all_entry < q10),]
  merged_sparse_region <- merge_adjacent_region(orig_sparse_region, chrlist)

  return(entry_num_list)
}

merge_adjacent_region <- function(region, chrlist = NULL, max_dis = 20){

  merged_region <- list()
  for (chr in chrlist){
    intra_region <- region[region[,1] == chr,]
    mid_pos <- (intra_region[,2] + intra_region[,3])/2
    jump_point   <- which(diff(mid_pos) > max_dis)
    pos_start    <- c(1,jump_point + 1)  # find jump point that distance > max_dis
    pos_end      <- c(jump_point, nrow(intra_region))       # define merged region's start and end position in original region

    chromosome           <- rep(chr, length(pos_start))
    merged_region_start  <- intra_region[pos_start,2]
    merged_region_end    <- intra_region[pos_end,3]
    merged_region[[chr]] <- data.frame(chromosome, merged_region_start, merged_region_end)
  }
  merged_region <- do.call(rbind, merged_region)
}
