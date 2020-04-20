# This is the source code for the manuscript "The biogenenis and function of nucleosome arrays.":

# The fastq files for MNase and ATAC samples were processed as described in the Methods section of the manuscript. All publicly available softwares
# and their version are stated in the Methods section.

# BigWig files for each sample are available from the GEO repository.





############################################################################################################################################################

# NRL and array regularity calculation

ocampo2 <- function (coverage = coverage, references = reference, beforeRef = 200, 
                     afterRef = 800, smoothingWindow = 75, spacing.low = 130, 
                     spacing.high = 220, shift.low = -70, shift.high = 60, mc.cores = 2,
                     sigma_scaled = FALSE) 
{
  
  gaussmf <- function(x, sigma, mean) {
    height <- 1
    mfVals = (height * exp(-(((x - mean)^2)/(2 * (sigma^2)))))
  }
  
  Pattern <- list()
  
  
  if(sigma_scaled){
    for (d in spacing.low:spacing.high) {
      Pattern[[d]] <- apply(sapply(0:9, function(x) {
        ### !!! scale sigma !!! ###
        gaussmf(seq(-beforeRef, afterRef), 40*(d/150), x * d)
      }), 1, sum)
    }
  } else {
    for (d in spacing.low:spacing.high) {
      Pattern[[d]] <- apply(sapply(0:9, function(x) {
        ### !!! scale sigma !!! ###
        gaussmf(seq(-beforeRef, afterRef), 40, x * d)
      }), 1, sum)
    }
  }
  
  windows <- data.frame(chr = references$chr, 
                        start = ifelse(references$strand == "+", 
                                       references$start - (beforeRef + (smoothingWindow/2)), 
                                       references$end - (afterRef + (smoothingWindow/2))), 
                        end = ifelse(references$strand == "+", 
                                     references$start + (afterRef + (-1 + smoothingWindow/2)), 
                                     references$end + (beforeRef + (-1 + smoothingWindow/2))), 
                        strand = references$strand)
  rownames(windows) <- rownames(references)
  mat <- coverageWindowsStranded(windows, coverage)
  res <- parallel::mclapply(1:nrow(mat), mc.cores = mc.cores, 
                            function(ridx) {
                              rx <- mat[ridx, ]
                              x <- zoo::rollmean(rx, smoothingWindow)
                              bestR <- 0
                              spacingV <- NA
                              shiftV <- NA
                              if (round(sd(x), 3) != 0) {
                                for (d in spacing.low:spacing.high) {
                                  
                    
                                  y <- Pattern[[d]]
                                  
                                  my_ccf <- ccf(x,y, lag.max = abs(shift.low), plot=FALSE)
                                  
                                  r <- max(my_ccf$acf)
                                  shift <- my_ccf$lag[my_ccf$acf == r]
                                  
                                  if (r > bestR) {
                                    bestR <- r
                                    shiftV <- shift
                                    spacingV <- d
                                    
                                    
                                  }
                                  
                                  
                                }
                              }
                              c(r = round(bestR, 2), space = spacingV, shift = shiftV)
                            })
  df <- t(as.data.frame(res))
  rownames(df) <- rownames(mat)
  df
}

############################################################################################################################################################
















############################################################################################################################################################

# Generate matrix for nucleosome dyad centers from coverage file

coverageWindowsStranded <- function(windows,  coverage) {
  
  #  cl <- makeCluster(getOption("cl.cores", 8))
  #  clusterExport(cl, list("centers","coverage","window.size") , envir=environment())
  
  windows <- windows[windows$chr %in% names(coverage),]
  
  #  result <- parSapply(cl, names(cov), function(x) {
  result <- lapply(names(coverage), function(x) {
    my.cov <- coverage[[x]]
    my.windows <- windows[windows$chr==x,]
    mw.views <- IRanges::Views(my.cov, start=my.windows$start, my.windows$end)
    ## remove out-of bounds views
    flt <- start(mw.views)>0 & end(mw.views) < length(my.cov)
    mw.views <- mw.views[flt,]
    my.windows <- my.windows[flt,]
    if (length(mw.views) > 0) {
      mat <- as.matrix(mw.views)
      rownames(mat) <- rownames(my.windows)
      return(mat)
    } else {
      return(NULL)
    }
  })
  #  stopCluster(cl)
  mat <- Reduce(rbind, result)
  windows <- windows[rownames(windows) %in% rownames(mat),]
  match(rownames(windows), rownames(mat)) -> o
  mat <- mat[o,]
  mat[windows$strand=="-",] <- t(apply(mat[windows$strand=="-",],1,rev))
  mat
}

############################################################################################################################################################












############################################################################################################################################################

# coverage from bam file

bam2dyad <- function(bam_file, 
                     type = c("PAIRED", "SINGLE"),
                     chromosomes = c("I","II", "III"),
                     average_fraglength = 150,
                     min_fraglength = 140,
                     max_fraglength = 160,
                     subsample_number = 10e6,
                     subsample_withreplacement = TRUE,
                     smooth_width = 50,
                     divide_bytotal = TRUE){
  
  
  if (type == "SINGLE"){
    # read single reads
    single_reads <- readGAlignments(file = bam_file)
    read_ranges <- granges(single_reads)
    
    # select chromosomes
    #read_ranges <- keepSeqlevels(read_ranges, value = chromosomes, pruning.mode = "coarse")
    
    
    # sub select by number of reads
    if(!(is.null(subsample_number))){
      set.seed(100)
      my_sub_index <- sample(length(read_ranges), subsample_number, replace = subsample_withreplacement)
      read_ranges <- read_ranges[my_sub_index]
    }
    
    # resize fragments to average fragments
    frag_ranges <- GenomicRanges::resize(read_ranges, width = average_fraglength, fix = "start")
    
    # resize fragments to smoothed dyads
    dyads <- GenomicRanges::resize(frag_ranges, width = smooth_width, fix = "center")
    
    # coverage
    my_cov <- GenomicRanges::coverage(dyads)
    
    # reads per million
    if(divide_bytotal){
      my_cov <- my_cov / sum(unlist(my_cov)) * 1e6
    }
    
    return(my_cov)
    
  } else {
    # read readpairs
    paired_reads <- readGAlignmentPairs(file = bam_file)
    frag_ranges <- granges(paired_reads, on.discordant.seqnames="drop")
    
    # select chromosomes
    #frag_ranges <- keepSeqlevels(frag_ranges, value = chromosomes, pruning.mode = "coarse")
    
    # select by fragment length
    my_frag_index <- width(frag_ranges) > min_fraglength & width(frag_ranges) < max_fraglength
    frag_ranges <- frag_ranges[my_frag_index]
    
    
    # sub select by number of reads
    if(!(is.null(subsample_number))){
      set.seed(100)
      my_sub_index <- sample(length(frag_ranges), subsample_number, replace = subsample_withreplacement)
      frag_ranges <- frag_ranges[my_sub_index]
    }
    
    # resize fragments to smoothed dyads
    dyads <- GenomicRanges::resize(frag_ranges, width = smooth_width, fix = "center")
    
    # coverage
    my_cov <- GenomicRanges::coverage(dyads)
    
    # reads per million
    if(divide_bytotal){
      my_cov <- my_cov / sum(unlist(my_cov)) * 1e6
    }
    
    return(my_cov)
  }
}

############################################################################################################################################################
