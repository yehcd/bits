## -------------------------------------------------------------------------------
# essentials for analysis
library(tidyverse)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg38)
library(vroom)

# for summary logo figure generation
library(ggseqlogo)
library(Biostrings)


## -------------------------------------------------------------------------------
# apply boolean based on if "target_seq" has a "CT" at position 8-9
check_CT_match <- function(string, positions = c(8,9)){
  # extract check region 
  # to upper for case-insensitivty input
  check_string <- stringr::str_sub(string = string,
                                   start = positions[1],
                                   end = positions[2]) %>% stringr::str_to_upper()

  # do check
  check_value <- check_string == "CT"
  
  # return
  return(check_value)
}


# parse casOffinder output to GRanges, 
make_casoffinder_granges<- function(tb_casoffinder){
  
  # basic conversion to GRanges type
  tb_casoffinder <- tb_casoffinder %>% mutate(
    target_start = target_position + 1,
    target_end = target_start + str_length(target_seq) - 1,
    gr_CT_start = ifelse(
      test = target_strand == "+",
      yes = target_start + 7,
      no = target_start + 5
    ),
    gr_CT_end = ifelse(
      test = target_strand == "+",
      yes = target_start + 8,
      no = target_start + 6
    ),
    CT_start = ifelse(
      test = target_strand == "+",
      yes = target_start + 7,
      no = target_start + 5
    )
  ) 
  
  # drop unnecssary column
  tb_casoffinder <- tb_casoffinder %>% select(-target_position)

  
  # make granges
  gr_casoffinder <- makeGRangesFromDataFrame(
    df = tb_casoffinder,
    keep.extra.columns = TRUE,
    ignore.strand = FALSE,
    seqnames.field = "target_chr",
    start.field = "gr_CT_start",
    end.field = "gr_CT_end",
    strand.field = "target_strand",
    starts.in.df.are.0based = FALSE,
  )
  
  # return
  return(gr_casoffinder)
}


# Helper function to loop through coverage bigWigs & combine with CasOFFinder sites

combine_bigwig_casoffinder <- function(input_bigwig, output_dir, gr_casoffinder) {
  # load datafile and decompress bamCoverage merging
  gr_bigwig <- rtracklayer::import(con = input_bigwig, format = "bigWig") %>% fillEmptyGrangesBins(tileWidth = 10) %>% sum_score_overlaps(
    gr_casoffinder = gr_casoffinder,
    gr_coverage = .,
    casoff_regionSize = 1000,
    cov_tileWidth = 10
  )
  
  #save background datafile
  rds_output <- paste0("./", output_dir, "/", basename(input_bigwig), ".Rds.gz")
  write_rds(x = gr_bigwig,
            file = rds_output,
            compress = "gz")
  
}



## -------------------------------------------------------------------------------
# expand bamCoverage sites
fillEmptyGrangesBins <- function(gr.data, tileWidth = 10) {

  ## remap/bin coverage using GRanges functions
  # make tiles of target size
  gr.tiles <-
    slidingWindows(x = gr.data, width = tileWidth, step = tileWidth) %>% unlist()


  # name of quantitative metric is "score" by default for coverage data
  rlel.data.cov <- coverage(gr.data, weight = "score")


  # make output
  gr.output <-
    binnedAverage(bins = gr.tiles,
                  numvar = rlel.data.cov,
                  varname = "score")

  # output
  return(gr.output)
}


# main function to get score
#   gr_casoffinder    Granges of CasOffinder defining the exact sites to score, CENTERED on CT
#   gr_coverage       Output of bamCoverage
#   casoff_regionSize   size of tiles in 'gr_casoffinder' (i.e., 1000bps)
#   cov_tileWidth       'binSize' parameter used in bamCoverage (i.e., 10bps)

sum_score_overlaps <- function(gr_casoffinder,
                               gr_coverage,
                               casoff_regionSize = 1000,
                               cov_tileWidth = 10) {
  # comput basic params
  casoff_regionWidth <- casoff_regionSize / 2
  min_overlap <- as.integer(cov_tileWidth / 2) + 1 # avoid double-count
  
  # define upstream and downstream regions
  gr_casoffinder_roi_upstream <- GenomicRanges::resize(x = gr_casoffinder, width = casoff_regionWidth, fix = "start")
  
  gr_casoffinder_roi_downstream <- GenomicRanges::resize(x = gr_casoffinder, width = casoff_regionWidth, fix = "end")
  
  
  # HELPER FUNCTION to compute ROI scores
  # set minimum overlaps score scoringing
  
  
  get_roi_score <- function(gr_coverage, gr_casoffinder_roi) {
    # get overlaps for scoring
    overlaps <- GenomicRanges::findOverlaps(query = gr_casoffinder_roi,
                                            subject =
                                              gr_coverage,
                                            minoverlap = min_overlap)
    
    # get score
    overlap_scores <- tapply(
      INDEX =  queryHits(overlaps),
      X = mcols(gr_coverage)$score[subjectHits(overlaps)],
      FUN = sum
    )
    
    return(overlap_scores)
  }
  
  # compute scores
  scores_upstream <- get_roi_score(gr_coverage = gr_coverage, gr_casoffinder_roi = gr_casoffinder_roi_upstream)
  scores_downstream <- get_roi_score(gr_coverage = gr_coverage, gr_casoffinder_roi = gr_casoffinder_roi_downstream)
  
  
  # attach score to sites
  gr_casoffinder$score_upstream <- 0 # initialize
  gr_casoffinder$score_downstream <- 0 # initialize
  
  mcols(gr_casoffinder)$score_upstream[as.integer(names(scores_upstream))] <- scores_upstream
  mcols(gr_casoffinder)$score_downstream[as.integer(names(scores_downstream))] <- scores_downstream
  
  # return
  return(gr_casoffinder)
  
}



## -------------------------------------------------------------------------------
# Get hits unfiltered list of hits
get_hits_unfiltered <- function(gr_signal, gr_background, guide_seq = NULL) {

  
  # convert to tibble
  tb_signal <- gr_signal %>% as_tibble()
  tb_background <- gr_background %>% as_tibble()
  
  # basic join & compute some basic values
  tb_merged <- left_join(
    x = tb_signal, # signal
    y = tb_background, # background
    by = c(
      "seqnames",
      "start",
      "end",
      "strand",
      "gRNA_query",
      "target_seq",
      "mm_count",
      "CT_match",
      "target_start",
      "target_end",
      "CT_start"
    ),
    suffix = c(".sample", ".bkgrnd")
  ) %>% mutate(
    log2FE_upstream = log2((score_upstream.sample + 1) / (score_upstream.bkgrnd + 1)),
    log2FE_downstream = log2((score_downstream.sample + 1) / (score_downstream.bkgrnd + 1)),
    bkgrndSub_upstream = score_upstream.sample - score_upstream.bkgrnd,
    bkgrndSub_downstream = score_downstream.sample - score_downstream.bkgrnd
  )
  
  # compute a mixed FE x CPM total upstream/downstream score value
  tb_hits <- tb_merged %>% mutate(
    sum_FExCPM_score = (abs(log2FE_upstream) * bkgrndSub_upstream) + (abs(log2FE_downstream) * bkgrndSub_downstream)
  )
  
  # guide_seq filtering
  if (!is.null(guide_seq)) {
    tb_hits <- tb_hits %>% filter(gRNA_query == guide_seq)
  }
  
  # basic sorting from highest to lowest
  tb_hits <- tb_hits %>% arrange(desc(sum_FExCPM_score))
  
  # return
  return(tb_hits)
  
}


# basic wrapper for the filtering based on the criteria of 
# a) minimum signal in ROI (guide target site +/-1kbps)  being 1CPM (up & downstream) AFTER background subtraction
# b) minimum fold enrichment of signal/background being (log2 = 1)
# c) up/downstream signal are within 20% of each other (i.e., minimal symmetry
# d) remove anonamous super over-represented sites (score >10000)

apply_filtering_basic <- function(tb_hits_unfiltered) {
  # filter for symmetry
  tb_filtered <- tb_hits_unfiltered %>% filter(
    abs(bkgrndSub_upstream - bkgrndSub_downstream) / pmax(bkgrndSub_upstream, bkgrndSub_downstream) < 0.2
  )
  
  # filter for minimum CPM of 1 in up/downstream ROI after background subtraction
  tb_filtered <- tb_filtered %>% filter(bkgrndSub_upstream >= 1, bkgrndSub_downstream >= 1)
  
  # filter for minimum log2fold enrichment signal/background being 1 for up/downstream ROI
  tb_filtered <- tb_filtered  %>% filter(log2FE_upstream >= 1, log2FE_downstream >= 1)
  
  # filter for excessively overrepresented & non-specific
  tb_filtered <- tb_filtered %>% filter(sum_FExCPM_score < 10000)
  
  # return
  return(tb_filtered)
}


# apply blacklist (CUT&RUN)
apply_blacklist <- function(tb_input, gr_blacklist){
  # make GRanges of tb_input for findOverlaps() usage
  gr_input <- tb_input %>% select(-width.sample) %>%  makeGRangesFromDataFrame(
    keep.extra.columns = TRUE,
    ignore.strand = FALSE,
    seqnames.field = "seqnames",
    start.field = "start",
    end.field = "end"
  )
  
  # get overlap between blacklist & hits file
  overlaps <- findOverlaps(query = gr_blacklist, subject = gr_input)
  overlaps
  
  # remove blacklist hits
  gr_final <- gr_input[-subjectHits(overlaps)]
  
  # as tibble
  tb_final <- gr_final %>% as_tibble()

  # return
  return(tb_final)
}


## -------------------------------------------------------------------------------
# basic function to compute 
do_analysis <- function(file_background_rds,
                        file_sample_rds,
                        guide_seq,
                        gr_blacklist,
                        input_dir,
                        output_dir,
                        score_threshold = 30,
                        CT_filter = TRUE) {
  # Load pre-computed Rds data
  path_data_bkgrd <- paste0(input_dir, "/", file_background_rds)
  path_data_sample <- paste0(input_dir, "/", file_sample_rds)
  
  gr_bkgrd <- read_rds(path_data_bkgrd)
  gr_treated <- read_rds(path_data_sample)
  
  # GET RAW HITS SCORING (UNFILTERED)
  tb_hits <- get_hits_unfiltered(gr_signal = gr_treated,
                                 gr_background = gr_bkgrd,
                                 guide_seq = guide_seq)
  
  # APPLY BLACKLIST
  tb_hits <- tb_hits %>% apply_blacklist(gr_blacklist = gr_blacklist_merged)
  
  # APPLY BASIC FILTERING
  tb_hits <- tb_hits %>% apply_filtering_basic()
  
  # APPLY MINIMUM SCORE FILTERING & IMPOSE CT REQUIREMENT
  tb_hits <- tb_hits %>% filter(sum_FExCPM_score >= score_threshold)
  if (CT_filter == TRUE) {
    tb_hits <- tb_hits %>% filter(CT_match == CT_filter)
  }
  
  
  # save output
  output_name <- file_sample_rds %>% basename() %>% sub("\\.cov\\.bw\\.Rds\\.gz$", "", .) %>% paste0(output_dir, "/", ., "_", guide_seq, ".4.thresholded.csv")
  cat(output_name, "\n")
  
  write_csv(x = tb_hits, file = output_name)
}

