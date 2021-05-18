#' Returns all context sequences that can derive from a given junction id
#'
#' @param junc_id A junction id
#' @param transcript_db a GRangesList with transcripts composed of exon ranges
#' @param genome_db a GRangesList with transcripts composed of exon ranges
#' @param window_size number of nucleotides left and right from the "breakpoint"
#'
#' @return A tibble with sorted columns as given above
#'
#' @examples
#'
#'
#'@import dplyr
#'
#'@export
juncid2context <-
  function(junc_id,
           transcript_db,
           genome_db,
           window_size = 200) {
    # junction table
    junc_df <- tibble(junc_id = junc_id) %>%
      separate(
        junc_id,
        sep = "_",
        into = c("chr", "pos1", "pos2", "strand"),
        remove = FALSE
      ) %>%
      mutate(pos1 = as.integer(pos1), pos2 = as.integer(pos2))
    # build GRanges object for both junction positions
    junc_pos1 <-
      GenomicRanges::GRanges(junc_df$chr, junc_df$pos1, strand = junc_df$strand)
    junc_pos2 <-
      GenomicRanges::GRanges(junc_df$chr, junc_df$pos2, strand = junc_df$strand)

    # identify transcripts that relate to junction id
    transcripts4junc <-
      map_junc_transcript(
        junc_id = junc_id,
        junc_pos1 = junc_pos1,
        junc_pos2 = junc_pos2,
        junc_df = junc_df,
        transcript_db = transcript_db
      )

    #get context sequences for all transcripts covering a junction
    context_seq <- lapply(transcripts4junc$enst, function(x) {
      get_context_sequence(
        transcript_id = x,
        junc_pos1 = junc_pos1,
        junc_pos2 = junc_pos2,
        transcript_db = transcript_db,
        strand = junc_df$strand,
        genome_db = genome_db,
        window_size = window_size
      )
    })

    return(context_seq)
  }


#' Returns transcripts in which alternative splicing events are predicted
#'
#' @param junc_id A junction id
#' @param junc_pos1 GRanges object for junction position 1
#' @param junc_pos2 GRanges object for junction position 2
#' @param junc_df A tibble with columns `junc_id`, `chr`, `pos1`, `pos2`, `strand`
#' @param transcript_db a GRangesList with transcripts composed of exon ranges
#'
#' @return A tibble with columns `junc_id`, `chr`, `pos1`, `pos2`, `strand`,
#'   `jidx`, `subjectHits`, `enst`. This tibble returns all transcripts that are
#'   affected by a given junction.
#'
#' @examples
#'
#'@import dplyr
#'
#'@export
# this function returns transcripts to which predicted splicing events represent alternative events
map_junc_transcript <-
  function(junc_id,
           junc_pos1,
           junc_pos2,
           junc_df,
           transcript_db) {

    junc_to_transcript <- tibble()

    # Get transcripts overlapping both junctions
    transcripts_pos1 <-
      GenomicRanges::findOverlaps(junc_pos1, transcript_db) %>% as.data.frame() %>% as_tibble()
    transcripts_pos2 <-
      GenomicRanges::findOverlaps(junc_pos2, transcript_db) %>% as.data.frame() %>% as_tibble()

    if (nrow(transcripts_pos1) > 0 & nrow(transcripts_pos2) > 0) {
      # both coordinates of junction event are located in exon sequence related to a transcript
      # this does not capture acceptor gain in intron sequence
      # exon skipping, 3'ASS acceptor gain in exon (+/-) + 5' ASS donor gain (+/-)
      junc_idx_to_txidx <-
        inner_join(transcripts_pos1, transcripts_pos2, by = c("queryHits", "subjectHits")) %>%
        mutate(enst = names(transcript_db)[subjectHits])
      junc_to_transcript <- junc_df %>%
        mutate(jidx = seq_along(junc_id)) %>%
        left_join(junc_idx_to_txidx, by = c("jidx" = "queryHits"))

    } else if (nrow(transcripts_pos1) > 0 & nrow(transcripts_pos2) == 0) {
      # only first coordinate of junction event is located in exon sequence
      # second coordinate is located in intron sequence
      # intron retention, 3' acceptor gain in intron sequence (+), 5' acceptor gain in intron sequence (+)
      junc_to_transcript <- junc_df %>%
        mutate(jidx = seq_along(junc_id)) %>%
        left_join(transcripts_pos1, by = c("jidx" = "queryHits")) %>%
        mutate(enst = names(transcript_db)[subjectHits])

    } else if (nrow(transcripts_pos1) == 0 & nrow(transcripts_pos2) > 0) {
      # first coordinate of junction event is located in intron sequence
      # only second coordinate is located in exon sequence
      # 5' donor gain in intron sequence (+), 3' acceptor gain in intron sequence (-)
      junc_to_transcript <- junc_df %>%
        mutate(jidx = seq_along(junc_id)) %>%
        left_join(transcripts_pos2, by = c("jidx" = "queryHits")) %>%
        mutate(enst = names(transcript_db)[subjectHits])

    }
    return(junc_to_transcript)
  }


#' Given the coordinates of junctions( junc1_gr, junc2_gr, strand.dir),
#' returns the context sequence
#'
#' @param transcript_id An ENSEMBL transcript id
#' @param transcript_db a GRangesList with transcripts composed of exon ranges
#' @param junc_pos1 GRanges object for junction position 1
#' @param junc_pos2 GRanges object for junction position 2
#' @param genome_db a GRangesList with transcripts composed of exon ranges
#' @param strand_direction strand direction. Shall be "+" or "-"
#' @param window_size number of nucleotides left and right from the "breakpoint"
#'
#' @return A tibble with columns `junc_id`, `chr`, `pos1`, `pos2`, `strand`,
#'   `jidx`, `subjectHits`, `enst`. This tibble returns all transcripts that are
#'   affected by a given junction.
#'
#'
#'@import dplyr
#'
get_context_sequence <-
  function(transcript_id,
           transcript_db,
           junc_pos1,
           junc_pos2,
           strand_direction,
           genome_db,
           window_size = 200) {

    # genomic ranges of wt transcript
    wt_transcript_range <- transcript_db[[transcript_id]]

    # identify exons which overlap with junction
    exon1_index <-
      findOverlaps(wt_transcript_range, junc_pos1)@from
    exon2_index <-
      findOverlaps(wt_transcript_range, junc_pos2)@from
    print(exon2_index - exon1_index)

    #  construct mutated ranges
    mutated_transcript_range <-
      construct_mutated_range(
        exon1_index = exon1_index,
        exon2_index = exon2_index,
        wt_transcript_range = wt_transcript_range,
        strand_direction = strand_direction,
        junction_start = start(junc_pos1),
        junction_end = start(junc_pos2)
      )
    print(mutated_transcript_range)
    if (is_empty(mutated_transcript_range)) {
      return("")
    } else{

      # Add transcript coordinates to genomic coordinates
      mutated_transcript_range <- add_transcript_coordinates(mutated_transcript_range)

      # get mutated transcript sequence
      mutated_transcript_sequence <- genome_db[mutated_transcript_range]
      mutated_transcript_sequence <- unlist(mutated_transcript_sequence)

      # extract window of defined size around junction
      context_info <-
        extract_sequence_window(
          mutated_sequence = mutated_transcript_sequence,
          mutated_ranges = mutated_transcript_range,
          junction_start_range = junc_pos1,
          window.size = window_size
        )

      return(
        paste(
          context_info$context_sequence,
          context_info$position_junction,
          sep = "_"
        )
      )
    }
  }

#' returns genomic range of mutated transcript according to the splice event
#'
#' @param exon1_index index of exon related to first position of junction in
#'   transcript ranges
#' @param exon2_index index of exon related to second position of junction in
#'   transcript ranges
#' @param wt_transcript_range Genomic range of wildtype transcript in which the
#'   alternative splicing event is taking place
#' @param strand_direction strand direction. Shall be "+" or "-"
#' @param junction_start junction start coordinate
#' @param junction_end junction end coordinate
#'
#' @return A tibble with columns `junc_id`, `chr`, `pos1`, `pos2`, `strand`,
#'   `jidx`, `subjectHits`, `enst`. This tibble returns all transcripts that are
#'   affected by a given junction.
#'
#'
#'@import dplyr
#'
construct_mutated_range <-
  function(exon1_index,
           exon2_index,
           wt_transcript_range,
           strand_direction,
           junction_start,
           junction_end) {

    mutated_transcript_range <- as.data.frame(wt_transcript_range)

    if (!(is_empty(exon1_index) | is_empty(exon2_index))) {
      # both junction coordinats are located within an exon
      mutated_transcript_range <- construct_range_both_in_exon(
        exon1_index,
        exon2_index,
        mutated_transcript_range,
        strand_direction,
        junction_start,
        junction_end
      )

    } else if (!is_empty(exon1_index) & is_empty(exon2_index)) {
      # first coordinate located in exon, second coordinate located in intron
      mutated_transcript_range <- construct_range_first_in_exon(
        exon1_index,
        mutated_transcript_range,
        strand_direction,
        junction_start,
        junction_end
      )
    } else if (is_empty(exon1_index) & !is_empty(exon2_index)) {

      # first coordinate located in intron, second coordinate located in exon
      mutated_transcript_range <- construct_range_second_in_exon(exon2_index,
                                     mutated_transcript_range,
                                     strand_direction,
                                     junction_start,
                                     junction_end)
    }
    # check for validity of ranges
    mutated_transcript_range <- GRanges(mutated_transcript_range)
    return(mutated_transcript_range)
  }


#' returns genomic range of mutated transcript for junction with both coordinates
#' in an exon
#'
#' @param exon1_index index of exon related to first position of junction in
#'   transcript ranges
#' @param exon2_index index of exon related to second position of junction in
#'   transcript ranges
#' @param transcript_range Genomic range of mutated transcript in which
#'   the alternative splicing event is taking place
#' @param strand_direction strand direction. Shall be "+" or "-"
#' @param junction_start junction start coordinate
#' @param junction_end junction end coordinate
#'
#' @return a mutated transcript range based on the presented parameters
#'
#'
#'@import dplyr
#'
construct_range_both_in_exon <- function(exon1_index,
                                         exon2_index,
                                         transcript_range,
                                         strand_direction,
                                         junction_start,
                                         junction_end) {
  if (exon2_index == exon1_index) {
    # both coordinates of these junction are within the same exon
    # true junctions should not be within one exon --> false predicted events
    transcript_range <- GRanges()

  } else if (abs(exon2_index - exon1_index) == 1) {
    # CASE 1
    # junction coordinates cover consecutive exons
    # insert junction coordinates into exon end / exon start
    # e.g alternative splice sites

    transcript_range$end[exon1_index] <- junction_start
    transcript_range$start[exon2_index] <-  junction_end

  } else if (abs(exon2_index - exon1_index) != 1) {
    # CASE2
    # junction coordinates cover non-consecutive exons
    # complete exons were removed in these events
    # insert junction coordinates into exon end / exon start
    # remove spliced out exons from granges object
    # e.g. exon skipping events

    transcript_range$end[exon1_index] <- junction_start
    transcript_range$start[exon2_index] <- junction_end
    if (strand_direction == "-") {
      removed_exons_index <-
        seq(from = exon2_index + 1,
            to = exon1_index - 1,
            by = 1)
    } else if (strand_direction == "+") {
      removed_exons_index <-
        seq(from = exon1_index + 1,
            to = exon2_index - 1,
            by = 1)
    }
    transcript_range <-
      transcript_range[-removed_exons_index, ]
  }
  return(transcript_range)
}

#' returns genomic range of mutated transcript for junctions with first
#' coordinate in an exon and seconde coordinate located in an intron
#'
#' @param exon1_index index of exon related to first position of junction in
#'   transcript ranges
#' @param transcript_range Genomic range of mutated transcript in which
#'   the alternative splicing event is taking place
#' @param strand_direction strand direction. Shall be "+" or "-"
#' @param junction_start junction start coordinate
#' @param junction_end junction end coordinate
#'
#' @return a mutated transcript range based on the presented parameters
#'
#'
#'@import dplyr
#'
construct_range_first_in_exon <- function(exon1_index,
                                          transcript_range,
                                          strand_direction,
                                          junction_start,
                                          junction_end) {
  # first coordinate located in exon, second coordinate located in intron
  # e.g intron retention, 3' ASS acceptor gain in intron (+)
  # CASE3
  if (abs(junction_end - junction_start) == 1) {
    # intron retention
    if (strand_direction == "+") {
      transcript_range$end[exon1_index] <-
        transcript_range$end[exon1_index + 1]
      transcript_range <-
        transcript_range[-(exon1_index + 1), ]
    } else if (strand_direction == "-") {
      transcript_range$end[exon1_index] <-
        transcript_range$end[exon1_index - 1]
      transcript_range <-
        transcript_range[-(exon1_index - 1),]
    }
  } else {
    # CASE 4
    if (strand_direction == "+") {
      # ASS 3' on positive strand
      transcript_range$start[(exon1_index + 1)] <-
        junction_end
    } else if (strand_direction == "-") {
      # ASS 5' on negative strand strand
      transcript_range$start[(exon1_index - 1)] <-
        junction_end
    }

  }
  return(transcript_range)
}


#' returns genomic range of mutated transcript for junctions with first
#' coordinate in an exon and seconde coordinate located in an intron
#'
#' @param exon2_index index of exon related to second position of junction in
#'   transcript ranges
#' @param transcript_range Genomic range of mutated transcript in which
#'   the alternative splicing event is taking place
#' @param strand_direction strand direction. Shall be "+" or "-"
#' @param junction_start junction start coordinate
#' @param junction_end junction end coordinate
#'
#' @return a mutated transcript range based on the presented parameters
#'
#'
#'@import dplyr
#'
construct_range_second_in_exon <- function(exon2_index,
                                           transcript_range,
                                           strand_direction,
                                           junction_start,
                                           junction_end){
  if (abs(junction_end - junction_start) == 1) {
    # intron retention
    # CASE 6
    if (strand_direction == "+") {
      transcript_range$start[exon2_index] <-
        transcript_range$start[exon2_index - 1]
      transcript_range <- mutated_transcript_range[-(exon2_index - 1),]
    } else if (strand_direction == "-") {
      transcript_range$start[exon2_index] <-
        transcript_range$start[exon2_index + 1]
      transcript_range <- transcript_range[-(exon2_index + 1), ]
    }
  } else{
    # CASE 5
    # e.g 5' ASS donor gain in intron (+)
    if (strand_direction == "+") {
      exon1_index <- which(!junction_start < transcript_range$start)[length(which(!junction_start < transcript_range$start))]
      transcript_range$end[exon.1] <- junction_start
      if (abs(exon1_index - exon2_index) != 1) {
        # there are some examples of junctions ids with splice donor in intron sequence with additional exon skipping
        exons_removed_index <-
          seq(from = exon1_index + 1,
              to = exon2_index - 1,
              by = 1)
        transcript_range <- transcript_range[-exons_removed_index, ]
      }
      # e.g 3' ASS acceptor gain in intron (-)
    } else if (strand_direction == "-") {
      exon1_index <- which(junction_start > transcript_range$end)[1]
      transcript_range$end[exon.1] <- junction_start
      if (abs(exon1_index - exon2_index) != 1) {
        exons_removed_index <-
          seq(from = exon2_index + 1,
              to = exon1_index - 1,
              by = 1)
        transcript_range <- transcript_range[-exons_removed_index, ]
      }
    }
  }
  return(transcript_range)
}
