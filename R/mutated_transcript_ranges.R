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
#'@import rlang
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
