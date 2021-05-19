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
    transcripts_covering_junction <-
      map_junc_transcript(
        junc_id = junc_id,
        junc_pos1 = junc_pos1,
        junc_pos2 = junc_pos2,
        junc_df = junc_df,
        transcript_db = transcript_db
      )

    #get context sequences for all transcripts covering a junction
    context_sequences <- lapply(transcripts_covering_junction$enst, function(x) {
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

    return(context_sequences)
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
          transcript_sequence = mutated_transcript_sequence,
          transcript_range = mutated_transcript_range,
          junc_pos1 = junc_pos1,
          window_size = window_size
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



#' adds transcript coordinates in the resulting transcript to the genomic
#' ranges of the exome
#'
#' @param transcript_range Genomic range of a transcript
#'
#' @return The Genomic range of a transcript with transcript coordinates
#'
#'
#'@import dplyr
#'
#
add_transcript_coordinates <- function(transcript_range) {
  transcript_data <- as.data.frame(transcript_range)
  transcript_data$start_transcript <-
    rep(0, length(transcript_range))
  transcript_data$end_transcript <-  rep(0, length(transcript_range))
  if (nrow(transcript_data) == 1) {
    # only one exon
    transcript_data$end_transcript <- width(transcript_range)
  } else {
    for (i in 2:length(transcript_range)) {
      transcript_data$start_transcript[i] =
        transcript_data$start_transcript[i -1] +
        width(transcript_range)[i - 1] + 1
    }
    for (i in 1:length(transcript_range)) {
      transcript_data$end_transcript[i] = transcript_data$start_transcript[i] +
        width(transcript_range)[i]
    }
  }
  extended_transcript_range <- GRanges(transcript_data)
  return(extended_transcript_range)
}


#' adds transcript coordinates in the resulting transcript to the genomic
#' ranges of the exome
#'
#' @param transcript_sequence Sequence of the full transcript
#' @param transcript_range Genomic range of a transcript
#' @param junc_pos1 Genomic range of the first junction coordinate
#' @param window_size number of nucleotides left and right from the "breakpoint"
#'
#'
#' @return The Genomic range of a transcript with transcript coordinates
#'
#'
#'@import dplyr
#'
#
extract_sequence_window <-
  function(transcript_sequence,
           transcript_range,
           junc_pos1,
           window_size = 200) {
    transcript_data <- as.data.frame(transcript_range)
    junction_exon_index <-
      findOverlaps(transcript_range, junc_pos1)@from
    junction_position_transcript <-
      transcript_data$end_transcript[junction_exon_index]
    if (junction_position_transcript - window_size > 0) {
      start_pos <- junction_position_transcript - window_size
      junction_position_context <- window_size
    } else{
      start_pos <- 1
      junction_position_context <- junction_position_transcript
    }
    if (junction_position_transcript + window_size < length(transcript_sequence)) {
      end_pos <- junction_position_transcript + window_size
    } else{
      end_pos <- length(transcript_sequence)
    }
    context.sequence <-
      as.character(transcript_sequence[start_pos:end_pos])
    res <-
      list(position_junction = junction_position_context, context_sequence = context.sequence)
    return(res)
  }
