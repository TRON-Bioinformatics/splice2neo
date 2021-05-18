#' Returns all context sequences that can derive from a given junction id
#'
#' @param junc_id A junction id
#' @param transcript_db a GRangesList with transcripts composed of exon ranges
#'
#' @return A tibble with sorted columns as given above
#'
#' @examples
#'
#' unsorted_junc_df
#' sorted_junc_df <- sort_columns(unsorted_junc_df)
#' sorted_junc_df
#'
#'@import dplyr
#'@import GenomicRanges
#'
#'@export
juncid2context <-
  function(junc_id,
           transcript_db,
           genome.seqs,
           window.size = 200) {
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
      GRanges(junc_df$chr, junc_df$pos1, strand = junc_df$strand)
    junc_pos2 <-
      GRanges(junc_df$chr, junc_df$pos2, strand = junc_df$strand)

    # identify transcripts that relate to junction id
    transcripts4junc <-
      map_junc_transcript(
        junc_id = junc_id,
        junc1_gr = junc_pos1,
        junc2_gr = junc_pos2,
        junc_df = junc_df,
        transcript_db = transcript_db
      )

    #get context sequences for all transcripts covering a junction
    context_seq <- lapply(transcripts4junc$enst, function(x) {
      get_context_sequence(
        transcript_id = x,
        junc1_gr = junc1_gr,
        junc2_gr = junc2_gr,
        transcript_db = transcript_db,
        strand.dir = junc_df$strand,
        genome.seqs = genome.seqs,
        window.size = window.size
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
#' unsorted_junc_df
#' sorted_junc_df <- sort_columns(unsorted_junc_df)
#' sorted_junc_df
#'
#'@import dplyr
#'@import GenomicRanges
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
      findOverlaps(junc_pos1, transcript_db) %>% as.data.frame() %>% as_tibble()
    transcripts_pos2 <-
      findOverlaps(junc_pos2, transcript_db) %>% as.data.frame() %>% as_tibble()

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
#' unsorted_junc_df
#' sorted_junc_df <- sort_columns(unsorted_junc_df)
#' sorted_junc_df
#'
#'@import dplyr
#'@import GenomicRanges
#'
#'@export
get_context_sequence <-
  function(transcript_id,
           transcript_db,
           junc1_gr,
           junc2_gr,
           strand.dir,
           genome.seqs,
           window.size) {
    # genomic ranges of wt transcript
    transcript.wt <- transcript_db[[transcript_id]]

    #=============================================================================
    # 2.) identify exons which overlap with junction
    #=============================================================================
    start.junc.ind <-
      GenomicRanges::findOverlaps(transcript.wt, junc1_gr)@from
    end.junc.ind <-
      GenomicRanges::findOverlaps(transcript.wt, junc2_gr)@from
    print(end.junc.ind - start.junc.ind)
    #=============================================================================
    # 2.) construct mutated ranges
    #=============================================================================
    transcript.mut <-
      construct_mutated_range(
        exon.index1 = start.junc.ind,
        exon.index2 = end.junc.ind,
        transcript_range.wt = transcript.wt,
        strand.dir = strand.dir,
        junction_start = start(junc1_gr),
        junction_end = start(junc2_gr)
      )
    print(transcript.mut)
    if (is_empty(transcript.mut)) {
      return("")
    } else{
      #=============================================================================
      # 3.) Add transcript coordinates to genomic coordinates
      #=============================================================================
      transcript.mut <- add_transcript_coordinates(transcript.mut)

      #=============================================================================
      # 4.) mutated transcript sequence
      #=============================================================================
      mutated.transcript.seq <-
        BSgenome::getSeq(genome.seqs, transcript.mut)
      mutated.transcript.seq <- unlist(mutated.transcript.seq)
      #print(mutated.transcript.seq)

      #=============================================================================
      # 5.) Calculate window of defined size around junction
      #=============================================================================
      # transcript sequence coordinates around junction
      context_info <-
        extract_sequence_window(
          mutated_sequence = mutated.transcript.seq,
          mutated_ranges = transcript.mut,
          junction_start_range = junc1_gr,
          window.size = window.size
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
