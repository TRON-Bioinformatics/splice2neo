skip("No automatical visualisation testing!")

plot_tracks <- function(trck1, trck2, main = NULL){
  tracks("canonical junction" = trck1, "alternative junctions" = trck2,
         heights = c(1, 4),label.text.cex = 0.5, main = main) +
    theme_bw()+
    theme( # remove the vertical grid lines
      panel.grid.major.y = element_blank() ,
      panel.grid.minor.y = element_blank() ,
    )
}

test_that("modify_tx works on toy example data on positive strand", {

  # example data
  transcripts <- GenomicRanges::GRangesList(list(
    tx1 = GenomicRanges::GRanges(
      c("1", "1", "1"),
      IRanges::IRanges(
        c(2, 15, 27),
        c(8, 23, 35),
      ),
      strand = c("+", "+", "+")
    )
  ))

  # junction examples
  junc_df = tibble::tribble(
    ~"chr", ~"pos1", ~"pos2", ~"strand", ~"comment",
    "1",  17, 20, "+", "exitron",
    "1",  6, 15, "+", "A5SS (exon)",
    "1",  10, 15, "+", "A5SS (intron)",
    "1",  8, 18, "+", "A3SS (exon)",
    "1",  8, 13, "+", "A3SS (intron)",
    "1",  8, 27, "+", "exon skipping",
    "1",  8, 9, "+", "intron retention",
  ) %>%
    mutate(
      junc_id = generate_junction_id(chr, pos1, pos2, strand)
    )
  junc_gr <- junc_to_gr(junc_df$junc_id)
  names(junc_gr) <- junc_df$comment



  txs <- rep(transcripts,length(junc_gr))
  names(txs) <- names(junc_gr)

  transcripts_mod <- modify_tx(txs, junc_gr)

  transcripts_mod <- c(transcripts, transcripts_mod)

  require(ggbio)

  p1 <- ggplot() +
    geom_alignment(transcripts_mod$tx1)#+
  #geom_point(aes(x=18, y=1))
  p2 <- ggplot()+
    geom_alignment(transcripts_mod[2:8], fill = "brown")

  tks <- plot_tracks(p1, p2, main = "positive strand")

  tks

  expect_true(length(transcripts_mod)-1 == length(junc_gr))

})

test_that("modify_tx works on toy example data with mutually exclusive exons", {

  # example data
  transcripts <- GenomicRanges::GRangesList(list(
    tx1 = GenomicRanges::GRanges(
      c("1", "1", "1", "1"),
      IRanges::IRanges(
        c(2, 15, 27, 40),
        c(8, 23, 35, 45),
      ),
      strand = c("+", "+", "+", "+")
    )
  ))

  # junction examples
  junc_df = tibble::tribble(
    ~"chr", ~"pos1", ~"pos2", ~"strand", ~"comment",
    "1",  8, 27, "+", "MXE1",
    "1",  23, 40, "+", "MXE2",
  ) %>%
    mutate(
      junc_id = generate_junction_id(chr, pos1, pos2, strand)
    )
  junc_gr <- junc_to_gr(junc_df$junc_id)
  names(junc_gr) <- junc_df$comment



  txs <- rep(transcripts,length(junc_gr))
  names(txs) <- names(junc_gr)

  transcripts_mod <- modify_tx(txs, junc_gr)

  transcripts_mod <- c(transcripts, transcripts_mod)

  require(ggbio)

  p1 <- ggplot() +
    geom_alignment(transcripts_mod$tx1)
  p2 <- ggplot()+
    geom_alignment(transcripts_mod[2:3], fill = "brown")

  tks <- plot_tracks(p1, p2, main = "positive strand")

  tks


  expect_true(length(transcripts_mod)-1 == length(junc_gr))

})





test_that("modify_tx works on toy example data on negative strand", {

  # example data
  transcripts <- GenomicRanges::GRangesList(list(
    tx1 = GenomicRanges::GRanges(
      c("1", "1", "1"),
      IRanges::IRanges(
        c(2, 15, 27),
        c(8, 23, 35),
      ),
      strand = c("-", "-", "-")
    )
  ))

  # junction examples
  junc_df = tibble::tribble(
    ~"chr", ~"pos1", ~"pos2", ~"strand", ~"comment",
    "1",  17, 20, "-", "exitron",
    "1",  6, 15, "-", "A3SS (exon)",
    "1",  10, 15, "-", "A3SS (intron)",
    "1",  8, 18, "-", "A5SS (exon)",
    "1",  8, 13, "-", "A5SS (intron)",
    "1",  8, 27, "-", "exon skipping",
    "1",  8, 9, "-", "intron retention",
  ) %>%
    mutate(
      junc_id = generate_junction_id(chr, pos1, pos2, strand)
    )
  junc_gr <- junc_to_gr(junc_df$junc_id)
  names(junc_gr) <- junc_df$comment



  txs <- rep(transcripts,length(junc_gr))
  names(txs) <- names(junc_gr)

  transcripts_mod <- modify_tx(txs, junc_gr)

  transcripts_mod <- c(transcripts, transcripts_mod)

  require(ggbio)

  p1 <- ggplot() +
    geom_alignment(transcripts_mod$tx1)
  p2 <- ggplot()+
    geom_alignment(transcripts_mod[2:8], fill = "brown")

  tks <- plot_tracks(p1, p2, main = "negative strand")

  tks


  expect_true(length(transcripts_mod)-1 == length(junc_gr))

})


