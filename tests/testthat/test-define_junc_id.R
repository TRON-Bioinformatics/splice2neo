test_that("generate_junction_id works", {

  junc_id <- generate_junction_id(chr = "chr1", start = 50, end = 100, strand = "+")

  junc_id_split <- str_split_fixed(junc_id, ":", n=3)

  expect_true(ncol(junc_id_split) == 3)
  expect_true(junc_id_split[,3] =="+")
})

test_that("generate_junction_id works on tibbles", {

  df <- tibble(chr = c("chr1", "chr2", "chr3"),
               start = c(50, 60, 80),
               end = c(100, 110, 190),
               strand = c("+", "-", "+")
               )
  df_id <- df %>%
    mutate(junc_id = generate_junction_id(chr, start, end, strand))

  expect_true(ncol(df_id) == 3)
  expect_true(nrow(df_id) == 3)
  expect_true(df_id$junc_id[1] == "chr1:50-100:+")
})


test_that("generate_junction_id works in case of numbers that are represented as scientific", {

  junc_id <- generate_junction_id(chr = "chr1", start = 500000, end = 1000000, strand = "+")

  expect_true(nchar(junc_id) == 21)

})
