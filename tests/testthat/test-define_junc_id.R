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

  expect_true(ncol(df_id) == 5)
  expect_true(nrow(df_id) == 3)
  expect_true(df_id$junc_id[1] == "chr1:50-100:+")
})


test_that("generate_junction_id works in case of numbers that are represented as scientific", {

  junc_id <- generate_junction_id(chr = "chr1", start = 500000, end = 1000000, strand = "+")

  expect_true(nchar(junc_id) == 21)

})



test_that("junc2breakpoint works", {

  breakpoint_ids <- junc2breakpoint(toy_junc_id)


  expect_true(length(breakpoint_ids) == length(toy_junc_id))
  expect_true(stringr::str_count(breakpoint_ids[1], "chr2") == 2)

})

test_that("breakpoint2junc works", {

  breakpoint_ids <- junc2breakpoint(toy_junc_id)
  strands <- stringr::str_split_fixed(toy_junc_id, ":", n = 3)[,3]

  junc_id <- breakpoint2junc(breakpoint_ids, strand = strands)

  expect_true(all(junc_id == toy_junc_id))
  expect_error(breakpoint2junc("chr2:152389996-chr5:152392205", strand = "+"))

})

test_that("breakpoint2junc works with toy example", {

  breakpoint_id <- "chr1:500000-chr1:1000000"
  junc_id <- breakpoint2junc(breakpoint_id = breakpoint_id, strand = "+")

  expect_equal(junc_id, "chr1:500000-1000000:+")

})

