# get_allele_freqs() turns the allele x marker data frames into the named
# list that pedtools wants. Mostly checking it doesn't mix up rows and
# columns, since that would silently give the wrong frequencies.

test_that("get_allele_freqs reshapes a small table correctly", {
  df <- data.frame(Allele = c(8, 9, 9.3, 10),
                   M1 = c(0.1, 0.2, 0.3, 0.4),
                   M2 = c(0, 0.5, 0, 0.5))
  out <- get_allele_freqs(df)

  expect_type(out, "list")
  expect_named(out, c("M1", "M2"))
  expect_identical(out$M1, c(`8` = 0.1, `9` = 0.2, `9.3` = 0.3, `10` = 0.4))
  expect_identical(out$M2, c(`8` = 0, `9` = 0.5, `9.3` = 0, `10` = 0.5))
})

test_that("get_allele_freqs keeps every marker of the bundled databases", {
  dbs <- list(Argentina = Argentina, Asia = Asia, Austria = Austria,
              BosniaHerz = BosniaHerz, China = China, Europe = Europe,
              Japan = Japan, USA = USA)
  for (nm in names(dbs)) {
    db <- dbs[[nm]]
    out <- get_allele_freqs(db)
    expect_identical(names(out), setdiff(names(db), "Allele"), info = nm)
    for (m in names(out)) {
      expect_identical(names(out[[m]]), as.character(db$Allele), info = nm)
    }
  }
})

test_that("get_allele_freqs picks the right cell for a given allele", {
  # look up D3S1358 allele 15 directly in the table, then via the list
  out <- get_allele_freqs(Argentina)
  row <- which(Argentina$Allele == 15)
  expect_length(row, 1)
  expect_equal(unname(out$D3S1358["15"]), Argentina$D3S1358[row])
  expect_gt(out$D3S1358[["15"]], 0)
})

test_that("get_allele_freqs output works with pedtools::setMarkers", {
  skip_if_not_installed("pedtools")
  freqs <- get_allele_freqs(Argentina)
  x <- pedtools::nuclearPed(1)
  x <- pedtools::setMarkers(x, locusAttributes = freqs["D3S1358"])

  af <- pedtools::afreq(x, "D3S1358")
  expect_equal(af[names(freqs$D3S1358)], freqs$D3S1358)
})
