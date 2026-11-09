suppressMessages({library(pedtools); library(forrel); library(mispitools); library(parallel)})
suppressMessages({library(pedtools)})
three <- function() {
  p <- nuclearPed(3)
  p <- addSon(p, 3, id = "c3"); p <- addSon(p, 4, id = "c4"); p <- addSon(p, 5, id = "MP")
  p
}
scen <- list(
  full_sibling        = list(label = "Full sibling",                ped = relabel(nuclearPed(2), new = c("F","M","R1","MP")), refs = "R1"),
  half_sibling        = list(label = "Half sibling",                ped = relabel(halfSibPed(), new = c("F","M1","M2","R1","MP")), refs = "R1"),
  grandparent         = list(label = "Grandparent",                 ped = { p <- nuclearPed(1); p <- addSon(p, 3, id = "MP"); p }, refs = "1"),
  uncle               = list(label = "Uncle/Aunt",                  ped = { p <- nuclearPed(2); p <- addSon(p, 4, id = "MP"); p }, refs = "3"),
  two_uncles          = list(label = "Two uncles/aunts",            ped = { p <- nuclearPed(3); p <- addSon(p, 5, id = "MP"); p }, refs = c("3","4")),
  first_cousin        = list(label = "First cousin",                ped = { p <- nuclearPed(2); p <- addSon(p, 3, id = "R1"); p <- addSon(p, 4, id = "MP"); p }, refs = "R1"),
  two_first_cousins   = list(label = "Two first cousins",           ped = three(), refs = c("c3","c4")),
  uncle_first_cousin  = list(label = "Uncle/aunt and first cousin", ped = three(), refs = c("3","c4"))
)
freqs <- lapply(get_allele_freqs(Argentina), function(x) { x <- x[!is.na(x) & x > 0]; x / sum(x) })
args <- commandArgs(TRUE); slug <- args[1]; N <- as.integer(args[2]); cores <- as.integer(args[3]); out_dir <- if (length(args) >= 4) args[4] else "examples/genetic"
s <- scen[[slug]]
la <- lapply(names(freqs), function(m) list(name = m, alleles = names(freqs[[m]]), afreq = unname(freqs[[m]])))
ped <- setMarkers(s$ped, locusAttributes = la)
sing <- setMarkers(singleton("POI"), locusAttributes = la)
one <- function(i) {
  set.seed(1e6 * match(slug, names(scen)) + i)
  sim <- profileSim(ped, N = 1, ids = c(s$refs, "MP"), verbose = FALSE)
  refped <- setGenotype(sim, ids = "MP", marker = seq_len(nMarkers(sim)), geno = "-/-")
  poi1 <- relabel(subset(sim, "MP"), old = "MP", new = "POI")
  un <- profileSim(sing, N = 1, verbose = FALSE)
  lr_rel <- missingPersonLR(refped, missing = "MP", poi = poi1, verbose = FALSE)$LRtotal[[1]]
  lr_un  <- missingPersonLR(refped, missing = "MP", poi = un, verbose = FALSE)$LRtotal[[1]]
  c(lr_un, lr_rel)
}
t0 <- Sys.time()
res <- do.call(rbind, mclapply(seq_len(N), one, mc.cores = cores))
message(slug, " N=", N, " ", format(Sys.time() - t0))
out <- file.path(out_dir, paste0(slug, ".txt"))
writeLines(c(paste("Missing person search", "Missing person search", sep = "\t"), paste("Reference Pedigree", s$label, sep = "\t"),
             paste(format(res[, 1], digits = 6, trim = TRUE), format(res[, 2], digits = 6, trim = TRUE), sep = "\t")), out)
