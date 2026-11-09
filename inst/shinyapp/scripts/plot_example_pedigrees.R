suppressMessages(library(pedtools))
src <- readLines("scripts/simulate_genetic_examples.R")
eval(parse(text = src[seq_len(grep("^scen <- list", src) + 9)]))
for (slug in names(scen)) {
  s <- scen[[slug]]
  rl <- if (length(s$refs) > 1) paste0("Ref ", seq_along(s$refs)) else "Ref"
  labs <- setNames(c(s$refs, "MP"), c(rl, "MP"))
  png(file.path("www/pedigrees", paste0(slug, ".png")), width = 900, height = 600, res = 150, bg = "white")
  plot(s$ped, labs = labs, hatched = s$refs, fill = list("#CC79A7" = "MP"),
       col = list("#CC79A7" = "MP"), lwd = 1.6, cex = 1.1, margins = c(1, 1, 1, 1))
  dev.off()
  cat(slug, "refs:", s$refs, "\n")
}
