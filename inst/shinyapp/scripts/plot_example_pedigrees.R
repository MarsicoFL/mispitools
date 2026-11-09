suppressMessages(library(pedtools))
src <- readLines("scripts/simulate_genetic_examples.R")
eval(parse(text = src[seq_len(grep("^scen <- list", src) + 9)]))
for (slug in names(scen)) {
  s <- scen[[slug]]
  rl <- if (length(s$refs) > 1) paste0("Ref ", seq_along(s$refs)) else "Ref"
  labs <- setNames(c(s$refs, "MP"), c(rl, "MP"))
  png(file.path("www/pedigrees", paste0(slug, ".png")), width = 900, height = 650, res = 150, bg = "white")
  plot(s$ped, labs = labs, hatched = s$refs, fill = list("#d9534f" = "MP"),
       col = list("#d9534f" = "MP"), lwd = 1.6, cex = 1.1, margins = c(4, 1, 1, 1))
  par(xpd = NA)
  legend("bottom", inset = c(0, -0.13), horiz = TRUE, bty = "n", cex = 0.9,
         legend = c("Genotyped reference", "Missing person (MP)", "Not genotyped"),
         fill = c("black", "#d9534f", "white"), density = c(30, -1, 0), angle = 45,
         border = c("black", "#d9534f", "black"))
  dev.off()
  cat(slug, "refs:", s$refs, "\n")
}
