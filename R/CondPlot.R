#' General plot for condiionted probabilities and LR combining variables
#'
#' @param CPT_POP Population conditioned probability table
#' @param CPT_MP Missing person conditioned probability table
#' @export
#' @import reshape2
#' @import patchwork
#' @return A value of Likelihood ratio based on preliminary investigation data. In this case, sex.
#' @examples
#' Cmodel()

CondPlot <- function(CPT_POP, CPT_MP) {

par(mfrow = c(2, 1), mar = c(2, 4, 4, 2))
POP <- reshape2::melt(CPT_POP)

p1 <- ggplot2::ggplot(POP, aes(x = Var2, y = Var1)) +
  ggplot2::geom_raster(aes(fill=value)) +
 scale_fill_gradient(low="grey90", high="blue") +
  labs(x="Hair colour (C)", y="Biologocial sex-Age", title="P(D|H2)", limits = c(0,1)) +
  theme_bw() + theme(axis.text.x=element_text(size=13, angle=0, vjust=0.3),
                     axis.text.y=element_text(size=13),
                     plot.title=element_text(size=13))

MP <- reshape2::melt(CPT_MP)
p2 <- ggplot2::ggplot(MP, aes(x = Var2, y = Var1)) +
  ggplot2::geom_raster(aes(fill=value)) +
  scale_fill_gradient(low="grey90", high="blue", limits = c(0,1)) +
  labs(x="Hair colour (C)", y="Biologocial sex-Age", title="P(D|H1)") +
  theme_bw() + theme(axis.text.x=element_text(size=13, angle=0, vjust=0.3),
                     axis.text.y=element_text(size=13),
                     plot.title=element_text(size=13))
LRtable <- log10(CPT_MP/CPT_POP)
LR <- reshape2::melt(LRtable)
p3 <- ggplot2::ggplot(LR, aes(x = Var2, y = Var1)) +
  ggplot2::geom_raster(aes(fill=value)) +
  scale_fill_gradient(low="grey90", high="blue",limits = c(-5,3)) +
  labs(x="Hair colour (C)", y="Biologocial sex-Age", title="Log10(LR)") +
  theme_bw() + theme(axis.text.x=element_text(size=13, angle=0, vjust=0.3),
                     axis.text.y=element_text(size=13),
                     plot.title=element_text(size=13))
p <- ((p1+p2+p3)) + patchwork::plot_annotation(tag_levels = 'A')
p
return(p)}

