# code to create figure 3 - cumulative frequency plots of first return and first breed
library(tidyverse)
library(ggplot2)
library(scales)
library(cowplot)

#### read data  ####
dat <- read.csv("data/firstbreedetc_reduced.csv") 

# table includes ringing occasion, occasion an individual was first seen at the colony, 
# the occasion an individual first bred at the colony, and the age at return and first breeding

return <- dat$firstreturn
return <- return[!is.na(return)]

breed <- dat$firstbreed
breed <- breed[!is.na(breed)]

####  plot age of first breeding  ####

# 5971 individuals
# 2239 are ever resighted at the colony
# 628 individuals seen breeding

breeders <- dat %>% drop_na()

range(breeders$breedage)

# cumulative relative frequency - plot as a proportion

# calculate proportion breeding at each age
p <- breeders$breedage
breaks <- seq(3, 20, by = 1)
p.cut <- cut(p, breaks, right = F)
p.freq <- table(p.cut)
p.cumfreq <- cumsum(p.freq)
p.cumrelfreq <- p.cumfreq/nrow(breeders)

cumrelfreq.breeders <- as.data.frame(cbind(breaks, p.cumrelfreq))

b1 <- function() {
  par(mar = c(2, 3, 0.5, 0.5))
  hist(breeders$breedage, breaks = seq(0, 20, by = 1), main = NA, xlab =" ", ylab = "", 
       col = alpha("#4B0055", 0.5), axes = F)
  axis(2, ylim=c(0, 175),las=1, tck = -0.01, cex.axis = 0.7, mgp = c(2, 0.3, 0))
  # time axis
  axis(1, xlim = c(0, 20), mgp=c(2, 0, 0), tck = -0.01, cex.axis = 0.7, at = c(0, 3, 6, 9, 12, 15, 18))
  mtext("Age at first breeding", side=1, line=0.8, cex = 0.8)
  par(new=TRUE)
  plot(cumrelfreq.breeders$breaks, cumrelfreq.breeders$p.cumrelfreq, pch=NA, axes=FALSE, 
       ylim=c(0, 1), xlim=c(0, 18), xlab="", ylab="", type="o", lty = 2, col = "black")
  axis(2, ylim=c(0, 1),las=1, tck = 0.01, labels = F, cex.axis = 0.7, col = "black")
  text(0, seq(0, 1, 0.2), labels = c(" ", 0.2, 0.4, 0.6, 0.8, 1), adj = c(0.5, 0.5), cex = 0.7)
  box(bty="l")
}

ggdraw(b1)

#### plot age of first return  ####

# format data
return <- dat %>% drop_na(returnage)

hist(return$returnage, breaks = seq(0, 15, by = 1))

# calculate proportion returning at each age
r <- return$returnage
breaks <- seq(1, 15, by = 1)
r.cut <- cut(r, breaks, right = F)
r.freq <- table(r.cut)
r.cumfreq <- cumsum(r.freq)
r.cumrelfreq <- r.cumfreq/nrow(return)

cumrelfreq.return <- as.data.frame(cbind(breaks, r.cumrelfreq))

r1 <- function() { 
  par(mar = c(2, 3, 0.5, 0.5))
  hist(return$returnage, breaks = seq(0, 15, by = 1), main = NA, xlab =" ", ylab = "", 
       col = alpha("#4B0055", 0.5), axes = F)
  axis(2, ylim=c(0, 1125),las=1, tck = -0.01, cex.axis = 0.7, mgp = c(2, 0.3, 0))
  mtext("Frequency", side = 2, line = 2, cex = 0.8)
  # time axis
  axis(1, xlim = c(0, 15), mgp=c(2, 0, 0), tck = -0.01, cex.axis = 0.7, at = c(0, 3, 6, 9, 12, 15))
  mtext("Age at first return", side=1, line=0.8, cex = 0.8)
  par(new=TRUE)
  plot(cumrelfreq.return$breaks, cumrelfreq.return$r.cumrelfreq, pch=NA, axes=FALSE, 
       ylim=c(0, 1), xlim=c(0, 13), xlab="", ylab="", type="o", lty = 2, col = "black")
  axis(2, ylim=c(0, 1),las=1, tck = 0.01, labels = F, cex.axis = 0.7, col = "black")
  text(0, seq(0, 1, 0.2), labels = c(" ", 0.2, 0.4, 0.6, 0.8, 1), adj = c(0.5, 0.5), cex = 0.7)
  mtext("Cumulative relative frequency", side = 2, line = 1.3, cex = 0.8, col = "black")
  box(bty="l")
}

ggdraw(r1)

####  combine into one plot with cowplot  ####

(plot_grid(r1, NULL, b1, rel_widths = c(1, -0.1,  1), nrow = 1))

(cumfreqs <- plot_grid(r1, NULL, b1, NULL, rel_widths = c(1, -0.08,  1, 0.05), nrow = 1, 
                       labels = c("A", " ", "B", " "), label_x = 0.05))

(plot_to_save <- ggdraw() +
    draw_plot(cumfreqs) +
    theme(panel.background = element_rect(fill = "white")) +
    panel_border(colour = "white"))

save_plot("cumfreq_return_breed_1995-2019_bw.png", plot_to_save, base_width = 7.4, dpi = 600)
save_plot("cumfreq_return_breed_1995-2019_bw.pdf", plot_to_save, base_width = 7.4, dpi = 600)
