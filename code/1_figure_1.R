# create guillemot non-breeding distribution map for multistate model chapter

library(tidyverse) # dplyr, ggplot2 etc.
library(sf)        # provides function to read shapefiles
library(ggstar)    # add in star shaped points
library(cowplot)   # better functionality for saving plots
library(magick)

# read data
shp <- st_read("data/landhi_wgs84.shp") # shapefile 
gls <- read.csv("guille_gls_overall_mcp.csv") # gls mcp
rr <- read.csv("guille_rrw_overall_mcp.csv") # ring recovery MCP

# star coords (location of Skomer)
slat <- 51.74
slong <- -5.3
stardata <- data.frame(slat, slong)

####  plot polygon outlines only  ####

# combine data into one dataframe
gls2 <- gls
gls2$id <- "b"
gls2$group <- "b.1"

df <- rbind(rr, gls2)

(lines <- ggplot() +
    geom_sf(data = shp, col = "grey80", fill = "grey80") +
    xlim(-14, 8) + ylim(35, 61) +
    geom_path(data = df, aes(x = long, y = lat, linetype = group), show.legend = F)+
    geom_rect(aes(xmin = -13.06, xmax = 7.15, ymin = 36.15, ymax = 60.38), color = "#D55E00", fill = NA) +
    geom_star(data = stardata, aes(x = slong, y = slat), show.legend = F, size = 2, fill = "#F0E442") +
    scale_linetype_manual(values = c("dashed", "dotted")) +
    theme_light() +
    xlab("Longitude") + ylab("Latitude") +
    theme(axis.title = element_text(size = 10),
          axis.text = element_text(size = 8),
          plot.margin = unit(c(3,2,2,2), "pt")))

(lines_final <- ggdraw() +
    draw_plot(lines) +
    draw_image(image = adult, x = -0.12, y =  -0.32, scale = 0.3) +
    draw_image(image = immature, x = 0.3, y =  -0.2, scale = 0.3) +
    theme(panel.background = element_rect(fill = "white")) +
    panel_border(colour = "white"))

save_plot("GU_mcp_lines_cov_extent_2026.pdf", lines_final, base_width = 2.5, dpi = 600)
save_plot("GU_mcps_lines_cov_extent_2026.png", lines_final, base_width = 2.5, dpi = 600)

# GLS dotted, ring recoveries dashed