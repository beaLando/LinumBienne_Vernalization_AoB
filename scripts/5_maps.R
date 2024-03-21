# LIBRARIES
library(geosphere)
library(tidyverse)
library(ggplot2)
library(patchwork)
library(ggrepel)

# DATA
dat <- read.csv("./output/exp_summary_pop.csv")
str(dat)

dat <- dat [,-c(1)]
dat$genotyped <- if_else(dat$Pop %in% c("6", "11", "Lla", "Vil", "IOW2", "Sut"), "yes", "no")
dat <- droplevels(subset(dat, species == "bienne"))
str(dat)


# PLOTs

## Get map
world_map <- map_data("world")

mapa <- ggplot(legend = FALSE) +
  geom_polygon(data = world_map, aes(x = long, y = lat, group = group), fill = "lightgray") +
  theme_minimal()

## Map all
m.all <- mapa + 
  coord_quickmap(xlim = c(-14, 20), ylim = c(34, 60)) +
  geom_point(data = dat, size = 4,  
             aes(x = Lon, y = Lat, 
                 fill = genotyped),
             shape = 21) +
  geom_text_repel(data = filter(dat, Lat > 40 | Lon > 10 | Pop == "El Bosque"), size = 5,
            aes(x = Lon, y = Lat, 
                label = Pop)) +
  scale_fill_manual(values = c("white", "cornflowerblue")) +
  labs(x = "Longitude", y = "Latitude", fill = "Genotyped") +
  theme_classic() +
  theme(legend.position = c(0.45, 0.95),
        legend.direction = "horizontal",
        legend.background = element_rect(fill=scales::alpha('white', alpha = 0)),
        text = element_text(size = 20))

## Map inset for southern Spain
dat.inset <- droplevels(filter(dat, Lat < 40 & Lon < 10))

m.inset <- mapa + 
  coord_quickmap(xlim = c(-7, -2.3), ylim = c(35.8, 38.5)) +
  geom_point(data = dat.inset, size = 4,  
             aes(x = Lon, y = Lat, 
                 fill = genotyped),
             shape = 21) +
  geom_text_repel(data = dat.inset, size = 5,
                  aes(x = Lon, y = Lat, 
                      label = Pop)) +
  scale_fill_manual(values = c("white", "cornflowerblue")) +
  labs(x = "Longitude", y = "Latitude") +
  theme_classic() +
  theme(legend.position = "none",
        axis.title = element_blank(),
        axis.text = element_blank(),
        text = element_text(size = 20))


## PC1 x Latitude
dat$genotyped <- factor(dat$genotyped, levels = c("yes", "no"))

dat <- dat %>% arrange(desc(genotyped)) %>% as.data.frame()

pclat <- ggplot() +
  geom_point(data = droplevels(subset(dat, genotyped == "no")),
             aes(x = PC1, y = Lat), 
             shape = 21,
             fill = "white",
             size = 4) +
  geom_point(data = droplevels(subset(dat, genotyped == "yes")),
             aes(x = PC1, y = Lat), 
             shape = 21,
             fill = "cornflowerblue",
             size = 4) +
  ylim(34, 60) +
  labs(x = "climatic PC1", y = "") +
  theme_classic() +
  theme(legend.position = "none",
        text = element_text(size = 20))
  
## Save
together <- m.all + pclat + plot_layout(widths = c(1.5, 1))
ggsave("./output/map&pc.pdf",
       together, 
       width = 15,
       height = 7.5,
       units = "in")

ggsave("./output/map&pc_inset.pdf",
       m.inset)
