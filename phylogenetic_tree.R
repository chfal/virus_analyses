# load libraries
library(phytools)
library(phyclust)
library(ggtree)
library(phangorn)
library(geiger)
library(treeio)
library(tidyverse)
library(ape)


# set working directory
setwd("C:/Users/chfal/OneDrive - Rutgers University/ATDV_stuff")

# read tree from astral
species_tree <- read.tree("C:/Users/chfal/OneDrive - Rutgers University/ATDV_stuff/inferred_species_tree.nw")

# this is a .csv file that i just made with the old names from astral/orthofinder and the new names that i want them to be
names <- read_csv("names.csv")

# we're going to root the tree at the outgroup, which is white sturgeon
rerooted <- phytools::reroot(species_tree, which(species_tree$tip.label=="white_sturgeon_adv"))

plot(rerooted)

rerooted$node.label
# tree is graphed


my_tree <- ggtree(rerooted) %<+% names +
  xlim(0, 6) +
  geom_tiplab(
    align = TRUE,
    linesize = 0.2,
    offset = 0.03
  )


my_tree


# saving it
# ggsave("species_tree.pdf", my_tree, width=8.47, height=7.15, units="in")

ggsave("species_tree.svg", my_tree, width=8.47, height=7.15, units="in")





# now let's read orthogroups; we get this file from the output of orthofinder, but we need to do a little fixing it because it comes out of orthofinder with counts and we just want it to be 0/1 if it has an orthogroup present or not

# read dataframe
orthogroups <- read_csv("OG000000-23.csv")
orthogroups <- as.data.frame(orthogroups)



orthogroups[] <- lapply(orthogroups, function(x) if (is.numeric(x)) pmin(x, 1) else x)

orthogroups_long <- orthogroups %>%
  pivot_longer(-ID, names_to = "sample", values_to = "presence")


# now we are going to set custom tip order
tip_order <- my_tree$data |>
  subset(isTip) |>
  arrange(y) |>
  pull(label)

tip_order <- rev(tip_order)

# this is our long format and we are reordering it
orthogroups_long$ID <- factor(orthogroups_long$ID, levels = tip_order)


# rename orthogroups

orthogroup_names <- read_csv("orthogroup_names.csv")


orthogroups_long_named <- orthogroups_long %>%
  left_join(orthogroup_names) %>%
  mutate(Name = if_else(is.na(Name), "OG0000023", Name))




# heatmap
heatmap <- ggplot(orthogroups_long_named, aes(x = ID, y = sample, fill = factor(presence))) +
  geom_tile(color = "white") +
  scale_fill_manual(values = c("0" = "white", "1" = "pink")) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  labs(fill = "Presence") +
  scale_y_discrete(labels = orthogroups_long_named$Name)
  

heatmap

ggsave("heatmap.svg", heatmap)




# reading tree once more
# species_tree <- read.tree("inferred_species_tree.nw")

# rooting tree again
# rerooted <- phytools::reroot(species_tree, which(species_tree$tip.label=="white_sturgeon_adv"))

# getting the values for the rerooted trees
rerooted$node.label <- signif(as.numeric(rerooted$node.label), 2)


# making the probability tree
probability_tree <- ggtree(rerooted) %<+% names +
  geom_tree() +
  geom_rootedge(rootedge = 0) +
  geom_tiplab(
    aes(label = new_name),
    align = TRUE,
    offset = 0.3,
    size = 3
  ) +
  geom_nodepoint(
    aes(
      subset = !isTip,
      fill = cut(
        replace(as.numeric(label), is.na(as.numeric(label)), 0),
        breaks = c(-Inf, 0.7, 0.95, Inf),
        labels = c("<.7", "0.7-0.94", ">0.95"),
        include.lowest = TRUE,
        right = FALSE
      )
    ),
    shape = 21,
    size = 4,
    color = "black"
  ) +
  scale_fill_manual(
    values = c(
      "<.7" = "white",
      "0.7-0.94" = "lightpink",
      ">0.95" = "#CC3366"
    ),
    name = "Support",
    na.value = "white"
  ) +
  xlim(0, 5) +
  theme(legend.position = "left")


probability_tree

# ggsave("probability.svg",probability_tree,height=8.5, width=10, units="in")
# ggsave("probability.pdf",probability_tree,height=8.5, width=10, units="in")
