library(ggplot2)
library(gggenes)
library(readr)
library(dplyr)
library(stringr)

# set working directory
setwd("C:/Users/chfal/OneDrive - Rutgers University/ATDV_stuff")

# read all data in
genes <- read_csv("gggenes_input.csv")
blast <- read_csv("blast_results_edited.csv")
names <- read_csv("names.csv")

# merge the data together
merged <- genes %>%
  left_join(blast)

merged[[2]] <- ifelse(
  !is.na(merged[[ncol(merged)]]) & merged[[ncol(merged)]] != "NA" & merged[[ncol(merged)]] != "",
  merged[[ncol(merged)]],
  merged[[2]]
)

# remove the ID
merged <- merged %>%
  select(-ID)

# now filter 8054 and 8158 only
genes_8054 <- merged %>%
  filter(molecule==8054 | molecule == 8158)

# remove prefixes
genes_8054$gene <- gsub("^8054_", "", genes_8054$gene)

# remove prefixes
genes_8054$gene <- gsub("^8158_", "", genes_8054$gene)



# i think we accidentally reverse complemented it compared to the others so let's try and undo that

only_8054 <- genes_8054 %>%
  filter(molecule == 8054)
  

only_8158 <- genes_8054 %>%
  filter(molecule == 8158)



genome_length_8054 <- 37517

genome_length_8158 <- 37535

rc_8054 <- only_8054 %>%
  mutate(
    start_new = genome_length_8054 - end + 1,
    end_new   = genome_length_8054 - start + 1,
    strand_new = ifelse(strand == "+", "-", "+")
  ) %>%
  select(molecule, gene, start_new, end_new, strand_new, orientation) %>%
  rename(start = start_new) %>%
  rename(end =end_new) %>%
  rename(strand=strand_new)

rc_8054$new_name <-"Anolis-Adv-2 (8054)"


rc_8158 <- only_8158 %>%
  mutate(
    start_new = genome_length_8158 - end + 1,
    end_new   = genome_length_8158 - start + 1,
    strand_new = ifelse(strand == "+", "-", "+")
  ) %>%
  select(molecule, gene, start_new, end_new, strand_new, orientation) %>%
  rename(start = start_new) %>%
  rename(end =end_new) %>%
  rename(strand=strand_new)

rc_8158$new_name <-"Anolis-Adv-2 (8158)"


rc_8054_8158 <- rbind(rc_8054,rc_8158)

# set the color to the gene order

gene_order <- unique(rc_8054_8158$gene[order(rc_8054_8158$start)])

rc_8054_8158$gene <- factor(rc_8054_8158$gene, levels = gene_order)


full_gggenes <- ggplot(rc_8054_8158, aes(xmin = start, xmax = end, y = molecule, fill = gene)) +
  geom_gene_arrow(aes(forward = strand == "+"),arrowhead_height = unit(0.5,"in")) + facet_wrap(~ molecule, scales = "free", ncol = 1) +
  theme_bw() +
  theme(axis.title.y = element_blank()) +
  scale_y_discrete(labels = c("Anolis-Adv-2 (8054)" = "", "Anolis-Adv-2 (8158)" = "")) +
  theme(
    axis.text.y = element_blank(),  # Removes the text numbers
    axis.ticks.y = element_blank()  # Removes the small tick lines
  )

full_gggenes



ggsave("full_gggenes.pdf", full_gggenes, width=9.5, height=5)
ggsave("full_gggenes.png", full_gggenes,width=9.5, height=5)


# filter out all species-specific unidentified ORFs
filtered_merged <- merged %>%
  filter(!str_detect(gene, fixed(molecule))) %>%
  filter(!grepl('Y', gene)) %>%
  filter(!grepl('ORF1', gene))
  

# join the new name we want it to be
filtered_merged_names <- left_join(filtered_merged, names, by=c("molecule"="old_name"))

# reverse complement filtered
filtered_merged_rc <- filtered_merged_names %>%
  filter(!molecule==8054) %>%
  filter(!molecule==8158)


filtered_merged_rc_final <- rbind(filtered_merged_rc,rc_8054_8158) %>%  filter(!grepl("^\\d+$", gene))

gene_order_filtered_merged <- unique(filtered_merged_rc_final$gene[order(filtered_merged_rc_final$start)])

filtered_merged_rc_final$gene <- factor(filtered_merged_rc_final$gene, levels = gene_order_filtered_merged)


filtered_merged_gggenes <- ggplot(filtered_merged_rc_final, aes(xmin=start, xmax=end, y=new_name, fill= gene)) +
  geom_gene_arrow(aes(forward = strand == "+"),arrowhead_height = unit(0.5,"in")) + facet_wrap(~ new_name, scales = "free", ncol = 1) +
  theme_bw() +
  theme(axis.title.y = element_blank()) +
  xlim(0,40000) +
  theme(
    axis.text.y = element_blank(),  # Removes the text numbers
    axis.ticks.y = element_blank()  # Removes the small tick lines
  )


ggsave("filtered_merged_gggenes.svg", filtered_merged_gggenes, height=50, width=20, units="in", limitsize=FALSE)


ggsave("filtered_merged_gggenes.pdf", filtered_merged_gggenes, height=50, width=20, units="in", limitsize=FALSE)


filtered_even_more <- filtered_merged_rc_final %>%
   filter(molecule=="8054" | molecule=="8158" | molecule =="lizard_adv_2" | molecule == "barthadenovirus_zootocae" | molecule == "bearded_dragon_adv_1" | molecule == "deer_adv_A" | molecule == "tern_adv_1" | molecule == "barthadenovirus_varani" | molecule =="psittacine_adv_3")



ggplot(filtered_even_more, aes(xmin=start, xmax=end, y=new_name, fill= gene)) +
  geom_gene_arrow(aes(forward = strand == "+"),arrowhead_height = unit(0.25,"in")) + facet_wrap(~ new_name, scales = "free", ncol = 1) +
  theme_bw() +
  theme(axis.title.y = element_blank()) +
  xlim(0,40000) +
  theme(
    axis.text.y = element_blank(),  # Removes the text numbers
    axis.ticks.y = element_blank()  # Removes the small tick lines
  )

ggsave("filtered_even_more_plot.pdf", width = 7, height = 11, units="in")

ggsave("filtered_even_more_plot.svg", width = 7, height = 11, units="in")



filtered_22_orthogroups <- filtered_even_more %>%
  filter(gene%in%c("SPIKE","CAP6","CAPSP","CAP3","PKG3","TERM","DPOL","SHUT","DNB2","PRO","CAPSH","CAP8","PKG1","E43","NP","RH1","E41","E42","LH3","L2MU","P32K","UXP","LH2"))

ggplot(filtered_22_orthogroups, aes(xmin=start, xmax=end, y=new_name, fill= gene)) +
  geom_gene_arrow(aes(forward = strand == "+"),arrowhead_height = unit(0.25,"in")) +
  facet_wrap(~ new_name, scales = "free",ncol=1) +
  theme(axis.title.y = element_blank()) +
  xlim(0,40000) +
  theme(
    axis.text.y = element_blank(),  # Removes the text numbers
    axis.ticks.y = element_blank()  # Removes the small tick lines
  )


ggsave("filtered_22_orthogroups.pdf", width = 7, height = 11, units="in")

ggsave("filtered_22_orthogroups.svg", width = 7, height = 11, units="in")


         
