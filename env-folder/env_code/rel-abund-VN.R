



# community data aggregation at a taxonomic level. e.g. genus
# take the sum of each phylum in each sample
taxa.agg <- aggregate(t(comm_rarfy),
                      by = list(taxo.sub[colnames(comm_rarfy),6]),
                      FUN = sum)
# clean up resulting object
rownames(taxa.agg) <- taxa.agg$Group.1
phylum_data <- t(taxa.agg[,-1])
# convert abundances to relative abundances
phylum_data <- phylum_data/rowSums(phylum_data)
# remove rare phyla
phylum_data <- phylum_data[,colSums(phylum_data)>0.01]
# now reshape phylum data to long format
phylum_data <- reshape2::melt(phylum_data)
# rename columns
colnames(phylum_data)[1:2] <- c('Samples','Bacteria_genus')

# Add column
phylum_data$env <- rep(c(rep("desert", 7), rep("urban", 4)), 5)

# now we can plot phylum relative abundance per sample
 ggplot(phylum_data,
        aes(Samples,
            weight = value,
            fill = Bacteria_genus)) +
  geom_bar(color = "black",
           width = .7,
           position = 'fill') +
  labs(y = "\nRelative abundance (%)") +
  scale_y_continuous(expand = c(0,0)) +
  scale_fill_viridis_d(direction = -1L) +
  theme_bw() +
  theme(axis.title.y = element_blank(),
        panel.spacing = unit(1, "cm",
                             data = NULL)) +
  coord_flip() +
  facet_grid(. ~ env)