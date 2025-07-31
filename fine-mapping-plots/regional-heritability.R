# =============================================================================
# Fine-mapping Regional Heritability Plot
# 
# Description:
#   Creates a stacked barplot of regional heritability (h²) per trait, 
#   high-impact loci (h² > 2%) are colored and the plot is faceted by trait groups.
#
# Output:
#   reg_h2_stacked.png
# =============================================================================


suppressPackageStartupMessages({
  library(ggplot2)
  library(ggh4x)
  library(RColorBrewer)
  library(dplyr)
})

#Read in fine-mapping results with regional heritability
finemap=read.delim("/Users/ottensma/Documents/second_paper/lipidome_disease/files/finemap_regh2.txt",header = T,sep = "\t",stringsAsFactors  = F)

#Define high-impact loci as those with regional heritability > 2%
high_impact_loci=unique(finemap_uv_comb[which(finemap_uv_comb$h2g>0.02),]$locus) 

#Label loci not among high-impact loci as other loci and create factor variable
finemap <- finemap %>%
  mutate(locus_grouped = ifelse(locus %in% high_impact_loci, locus, "other"),
         locus_grouped = factor(locus_grouped, levels = c(high_impact_loci, "other")))

#Assign unique colors for high-impact loci, remaining loci are colored grey
n_colors <- length(high_impact_loci)
colors <- c(brewer.pal(n_colors, "Set3"), "grey") #NOTE: for more than 12 high-impact loci, unique colors need to be assigned in a different way
names(colors) <- c(high_impact_loci, "other")

#Get panel size weights
trait_group_counts <- finemap %>%
  distinct(trait, group) %>%
  count(group, name = "n_traits")
finemap <- finemap %>%
  left_join(trait_group_counts, by = "group")

#Produce stacked barplot of regional heritability with facets for different trait groups such as lipid classes
options(bitmapType='cairo')
plot=ggplot(finemap, aes( x=trait, y=h2g, fill=locus_grouped)) + 
  geom_bar(position="stack", stat="identity",colour="white") + 
  scale_fill_manual(values=colors,name="locus") +
  scale_y_continuous(breaks = seq(0, 0.25, by = 0.05),limits=c(0,0.26)) + #y-axis scale can be modified according heritability values
  facet_wrap(~group,scales="free_x",nrow=1) + 
  force_panelsizes(cols = 0.2*trait_group_counts$n_traits) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, size = 6),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "right",
    plot.background = element_rect(fill = "white", color = NA),
    panel.spacing = unit(1, "lines")
    ) + 
  labs(x = "trait",y=bquote(''~h^2),title = "Regional heritability by trait") 
ggsave(plot,filename = "reg_h2_stacked.png",dpi = 300,units="in", height=7, width=15)

