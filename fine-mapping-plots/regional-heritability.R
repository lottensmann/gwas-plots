library(ggh4x)
library(RColorBrewer)

#read in fine-mapping results with regional heritability
finemap=read.delim("/Users/ottensma/Documents/second_paper/lipidome_disease/files/finemap_regh2.txt",header = T,sep = "\t",stringsAsFactors  = F)

#define high impact loci as those with regional heritability > 2%
high_impact_loci=unique(finemap_uv_comb[which(finemap_uv_comb$h2g>0.02),]$locus) 

#unique colors for high-impact loci, remaining loci are colored grey
colors=brewer.pal(n = length(high_impact_loci), name = "Set3")

#create factor variable for loci
finemap$locus=factor(finemap$locus,levels=c(high_impact_loci,"other"))

#get number of unique traits per group of trait, needed for size of plot facets
finemap_nod=finemap[!duplicated(finemap$trait),]
group_size=table(finemap_nod$group)

#produce stacked barplot of regional heritability with facets for different trait groups such as lipid classes
options(bitmapType='cairo')
plot=ggplot(finemap, aes(fill=locus, y=h2g, x=trait)) + 
  geom_bar(position="stack", stat="identity",colour="white") + scale_fill_manual(values=colors,name="locus") +
  scale_y_continuous(breaks=c(0,0.05,0.1,0.15,0.2,0.25),limits=c(0,0.26)) + #y-axis scale can be modified according heritability values
  theme(axis.text.x = element_text(angle = 90,size=6),panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),panel.background = element_rect(fill = "white"), 
        plot.background = element_rect(fill = "white", color = NA),legend.position = "n") + labs(x = "trait",y=bquote(''~h^2)) + 
  facet_wrap(~group,scales="free_x",nrow=1) +  force_panelsizes(cols = 0.2*class_siz)
ggsave(plot,filename = "reg_h2_stacked.png",dpi = 300,units="in", height=7, width=15)

