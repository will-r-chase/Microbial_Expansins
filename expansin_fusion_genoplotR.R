library(genoPlotR)
library(plyr)
library(dplyr)
library(ape)
library(phytools)
library(ade4)


#read domain info dataframe
domain_data<-read.table("domain_positions_genoplotR.txt", sep="\t", stringsAsFactors = F, header=T)

#add strand column
domain_data$strand<-"+"

#reorder columns for dna_segs function
domain_data<-data.frame(name=domain_data$domain, 
                        start=domain_data$start, 
                        end=domain_data$end,
                        strand=domain_data$strand,
                        species=domain_data$name, 
                        accession=domain_data$id
                        )

#add colors, fill, and gene_type columns
domain_data$col<-NA
domain_data$col[which(domain_data$name=="Full Protein")]<-"black"
domain_data$col[which(domain_data$name=="Full Protein 1")]<-"black"
domain_data$col[which(domain_data$name=="Full Protein 2")]<-"black"
domain_data$col[which(domain_data$name=="Expansin")]<-"steelblue"
domain_data$col[which(domain_data$name=="CBM1")]<-"firebrick1"
domain_data$col[which(domain_data$name=="CBM2")]<-"mediumorchid1"
domain_data$col[which(domain_data$name=="CBM32")]<-"lightpink2"
domain_data$col[which(domain_data$name=="GH9")]<-"yellow3"
domain_data$col[which(domain_data$name=="GH5")]<-"seagreen2"
domain_data$fill<-domain_data$col
domain_data$fill[which(domain_data$name=="Full Protein")]<-"transparent"
domain_data$fill[which(domain_data$name=="Full Protein 1")]<-"transparent"
domain_data$fill[which(domain_data$name=="Full Protein 2")]<-"transparent"
domain_data$gene_type<-"blocks"
domain_data$gene_type[which(domain_data$name=="Full Protein")]<-"arrows"
domain_data$gene_type[which(domain_data$name=="Full Protein 1")]<-"arrows"
domain_data$gene_type[which(domain_data$name=="Full Protein 2")]<-"arrows"


#reorder rows so "full protein" line is always in front of domain boxes in plot
order<-c("Expansin", "GH5", "GH9", "CBM1", "CBM2", "CBM32", "Full Protein", "Full Protein 1", "Full Protein 2")
domain_data<-domain_data %>% arrange(species, match(name, order))

#split plots by species, make into dna_segs object, and plot gene map
domains_split<-split(domain_data, domain_data$species)
dna_segs<-lapply(domains_split, function(x) dna_seg(x))

#read tree and rename
tree <- read.tree("newseqs_combined_treegraph_rooted_fixed.nwk")
ult.tree <- force.ultrametric(tree)
ult.tree<-read.tree("ultrametric_tree.nwk")
#read from file with replacement tip names to match domain data
new_tip_names<-readClipboard()
ult.tree$tip.label<-new_tip_names
ult.tree<-drop.tip(ult.tree, "Streptomyces acidiscabies2")
tr_hclust<-as.hclust.phylo(ult.tree)
phylog<-hclust2phylog(tr_hclust)

#read new tip names from clipboard
new_phylog_names <- readClipboard()
names(phylog$leaves) <- new_phylog_names
dna_segs2 <- dna_segs[names(phylog$leaves)]

tr_hclust<-as.hclust.phylo(ult.tree)
phylog<-hclust2phylog(tr_hclust)
phylog_names<-names(phylog$leaves)
names(dna_segs2)<-phylog_names

plot_gene_map(dna_segs=dna_segs2, tree=phylog, tree_width = 8, legend = T)




###
