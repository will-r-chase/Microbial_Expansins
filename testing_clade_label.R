tr <- read.tree(text = "((a,((b,c),(d,e))),f);")
ggtree(tr) + geom_tiplab(size = 3) + geom_text2(aes(subset = !isTip, label = node), size = 3, hjust = -0.3)
p<-ggtree(tr) + geom_tiplab(size = 3)
nodes_to_collapse_2<-c(10, 11)
p_collapsed<-collapse(p, 11, clade_name = "Clade 11")
p_collapsed
labels<-c("Clade 10", "Clade 11")

for(i in 1:length(nodes_to_collapse_2)){
  p <- collapse(p, node = nodes_to_collapse_2[i], clade_label=labels[i])
}
p + geom_tiplab(size = 3) + geom_text2(aes(subset = !isTip, label = node), size = 3, hjust = -0.3)

p+geom_point2(aes(subset = (node %in% nodes_to_collapse_2)), size = 3, shape = 23, fill = "blue")

p+geom_nodelab(aes(subset = node==10 & node==11))

p_collapsed+geom_nodelab(aes(subset=node==11))

nodes<-paste(paste0("node==", nodes_to_collapse_2), collapse = "|")
