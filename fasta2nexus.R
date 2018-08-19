library(seqinr)
library(ape)

fasta <- read.fasta("microbe_expansin_alignment_removeseqs.fas")
write.nexus.data(fasta, file = "microbe_expansin_alignment_removeseqs_nex.nexus", format = "protein")
