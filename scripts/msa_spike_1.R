library(BiocManager)
library(tidyverse)
library(msa)
library(tinytex)
library(Biostrings)
library(ape)
library(seqinr)


# Simple MSA test on Sars-CoV-2 spike proteins.
# 7/1/2025
# brydwat

sars2_spike_sequences <- readAAStringSet("spike_proteins.fasta")

variant_names <- c(
  "YP_009724390.1 surface glycoprotein [Severe acute respiratory syndrome coronavirus 2]" = "Wuhan",
  "ULQ63197.1 surface glycoprotein [Severe acute respiratory syndrome coronavirus 2]" = "Alpha", 
  "QWW93436.1 surface glycoprotein [Severe acute respiratory syndrome coronavirus 2]" = "Beta",  
  "WPW42754.1 surface glycoprotein [Severe acute respiratory syndrome coronavirus 2]" = "Delta", 
  "WPW42837.1 surface glycoprotein [Severe acute respiratory syndrome coronavirus 2]" = "Omicron BA.1",
  "WZD59850.1 surface glycoprotein [Severe acute respiratory syndrome coronavirus 2]" = "Omicron JN.1"
)

names(sars2_spike_sequences) <- variant_names[names(sars2_spike_sequences)]

spike_alignment <- msa(sars2_spike_sequences)
spike_alignment


# Print string results to PDF
msaPrettyPrint(spike_alignment, 
               output="pdf", 
               showNames="left",
               showLogo="none", 
               askForOverwrite=FALSE, 
               verbose=FALSE)

# Compute distance matrix
spike_alignment_sequence <- msaConvert(spike_alignment, type = "seqinr::alignment")

distance_alignment <- dist.alignment(spike_alignment_sequence)

# Compute phylogenetic tree using neighbor joining

spike_tree <- bionj(distance_alignment)

plot(spike_tree)
