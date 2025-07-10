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
  "QWE88920.1 surface glycoprotein [Severe acute respiratory syndrome coronavirus 2]" = "Alpha", 
  "QWK65230.1 surface glycoprotein [Severe acute respiratory syndrome coronavirus 2]" = "Beta",  
  "UFO69279.1 surface glycoprotein [Severe acute respiratory syndrome coronavirus 2]" = "Delta", 
  "UVC22721.1 surface glycoprotein [Severe acute respiratory syndrome coronavirus 2]" = "Omicron",
  "WZD59850.1 surface glycoprotein [Severe acute respiratory syndrome coronavirus 2]" = "Omicron JN.1"
)

names(sars2_spike_sequences) <- variant_names[names(sars2_spike_sequences)]

# Optional reorder to match name order above.
desired_order <- c("Wuhan", "Alpha", "Beta", "Delta", "Omicron", "Omicron JN.1")
sars2_spike_sequences <- sars2_spike_sequences[desired_order]

# Optional naming
names(sars2_spike_sequences) <- c(
  "Omicron JN.1 | WZD59850.1",
  "Omicron | UVC22721.1",
  "Delta | UFO69279.1",
  "Alpha | QWE88920.1",
  "Beta | QWK65230.1",
  "Wuhan | YP_009724390.1"
)

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
