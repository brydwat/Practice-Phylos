library(BiocManager)
library(tidyverse)
library(msa)
library(tinytex)
library(Biostrings)


# Simple MSA test on Sars-CoV-2 spike proteins.
# 7/1/2025
# brydwat

sars2_spike_sequences <- readAAStringSet("spike_proteins.fasta")

names(sars2_spike_sequences) <- c(
  "Alpha | QWE88920.1",
  "Beta | QWK65230.1",
  "Delta | UFO69279.1",
  "Wuhan | YP_009724390.1"
)

spike_alignment <- msa(sars2_spike_sequences)
spike_alignment



msaPrettyPrint(spike_alignment, output="pdf", showNames="left",
               showLogo="none", askForOverwrite=FALSE, verbose=FALSE)
