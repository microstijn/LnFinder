# LnFinder

LnFinder.jl is a Julia package designed to identify and classify Lanthanide-dependent enzymes from large-scale metagenomic datasets.

It operationalizes the biochemical "Lanthanide Switch" described by Voutsinos et al. (2025), allowing for the (eventual) rapid screening of PQQ-dependent dehydrogenases in diverse ecosystems (soil, freshwater, marine) to determine if they encode for Lanthanide (Ln) or Calcium (Ca) as a cofactor.

A recent paper by Voutsinos et al. (2025) proposes that Lanthanide-dependent metabolism is far more ubiquitous than previously thought, particularly in the ocean. The distinction between Lanthanide-dependent and Calcium-dependent PQQ-alcohol dehydrogenases (PQQ-ADHs) can be predicted by a specific amino acid motif in the active site:

> Lanthanide-dependent: Contains an additional Aspartate residue essential for Ln coordination (D-x-D).
> Calcium-dependent: Contains an Alanine, Threonine, or Serine in the equivalent position (D-x-[ATS]).

LnFinder automates the detection of these motifs in aligned sequences.

## Features

**Motif classification**: Distinguishes between Lanthanide, Calcium, and Unknown motifs using PROSITE notation.

**Streaming analysis**: Optimized analyze_alignment function handles massive FASTA files by streaming data rather than loading it entirely into RAM.

**Heuristic column detection**: Automatically identifies the consensus active-site column by using ASLIGNED sequences, sampling the initial subset of the alignment.
