This repository contains two Python scripts developed during my internship in MAIAGE, INRAE, for *Characterization of TraB in Streptomyces* to analyze AICE13_RLB1-9_QDN99831.1(putative TraB) homologs in Streptomyces and related taxa. The tools work with WHOOPER metadata (TSV output) and provide two main functions:

Pairing AICE13_RLB1-9_QDN99831.1 hits with nearby integrases using genomic coordinates from NCBI RefSeq.

Scanning genomic flanks around each hit for user-defined DNA motifs (with mismatches).

## Input: WHOOPER metadata (metadata.tsv)

WHOOPER identifies homologs of a protein family across RefSeq proteomes. The metadata file contains both hit information and genomic context. Example columns:

| Column                                         | Description                                                       |
| ---------------------------------------------- | ----------------------------------------------------------------- |
| `Target Proteome Accession`                    | RefSeq assembly ID (e.g., `GCF_017639205.1`)                      |
| `Target Proteome`                              | Species and strain name                                           |
| `Query`                                        | Query sequence (e.g., `AICE13_QDN99831.1`)                        |
| `Hit`                                          | RefSeq protein accession (e.g., `WP_189734947.1`)                 |
| `Hit desc`                                     | Functional annotation (e.g., `integrase`, `hypothetical protein`) |
| `Target chr`                                   | Chromosome/contig name                                            |
| `Target start`, `Target stop`, `Target strand` | Genomic coordinates of the hit                                    |
| Taxonomy (`phylum`, `class`, `order`, …)       | Classification info                                               |

## Use:
### Script 1: find_colocated_gene.py

#### Purpose:
Identifies AICE13_RLB1-9_QDN99831.1 proteins and searches for colocated integrase genes in the same genome.

#### How it works:

Reads metadata.tsv (from WHOOPER).

#### For each AICE13_RLB1-9_QDN99831.1 hit:

  Finds integrase hits (from metadata, if present).

  If missing, downloads the corresponding RefSeq genome + annotation (GFF).

  Searches for CDSs annotated as “integrase” within a configurable distance window (default: 10 kb).

  Computes distances between AICE13_RLB1-9_QDN99831.1 and integrase genes (strand-specific if requested).

#### Produces:

  *stats.tsv:* table of AICE13–integrase pairs, with distance and annotation.

  *stats.txt:* summary counts.

  *genomes/:* cached .fna.gz and .gff.gz files from NCBI.

#### Example usage:

     python find_colocated_gene.py \
     --input metadata.tsv \
     --output results.tsv \
     --stats stats.tsv \
     --distance 10000 \
     --workdir genomes \
     --strand-specific

### Script 2: find_motif_downstream.py

#### Purpose:
Scans upstream and downstream regions of each AICE13_RLB1-9_QDN99831.1 hit for a user-defined DNA motif, allowing mismatches.

#### How it works:

Uses stats.tsv from step 1 to locate each AICE13_RLB1-9_QDN99831.1 gene.

#### Extracts the coding sequence and flanks (--upstream-size, --downstream-size).

Applies a sliding window with Hamming distance to count motif matches with up to --max-mismatches.

#### Reports:

  Motif counts upstream and downstream.
  Extracted sequences and coordinates.
  Total motif counts per element.

#### Example usage:
    python find_motif_downstream.py \
    --stats stats.tsv \
    --genome-dir genomes \
    --upstream-size 200 \
    --downstream-size 200 \
    --motif ATGCGTN \
    --max-mismatches 1 \
    --output motif_scan.tsv

## Typical workflow
1. Run WHOOPER to get metadata.tsv with AICE13 homologs.

2. Use find_colocated_gene.py to identify colocated integrases → stats.tsv.

3. Use find_motif_downstream.py to scan flanks for motifs → motif_scan.tsv.

4. Perform downstream comparative or phylogenetic analysis.
