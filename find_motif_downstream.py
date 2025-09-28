#!/usr/bin/env python3
"""
Script: find_motif_downstream.py

Given a stats.tsv of AICE13 hits and integrase annotations, extract both upstream
and downstream regions of each AICE13 hit (plus the AICE13 sequence itself) from its
genome assembly, search for a user-specified DNA motif allowing a given number of
mismatches, and write out a result.tsv with counts, sequences, and coordinates.

Usage:
    python find_motif_flanks.py \
        --stats stats.tsv \
        --genome-dir genomes \
        --upstream-size 100 \
        --downstream-size 100 \
        --motif ATGCGTN \
        --max-mismatches 1 \
        --output result.tsv

Output:
    result.tsv (tab-delimited) with original columns + AICE13_Seq + \
    Upstream_Motif_Count + Downstream_Motif_Count + Total_Motif_Count + \
    Upstream_Seq + Downstream_Seq + Upstream_Coords + Downstream_Coords
"""
import os
import sys
import csv
import gzip
import argparse
from Bio import SeqIO


def parse_args():
    parser = argparse.ArgumentParser(
        description="Search for a DNA motif upstream and downstream of AICE13 hits")
    parser.add_argument(
        "--stats", required=True,
        help="Path to stats.tsv from previous integrase detection script")
    parser.add_argument(
        "--genome-dir", dest="genome_dir", default="genomes",
        help="Directory containing genome FASTA (.fna/.fna.gz) and GFF (.gff.gz) files named by assembly accession")
    parser.add_argument(
        "--upstream-size", type=int, required=True,
        help="Number of bases upstream of each AICE13 hit to search")
    parser.add_argument(
        "--downstream-size", type=int, required=True,
        help="Number of bases downstream of each AICE13 hit to search")
    parser.add_argument(
        "--motif", required=True,
        help="DNA motif to search for (e.g. ATGCGT)")
    parser.add_argument(
        "--max-mismatches", type=int, default=0,
        help="Maximum number of mismatches allowed when matching the motif")
    parser.add_argument(
        "--output", default="result.tsv",
        help="Path to output TSV file")
    return parser.parse_args()


def load_fasta_for_assembly(assembly, genome_dir):
    # Build paths for gzipped and plain FASTA files
    gz = os.path.join(genome_dir, f"{assembly}.fna.gz")
    fa = os.path.join(genome_dir, f"{assembly}.fna")
    # Choose the existing file (gzipped preferred), else bail out
    if os.path.isfile(gz):
        opener = lambda p: gzip.open(p, 'rt')
        path = gz
    elif os.path.isfile(fa):
        opener = lambda p: open(p, 'r')
        path = fa
    else:
        # No FASTA found for this assembly
        return None, None
    # Parse the FASTA into a dict of SeqRecord objects
    # Key = full header ID string (everything after '>')
    with opener(path) as handle:
        records = SeqIO.to_dict(SeqIO.parse(handle, 'fasta'))
    # Build an alias map: strip off any " .version" and description text
    alias = {rid.split()[0].split('.')[0]: rid for rid in records}
    return records, alias

#find the protein in gff
def find_cds_by_protein_id(gff_path, protein_id):
    with gzip.open(gff_path, 'rt') as gff:
        for line in gff:
            if line.startswith('#'): continue
            parts = line.strip().split('\t')
            if len(parts) != 9: continue
            seqid, _, ftype, start, end, _, strand, _, attr = parts
            if ftype == 'CDS' and f"protein_id={protein_id}" in attr:
                return seqid, int(start), int(end), strand
    return None


def count_motif_hits(seq, motif, max_mismatches):
    s = str(seq).upper()
    m = motif.upper()
    #Window width is the length of the motif
    w = len(m)
    count = 0
    #Compare the window to the motif position by position
    #count how many positions differ (Hamming distance).
    #If that mismatch count is ≤ allowed, it’s a hit.
    
    for i in range(len(s) - w + 1):
        window = s[i:i+w]
        if sum(a!=b for a, b in zip(window, m)) <= max_mismatches:
            count += 1
    return count


def extract_flank_coords(seq_len, start, end, strand, size, direction):
    if direction == 'downstream':
        if strand == '+':
            a, b = end + 1, end + size
        else:
            a, b = start - size, start - 1
    else:
        if strand == '+':
            a, b = start - size, start - 1
        else:
            a, b = end + 1, end + size
    a = max(1, a)
    b = min(seq_len, b)
    return a, b


def extract_seq(seq, a, b, strand):
    sub = seq[a-1:b]
    return sub.reverse_complement() if strand == '-' else sub


def main():
    args = parse_args()
    with open(args.stats) as f:
        reader = csv.DictReader(f, delimiter='\t')
        rows = list(reader)
    new_cols = [
        'AICE13_Seq',
        'Upstream_Motif_Count','Downstream_Motif_Count','Total_Motif_Count',
        'Upstream_Seq','Downstream_Seq',
        'Upstream_Coords','Downstream_Coords'
    ]
    out_fields = reader.fieldnames + new_cols
    #Two caches so you only load each genome FASTA or GFF once per assembly
    fasta_cache = {}
    gff_cache = {}
    results = []
    for row in rows:
		#assembly accession string from current row
        asm = row['Target_Proteome']
        #Avoid re-loading the same genome file over and over when multiple hits belong to the same assembly.
        if asm not in fasta_cache:
            fasta_cache[asm] = load_fasta_for_assembly(asm, args.genome_dir)
        records, alias = fasta_cache[asm]
        # If no FASTA was found, mark all new columns as NA and skip
        if not records:
            for col in new_cols:
                row[col] = 'NA'
            results.append(row)
            continue
        # Cache the GFF path (or None if missing) to avoid repeated filesystem checks
        if asm not in gff_cache:
            p = os.path.join(args.genome_dir, f"{asm}.gff.gz")
            gff_cache[asm] = p if os.path.isfile(p) else None
        gff_path = gff_cache[asm]
        if gff_path:
			# Lookup CDS coords by protein_id in the GFF
            info = find_cds_by_protein_id(gff_path, row['AICE13_Hit'])
            if info:
                contig, a_start, a_end, strand = info
            else:
				# protein_id not in GFF > cannot locate element
                for col in new_cols:
                    row[col] = 'NA'
                results.append(row)
                continue
        else:
            # No GFF available → parse locus string like "contigX:100-500"
            try:
                contig, coords = row['AICE13_Locus'].split(':')
                a_start, a_end = map(int, coords.split('-'))
                strand = row.get('AICE13_Strand','+')
            except:
                for col in new_cols:
                    row[col] = 'NA'
                results.append(row)
                continue
        # If contig name needs version-stripped lookup, use alias map
        if contig not in records:
            contig = alias.get(contig)
        # If still not found, mark as NA and skip
        if not contig or contig not in records:
            for col in new_cols:
                row[col] = 'NA'
            results.append(row)
            continue
        seq = records[contig].seq
        aice_seq = extract_seq(seq, a_start, a_end, strand)
        up_a, up_b = extract_flank_coords(len(seq), a_start, a_end, strand, args.upstream_size, 'upstream')
        down_a, down_b = extract_flank_coords(len(seq), a_start, a_end, strand, args.downstream_size, 'downstream')
        up_seq = extract_seq(seq, up_a, up_b, strand)
        down_seq = extract_seq(seq, down_a, down_b, strand)
        up_count = count_motif_hits(up_seq, args.motif, args.max_mismatches)
        down_count = count_motif_hits(down_seq, args.motif, args.max_mismatches)
        total = up_count + down_count
        row['AICE13_Seq'] = str(aice_seq)
        row['Upstream_Motif_Count'] = str(up_count)
        row['Downstream_Motif_Count'] = str(down_count)
        row['Total_Motif_Count'] = str(total)
        row['Upstream_Seq'] = str(up_seq)
        row['Downstream_Seq'] = str(down_seq)
        row['Upstream_Coords'] = f"{contig}:{up_a}-{up_b}"
        row['Downstream_Coords'] = f"{contig}:{down_a}-{down_b}"
        results.append(row)
    with open(args.output, 'w') as out_f:
        writer = csv.DictWriter(out_f, out_fields, delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(results)
    print(f"Results written to {args.output}")

if __name__ == '__main__':
    main()
