import os
import sys
import csv
import gzip
import json
import time
import shutil
import argparse
import requests
from io import StringIO
from collections import defaultdict
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# ------------------ PARAMETERS ------------------
REFSEQ_SUMMARY_URL = "https://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_refseq.txt"
MAX_RETRIES = 3
USER_AGENT = "AICE13Analyzer/1.0 (lzlodhi@inrae.com)"
DEFAULT_MAX_DISTANCE = 10000

# ------------------ UTILITY FUNCTIONS ------------------

# ftp_paths : a dictionary mapping assembly_id -> ftp_path

def fetch_refseq_genome(assembly_id, output_dir, ftp_paths):
    """Download .fna and .gff files for a given assembly from RefSeq."""
    if assembly_id not in ftp_paths:
        return None, None
    # Construct URL and path
    ftp_path = ftp_paths[assembly_id]
    basename = os.path.basename(ftp_path)
    fasta_url = f"{ftp_path}/{basename}_genomic.fna.gz"
    gff_url   = f"{ftp_path}/{basename}_genomic.gff.gz"
    # local save path
    os.makedirs(output_dir, exist_ok=True)
    fasta_path = os.path.join(output_dir, f"{assembly_id}.fna.gz")
    gff_path   = os.path.join(output_dir, f"{assembly_id}.gff.gz")
    # download files
    for url, file_path in [(fasta_url, fasta_path), (gff_url, gff_path)]:
        if not os.path.exists(file_path):
            for _ in range(MAX_RETRIES):
                try:
                    r = requests.get(url, headers={"User-Agent": USER_AGENT}, stream=True)
                    if r.status_code == 200:
                        with open(file_path, 'wb') as f:
                            f.write(r.content)
                        break
                    else:
                        time.sleep(2)
                except Exception:
                    time.sleep(2)
    return fasta_path, gff_path


def find_cds_by_protein_id(gff_path, protein_id):
    """Find a CDS feature in the GFF file by its protein_id."""
    # Opens the .gff.gz file in text mode.
    with gzip.open(gff_path, 'rt') as gff:
        for line in gff:
            # Skips comment lines
            if line.startswith("#"):
                continue
            # Splits the line into 9 tab-separated fields.
            parts = line.strip().split('	')
            if len(parts) != 9:
                continue
            seqid, source, feature_type, start, end, score, strand, phase, attributes = parts
            # Only CDS features are considered
            if feature_type != "CDS":
                continue
            if f"protein_id={protein_id}" in attributes:
                return {
                    "contig": seqid,
                    "start": int(start),
                    "end": int(end),
                    "strand": strand,
                    "attributes": attributes
                }
    return None


def find_all_integrases(gff_path, contig, strand, aice13_start, aice13_end, max_distance, strand_specific):
    """Find all nearby integrase CDS features near AICE13 hit in the GFF."""
    integrases = []
    # Defines the genomic search window
    aice13_window_start = aice13_start - max_distance
    aice13_window_end   = aice13_start + max_distance
    # Opens the .gff.gz file in text mode and skips all comment lines.
    with gzip.open(gff_path, 'rt') as gff:
        for line in gff:
            if line.startswith("#"):
                continue
            # Parses each GFF line and skips malformed ones.
            parts = line.strip().split('	')
            if len(parts) != 9:
                continue
            seqid, source, feature_type, start, end, score, fstrand, phase, attributes = parts
            # We only care about CDS features and only on the same contig as the AICE13 hit.
            if feature_type != "CDS" or seqid != contig:
                continue
            start, end = int(start), int(end)
            # If strand_specific is True, ignore CDSs that are not on the same strand as AICE13.
            if strand_specific and fstrand != strand:
                continue
            # Checks whether the CDS overlaps the AICE13 window in either start or end position.
            if (start >= aice13_window_start and start <= aice13_window_end) or \
               (end   >= aice13_window_start and end   <= aice13_window_end):
                # If "integrase" appears in the attribute field, it's recorded
                if "integrase" in attributes.lower():
                    integrases.append({
                        "start": start,
                        "end": end,
                        "strand": fstrand,
                        "attributes": attributes
                    })
    return integrases


def fetch_assembly_summary():
    """Download and parse RefSeq assembly summary file for FTP paths."""
    print("Downloading RefSeq assembly summary...")
    # Downloads a summary file which is a big tab-separated text file with a line for each RefSeq genome assembly.
    response = requests.get(REFSEQ_SUMMARY_URL, headers={"User-Agent": USER_AGENT})
    lines = response.text.splitlines()
    header_index = next(i for i, line in enumerate(lines) if 'assembly_accession' in line and 'ftp_path' in line)
    # Cleans the header line
    header_line = lines[header_index].lstrip("# ").strip()
    # Convert the data lines into a list of dictionaries, where each dictionary represents one genome assembly row
    data_lines = lines[header_index + 1:]
    reader = csv.DictReader([header_line] + data_lines, delimiter='	')
    # Returns a dictionary where key is assembly accession and value is ftp path to that assembly's files.
    return {row['assembly_accession']: row['ftp_path'] for row in reader}

# ------------------ MAIN FUNCTION ------------------

def main():
    parser = argparse.ArgumentParser(description="Analyze AICE13 hits with integrase detection")
    parser.add_argument("--input", required=True, help="Path to metadata.tsv")
    parser.add_argument("--output", default="results.tsv", help="Path to output TSV file")
    parser.add_argument("--stats", default="stats.tsv", help="Path to output stats file")
    parser.add_argument("--distance", type=int, default=DEFAULT_MAX_DISTANCE, help="Max distance to search for integrase")
    parser.add_argument("--workdir", default="genomes", help="Directory to download genomes")
    parser.add_argument("--strand-specific", action="store_true", help="Limit integrase search to same strand only")
    args = parser.parse_args()

    # loads the metadata
    entries = []
    with open(args.input, 'r') as f:
        reader = csv.DictReader(f, delimiter='	')
        for row in reader:
            row["coverage"] = float(row.get("Query coverage", "0").strip() or 0.0)
            entries.append(row)

    # download assembly summary
    ftp_paths = fetch_assembly_summary()

    stats = []
    found_metadata_integrase = 0
    found_manual_integrase = 0
    total_processed = 0
    no_integrase_found = 0

    # Iterate over each AICE13 hit in metadata
    for row in entries:
        query_type = row.get("Query", "").strip()
        # skip rows that are not AICE13 homologs or are integrase rows
        if not query_type.startswith("AICE13") or query_type == "AICE13_integrase":
            continue

        total_processed += 1
        aice13_hit = row.get("Hit", "").strip()
        target_proteome = row.get("Target Proteome Accession", "").strip()
        # Coordinates & strand from metadata for AICE13
        a_start  = int(row.get("Target start", 0))
        a_end    = int(row.get("Target stop", 0))
        a_strand = row.get("Target strand", "+")

        # Phase 1: metadata integrase pairing
        meta_integrases = [cand for cand in entries
                           if cand.get("Query", "").strip() == "AICE13_integrase"
                              and cand.get("Target Proteome Accession", "").strip() == target_proteome]
        if meta_integrases:
            # For each metadata integrase, compute distance using metadata coords
            for cand in meta_integrases:
                i_start  = int(cand.get("Target start", 0))
                i_end    = int(cand.get("Target stop", 0))
                i_strand = cand.get("Target strand", "+")
                # Compute distance
                if args.strand_specific and a_strand == i_strand:
                    if a_end < i_start:
                        dist = i_start - a_end - 1
                    elif i_end < a_start:
                        dist = a_start - i_end - 1
                    else:
                        dist = 0
                else:
                    a_mid = (a_start + a_end) // 2
                    i_mid = (i_start + i_end) // 2
                    dist = abs(a_mid - i_mid)
                stats.append({
                    "AICE13_Hit": aice13_hit,
                    "AICE13_Locus": f"{row['Target chr']}:{a_start}-{a_end}",
                    "AICE13_Strand": a_strand,
                    "Integrase_Hit": cand.get("Hit", ""),
                    "Integrase_Locus": f"{cand['Target chr']}:{i_start}-{i_end}",
                    "Integrase_Strand": i_strand,
                    "Source": "metadata",
                    "Distance": dist,
                    "Target_Proteome": target_proteome
                })
                found_metadata_integrase += 1
            continue  # skip manual search for this hit

        # Phase 2: manual GFF search for integrases
        # download genome & annotation
        assembly = target_proteome
        fasta_path, gff_path = fetch_refseq_genome(assembly, args.workdir, ftp_paths)
        if not (fasta_path and gff_path):
            continue
        # Locate AICE13 hit in GFF to double-check coords (optional)
        aice13_info = find_cds_by_protein_id(gff_path, aice13_hit)
        if not aice13_info:
            continue
        # Search for all integrases in window
        manual_integrases = find_all_integrases(
            gff_path,
            aice13_info['contig'],
            aice13_info['strand'],
            aice13_info['start'],
            aice13_info['end'],
            args.distance,
            args.strand_specific
        )
        if not manual_integrases:
            no_integrase_found += 1
            stats.append({
                "AICE13_Hit": aice13_hit,
                "AICE13_Locus": f"{aice13_info['contig']}:{aice13_info['start']}-{aice13_info['end']}",
                "AICE13_Strand": aice13_info['strand'],
                "Integrase_Hit": "NA",
                "Integrase_Locus": "NA",
                "Integrase_Strand": "NA",
                "Source": "none",
                "Distance": "NA",
                "Target_Proteome": target_proteome
            })
            continue
        # For each manual integrase, compute distance and append
        for integ in manual_integrases:
            # Extract integrase protein_id from GFF attributes
            attr = integ['attributes']
            prot_match = None
            for field in attr.split(';'):
                if field.startswith('protein_id='):
                    prot_match = field.split('=')[1]
                    break
            integrase_id = prot_match or 'NA'
            i_start = integ['start']
            i_end   = integ['end']
            i_strand = integ['strand']
            # Compute distance
            if args.strand_specific and aice13_info['strand'] == i_strand:
                if a_end < i_start:
                    dist = i_start - a_end - 1
                elif i_end < a_start:
                    dist = a_start - i_end - 1
                else:
                    dist = 0
            else:
                a_mid = (a_start + a_end) // 2
                i_mid = (i_start + i_end) // 2
                dist = abs(a_mid - i_mid)
            # Append with actual integrase ID
            stats.append({
                "AICE13_Hit": aice13_hit,
                "AICE13_Locus": f"{aice13_info['contig']}:{aice13_info['start']}-{aice13_info['end']}",
                "AICE13_Strand": aice13_info['strand'],
                "Integrase_Hit": integrase_id,
                "Integrase_Locus": f"{i_start}-{i_end}",
                "Integrase_Strand": i_strand,
                "Source": "manual",
                "Distance": dist,
                "Target_Proteome": target_proteome
            })
            found_manual_integrase += 1

    # Write stats TSV
    with open(args.stats, 'w') as out:
        writer = csv.DictWriter(
            out,
            fieldnames=["AICE13_Hit","AICE13_Locus","AICE13_Strand",
                       "Integrase_Hit","Integrase_Locus","Integrase_Strand",
                       "Source","Distance","Target_Proteome"],
            delimiter='\t'
        )
        writer.writeheader()
        writer.writerows(stats)

    # Summary log
    with open("stats.txt", "w") as f:
        f.write("Summary of Integrase Detection for AICE13 Homologs\n")
        f.write("==============================================\n")
        f.write(f" Total AICE13_homolog entries processed: {total_processed}\n")
        f.write(f" AICE13_homolog and AICE13_integrase_homolog pairs found in metadata: {found_metadata_integrase}\n")
        f.write(f" AICE13_homologs with integrase found manually from GFF: {found_manual_integrase}\n")
        f.write(f" AICE13_homologs with no integrase found at all: {no_integrase_found}\n")

if __name__ == "__main__":
    main()
