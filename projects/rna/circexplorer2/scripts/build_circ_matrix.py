#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Build circRNA locus (gene + back-splice site) x cell UMI count matrix
Row: gene_id|chr:start-end  (e.g. ENSMUSG00000064337|MT:69-69)
Column: all CeleScope filtered barcodes

Additional output: a barcode-filtered copy of the raw UMI table
  ({name}_circularRNA_umi_barcode_filtered.tsv) - keeps every original
  column/row whose CB is in the CeleScope filtered barcode set.
"""
import argparse
import gzip
import statistics
from collections import defaultdict

TOP_N = 10


def is_valid_gene_id(gene_id):
    """Validate gene_id format: retain only valid Ensembl gene IDs starting with ENSMUSG"""
    if not gene_id or gene_id in ("NA", "None", ""):
        return False
    gene_clean = gene_id.strip()
    return gene_clean.startswith("ENSMUSG")


def make_locus_key(chr_, start, end, gene_id):
    """Generate unique locus key: gene_id|chr:start-end"""
    return f"{gene_id}|{chr_}:{start}-{end}"


def build_circ_matrix(name, barcode):
    input = f'circexplorer2/{name}/04.matrix/{name}_circularRNA_umi.tsv'
    out_mtx = f'circexplorer2/{name}/04.matrix/{name}_circularRNA_mtx.tsv'
    out_stat = f'circexplorer2/{name}/04.matrix/{name}_circularRNA_stats.txt'
    # NEW: barcode-filtered copy of the raw UMI table
    out_filtered = f'circexplorer2/{name}/04.matrix/{name}_circularRNA_umi_barcode_filtered.tsv'

    # Column index mapping
    COL = {
        "chr": 0, "bed_start": 1, "bed_end": 2,
        "FUSIONJUNC_ID": 3, "circRNA_id": 4,
        "CB": 5, "UMI": 6,
        "gene_id": 7, "transcript_id": 8, "exon_info": 9,
    }

    # 1. Load all CeleScope filtered barcodes
    all_cbs = []
    valid_cb_set = set()
    with gzip.open(barcode, "rt") as f:
        for line in f:
            b = line.strip()
            if b:
                all_cbs.append(b)
                valid_cb_set.add(b)
    print(f"CeleScope filtered barcode: {len(all_cbs)}")

    # cell_locus_umi[CB][locus_key] = set(UMIs)
    cell_locus_umi = defaultdict(lambda: defaultdict(set))
    positive_cbs = set()
    discarded_dirty_cb = 0
    discarded_invalid_gene = 0
    total_records = 0
    written_filtered_records = 0

    # 2. Read merged circ table
    with open(input, "r") as f, open(out_filtered, "w") as fout_filtered:
        header = f.readline()
        # write header to the barcode-filtered copy
        fout_filtered.write(header)

        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            total_records += 1
            parts = line.split("\t")
            if len(parts) < 10:
                continue

            chr_ = parts[COL["chr"]]
            s_bed = parts[COL["bed_start"]]
            e_bed = parts[COL["bed_end"]]
            cb   = parts[COL["CB"]]
            umi  = parts[COL["UMI"]]
            gene = parts[COL["gene_id"]]

            # --- Barcode-based filter: decide first, before gene/UMI checks ---
            # Rows without a usable CB cannot be classified -> skip (also not written)
            if not cb or cb == "NA":
                continue
            if cb not in valid_cb_set:
                discarded_dirty_cb += 1
                continue
            # CB passes the filtered set -> keep the original row as-is in the filtered copy
            fout_filtered.write(line + "\n")
            written_filtered_records += 1

            # --- Downstream matrix field validation ---
            if not umi or umi == "NA" or not gene or gene == "NA":
                continue
            if not is_valid_gene_id(gene):
                discarded_invalid_gene += 1
                continue

            # Build locus key: gene_id|chr:start-end
            locus_key = make_locus_key(chr_, s_bed, e_bed, gene)
            cell_locus_umi[cb][locus_key].add(umi)
            positive_cbs.add(cb)

    # Collect all unique loci
    all_loci = set()
    for cb, ldict in cell_locus_umi.items():
        all_loci.update(ldict.keys())
    all_loci = sorted(all_loci)

    # ========== Locus-level statistics ==========
    locus_stats = {}
    for loc in all_loci:
        total_umi = 0
        cells_detected = 0
        for cb in positive_cbs:
            umi_set = cell_locus_umi[cb].get(loc, set())
            if len(umi_set) > 0:
                cells_detected += 1
                total_umi += len(umi_set)
        locus_stats[loc] = {
            "total_umi": total_umi,
            "cells_detected": cells_detected,
            "detection_rate": cells_detected / len(positive_cbs) if positive_cbs else 0
        }

    sorted_loci = sorted(all_loci, key=lambda l: locus_stats[l]["total_umi"], reverse=True)
    top_loci = sorted_loci[:TOP_N]

    # ========== Write expression matrix ==========
    with open(out_mtx, "w") as fout:
        header_line = "locus\t" + "\t".join(all_cbs) + "\n"
        fout.write(header_line)
        for loc in sorted_loci:
            row_parts = [loc]
            for cb in all_cbs:
                umi_set = cell_locus_umi[cb].get(loc, set())
                row_parts.append(str(len(umi_set)))
            fout.write("\t".join(row_parts) + "\n")

    # ========== Calculate summary metrics ==========
    n_total_cells = len(all_cbs)
    n_cells_with_circ = len(positive_cbs)
    n_no_circ = n_total_cells - n_cells_with_circ

    loci_per_cell = []
    umis_per_cell  = []
    for cb in positive_cbs:
        ldict = cell_locus_umi[cb]
        loci_per_cell.append(len(ldict))
        umis_per_cell.append(sum(len(umis) for umis in ldict.values()))

    median_loci = statistics.median(loci_per_cell) if loci_per_cell else 0
    median_umi  = statistics.median(umis_per_cell)  if umis_per_cell else 0
    mean_loci   = statistics.mean(loci_per_cell) if loci_per_cell else 0
    mean_umi    = statistics.mean(umis_per_cell) if umis_per_cell else 0
    max_umi     = max(umis_per_cell) if umis_per_cell else 0
    max_loci    = max(loci_per_cell) if loci_per_cell else 0

    total_umi_all = sum(locus_stats[l]["total_umi"] for l in all_loci)

    # ========== Output statistics summary ==========
    summary_lines = [
        f"Input circ file: {input}",
        f"Barcode source: {barcode}",
        f"Barcode-filtered UMI table: {out_filtered}",
        f"Records written to barcode-filtered table: {written_filtered_records}",
        f"Total records (merged table): {total_records}",
        f"Matrix dimension: {len(all_loci)} circRNA loci x {n_total_cells} CeleScope filtered cells",
        f"Matrix file: {out_mtx}",
        f"Row format: gene_id|chr:start-end",
        f"----------------------------------------",
        f"CeleScope filtered total cells: {n_total_cells}",
        f"Cells detected with circRNA: {n_cells_with_circ}",
        f"Cells without circ signal: {n_no_circ}",
        f"Cell detection rate: {n_cells_with_circ/n_total_cells*100:.2f}%",
        f"----------------------------------------",
        f"Discarded records with barcodes not in filtered set: {discarded_dirty_cb}",
        f"Discarded records with invalid gene_id (non-ENSMUSG format): {discarded_invalid_gene}",
        f"----------------------------------------",
        f"[Positive cell statistics]",
        f"  Detected circRNA loci per cell: median={median_loci:.1f}, mean={mean_loci:.2f}, max={max_loci}",
        f"  Total circ UMI per cell: median={median_umi:.1f}, mean={mean_umi:.2f}, max={max_umi}",
        f"----------------------------------------",
        f"[Locus-level statistics] TOP{TOP_N} (sorted by total UMI, descending)",
        f"  Total unique circRNA loci: {len(all_loci)}",
        f"  Total circ UMI across all loci: {total_umi_all}",
        f"",
        f"  {'locus':<45} {'Total_UMI':>10} {'Detected_Cells':>12} {'Detection_Rate(%)':>12}",
        f"  {'-'*45} {'-'*10} {'-'*12} {'-'*12}",
    ]
    for loc in top_loci:
        s = locus_stats[loc]
        summary_lines.append(
            f"  {loc:<45} {s['total_umi']:>10} {s['cells_detected']:>12} {s['detection_rate']*100:>11.2f}%"
        )
    echo_summary = "\n".join(summary_lines) + "\n"
    print(echo_summary)

    for loc in sorted_loci:
        s = locus_stats[loc]
        summary_lines.append(
            f"  {loc:<45} {s['total_umi']:>10} {s['cells_detected']:>12} {s['detection_rate']*100:>11.2f}%"
        )
    summary = "\n".join(summary_lines) + "\n"
    with open(out_stat, "w") as f:
        f.write(summary)

    print(f"Matrix output completed: {out_mtx}")
    print(f"Barcode-filtered UMI table output completed: {out_filtered}")
    print(f"Statistics summary output completed (TOP{TOP_N}): {out_stat}")


def parse_args():
    parser = argparse.ArgumentParser(description="Build circRNA locus x cell UMI expression matrix")
    parser.add_argument("--name", required=True, help="name")
    parser.add_argument("--barcode", required=True, help="Path to CeleScope filtered barcodes.tsv.gz")
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    build_circ_matrix(args.name, args.barcode)
