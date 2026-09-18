#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Parse STARsolo chimeric junction output, match with circRNA BED and annotation file,
generate merged circ UMI annotation table (CLI version)
"""
import argparse
from collections import defaultdict


def star_parse_convert(chr1, site1, strand1, chr2, site2, strand2):
    """
    Convert STAR chimeric junction coordinates to back-spliced BED coordinates
    Return tuple (chrom, bed_start, bed_end) or None if invalid
    """
    if chr1 != chr2 or strand1 != strand2:
        return None
    s1 = int(site1)
    s2 = int(site2)
    if strand1 == "+":
        bed_start = s2
        bed_end = s1 - 1
    else:
        bed_start = s1
        bed_end = s2 - 1
    if bed_start > bed_end:
        return None
    return (chr1, bed_start, bed_end)


def merge_circ_annotation(name, match_window):
    chimeric = f'circexplorer2/{name}/01.chimeric/{name}_valid_Chimeric.out.junction'
    bed = f'circexplorer2/{name}/02.parse/{name}_back_spliced_junction.bed'
    anno = f'circexplorer2/{name}/03.annotate/{name}_circularRNA_known.txt'
    out = f'circexplorer2/{name}/04.matrix/{name}_circularRNA_umi.tsv'

    # Statistics counter
    stat = {
        "total_chimeric": 0,
        "valid_circ_coord": 0,
        "match_bed": 0,
        "anno_hit": 0,
        "anno_miss": 0,
        "final_unique_mol": 0
    }

    # 1. Read BED file: map (chr,start,end) to FUSIONJUNC_ID
    bed_key2fj = dict()
    with open(bed, "r") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            c, s, e, fjid, score, strand = line.split("\t")
            ck = (c, int(s), int(e))
            bed_key2fj[ck] = fjid

    # 2. Read annotation file, group records by chromosome
    anno_chr_dict = defaultdict(list)
    with open(anno, "r") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            c = parts[0]
            s_anno = int(parts[1])
            e_anno = int(parts[2])
            circ_id = parts[3]
            # parts[14] = gene_id(ENSMUSG), parts[15] = transcript_id
            gene = parts[14]
            tx = parts[15]
            exon_info = parts[-1]
            anno_chr_dict[c].append({
                "s": s_anno,
                "e": e_anno,
                "circRNA_id": circ_id,
                "gene_id": gene,
                "transcript_id": tx,
                "exon_info": exon_info
            })

    # 3. Annotation matching function: find nearest circRNA annotation within window
    def find_anno(chr_bed, s_bed, e_bed):
        if chr_bed not in anno_chr_dict:
            return None
        candidates = []
        for ann in anno_chr_dict[chr_bed]:
            if abs(ann["s"] - s_bed) <= match_window and abs(ann["e"] - e_bed) <= match_window:
                dist = abs(ann["s"] - s_bed) + abs(ann["e"] - e_bed)
                candidates.append((dist, ann))
        if not candidates:
            return None
        candidates.sort(key=lambda x: x[0])
        return candidates[0][1]

    # 4. Read chimeric junction, convert coordinates, match annotation, deduplicate by FUSIONJUNC+CB+UMI
    circ_umi_set = set()
    out_records = []
    with open(chimeric, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            stat["total_chimeric"] += 1
            parts = line.split("\t")
            flag = int(parts[6])
            if flag < 0:
                continue
            chr1, site1, strand1, chr2, site2, strand2 = parts[:6]
            cb = parts[-2]
            umi = parts[-1]
            circ_key = star_parse_convert(chr1, site1, strand1, chr2, site2, strand2)
            if circ_key is None:
                continue
            stat["valid_circ_coord"] += 1
            c_bed, s_bed, e_bed = circ_key
            if circ_key not in bed_key2fj:
                continue
            stat["match_bed"] += 1
            fj_id = bed_key2fj[circ_key]
            ann_info = find_anno(c_bed, s_bed, e_bed)
            if ann_info is None:
                stat["anno_miss"] += 1
                continue
            stat["anno_hit"] += 1
            # Filter empty gene_id
            gene_id = ann_info["gene_id"].strip()
            if gene_id == "" or gene_id == "None":
                continue
            # Deduplication key: FUSIONJUNC_ID + CB + UMI
            unique_key = (fj_id, cb, umi)
            if unique_key not in circ_umi_set:
                circ_umi_set.add(unique_key)
                out_records.append({
                    "chr": c_bed,
                    "bed_start": s_bed,
                    "bed_end": e_bed,
                    "FUSIONJUNC_ID": fj_id,
                    "circRNA_id": ann_info["circRNA_id"],
                    "CB": cb,
                    "UMI": umi,
                    "gene_id": gene_id,
                    "transcript_id": ann_info["transcript_id"],
                    "exon_info": ann_info["exon_info"]
                })
    stat["final_unique_mol"] = len(circ_umi_set)

    # 5. Write output merged TSV
    with open(out, "w") as fout:
        header = "chr\tbed_start\tbed_end\tFUSIONJUNC_ID\tcircRNA_id\tCB\tUMI\tgene_id\ttranscript_id\texon_info\n"
        fout.write(header)
        for rec in out_records:
            row = (f"{rec['chr']}\t{rec['bed_start']}\t{rec['bed_end']}\t{rec['FUSIONJUNC_ID']}\t{rec['circRNA_id']}\t"
                   f"{rec['CB']}\t{rec['UMI']}\t{rec['gene_id']}\t{rec['transcript_id']}\t{rec['exon_info']}\n")
            fout.write(row)

    print(f"Output file: {out}")
    print("===== Summary Statistics =====")
    for k, v in stat.items():
        print(f"{k}: {v}")
    print(f"Retained unique molecular records: {len(circ_umi_set)}")



def parse_args():
    parser = argparse.ArgumentParser(
        description="Merge STARsolo chimeric junctions with circRNA BED and annotation to produce circ_umi_anno_merged.tsv"
    )
    parser.add_argument("--name", required=True, help="name")
    parser.add_argument("--match_window", type=int, default=10, help="Tolerance window (bp) for coordinate matching, default=10")
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    merge_circ_annotation(args.name, args.match_window)
