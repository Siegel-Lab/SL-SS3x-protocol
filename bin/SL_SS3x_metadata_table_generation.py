#!/usr/bin/env python3
import argparse
import os
import pandas as pd
from collections import defaultdict
from Bio import SeqIO
from Bio.Seq import Seq
import string
import sys

ROWS = list(string.ascii_uppercase[:16])  # A..P
COLS = list(range(1, 25))                  # 1..24


def revcomp(seq: str) -> str:
    return str(Seq(seq).reverse_complement())


def parse_well_spec(spec: str):
    """
    Parse a 'Well(s)' specification:
    - 'A3,C4,F12' -> ['A3','C4','F12']
    - 'A-P:1-12'   -> list of wells rows A..P and cols 1..12
    - 'Rest'       -> return 'Rest' sentinel
    """
    if pd.isna(spec):
        return []
    s = str(spec).strip()
    if s.lower() == "rest":
        return "Rest"

    wells = []
    # comma separated parts
    parts = [p.strip() for p in s.split(",") if p.strip()]
    for p in parts:
        if ":" in p:
            row_part, col_part = p.split(":", 1)
            if "-" not in row_part or "-" not in col_part:
                raise ValueError(f"Bad range spec: {p}")
            row_start, row_end = row_part.split("-", 1)
            row_start = row_start.strip().upper()
            row_end = row_end.strip().upper()
            col_start, col_end = map(int, col_part.split("-", 1))
            if row_start not in ROWS or row_end not in ROWS:
                raise ValueError(f"Row range out of A-P: {row_start}-{row_end}")
            if not (1 <= col_start <= 24 and 1 <= col_end <= 24):
                raise ValueError(f"Column range out of 1-24: {col_start}-{col_end}")
            for r_ord in range(ord(row_start), ord(row_end) + 1):
                r = chr(r_ord)
                for c in range(col_start, col_end + 1):
                    wells.append(f"{r}{c}")
        else:
            # single well token like A3
            token = p.upper()
            # very permissive parse: first letter(s) = row, trailing digits = col
            import re
            m = re.match(r"^([A-P])([1-9]|1[0-9]|2[0-4])$", token)
            if not m:
                raise ValueError(f"Invalid well token: '{p}'")
            wells.append(token)
    return wells


def load_index_records(index_name: str, folder: str):
    """
    Look for index_name.fasta in folder and return list(SeqRecord).
    """
    path = os.path.join(folder, f"{index_name}.fasta")
    if not os.path.exists(path):
        raise FileNotFoundError(f"Index FASTA not found: {path}")
    recs = list(SeqIO.parse(path, "fasta"))
    if not recs:
        raise ValueError(f"No records parsed from {path}")
    return recs


def main():
    p = argparse.ArgumentParser(description="Generate metadata table for multiplexed 384-well plates")
    p.add_argument("--config", required=True, help="Config file (Plate, i7_set, i5_set, tag) TSV")
    p.add_argument("--assignments", required=True, help="Assignments TSV (Plate, Well(s), Feature, Value)")
    p.add_argument("--tag_fasta", required=True, help="TAG FASTA (IDs like T000_TAG; sequences 11 nt)")
    p.add_argument("--index_sets_folder", default=".", help="Folder containing index set FASTA files (default current dir)")
    p.add_argument("--output", required=True, help="Output metadata TSV")
    p.add_argument("--tag_output", required=True, help="Output FASTA of used TAGs (with AGAGACAG prefix + revcomp)")
    args = p.parse_args()

    # --- read inputs ---
    config_df = pd.read_csv(args.config, sep="\t", dtype=str).fillna("")
    assignments_df = pd.read_csv(args.assignments, sep="\t", dtype=str).fillna("")
    tag_records = {rec.id: str(rec.seq) for rec in SeqIO.parse(args.tag_fasta, "fasta")}

    # Basic validation
    required_cfg_cols = {"Plate", "i7_set", "i5_set", "tag"}
    if not required_cfg_cols.issubset(set(config_df.columns)):
        raise ValueError(f"Config file must contain columns: {required_cfg_cols}")

    required_assign_cols = {"Plate", "Well(s)", "Feature", "Value"}
    if not required_assign_cols.issubset(set(assignments_df.columns)):
        raise ValueError(f"Assignments file must contain columns: {required_assign_cols}")

    # Map each plate -> well -> {feature: value}
    all_wells = [f"{r}{c}" for r in ROWS for c in COLS]
    plate_well_features = defaultdict(lambda: defaultdict(dict))
    # track which wells were explicitly assigned for each plate+feature (to handle Rest)
    assigned_by_plate_feature = defaultdict(lambda: defaultdict(set))
    rest_rows = []

    # First pass: explicit assignments (not Rest)
    for _, r in assignments_df.iterrows():
        plate = r["Plate"]
        spec = str(r["Well(s)"]).strip()
        feature = r["Feature"]
        value = r["Value"]

        parsed = parse_well_spec(spec)
        if parsed == "Rest":
            rest_rows.append((plate, feature, value))
            continue

        for well in parsed:
            if feature in plate_well_features[plate][well]:
                raise ValueError(f"Conflict: feature '{feature}' for {plate} {well} already assigned")
            plate_well_features[plate][well][feature] = value
            assigned_by_plate_feature[plate][feature].add(well)

    # Second pass: handle Rest per plate+feature
    for plate, feature, value in rest_rows:
        already = assigned_by_plate_feature[plate][feature]
        remaining = set(all_wells) - already
        for well in remaining:
            # do not overwrite an already existing feature (shouldn't happen given above)
            if feature in plate_well_features[plate][well]:
                continue
            plate_well_features[plate][well][feature] = value

    # --- Iterate plates from config (only plates listed in config) ---
    rows_out = []
    used_tags = set()
    for _, cfg in config_df.iterrows():
        plate = cfg["Plate"]
        i7_set_name = cfg["i7_set"]
        i5_set_name = cfg["i5_set"]
        tag_name = cfg["tag"]

        if tag_name not in tag_records:
            raise ValueError(f"TAG '{tag_name}' not found in tag_fasta")

        # load index sets for this plate
        i7_recs = load_index_records(i7_set_name, args.index_sets_folder)
        i5_recs = load_index_records(i5_set_name, args.index_sets_folder)

        if len(i7_recs) != 24:
            raise ValueError(f"i7 set '{i7_set_name}' must contain 24 records; found {len(i7_recs)}")
        if len(i5_recs) != 16:
            raise ValueError(f"i5 set '{i5_set_name}' must contain 16 records; found {len(i5_recs)}")

        # full tag sequence used in BC_TAG (prefix 'AGAGACAG' + 11-nt)
        tag11 = tag_records[tag_name]
        if len(tag11) != 11:
            raise ValueError(f"TAG '{tag_name}' sequence must be 11 nt (found {len(tag11)})")
#        tag_full = "AGAGACAG" + tag11
        tag_full_for_fasta = "AGAGACAG" + tag11  # only for FASTA output
        used_tags.add(tag_name)

        # For all wells (iterate in standard plate order A1..P24)
        for row_letter in ROWS:
            row_idx = ROWS.index(row_letter)  # 0-based -> i5 index
            for col_num in COLS:
                well = f"{row_letter}{col_num}"
                # pick sequences by column (i7) and row (i5)
                i7_rec = i7_recs[col_num - 1]   # column 1 -> index 0
                i5_rec = i5_recs[row_idx]       # row A -> index 0

                i7_seq = str(i7_rec.seq)
                i5_seq = str(i5_rec.seq)
                i7_rc = revcomp(i7_seq)
                i5_rc = revcomp(i5_seq)
                BC = i7_rc + i5_rc
                BC_TAG = BC + tag11

                # features (may be empty)
                feats = plate_well_features.get(plate, {}).get(well, {})
                # Build Library_ID by sorting feature names (stable)
                sorted_feat_names = sorted(feats.keys())
                library_id_parts = [feats[name] for name in sorted_feat_names]
                library_id = "_".join([p for p in library_id_parts if p])  # skip empty
                well_id = f"{plate}_{well}"
                if library_id:
                    well_id = f"{well_id}_{library_id}"

                rowdata = {
                    "Plate": plate,
                    "Well": well,
                    "i7": i7_seq,
                    "i7_rc": i7_rc,
                    "i5": i5_seq,
                    "i5_rc": i5_rc,
                    "BC": BC,
                    "TAG": tag11,
                    "BC_TAG": BC_TAG,
                    # include features as separate columns
                }
                # add feature columns (consistent set across plate later)
                for fn in sorted_feat_names:
                    rowdata[fn] = feats.get(fn, "")

                # identifier columns
                rowdata["Library_ID"] = library_id
                rowdata["Well_ID"] = well_id
                # add index record ids for traceability if desired
                rowdata["i7_name"] = i7_rec.id
                rowdata["i5_name"] = i5_rec.id

                rows_out.append(rowdata)

    # Build DataFrame
    df = pd.DataFrame(rows_out)

    # ensure Library_ID and Well_ID are last two columns
    cols = [c for c in df.columns if c not in ("Library_ID", "Well_ID")]
    cols = cols + ["Library_ID", "Well_ID"]
    df = df[cols]

    df.to_csv(args.output, sep="\t", index=False)

    # Write used TAGs FASTA (with AGAGACAG prefix) and their reverse complement
    with open(args.tag_output, "w") as outfa:
        for tag in sorted(used_tags):
            tag11 = tag_records[tag]
            seq_full = "AGAGACAG" + tag11  # prefix only here
#            seq_full = "AGAGACAG" + tag_records[tag]
            outfa.write(f">partialTn5_plus_{tag}\n{seq_full}\n")
            outfa.write(f">partialTn5_plus_{tag}_revcomp\n{revcomp(seq_full)}\n")

    print(f"[OK] wrote {args.output} and {args.tag_output}")


if __name__ == "__main__":
    main()

