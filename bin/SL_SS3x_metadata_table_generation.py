
import pandas as pd
import argparse
from pathlib import Path
from collections import defaultdict
import re

def parse_range(spec):
    if spec == "Rest":
        return "Rest"
    wells = []
    for part in spec.split(','):
        if ':' in part:
            row_range, col_range = part.split(':')
            row_start, row_end = row_range.split('-')
            col_start, col_end = map(int, col_range.split('-'))
            for row in range(ord(row_start), ord(row_end)+1):
                for col in range(col_start, col_end+1):
                    wells.append(f"{chr(row)}{col}")
        else:
            wells.append(part)
    return wells

def reverse_complement(seq):
    complement = str.maketrans("ATCGN", "TAGCN")
    return seq.translate(complement)[::-1]

def generate_all_wells(rows="A-P", cols=range(1, 25)):
    all_wells = []
    for row in range(ord(rows[0]), ord(rows[-1])+1):
        for col in cols:
            all_wells.append(f"{chr(row)}{col}")
    return all_wells

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--assignments", required=True)
    parser.add_argument("--config", required=True)
    parser.add_argument("--i7", required=True)
    parser.add_argument("--i5", required=True)
    parser.add_argument("--tags", required=True)
    parser.add_argument("--output_table", required=True)
    parser.add_argument("--used_tags_fasta", required=True)
    args = parser.parse_args()

    assignments_df = pd.read_csv(args.assignments, sep='\t')
    config_df = pd.read_csv(args.config, sep='\t')
    i7_df = pd.read_csv(args.i7, sep='\t', index_col=0)
    i5_df = pd.read_csv(args.i5, sep='\t', index_col=0)
    tag_df = pd.read_csv(args.tags, sep='\t', index_col=0)

    all_possible_wells = generate_all_wells()

    # Step 1: Collect feature assignments per plate
    feature_assignments = defaultdict(lambda: defaultdict(dict))
    assigned_features = defaultdict(lambda: defaultdict(set))

    for _, row in assignments_df.iterrows():
        plate = row["Plate"]
        wells_spec = row["Well(s)"]
        feature = row["Feature"]
        value = row["Value"]

        if wells_spec == "Rest":
            continue  # Process after
        wells = parse_range(wells_spec)
        for well in wells:
            feature_assignments[plate][well][feature] = value
            assigned_features[plate][feature].add(well)

    # Step 2: Handle "Rest"
    for _, row in assignments_df.iterrows():
        plate = row["Plate"]
        wells_spec = row["Well(s)"]
        feature = row["Feature"]
        value = row["Value"]

        if wells_spec != "Rest":
            continue

        # Find unassigned wells for this feature
        assigned = assigned_features[plate][feature]
        unassigned = set(all_possible_wells) - assigned
        for well in unassigned:
            feature_assignments[plate][well][feature] = value

    # Step 3: Assign barcodes
    results = []
    used_tags = set()

    for _, config in config_df.iterrows():
        plate = config["Plate"]
        i7_set = config["i7_set"]
        i5_set = config["i5_set"]
        tag_id = config["tag"]

        tag_seq = tag_df.loc[tag_id].values[0]
        tag_full = f"AGAGACAG{tag_seq}"
        tag_rc = reverse_complement(tag_full)
        used_tags.add((f"partialTn5_plus_{tag_id}", tag_full))
        used_tags.add((f"partialTn5_plus_{tag_id}_revcomp", tag_rc))

        i7_indexes = i7_df.loc[i7_set]
        i5_indexes = i5_df.loc[i5_set]

        for well, i7 in zip(all_possible_wells, i7_indexes.index):
            i5 = i5_indexes.index[list(all_possible_wells).index(well) % len(i5_indexes)]
            i7_seq = i7_indexes.loc[i7]
            i5_seq = i5_indexes.loc[i5]
            i7_rc = reverse_complement(i7_seq)
            i5_rc = reverse_complement(i5_seq)
            BC = i7_rc + i5_seq
            BC_TAG = BC + tag_full

            well_features = feature_assignments[plate].get(well, {})
            library_id_parts = [plate, well] + [well_features.get(f, "") for f in sorted(well_features)]
            library_id = "_".join(filter(None, library_id_parts))
            row_data = {
                "Plate": plate,
                "Well": well,
                "i7": i7_seq,
                "i7_rc": i7_rc,
                "i5": i5_seq,
                "i5_rc": i5_rc,
                "TAG": tag_seq,
                "BC": BC,
                "BC_TAG": BC_TAG,
                "Well_ID": f"{plate}_{well}",
                "Library_ID": library_id
            }
            row_data.update(well_features)
            results.append(row_data)

    df = pd.DataFrame(results)

    # Move Well_ID and Library_ID to end
    cols = [c for c in df.columns if c not in ["Well_ID", "Library_ID"]] + ["Well_ID", "Library_ID"]
    df = df[cols]
    df.to_csv(args.output_table, sep='\t', index=False)

    # Output used TAGs fasta
    with open(args.used_tags_fasta, 'w') as f:
        for name, seq in sorted(used_tags):
            f.write(f">{name}\n{seq}\n")

if __name__ == "__main__":
    main()
