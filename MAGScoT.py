#!/usr/bin/env python3
"""
MAGScoT - Bin Scoring and Refinement Tool for Metagenome Assemblies.
Python port of https://github.com/ikmb/MAGScoT/blob/main/MAGScoT.R
"""

import argparse
import hashlib
import os
import re
import sys
from collections import defaultdict

import pandas as pd


def parse_args():
    p = argparse.ArgumentParser(
        description="MAGScoT: Bin Scoring and Refinement Tool"
    )
    p.add_argument("-i", "--input", required=True,
                   help="Tab-separated input: bin, contig, set (no header)")
    p.add_argument("--hmm", required=True,
                   help="Tab-separated HMM markers: gene_id, marker, e-value (no header)")
    p.add_argument("-p", "--profile", default="default",
                   choices=["default", "bac120", "ar53"],
                   help="Marker profile from GTDB rel207 [default]")
    p.add_argument("--profile_dir", default=None,
                   help="Directory containing marker profile TSVs "
                        "(default: profiles/ next to this script)")
    p.add_argument("-o", "--out", default="MAGScoT",
                   help="Output file name base [MAGScoT]")
    p.add_argument("--bin_separator", default="cleanbin",
                   help="Separator for refined output bin names [cleanbin]")
    p.add_argument("-a", "--score_a", type=float, default=1.0,
                   help="Scoring parameter a [1.0]")
    p.add_argument("-b", "--score_b", type=float, default=0.5,
                   help="Scoring parameter b [0.5]")
    p.add_argument("-c", "--score_c", type=float, default=0.5,
                   help="Scoring parameter c [0.5]")
    p.add_argument("--max_cont", type=float, default=1.0,
                   help="Maximum contamination threshold [1.0]")
    p.add_argument("-t", "--threshold", type=float, default=0.5,
                   help="Minimum completeness threshold [0.5]")
    p.add_argument("--score_only", action="store_true",
                   help="Only do scoring, no refinement")
    p.add_argument("--skip_merge_bins", action="store_true",
                   help="Skip bin merging")
    p.add_argument("-m", "--min_markers", type=int, default=25,
                   help="Min unique markers for merge seed bins [25]")
    p.add_argument("-s", "--min_sharing", type=float, default=0.8,
                   help="Min fraction of shared markers for merging [0.8]")
    p.add_argument("-n", "--n_iterations", type=int, default=2,
                   help="Number of merging iterations [2]")
    return p.parse_args()


def md5_hash(s):
    return hashlib.md5(s.encode()).hexdigest()


def compute_binhashes(ctb, contig_idx_map):
    """Compute MD5 hash per bin based on sorted contig indices."""
    hashes = {}
    for bin_name, grp in ctb.groupby("bin"):
        indices = sorted(contig_idx_map[c] for c in grp["contig"] if c in contig_idx_map)
        hashes[bin_name] = md5_hash("".join(str(i) for i in indices))
    return hashes


def score_bins(contig_to_bin_remain, gene_to_contig, markers, zz_prev,
               recalc_bins, mag_exclude, a, b, c):
    """Score all remaining bins. Returns (zz, df_list) DataFrames."""
    # Join remaining contigs to their marker genes
    ctb_with_genes = contig_to_bin_remain[
        contig_to_bin_remain["contig"].isin(gene_to_contig["contig"])
    ].merge(gene_to_contig, on="contig")

    if recalc_bins is not None and zz_prev is not None:
        # Keep previously computed stats for unchanged bins
        keep = zz_prev[
            zz_prev["bin"].isin(contig_to_bin_remain["bin"].unique()) &
            ~zz_prev["bin"].isin(recalc_bins) &
            ~zz_prev["bin"].isin(mag_exclude)
        ]
        # Recompute for changed bins
        new_data = ctb_with_genes[
            ctb_with_genes["bin"].isin(recalc_bins) &
            ~ctb_with_genes["bin"].isin(mag_exclude)
        ]
        if len(new_data) > 0:
            new_zz = (new_data[["bin", "gene"]]
                      .groupby(["bin", "gene"]).size()
                      .reset_index(name="count"))
            new_zz = new_zz.merge(markers, left_on="gene", right_on="marker",
                                  how="left").drop(columns=["marker"])
            zz = pd.concat([keep, new_zz], ignore_index=True)
        else:
            zz = keep
    else:
        # First iteration: compute all
        zz = (ctb_with_genes[["bin", "gene"]]
              .groupby(["bin", "gene"]).size()
              .reset_index(name="count"))
        zz = zz.merge(markers, left_on="gene", right_on="marker",
                       how="left").drop(columns=["marker"])

    if len(zz) == 0:
        return zz, pd.DataFrame()

    # Group by (bin, set) to get completeness/contamination per marker set
    agg = zz.groupby(["bin", "set"]).agg(
        uniqueSCGs=("gene", "size"),
        multipleSCGs=("count", lambda x: (x > 1).sum()),
        sumSCGs=("count", "sum")
    ).reset_index()

    # Number of markers per set
    n_markers = markers.groupby("set").size().reset_index(name="n_markers")
    agg = agg.merge(n_markers, on="set", how="left")

    # Number of contigs per bin
    n_contigs = (contig_to_bin_remain.groupby("bin").size()
                 .reset_index(name="n_contigs"))
    agg = agg.merge(n_contigs, on="bin", how="left")

    # Compute scores
    agg["Completeness"] = agg["uniqueSCGs"] / agg["n_markers"]
    agg["additionalSCGs"] = agg["sumSCGs"] - agg["uniqueSCGs"] - agg["multipleSCGs"]
    agg["Contamination"] = agg["multipleSCGs"] / agg["uniqueSCGs"]
    agg["MultiContam"] = agg["additionalSCGs"] / agg["n_markers"]
    agg["Contamination_Score"] = b * agg["Contamination"] + c * agg["MultiContam"]
    agg["score"] = a * agg["Completeness"] - agg["Contamination_Score"]

    return zz, agg


def main():
    args = parse_args()

    # Resolve profile directory
    if args.profile_dir:
        profile_dir = args.profile_dir
    else:
        profile_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                   "profiles")

    profile_map = {
        "default": "gtdb_rel207_default_markers.tsv",
        "bac120": "gtdb_rel207_bac120_markers.tsv",
        "ar53": "gtdb_rel207_ar53_markers.tsv",
    }
    profile_path = os.path.join(profile_dir, profile_map[args.profile])
    if not os.path.exists(profile_path):
        print(f"ERROR: Profile file not found: {profile_path}", file=sys.stderr)
        sys.exit(1)

    # ---- Read inputs ----
    print("Loading data...", file=sys.stderr)

    contig_to_bin = pd.read_csv(args.input, sep="\t", header=None,
                                names=["bin", "contig", "set"])
    contig_to_bin = contig_to_bin.dropna(subset=["set"])

    if contig_to_bin.shape[1] != 3:
        print(f"ERROR: Input file should have 3 tab-separated columns, "
              f"got {contig_to_bin.shape[1]}", file=sys.stderr)
        sys.exit(1)

    n_bins = contig_to_bin["bin"].nunique()
    n_contigs = contig_to_bin["contig"].nunique()
    if n_bins > n_contigs:
        print("ERROR: More unique bins than contigs. "
              "Expected format: bin, contig, set (no header)", file=sys.stderr)
        sys.exit(1)

    # Build contig index map (sorted unique contigs → integer index)
    all_contigs_sorted = sorted(contig_to_bin["contig"].unique())
    contig_idx_map = {c: i + 1 for i, c in enumerate(all_contigs_sorted)}

    num_binsets = contig_to_bin["set"].nunique()
    if num_binsets == 1:
        print("Only one binning algorithm provided. Scoring bins only.",
              file=sys.stderr)
        args.score_only = True
        args.n_iterations = 0

    if args.n_iterations >= num_binsets and not args.skip_merge_bins:
        print(f"Merging iterations ({args.n_iterations}) >= number of binners "
              f"({num_binsets}). Setting to {num_binsets - 1}.", file=sys.stderr)
        args.n_iterations = num_binsets - 1

    # Read marker profiles
    markers = pd.read_csv(profile_path, sep="\t")
    markers["marker"] = markers["marker"].str.replace(r"\.HMM$|\.hmm$", "",
                                                       regex=True)

    # Read HMM hits: gene_id, marker, e-value
    hmm = pd.read_csv(args.hmm, sep="\t", header=None,
                       names=["gene_id", "marker", "evalue"])
    if hmm.shape[1] != 3:
        print(f"ERROR: HMM file should have 3 columns, got {hmm.shape[1]}",
              file=sys.stderr)
        sys.exit(1)

    # Best hit per gene (lowest e-value), strip trailing _N for contig name
    hmm = hmm.sort_values("evalue").drop_duplicates("gene_id", keep="first")
    hmm["contig"] = hmm["gene_id"].str.replace(r"_\d+$", "", regex=True)
    hmm["gene"] = hmm["marker"]
    gene_to_contig = hmm[["contig", "gene"]].drop_duplicates()

    # Compute initial bin hashes
    binhashes = compute_binhashes(contig_to_bin, contig_idx_map)

    # ---- Bin merging ----
    if args.skip_merge_bins or args.score_only or args.n_iterations == 0:
        print("Skipping bin merging...", file=sys.stderr)
    else:
        contig_to_bin_bisco = contig_to_bin.copy()
        contig_to_bin_bisco_it = contig_to_bin.copy()
        bisco_out_all = []

        for it in range(1, args.n_iterations + 1):
            print(f"Bin-Merging iteration: {it}", file=sys.stderr)
            if it == 1:
                contig_to_bin_bisco_it = contig_to_bin.copy()
                contig_to_bin_bisco = contig_to_bin.copy()

            # Join contigs with genes for current iteration bins
            valid_bins = set(binhashes.keys())
            ctb_genes = contig_to_bin_bisco_it[
                contig_to_bin_bisco_it["contig"].isin(gene_to_contig["contig"]) &
                contig_to_bin_bisco_it["bin"].isin(valid_bins)
            ].merge(gene_to_contig, on="contig")

            # Find merge candidates: bins with min_markers <= unique markers <= 150
            gene_counts = (ctb_genes[["bin", "gene"]]
                          .drop_duplicates()
                          .groupby("bin").size()
                          .reset_index(name="n_unique"))
            total_counts = (ctb_genes.groupby("bin")["gene"].size()
                           .reset_index(name="n_total"))
            cand_df = gene_counts.merge(total_counts, on="bin")
            merge_cand = cand_df[
                (cand_df["n_unique"] >= args.min_markers) &
                (cand_df["n_total"] <= 150)
            ]["bin"].tolist()

            if not merge_cand:
                break

            # Compare marker-carrying contigs across all candidate bins
            bisco_results = []
            prev_bisco_out = (pd.DataFrame(bisco_out_all)
                              if bisco_out_all else pd.DataFrame())

            for idx, thisbin in enumerate(merge_cand):
                print(f"\rProcessing candidate bin {idx + 1} of "
                      f"{len(merge_cand)}", end="", file=sys.stderr)

                thisbin_contigs = set(
                    contig_to_bin_bisco_it[
                        (contig_to_bin_bisco_it["bin"] == thisbin) &
                        contig_to_bin_bisco_it["contig"].isin(
                            gene_to_contig["contig"])
                    ]["contig"]
                )
                thisbin_genes = gene_to_contig[
                    gene_to_contig["contig"].isin(thisbin_contigs)
                ]
                n_unique_genes = thisbin_genes["gene"].nunique()
                if n_unique_genes == 0:
                    continue

                # Blacklist: bins already merged with thisbin in previous iters
                blacklist = set()
                if it > 1 and len(prev_bisco_out) > 0:
                    mask = prev_bisco_out["bin"] == thisbin
                    if mask.any():
                        bl = prev_bisco_out.loc[mask, ["bin_a", "bin_b"]]
                        blacklist = set(bl.values.flatten())

                # Find other bins sharing contigs with thisbin's marker genes
                shared = contig_to_bin[
                    contig_to_bin["contig"].isin(thisbin_genes["contig"])
                ].merge(thisbin_genes, on="contig")

                # Count shared markers per bin (distinct gene per bin)
                shared_counts = (shared[["bin", "gene"]]
                                .drop_duplicates()
                                .groupby("bin").size()
                                .reset_index(name="shared_markers"))
                shared_counts["shared_markers_rel"] = (
                    shared_counts["shared_markers"] / n_unique_genes
                )

                # Filter: different bin, not blacklisted, sufficient sharing
                hits = shared_counts[
                    (shared_counts["bin"] != thisbin) &
                    (~shared_counts["bin"].isin(blacklist)) &
                    (shared_counts["shared_markers_rel"] >= args.min_sharing)
                ].copy()

                for _, row in hits.iterrows():
                    other = row["bin"]
                    bin_a = max(thisbin, other)
                    bin_b = min(thisbin, other)
                    bisco_results.append({
                        "seed": thisbin, "bin_a": bin_a, "bin_b": bin_b,
                        "shared_markers": int(row["shared_markers"]),
                        "shared_markers_rel": row["shared_markers_rel"],
                    })

            print("", file=sys.stderr)

            if not bisco_results:
                print("Exiting bin merging, no (more) overlaps found",
                      file=sys.stderr)
                break

            bisco_out = pd.DataFrame(bisco_results)
            bisco_out = bisco_out.drop_duplicates(subset=["bin_a", "bin_b"],
                                                  keep="first")
            bisco_out = bisco_out[
                bisco_out["shared_markers"] >= args.min_markers
            ].reset_index(drop=True)

            if len(bisco_out) == 0:
                break

            bisco_out["bin"] = [
                f"MAGScoT_{it}_{j + 1}" for j in range(len(bisco_out))
            ]

            # Build contig-to-bin for merged bins
            merged_rows = []
            for _, row in bisco_out.iterrows():
                pair_contigs = contig_to_bin[
                    contig_to_bin["bin"].isin([row["bin_a"], row["bin_b"]])
                ][["contig"]].drop_duplicates()
                for contig in pair_contigs["contig"]:
                    merged_rows.append({
                        "bin": row["bin"],
                        "contig": contig,
                        "set": f"MAGScoT_{it}",
                    })

            contig_to_bin_bisco_it = pd.DataFrame(merged_rows)

            # Hash merged bins; only keep those distinct from existing
            merged_hashes = compute_binhashes(contig_to_bin_bisco_it,
                                              contig_idx_map)
            existing_hash_set = set(binhashes.values())
            new_bins = [b for b, h in merged_hashes.items()
                        if h not in existing_hash_set]

            contig_to_bin_bisco = pd.concat([
                contig_to_bin_bisco,
                contig_to_bin_bisco_it[
                    contig_to_bin_bisco_it["bin"].isin(new_bins)
                ]
            ], ignore_index=True)

            # Recompute all hashes
            binhashes = compute_binhashes(contig_to_bin_bisco, contig_idx_map)
            bisco_out_all.extend(bisco_out.to_dict("records"))

            if len(contig_to_bin_bisco_it) == 0:
                print("Exiting bin merging, no (more) overlaps found",
                      file=sys.stderr)
                break

        contig_to_bin = contig_to_bin_bisco

    # ---- Scoring and refinement ----
    contig_to_bin_remain = contig_to_bin.copy()
    contig_to_bin_out = pd.DataFrame(columns=contig_to_bin.columns)
    mag_exclude = set()
    scoreframe_out = []

    # Get binner tool per bin (from original contig_to_bin)
    bin_tool_map = (contig_to_bin[["bin", "set"]]
                    .drop_duplicates("bin")
                    .set_index("bin")["set"]
                    .to_dict())

    zz = None
    recalc = set(contig_to_bin_remain["bin"].unique())
    iteration = 1

    while len(contig_to_bin_remain) > 1:
        print(f"Refining bins, iteration: {iteration}", file=sys.stderr)
        print("Extracting SCG information for bins...", file=sys.stderr)

        zz, df_list = score_bins(
            contig_to_bin_remain, gene_to_contig, markers, zz,
            recalc if recalc else None, mag_exclude,
            args.score_a, args.score_b, args.score_c
        )

        if len(df_list) == 0:
            break

        # Build scoreframe: best score per bin, sorted descending
        # Prefer non-MAGScoT bins at ties (grepl("MAGScoT",bin) → sort True last)
        df_list["is_magscot"] = df_list["bin"].str.startswith("MAGScoT")
        scoreframe = (df_list
                      .sort_values(["score", "n_contigs", "is_magscot"],
                                   ascending=[False, True, True])
                      .drop_duplicates("bin", keep="first")
                      .reset_index(drop=True))

        # Add bintool column
        scoreframe["bintool"] = scoreframe["bin"].map(bin_tool_map)

        if iteration == 1:
            # Write initial scores
            scores_file = f"{args.out}.scores.out"
            scoreframe.to_csv(scores_file, sep="\t", index=False)
            print(f"Scores for all initial bins written to: {scores_file}",
                  file=sys.stderr)
            if args.score_only:
                sys.exit(0)

        if len(scoreframe) == 0 or scoreframe.iloc[0]["score"] < args.threshold:
            print(f"No bin surpasses the cutoff of: {args.threshold}",
                  file=sys.stderr)
            break

        # Exclude bins below completeness threshold (they'll never pass)
        low_comp = set(
            scoreframe[scoreframe["Completeness"] < args.threshold]["bin"]
        )
        mag_exclude.update(low_comp)
        scoreframe = scoreframe[
            (~scoreframe["bin"].isin(mag_exclude)) &
            (scoreframe["Contamination"] <= args.max_cont)
        ].reset_index(drop=True)

        if len(scoreframe) == 0:
            break

        # Winner streak: consecutive bins from same binner don't affect each other
        winnerset = scoreframe.iloc[0]["bintool"]
        row_idx = 0

        while (row_idx < len(scoreframe) and
               scoreframe.iloc[row_idx]["bintool"] == winnerset):
            winnerbin = scoreframe.iloc[row_idx]["bin"]
            scoreframe_out.append(scoreframe.iloc[row_idx].to_dict())
            print(f"The highest scoring bin is: {winnerbin}", file=sys.stderr)

            # Add winner's contigs to output
            winner_contigs = contig_to_bin_remain[
                contig_to_bin_remain["bin"] == winnerbin
            ]
            contig_to_bin_out = pd.concat([contig_to_bin_out, winner_contigs],
                                          ignore_index=True)

            print(f"Remaining: {contig_to_bin_remain['contig'].nunique()} "
                  f"contigs and {len(scoreframe)} candidate bins",
                  file=sys.stderr)

            row_idx += 1
            # Break after picking MAGScoT-merged bins (contigs may overlap)
            if str(winnerset).startswith("MAGScoT"):
                break

        # Remove picked bins from scoreframe
        scoreframe = scoreframe.iloc[row_idx:].reset_index(drop=True)

        # Find bins affected by removed contigs → need rescoring
        assigned_contigs = set(contig_to_bin_out["contig"])
        recalc = set(
            contig_to_bin_remain[
                contig_to_bin_remain["contig"].isin(assigned_contigs)
            ]["bin"]
        )

        # Remove assigned contigs and eliminated bins from remaining
        remaining_bins = set(scoreframe["bin"])
        contig_to_bin_remain = contig_to_bin_remain[
            contig_to_bin_remain["bin"].isin(remaining_bins) &
            ~contig_to_bin_remain["contig"].isin(assigned_contigs)
        ].reset_index(drop=True)

        iteration += 1

    # ---- Write output ----
    if scoreframe_out:
        result = pd.DataFrame(scoreframe_out)
        result["refined_id"] = [
            f"{args.out}_{args.bin_separator}_{i + 1:06d}"
            for i in range(len(result))
        ]

        # Map old bin names to new refined IDs in contig output
        bin_rename = dict(zip(result["bin"], result["refined_id"]))
        contig_to_bin_out["binnew"] = contig_to_bin_out["bin"].map(bin_rename)

        n_refined = len(result)
        print(f"Refinement led to a total of {n_refined} bins with a score "
              f">= {args.threshold}", file=sys.stderr)

        stats_file = f"{args.out}.refined.out"
        binning_file = f"{args.out}.refined.contig_to_bin.out"

        print(f"Refinement stats are written to: {stats_file}", file=sys.stderr)
        print(f"Contig-to-refined-bin mapping is written to: {binning_file}",
              file=sys.stderr)

        result.to_csv(stats_file, sep="\t", index=False)
        contig_to_bin_out[["binnew", "contig"]].to_csv(
            binning_file, sep="\t", index=False
        )
    else:
        print(f"No bins with score >= {args.threshold} were found.",
              file=sys.stderr)


if __name__ == "__main__":
    main()
