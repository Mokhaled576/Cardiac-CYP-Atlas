"""Summarize raw-count nucleus QC from explicit cohort/stage inputs."""

import argparse
import csv
import gzip
from pathlib import Path

import numpy as np


def summarize(record, base, chunk_size):
    import h5py
    from anndata.experimental import read_elem
    from anndata.io import sparse_dataset
    from scipy import sparse

    path = base / record["h5ad"]
    with h5py.File(path, "r") as handle:
        matrix = handle[record["count_key"]]
        if isinstance(matrix, h5py.Group):
            matrix = sparse_dataset(matrix)
        var = read_elem(handle["raw/var" if record["count_key"] == "raw/X" else "var"])
        symbol_col = record.get("symbol_column", "")
        symbols = var[symbol_col] if symbol_col else var.index
        mito = np.asarray([str(s).upper().startswith("MT-") for s in symbols])
        if matrix.shape[1] != len(mito):
            raise ValueError(f"{path}: count matrix and feature metadata differ")
        mask = np.ones(matrix.shape[0], dtype=bool)
        flag = record.get("pass_column", "")
        if flag:
            obs = read_elem(handle["obs"])
            if obs[flag].dtype.kind != "b" or obs[flag].isna().any():
                raise ValueError(f"{path}: {flag} must contain only boolean flags")
            mask = obs[flag].to_numpy(dtype=bool)
        n = int(mask.sum())
        if n == 0:
            raise ValueError(f"{path}: selection contains zero nuclei")
        metrics = np.empty((3, n), dtype=np.float64)
        pos = 0
        for start in range(0, matrix.shape[0], chunk_size):
            stop = min(start + chunk_size, matrix.shape[0])
            selected = mask[start:stop]
            if not selected.any():
                continue
            block = matrix[start:stop][selected]
            values = block.data if sparse.issparse(block) else block
            if np.any(~np.isfinite(values) | (values < 0) | (values != np.floor(values))):
                raise ValueError(f"{path}: {record['count_key']} contains non-count values")
            if sparse.issparse(block):
                block = block.tocsr(copy=True)
                block.sum_duplicates()
                block.eliminate_zeros()
                detected = np.diff(block.indptr)
            else:
                detected = np.count_nonzero(block, axis=1)
            totals = np.asarray(block.sum(axis=1, dtype=np.float64)).ravel()
            mt_counts = np.asarray(block[:, mito].sum(axis=1, dtype=np.float64)).ravel()
            mt_pct = np.full(totals.shape, np.nan)
            if mito.any():
                np.divide(mt_counts * 100, totals, out=mt_pct, where=totals > 0)
            end = pos + len(totals)
            metrics[:, pos:end] = (totals, detected, mt_pct)
            pos = end
        assert pos == n
    return summarize_metrics(record, base, path, record["count_key"], len(mask), metrics, int(mito.sum()))


def summarize_export(record, base):
    path = base / record["nucleus_qc"]
    opener = gzip.open if path.suffix == ".gz" else open
    columns = ("total_counts", "n_genes_by_counts", "pct_counts_mt")
    flag = record.get("pass_column", "")
    with opener(path, "rt", encoding="utf-8-sig", newline="") as file:
        reader = csv.DictReader(file, delimiter="\t")
        if not set((*columns, *([flag] if flag else []))).issubset(reader.fieldnames or ()):
            raise ValueError(f"{path}: missing expected per-nucleus QC columns")
        n_total = n_flags = mismatches = 0
        def values():
            nonlocal n_total, n_flags, mismatches
            for row in reader:
                n_total += 1
                passed = row[flag].strip().lower() if flag else "true"
                if passed not in {"true", "false", "1", "0"}:
                    raise ValueError(f"{path}: invalid QC flag at row {n_total}")
                metrics = tuple(float(row[c]) if row[c].strip().lower() not in {"", "na", "nan"}
                                else np.nan for c in columns)
                if "provisional_qc_flag" in row:
                    observed = row["provisional_qc_flag"].strip().lower()
                    if observed not in {"true", "false", "1", "0"}:
                        raise ValueError(f"{path}: invalid provisional flag at row {n_total}")
                    umi, genes, mito = metrics
                    expected = umi > 0 and genes >= 300 and (np.isnan(mito) or mito <= 15)
                    n_flags += 1
                    mismatches += expected != (observed in {"true", "1"})
                if passed in {"true", "1"}:
                    yield metrics
        metrics = np.fromiter(values(), dtype=np.dtype((np.float64, 3))).T
    if not metrics.shape[1]:
        raise ValueError(f"{path}: selection contains zero nuclei")
    if np.any(metrics < 0) or np.any(metrics[2] > 100):
        raise ValueError(f"{path}: invalid counts or mitochondrial percentage")
    result = summarize_metrics(record, base, path, "archived_nucleus_qc_export", n_total, metrics, "NA")
    result["n_flags_compared"] = n_flags
    result["n_provisional_rule_mismatches"] = mismatches if n_flags else "NA"
    if n_flags and mismatches == 0:
        result["qc_thresholds"] = "Archived flags match: UMIs >0, genes >=300, MT <=15% or missing"
    return result


def summarize_metrics(record, base, path, source, n_total, metrics, n_mt_features):
    evidence = record.get("evidence", "").strip()
    documented = any(record.get(key, "").strip() for key in ("qc_thresholds", "doublet_procedure"))
    if documented and (not evidence or not (base / evidence).is_file()):
        raise ValueError(f"{path}: documented thresholds/procedures require a readable evidence file")
    summary = {
        "cohort": record["cohort"], "stage": record["stage"], "input": str(path),
        "source": source, "pass_column": record.get("pass_column", ""), "n_total": n_total,
        "n_selected": metrics.shape[1], "n_mt_features": n_mt_features,
        "qc_thresholds": record.get("qc_thresholds") or "not documented",
        "doublet_procedure": record.get("doublet_procedure") or "not documented",
        "evidence": evidence,
        "n_flags_compared": "NA", "n_provisional_rule_mismatches": "NA",
    }
    for name, values in zip(("n_counts", "n_genes", "pct_counts_mt"), metrics):
        valid = values[np.isfinite(values)]
        stats = np.quantile(valid, [0, .05, .5, 1], method="linear") if len(valid) else ["NA"] * 4
        summary[f"{name}_n_valid"] = len(valid)
        summary.update({f"{name}_{label}": value for label, value in
                        zip(("min", "p05", "median", "max"), stats)})
    return summary


def audit_sample_export(record, base):
    path = base / record["sample_qc"]
    with path.open(encoding="utf-8-sig", newline="") as file:
        reader = csv.DictReader(file, delimiter="\t")
        required = {"nuclei", "provisional_pass"}
        if not required.issubset(reader.fieldnames or ()):
            raise ValueError(f"{path}: missing expected project Sample_QC.tsv columns")
        has_medians = {"median_umi", "median_genes", "median_mito_pct"}.issubset(reader.fieldnames)
        rows = list(reader)
    totals = [(int(row["nuclei"]), int(row["provisional_pass"])) for row in rows]
    if not totals or any(not 0 <= passed <= count for count, passed in totals):
        raise ValueError(f"{path}: empty or invalid sample QC export")
    return {
        "cohort": record["cohort"], "stage": record["stage"], "sample_qc": str(path),
        "n_sample_rows": len(rows), "nuclei_reported": sum(n for n, _ in totals),
        "provisional_pass_reported": sum(n for _, n in totals),
        "has_sample_medians": has_medians,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("manifest", type=Path)
    parser.add_argument("output_dir", type=Path)
    parser.add_argument("--chunk-size", type=int, default=512)
    args = parser.parse_args()
    if args.chunk_size < 1:
        parser.error("--chunk-size must be positive")
    with args.manifest.open(encoding="utf-8-sig", newline="") as file:
        reader = csv.DictReader(file, delimiter="\t")
        if not {"cohort", "stage", "h5ad", "count_key", "sample_qc"}.issubset(reader.fieldnames or ()):
            parser.error("manifest must define cohort, stage, h5ad, count_key and sample_qc")
        rows = list(reader)
    keys = [(row["cohort"], row["stage"]) for row in rows]
    if not keys or len(set(keys)) != len(keys):
        parser.error("provide one row per unique cohort and stage")
    base = args.manifest.resolve().parent
    distributions, sample_audits = [], []
    for row in rows:
        if not any(row.get(key) for key in ("h5ad", "sample_qc", "nucleus_qc")):
            parser.error(f"no input for {row['cohort']} / {row['stage']}")
        if row["h5ad"] and row.get("nucleus_qc"):
            parser.error("use separate stages for h5ad and nucleus_qc inputs")
        if row["h5ad"]:
            distributions.append(summarize(row, base, args.chunk_size))
        if row.get("nucleus_qc"):
            distributions.append(summarize_export(row, base))
        if row["sample_qc"]:
            sample_audits.append(audit_sample_export(row, base))
    args.output_dir.mkdir()
    for name, results in [("Nucleus_QC_Distributions.tsv", distributions), ("Sample_QC_Audit.tsv", sample_audits)]:
        if results:
            with (args.output_dir / name).open("x", encoding="utf-8", newline="") as file:
                writer = csv.DictWriter(file, fieldnames=list(results[0]), delimiter="\t")
                writer.writeheader()
                writer.writerows(results)


if __name__ == "__main__":
    main()
