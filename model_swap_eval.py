#!/usr/bin/env python3
"""
model_swap_eval.py

Evaluates whether replacing the 6 qfix pilot sequences' (bind_022/019/020_
binder, nonb_006/008/009_nb) feature values with their bond-order-fixed
reparameterization changes model performance, versus the current baseline
(original parameterization, exactly as ML_classification.ipynb trains on
today).

Both the baseline and the qfix-swapped run use IDENTICAL StratifiedGroupKFold
splits (same y and seq_groups -- the label and the amino-acid sequence for
these 6 rows are unaffected by the ligand reparameterization, only their MD-
derived feature values change), so this is a paired before/after comparison
on the exact same folds, not two independently-resampled runs.

Mirrors ML_classification.ipynb's feature-merge (cell 1), GroupAwareSelector
(cell 2), and CV pipeline (cell 5) exactly -- see that notebook for the
rationale behind each feature family / hyperparameter choice.

Usage:
    python model_swap_eval.py --qfix_table qfix_pilot_feat_table.csv
"""
import argparse
import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import StratifiedGroupKFold, cross_validate
from sklearn.metrics import make_scorer, f1_score, precision_score, recall_score, roc_auc_score
from sklearn.pipeline import Pipeline
from sklearn.base import BaseEstimator, TransformerMixin
from sklearn.preprocessing import MinMaxScaler
from scipy.stats import spearmanr
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform

from seq_utils import sequence_similarity_groups

POCKET_RESIDUES = [58, 59, 62, 83, 87, 88, 89, 92, 110, 115, 116, 117, 120, 122, 159, 160, 163, 164]
RMSD_REGIONS = ["Gate", "Latch"]

FLIPPED_LIGAND_NAMES = [
    "pair_3064_binder", "bind_033_binder", "bind_109_binder",
    "pair_1708_low_pkt", "nonb_012_nb", "nonb_020_nb", "nonb_055_nb",
]

TARGET_SEQ_IDS = ["bind_022_binder", "bind_019_binder", "bind_020_binder",
                   "nonb_006_nb", "nonb_008_nb", "nonb_009_nb"]

N_SPLITS = 7
RANDOM_STATE = 42


# ── GroupAwareSelector -- verbatim copy of ML_classification.ipynb cell 2 ────
def _feature_score(x_col, y):
    classes = np.unique(y)
    if len(classes) <= 2:
        y_bin = (y == classes[-1]).astype(int)
        if y_bin.sum() in (0, len(y_bin)):
            return 0.0
        return abs(roc_auc_score(y_bin, x_col) - 0.5)
    best = 0.0
    for c in classes:
        y_bin = (y == c).astype(int)
        if y_bin.sum() in (0, len(y_bin)):
            continue
        best = max(best, abs(roc_auc_score(y_bin, x_col) - 0.5))
    return best


class GroupAwareSelector(BaseEstimator, TransformerMixin):
    def __init__(self, groups, corr_prune_groups=frozenset(), k_per_group=None,
                 corr_threshold=0.8):
        self.groups = groups
        self.corr_prune_groups = corr_prune_groups
        self.k_per_group = k_per_group or {}
        self.corr_threshold = corr_threshold

    def fit(self, X, y):
        X = np.asarray(X)
        y = np.asarray(y)
        selected = []
        self.group_report_ = {}

        for gname, idx in self.groups.items():
            idx = np.array(idx)
            variances = X[:, idx].var(axis=0)
            valid = idx[variances > 1e-12]
            if len(valid) == 0:
                self.group_report_[gname] = {"input": len(idx), "after_corr": 0, "kept": 0}
                continue

            if gname in self.corr_prune_groups and len(valid) > 1:
                Xg = X[:, valid]
                corr, _ = spearmanr(Xg)
                corr = np.atleast_2d(corr)
                if corr.shape[0] == 1:
                    corr = np.array([[1.0, corr.item()], [corr.item(), 1.0]])
                dist = 1.0 - np.abs(np.nan_to_num(corr, nan=0.0))
                np.fill_diagonal(dist, 0.0)
                dist = (dist + dist.T) / 2
                Z = linkage(squareform(dist, checks=False), method="average")
                cluster_ids = fcluster(Z, t=(1 - self.corr_threshold), criterion="distance")
                reps = []
                for cid in np.unique(cluster_ids):
                    members = valid[cluster_ids == cid]
                    scores = [_feature_score(X[:, m], y) for m in members]
                    reps.append(members[int(np.argmax(scores))])
                pruned = np.array(reps)
            else:
                pruned = valid

            k = self.k_per_group.get(gname)
            if k is not None and len(pruned) > k:
                scores = np.array([_feature_score(X[:, m], y) for m in pruned])
                top_k = pruned[np.argsort(scores)[-k:]]
            else:
                top_k = pruned

            self.group_report_[gname] = {
                "input": len(idx), "after_corr": len(pruned), "kept": len(top_k)
            }
            selected.extend(top_k.tolist())

        self.selected_indices_ = np.array(sorted(selected))
        return self

    def transform(self, X):
        return np.asarray(X)[:, self.selected_indices_]


def load_baseline_df():
    """Reproduces ML_classification.ipynb cell 1's merged df exactly
    (SEQ_SOURCE='all'), reading the same repo-relative paths it does."""
    df = pd.read_excel("feat_table_500ns.xlsx", sheet_name="all_feats_500ns")
    df = df[~df["name"].isin(FLIPPED_LIGAND_NAMES)].reset_index(drop=True)
    n_base = len(df)

    gl_cols = []
    for r in RMSD_REGIONS:
        gl_cols += [f"{r} RMSD mean (A)", f"{r} RMSD SD (A)", f"{r} drift100 (A)", f"{r} slope (A/ns)"]
    gl = pd.read_csv("analysis/gate_latch_rmsd_to_ref_summary_500ns.csv")
    gl = gl.rename(columns={"Sequence": "seq_id"})[["seq_id"] + gl_cols]
    df = df.merge(gl, left_on="name", right_on="seq_id", how="left").drop(columns=["seq_id"])

    sb_cols = ["max_saltbridge_occupancy_pct", "n_saltbridges_gt50pct", "mean_top3_occupancy_pct"]
    sb = pd.read_csv("salt_bridge/saltbridge_features_all_seqs.csv")[["seq_id"] + sb_cols]
    df = df.merge(sb, left_on="name", right_on="seq_id", how="left").drop(columns=["seq_id"])
    df[sb_cols] = df[sb_cols].fillna(0.0)

    core = pd.read_csv("water_analysis/agg_out/dw_scores_all_sequences_40_500ns_core.csv")
    tail = pd.read_csv("water_analysis/agg_out/dw_scores_all_sequences_40_500ns_tail.csv")
    dw_names = [f"D_{r}" for r in POCKET_RESIDUES] + [f"W_{r}" for r in POCKET_RESIDUES]
    core_sub = core[["seq_id"] + dw_names].rename(columns={c: f"{c}_core" for c in dw_names})
    tail_sub = tail[["seq_id"] + dw_names].rename(columns={c: f"{c}_tail" for c in dw_names})
    ct = core_sub.merge(tail_sub, on="seq_id", how="inner")
    ct_cols = []
    for r in POCKET_RESIDUES:
        for letter in ("D", "W"):
            col = f"delta_{letter}_{r}"
            ct[col] = ct[f"{letter}_{r}_core"] - ct[f"{letter}_{r}_tail"]
            ct_cols.append(col)
    df = df.merge(ct[["seq_id"] + ct_cols], left_on="name", right_on="seq_id", how="left").drop(columns=["seq_id"])
    df[ct_cols] = df[ct_cols].fillna(0.0)

    hyd_cols = ["hydration_count_pocket_4A_mean", "hydration_count_pocket_4A_std",
                "hydration_count_pocket_4A_min", "hydration_count_pocket_4A_max",
                "hydration_count_pocket_4A_early20_mean", "hydration_count_pocket_4A_late20_mean",
                "hydration_count_pocket_4A_drift20", "hydration_count_pocket_4A_slope",
                "hydration_count_pocket_4A_slope_per_ns"]
    hyd = pd.read_csv("water_spatial/water_density_feats_pocket.csv")[["seq_id"] + hyd_cols]
    df = df.merge(hyd, left_on="name", right_on="seq_id", how="left").drop(columns=["seq_id"])

    wb_cols = ["gate_bridge_occupancy", "triple_bridge_occupancy",
               "mean_n_triple_bridge_waters", "co_occurrence_occupancy",
               "mean_run_duration_ns"]
    wb = pd.read_csv("water_analysis/gate_latch_bridge_all_0_500ns_core.csv")[["seq_id"] + wb_cols]
    df = df.merge(wb, left_on="name", right_on="seq_id", how="left").drop(columns=["seq_id"])
    df[wb_cols] = df[wb_cols].fillna(0.0)

    dw_cols = []
    for r in POCKET_RESIDUES:
        dw_cols += [f"D_{r}", f"W_{r}"]

    feature_group_cols = {
        "dw_pocket":        dw_cols,
        "contact_type":     ["mean_n_neg_charged", "std_n_neg_charged", "occ_n_neg_charged_gt0"],
        "rmsf":             ["Y23 RMSF (A)", "R79 RMSF (A)", "I110 RMSF (A)", "G163 RMSF (A)",
                              "Gate (r84-90) mean (A)", "Gate (r84-90) SD (A)"],
        "gate_latch_rmsd":  gl_cols,
        "salt_bridge":      sb_cols,
        "core_tail_delta":  ct_cols,
        "hydration_pocket": hyd_cols,
        "water_bridge":     wb_cols,
    }

    print(f"Baseline cohort: {n_base} sequences ({len(FLIPPED_LIGAND_NAMES)} flipped-ligand QC exclusions removed)")
    return df, feature_group_cols


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--qfix_table", default="qfix_pilot_feat_table.csv")
    p.add_argument("--out", default="model_swap_eval_results.csv")
    args = p.parse_args()

    df, feature_group_cols = load_baseline_df()
    feature_cols = [c for cols in feature_group_cols.values() for c in cols]
    col_index = {c: i for i, c in enumerate(feature_cols)}
    feature_groups = {g: [col_index[c] for c in cols] for g, cols in feature_group_cols.items()}
    corr_prune_groups = {"dw_pocket", "water_bridge"}
    k_per_group = {"dw_pocket": 12, "core_tail_delta": 16}

    missing_targets = set(TARGET_SEQ_IDS) - set(df["name"])
    if missing_targets:
        raise ValueError(f"Target sequences not found in baseline cohort: {missing_targets}")

    y = df["Label"].values
    sequences = df["Sequence"].astype(str).str.strip().tolist()
    seq_groups = sequence_similarity_groups(sequences, identity_threshold=0.95)

    X_base = df[feature_cols].values

    qfix = pd.read_csv(args.qfix_table).set_index("name")
    missing_qfix_cols = set(feature_cols) - set(qfix.columns)
    if missing_qfix_cols:
        raise ValueError(f"{args.qfix_table} is missing feature columns: {missing_qfix_cols}")

    df_qfix = df.copy()
    for seq_id in TARGET_SEQ_IDS:
        row_idx = df_qfix.index[df_qfix["name"] == seq_id]
        df_qfix.loc[row_idx, feature_cols] = qfix.loc[seq_id, feature_cols].values
    X_qfix = df_qfix[feature_cols].values

    # np.isclose(..., equal_nan=True) rather than a plain != comparison --
    # several baseline rows have NaN in the gate_latch_rmsd columns (that
    # source CSV doesn't cover the full cohort), and NaN != NaN is always
    # True, which would otherwise flag every NaN-containing row as "changed"
    # even though the swap never touched it.
    n_changed = (~np.isclose(X_base.astype(float), X_qfix.astype(float),
                              equal_nan=True)).any(axis=1).sum()
    print(f"Rows with changed features in the qfix-swapped table: {n_changed} "
          f"(expected {len(TARGET_SEQ_IDS)})")

    cv = StratifiedGroupKFold(n_splits=N_SPLITS, shuffle=True, random_state=RANDOM_STATE)
    pipe = Pipeline([
        ("scale",  MinMaxScaler()),
        ("select", GroupAwareSelector(groups=feature_groups, corr_prune_groups=corr_prune_groups,
                                       k_per_group=k_per_group)),
        ("rf",     RandomForestClassifier(n_estimators=500, max_features="sqrt",
                                           class_weight="balanced", random_state=RANDOM_STATE)),
    ])
    scoring = {
        "accuracy":          "accuracy",
        "balanced_accuracy": "balanced_accuracy",
        "roc_auc":           "roc_auc",
        "f1":                make_scorer(f1_score,        zero_division=0),
        "precision":         make_scorer(precision_score, zero_division=0),
        "recall":            make_scorer(recall_score,    zero_division=0),
    }

    # Fold assignment depends only on y/groups (both identical between the
    # two runs, unaffected by the ligand reparameterization), not on X's
    # values -- so both cv.split() calls below yield IDENTICAL folds, making
    # this a paired before/after comparison rather than two independent runs.
    print(f"\nRunning baseline CV ({N_SPLITS}-fold grouped stratified)...")
    res_base = cross_validate(pipe, X_base, y, cv=cv, groups=seq_groups, scoring=scoring)

    print(f"Running qfix-swapped CV ({N_SPLITS}-fold grouped stratified, same folds)...")
    res_qfix = cross_validate(pipe, X_qfix, y, cv=cv, groups=seq_groups, scoring=scoring)

    rows = []
    print(f"\n{'metric':<20}{'baseline':>18}{'qfix-swapped':>18}{'mean paired diff':>20}")
    for metric in scoring:
        b = res_base[f"test_{metric}"]
        q = res_qfix[f"test_{metric}"]
        diff = q - b
        print(f"{metric:<20}{b.mean():>10.4f} ± {b.std():<5.4f}"
              f"{q.mean():>10.4f} ± {q.std():<5.4f}{diff.mean():>+20.4f}")
        rows.append({
            "metric": metric,
            "baseline_mean": b.mean(), "baseline_std": b.std(),
            "qfix_mean": q.mean(), "qfix_std": q.std(),
            "mean_paired_diff": diff.mean(), "std_paired_diff": diff.std(),
        })

    pd.DataFrame(rows).to_csv(args.out, index=False)
    print(f"\nSaved -> {args.out}")


if __name__ == "__main__":
    main()
