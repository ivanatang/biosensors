"""
ml_utils.py

Shared CV/model-selection pieces for comparing the binder/nonbinder RF
classifier across feature tables built at different MD trajectory durations
(feat_table_100ns.xlsx / feat_table_250ns.xlsx / feat_table_500ns.xlsx).

This is a verbatim copy (not a refactor) of the GroupAwareSelector,
_feature_score, grouped_permutation_test, CORR_PRUNE_GROUPS, and K_PER_GROUP
definitions from ML_classification.ipynb's cell 2, plus the feature-group
column layout from cell 1 and the CV/pipeline setup from cell 5. Kept as a
copy rather than importing ML_classification.ipynb so the primary,
already-working notebook is never modified or re-pointed at this module.

Used by ML_duration_comparison.ipynb.
"""
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

# ── Feature-group layout (mirrors ML_classification.ipynb cell 1) ────────────
POCKET_RESIDUES = [58, 59, 62, 83, 87, 88, 89, 92, 110, 115, 116, 117, 120, 122, 159, 160, 163, 164]
RMSD_REGIONS = ["Gate", "Latch"]

CORR_PRUNE_GROUPS = {"dw_pocket", "water_bridge"}
K_PER_GROUP = {"dw_pocket": 12, "core_tail_delta": 16}

N_SPLITS = 7
RANDOM_STATE = 42

SCORING = {
    "accuracy":          "accuracy",
    "balanced_accuracy": "balanced_accuracy",
    "roc_auc":           "roc_auc",
    "f1":                make_scorer(f1_score,        zero_division=0),
    "precision":         make_scorer(precision_score, zero_division=0),
    "recall":            make_scorer(recall_score,    zero_division=0),
}


def build_feature_group_cols():
    """Same 8-family column layout as ML_classification.ipynb's
    FEATURE_GROUP_COLS, valid for any duration since build_feat_table.py
    reproduces this exact column naming for every window."""
    dw_cols = []
    for r in POCKET_RESIDUES:
        dw_cols += [f"D_{r}", f"W_{r}"]

    gl_cols = []
    for r in RMSD_REGIONS:
        gl_cols += [f"{r} RMSD mean (A)", f"{r} RMSD SD (A)", f"{r} drift100 (A)", f"{r} slope (A/ns)"]

    ct_cols = []
    for r in POCKET_RESIDUES:
        for letter in ("D", "W"):
            ct_cols.append(f"delta_{letter}_{r}")

    return {
        "dw_pocket":        dw_cols,
        "contact_type":     ["mean_n_neg_charged", "std_n_neg_charged", "occ_n_neg_charged_gt0"],
        "rmsf":             ["Y23 RMSF (A)", "R79 RMSF (A)", "I110 RMSF (A)", "G163 RMSF (A)",
                              "Gate (r84-90) mean (A)", "Gate (r84-90) SD (A)"],
        "gate_latch_rmsd":  gl_cols,
        "salt_bridge":      ["max_saltbridge_occupancy_pct", "n_saltbridges_gt50pct", "mean_top3_occupancy_pct"],
        "core_tail_delta":  ct_cols,
        "hydration_pocket": ["hydration_count_pocket_4A_mean", "hydration_count_pocket_4A_std",
                              "hydration_count_pocket_4A_min", "hydration_count_pocket_4A_max",
                              "hydration_count_pocket_4A_early20_mean", "hydration_count_pocket_4A_late20_mean",
                              "hydration_count_pocket_4A_drift20", "hydration_count_pocket_4A_slope",
                              "hydration_count_pocket_4A_slope_per_ns"],
        "water_bridge":     ["gate_bridge_occupancy", "triple_bridge_occupancy",
                              "mean_n_triple_bridge_waters", "co_occurrence_occupancy",
                              "mean_run_duration_ns"],
    }


def load_feat_table(path, sheet_name):
    """Load a self-contained feat_table_{N}ns.xlsx built by build_feat_table.py."""
    df = pd.read_excel(path, sheet_name=sheet_name)
    feature_group_cols = build_feature_group_cols()
    feature_cols = [c for cols in feature_group_cols.values() for c in cols]
    return df, feature_cols, feature_group_cols


# ── GroupAwareSelector: fold-safe transformer (verbatim from cell 2) ─────────
def _feature_score(x_col, y):
    """|AUC - 0.5| of a single feature vs y; for >2 classes, best one-vs-rest AUC."""
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
    """
    groups            : {group_name: [column indices into X]}
    corr_prune_groups : group names to cluster-prune by |Spearman r| before budgeting
    k_per_group       : {group_name: k or None} cap after pruning
    corr_threshold    : |corr| above this -> same cluster (one representative kept)
    """
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
                if corr.shape[0] == 1:   # spearmanr collapses to scalar for exactly 2 cols
                    corr = np.array([[1.0, corr.item()], [corr.item(), 1.0]])
                dist = 1.0 - np.abs(np.nan_to_num(corr, nan=0.0))
                np.fill_diagonal(dist, 0.0)
                dist = (dist + dist.T) / 2   # enforce symmetry against fp noise
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


def grouped_permutation_test(estimator, X, y, groups, cv, scoring="balanced_accuracy",
                              n_permutations=100, random_state=None):
    """Shuffles y completely freely each permutation (not just within-group,
    which would be close to a no-op given how many singleton sequence-
    similarity groups exist), while still using `groups` for fold-safe CV
    splitting on every permuted run."""
    rng = np.random.default_rng(random_state)

    observed = cross_validate(estimator, X, y, cv=cv, groups=groups,
                               scoring=scoring)["test_score"].mean()

    perm_scores = np.empty(n_permutations)
    for i in range(n_permutations):
        y_perm = rng.permutation(y)
        perm_scores[i] = cross_validate(estimator, X, y_perm, cv=cv, groups=groups,
                                         scoring=scoring)["test_score"].mean()

    pval = (1 + np.sum(perm_scores >= observed)) / (1 + n_permutations)
    return observed, perm_scores, pval


def make_pipeline(feature_group_cols, feature_cols):
    col_index = {c: i for i, c in enumerate(feature_cols)}
    feature_groups = {g: [col_index[c] for c in cols] for g, cols in feature_group_cols.items()}
    return Pipeline([
        ("scale",  MinMaxScaler()),
        ("select", GroupAwareSelector(groups=feature_groups, corr_prune_groups=CORR_PRUNE_GROUPS,
                                       k_per_group=K_PER_GROUP)),
        ("rf",     RandomForestClassifier(n_estimators=500, max_features="sqrt",
                                           class_weight="balanced", random_state=RANDOM_STATE)),
    ])


def run_cv(df, feature_cols, feature_group_cols, n_splits=N_SPLITS, random_state=RANDOM_STATE):
    """Runs the same StratifiedGroupKFold CV as ML_classification.ipynb cell 5
    and returns (metrics_df, cv_results, seq_groups)."""
    X = df[feature_cols].values
    y = df["Label"].values
    sequences = df["Sequence"].astype(str).str.strip().tolist()
    seq_groups = sequence_similarity_groups(sequences, identity_threshold=0.95)

    cv = StratifiedGroupKFold(n_splits=n_splits, shuffle=True, random_state=random_state)
    pipe = make_pipeline(feature_group_cols, feature_cols)

    cv_results = cross_validate(pipe, X, y, cv=cv, groups=seq_groups,
                                 scoring=SCORING, return_train_score=False)
    metrics = {k.replace("test_", ""): v for k, v in cv_results.items() if k.startswith("test_")}

    metrics_df = pd.DataFrame({
        k: {"Mean": round(v.mean(), 3), "SD": round(v.std(), 3)}
        for k, v in metrics.items()
    }).T

    return metrics_df, metrics, seq_groups
