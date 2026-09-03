"""Shared sequence-comparison utilities (amino acid distance/similarity, variable positions).

Extracted from sequence_similarity.ipynb so ML_classification.ipynb can reuse the same
pairwise sequence-distance logic to build sequence-similarity groups for GroupKFold CV.
"""
import numpy as np

AA_ORDER = list("ACDEFGHIKLMNPQRSTVWY")
AA_INDEX = {aa: i for i, aa in enumerate(AA_ORDER)}

# BLOSUM62 substitution matrix (symmetric, 20x20, indexed by AA_ORDER)
# Source: Henikoff & Henikoff (1992) PNAS 89:10915
_BLOSUM62_RAW = [
# A   C   D   E   F   G   H   I   K   L   M   N   P   Q   R   S   T   V   W   Y
[ 4,  0, -2, -1, -2,  0, -2, -1, -1, -1, -1, -2, -1, -1, -1,  1,  0,  0, -3, -2],  # A
[ 0,  9, -3, -4, -2, -3, -3, -1, -3, -1, -1, -3, -3, -3, -3, -1, -1, -1, -2, -2],  # C
[-2, -3,  6,  2, -3, -1, -1, -3, -1, -4, -3,  1, -1,  0, -2,  0, -1, -3, -4, -3],  # D
[-1, -4,  2,  5, -3, -2,  0, -3,  1, -3, -2,  0, -1,  2,  0,  0, -1, -2, -3, -2],  # E
[-2, -2, -3, -3,  6, -3, -1,  0, -3,  0,  0, -3, -4, -3, -3, -2, -2, -1,  1,  3],  # F
[ 0, -3, -1, -2, -3,  6, -2, -4, -2, -4, -3,  0, -2, -2, -2,  0, -2, -3, -2, -3],  # G
[-2, -3, -1,  0, -1, -2,  8, -3, -1, -3, -2,  1, -2,  0,  0, -1, -2, -3, -2,  2],  # H
[-1, -1, -3, -3,  0, -4, -3,  4, -3,  2,  1, -3, -3, -3, -3, -2, -1,  3, -3, -1],  # I
[-1, -3, -1,  1, -3, -2, -1, -3,  5, -2, -1,  0, -1,  1,  2,  0, -1, -2, -3, -2],  # K
[-1, -1, -4, -3,  0, -4, -3,  2, -2,  4,  2, -3, -3, -2, -2, -2, -1,  1, -2, -1],  # L
[-1, -1, -3, -2,  0, -3, -2,  1, -1,  2,  5, -2, -2,  0, -1, -1, -1,  1, -1, -1],  # M
[-2, -3,  1,  0, -3,  0,  1, -3,  0, -3, -2,  6, -2,  0,  0,  1,  0, -3, -4, -2],  # N
[-1, -3, -1, -1, -4, -2, -2, -3, -1, -3, -2, -2,  7, -1, -2, -1, -1, -2, -4, -3],  # P
[-1, -3,  0,  2, -3, -2,  0, -3,  1, -2,  0,  0, -1,  5,  1,  0, -1, -2, -2, -1],  # Q
[-1, -3, -2,  0, -3, -2,  0, -3,  2, -2, -1,  0, -2,  1,  5, -1, -1, -3, -3, -2],  # R
[ 1, -1,  0,  0, -2,  0, -1, -2,  0, -2, -1,  1, -1,  0, -1,  4,  1, -2, -3, -2],  # S
[ 0, -1, -1, -1, -2, -2, -2, -1, -1, -1, -1,  0, -1, -1, -1,  1,  5,  0, -2, -2],  # T
[ 0, -1, -3, -2, -1, -3, -3,  3, -2,  1,  1, -3, -2, -2, -3, -2,  0,  4, -3, -1],  # V
[-3, -2, -4, -3,  1, -2, -2, -3, -3, -2, -1, -4, -4, -2, -3, -3, -2, -3, 11,  2],  # W
[-2, -2, -3, -2,  3, -3,  2, -1, -2, -1, -1, -2, -3, -1, -2, -2, -2, -1,  2,  7],  # Y
]
BLOSUM62 = np.array(_BLOSUM62_RAW, dtype=np.float32)
# Max self-score per AA for normalisation
BLOSUM62_SELFSCORES = np.array([BLOSUM62[i, i] for i in range(20)])


def variable_positions(sequences, min_entropy=0.1):
    """Finds sequence positions with enough amino-acid diversity to matter.

    Args:
        sequences (list[str]): Equal-length amino-acid sequences.
        min_entropy (float): Minimum Shannon entropy (bits) at a position
            to be kept (default: 0.1).

    Returns:
        numpy.ndarray: 0-indexed positions with entropy >= min_entropy.
    """
    L = len(sequences[0])
    entropies = []
    for pos in range(L):
        counts = np.zeros(20)
        for seq in sequences:
            idx = AA_INDEX.get(seq[pos].upper())
            if idx is not None: counts[idx] += 1
        freq = counts / counts.sum() if counts.sum() > 0 else counts
        freq = freq[freq > 0]
        entropies.append(-np.sum(freq * np.log2(freq)) if len(freq) else 0)
    return np.where(np.array(entropies) >= min_entropy)[0]


def hamming_distance_matrix(sequences):
    """Computes pairwise fractional Hamming distance between sequences.

    Args:
        sequences (list[str]): Equal-length amino-acid sequences.

    Returns:
        numpy.ndarray: (n, n) matrix with values in [0, 1].
    """
    n, L = len(sequences), len(sequences[0])
    S = np.array([[AA_INDEX.get(aa.upper(), -1) for aa in seq]
                  for seq in sequences], dtype=np.int8)
    D = np.zeros((n, n), dtype=np.float32)
    for i in range(n):
        for j in range(i+1, n):
            d = np.mean(S[i] != S[j])
            D[i, j] = D[j, i] = d
    return D


def blosum62_similarity_matrix(sequences):
    """Computes pairwise BLOSUM62-based similarity, normalized to [0, 1].

    Normalization: score(i,j) / sqrt(score(i,i) * score(j,j)), so identical
    sequences score 1.0 and divergent ones score lower.

    Args:
        sequences (list[str]): Equal-length amino-acid sequences.

    Returns:
        numpy.ndarray: (n, n) similarity matrix, values in [0, 1].
    """
    n, L = len(sequences), len(sequences[0])
    S = np.array([[AA_INDEX.get(aa.upper(), -1) for aa in seq]
                  for seq in sequences], dtype=np.int32)

    def seq_score(a, b):
        total = 0
        for pos in range(L):
            ia, ib = a[pos], b[pos]
            if ia >= 0 and ib >= 0:
                total += BLOSUM62[ia, ib]
        return total

    self_scores = np.array([seq_score(S[i], S[i]) for i in range(n)])

    sim = np.zeros((n, n), dtype=np.float32)
    for i in range(n):
        for j in range(i, n):
            raw = seq_score(S[i], S[j])
            denom = np.sqrt(self_scores[i] * self_scores[j])
            val = raw / denom if denom > 0 else 0.0
            sim[i, j] = sim[j, i] = val
    return sim


def sequence_similarity_groups(sequences, identity_threshold=0.95, min_entropy=0.1):
    """Clusters sequences into groups of near-identical variants.

    Sequences >= identity_threshold identical (fraction of matching
    residues on the variable positions) fall into the same group. Used to
    build GroupKFold groups so near-identical designed variants don't
    split across train/test folds.

    Args:
        sequences (list[str]): Equal-length amino-acid sequences.
        identity_threshold (float): Minimum fractional identity (on
            variable positions) for two sequences to share a group
            (default: 0.95).
        min_entropy (float): Passed to variable_positions (default: 0.1).

    Returns:
        numpy.ndarray: (n,) array of integer group ids.
    """
    from scipy.cluster.hierarchy import fcluster, linkage
    from scipy.spatial.distance import squareform

    var_pos = variable_positions(sequences, min_entropy=min_entropy)
    var_sequences = ["".join(seq[p] for p in var_pos) for seq in sequences]
    D = hamming_distance_matrix(var_sequences)
    np.fill_diagonal(D, 0.0)
    Z = linkage(squareform(D, checks=False), method="average")
    max_dist = 1.0 - identity_threshold
    return fcluster(Z, t=max_dist, criterion="distance")
