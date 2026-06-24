#!/usr/bin/env python3
import argparse
import numpy as np

def main():
    p = argparse.ArgumentParser(description="Find medoid frame (min sum of squared distances) after aligning frames.")
    p.add_argument("-s", "--tpr", required=True, help="GROMACS .tpr topology")
    p.add_argument("-f", "--xtc", required=True, help="Trajectory (.xtc)")
    p.add_argument("--select", default="protein and name CA",
                   help='Atom selection for alignment + medoid (default: "protein and name CA")')
    p.add_argument("--stride", type=int, default=1, help="Stride frames (default: 1)")
    p.add_argument("--out", default="medoid_CA.txt", help="Output text file")
    args = p.parse_args()

    import MDAnalysis as mda
    from MDAnalysis.analysis import align

    u = mda.Universe(args.tpr, args.xtc)
    ag = u.select_atoms(args.select)
    if ag.n_atoms == 0:
        raise ValueError(f"Selection matched 0 atoms: {args.select}")

    # Align all frames to the first frame using the same atom selection
    # (in_memory=True is important because we need coords post-alignment)
    align.AlignTraj(u, u, select=args.select, in_memory=True).run()

    coords = []
    times_ps = []

    for ts in u.trajectory[::args.stride]:
        coords.append(ag.positions.astype(np.float64))  # (n_atoms, 3)
        times_ps.append(float(ts.time))                 # usually ps in GROMACS trajectories

    coords = np.asarray(coords)  # (N, n_atoms, 3)
    times_ps = np.asarray(times_ps)

    N = coords.shape[0]
    n_atoms = coords.shape[1]

    # Flatten coordinates so each frame is a vector in R^(3*n_atoms)
    X = coords.reshape(N, 3 * n_atoms)

    # Compute the medoid using the O(N) trick on squared Euclidean distances:
    # S_i = sum_j ||x_i - x_j||^2 = N||x_i||^2 + sum||x_j||^2 - 2 x_i · sum(x_j)
    X_sqnorm = np.sum(X * X, axis=1)          # ||x_i||^2
    sum_sqnorm = float(np.sum(X_sqnorm))      # sum_j ||x_j||^2
    sum_vec = np.sum(X, axis=0)               # sum_j x_j

    S = N * X_sqnorm + sum_sqnorm - 2.0 * (X @ sum_vec)  # shape (N,)

    medoid_idx = int(np.argmin(S))
    medoid_time_ps = float(times_ps[medoid_idx])

    # Convert to RMSD^2 units (nm^2) if coords are in nm (MDAnalysis uses same units as input;
    # GROMACS typically stores nm). This is still "mean squared displacement per atom" after alignment.
    # RMSD^2_i = (1/n_atoms) * (1) * ||x_i - x_j||^2 averaged over j isn't needed here;
    # we only need argmin. But we can report an interpretable scale:
    # average squared distance to others (per atom) for medoid:
    avg_sq_dist_per_atom = float(S[medoid_idx] / (N * n_atoms))

    with open(args.out, "w") as f:
        f.write(f"selection: {args.select}\n")
        f.write(f"stride: {args.stride}\n")
        f.write(f"n_frames_used: {N}\n")
        f.write(f"n_atoms: {n_atoms}\n")
        f.write(f"medoid_index_in_strided_traj: {medoid_idx}\n")
        f.write(f"medoid_time_ps: {medoid_time_ps:.3f}\n")
        f.write(f"medoid_time_ns: {medoid_time_ps/1000.0:.6f}\n")
        f.write(f"avg_sq_dist_per_atom (nm^2, aligned): {avg_sq_dist_per_atom:.6e}\n")

    print(f"Medoid frame (strided index): {medoid_idx}")
    print(f"Medoid time: {medoid_time_ps/1000.0:.6f} ns")
    print(f"Wrote: {args.out}")

if __name__ == "__main__":
    main()
