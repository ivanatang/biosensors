#!/usr/bin/env python3
import argparse
import numpy as np

def main():
    p = argparse.ArgumentParser(
        description="Find medoid frame (min sum of squared distances) after aligning frames."
    )
    p.add_argument("-s", "--tpr",      required=True, help="GROMACS .tpr topology")
    p.add_argument("-f", "--xtc",      required=True, help="Trajectory (.xtc)")
    p.add_argument("--select", default="protein and name CA",
                   help='Atom selection for alignment + medoid (default: "protein and name CA")')
    p.add_argument("--stride",   type=int,   default=1,    help="Stride frames (default: 1)")
    p.add_argument("--start-ns", type=float, default=None,
                   help="Start of analysis window in ns (default: beginning of trajectory)")
    p.add_argument("--end-ns",   type=float, default=None,
                   help="End of analysis window in ns (default: end of trajectory)")
    p.add_argument("--out", default="medoid_CA.txt", help="Output text file")
    args = p.parse_args()

    import MDAnalysis as mda
    from MDAnalysis.analysis import align

    u  = mda.Universe(args.tpr, args.xtc)
    ag = u.select_atoms(args.select)
    if ag.n_atoms == 0:
        raise ValueError(f"Selection matched 0 atoms: {args.select}")

    # ── Read trajectory metadata BEFORE alignment resets timestamps ──
    # AlignTraj with in_memory=True replaces the trajectory with a
    # MemoryReader that resets ts.time to start from 0. We therefore
    # capture the actual simulation times from frame indices now.
    dt            = u.trajectory.dt        # ps per frame
    n_total       = len(u.trajectory)
    first_time_ps = u.trajectory[0].time   # actual simulation time of frame 0 (e.g. 40012.5 ps)

    # ── Resolve time window → frame indices ───────────────────────
    start_frame = (max(0, round((args.start_ns * 1000.0 - first_time_ps) / dt))
                   if args.start_ns is not None else 0)
    end_frame   = (min(n_total, round((args.end_ns * 1000.0 - first_time_ps) / dt))
                   if args.end_ns   is not None else n_total)

    if end_frame <= start_frame:
        raise ValueError(
            f"Empty window: start_frame={start_frame}, end_frame={end_frame} "
            f"(start_ns={args.start_ns}, end_ns={args.end_ns}, "
            f"first_time_ps={first_time_ps}, dt={dt} ps)"
        )

    # Pre-compute actual simulation times for every window frame
    # (using frame index arithmetic, not ts.time, to avoid MemoryReader reset)
    window_frame_indices = np.arange(start_frame, end_frame, args.stride)
    actual_times_ps      = first_time_ps + window_frame_indices * dt

    actual_start_ns = actual_times_ps[0]  / 1000.0
    actual_end_ns   = actual_times_ps[-1] / 1000.0

    print(f"Trajectory first frame : {first_time_ps/1000.0:.4f} ns  "
          f"({n_total} frames, dt={dt:.2f} ps)")
    print(f"Analysis window        : {actual_start_ns:.4f}–{actual_end_ns:.4f} ns  "
          f"(frames {start_frame}–{end_frame})")

    # ── Align ALL frames to first frame ───────────────────────────
    # Full trajectory alignment for a consistent reference frame.
    # Timestamps are reset to 0-based in the resulting MemoryReader;
    # we use actual_times_ps (computed above) for all time reporting.
    align.AlignTraj(u, u, select=args.select, in_memory=True).run()

    # ── Collect aligned coordinates for the window only ───────────
    coords = []
    for ts in u.trajectory[start_frame:end_frame:args.stride]:
        coords.append(ag.positions.astype(np.float64))   # (n_atoms, 3)

    coords = np.asarray(coords)   # (N, n_atoms, 3)
    N      = coords.shape[0]
    n_atoms = coords.shape[1]

    print(f"Frames collected for medoid: {N}")

    # ── Medoid via O(N) squared-distance trick ─────────────────────
    # S_i = sum_j ||x_i - x_j||^2 = N||x_i||^2 + sum||x_j||^2 - 2 x_i · sum(x_j)
    X          = coords.reshape(N, 3 * n_atoms)
    X_sqnorm   = np.sum(X * X, axis=1)
    sum_sqnorm = float(np.sum(X_sqnorm))
    sum_vec    = np.sum(X, axis=0)
    S          = N * X_sqnorm + sum_sqnorm - 2.0 * (X @ sum_vec)

    medoid_idx           = int(np.argmin(S))
    medoid_time_ps       = float(actual_times_ps[medoid_idx])   # actual simulation time
    avg_sq_dist_per_atom = float(S[medoid_idx] / (N * n_atoms))

    # ── Write output ───────────────────────────────────────────────
    with open(args.out, "w") as f:
        f.write(f"selection: {args.select}\n")
        f.write(f"stride: {args.stride}\n")
        f.write(f"traj_first_time_ps: {first_time_ps:.3f}\n")
        f.write(f"window_start_ns: {actual_start_ns:.6f}\n")
        f.write(f"window_end_ns: {actual_end_ns:.6f}\n")
        f.write(f"n_frames_used: {N}\n")
        f.write(f"n_atoms: {n_atoms}\n")
        f.write(f"medoid_index_in_strided_traj: {medoid_idx}\n")
        f.write(f"medoid_time_ps: {medoid_time_ps:.3f}\n")
        f.write(f"medoid_time_ns: {medoid_time_ps/1000.0:.6f}\n")
        f.write(f"avg_sq_dist_per_atom (nm^2, aligned): {avg_sq_dist_per_atom:.6e}\n")

    print(f"Medoid frame (strided index): {medoid_idx}")
    print(f"Medoid time: {medoid_time_ps/1000.0:.6f} ns  ({medoid_time_ps:.3f} ps)")
    print(f"Wrote: {args.out}")

if __name__ == "__main__":
    main()
