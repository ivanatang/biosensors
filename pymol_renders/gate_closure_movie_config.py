"""Per-movie settings shared by gate_latch_closure_ghost_movie.py (PyMOL, 3D
frames) and gate_latch_trace_panel.py (matplotlib, trace strip + stitch).

Pick a movie with GATE_MOVIE=<key> (default DEFAULT_MOVIE). Both scripts must
see the same key so frame times, output names, and the plotted distance agree.

Window trajectories and full-resolution traces live in the gitignored
pymol_renders/scratch_traj/; the commands that made them are listed per entry.

Frame timing: each window xtc's raw frame 1 is at t0_ns, spaced 37.5 ps.
PyMOL's load_traj(start=1, interval=STRIDE) keeps raw frames STRIDE, 2*STRIDE,
... (it skips frame 1), so state s is at t0_ns + (s*STRIDE - 1) * 0.0375 ns.
"""

import os

REPO_ROOT = "/Users/ivanatang/Developer/biosensors"
SCRATCH = os.path.join(REPO_ROOT, "pymol_renders", "scratch_traj")
BINDERS = "/Users/ivanatang/Library/CloudStorage/OneDrive-UCB-O365/Shirts Lab/LCA_boltz_models/binders"
RUN_SUBDIR = "prod_md_0p9_cutoff_3dt_64x1_16PME_642dd"
TIME_SPACING_NS = 0.0375

MOVIES = {
    # Ca88 -> Ca116 (run_gate_latch.sh's gate-latch metric). Largest sustained
    # Ca88-Ca116 drop among binders: ~1.45 -> ~0.80 nm over 444-450 ns, closed
    # to 500 ns. A re-closure (gate starts closed, drifts open ~125-400 ns).
    # xtc: gmx trjconv -f PL_only_40_500ns.xtc -b 440000 -e 460000
    "pair3088_ca116": dict(
        seq="pair_3088_binder",
        xtc="gate_closure_pair3088_440_460ns_PL.xtc",
        t0_ns=440.025,
        n_raw_frames=533,
        partner_sel="protein_sel and resi 116 and name CA",
        partner_label="Cα116",
        trace_xvg=os.path.join(BINDERS, "pair_3088_binder", RUN_SUBDIR, "gate_latch116_timeseries.xvg"),
        ylim=(0.6, 1.8),
        out_name="gate_latch_closure_ghost_pair3088_binder_440-460ns",
    ),
    # Ca88 -> ligand 3-OH oxygen: the gate closing onto the ligand itself.
    # Ranked first of 39 binders (40-500 ns PL-only trajectories) for the
    # largest drop that stays down >= 20 ns: open ~1.2 nm from 40-115 ns, two-
    # step closure 115-119 ns to ~0.6 nm, closed (median 0.59 nm) to 500 ns.
    # Atom 2888 is the 3-OH O: of LIG's three O atoms it is the only one with
    # no other O within 3 A (2886/2887 are the carboxylate pair, 2.2 A apart).
    # Ligand atom names are anonymized ("O"), so it is selected by PDB serial.
    # xtc: gmx trjconv -f PL_only_40_500ns.xtc -b 108000 -e 128000
    # xvg: gmx distance -s medoid_PL.pdb -f <xtc above> -select "atomnr 1383 2888"
    "pair3087_lig3oh": dict(
        seq="pair_3087_binder",
        xtc="gate_closure_pair3087_108_128ns_PL.xtc",
        t0_ns=108.0,
        n_raw_frames=534,
        partner_sel="ligand_sel and id 2888",
        partner_label="LCA 3-OH O",
        trace_xvg=os.path.join(SCRATCH, "gate_closure_pair3087_108_128ns_ca88_lig3OH.xvg"),
        ylim=(0.4, 1.6),
        out_name="gate_closure_ghost_pair3087_binder_108-128ns_ca88_lig3OH",
    ),
}
DEFAULT_MOVIE = "pair3087_lig3oh"


def get_movie() -> dict:
    """Returns the selected movie's settings, with derived paths filled in.

    Returns:
        The MOVIES entry for $GATE_MOVIE (or DEFAULT_MOVIE), plus "key",
        "pdb", "xtc" (absolute), "t_end_ns", and "open_ns" (time of the
        first loaded state, i.e. the ghost's time).

    Raises:
        KeyError: $GATE_MOVIE is not a MOVIES key.
    """
    key = os.environ.get("GATE_MOVIE", DEFAULT_MOVIE)
    if key not in MOVIES:
        raise KeyError(f"GATE_MOVIE={key!r}; choose from {sorted(MOVIES)}")
    m = dict(MOVIES[key], key=key)
    stride = int(os.environ.get("GATE_LATCH_STRIDE", "2"))
    m["pdb"] = os.path.join(BINDERS, m["seq"], RUN_SUBDIR, "medoid_PL.pdb")
    m["xtc"] = os.path.join(SCRATCH, m["xtc"])
    m["t_end_ns"] = m["t0_ns"] + (m["n_raw_frames"] - 1) * TIME_SPACING_NS
    m["open_ns"] = m["t0_ns"] + (stride - 1) * TIME_SPACING_NS
    return m
