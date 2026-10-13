#!/usr/bin/env python3
"""Assemble certificates/selected-results.json, certificates/SHA256SUMS and commitment/*.

Standard library only.
  python tools/build_certificates.py --write    regenerate the three outputs from results/*/certificate.json
  python tools/build_certificates.py --check    fail if any output differs from what would be generated,
                                                or if a hash in SHA256SUMS or the commitment does not match
selected-results.json has one record per result: statement, value, evidence label, verifier command,
input hashes (copied from the certificate), the certificate hash and the prior-art status.
Evidence labels: PROVEN (written proof, in PROOF.md), VERIFIED (finite computation, not a proof),
KNOWN (published; cited), NEW-COMPUTATION (our own, not independently reproduced beyond the listed checkers).
"""
import json
import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import certlib  # noqa: E402

ROOT = certlib.ROOT
RESULTS = [
    {
        "id": "rigidity-nd",
        "statement": "The contact framework of the truncated n-cubic honeycomb (2n balls per cell, 2n^2-n bars of one length, dots at the pyramid incentres) has, at a generic wavevector, exactly n(n-1)/2 zero modes (the cluster rotations) and n(n-3)/2 self-stresses per cell, its cross-polytope cell is infinitesimally rigid for n >= 3, and with unit springs and full internal relaxation C11 = 4(n-1)^2/(n^2 2^n), C12 = 4/(n^2 2^n), all shear moduli vanish, nu_100 = 1/(n^2-n-1).",
        "value": {"nu_100": {"n=3": "1/5", "n=4": "1/11", "n=5": "1/19", "n=6": "1/29", "n=7": "1/41", "n=8": "1/55"},
                  "generic_zero_modes_per_cell": "n(n-1)/2", "generic_self_stresses_per_cell": "n(n-3)/2"},
        "evidence_label": "PROVEN for all n >= 3 (written proof, results/rigidity-nd/PROOF.md); VERIFIED (finite check) for n = 2..8",
        "verify_command": "make verify-rigidity",
        "prior_art": "3D case KNOWN in substance (Goodwin 2006, rigid unit modes of linearly bridged frameworks; Power 2014). The n-dimensional closed forms were NOT FOUND in the literature (absence of evidence only).",
    },
    {
        "id": "radius-ratio",
        "statement": "For a convex polytope P (origin inside) whose pyramids over the facets each have an inscribed ball (the dots), if every wall is a mirror of P (H1), or P has an insphere centred at the origin and all dihedral angles are equal (H2), then all dots have equal radius and equal distance D from the origin, adjacent balls touch, and hole/ball = (D-r)/r equals the classical radius-ratio limit of the dot polytope.",
        "value": {"example_hole_over_ball": {"cube->octahedron": "sqrt(2)-1", "octahedron->cube": "sqrt(3)-1", "rhombic dodecahedron->cuboctahedron": 1},
                  "counterexample": "tetragonal trapezohedra polar to square antiprisms of height 0.7x and 1.3x the uniform one (facet-transitive, insphere, but contacts open)"},
        "evidence_label": "PROVEN under H1 or H2 (written proof); VERIFIED (finite check) on 20 polytopes plus 4 counterexample controls",
        "verify_command": "make verify-radius-ratio",
        "prior_art": "NOT FOUND as a statement; the ingredients (Wythoff construction, classical radius-ratio rule, polar duality) are KNOWN and the proof is elementary.",
    },
    {
        "id": "cab6-packing",
        "statement": "The maximal hard-sphere packing fraction of the CaB6 structure type at size ratio gamma = rS/rL is (pi/6)(1+6 gamma^3) for gamma <= sqrt2-1 and (pi/6)(sqrt2-1)^3 (6+gamma^-3) for gamma >= sqrt2-1; the global maximum (pi/6)(30 sqrt2 - 41) = 0.746865 is attained at gamma = sqrt2-1, where Ca is held by Ca-Ca contact and Ca-B do not touch; the packing beats fcc exactly for 0.41023 < gamma < 0.41591.",
        "value": {"phi_max": "(pi/6)(30*sqrt(2) - 41)", "phi_max_decimal": 0.746864891262, "window": [0.410227, 0.415911]},
        "evidence_label": "PROVEN (elementary upper-bound argument, PROOF.md); attainment and all-pairs non-overlap VERIFIED (finite check, 5x5x5 block, 15 size ratios)",
        "verify_command": "make verify-cab6",
        "prior_art": "PARTIALLY KNOWN: the curve is plotted for this and 23 other lattices by Horst and Travesset (2016); the closed form was not found in text. CaB6 is absent from the densest binary packings of Hopkins, Stillinger and Torquato (2012).",
    },
    {
        "id": "fold-2d",
        "statement": "The 4.8.8 hinged-squares linkage read as a ball framework: along the symmetric fold the period is L = 2R + 2c cos(theta), the first new contact is at theta = pi/4 for the incentre c = 2-sqrt2, where L = 4R, the linear strain is -(sqrt2-1)^2 and the density rises from pi(sqrt2-1)^2 = 0.5390 to pi/4.",
        "value": {"density_start": 0.5390120845, "density_end": "pi/4", "lock_angle_deg": 45},
        "evidence_label": "KNOWN linkage; corollary numbers PROVEN (exact, sympy); uniqueness of the incentre among c VERIFIED on a scan only",
        "verify_command": "make verify-fold2d",
        "prior_art": "KNOWN: Flores (Hinged tilings) crediting Stuart (1963); rotating squares, Grima and Evans (2000). The ball-packing numbers are trivial corollaries.",
    },
    {
        "id": "fold-3d-cubic",
        "statement": "In the 3D framework, rotating the octahedra alternately by +-theta about [111] with the lattice kept cubic reaches a=8c/3 at theta = 60 degrees with six contacts per ball and density 27 pi sqrt2/256 = 0.4686; this is not a jam, because trigonal lattices keep every bar and reach density 0.6477 on a grid scan.",
        "value": {"density_at_60_degrees": 0.468585, "trigonal_grid_best_density": 0.6477},
        "evidence_label": "VERIFIED (finite check; exact algebra at 60 degrees, grid scan elsewhere); not a proof of any path",
        "verify_command": "make verify-fold3d",
        "prior_art": "PARTIALLY KNOWN: Prussian-blue topology and a-a-a- tilts known; the 60-degree end net matches RCSR sms in its first 10 coordination-sequence terms (not an isomorphism proof).",
    },
    {
        "id": "fold-close-packing",
        "statement": "In the 2x2x2 supercell (48 balls, 120 bars of length 2R, free lattice) there is a sampled path of 793 states from the open cab state, keeping every bar length and with no overlap, to an hcp packing of density pi/sqrt18; two further rigid layer slides (400 states) reach pure fcc. Both paths are accepted by two independent checkers.",
        "value": {"hcp_states": 793, "fcc_states": 1193, "end_density": "pi/sqrt(18)", "end_contacts_per_ball": 12,
                  "fcc_end_cna": {"(4,2,1)": 576}},
        "evidence_label": "VERIFIED (finite check): continuity certified only by sampling at step <= 0.01 with midpoint / 10-substate interpolation; hcp path checked by two checkers; the fcc extension is NEW-COMPUTATION, checked by the same two checkers",
        "verify_command": "make verify-fold-cp (and make verify-fold-tests for the checker's own tests)",
        "prior_art": "NOT FOUND as a statement. Closest: ReO3 to VF3 a-a-a- tilt towards hcp anions (O'Keeffe and Hyde 1977; Glazer 1972); sequential-LP packing methods (Torquato and Jiao 2010). See autoresearch-fold/PRIOR-ART-FOLD.md.",
    },
    {
        "id": "fold-end-exact",
        "statement": "The m = 2 hcp and fcc closed forms have all 120 bars exactly 2R, no overlap, 12 contacts per ball and density pi/sqrt18, with the hcp and fcc common-neighbour signatures. The final states of the stored paths agree with these forms to 1e-7 after a rigid motion.",
        "value": {"balls": 48, "bars": 120, "end_density": "pi/sqrt(18)", "end_contacts_per_ball": 12,
                  "hcp_cna": {"(4,2,1)": 288, "(4,2,2)": 288}, "fcc_cna": {"(4,2,1)": 576}},
        "evidence_label": "PROVEN (exact integer arithmetic) for the closed forms; VERIFIED (finite check, 1e-7) for agreement with stored endpoints; no proof of continuous paths",
        "verify_command": "make verify-fold-exact",
        "prior_art": "Exact forms of the endpoints from results/fold-close-packing; no separate prior-art search. See autoresearch-fold/PRIOR-ART-FOLD.md. One internal implementation, no outside review.",
    },
    {
        "id": "fold-m1",
        "statement": "In the m = 1 cell (6 balls, 15 bars, free lattice), a 333-state sampled path reaches hcp and a 268-state sampled path reaches fcc from the open cab state. Both pass the locked and independent checkers. The fcc density rises at every sample; both closed-form endpoints have exact bars, no overlap, 12 contacts per ball and density pi/sqrt18.",
        "value": {"balls": 6, "bars": 15, "hcp_states": 333, "fcc_states": 268,
                  "fcc_density_increasing_at_all_sampled_steps": 267, "end_density": "pi/sqrt(18)",
                  "end_contacts_per_ball": 12, "hcp_cna": {"(4,2,1)": 36, "(4,2,2)": 36},
                  "fcc_cna": {"(4,2,1)": 72}},
        "evidence_label": "VERIFIED (finite check, sampled continuity), NEW-COMPUTATION; monotonicity at the fcc samples only; endpoints PROVEN (integer arithmetic); interpolation is not exactly bar-preserving",
        "verify_command": "make verify-fold-m1",
        "prior_art": "No separate search for m = 1; see autoresearch-fold/PRIOR-ART-FOLD.md for the fold search. Checkers and endpoint implementation are internal, with no outside review.",
    },
]

RESULTS.append({
    "id": "fold-continuous-fcc",
    "statement": "An explicit smooth m = 1 motion preserves all 15 bars of the cab contact framework, keeps all spheres non-overlapping for every parameter value and periodic image, and increases density strictly from the open cubic state to fcc at pi/sqrt18. The construction lifts to every integer supercell size.",
    "value": {"balls": 6, "bars": 15, "end_density": "pi/sqrt(18)",
              "distance_polynomials": 3248, "nonnegative_bernstein_coefficients": 16240},
    "evidence_label": "PROVEN (written explicit construction, analytic image bound and monotonicity, exact rational Bernstein certificate); internal AI review only, independent human review pending",
    "verify_command": "make verify-fold-continuous",
    "prior_art": "Bridged-unit flexibility and cross-polytope deployments are known (Goodwin 2006; Borcea and Streinu 2021). The endpoint is the project's existing m = 1 fcc endpoint. Priority of this explicit equal-bar, non-overlapping motion is unresolved; no novelty claim.",
})

SUM_FILES = ["autoresearch-fold/best/path.npz", "results/fold-close-packing/path_fcc.npz",
             "results/fold-m1/path_hcp_m1.npz", "results/fold-m1/path_fcc_m1.npz"]


def cert_path(rid):
    return f"results/{rid}/certificate.json"


def validate_registry():
    """Every result certificate needs exactly one reviewed metadata record."""
    ids = [r["id"] for r in RESULTS]
    duplicates = sorted({rid for rid in ids if ids.count(rid) > 1})
    discovered = {p.parent.name for p in (ROOT / "results").glob("*/certificate.json")}
    unregistered = sorted(discovered - set(ids))
    missing = sorted(set(ids) - discovered)
    if duplicates or unregistered or missing:
        raise SystemExit(
            f"result registry mismatch: duplicate IDs={duplicates}; "
            f"unregistered certificates={unregistered}; missing certificates={missing}"
        )


def build_selected():
    validate_registry()
    recs = []
    path_inputs = set()
    for r in RESULTS:
        cp = ROOT / cert_path(r["id"])
        if not cp.exists():
            raise SystemExit(f"missing {cp}; run the verifier first")
        cert = json.loads(cp.read_text())
        path_inputs.update(p for p in cert["inputs"] if p.endswith(".npz"))
        rec = dict(r)
        rec["certificate"] = cert_path(r["id"])
        rec["certificate_sha256"] = certlib.sha256_file(cp)
        rec["input_hashes"] = cert["inputs"]
        recs.append(rec)
    unregistered_paths = sorted(path_inputs - set(SUM_FILES))
    if unregistered_paths:
        raise SystemExit(f"unregistered path inputs: {unregistered_paths}; add to SUM_FILES")
    return {"note": "Generated by tools/build_certificates.py from results/*/certificate.json. Labels: PROVEN = written proof; VERIFIED = finite computation, not a proof.",
            "results": recs}


def sums_text(paths):
    return "".join(f"{certlib.sha256_file(ROOT / p)}  {p}\n" for p in paths)


def outputs():
    sel = json.dumps(build_selected(), indent=2, sort_keys=True) + "\n"
    (ROOT / "certificates").mkdir(exist_ok=True)
    return sel


def main(mode):
    sel_text = outputs()
    sel_path = ROOT / "certificates" / "selected-results.json"
    plan = {sel_path: sel_text}
    if mode == "--write":
        sel_path.write_text(sel_text)
    elif sel_path.read_text() != sel_text:
        raise SystemExit("certificates/selected-results.json is out of date (run make certificates)")
    cert_files = [cert_path(r["id"]) for r in RESULTS]
    sums_files = cert_files + ["certificates/selected-results.json"] + SUM_FILES
    sums = sums_text(sums_files)
    # commitment: results/*/certificate.json, certificates/* (incl. SHA256SUMS), path files
    sums_path = ROOT / "certificates" / "SHA256SUMS"
    if mode == "--write":
        sums_path.write_text(sums)
    elif sums_path.read_text() != sums:
        raise SystemExit("certificates/SHA256SUMS does not match the files")
    manifest_files = cert_files + ["certificates/SHA256SUMS", "certificates/selected-results.json"] + SUM_FILES
    manifest = sums_text(sorted(manifest_files))
    root = certlib.hashlib.sha256(manifest.encode()).hexdigest()
    cm = ROOT / "commitment"
    cm.mkdir(exist_ok=True)
    if mode == "--write":
        (cm / "SHA256SUMS").write_text(manifest)
        (cm / "ROOT").write_text(root + "\n")
    else:
        if (cm / "SHA256SUMS").read_text() != manifest or (cm / "ROOT").read_text().strip() != root:
            raise SystemExit("commitment/ does not match the files (run make certificates)")
    print(("wrote" if mode == "--write" else "checked") + f" selected-results.json, certificates/SHA256SUMS, commitment (ROOT {root})")


if __name__ == "__main__":
    if len(sys.argv) != 2 or sys.argv[1] not in ("--write", "--check"):
        raise SystemExit(__doc__)
    main(sys.argv[1])
