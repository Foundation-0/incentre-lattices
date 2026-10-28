# 2 · Chemistry, crystals and materials: what the "dot in each pyramid" rule does

Code: `research/chemistry/` (install with `pip install -r research/chemistry/requirements.txt`: numpy, scipy, ase,
pandas; certifi only for `fetch_borides.py`; run each script from inside `research/chemistry/`). Scripts:
`polytopes.py` (rule on every polytope, linear-programming insphere, brute-force contact check),
`isohedral.py` (bipyramids, trapezohedra, Catalan solids), `borides_ideal.py` (tile space with the rule),
`fetch_borides.py` + `borides_real.py` (45 COD crystal structures), `bnsl.py` (binary hard-sphere packings),
`perovskite_mu.py` (Bartel et al. 2019 data). Computed numbers below are printed by one of these or by the
`results/` folder cited beside them, or derived by hand from that output where stated. The source articles of the
IUCr-provided COD files are listed under "Data sources" at the end.

Labels: **PROVEN** (argument given), **VERIFIED NUMERICALLY**, **CONJECTURE**, **KNOWN** (with citation).

---

## Finding 1. Tiled, the drawing is Pauling's 1934 calcium hexaboride (CaB6) structure, exactly

Put the cube's corners on a lattice and apply the rule in every cube. The dots are the boron atoms of CaB6,
and the corners are the calcium atoms.

- Boron sits at (x, ½, ½) with Ca at the origin. The rule gives x = ½ − (2−√2)/2 = **(√2−1)/2 = 0.20711**.
  This is the unique x at which the B–B distance inside an octahedron equals the B–B distance between
  octahedra (both 0.41421 a). "Neighbouring balls touch" means "all five B–B bonds equal". (PROVEN, VERIFIED NUMERICALLY)
- Lattice constant a = (1+√2)·d(B–B). Real MB6: mean a/d̄(B–B) = 2.399 vs 2.414 ideal. (VERIFIED NUMERICALLY on COD data)
- Each metal sits in a 24-boron cage that is a **truncated cube**, at M–B = 1.7788 d(B–B). The cage's
  hole/ball ratio is 2.5576, which is the radius-ratio value for the truncated cube. Real: 1.765. (VERIFIED NUMERICALLY)
- **History (KNOWN).** Pauling & Weinbaum, *Z. Kristallogr.* 87, 181 (1934), COD 1011310, give B at
  z = 0.293(1) from the octahedron centre, i.e. x = 0.207. In a 1955 letter Pauling wrote that the B–B
  distances between octahedra are "identical, to within 0.01 Å" with those inside an octahedron
  (Oregon State Pauling archive). So the touching-balls picture is Pauling's own CaB6 model. The vertex net of
  the truncated cubic honeycomb is the RCSR net **cab** ("CaB6 boron net", = pcu-a, a tiling by octahedra
  and truncated cubes; OEIS A299266 gives coordination sequence 1, 5, 9, 22, …). Known too.
- **What real crystals do (VERIFIED NUMERICALLY, 45 COD entries).** 14 of the 31 MB6 entries carry x = 0.207
  exactly (Ba, Ca, Ce, Er, Gd, Nd, Pr, Pu, Sr, Th, Y, Yb; these are records from 1934 to 1960 and handbook
  re-entries; COD 1009052, 2311004 and 5910016 were checked by hand). These entries are the equal-bond
  assumption copied forward, not measurements. The 17 entries with x ≠ 0.207 (16 refinements for Ce, Eu, La, Sm
  and Th, dated 1981–2010 except the undated COD 1000055, plus the 1932 von Stackelberg–Neumann CaB6 entry,
  COD 1010302, x = 0.20(2)) give x = 0.197–0.204. Over those 17 entries inter/intra B–B = **0.947 on average**
  (range 0.920–0.977; also 0.947 over the 16 without COD 1010302); `borides_real.py` itself prints the mean
  over all 31 MB6 entries, 0.9704 (range 0.920–0.999), because the 14 copied-forward entries sit at 0.999.
  So in modern refinements the bond between octahedra is ~5% shorter than the "touching" ideal. This is
  consistent with the known bonding picture:
  each B has one 2-centre exo bond and multicentre bonds inside the cage (Longuet-Higgins & Roberts 1954).
  **The ideal is a geometric zero, and the deviation from it measures the bonding asymmetry.**
- The perovskite reading (KNOWN in outline, framing ours): Ca is at the A site and the empty octahedron centre
  is at the B site. The touching pair that straddles each face centre, the "split oxygen" pair, is the
  inter-octahedral B–B bond, whose midpoint lies exactly on the perovskite O site. CaB6 is often called
  "CsCl-like Ca + B6", and B6 units are compared with corner-linked octahedra. We found no source that says
  in words that "the exo B–B bond midpoint is the perovskite O site". This is a framing, not a discovery.
- Negative result: real perovskite O displacements are transverse (Glazer tilts), not split along the B–O–B
  axis. Breathing distortions (RNiO3, BaBiO3) shift O along the axis, but they move it as a whole and do not
  split it into a touching pair. No oxide in this check realised the split-O picture. Hexaborides do.

## Finding 2. The same rule on the fcc cell gives the dodecaboride (UB12) structure, and the metal cages come out as Archimedean radius ratios

Use the fcc Voronoi cell (rhombic dodecahedron) instead of the cube. It is face-transitive, all its walls
are mirrors, and it has 12 rhombic pyramids:
- The dots form a B12 **cuboctahedron** in each cell, and each dot touches 4 + 1 = 5 others, as in CaB6.
  B sits at 48i (½, x, x) with **x = 1/6**, and a = 3√2·d(B–B). This is the UB12 (D2f, cF52) structure.
  (VERIFIED NUMERICALLY)
- The metal sits at the 4-valent cell vertex, in a 24-B **truncated octahedron**, at M–B = √(5/2)·d =
  1.5811 d. Its hole/ball ratio is 2.1623, which is the truncated-octahedron radius ratio. Real (Sc, U, Y, Yb, Zr):
  M–B/d = 1.567; inter/intra = 0.9559 on average over the 14 entries (range 0.931–1.012), as `borides_real.py`
  prints. That is a shortening of about 4% on average, not a uniform ~5%: U (COD 1008115) has 1.012, a
  lengthening, and the other four metals lie between 0.931 and 0.967. (The earlier figure 0.949 is a mean of
  per-metal means with U left out.) The empty tetrahedral hole is a
  12-B truncated tetrahedron with ratio 1.3452. (VERIFIED NUMERICALLY)
- So the rule, tiled over a lattice whose Voronoi cell is face-transitive with mirror walls, produces the
  equal-edge augmented net (pcu → CaB6, fcu → UB12). The cage around each lattice vertex is an Archimedean
  polyhedron whose radius ratio is the cage's hole/ball ratio. The bcc cell (truncated octahedron) **fails**:
  it has two pyramid types (R = 0.500 and 0.634) and no contacts across the square/hexagon walls. Real
  bcc-based borides are not of this type. (VERIFIED NUMERICALLY)
- Using the cage ratios, MB6 cages hold M–B ≈ 3.07 Å and MB12 cages ≈ 2.79 Å (worked out by hand: the ideal
  M–B/d(B–B) of 1.7788 and 1.5811 times the real mean B–B 1.726 and 1.767 Å printed by `borides_real.py`).
  That fits the known trend that MB12 forms with smaller metals (Zr, Lu, Y, U) and MB6 with larger ones
  (Ca, Sr, Ba, La, Eu). (CONJECTURE as a predictor. We did not build a classifier with a vetted list of formers,
  so this is qualitative only, and the size trend itself is KNOWN in boride literature.)
- Novelty: UB12 = B12 cuboctahedra in a NaCl arrangement, with M in truncated B24 octahedra, is KNOWN (AFLOW
  prototype A12B_cF52_225_h_b; arXiv 2010.16239). That the equal-bond x is 1/6, and that both boride families
  come from one "pyramid insphere in the Voronoi cell" rule: **we found no prior reference** after
  searching for augmented nets fcu-a/pcu-a with borides and for Voronoi/insphere constructions. RCSR's raw data file
  (http://rcsr.net/data/3dall.txt) is reachable. Its `ubt` entry ("UB12", Fm-3m, a = 3√2 for unit edge, vertex at
  (1/6, 1/6, 0), 5-coordinated) and its `cab` entry (Pm-3m, a = 1 + √2, vertex at (0.2929, 0, 0)) list the
  equal-edge geometries derived here, so x = 1/6 is stated there implicitly; the one-rule link between the two
  families is still not found in print.

## Finding 3. CaB6-type binary sphere packing peaks exactly at the dot geometry: γ = √2−1, φ = (π/6)(30√2−41) = 0.746865

Treat Ca as a large sphere and B as small spheres, then maximise the packing fraction over x and a for each
size ratio γ = r_S/r_L (`bnsl.py`):
- The maximum of the CaB6 type is at **γ = √2−1 = 0.4142, x = 0.2071, φ = (π/6)(30√2−41) = 0.746865**
  (`bnsl.py` prints φ = 0.7469 for the dot-ideal configuration; its 0.005 γ grid misses √2−1, and its grid
  maximum is 0.7439 at γ = 0.415; exact value and proof in `results/cab6-packing`). At that point the large spheres
  form a touching simple-cubic array, and each small sphere touches five other small spheres. Ca and B do **not**
  touch: d(Ca–B) = 0.7368 a against r_L + r_S = 0.7071 a, a gap of 0.0297 a, so Ca is held by Ca–Ca contact
  (`results/cab6-packing`). φ = (4π/3)[1/8 + 6((√2−1)/2)³] = (π/6)(30√2 − 41). (PROVEN: upper bound and exact
  all-pairs attainment, `results/cab6-packing/PROOF.md` and `ATTAINMENT.md`, internal review only; `bnsl.py`
  corroborates numerically)
- This is slightly above fcc (0.7405) but below NaCl at the same γ (0.793; both values from
  `results/cab6-packing`), so hard spheres alone do not favour it. That is consistent with CaB6 nanocrystal
  superlattices being described as low-density, kinetically or ligand-stabilised phases (KNOWN: Shevchenko
  et al., *Nature* 439, 55 (2006) reports CaB6-type BNSLs; ligand-exchange study, ACS figshare 2105017). One
  supplementary-table example gives core diameters Au 5 nm / γ-Fe2O3 13.4 nm (γ_core ≈ 0.37). Ligand shells
  push the effective γ up toward ~0.4, but we could not pin down effective sizes.
  (CONJECTURE that observed CaB6 BNSLs cluster near γ_eff ≈ 0.41)
- The UB12 type's optimum is **not** at its dot point (γ = 0.4625 gives φ = 0.606; the best is a degenerate
  small-γ limit), so the coincidence is special to the cube. (VERIFIED NUMERICALLY)
- Novelty (PARTIALLY KNOWN): the curve φ_max(γ) for this lattice was plotted by Horst and Travesset,
  *J. Chem. Phys.* 144, 014502 (2016), Fig. 3, so the peak can be read off it; Travesset, *ACS Nano* 11,
  5375 (2017) may state its location (not read). No explicit statement of the closed form, or of its identity
  with the √2−1 octahedral limit, was found (`results/cab6-packing/PROOF.md`, Prior art).

## Finding 4. The general rule reproduces exactly the radius-ratio limits of equal-edged, vertex-transitive coordination polyhedra, and nothing else

**Theorem (PROVEN under H1 or H2; clause (iv) under H2 only).** Take a polytope P and put the dot of each
facet-pyramid at the point equidistant from its walls, with radius R = that distance. H1: every wall (the plane
through the centre and a ridge) is a mirror of P. H2: P has an insphere centred at the centre and the dihedral
angle is the same at every ridge. Under H1 or H2:
(i) all dots are equivalent, at distance t;
(ii) a ball and its mirror image touch at the wall, so every pair of ridge-adjacent dots is at 2R: the
dots form an equal-edged, vertex-transitive polytope Q combinatorially dual to P;
(iii) hole/ball = (t−R)/R = t/(e/2) − 1, which is the classical radius-ratio limit of Q;
(iv) under H2 only, hole/ball = 1/sin θ − 1 with 2θ the angle between adjacent facet normals. Mirror walls alone do
not give (iv): the triangular bipyramid satisfies H1 at every height and its true ratio is 0.5275 throughout, but
1/sin θ − 1 gives 1.236 and 0.291 at height 1 (the two kinds of adjacent normals) and 0.633 and 0.414 at height 0.5
(`results/radius-ratio`). **The base facet never enters the ratio.** It only sets the scale, and it gives the contacts across faces when tiled.
In kaleidoscope language the dot is a Wythoff point, and the walls being mirrors is exactly what makes
the balls touch (H1).

Results (VERIFIED NUMERICALLY, LP insphere plus brute-force contacts):

| P (cut into pyramids) | balls | all touch | hole/ball | coordination polyhedron Q, classical limit |
|---|---|---|---|---|
| triangle (2D) | 3 | yes | 0.1547 | triangle CN3, 0.155 ✓ |
| square (2D) | 4 | yes | 0.4142 | square planar CN4, 0.414 ✓ |
| pentagon / hexagon (2D) | 5 / 6 | yes | 0.7013 / 1.0000 | 2D kissing at 6 |
| tetrahedron | 4 | 6/6 | 0.2247 | tetrahedron CN4 ✓ |
| cube | 6 | 12/12 | 0.4142 | octahedron CN6 ✓ |
| octahedron | 8 | 12/12 | 0.7321 | cube CN8 ✓ |
| triangular bipyramid, **any height** | 6 | 9/9 | 0.5275 | trigonal prism CN6 ✓ (self-tuning) |
| n-bipyramid, any height (n=4,5,6) | 2n | all | 0.7321, 0.9734, 1.2361 | uniform n-prism (hex prism CN12) ✓ |
| dodecahedron | 12 | 30/30 | 0.9021 | icosahedron CN12 ✓ |
| rhombic dodecahedron | 12 | 24/24 | 1.0000 | cuboctahedron CN12 ✓ (3D kissing) |
| icosahedron | 20 | 30/30 | 1.8025 | dodecahedron CN20 |
| Catalan duals | 12–30 | all | 1.3452, 2.1623, 2.5576, 1.7979, 2.2361 | truncated tetra., truncated octa., truncated cube, rhombicuboct., icosidodec. ✓ |
| tetragonal trapezohedron | 8 | **only at one shape** (16/16), otherwise 8/16 | 0.6453 | square antiprism CN8 ✓ only by tuning |
| 5-cell / tesseract / 16-cell / 24-cell (4D) | 5/8/16/24 | all | 0.2649 / 0.4142 / 1.0000 / 1.0000 | 4D: 24-cell = 4D kissing |
| cuboctahedron, truncated octahedron, prisms (not facet-transitive) | n/a | **no** (unequal R) | n/a | fails |

Why it fails where it does (PROVEN):
- **Antiprisms** (square antiprism CN8, 0.645) are snub-type and not Wythoffian. The trapezohedron's walls are
  not mirrors, so touching is not forced. It occurs only on a one-parameter coincidence (the exact dual of the
  uniform antiprism). Off that shape, half the contacts open up (ratio 0.6489–0.6490 at stretch 0.8 and 1.3).
- **CN7, CN9 (tricapped trigonal prism 0.732), dodecahedral CN8 (0.668), capped octahedron (0.592)**: these
  anion polyhedra are not vertex-transitive. The rule's dots are always one orbit, so it cannot produce them.
- **bcc cell / non-isohedral P**: there are several pyramid types with different R and no contacts.
- **Self-tuning (VERIFIED NUMERICALLY, and a corollary of the theorem):** a bipyramid of any height gives the
  *uniform* prism, because mirror walls force all edges to 2R. The rule finds the equal-edge geometry without
  being told.

Novelty: the radius ratios themselves are KNOWN (Pauling 1929; Wikipedia's "Pauling's rules" table lists
0.155–1.00). The vertex/facet duality is trivial. The statement "the limits reachable this way are exactly the
Wythoffian (mirror-generated) equal-edge polytopes; antiprisms are reachable only by fine-tuning; CN7, CN9 and
dodecahedral CN8 are unreachable" is one we found **no prior reference** for, after searching radius ratio +
Wythoff / kaleidoscope / uniform polyhedra. It is mathematically elementary, so it may well appear in some
textbook.

## Finding 5. Perovskites: √2−1 is the near-optimal octahedral-factor cut, but the geometry gives no sharper predictor

Data: Bartel et al., *Sci. Adv.* 5, eaav0693 (2019), Table S1 (576 ABX3, 313 perovskites; GitHub CJBartel/perovskite-stability).
- μ = r_B/r_X ≥ √2−1 alone: accuracy 0.752. The best fitted single μ cut is 0.426 (0.769), within 3% of √2−1.
  Goldschmidt window 0.825<t<1.059: 0.736. Both together: **0.835**. Fitting all three cuts gives 0.847.
  Bartel's τ: **0.915**. (VERIFIED NUMERICALLY)
- By anion: oxides tolerate μ down to 0.29 (28 perovskites violate √2−1, mostly oxides and fluorides with
  small B by Shannon radii). Halides (Cl, Br, I) have no perovskite below 0.43. That matches Li et al.'s
  empirical halide bound ~0.44 (KNOWN; Li et al., *Acta Cryst. B* 64, 702 (2008)).
- PROVEN identity: in the ideal cubic perovskite (t = 1), μ = √2−1 ⟺ r_A = r_X. The B-site octahedral limit
  and the A-site cuboctahedral limit (ratio 1.000, rhombic-dodecahedron row above) coincide only at t = 1.
  This is the geometric content of Goldschmidt's √2 (A–X = √2·B–X). It is a restatement, not a new predictor.
- Conclusion: the dot geometry *derives* the 0.414 bound (KNOWN as the octahedral radius ratio) but does not
  beat τ. We do not claim more.

## Radius-ratio failures (George et al. 2020)

Verified: J. George, D. Waroquiers, …, G. Hautier, "The Limited Predictive Power of the Pauling Rules",
*Angew. Chem. Int. Ed.* 59, 7569 (2020), doi:10.1002/anie.202000829 (PMC7217010). Pauling's first rule holds
for ~66% of ~5000 oxide environments, and only 13% satisfy rules 2–5 together. The kaleidoscope view explains
one structural reason for failure (Finding 4): limits exist only for equal-edged vertex-transitive polyhedra,
while real environments are distorted and include the non-uniform CN7–CN9 shapes, which have no "touching"
limit of this kind. It does not explain the softness and size failures, which dominate. (CONJECTURE that the
unreachable shapes are over-represented among violations; not tested.)

## What was not done

- RCSR's raw data file (http://rcsr.net/data/3dall.txt) was reached late: it confirms the `cab` and `ubt` geometries, but the net symbols elsewhere in these notes still come from OEIS A299266 and the Yaghi group's MOK-net notes.
- No DFT. The shortening of the inter-octahedron bond in modern refinements (about 5% for MB6, about 4% on average for MB12, with U lengthened) is from the literature data only.
- Glazer tilts as rigid-unit modes were not explored.

## Best open question (PhD-sized)

**Is "Voronoi-cell Wythoff points → equal-bond cluster frameworks" a generative rule for boron-rich (and other
cluster) solids?** Cube → CaB6 and rhombic dodecahedron → UB12 both work, with the metal cage giving an
Archimedean radius ratio and the real MB6 structures about 5% off ideal and the MB12 ones about 4% off on average (U goes the other way). Next steps:
(a) enumerate all face-transitive, mirror-walled cells of space groups (and their sub-tilings) and list the
resulting augmented nets as candidate MBn frameworks;
(b) explain the dodecahedron case: icosahedral B12 (5+1 bonding, α-boron) is the frustrated version that cannot
tile, a Frank–Kasper-like story;
(c) test whether the measured inter/intra ratio (0.92–1.01) tracks the metal's electron count, cage filling
or temperature with DFT and modern refinements. That would turn the geometric zero into a quantitative bond
descriptor.

## Sources
- Pauling & Weinbaum 1934, COD 1011310; Pauling 1955 letter: https://scarc.library.oregonstate.edu/coll/pauling/calendar/1955/07/11.html
- AFLOW CaB6 prototype: https://aflowlib.duke.edu/prototype-encyclopedia/A6B_cP7_221_e_b-001 ; UB12: https://aflowlib.duke.edu/prototype-encyclopedia/A12B_cF52_225_h_b-001
- Dodecaboride review: https://arxiv.org/pdf/2010.16239
- cab net tiling: https://oeis.org/A299266 ; https://yaghi.cchem.berkeley.edu/research-news/MOK-nets2.pdf
- George et al. 2020: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7217010/
- Shevchenko et al. 2006: https://ideas.repec.org/a/nat/nature/v439y2006i7072d10.1038_nature04414.html
- Horst & Travesset 2016, *J. Chem. Phys.* 144, 014502, Fig. 3: https://doi.org/10.1063/1.4939238 ; Travesset 2017,
  *ACS Nano* 11, 5375 (not read): https://doi.org/10.1021/acsnano.7b02219
- Bartel et al. 2019 data: https://github.com/CJBartel/perovskite-stability
- Hexaboride B–B bond lengths: https://www.jmst.org/EN/10.1016/j.jmst.2017.01.022

## Data sources

25 of the 45 COD files in `research/chemistry/cod/` (all the 2xxxxxx entries: 15 MB6 and 10 MB12) state that
their data were provided by IUCr Journals and may be used "so long as proper attribution is given to the journal
article from which the data were obtained" (see `NOTICE`). All 25 enter the statistics of Findings 1 and 2. Their
articles, as each file names them (`_publ_author_name`, `_journal_name_full`, `_journal_volume`, pages,
`_journal_year`, `_journal_paper_doi`):

| COD | compound | authors | article | DOI |
|---|---|---|---|---|
| 2101693 | CeB6 | K. Tanaka, Y. Kato, Y. Ōnuki | *Acta Cryst. B* 53, 143–152 (1997) | 10.1107/S0108768196009019 |
| 2102953 | CeB6 | R. Makita, K. Tanaka, Y. Ōnuki, H. Tatewaki | *Acta Cryst. B* 63, 683–692 (2007) | 10.1107/S0108768107035021 |
| 2103133, 2103136 | CeB6 | K. Tanaka, Y. Ōnuki | *Acta Cryst. B* 58, 423–436 (2002) | 10.1107/S010876810102167X |
| 2104215, 2104216 | CeB6 | R. Makita, K. Tanaka, Y. Ōnuki | *Acta Cryst. B* 64, 534–549 (2008) | 10.1107/S0108768108026542 |
| 2104736, 2104747 | LaB6 (a phase in a powder study of pyrrhotite) | A. D. Elliot | *Acta Cryst. B* 66, 271–279 (2010) | 10.1107/S0108768110011845 |
| 2104762–2104765 | SmB6 | S. Funahashi, K. Tanaka, F. Iga | *Acta Cryst. B* 66, 292–306 (2010) | 10.1107/S0108768110009250 |
| 2108722–2108730 | YbB12 | N. Bolotina, O. Khrykina, A. Azarevich, S. Gavrilkin, N. Sluchanko | *Acta Cryst. B* 76, no. 6 (2020; the files give no page) | 10.1107/S2052520620013566 |
| 2310102 | YB12 | R. Barrett, V. I. Matkovic, R. F. Giese Jr., J. Economy | *Acta Cryst.* 19, 1056–1058 (1965) | 10.1107/S0365110X65004954 |
| 2310112 | PuB6 | W. I. Stuart, B. J. McDonald | *Acta Cryst.* 13, 447–448 (1960) | 10.1107/S0365110X60001059 |
| 2311004, 2311005 | NdB6, YbB6 | P. Blum, F. Bertaut | *Acta Cryst.* 7, 81–86 (1954) | none in these files; COD 1008115 and 1009032 give 10.1107/S0365110X54000151 for the same article |

The other 20 files carry the COD public-domain header and are not listed here.
