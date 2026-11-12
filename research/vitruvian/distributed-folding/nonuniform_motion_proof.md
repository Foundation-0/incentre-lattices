# A compatible family of nonuniform folding patterns

Derived in this investigation, 2026-11-12. Reviewed independently by gpt-6-astra. Historical priority is unestablished. This is a kinematic statement for the ideal sphere-and-bar model, not a material dynamics result.

## Statement

Let P={±e1,±e2,±e3}, with sphere radius 1/sqrt(2), and let C have rows (0,-1,-1), (1,0,-1), (1,-1,0). For each direction i choose any sequence theta_i:Z→[40°,50°]. The same sequence is shared across lines parallel to that direction: theta_i depends on n_i only.

Define

    L_i(k) = 2e_i + sqrt(2)cos(theta_i(k))e_i + sin(theta_i(k))C_i,
    S_i(0)=0,  S_i(k+1)-S_i(k)=L_i(k),
    X(n)=S_1(n_1)+S_2(n_2)+S_3(n_3).

Put an identical, fixed-orientation octahedral motif X(n)+P at every integer cell n. Then all internal bars and neighboring bridges have length sqrt(2), and no pair of distinct spheres overlaps. The statement holds for arbitrary sequences in the specified intervals, including continuous time-dependent sequences staying in those intervals.

## Closure and prescribed lengths

The position formula makes loop closure automatic. Moving in direction i changes position by L_i(n_i), independently of the other coordinates. Each bridge joins X(n)+e_i to X(n+e_i)-e_i, with vector L_i(n_i)-2e_i. Its squared length is 2cos²(theta_i)+2sin²(theta_i)=2 because C_i is orthogonal to e_i and has squared length 2. The motif is fixed, so its internal edges are unchanged.

## Clearance bound

Use the uniform reference theta=45°, L0=3I+C/sqrt(2). A row changes by a circular chord of length at most 2sqrt(2)sin(2.5°)<1/8. The inequality follows from sin x≤x and pi<22/7. The final comparison can be squared and checked with integers.

For any two cells separated in material coordinates by m, telescoping positive or negative steps yields

    X(n+m)-X(n)=mL0+E,
    |E|≤|m|_1/8≤sqrt(3)|m|/8.

The symmetric part of L0 has smallest eigenvalue 3-1/sqrt(2)>9/4, so |mL0|>(9/4)|m|. Any motif difference v has |v|≤2. For |m|≥2,

    |X(n+m)-X(n)+v|
      ≥|mL0|-|E|-|v|
      >(9/4-1/4)|m|-2
      ≥2>sqrt(2).

Thus all distant pairs are safe uniformly over all choices of sequences and time. The same-cell pairs are fixed and safe. It remains to consider the 26 nonzero integer m with |m|<2 and 19 motif differences v. Six are prescribed bridges, handled exactly above.

For each of the other 488 pairs it suffices that

    |mL0+v| > sqrt(2)+|m|_1/8.

Write z=mC, u=3m+v, and q=|m|_1. Sixty-four times the difference of the squares of these two sides equals A+B sqrt(2), with

    A=64|u|²+32|z|²-128-q²,
    B=64 u·z-16q.

`nonuniform_motion.py` enumerates exactly those pairs and certifies A+B sqrt(2)>0 using integer arithmetic. This completes the non-overlap argument. The code is a finite trusted arithmetic component, not a proof-assistant formalization. The independent review checked the enumeration, exclusions, sign routine, telescoping, and far-pair bounds.

## Traveling example

For any positive integer wavelength m, choose

    theta_1(k,t)=45°+5°sin(2pi k/m-omega t),
    theta_2=theta_3=45°.

This gives a traveling modulation of local link angles. Integrating the increments gives a compatible, collision-free placement at every time. It avoids the earlier failed prescription of one common phase varying across all three directions at each position. It needs no local motif rotations within this restricted separated-coordinate construction.

The global periodic cell may vary with time, and no fixed-boundary condition is asserted. Forces, speed limits, actuation, finite connector thickness, elastic energy, and stability have not been modeled. The statement does not establish a self-propagating physical wave.

## Research significance and boundary

The result establishes finite-amplitude nonuniform configurations rather than merely an infinitesimal flexibility count or a sample animation. General rigid-unit flexibility is established prior art. The contribution here is an explicit compatible family with a uniform all-pairs clearance certificate. Whether that family or an equivalent construction is already published remains to be checked.

The useful next question is how far the certified amplitude range can be expanded, particularly as a background state approaches a transport threshold, and what controls can realize it with finite boundaries and physical connectors. A general physical claim requires those additional tests.
