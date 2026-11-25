# Fixed angle springs along the octahedral fold

A fixed angular energy gives a large but finite internal-motion response near the open fold. At the fully folded endpoint, an alternating displacement proves instability for every positive angular stiffness under the specified distributed holding loads.

The author calculation and a separate Astra review agree on the response asymptote and loading convention. This is internal project review, not external peer review or a prediction for a measured material.

- [Results and finite stability checks](RESULTS.md), with executable code and source hashes.
- [Fixed physical protocol](PROTOCOL.md) and [response derivation](DERIVATION.md).
- [Independent review](independent-review/REVIEW.md), including an independent geometry and energy implementation.
- [Exact endpoint instability proof](independent-review/ENDPOINT-PROOF.md), stronger than the original sampled result.

The author files remain unchanged at their reviewed hashes. Their original “independent review pending” status is superseded by the review linked above. The negative-curvature proof does not establish stability elsewhere. Positive spatial results remain limited to the sampled wavevectors and directions.

Run the author commands in the results note, then `python3 independent-review/independent_check.py`. The independent check uses SymPy and mpmath, which is a dependency of SymPy. It does not import the author's mechanics model.

The ideal-model question has reached a useful stopping point. A next physical calculation needs a specified hinge and a specimen held at its boundaries. Further arbitrary spring variants would not establish that realization.

This continues the constrained-response comparison in [gravity PR #11](https://github.com/Foundation-0/gravity/pull/11). The connection concerns what an effective response law omits. No equivalence between this lattice and spacetime is established.
