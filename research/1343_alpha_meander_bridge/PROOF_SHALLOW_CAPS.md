# HC: shallow connected caps for a two-color noose

## Statement

Let a disc contain portions of two simple closed curves (red and blue), with transverse intersection points in its interior. Suppose their intersections with the boundary are transverse and finite, and each color has at least two boundary ports. For each color separately, the interior components determine a noncrossing perfect matching of its boundary ports. Write `r` and `b` for the respective port counts and `w=r+b`.

There are noncrossing matchings of the red and blue ports in the complementary disc such that each monochromatic curve is completed to **one** simple closed curve and their complementary-disc arcs have at most `2r` crossings. Consequently this topological closure creates at most `2w` new red-blue crossings, with no changes to the old interior crossings. This is a linear-size closure statement. It does **not** by itself bound the inverse multiplicity of the construction: the separating noose must still be encoded or recovered from the output.

## A shallow connected complement for one color

Label the `2m` ports of a single color in boundary order `1,…,2m`. Let `P` be their matching by interior arcs. Begin with the adjacent exterior matching

`Q_0={(1,2),(3,4),…,(2m−1,2m)}`.

Process the `m−1` gaps between consecutive adjacent-pair blocks, from left to right. At gap `k`, compare the loops of `P∪Q` containing ports `2k` and `2k+1`. If they are the same, leave `Q` unchanged. If they are distinct, the two exterior `Q` edges at those ports have the form `(a,2k)` and `(2k+1,d)` with `a≤2k<2k+1≤d`. Replace them by `(a,d)` and `(2k,2k+1)`.

Inductively, `Q` consists of disjoint consecutive runs of original blocks. Within each run its two extreme ports are joined by one outer arc, and the remaining ports are joined by consecutive short arcs across the internal gaps. A one-block run has its original adjacent arc. The stated swap joins the two runs on either side of the current gap, preserving noncrossing. Each boundary point lies beneath at most two arcs of `Q`: the outer arc of its run and possibly one short arc.

Each swap acts on edges of distinct alternating cycles of `P∪Q`, and therefore joins those two cycles into one. It never splits a cycle. The two ports of each original block begin in the same cycle and remain so. Once gap `k` has been processed, the cycles containing original blocks `k` and `k+1` agree; subsequent swaps preserve this. After the last gap, all blocks belong to one cycle. Thus `P∪Q` is one cycle and `Q` has nesting depth at most two. The procedure works for any interior matching `P`; its noncrossing property is needed for the geometric interior drawing, not for the cycle argument.

## Two colors and crossing count

Apply the preceding construction to the red and blue port orders independently. Draw each exterior matching as disjoint arcs in the complementary disc, with arcs of different colors transverse and in general position. Two such arcs with endpoints on the boundary must cross only when their endpoints alternate; choose chord-style representatives so each alternating pair crosses exactly once and every other pair is disjoint.

If a red exterior arc crosses a blue exterior arc, exactly one of its red endpoints lies between the blue arc's endpoints along the chosen linear boundary order. Charge the crossing to that red endpoint and the blue arc. At each red endpoint, at most two blue arcs cover it because the blue matching has nesting depth at most two even after interleaving the red ports. Hence the number of exterior crossings is at most `2r≤2w`.

For either color, the alternating union of its interior and exterior matchings is one cycle. The resulting red and blue curves are therefore simple and connected. Smoothing their boundary joins and perturbing intersections generically realizes a meander on the sphere with the old interior crossings plus at most `2w` new crossings. The cyclic order of all ports, and the exact cap matching, are recoverable from the source and deterministic rule, but locating the old disc in the unmarked output remains the separate inverse problem.

## Relevance to a balanced sphere-cut noose

The single balanced noose in `1340_alpha_campaign/PROOF_SINGLE_NOOSE.md` has `w=O(√n)` edge crossings and balances the old intersection vertices. This lemma shows that both disc pieces can be separately completed to meanders with an additive `O(√n)` number of intersections. The old `n^{O(w)}` localization cost is **not** removed until a bounded-alphabet record identifies the separating collar or noose from an output root; the proposed primal ring construction is a candidate.
