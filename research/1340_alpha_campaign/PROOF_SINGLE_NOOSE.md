# A single balanced disc separator for meanders

Status: proved, independently reviewed in `1342_alpha_tree_entropy/REFEREE_SINGLE_NOOSE.md`. This does not yet improve the counting bound.

## Statement

For every sufficiently large order N meander, viewed as two transverse Jordan curves on the sphere with 2N intersections, there is one Jordan curve beta disjoint from those intersections, transverse to both curves, with at most 4 sqrt(27N) intersections in total. Each of its two complementary discs contains between N/2 and 3N/2 original intersections. More precisely the counts lie between 2N/3 - sqrt(27N) and 4N/3 + sqrt(27N).

Unlike the separator used in 628, this separator has exactly two disc sides, with no holes. It need not be a river gap, and it need not intersect the river in only two points.

## Published inputs

Dorn, Penninkx, Bodlaender and Fomin, *Efficient Exact Algorithms on Planar Graphs: Exploiting Sphere Cut Decompositions*, technical report UU-CS-2006-006, Proposition 1 and Theorem 1, pp. 4-5:

https://webdoc.sub.gwdg.de/ebook/serien/ah/UU-CS/2006-006.pdf

These give branchwidth at most sqrt(4.5v) for a planar v-vertex graph and a sphere-cut decomposition of that width for a connected plane graph without degree-one vertices. Each decomposition edge is realized by a noose meeting only its middle-set vertices, with the two edge sets on opposite disc sides. The theorem is attributed there to Seymour--Thomas and Gu--Tamaki. No algorithmic running-time assertion is needed here.

## Application and proof

The union of river and road is a plane multigraph G with 2N vertices, 4N edges and degree four. For N>=2 it has no loops (parallel edges are allowed). Subdivide each edge once to obtain a simple connected plane graph H with 6N vertices and minimum degree two. This avoids any convention about applying the cited results to parallel edges. Put k=sqrt(27N). The published inputs supply a sphere-cut decomposition of H of width at most k.

For each original vertex v choose one incident half-edge h(v) in H. These chosen half-edges are distinct because every half-edge has exactly one original endpoint. Give the corresponding decomposition leaves weight one and all other leaves weight zero. Their total weight is 2N.

A finite ternary tree with leaf weights zero or one and total weight W>=3 has an edge whose two components each have weight in [W/3,2W/3]. To see this, select a weighted centroid, so every component left by deleting it has weight at most W/2. A centroid cannot be a leaf when W>2. Its three branch weights sum to W, so one lies in [W/3,W/2]; its incident edge works.

Apply this with W=2N, and let O be the noose of the selected edge. At most k vertices of H lie on O. For an original vertex not on O all its incident half-edges lie on the same side as that vertex. Consequently the number of original vertices strictly on either side differs from the assigned leaf weight on that side by at most the number of original vertices on O, hence by at most k. The strict-side counts are therefore between 2N/3-k and 4N/3. This also follows directly from the defining containment of each edge subgraph in its closed disc side.

Perturb O in pairwise disjoint small neighborhoods of its boundary vertices. In each neighborhood replace the passage through the vertex by a simple detour avoiding it. The detour crosses each incident edge at most once, hence at most four edges at an original vertex and at most two at a subdivision vertex. These replacements remain disjoint and preserve the Jordan property. Outside these neighborhoods the curve is unchanged and disjoint from H. After removing the auxiliary subdivisions we obtain beta, transverse to the original curves, with at most 4k crossings. Original vertices off O have not changed sides; the at most k original boundary vertices have each been assigned to one side. Thus each disc has between 2N/3-k and 4N/3+k original intersections. For k<=N/6 this implies the asserted [N/2,3N/2] balance. Both colors meet beta in a positive even number, since each original Jordan curve visits original vertices in both discs.

## What this removes and what it does not

This replaces the potentially many-holed region in 628 by two discs. In a subsequent counting proof, there are only two interior sizes to sum over, costing a polynomial in N rather than N raised to the number of boundary components. All dangling road and river endpoints occur on one cyclic boundary on each side.

It does not supply a completion injection. A disc side contains several river strands, and their attachment locations are not a contiguous suffix in the completed meander. The remaining target is a reversible completion with O(w) added intersections and exp(O(w)) records, including recognition of all boundary ports. Standard sphere-cut dynamic programming counts boundary connectivity states; it does not count the possible interiors by M_(m+O(w)). That distinction is essential.

Ledger search: `python3 ledger/find.py text=noose` returned no findings on 2026-10-25. This is an application of published separator machinery, not a claim of a new general planar separator theorem.
