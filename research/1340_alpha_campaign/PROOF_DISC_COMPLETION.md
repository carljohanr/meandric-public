# Disc completion with a recognizable primal cycle

Status: proved, independently reviewed in `1342_alpha_tree_entropy/REFEREE_DISC_COMPLETION.md` and `1341_alpha_tail/REFEREE_DISC_COMPLETION.md`. The shallow-cap input in 1343 is independently reviewed in `1342_alpha_tree_entropy/REFEREE_SHALLOW_CAPS_COLLAR.md`.

## Lemma

Let a disc patch be obtained by cutting two transverse Jordan curves on the sphere along a transverse Jordan separator. Suppose both sides contain original crossings. Write v for the number of original crossings in the patch and w for its boundary ports (intersections with the separator). Each color has a positive even number of ports. There is a completion to two transverse Jordan curves with at most v+8w crossings, containing a distinguished simple primal cycle of length at most 2w. The original patch, with its cyclic port order, is recovered from the rooted completed meander with at most A^w inverse records, for an absolute constant A. No bulk edge labels are needed.

## Construction of the protective cycle

Extend the cut strands a short distance through an annular collar. Draw a simple cycle C following the boundary, crossing each extended port strand once. At a red port C is blue and goes straight through the new crossing; at a blue port C is red. Between consecutive ports with different colors, insert one switch crossing, at which C changes color using two adjacent half-edges. Place both unused half-edges in the exterior of C, so that the interior sector between the two cycle edges is empty. Between consecutive ports of the same color no switch is needed. This is a local construction in disjoint collar sectors; curving the edges freely realizes the indicated cyclic orders. At a switch the red and blue pairs are opposite in the local rotation, as required for transverse crossings. There is no angle-sum restriction because edges are curved.

If t is the number of color changes around the boundary, C has w+t<=2w vertices. Its port vertices each have one free exterior half-edge of the port color. Its switch vertices each have two free exterior half-edges, one of each color. These are distinct terminals after extending a short distance outward to another parallel boundary curve; terminals from the same switch retain their local cyclic order. Thus the total number of exterior terminals is w+2t<=3w.

Inside this outer boundary, each monochromatic component is a simple arc between two terminals. Indeed, original strands are disjoint arcs of a Jordan curve cut by a separator; they cannot contain a closed component because the original curve also visits the other disc. They pass straight through port vertices. The other monochromatic paths are the segments of C between successive switches: at each port on such a segment C has the opposite color to the port strand, so the crossing does not connect the two. Because both port colors occur, t>=2, so none of these segments closes into a monochromatic cycle. This verifies the no-closed-component hypothesis rather than assuming it from connectivity of the whole patch. For either color, its path endpoint pairing is therefore noncrossing.

## Exterior caps

Apply the shallow-cap lemma to each monochromatic pairing, with a common gap defining the linear order of all exterior terminals. It supplies a complementary noncrossing pairing of depth at most two that joins all paths of that color into one loop. Draw both complementary pairings in the exterior disc so same-color arcs are disjoint and two opposite-color arcs cross once exactly when their endpoints alternate. One may use straight chords after mapping that disc to a round disc, in general position. The depth-two crossing bound in the shallow-cap proof gives at most twice the number of red terminals, hence at most 6w added crossings. Small perturbations eliminate accidental multiple intersections without changing the alternation count.

These cap arcs lie outside the protective cycle and meet its existing strands only at the intended terminals. Each monochromatic curve is now one simple closed curve. Together with the at most 2w cycle vertices the construction adds at most 8w crossings. The total number of crossings is automatically even, by the intersection parity of two Jordan curves on a sphere.

## Inverse records

Root the completed meander at a cycle vertex, with the incident river direction specifying its usual rooted-meander representation. Allow a bounded extra choice for local orientation and the first directed cycle edge. Record the length L<=2w, its nonbacktracking step word, and which side of the cycle contains the original patch. There are at most three possible continuations at each vertex of a degree-four primal graph, even when the surrounding faces have large degree. Thus the number of such records is at most a constant times sum_(L<=2w) 3^L, which is A^w after increasing A. Extra port/switch flags, if desired, also cost only 2^L.

Starting at the root, the word reconstructs C. At its straight vertices the strand entering the interior is an original port; at its turning vertices no strand enters the interior. Delete C, the exterior, and the short extended portions of the original port strands. This recovers exactly the original patch and its cyclic order of ports. The deletion is an inverse reconstruction, not an operation on the live project files. The root of the original global meander need not be retained in each piece: the global counting argument records that one root separately at polynomial cost.

The recovery argument is uniform over v and w. It does not require that a boundary face have small degree, that ports occur consecutively along the river, or that the protective cycle be monochromatic.
