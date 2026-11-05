# 677: independent referee audit of 637 Proposition T

Scope: sections 0–1 of 637_rung2_matrix_product/PROOF.md against section 2, Theorem M, of 630_referee_band_growth/REPAIRS.md. No optimization, certificate generation, or numerical count was used. Verdict: Proposition T and its counting consequence are correct. The proof needs minor convention clarifications listed below, but no mathematical repair changes its statement or counts.

## 1. Source machine and pruning

637's label ranges agree with 630: a lower push has h' in [Lo,Hi-1], an upper push in [Lo+1,Hi]. In 630 the stated one-sided push inequalities plus h' in the band imply the same ranges. The step and pop conditions coincide exactly.

630 allows states with no successor and observes that forbidding entry changes no walk ending at t. In the reachable invariant class used by 637, both tops marked means no permitted direction; therefore such a state cannot occur on a return walk. The tape machine represents the pruned machine. This should be stated explicitly when invoking Theorem M: reachable dead states in the unpruned machine are discarded, and all reachable states below mean reachable in this pruned machine.

## 2. Invariants and tape recovery

Consecutive labels in a nonempty stack differ by one. A push preserves this because its new label is h' and the previous top differs by one; a pop exposes a label h' or h'-2 below on the lower side, and h' or h'+2 above. The empty-stack case is separate: the only disappearing bottom label is 0. Lower emptiness is created at h'=1 and the next permitted step is to h'=0, where the sole lower action is to push 0. Upper emptiness is its reflected case at h'=-1. This also verifies that no adjacent pair of empty stacks is reachable.

Thus the two sentinel levels are legitimate. With an empty lower stack set l=-1; with an empty upper stack set u=1. The head difference is (u-1)-l. It equals -1 when both tops are unmarked, and +1 when exactly one is marked. Both-marked states are excluded. The end steps and every internal stack step have magnitude one; all vertices lie in [Lo,Hi-1].

The head position uniquely cuts the path into the lower labels and reversed shifted upper labels. Hence a tape determines both stacks, including their lengths. For head '-', h=l is unique. For head '+', the only possibilities are h=l or h=l+2. Both possibilities need not belong to the band or be reachable; the fiber statement is a restriction on whichever preimages exist. The phrase 'both twins' in the local proof can safely be read as 'either existing twin'.

## 3. Exhaustive outgoing-move audit

For a '-' head both tops equal h and both stacks are nonempty. A positive step must push above and either push or pop below. The two outputs are respectively the split '- to +[-]-' and a left shift. The upper push is always permitted because h<=Hi-1, while the lower push requires h+1<=Hi-1, exactly the extra strip vertex in that split. The negative step gives the reflected split '- to -[-]+' or right shift, with the sole extra requirement h-2>=Lo. A shift needing an absent neighboring letter is excluded precisely by the absent pop.

For a '+' head u=l+2. Either existing twin has the unique permitted next winding c=l+1. Both pushes are allowed: c<=Hi-1 follows from the upper label or its empty sentinel, and c>=Lo+1 follows from the lower label or its empty sentinel. The four choices are: two pushes give '+ to +[-]+'; lower pop gives left shift; upper pop gives right shift; two pops give merge. A missing stack rules out exactly the corresponding pop, matching the tape boundary.

For two pops, let r be the exposed lower label and r'' the shifted exposed upper label, using sentinels when a pop empties a stack. The new state is dead exactly when r=c-2 and r''=c+1. Then r''-r=3, precisely the local word '+++'. In all three other neighbor-sign choices r''-r is +1 or -1 and the merge is legal. This includes pops that empty one stack: the sentinel calculation still applies.

All outputs are distinct as tapes. Shifts retain length and have different head positions. Splits increase length by two; at a '-' head their different outer letters distinguish the two possible splits. A merge decreases length by two. Therefore this is a bijection on individual outgoing moves as well as successor states. No hidden multiplicity is lost.

## 4. Quotient path lifting and terminal fiber

When the head is '+', all source twins move to the same winding c and have the same stacks. Their actual label successor states, not merely their tape images, coincide. Each individual source still has exactly one successor above each allowed tape successor, by the preceding case audit.

A tape walk starting at the root consequently lifts uniquely after fixing the initial state t: inductively the current label state fixes its unique next state over the next tape. Conversely projecting a label walk gives that tape walk. No arbitrary factor two is introduced at any twin visit.

The root tape has stacks L=U=(0) and head '-'. Its winding is uniquely h=0, so its fiber is exactly {t}. Thus a lift ends at t if and only if its tape ends at the root. This proves exact preservation of return-walk counts, rather than only existence of a state quotient. Theorem M then gives the claimed band meander count with 2m-2 moves. The case m=1 is the empty return walk and is included.

## 5. Uncapped exhaustion

The uncapped invariant and local audit hold with all label bounds removed. Alternatively take any fixed finite tape walk: its tape words and vertices form a finite set, and choose a band with margin at least two around every visited vertex as well as Lo<=-1, Hi>=1. Its local lifts then satisfy the capped proof. The root fiber argument is unchanged. Conversely every uncapped label return walk is finite and fits a finite band. Hence the uncapped tape and label return walks are in bijection.

Every order-m meander has winding obtained from 2m-2 signs starting at zero, so it lies in [-2m,2m]. Applying Theorem M to this band shows that all meanders occur in the uncapped correspondence. This makes explicit the exhaustion argument compressed into 637's last sentence.

## Referee disposition

Accept Proposition T and Corollary T1. The only editorial amendments recommended are: specify pruning when comparing to 630; say 'either existing twin' in bands; and insert the finite-band exhaustion argument. They do not require changed transition rules or any correction factor. The later numerical or local-window certificates were outside this review.
