# CityShield - viva notes

One page per tab: **idea - correctness - time - space - CityShield use - likely questions**. Complexities are those documented in the Roxygen header of each function.

## How to open the viva (60 seconds)
*"CityShield plans emergency response on the real Chennai road network (OpenStreetMap, ~103,000 junctions, ~253,000 directed roads). Every algorithm the dashboard uses is written from scratch in R with explicit loops and hand-made queues, heaps, stacks and trees; each has a brute-force or library reference in the tests. A Shiny dashboard with 6 tabs shows the algorithms working on the city."*

Three design points to mention: (1) CSR graph + separate residual-edge flow graph; (2) mutable structures live in environments because R copies vectors on modification (two accidental O(n)-per-update bugs were found with benchmarks and fixed); (3) a test scans the source and fails if `which/order/outer/*apply/sort` is called.

---

## Dispatch
**Push-relabel max flow, used for ambulance matching.** Build a unit-capacity network (source -> ambulance -> incident -> sink); an integral max flow of value k is a matching of size k (flow integrality). Push-relabel keeps a *preflow* (nodes may hold excess) and heights; push along admissible edges (residual > 0, height[u] = height[v] + 1), relabel when stuck. Global relabelling (reverse BFS) gives exact starting heights, the gap heuristic lifts cut-off nodes. Phase 1 already gives the flow *value*; `complete = TRUE` returns stranded excess to the source so the pairs can be read off. FIFO selection O(V^3). *Q: why does the residual graph have reverse edges?* They let a later push undo earlier flow.

**Dijkstra with early exit (routes).** Greedy: settle the unsettled node with the smallest tentative distance (binary heap, lazy deletion); correct for non-negative weights. Stopping when the target is settled saves most of the work. O((V+E) log V), O(V+E). *Q: what breaks with negative weights?* A settled node could later be improved, so the greedy choice is no longer safe.

**Quicksort (ranking incidents).** Partition around a pivot (Lomuto), smaller side first with an explicit stack so the stack is O(log n) even on adversarial input. Random pivot gives 2 n ln n expected comparisons on *every* input; the last-element pivot is Theta(n^2) on sorted input. The dashboard uses it silently to pick the ten most urgent incidents. *Q: Monte Carlo or Las Vegas?* Las Vegas: always correct, random running time.

## Evacuation
**Push-relabel max flow** on the road network: super-source -> people at risk, road edges with capacity (vehicles per hour x 3 people), shelters -> super-sink. Max flow = min cut = how many people can leave in an hour.

**Sweep line (Bentley-Ottmann)** for the roads crossing the flood boundary. Events (segment start, end, discovered crossing) in x-order from a heap; the *status* is a treap ordered by y; two segments can only cross after they become neighbours, so only neighbours are tested. O((n+k) log n). *Q: what broke it on real data?* Touching segments and 1e-14 coordinate noise - fixed by snapping to a 1e-6 grid; regression tests keep it fixed.

**Graham scan** for the hull of the affected incidents: sort by angle around the lowest point, keep a stack of left turns. O(n log n).

## Resilience
**Karger's min cut.** Contract a uniformly random edge (union-find) until two super-nodes remain; a fixed minimum cut survives with probability >= 2/(n(n-1)), so n(n-1)/2 ln(1/delta) trials fail with probability <= delta. The dashboard fixes 100 trials and a fixed seed. Monte Carlo (may be wrong, always fast). Validated against igraph and an exhaustive reference *in the tests*.

**Divide-and-conquer maximum subarray** for the worst load window: best block is in the left half, the right half, or crosses the middle (best suffix + best prefix). T(n) = 2T(n/2) + Theta(n) = Theta(n log n).

## Placement
**Set cover (stations)** and **vertex cover (sensors)** are NP-hard. The dashboard uses the *exact* branch and bound whenever it finishes within its node budget and the *approximation* otherwise.
- Exact set cover: branch on the sets covering the rarest uncovered place; bound = chosen + ceil(uncovered / largest set).
- Greedy set cover: pick the set covering most uncovered places; |greedy| <= H(d) OPT <= (ln d + 1) OPT (charging argument).
- Exact vertex cover: branch on an uncovered edge (one endpoint must be in the cover), degree-1 rule, bound uncovered / max degree.
- 2-approximation: both endpoints of a maximal matching; the matched edges are disjoint so OPT >= |M| while the cover has 2|M| vertices.
*Q: when is B&B correct?* If the bound never overestimates the best completion. *Q: is the 2 tight?* Yes - a perfect matching forces ratio exactly 2.

## Logistics
**Fractional knapsack (supplies)** - rank by value/weight, take whole items, then a fraction. Correct by exchange argument. O(n log n). *Q: why does greedy fail for 0/1 knapsack?* The last item cannot be split, so the greedy-choice property breaks.

**0/1 knapsack DP (upgrades)** - K[i,c] = max(K[i-1,c], K[i-1,c-w_i] + v_i); Theta(nW), *pseudo-polynomial* (W is a number, not a size).

**TSP (inspection route)** - exact branch and bound with a reduced-cost-matrix bound (row/column minima must still be paid) when it finishes, otherwise the MST 2-approximation: MST <= OPT, a preorder walk with shortcuts costs <= 2 MST <= 2 OPT (needs the triangle inequality, which road travel times satisfy). The approximation's tour also gives the exact search its first upper bound.

## Intel
**KMP** - the prefix function says how much of the pattern is still matched after a mismatch; the matched length rises at most once per text character so total fallbacks <= n. O(n+m). *Q: why never re-read the text?* A mismatch at pattern position q tells us the last q text characters, so the longest border is the only possible restart.

**LCS** - L[i,j] = L[i-1,j-1] + 1 or max(L[i-1,j], L[i,j-1]); Theta(mn); two rows suffice for the length. Similarity = LCS length / longer report; reports above the threshold are duplicates. *Q: precision vs recall?* A low threshold flags template-similar reports (more correct duplicates found, more wrong flags); a high one misses noisy copies.

---

## Demo script (3 minutes)
1. **Dispatch** - pick the top incident, drag the trip limit: watch the matching and the route.
2. **Evacuation** - scenario 1: flood circle, flows, red exit roads, affected area.
3. **Resilience** - the red road is the weakest link; the shaded bars are the worst load period.
4. **Placement** - stations: radius 1 is too big for the exact search (approximation used), radius 2 and 3 use the exact search; sensors use the exact search.
5. **Logistics** - supplies slider, upgrades budget, inspection sites 4 -> 14 (exact up to 10 sites, then the approximation).
6. **Intel** - search "flood", then slide the duplicate threshold.
