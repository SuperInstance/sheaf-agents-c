# sheaf-agents-c

**Sheaf cohomology for multi-agent coordination** — C11 implementation.

## What is this?

This library models multi-agent systems using **cellular sheaves on graphs**. The key insight (from Step-3.5-Flash): **H¹ measures where agents structurally cannot agree.**

### Sheaf Cohomology in One Paragraph

A *cellular sheaf* assigns vector spaces (stalks) to graph vertices and linear maps (restriction maps) to edges. The **sheaf Laplacian** generalizes the graph Laplacian. Its kernel gives **H⁰** (global sections = possible agreements), and the cokernel gives **H¹** (obstructions = disagreements that *no amount of communication can resolve*).

## Key Concepts

| Concept | Meaning in Agent Terms |
|---------|----------------------|
| Stalk at vertex v | Agent v's local belief space |
| Restriction map on edge (v,w) | How agent v's beliefs project to what agent w can observe |
| H⁰ (global sections) | Space of beliefs where all agents agree |
| H¹ (obstructions) | **Structural disagreements — agents CANNOT agree, regardless of communication** |
| Sheaf Laplacian | Disagreement energy: L_sheaf measures total local inconsistency |

## H¹ as Structural Disagreement

When H¹ > 0, agents face a fundamental impossibility: their restriction maps create conflicting constraints. This isn't a communication problem — it's a **mathematical** one. No protocol, no amount of messaging, and no amount of patience will produce agreement.

Example: Three agents on a triangle where two edges use identity maps but one uses a 90° rotation. The rotation creates an obstruction — there's no global belief compatible with all three edges simultaneously.

## Building

```bash
make test
```

No dependencies beyond C11 + libm.

## API

### Cellular Sheaf

- `sheaf_create(n_verts, stalk_dims, n_edges)` — create a sheaf on a graph
- `sheaf_set_edge(s, idx, v1, v2, r1, r2)` — set restriction maps for an edge
- `sheaf_laplacian(s)` — compute the sheaf Laplacian
- `sheaf_cohomology(s, tol)` — compute H⁰ and H¹

### Agent Layer

- `agent_network_create(s)` — create agents from sheaf
- `agent_set_belief(net, agent, belief)` — set agent's local belief
- `agent_synchronize(net, step, converged)` — one synchronization step
- `agent_converge(net, step, max_iter, tol)` — iterate until convergence
- `agent_disagreement(net)` — measure total disagreement

### Consensus

- `can_agree(s, tol)` — check if H¹ = 0 (agreement is structurally possible)
- `forced_disagreement(s, tol)` — compute H¹ basis (what agents CANNOT agree on)
- `quality_metric(net, tol)` — score consensus quality

## License

MIT
