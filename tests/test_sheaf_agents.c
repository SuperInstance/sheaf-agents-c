#include "sheaf_agents.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

static int tests_passed = 0;
static int tests_failed = 0;

#define ASSERT(cond, msg) do { \
    if (cond) { tests_passed++; } \
    else { tests_failed++; printf("FAIL: %s\n", msg); } \
} while(0)

#define ASSERT_EQ(a, b, msg) ASSERT((a) == (b), msg)
#define ASSERT_APPROX(a, b, tol, msg) ASSERT(fabs((a) - (b)) < (tol), msg)

/* ── Test 1: Sheaf on single edge — H⁰ = dim(stalk), H¹ = 0 ── */
static void test_single_edge(void) {
    int dims[] = {2, 2};
    CellularSheaf s = sheaf_create(2, dims, 1);

    Matrix id2 = matrix_identity(2);
    sheaf_set_edge(&s, 0, 0, 1, &id2, &id2);

    Cohomology c = sheaf_cohomology(&s, 1e-8);
    ASSERT_EQ(c.h0_dim, 2, "single edge: H0 = dim(stalk) = 2");
    ASSERT_EQ(c.h1_dim, 0, "single edge: H1 = 0");

    cohomology_free(&c);
    matrix_free(&id2);
    sheaf_free(&s);
}

/* ── Test 2: Triangle with identity restrictions — perfect agreement ── */
static void test_triangle_identity(void) {
    int dims[] = {2, 2, 2};
    CellularSheaf s = sheaf_create(3, dims, 3);

    Matrix id2 = matrix_identity(2);
    sheaf_set_edge(&s, 0, 0, 1, &id2, &id2);
    sheaf_set_edge(&s, 1, 1, 2, &id2, &id2);
    sheaf_set_edge(&s, 2, 0, 2, &id2, &id2);

    Cohomology c = sheaf_cohomology(&s, 1e-8);
    ASSERT_EQ(c.h0_dim, 2, "triangle identity: H0 = 2");
    ASSERT_EQ(c.h1_dim, 0, "triangle identity: H1 = 0 (compatible restrictions)");

    cohomology_free(&c);
    matrix_free(&id2);
    sheaf_free(&s);
}

/* ── Test 3: Triangle with orthogonal restrictions — H¹ > 0 ── */
static void test_triangle_orthogonal(void) {
    int dims[] = {2, 2, 2};
    CellularSheaf s = sheaf_create(3, dims, 3);

    Matrix r1 = matrix_identity(2);
    Matrix r2 = matrix_zeros(2, 2);
    r2.data[0 * 2 + 1] = 1.0;  /* [0 1] */
    r2.data[1 * 2 + 0] = -1.0; /* [-1 0] */

    /* Edge 0: v0→v1 with identity/rotation
     * Edge 1: v1→v2 with identity/rotation
     * Edge 2: v0→v2 with identity/identity (creates cycle inconsistency) */
    sheaf_set_edge(&s, 0, 0, 1, &r1, &r2);
    sheaf_set_edge(&s, 1, 1, 2, &r1, &r2);
    sheaf_set_edge(&s, 2, 0, 2, &r1, &r1);

    Cohomology c = sheaf_cohomology(&s, 1e-6);
    /* Going around the cycle: v0 → v1 (via rot90) → v2 (via rot90) → v0 (via identity)
     * Composition = rot90 * rot90 * identity = rot180 = -identity ≠ identity
     * So H⁰ = 0 (no global sections) and H¹ = 2 > 0 */
    ASSERT_EQ(c.h0_dim, 0, "triangle orthogonal: H0 = 0 (no global sections)");
    ASSERT(c.h1_dim > 0, "triangle orthogonal: H1 > 0 (structural obstruction)");

    cohomology_free(&c);
    matrix_free(&r1); matrix_free(&r2);
    sheaf_free(&s);
}

/* ── Test 4: Star topology — central sees all ── */
static void test_star_topology(void) {
    int dims[] = {3, 1, 1, 1};
    CellularSheaf s = sheaf_create(4, dims, 3);

    for (int i = 0; i < 3; i++) {
        Matrix r_center = matrix_zeros(1, 3);
        r_center.data[i] = 1.0;
        Matrix r_leaf = matrix_identity(1);
        sheaf_set_edge(&s, i, 0, i + 1, &r_center, &r_leaf);
        matrix_free(&r_center);
        matrix_free(&r_leaf);
    }

    Cohomology c = sheaf_cohomology(&s, 1e-8);
    ASSERT(c.h0_dim >= 1, "star: H0 >= 1 (can agree on center's view)");
    ASSERT(c.h1_dim >= 0, "star: H1 >= 0");

    cohomology_free(&c);
    sheaf_free(&s);
}

/* ── Test 5: Agent sync converges when H¹ = 0 ── */
static void test_agent_converge_h1_zero(void) {
    int dims[] = {2, 2};
    CellularSheaf s = sheaf_create(2, dims, 1);
    Matrix id2 = matrix_identity(2);
    sheaf_set_edge(&s, 0, 0, 1, &id2, &id2);

    AgentNetwork net = agent_network_create(&s);
    Matrix b0 = matrix_zeros(2, 1);
    b0.data[0] = 1.0; b0.data[1] = 2.0;
    Matrix b1 = matrix_zeros(2, 1);
    b1.data[0] = 3.0; b1.data[1] = 4.0;
    agent_set_belief(&net, 0, &b0);
    agent_set_belief(&net, 1, &b1);

    ConvergenceResult cr = agent_converge(&net, 0.1, 500, 1e-6);
    ASSERT(cr.converged, "agent converge: converges when H1=0");
    ASSERT_APPROX(cr.final_disagreement, 0.0, 1e-4, "agent converge: disagreement -> 0");

    matrix_free(&b0); matrix_free(&b1);
    agent_network_free(&net);
    matrix_free(&id2);
    sheaf_free(&s);
}

/* ── Test 6: Triangle obstruction detected ── */
static void test_agent_obstruction(void) {
    int dims[] = {2, 2, 2};
    CellularSheaf s = sheaf_create(3, dims, 3);

    Matrix r1 = matrix_identity(2);
    Matrix r2 = matrix_zeros(2, 2);
    r2.data[0 * 2 + 1] = 1.0;
    r2.data[1 * 2 + 0] = -1.0;

    sheaf_set_edge(&s, 0, 0, 1, &r1, &r2);
    sheaf_set_edge(&s, 1, 1, 2, &r1, &r2);
    sheaf_set_edge(&s, 2, 0, 2, &r1, &r1);

    bool agree = can_agree(&s, 1e-6);
    ASSERT(!agree, "agent obstruction: can_agree returns false for inconsistent cycle");

    matrix_free(&r1); matrix_free(&r2);
    sheaf_free(&s);
}

/* ── Test 7: Quality metric decreases with H¹ ── */
static void test_quality_vs_h1(void) {
    /* H¹ = 0 case: triangle with identity */
    int dims1[] = {2, 2, 2};
    CellularSheaf s1 = sheaf_create(3, dims1, 3);
    Matrix id2 = matrix_identity(2);
    sheaf_set_edge(&s1, 0, 0, 1, &id2, &id2);
    sheaf_set_edge(&s1, 1, 1, 2, &id2, &id2);
    sheaf_set_edge(&s1, 2, 0, 2, &id2, &id2);

    AgentNetwork net1 = agent_network_create(&s1);
    Matrix b = matrix_zeros(2, 1); b.data[0] = 1.0;
    for (int v = 0; v < 3; v++) agent_set_belief(&net1, v, &b);

    ConsensusQuality q1 = quality_metric(&net1, 1e-8);

    /* H¹ > 0 case: triangle with orthogonal */
    int dims2[] = {2, 2, 2};
    CellularSheaf s2 = sheaf_create(3, dims2, 3);
    Matrix r1 = matrix_identity(2);
    Matrix r2 = matrix_zeros(2, 2);
    r2.data[1] = 1.0; r2.data[2] = -1.0;
    sheaf_set_edge(&s2, 0, 0, 1, &r1, &r2);
    sheaf_set_edge(&s2, 1, 1, 2, &r1, &r2);
    sheaf_set_edge(&s2, 2, 0, 2, &r1, &r1);

    AgentNetwork net2 = agent_network_create(&s2);
    for (int v = 0; v < 3; v++) agent_set_belief(&net2, v, &b);

    ConsensusQuality q2 = quality_metric(&net2, 1e-6);

    ASSERT(q1.h1_dim == 0, "quality: H1=0 case has h1_dim=0");
    ASSERT(q2.h1_dim > 0, "quality: H1>0 case has h1_dim>0");
    /* Quality of agreement space is worse when H¹ > 0 */
    ASSERT(q1.h0_dim > q2.h0_dim, "quality: more global sections when H1=0");

    matrix_free(&b);
    agent_network_free(&net1);
    agent_network_free(&net2);
    matrix_free(&id2); matrix_free(&r1); matrix_free(&r2);
    sheaf_free(&s1); sheaf_free(&s2);
}

/* ── Test 8: Empty sheaf ── */
static void test_empty_sheaf(void) {
    CellularSheaf s = sheaf_create(0, NULL, 0);
    Cohomology c = sheaf_cohomology(&s, 1e-8);
    ASSERT_EQ(c.h0_dim, 0, "empty sheaf: H0 = 0");
    ASSERT_EQ(c.h1_dim, 0, "empty sheaf: H1 = 0");
    cohomology_free(&c);
    sheaf_free(&s);
}

/* ── Test 9: Single vertex ── */
static void test_single_vertex(void) {
    int dims[] = {3};
    CellularSheaf s = sheaf_create(1, dims, 0);
    Cohomology c = sheaf_cohomology(&s, 1e-8);
    ASSERT_EQ(c.h0_dim, 3, "single vertex: H0 = 3 (whole stalk)");
    ASSERT_EQ(c.h1_dim, 0, "single vertex: H1 = 0");
    cohomology_free(&c);
    sheaf_free(&s);
}

/* ── Test 10: Disconnected graph ── */
static void test_disconnected(void) {
    int dims[] = {2, 2};
    CellularSheaf s = sheaf_create(2, dims, 0);
    Cohomology c = sheaf_cohomology(&s, 1e-8);
    ASSERT_EQ(c.h0_dim, 4, "disconnected: H0 = 2+2 = 4");
    ASSERT_EQ(c.h1_dim, 0, "disconnected: H1 = 0");
    cohomology_free(&c);
    sheaf_free(&s);
}

/* ── Test 11: Sheaf Laplacian is symmetric ── */
static void test_laplacian_psd(void) {
    int dims[] = {2, 2};
    CellularSheaf s = sheaf_create(2, dims, 1);
    Matrix id2 = matrix_identity(2);
    sheaf_set_edge(&s, 0, 0, 1, &id2, &id2);

    Matrix L = sheaf_laplacian(&s);

    bool symmetric = true;
    for (int i = 0; i < L.rows; i++)
        for (int j = 0; j < L.cols; j++)
            if (fabs(L.data[i * L.cols + j] - L.data[j * L.cols + i]) > 1e-10)
                symmetric = false;
    ASSERT(symmetric, "laplacian: symmetric");

    matrix_free(&L);
    matrix_free(&id2);
    sheaf_free(&s);
}

/* ── Test 12: Laplacian structure for identity sheaf ── */
static void test_laplacian_h0(void) {
    int dims[] = {2, 2};
    CellularSheaf s = sheaf_create(2, dims, 1);
    Matrix id2 = matrix_identity(2);
    sheaf_set_edge(&s, 0, 0, 1, &id2, &id2);

    Matrix L = sheaf_laplacian(&s);
    ASSERT_APPROX(L.data[0], 1.0, 1e-10, "laplacian h0: L[0,0] = 1");
    ASSERT_APPROX(L.data[5], 1.0, 1e-10, "laplacian h0: L[3,3] = 1");
    ASSERT_APPROX(L.data[2], -1.0, 1e-10, "laplacian h0: L[0,2] = -1");

    matrix_free(&L);
    matrix_free(&id2);
    sheaf_free(&s);
}

/* ── Test 13: can_agree true for identity sheaf ── */
static void test_can_agree_true(void) {
    int dims[] = {2, 2, 2};
    CellularSheaf s = sheaf_create(3, dims, 3);
    Matrix id2 = matrix_identity(2);
    sheaf_set_edge(&s, 0, 0, 1, &id2, &id2);
    sheaf_set_edge(&s, 1, 1, 2, &id2, &id2);
    sheaf_set_edge(&s, 2, 0, 2, &id2, &id2);

    ASSERT(can_agree(&s, 1e-8), "can_agree: true for identity sheaf");

    matrix_free(&id2);
    sheaf_free(&s);
}

/* ── Test 14: can_agree false for inconsistent cycle ── */
static void test_can_agree_false(void) {
    int dims[] = {2, 2, 2};
    CellularSheaf s = sheaf_create(3, dims, 3);
    Matrix r1 = matrix_identity(2);
    Matrix r2 = matrix_zeros(2, 2);
    r2.data[1] = 1.0; r2.data[2] = -1.0;
    sheaf_set_edge(&s, 0, 0, 1, &r1, &r2);
    sheaf_set_edge(&s, 1, 1, 2, &r1, &r2);
    sheaf_set_edge(&s, 2, 0, 2, &r1, &r1);

    ASSERT(!can_agree(&s, 1e-6), "can_agree: false for inconsistent cycle");

    matrix_free(&r1); matrix_free(&r2);
    sheaf_free(&s);
}

/* ── Test 15: Agent disagreement measure ── */
static void test_agent_disagreement(void) {
    int dims[] = {2, 2};
    CellularSheaf s = sheaf_create(2, dims, 1);
    Matrix id2 = matrix_identity(2);
    sheaf_set_edge(&s, 0, 0, 1, &id2, &id2);

    AgentNetwork net = agent_network_create(&s);
    Matrix b0 = matrix_zeros(2, 1); b0.data[0] = 1.0;
    Matrix b1 = matrix_zeros(2, 1); b1.data[0] = 2.0;
    agent_set_belief(&net, 0, &b0);
    agent_set_belief(&net, 1, &b1);

    double dis = agent_disagreement(&net);
    ASSERT(dis > 0, "agent disagreement: positive when beliefs differ");

    matrix_free(&b0); matrix_free(&b1);
    agent_network_free(&net);
    matrix_free(&id2);
    sheaf_free(&s);
}

/* ── Test 16: Agent disagreement zero when beliefs agree ── */
static void test_agent_agreement(void) {
    int dims[] = {2, 2};
    CellularSheaf s = sheaf_create(2, dims, 1);
    Matrix id2 = matrix_identity(2);
    sheaf_set_edge(&s, 0, 0, 1, &id2, &id2);

    AgentNetwork net = agent_network_create(&s);
    Matrix b = matrix_zeros(2, 1); b.data[0] = 1.0; b.data[1] = 2.0;
    agent_set_belief(&net, 0, &b);
    agent_set_belief(&net, 1, &b);

    double dis = agent_disagreement(&net);
    ASSERT_APPROX(dis, 0.0, 1e-10, "agent agreement: zero disagreement when beliefs match");

    matrix_free(&b);
    agent_network_free(&net);
    matrix_free(&id2);
    sheaf_free(&s);
}

/* ── Test 17: forced_disagreement returns non-empty for H¹ > 0 ── */
static void test_forced_disagreement(void) {
    int dims[] = {2, 2, 2};
    CellularSheaf s = sheaf_create(3, dims, 3);
    Matrix r1 = matrix_identity(2);
    Matrix r2 = matrix_zeros(2, 2);
    r2.data[1] = 1.0; r2.data[2] = -1.0;
    sheaf_set_edge(&s, 0, 0, 1, &r1, &r2);
    sheaf_set_edge(&s, 1, 1, 2, &r1, &r2);
    sheaf_set_edge(&s, 2, 0, 2, &r1, &r1);

    Matrix fd = forced_disagreement(&s, 1e-6);
    ASSERT(fd.cols > 0, "forced disagreement: returns basis vectors for H1>0");
    ASSERT(fd.rows == 6, "forced disagreement: correct total dimension");

    matrix_free(&fd);
    matrix_free(&r1); matrix_free(&r2);
    sheaf_free(&s);
}

/* ── Test 18: forced_disagreement empty for H¹ = 0 ── */
static void test_no_forced_disagreement(void) {
    int dims[] = {2, 2};
    CellularSheaf s = sheaf_create(2, dims, 1);
    Matrix id2 = matrix_identity(2);
    sheaf_set_edge(&s, 0, 0, 1, &id2, &id2);

    Matrix fd = forced_disagreement(&s, 1e-8);
    ASSERT_EQ(fd.cols, 0, "no forced disagreement: empty when H1=0");

    matrix_free(&fd);
    matrix_free(&id2);
    sheaf_free(&s);
}

/* ── Test 19: Agent network creation/free ── */
static void test_agent_network_lifecycle(void) {
    int dims[] = {3, 2, 1};
    CellularSheaf s = sheaf_create(3, dims, 0);
    AgentNetwork net = agent_network_create(&s);

    ASSERT_EQ(net.n_agents, 3, "agent lifecycle: 3 agents created");

    agent_network_free(&net);
    sheaf_free(&s);
}

/* ── Test 20: Synchronization reduces disagreement ── */
static void test_sync_reduces(void) {
    int dims[] = {2, 2};
    CellularSheaf s = sheaf_create(2, dims, 1);
    Matrix id2 = matrix_identity(2);
    sheaf_set_edge(&s, 0, 0, 1, &id2, &id2);

    AgentNetwork net = agent_network_create(&s);
    Matrix b0 = matrix_zeros(2, 1); b0.data[0] = 5.0; b0.data[1] = 3.0;
    Matrix b1 = matrix_zeros(2, 1); b1.data[0] = 1.0; b1.data[1] = 7.0;
    agent_set_belief(&net, 0, &b0);
    agent_set_belief(&net, 1, &b1);

    double dis_before = agent_disagreement(&net);
    int conv;
    Matrix x = agent_synchronize(&net, 0.1, &conv);
    double dis_after = agent_disagreement(&net);

    ASSERT(dis_after < dis_before, "sync: disagreement decreases after step");

    matrix_free(&x);
    matrix_free(&b0); matrix_free(&b1);
    agent_network_free(&net);
    matrix_free(&id2);
    sheaf_free(&s);
}

/* ── Test 21: H⁰ basis vectors are actual global sections ── */
static void test_h0_basis_global(void) {
    int dims[] = {2, 2};
    CellularSheaf s = sheaf_create(2, dims, 1);
    Matrix id2 = matrix_identity(2);
    sheaf_set_edge(&s, 0, 0, 1, &id2, &id2);

    Cohomology c = sheaf_cohomology(&s, 1e-8);
    ASSERT_EQ(c.h0_basis.cols, 2, "h0 basis: 2 global sections");

    /* For identity sheaf, global section: x0 = x1 */
    /* First basis vector should have x0[0] proportional to x1[0] */
    double x0_0 = c.h0_basis.data[0 * 2 + 0]; /* v0, component 0 of basis 0 */
    double x1_0 = c.h0_basis.data[2 * 2 + 0]; /* v1, component 0 of basis 0 */
    /* They should be equal (or proportional) */
    ASSERT_APPROX(x0_0, x1_0, 1e-6, "h0 basis: global section has matching v0=v1");

    cohomology_free(&c);
    matrix_free(&id2);
    sheaf_free(&s);
}

/* ── Test 22: Single vertex agent network ── */
static void test_single_agent(void) {
    int dims[] = {3};
    CellularSheaf s = sheaf_create(1, dims, 0);

    AgentNetwork net = agent_network_create(&s);
    Matrix b = matrix_zeros(3, 1);
    b.data[0] = 1.0; b.data[1] = 2.0; b.data[2] = 3.0;
    agent_set_belief(&net, 0, &b);

    double dis = agent_disagreement(&net);
    ASSERT_APPROX(dis, 0.0, 1e-10, "single agent: zero disagreement");

    matrix_free(&b);
    agent_network_free(&net);
    sheaf_free(&s);
}

/* ── Test 23: Matrix operations ── */
static void test_matrix_ops(void) {
    Matrix a = matrix_zeros(2, 3);
    a.data[0] = 1.0; a.data[1] = 2.0; a.data[2] = 3.0;
    a.data[3] = 4.0; a.data[4] = 5.0; a.data[5] = 6.0;

    Matrix at = matrix_transpose(&a);
    ASSERT_EQ(at.rows, 3, "transpose: rows");
    ASSERT_EQ(at.cols, 2, "transpose: cols");
    ASSERT_APPROX(at.data[0], 1.0, 1e-10, "transpose: value check");
    ASSERT_APPROX(at.data[1], 4.0, 1e-10, "transpose: value check");

    Matrix id3 = matrix_identity(3);
    Matrix prod = matrix_mul(&a, &id3);
    for (int i = 0; i < 6; i++)
        ASSERT_APPROX(prod.data[i], a.data[i], 1e-10, "mul by identity");

    matrix_free(&a); matrix_free(&at);
    matrix_free(&id3); matrix_free(&prod);
}

/* ── Test 24: Convergence result fields ── */
static void test_convergence_result(void) {
    int dims[] = {2, 2};
    CellularSheaf s = sheaf_create(2, dims, 1);
    Matrix id2 = matrix_identity(2);
    sheaf_set_edge(&s, 0, 0, 1, &id2, &id2);

    AgentNetwork net = agent_network_create(&s);
    Matrix b0 = matrix_zeros(2, 1); b0.data[0] = 1.0;
    Matrix b1 = matrix_zeros(2, 1); b1.data[0] = -1.0;
    agent_set_belief(&net, 0, &b0);
    agent_set_belief(&net, 1, &b1);

    ConvergenceResult cr = agent_converge(&net, 0.1, 500, 1e-6);
    ASSERT(cr.converged, "convergence result: converged");
    ASSERT(!cr.obstruction_detected, "convergence result: no obstruction");
    ASSERT(cr.iterations > 0, "convergence result: iterations > 0");

    matrix_free(&b0); matrix_free(&b1);
    agent_network_free(&net);
    matrix_free(&id2);
    sheaf_free(&s);
}

/* ── Test 25: Quality metric for identity sheaf ── */
static void test_quality_metric_identity(void) {
    int dims[] = {2, 2};
    CellularSheaf s = sheaf_create(2, dims, 1);
    Matrix id2 = matrix_identity(2);
    sheaf_set_edge(&s, 0, 0, 1, &id2, &id2);

    AgentNetwork net = agent_network_create(&s);
    Matrix b = matrix_zeros(2, 1); b.data[0] = 1.0;
    agent_set_belief(&net, 0, &b);
    agent_set_belief(&net, 1, &b);

    ConsensusQuality q = quality_metric(&net, 1e-8);
    ASSERT_APPROX(q.agreement_score, 1.0, 1e-6, "quality: agreement = 1 when beliefs match");
    ASSERT(q.h0_dim > 0, "quality: h0_dim > 0 for identity sheaf");
    ASSERT_EQ(q.h1_dim, 0, "quality: h1_dim = 0 for identity sheaf");

    matrix_free(&b);
    agent_network_free(&net);
    matrix_free(&id2);
    sheaf_free(&s);
}

/* ── Test 26: Larger graph — path of 4 vertices ── */
static void test_path_graph(void) {
    int dims[] = {2, 2, 2, 2};
    CellularSheaf s = sheaf_create(4, dims, 3);
    Matrix id2 = matrix_identity(2);
    sheaf_set_edge(&s, 0, 0, 1, &id2, &id2);
    sheaf_set_edge(&s, 1, 1, 2, &id2, &id2);
    sheaf_set_edge(&s, 2, 2, 3, &id2, &id2);

    Cohomology c = sheaf_cohomology(&s, 1e-8);
    ASSERT_EQ(c.h0_dim, 2, "path graph: H0 = 2");
    ASSERT_EQ(c.h1_dim, 0, "path graph: H1 = 0 (tree, no cycle obstruction)");

    cohomology_free(&c);
    matrix_free(&id2);
    sheaf_free(&s);
}

/* ── Test 27: Frobenius norm ── */
static void test_frobenius(void) {
    Matrix m = matrix_zeros(2, 2);
    m.data[0] = 3.0; m.data[3] = 4.0;
    ASSERT_APPROX(matrix_frobenius(&m), 5.0, 1e-10, "frobenius: 3,4 diagonal = 5");
    matrix_free(&m);
}

/* ── Test 28: Single edge always has H¹ = 0 (no cycles) ── */
static void test_single_edge_h1_zero(void) {
    int dims[] = {2, 2};
    CellularSheaf s = sheaf_create(2, dims, 1);

    /* Even with mismatched restrictions, single edge has H¹ = 0 */
    Matrix r1 = matrix_identity(2);
    Matrix r2 = matrix_zeros(2, 2);
    r2.data[1] = 1.0; r2.data[2] = -1.0;
    sheaf_set_edge(&s, 0, 0, 1, &r1, &r2);

    Cohomology c = sheaf_cohomology(&s, 1e-8);
    ASSERT_EQ(c.h0_dim, 2, "single edge mismatched: H0 = 2");
    ASSERT_EQ(c.h1_dim, 0, "single edge mismatched: H1 = 0 (no cycle)");

    cohomology_free(&c);
    matrix_free(&r1); matrix_free(&r2);
    sheaf_free(&s);
}

/* ── Test 29: Triangle with one bad edge, H¹ = 2 ── */
static void test_triangle_h1_dimension(void) {
    int dims[] = {2, 2, 2};
    CellularSheaf s = sheaf_create(3, dims, 3);
    Matrix r1 = matrix_identity(2);
    Matrix r2 = matrix_zeros(2, 2);
    r2.data[1] = 1.0; r2.data[2] = -1.0;
    sheaf_set_edge(&s, 0, 0, 1, &r1, &r2);
    sheaf_set_edge(&s, 1, 1, 2, &r1, &r2);
    sheaf_set_edge(&s, 2, 0, 2, &r1, &r1);

    Cohomology c = sheaf_cohomology(&s, 1e-6);
    ASSERT_EQ(c.h1_dim, 2, "triangle orthogonal: H1 = 2 (full obstruction)");

    cohomology_free(&c);
    matrix_free(&r1); matrix_free(&r2);
    sheaf_free(&s);
}

/* ── Main ── */
int main(void) {
    printf("Running sheaf-agents-c tests...\n\n");

    test_matrix_ops();
    test_frobenius();
    test_single_edge();
    test_triangle_identity();
    test_triangle_orthogonal();
    test_star_topology();
    test_agent_converge_h1_zero();
    test_agent_obstruction();
    test_quality_vs_h1();
    test_empty_sheaf();
    test_single_vertex();
    test_disconnected();
    test_laplacian_psd();
    test_laplacian_h0();
    test_can_agree_true();
    test_can_agree_false();
    test_agent_disagreement();
    test_agent_agreement();
    test_forced_disagreement();
    test_no_forced_disagreement();
    test_agent_network_lifecycle();
    test_sync_reduces();
    test_h0_basis_global();
    test_single_agent();
    test_convergence_result();
    test_quality_metric_identity();
    test_path_graph();
    test_single_edge_h1_zero();
    test_triangle_h1_dimension();

    printf("\n%d passed, %d failed, %d total\n",
           tests_passed, tests_failed, tests_passed + tests_failed);

    return tests_failed > 0 ? 1 : 0;
}
