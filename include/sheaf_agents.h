#ifndef SHEAF_AGENTS_H
#define SHEAF_AGENTS_H

#include <stddef.h>
#include <stdbool.h>

#ifdef __cplusplus
extern "C" {
#endif

/* ── Linear algebra helpers ── */

typedef struct {
    double *data;
    int rows, cols;
} Matrix;

Matrix matrix_alloc(int rows, int cols);
void   matrix_free(Matrix *m);
Matrix matrix_zeros(int rows, int cols);
Matrix matrix_identity(int n);
Matrix matrix_transpose(const Matrix *m);
Matrix matrix_mul(const Matrix *a, const Matrix *b);
Matrix matrix_add(const Matrix *a, const Matrix *b);
Matrix matrix_sub(const Matrix *a, const Matrix *b);
double matrix_frobenius(const Matrix *m);
void   matrix_print(const Matrix *m);

/* ── Cellular Sheaf on Graphs ── */

typedef struct {
    int v1, v2;               /* edge endpoints */
    Matrix restriction;       /* maps stalk at v1 → common edge space */
    Matrix restriction_t;     /* maps stalk at v2 → common edge space */
} SheafEdge;

typedef struct {
    int        n_verts;
    int       *stalk_dims;    /* stalk_dims[v] = dimension of stalk at v */
    int        n_edges;
    SheafEdge *edges;
} CellularSheaf;

CellularSheaf sheaf_create(int n_verts, const int *stalk_dims, int n_edges);
void          sheaf_free(CellularSheaf *s);
void          sheaf_set_edge(CellularSheaf *s, int idx,
                             int v1, int v2,
                             const Matrix *r1, const Matrix *r2);

/* Sheaf Laplacian: L = B^T diag(F_e) B, generalizes graph Laplacian */
Matrix sheaf_laplacian(const CellularSheaf *s);

/* Cohomology: H⁰ = ker(d₀), H¹ = coker(d₀) */
typedef struct {
    int h0_dim;     /* dimension of global sections (agreements) */
    int h1_dim;     /* dimension of obstructions (structural disagreements) */
    Matrix h0_basis; /* basis for global sections, columns are basis vectors */
    Matrix h1_basis; /* basis for obstructions */
} Cohomology;

Cohomology sheaf_cohomology(const CellularSheaf *s, double tol);
void       cohomology_free(Cohomology *c);

/* ── Agent Layer ── */

typedef struct {
    int    vertex;       /* which vertex this agent sits on */
    Matrix belief;       /* column vector in stalk space */
} AgentState;

typedef struct {
    CellularSheaf sheaf;
    int           n_agents;
    AgentState   *agents;
} AgentNetwork;

AgentNetwork agent_network_create(const CellularSheaf *s);
void         agent_network_free(AgentNetwork *net);
void         agent_set_belief(AgentNetwork *net, int agent, const Matrix *belief);

/* Synchronize: push beliefs through restriction maps */
Matrix agent_synchronize(const AgentNetwork *net, double step, int *converged);

/* Measure disagreement (H¹ obstruction in current state) */
double agent_disagreement(const AgentNetwork *net);

/* Iterate until convergence or detect obstruction */
typedef struct {
    int      iterations;
    double   final_disagreement;
    bool     converged;
    bool     obstruction_detected;
} ConvergenceResult;

ConvergenceResult agent_converge(AgentNetwork *net, double step, int max_iter, double tol);

/* ── Consensus ── */

typedef struct {
    double agreement_score;    /* 1.0 = perfect agreement */
    double h0_quality;         /* quality of global section */
    double h1_obstruction;     /* measure of structural disagreement */
    int    h0_dim;
    int    h1_dim;
} ConsensusQuality;

/* Can agents agree? H¹ = 0 means yes */
bool can_agree(const CellularSheaf *s, double tol);

/* Compute H¹ basis — what agents CANNOT agree on */
Matrix forced_disagreement(const CellularSheaf *s, double tol);

/* Quality metric for consensus */
ConsensusQuality quality_metric(const AgentNetwork *net, double tol);

#ifdef __cplusplus
}
#endif

#endif /* SHEAF_AGENTS_H */
