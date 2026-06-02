#include "sheaf_agents.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

/* ── Matrix operations ── */

Matrix matrix_alloc(int rows, int cols) {
    Matrix m;
    m.rows = rows;
    m.cols = cols;
    m.data = (rows && cols) ? calloc((size_t)rows * cols, sizeof(double)) : NULL;
    return m;
}

void matrix_free(Matrix *m) {
    free(m->data);
    m->data = NULL;
    m->rows = m->cols = 0;
}

Matrix matrix_zeros(int rows, int cols) {
    return matrix_alloc(rows, cols);
}

Matrix matrix_identity(int n) {
    Matrix m = matrix_zeros(n, n);
    for (int i = 0; i < n; i++) m.data[i * n + i] = 1.0;
    return m;
}

Matrix matrix_transpose(const Matrix *m) {
    Matrix t = matrix_alloc(m->cols, m->rows);
    for (int i = 0; i < m->rows; i++)
        for (int j = 0; j < m->cols; j++)
            t.data[j * m->rows + i] = m->data[i * m->cols + j];
    return t;
}

Matrix matrix_mul(const Matrix *a, const Matrix *b) {
    if (a->cols != b->rows) return matrix_zeros(0, 0);
    Matrix c = matrix_zeros(a->rows, b->cols);
    for (int i = 0; i < a->rows; i++)
        for (int k = 0; k < a->cols; k++) {
            double aik = a->data[i * a->cols + k];
            for (int j = 0; j < b->cols; j++)
                c.data[i * c.cols + j] += aik * b->data[k * b->cols + j];
        }
    return c;
}

Matrix matrix_add(const Matrix *a, const Matrix *b) {
    Matrix c = matrix_alloc(a->rows, a->cols);
    for (int i = 0; i < a->rows * a->cols; i++)
        c.data[i] = a->data[i] + b->data[i];
    return c;
}

Matrix matrix_sub(const Matrix *a, const Matrix *b) {
    Matrix c = matrix_alloc(a->rows, a->cols);
    for (int i = 0; i < a->rows * a->cols; i++)
        c.data[i] = a->data[i] - b->data[i];
    return c;
}

double matrix_frobenius(const Matrix *m) {
    double sum = 0;
    for (int i = 0; i < m->rows * m->cols; i++)
        sum += m->data[i] * m->data[i];
    return sqrt(sum);
}

void matrix_print(const Matrix *m) {
    for (int i = 0; i < m->rows; i++) {
        for (int j = 0; j < m->cols; j++)
            printf("%8.4f ", m->data[i * m->cols + j]);
        printf("\n");
    }
}

/* ── Null space computation via row echelon form ── */

/* Compute rank and null space of an m×n matrix.
 * Returns rank, writes null space basis to null_basis (n × nullity). */
static int compute_null_space(const double *A, int m, int n, double *null_basis, double tol) {
    /* Augmented: work on a copy */
    double *W = malloc((size_t)m * n * sizeof(double));
    memcpy(W, A, (size_t)m * n * sizeof(double));

    int *pivot_col = malloc((size_t)m * sizeof(int));
    int *pivot_row = malloc((size_t)n * sizeof(int));
    for (int i = 0; i < m; i++) pivot_col[i] = -1;
    for (int j = 0; j < n; j++) pivot_row[j] = -1;

    /* Gaussian elimination with partial pivoting */
    int rank = 0;
    for (int col = 0; col < n && rank < m; col++) {
        /* Find pivot */
        int max_row = -1;
        double max_val = tol;
        for (int row = rank; row < m; row++) {
            double val = fabs(W[row * n + col]);
            if (val > max_val) { max_val = val; max_row = row; }
        }
        if (max_row < 0) continue;

        /* Swap rows */
        if (max_row != rank) {
            for (int j = 0; j < n; j++) {
                double tmp = W[rank * n + j];
                W[rank * n + j] = W[max_row * n + j];
                W[max_row * n + j] = tmp;
            }
        }

        pivot_col[rank] = col;
        pivot_row[col] = rank;

        /* Eliminate below */
        double piv = W[rank * n + col];
        for (int row = rank + 1; row < m; row++) {
            double factor = W[row * n + col] / piv;
            W[row * n + col] = 0;
            for (int j = col + 1; j < n; j++)
                W[row * n + j] -= factor * W[rank * n + j];
        }
        rank++;
    }

    /* Back-substitute to get reduced row echelon form */
    for (int i = rank - 1; i >= 0; i--) {
        int pc = pivot_col[i];
        double piv = W[i * n + pc];
        /* Normalize pivot row */
        for (int j = pc; j < n; j++)
            W[i * n + j] /= piv;
        /* Eliminate above */
        for (int row = 0; row < i; row++) {
            double factor = W[row * n + pc];
            for (int j = pc; j < n; j++)
                W[row * n + j] -= factor * W[i * n + j];
        }
    }

    /* Extract null space basis */
    int nullity = n - rank;
    if (null_basis && nullity > 0) {
        /* Free columns are those not in pivot_col */
        bool *is_pivot = calloc((size_t)n, sizeof(bool));
        for (int i = 0; i < rank; i++) is_pivot[pivot_col[i]] = true;

        int free_idx = 0;
        for (int j = 0; j < n; j++) {
            if (!is_pivot[j]) {
                /* Set null space vector: free variable j = 1, others = 0, solve for pivots */
                for (int i = 0; i < n; i++)
                    null_basis[i * nullity + free_idx] = 0.0;

                for (int i = 0; i < rank; i++) {
                    int pc = pivot_col[i];
                    /* x[pc] = -sum over free cols of W[i][free] * x[free] */
                    null_basis[pc * nullity + free_idx] = -W[i * n + j];
                }
                null_basis[j * nullity + free_idx] = 1.0;
                free_idx++;
            }
        }
        free(is_pivot);
    }

    free(W);
    free(pivot_col);
    free(pivot_row);
    return rank;
}

/* ── Find connected components ── */

static int *find_components(int n_verts, int n_edges, const SheafEdge *edges, int *n_components) {
    int *comp = malloc((size_t)n_verts * sizeof(int));
    for (int i = 0; i < n_verts; i++) comp[i] = i;

    for (int e = 0; e < n_edges; e++) {
        int a = edges[e].v1, b = edges[e].v2;
        int ra = a; while (comp[ra] != ra) ra = comp[ra];
        int rb = b; while (comp[rb] != rb) rb = comp[rb];
        if (ra != rb) comp[rb] = ra;
    }
    for (int i = 0; i < n_verts; i++)
        while (comp[i] != comp[comp[i]]) comp[i] = comp[comp[i]];

    int nc = 0;
    bool *seen = calloc((size_t)n_verts, sizeof(bool));
    for (int i = 0; i < n_verts; i++) {
        if (!seen[comp[i]]) { seen[comp[i]] = true; nc++; }
    }
    free(seen);
    *n_components = nc;
    return comp;
}

/* ── Cellular Sheaf ── */

CellularSheaf sheaf_create(int n_verts, const int *stalk_dims, int n_edges) {
    CellularSheaf s;
    s.n_verts = n_verts;
    s.n_edges = n_edges;
    s.stalk_dims = malloc((size_t)n_verts * sizeof(int));
    if (n_verts > 0 && stalk_dims) memcpy(s.stalk_dims, stalk_dims, (size_t)n_verts * sizeof(int));
    s.edges = n_edges ? malloc((size_t)n_edges * sizeof(SheafEdge)) : NULL;
    return s;
}

void sheaf_free(CellularSheaf *s) {
    free(s->stalk_dims);
    for (int i = 0; i < s->n_edges; i++) {
        matrix_free(&s->edges[i].restriction);
        matrix_free(&s->edges[i].restriction_t);
    }
    free(s->edges);
    s->stalk_dims = NULL;
    s->edges = NULL;
}

void sheaf_set_edge(CellularSheaf *s, int idx,
                     int v1, int v2,
                     const Matrix *r1, const Matrix *r2) {
    s->edges[idx].v1 = v1;
    s->edges[idx].v2 = v2;
    s->edges[idx].restriction = matrix_alloc(r1->rows, r1->cols);
    memcpy(s->edges[idx].restriction.data, r1->data,
           (size_t)r1->rows * r1->cols * sizeof(double));
    s->edges[idx].restriction_t = matrix_alloc(r2->rows, r2->cols);
    memcpy(s->edges[idx].restriction_t.data, r2->data,
           (size_t)r2->rows * r2->cols * sizeof(double));
}

Matrix sheaf_laplacian(const CellularSheaf *s) {
    int total_dim = 0;
    for (int v = 0; v < s->n_verts; v++) total_dim += s->stalk_dims[v];

    Matrix L = matrix_zeros(total_dim, total_dim);

    for (int e = 0; e < s->n_edges; e++) {
        SheafEdge *edge = &s->edges[e];
        int v1 = edge->v1, v2 = edge->v2;

        int off1 = 0, off2 = 0;
        for (int v = 0; v < v1; v++) off1 += s->stalk_dims[v];
        for (int v = 0; v < v2; v++) off2 += s->stalk_dims[v];

        int d1 = s->stalk_dims[v1], d2 = s->stalk_dims[v2];

        Matrix r1t = matrix_transpose(&edge->restriction);
        Matrix r1t_r1 = matrix_mul(&r1t, &edge->restriction);

        Matrix r2t = matrix_transpose(&edge->restriction_t);
        Matrix r2t_r2 = matrix_mul(&r2t, &edge->restriction_t);

        Matrix r1t_r2 = matrix_mul(&r1t, &edge->restriction_t);
        Matrix r2t_r1 = matrix_mul(&r2t, &edge->restriction);

        for (int i = 0; i < d1; i++)
            for (int j = 0; j < d1; j++)
                L.data[(off1 + i) * total_dim + (off1 + j)] += r1t_r1.data[i * d1 + j];

        for (int i = 0; i < d2; i++)
            for (int j = 0; j < d2; j++)
                L.data[(off2 + i) * total_dim + (off2 + j)] += r2t_r2.data[i * d2 + j];

        for (int i = 0; i < d1; i++)
            for (int j = 0; j < d2; j++) {
                L.data[(off1 + i) * total_dim + (off2 + j)] -= r1t_r2.data[i * d2 + j];
                L.data[(off2 + j) * total_dim + (off1 + i)] -= r2t_r1.data[j * d1 + i];
            }

        matrix_free(&r1t); matrix_free(&r1t_r1);
        matrix_free(&r2t); matrix_free(&r2t_r2);
        matrix_free(&r1t_r2); matrix_free(&r2t_r1);
    }

    return L;
}

Cohomology sheaf_cohomology(const CellularSheaf *s, double tol) {
    Cohomology c;
    memset(&c, 0, sizeof(c));

    int total_vertex_dim = 0;
    for (int v = 0; v < s->n_verts; v++) total_vertex_dim += s->stalk_dims[v];

    int total_edge_dim = 0;
    for (int e = 0; e < s->n_edges; e++)
        total_edge_dim += s->edges[e].restriction.rows;

    /* Trivial case: empty sheaf */
    if (total_vertex_dim == 0 && total_edge_dim == 0) {
        c.h0_dim = 0;
        c.h1_dim = 0;
        c.h0_basis = matrix_zeros(0, 0);
        c.h1_basis = matrix_zeros(0, 0);
        return c;
    }

    /* No edges: H⁰ = full vertex stalk space, H¹ = 0 */
    if (s->n_edges == 0) {
        c.h0_dim = total_vertex_dim;
        c.h0_basis = (total_vertex_dim > 0) ? matrix_identity(total_vertex_dim)
                                             : matrix_zeros(0, 0);
        c.h1_dim = 0;
        c.h1_basis = matrix_zeros(0, 0);
        return c;
    }

    /* Build coboundary operator d₀: (total_edge_dim × total_vertex_dim)
     * For each edge e=(v1,v2) with restriction maps r1 (v1→edge), r2 (v2→edge):
     *   (d₀ s)_e = r1(s_{v1}) - r2(s_{v2})
     * H⁰ = ker(d₀)  (global sections)
     * H¹ = coker(d₀) = C¹/im(d₀)  (obstructions) */
    double *d0 = calloc((size_t)total_edge_dim * total_vertex_dim, sizeof(double));

    int edge_off = 0;
    for (int e = 0; e < s->n_edges; e++) {
        SheafEdge *edge = &s->edges[e];
        int v1 = edge->v1, v2 = edge->v2;
        int de = edge->restriction.rows;
        int d1 = edge->restriction.cols;
        int d2 = edge->restriction_t.cols;

        int off1 = 0, off2 = 0;
        for (int v = 0; v < v1; v++) off1 += s->stalk_dims[v];
        for (int v = 0; v < v2; v++) off2 += s->stalk_dims[v];

        /* Place +r1 block */
        for (int i = 0; i < de; i++)
            for (int j = 0; j < d1; j++)
                d0[(edge_off + i) * total_vertex_dim + (off1 + j)]
                    += edge->restriction.data[i * d1 + j];

        /* Place -r2 block */
        for (int i = 0; i < de; i++)
            for (int j = 0; j < d2; j++)
                d0[(edge_off + i) * total_vertex_dim + (off2 + j)]
                    -= edge->restriction_t.data[i * d2 + j];

        edge_off += de;
    }

    /* H⁰ = ker(d₀): null space of d₀ */
    double *h0_data = calloc((size_t)total_vertex_dim * total_vertex_dim, sizeof(double));
    int rank_d0 = compute_null_space(d0, total_edge_dim, total_vertex_dim,
                                     h0_data, tol);
    int h0 = total_vertex_dim - rank_d0;

    c.h0_dim = h0;
    if (h0 > 0) {
        c.h0_basis = matrix_alloc(total_vertex_dim, h0);
        memcpy(c.h0_basis.data, h0_data, (size_t)total_vertex_dim * h0 * sizeof(double));
    } else {
        c.h0_basis = matrix_zeros(total_vertex_dim, 0);
    }
    free(h0_data);

    /* H¹ = coker(d₀): null space of d₀^T
     * d₀^T is (total_vertex_dim × total_edge_dim)
     * Elements of ker(d₀^T) are 1-cochains orthogonal to im(d₀) */
    double *d0T = calloc((size_t)total_vertex_dim * total_edge_dim, sizeof(double));
    for (int i = 0; i < total_edge_dim; i++)
        for (int j = 0; j < total_vertex_dim; j++)
            d0T[j * total_edge_dim + i] = d0[i * total_vertex_dim + j];

    double *h1_data = calloc((size_t)total_edge_dim * total_edge_dim, sizeof(double));
    int rank_d0T = compute_null_space(d0T, total_vertex_dim, total_edge_dim,
                                      h1_data, tol);
    int h1 = total_edge_dim - rank_d0T;

    c.h1_dim = h1;
    if (h1 > 0) {
        c.h1_basis = matrix_alloc(total_edge_dim, h1);
        memcpy(c.h1_basis.data, h1_data, (size_t)total_edge_dim * h1 * sizeof(double));
    } else {
        c.h1_basis = matrix_zeros(0, 0);
    }

    free(h1_data);
    free(d0T);
    free(d0);

    return c;
}

void cohomology_free(Cohomology *c) {
    matrix_free(&c->h0_basis);
    matrix_free(&c->h1_basis);
}
