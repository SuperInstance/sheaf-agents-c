#include "sheaf_agents.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

/* ── Agent Network ── */

AgentNetwork agent_network_create(const CellularSheaf *s) {
    AgentNetwork net;
    net.sheaf = sheaf_create(s->n_verts, s->stalk_dims, s->n_edges);
    for (int e = 0; e < s->n_edges; e++) {
        net.sheaf.edges[e].v1 = s->edges[e].v1;
        net.sheaf.edges[e].v2 = s->edges[e].v2;
        net.sheaf.edges[e].restriction = matrix_alloc(s->edges[e].restriction.rows, s->edges[e].restriction.cols);
        memcpy(net.sheaf.edges[e].restriction.data, s->edges[e].restriction.data,
               (size_t)s->edges[e].restriction.rows * s->edges[e].restriction.cols * sizeof(double));
        net.sheaf.edges[e].restriction_t = matrix_alloc(s->edges[e].restriction_t.rows, s->edges[e].restriction_t.cols);
        memcpy(net.sheaf.edges[e].restriction_t.data, s->edges[e].restriction_t.data,
               (size_t)s->edges[e].restriction_t.rows * s->edges[e].restriction_t.cols * sizeof(double));
    }
    net.n_agents = s->n_verts;
    net.agents = malloc((size_t)s->n_verts * sizeof(AgentState));
    for (int v = 0; v < s->n_verts; v++) {
        net.agents[v].vertex = v;
        net.agents[v].belief = matrix_zeros(s->stalk_dims[v], 1);
    }
    return net;
}

void agent_network_free(AgentNetwork *net) {
    for (int i = 0; i < net->n_agents; i++)
        matrix_free(&net->agents[i].belief);
    free(net->agents);
    net->agents = NULL;
    sheaf_free(&net->sheaf);
}

void agent_set_belief(AgentNetwork *net, int agent, const Matrix *belief) {
    int d = net->sheaf.stalk_dims[agent];
    memcpy(net->agents[agent].belief.data, belief->data, (size_t)d * sizeof(double));
}

double agent_disagreement(const AgentNetwork *net) {
    double total = 0;
    for (int e = 0; e < net->sheaf.n_edges; e++) {
        SheafEdge *edge = &net->sheaf.edges[e];
        Matrix x1 = net->agents[edge->v1].belief;
        Matrix x2 = net->agents[edge->v2].belief;

        Matrix r1x = matrix_mul(&edge->restriction, &x1);
        Matrix r2x = matrix_mul(&edge->restriction_t, &x2);
        Matrix diff = matrix_sub(&r1x, &r2x);
        total += matrix_frobenius(&diff);

        matrix_free(&r1x);
        matrix_free(&r2x);
        matrix_free(&diff);
    }
    return total;
}

Matrix agent_synchronize(const AgentNetwork *net, double step, int *converged) {
    int total_dim = 0;
    for (int v = 0; v < net->sheaf.n_verts; v++)
        total_dim += net->sheaf.stalk_dims[v];

    double *gradient = calloc((size_t)total_dim, sizeof(double));

    for (int e = 0; e < net->sheaf.n_edges; e++) {
        SheafEdge *edge = &net->sheaf.edges[e];
        int v1 = edge->v1, v2 = edge->v2;
        int off1 = 0, off2 = 0;
        for (int v = 0; v < v1; v++) off1 += net->sheaf.stalk_dims[v];
        for (int v = 0; v < v2; v++) off2 += net->sheaf.stalk_dims[v];

        Matrix x1 = net->agents[v1].belief;
        Matrix x2 = net->agents[v2].belief;

        Matrix r1x = matrix_mul(&edge->restriction, &x1);
        Matrix r2x = matrix_mul(&edge->restriction_t, &x2);
        Matrix diff = matrix_sub(&r1x, &r2x);

        Matrix r1t = matrix_transpose(&edge->restriction);
        Matrix g1 = matrix_mul(&r1t, &diff);
        for (int i = 0; i < g1.rows; i++)
            gradient[off1 + i] += g1.data[i];

        Matrix r2t = matrix_transpose(&edge->restriction_t);
        Matrix g2 = matrix_mul(&r2t, &diff);
        for (int i = 0; i < g2.rows; i++)
            gradient[off2 + i] -= g2.data[i];

        matrix_free(&r1x); matrix_free(&r2x); matrix_free(&diff);
        matrix_free(&r1t); matrix_free(&g1);
        matrix_free(&r2t); matrix_free(&g2);
    }

    double grad_norm = 0;
    for (int i = 0; i < total_dim; i++) grad_norm += gradient[i] * gradient[i];
    grad_norm = sqrt(grad_norm);

    if (converged) *converged = (grad_norm < 1e-10);

    /* Apply gradient descent to a mutable copy */
    AgentNetwork *mut_net = (AgentNetwork *)net;
    int off = 0;
    for (int v = 0; v < mut_net->sheaf.n_verts; v++) {
        int d = mut_net->sheaf.stalk_dims[v];
        for (int i = 0; i < d; i++)
            mut_net->agents[v].belief.data[i] -= step * gradient[off + i];
        off += d;
    }

    Matrix combined = matrix_zeros(total_dim, 1);
    off = 0;
    for (int v = 0; v < net->sheaf.n_verts; v++) {
        int d = net->sheaf.stalk_dims[v];
        for (int i = 0; i < d; i++)
            combined.data[off + i] = net->agents[v].belief.data[i];
        off += d;
    }

    free(gradient);
    return combined;
}

ConvergenceResult agent_converge(AgentNetwork *net, double step, int max_iter, double tol) {
    ConvergenceResult res = {0, 0, false, false};

    Cohomology coh = sheaf_cohomology(&net->sheaf, tol);
    bool has_obstruction = (coh.h1_dim > 0);
    cohomology_free(&coh);

    for (int i = 0; i < max_iter; i++) {
        double dis = agent_disagreement(net);
        res.iterations = i + 1;
        res.final_disagreement = dis;

        if (dis < tol) {
            res.converged = true;
            return res;
        }

        int conv;
        Matrix x = agent_synchronize(net, step, &conv);
        matrix_free(&x);

        if (conv) {
            res.converged = true;
            return res;
        }
    }

    res.obstruction_detected = has_obstruction;
    res.final_disagreement = agent_disagreement(net);
    return res;
}

/* ── Consensus ── */

bool can_agree(const CellularSheaf *s, double tol) {
    Cohomology c = sheaf_cohomology(s, tol);
    bool result = (c.h1_dim == 0);
    cohomology_free(&c);
    return result;
}

Matrix forced_disagreement(const CellularSheaf *s, double tol) {
    Cohomology c = sheaf_cohomology(s, tol);
    Matrix result;
    if (c.h1_dim > 0) {
        result = matrix_alloc(c.h1_basis.rows, c.h1_basis.cols);
        memcpy(result.data, c.h1_basis.data,
               (size_t)result.rows * result.cols * sizeof(double));
    } else {
        result = matrix_zeros(0, 0);
    }
    cohomology_free(&c);
    return result;
}

ConsensusQuality quality_metric(const AgentNetwork *net, double tol) {
    ConsensusQuality q;
    memset(&q, 0, sizeof(q));

    Cohomology c = sheaf_cohomology(&net->sheaf, tol);
    q.h0_dim = c.h0_dim;
    q.h1_dim = c.h1_dim;

    double dis = agent_disagreement(net);

    int total_dim = 0;
    for (int v = 0; v < net->sheaf.n_verts; v++)
        total_dim += net->sheaf.stalk_dims[v];

    q.agreement_score = 1.0 - (dis / (dis + 1.0));
    q.h0_quality = (total_dim > 0) ? (double)c.h0_dim / total_dim : 1.0;
    q.h1_obstruction = (double)c.h1_dim;

    cohomology_free(&c);
    return q;
}
