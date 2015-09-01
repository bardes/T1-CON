#include <stdio.h>
#include <math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>

// Matriz constante (pode ser compartilhada entre threads)
typedef struct {
    const gsl_matrix *A;        // Coeficientes lineares com a diagonal zerada
    const gsl_vector *diag;     // Diagonal dos coeficientes
    const gsl_vector *B;        // Vetor de resultados
    const gsl_vector *x;        // Soluções do sistema
} jacobi_matrix;

auto double sum_elements(const gsl_vector *v)
{
    double sum = 0;
    for(size_t i = 0; i < v->size; i++)
        sum += v->data[i];
    return sum;
}

double compute_line(size_t line, const jacobi_matrix *m,
                    gsl_vector *work_area)
{
    gsl_matrix_get_row(work_area, m->A, line);
    gsl_vector_mul(work_area, m->x);
    return (gsl_vector_get(m->B, line) - sum_elements(work_area)) /
            gsl_vector_get(m->diag, line); 
}

double check_error(size_t line, const jacobi_matrix *m,
                   gsl_vector *work_area)
{
    gsl_matrix_get_row(work_area, m->A, line);
    gsl_vector_mul(work_area, m->x);
    return gsl_vector_get(m->B, line) - sum_elements(work_area)
           - gsl_vector_get(m->diag, line) * gsl_vector_get(m->x, line);
}

int main()
{
    // Lê os "metadaos"
    double max_err;
    size_t order, test_row, max_iter;
    scanf("%zu %zu %lf %zu", &order, &test_row, &max_err, &max_iter);

    // Cria a matriz e lê do arquivo de entrada
    gsl_matrix *A = gsl_matrix_alloc(order, order);
    gsl_vector_view diag_view = gsl_matrix_diagonal(A);
    gsl_matrix_fscanf(stdin, A);

    // Aloca os vetores
    gsl_vector *diag = gsl_vector_alloc(order);
    gsl_vector    *B = gsl_vector_alloc(order);
    gsl_vector    *x = gsl_vector_alloc(order);

    // Copia a diagonal para o vetor e zera a diagonal da matriz
    gsl_vector_memcpy(diag, &diag_view.vector);
    gsl_vector_set_zero(&diag_view.vector);

    // Lê o vetor de resultados 
    gsl_vector_fscanf(stdin, B);

    // Começa com a primeira estimativa de solução [0, 0, 0, ..., 0]
    gsl_vector_set_zero(x);
    
    // Struct usada como entrada para iterar sobre uma linha.
    jacobi_matrix m = {
        .A = A,
        .diag = diag,
        .B = B,
        .x = x
    };

    gsl_vector *tmp = gsl_vector_alloc(order);
    gsl_vector *wa = gsl_vector_alloc(order);
    
    double err;
    size_t itrs;
    for(itrs = 0; itrs < max_iter; ++itrs) {
        for(size_t line = 0; line < order; ++line) {
            gsl_vector_set(tmp, line, compute_line(line, &m, wa));
        }
        gsl_vector_memcpy(x, tmp);

        err = check_error(test_row, &m, wa);
        if(fabs(err) < max_err) break;
    }

    printf("Iterações: %zu, Erro: %.3g\n", itrs, err);

    gsl_matrix_free(A);
    gsl_vector_free(diag);
    gsl_vector_free(B);
    gsl_vector_free(x);
    gsl_vector_free(tmp);
    return 0;
}
