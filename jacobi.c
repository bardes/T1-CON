#include <pthread.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <fcntl.h>
#include <semaphore.h>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>

#include <sys/stat.h>

#include "utils.h"

typedef struct {
    gsl_matrix *A;        // Coeficientes lineares com a diagonal zerada
    gsl_vector *diag;     // Diagonal dos coeficientes
    gsl_vector *B;        // Vetor de resultados
    gsl_vector *x;        // Soluções do sistema
} jacobi_matrix;

// Matriz constante (pode ser compartilhada entre threads)
typedef struct {
    jacobi_matrix m;
    size_t from;                // Linha inicial (inclusive)
    size_t to;                  // Linha final (não inclusa)
    gsl_vector *result;         // Resultado da estimativa
    gsl_vector *work_area;      // Area de memória usada pra fazer as contas
} thread_data;

static double sum_elements(const gsl_vector *v)
{
    double sum = 0;
    for(size_t i = 0; i < v->size; i++)
        sum += v->data[i];
    return sum;
}

static void mul_elements(gsl_vector *result, const gsl_vector *a,
                         const gsl_vector *b)
                           
{
    if(a->size != b->size) return;
    if(a->size != result->size) return;

    for(size_t i = 0; i < a->size; ++i)
        result->data[i] = a->data[i] * b->data[i];
}

static double compute_line(size_t line, const jacobi_matrix *m,
                           gsl_vector *work_area)
{
    gsl_vector_const_view l = gsl_matrix_const_row(m->A, line);
    mul_elements(work_area, &l.vector, m->x);
    return (gsl_vector_get(m->B, line) - sum_elements(work_area)) /
            gsl_vector_get(m->diag, line); 
}

static double check_error(size_t line, const jacobi_matrix *m,
                          gsl_vector *work_area)
{
    gsl_vector_const_view l = gsl_matrix_const_row(m->A, line);
    mul_elements(work_area, &l.vector, m->x);
    return gsl_vector_get(m->B, line) - sum_elements(work_area)
           - gsl_vector_get(m->diag, line) * gsl_vector_get(m->x, line);
}

static void *compute_block(void *tdata)
{
    thread_data *d = tdata;
    for(size_t line = d->from; line < d->to; ++line)
        d->result->data[line] = compute_line(line, &d->m, d->work_area);
    pthread_exit(NULL);
}

static void read_jacobi_matrix(FILE *input, size_t order, jacobi_matrix *dest)
{
    // Cria a matriz e lê do arquivo de entrada
    dest->A = gsl_matrix_alloc(order, order);
    FATAL(dest->A, EXIT_FAILURE);
    gsl_vector_view diag_view = gsl_matrix_diagonal(dest->A);
    gsl_matrix_fscanf(input, dest->A);

    // Aloca os vetores
    dest->diag = gsl_vector_alloc(order);
    dest->B    = gsl_vector_alloc(order);
    dest->x    = gsl_vector_alloc(order);
    FATAL(dest->diag && dest->B && dest->x, EXIT_FAILURE);

    // Copia a diagonal para o vetor e zera a diagonal da matriz
    gsl_vector_memcpy(dest->diag, &diag_view.vector);
    gsl_vector_set_zero(&diag_view.vector);

    // Lê o vetor de resultados 
    gsl_vector_fscanf(stdin, dest->B);

    // Começa com a primeira estimativa de solução [0, 0, 0, ..., 0]
    gsl_vector_set_zero(dest->x);

}

static void free_jacobi_matrix(jacobi_matrix *m)
{
    gsl_matrix_free(m->A);
    gsl_vector_free(m->diag);
    gsl_vector_free(m->B);
    gsl_vector_free(m->x);
}

int main(int argc, char *argv[])
{
    // Determina a qtd. de threads a serem executadas
    unsigned long n_threads;
    if(argc == 1) {
        n_threads = 1;
    } else if (argc == 2) {
        n_threads = strtoul(argv[1], NULL, 10);
        if(errno == EINVAL || errno == ERANGE)
            exit(EXIT_FAILURE);
    } else {
        fprintf(stderr, "Usage: %s [n_threads=1]\n", argv[0]);
        exit(EXIT_FAILURE);
    }

    // Lê os "metadaos"
    double max_err;
    size_t order, test_row, max_iter;
    scanf("%zu %zu %lf %zu", &order, &test_row, &max_err, &max_iter);

    // Lê a matriz
    jacobi_matrix m;
    read_jacobi_matrix(stdin, order, &m);
    
    // Array de threads
    pthread_t *threads = malloc(sizeof(pthread_t) * n_threads);
    FATAL(threads, EXIT_FAILURE);

    // Vetor para guardar os resultados parciais de cada iteração
    gsl_vector *tmp = gsl_vector_alloc(order);
    FATAL(tmp, EXIT_FAILURE);


    // Calculando a divisão de linhas por threads
    size_t step = order / n_threads;
    thread_data *td = malloc(sizeof(thread_data) * n_threads);
    FATAL(td, EXIT_FAILURE);
    for(size_t i = 0; i < n_threads; ++i) {
        td[i].m = m;
        td[i].from = i * step;
        td[i].to = (i+1) * step;
        td[i].result = tmp;
        td[i].work_area = gsl_vector_alloc(order);
        FATAL(td[i].work_area, EXIT_FAILURE);
    }

    td[n_threads - 1].to += order % n_threads;
    
    double err = 0; size_t itrs;
    for(itrs = 0; itrs < max_iter; ++itrs) {

        // Dispara cada thread
        for(size_t i = 0; i < n_threads; ++i)
            if(pthread_create(threads + i, NULL, compute_block, td + i))
                FATAL_MSG(0, EXIT_FAILURE, "Failed to create thread!");

        // Espera os resultados
        for(size_t i = 0; i < n_threads; ++i)
            pthread_join(threads[i], NULL);
        
        // Copia para o vetor de resultados
        gsl_vector_memcpy(m.x, tmp);

        // Verifia se o erro esta dentro do limite
        err = check_error(test_row, &td[0].m, tmp);
        if(fabs(err) < max_err) break;
    }

    printf("Iterações: %zu, Erro: %.3g\n", itrs, err);

    free_jacobi_matrix(&m);
    gsl_vector_free(tmp);
    for(size_t i = 0; i < n_threads; ++i)
        gsl_vector_free(td[i].work_area);
    free(td);
    free(threads);
    return 0;
}
