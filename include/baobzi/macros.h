#ifndef MACROS_H
#define MACROS_H

#define BAOBZI_JOIN(A, B) (((A) << 16) | (B))

#define BAOBZI_CASE(DIM, ORDER)                                                                                        \
    case BAOBZI_JOIN(DIM, ORDER): {                                                                                    \
        res.obj = baobzi_init_##DIM##d_##ORDER(fin, center, half_length, tol);                                         \
        res.eval = &baobzi_eval_##DIM##d_##ORDER;                                                                      \
        res.free = &baobzi_free_##DIM##d_##ORDER;                                                                      \
        break;                                                                                                         \
    }

#define BAOBZI_DEFS(DIM, ORDER)                                                                                        \
    double baobzi_eval_##DIM##d_##ORDER(const void *f, const double *x) {                                              \
        return (*(baobzi::Function<DIM, ORDER> *)f).eval(x);                                                           \
    }                                                                                                                  \
    void baobzi_free_##DIM##d_##ORDER(void *f) { delete (baobzi::Function<DIM, ORDER> *)f; }                           \
    void *baobzi_init_##DIM##d_##ORDER(double (*fin)(const double *), const double *center, const double *half_length, \
                                       const double tol) {                                                             \
        return (void *)new baobzi::Function<DIM, ORDER>(fin, center, half_length, tol);                                \
    }

#define BAOBZI_DECLS(DIM, ORDER)                                                                                       \
    double baobzi_eval_##DIM##d_##ORDER(const void *f, const double *x);                                               \
    void baobzi_free_##DIM##d_##ORDER(void *f);                                                                        \
    void *baobzi_init_##DIM##d_##ORDER(double (*)(const double *), const double *center, const double *half_length,    \
                                       const double tol);


#endif
