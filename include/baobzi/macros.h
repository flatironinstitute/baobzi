#ifndef MACROS_H
#define MACROS_H

#define BAOBZI_JOIN(A, B, C) (((A) << 16) | (B) | ((C) << 24))

#define BAOBZI_CASE(DIM, ORDER, ISET)                                                                                  \
    case BAOBZI_JOIN(DIM, ORDER, ISET): {                                                                              \
        res->obj = baobzi_init_##DIM##d_##ORDER##_##ISET(input, center, half_length);                                  \
        res->save = &baobzi_save_##DIM##d_##ORDER##_##ISET;                                                            \
        res->eval = &baobzi_eval_##DIM##d_##ORDER##_##ISET;                                                            \
        res->eval_multi = &baobzi_eval_multi_##DIM##d_##ORDER##_##ISET;                                                \
        res->stats = &baobzi_stats_##DIM##d_##ORDER##_##ISET;                                                          \
        res->free = &baobzi_free_##DIM##d_##ORDER##_##ISET;                                                            \
        break;                                                                                                         \
    }

#define BAOBZI_CASE_RESTORE(DIM, ORDER, ISET)                                                                          \
    case BAOBZI_JOIN(DIM, ORDER, ISET): {                                                                              \
        res->obj = baobzi_restore_##DIM##d_##ORDER##_##ISET(&obj);                                                     \
        res->save = &baobzi_save_##DIM##d_##ORDER##_##ISET;                                                            \
        res->eval = &baobzi_eval_##DIM##d_##ORDER##_##ISET;                                                            \
        res->eval_multi = &baobzi_eval_multi_##DIM##d_##ORDER##_##ISET;                                                \
        res->stats = &baobzi_stats_##DIM##d_##ORDER##_##ISET;                                                          \
        res->free = &baobzi_free_##DIM##d_##ORDER##_##ISET;                                                            \
        break;                                                                                                         \
    }

#define BAOBZI_DEFS(DIM, ORDER, ISET)                                                                                  \
    void baobzi_eval_##DIM##d_##ORDER##_##ISET(const void *f, const double *x, double *y) {                            \
        (*(baobzi::Function<DIM, ORDER, ISET> *)f).eval(x, y);                                                         \
    }                                                                                                                  \
    void baobzi_eval_multi_##DIM##d_##ORDER##_##ISET(const void *f, const double *x, double *res, int ntrg) {          \
        (*(baobzi::Function<DIM, ORDER, ISET> *)f).eval(x, res, ntrg);                                                 \
    }                                                                                                                  \
    void baobzi_save_##DIM##d_##ORDER##_##ISET(const void *f, const char *filename) {                                  \
        (*(baobzi::Function<DIM, ORDER, ISET> *)f).save(filename);                                                     \
    }                                                                                                                  \
    void baobzi_stats_##DIM##d_##ORDER##_##ISET(void *f) { (*(baobzi::Function<DIM, ORDER, ISET> *)f).print_stats(); } \
    void baobzi_free_##DIM##d_##ORDER##_##ISET(void *f) { delete (baobzi::Function<DIM, ORDER, ISET> *)f; }            \
    void *baobzi_init_##DIM##d_##ORDER##_##ISET(const baobzi_input_t *input, const double *center,                     \
                                                const double *half_length) {                                           \
        return (void *)new (std::nothrow) baobzi::Function<DIM, ORDER, ISET>(input, center, half_length);              \
    }                                                                                                                  \
    void *baobzi_restore_##DIM##d_##ORDER##_##ISET(const void *obj_) {                                                 \
        const msgpack::object *obj = (msgpack::object *)obj_;                                                          \
        baobzi::Function<DIM, ORDER, ISET> *f_ptr = new baobzi::Function<DIM, ORDER, ISET>();                          \
        *f_ptr = obj->as<baobzi::Function<DIM, ORDER, ISET>>();                                                        \
        f_ptr->build_cache();                                                                                          \
        return f_ptr;                                                                                                  \
    }

#define BAOBZI_DECLS(DIM, ORDER, ISET)                                                                                 \
    void baobzi_eval_##DIM##d_##ORDER##_##ISET(const void *f, const double *x, double *y);                             \
    void baobzi_eval_multi_##DIM##d_##ORDER##_##ISET(const void *f, const double *x, double *res, int ntrg);           \
    void baobzi_save_##DIM##d_##ORDER##_##ISET(const void *f, const char *filename);                                   \
    void baobzi_stats_##DIM##d_##ORDER##_##ISET(void *f);                                                              \
    void baobzi_free_##DIM##d_##ORDER##_##ISET(void *f);                                                               \
    void *baobzi_init_##DIM##d_##ORDER##_##ISET(const baobzi_input_t *input, const double *center,                     \
                                                const double *half_length);                                            \
    void *baobzi_restore_##DIM##d_##ORDER##_##ISET(const void *obj);

#endif
