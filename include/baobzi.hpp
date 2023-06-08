#ifndef BAOBZI_HPP
#define BAOBZI_HPP

#include <baobzi.h>
#include <string>
#include <memory>

namespace baobzi {
/// Wrapper class for C library (which is a wrapper for the template library. oof)
class Baobzi {
  private:
    std::shared_ptr<baobzi_struct> obj_ = nullptr; ///< Pointer to C baobzi struct this class wraps
  public:
    /// @brief Construct Baobzi object from input function
    /// @param[in] input pointer to baobzi_input_t object
    /// @param[in] center [dim] center of the domain
    /// @param[in] half_length [dim] half the size of the domain in each dimension
    Baobzi(const baobzi_input_t *input, const double *center, const double *half_length)
        : obj_(baobzi_init(input, center, half_length), baobzi_free) {}

    /// @brief Restore baobzi object from serialized version
    /// @param[in] input_file path to file of serialized function
    Baobzi(const std::string &input_file) : obj_(baobzi_restore(input_file.c_str()), baobzi_free) {}

    /// @brief Save Baobzi object to file
    /// @param[in] output_file path of file to serialize object to
    inline void save(const std::string &output_file) { baobzi_save(obj_.get(), output_file.c_str()); }

    /// @brief Print various stats about Baobzi object
    inline void stats(const std::string &output_file) { baobzi_stats(obj_.get()); }

    /// @brief Evaluate Baobzi object at point
    /// @param[in] x [dim]
    /// @return approximate value of function at point x
    inline void operator()(const double *x, double *y) const { baobzi_eval(obj_.get(), x, y); }

    /// @brief eval function approximation at n_trg points
    /// @param[in] xp [DIM * n_trg] array of points to evaluate function at
    /// @param[out] res [DIM * n_trg] array of results
    /// @param[in] n_trg number of points to evaluate
    inline void operator()(const double *xp, double *res, int n_trg) const {
        baobzi_eval_multi(obj_.get(), xp, res, n_trg);
    }
};
} // namespace baobzi
#endif
