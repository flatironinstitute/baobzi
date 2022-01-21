#ifndef BAOBZI_HPP
#define BAOBZI_HPP

#include <baobzi.h>
#include <string>

namespace baobzi {
/// Wrapper class for C library (which is a wrapper for the template library. oof)
class Baobzi {
  private:
    baobzi_t obj_ = nullptr; ///< Pointer to C baobzi struct this class wraps
  public:
    /// @brief Construct Baobzi object from input function
    /// @param[in] fin pointer to function to fit
    /// @param[in] dim input dimension of the function
    /// @param[in] order order of polynomial to fit
    /// @param[in] center [dim] center of the domain
    /// @param[in] half_length [dim] half the size of the domain in each dimension
    /// @param[in] tol desired relative tolerance
    Baobzi(double (*fin)(const double *), uint16_t dim, uint16_t order, const double *center, const double *half_length,
           const double tol)
        : obj_(baobzi_init(fin, dim, order, center, half_length, tol)) {}

    /// @brief Restore baobzi object from serialized version
    /// @param[in] input_file path to file of serialized function
    Baobzi(const std::string &input_file) : obj_(baobzi_restore(input_file.c_str())) {}

    /// @brief Destroy Baobzi object and associated memory
    ~Baobzi() { obj_ = baobzi_free(obj_); }

    /// @brief Save Baobzi object to file
    /// @param[in] output_file path of file to serialize object to
    void save(const std::string &output_file) { baobzi_save(obj_, output_file.c_str()); }

    /// @brief Evaluate Baobzi object at point
    /// @param[in] x [dim]
    /// @return approximate value of function at point x
    double operator()(const double *x) const { return baobzi_eval(obj_, x); }
};
} // namespace baobzi
#endif
