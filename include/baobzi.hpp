#ifndef BAOBZI_HPP
#define BAOBZI_HPP

#include <baobzi.h>
#include <string>

namespace baobzi {
class Baobzi {
private:
  baobzi_t obj_ = nullptr;
public:
  Baobzi(double (*fin)(const double *), uint16_t dim, uint16_t order, const double *center, const double *half_length,
         const double tol)
      : obj_(baobzi_init(fin, dim, order, center, half_length, tol)) {}
  Baobzi(const std::string &input_file) : obj_(baobzi_restore(input_file.c_str())) {}
  ~Baobzi() { obj_ = baobzi_free(obj_); }

  void save(const std::string &output_file) { baobzi_save(obj_, output_file.c_str()); }
  double operator()(const double *x) const { return baobzi_eval(obj_, x); }
};
} // namespace baobzi
#endif
