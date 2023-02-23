#include "baobzi.h"
#include "class_handle.hpp"
#include "mex.h"

#include <cstdint>
#include <list>
#include <mutex>
#include <set>
#include <string>

struct wrapper_data_t {
    std::string fit_fcn;
    int dim;
    int order;
};

double matfun_wrapper(const double *x, const void *data_) {
    mxArray *prhs[1];
    mxArray *plhs[1];
    const wrapper_data_t *data = (wrapper_data_t *)data_;
    plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
    prhs[0] = mxCreateDoubleMatrix(data->dim, 1, mxREAL);

    double *xptr = mxGetPr(prhs[0]);
    for (int i = 0; i < data->dim; ++i)
        xptr[i] = x[i];

    int status = mexCallMATLAB(1, plhs, 1, prhs, data->fit_fcn.c_str());
    if (status != 0) {
        mexErrMsgIdAndTxt("MATLAB:mexfeval:mxCallMATLAB", "Failed to execute MATLAB function.");
        return 0.0;
    }

    return mxGetPr(plhs[0])[0];
}

class baobzi {
  public:
    baobzi(const std::string &matfun, int dim, int order, const double *center, const double *half_length,
           const double tol, const double minimum_leaf_fraction, const int split_multi_eval, const int max_depth)
        : data_({matfun, dim, order}) {

        baobzi_input_t input = {.func = matfun_wrapper,
                                .data = &data_,
                                .dim = dim,
                                .order = order,
                                .tol = tol,
                                .minimum_leaf_fraction = minimum_leaf_fraction,
                                .split_multi_eval = split_multi_eval,
                                .max_depth = max_depth};

        // to prevent matlab from segfaulting if your function doesn't exist, try to call it once
        matfun_wrapper(center, &data_);
        obj_ = baobzi_init(&input, center, half_length);
    }
    baobzi(const char *infile) {
        obj_ = baobzi_restore(infile);
        data_.dim = obj_->DIM;
        data_.order = obj_->ORDER;
    }

    ~baobzi() {
        if (obj_)
            obj_ = baobzi_free(obj_);
    }
    double eval(const double *x) { return baobzi_eval(obj_, x); };
    void eval(const double *x, double *res, int ntrg) { baobzi_eval_multi(obj_, x, res, ntrg); };
    void save(const std::string &fname) { baobzi_save(obj_, fname.c_str()); }
    void stats() const { baobzi_stats(obj_); }

    wrapper_data_t data_;
    baobzi_t obj_ = nullptr;
};

std::string to_string(const mxArray *arr) {
    size_t buflen = mxGetN(arr) + 1;
    std::string res(buflen, '\0');
    int status = mxGetString(arr, &res[0], buflen);

    if (status != 0)
        mexErrMsgIdAndTxt("MATLAB:mexfeval:mxGetString", "Failed to copy string.");

    return res;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    // Get the command string
    char cmd[64];
    if (nrhs < 1 || mxGetString(prhs[0], cmd, sizeof(cmd)))
        mexErrMsgTxt("First input should be a command string less than 64 characters long.");

    // New
    if (!strcmp("new", cmd)) {
        // Check parameters
        if (nlhs != 1)
            mexErrMsgTxt("New: One output expected.");
        if (nrhs != 10)
            mexErrMsgTxt("New: Ten inputs expected.");

        std::string baobzi_fit_fcn = to_string(prhs[1]);

        int dim = mxGetPr(prhs[2])[0];
        int order = mxGetPr(prhs[3])[0];
        double *center = mxGetPr(prhs[4]);
        double *half_length = mxGetPr(prhs[5]);
        double tol = mxGetPr(prhs[6])[0];
        double minimum_leaf_fraction = mxGetPr(prhs[7])[0];
        int split_multi_eval = mxGetPr(prhs[8])[0];
        int max_depth = mxGetPr(prhs[9])[0];

        // Return a handle to a new C++ instance
        plhs[0] = convertPtr2Mat<baobzi>(new baobzi(baobzi_fit_fcn, dim, order, center, half_length, tol,
                                                    minimum_leaf_fraction, split_multi_eval, max_depth));
        return;
    }

    if (!strcmp("restore", cmd)) {
        // Check parameters
        if (nlhs != 1)
            mexErrMsgTxt("New: One output expected.");
        if (nrhs != 2)
            mexErrMsgTxt("New: One input expected.");

        std::string infile = to_string(prhs[1]);

        // Return a handle to a new C++ instance
        plhs[0] = convertPtr2Mat<baobzi>(new baobzi(infile.c_str()));

        return;
    }

    // Check there is a second input, which should be the class instance handle
    if (nrhs < 2)
        mexErrMsgTxt("Second input should be a class instance handle.");

    // Delete
    if (!strcmp("free", cmd)) {
        // Destroy the C++ object
        destroyObject<baobzi>(prhs[1]);
        // Warn if other commands were ignored
        if (nlhs != 0 || nrhs != 2)
            mexWarnMsgTxt("Delete: Unexpected arguments ignored.");
        return;
    }

    // Get the class instance pointer from the second input
    baobzi *baobzi_instance = convertMat2Ptr<baobzi>(prhs[1]);

    if (!strcmp("stats", cmd)) {
        baobzi_instance->stats();
        return;
    }

    if (!strcmp("eval", cmd)) {
        // Check parameters
        if (nlhs > 1 || nrhs != 3)
            mexErrMsgTxt("eval: Unexpected arguments.");
        // Call the method
        int n_points = mxGetNumberOfElements(prhs[2]) / baobzi_instance->data_.dim;
        plhs[0] = mxCreateDoubleMatrix(n_points, 1, mxREAL);
        if (n_points == 1)
            mxGetPr(plhs[0])[0] = baobzi_instance->eval(mxGetPr(prhs[2]));
        else
            baobzi_instance->eval(mxGetPr(prhs[2]), mxGetPr(plhs[0]), n_points);

        return;
    }

    if (!strcmp("save", cmd)) {
        if (nlhs > 0 || nrhs != 3)
            mexErrMsgTxt("save: Unexpected arguments.");

        std::string outfile = to_string(prhs[2]);

        baobzi_instance->save(outfile);
        return;
    }

    // Got here, so command not recognized
    mexErrMsgTxt("Command not recognized.");
}
