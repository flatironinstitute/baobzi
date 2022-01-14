#include "baobzi.h"
#include "class_handle.hpp"
#include "mex.h"

#include <string>

static std::string baobzi_fit_fcn;
static uint16_t baobzi_dim;
static uint16_t baobzi_order;

double matfun_wrapper(const double *x) {
    mxArray *prhs[1];
    mxArray *plhs[1];
    plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
    prhs[0] = mxCreateDoubleMatrix(baobzi_dim, 1, mxREAL);

    double *xptr = mxGetPr(prhs[0]);
    for (int i = 0; i < baobzi_dim; ++i)
        xptr[i] = x[i];

    int status = mexCallMATLAB(1, plhs, 1, prhs, baobzi_fit_fcn.c_str());
    if (status != 0)
        mexErrMsgIdAndTxt("MATLAB:mexfeval:mxCallMATLAB", "Failed to execute MATLAB function.");

    return mxGetPr(plhs[0])[0];
}

class baobzi {
  public:
    baobzi(const std::string &matfun, int dim, int order, const double *center, const double *half_length,
           const double tol)
        : funname_(matfun), dim_(dim), order_(order) {
        set_all_params();
        obj_ = baobzi_init(matfun_wrapper, dim_, order_, center, half_length, tol);
    }
    baobzi(const std::string &matfun, const char *infile) : funname_(matfun) {
        baobzi_header_t header = baobzi_read_header_from_file(infile);
        dim_ = header.dim;
        order_ = header.order;
        set_all_params();
        obj_ = baobzi_restore(matfun_wrapper, infile);
    }

    ~baobzi() { obj_ = baobzi_free(obj_); }
    double eval(const double *x) {
        set_other_params();
        return baobzi_eval(obj_, x);
    };
    void save(const std::string &fname) { baobzi_save(obj_, fname.c_str()); }

    void set_all_params() {
        baobzi_fit_fcn = funname_;
        set_other_params();
    }

    void set_other_params() {
        baobzi_dim = dim_;
        baobzi_order = order_;
    }

  private:
    std::string funname_;
    int dim_;
    int order_;
    baobzi_t obj_;
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
        if (nrhs != 7)
            mexErrMsgTxt("New: Seven inputs expected.");

        baobzi_fit_fcn = to_string(prhs[1]);

        int dim = mxGetPr(prhs[2])[0];
        int order = mxGetPr(prhs[3])[0];
        double *center = mxGetPr(prhs[4]);
        double *half_length = mxGetPr(prhs[5]);
        double tol = mxGetPr(prhs[6])[0];

        // Return a handle to a new C++ instance
        plhs[0] = convertPtr2Mat<baobzi>(new baobzi(baobzi_fit_fcn, dim, order, center, half_length, tol));
        return;
    }

    if (!strcmp("restore", cmd)) {
        // Check parameters
        if (nlhs != 1)
            mexErrMsgTxt("New: One output expected.");
        if (nrhs != 3)
            mexErrMsgTxt("New: Two inputs expected.");

        baobzi_fit_fcn = to_string(prhs[1]);
        std::string infile = to_string(prhs[2]);

        // Return a handle to a new C++ instance
        plhs[0] = convertPtr2Mat<baobzi>(new baobzi(baobzi_fit_fcn, infile.c_str()));

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

    if (!strcmp("eval", cmd)) {
        // Check parameters
        if (nlhs > 1 || nrhs != 3)
            mexErrMsgTxt("eval: Unexpected arguments.");
        // Call the method
        plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
        mxGetPr(plhs[0])[0] = baobzi_instance->eval(mxGetPr(prhs[2]));

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
