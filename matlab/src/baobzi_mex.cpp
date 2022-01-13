#include "baobzi.h"
#include "class_handle.hpp"
#include "mex.h"

#include <string>

static char baobzi_fit_fcn[64];
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

    int status = mexCallMATLAB(1, plhs, 1, prhs, baobzi_fit_fcn);
    return mxGetPr(plhs[0])[0];
}

// The class that we are interfacing to
class baobzi {
  public:
    baobzi(const char *matfun, int dim, int order, const double *center, const double *half_length, const double tol)
        : funname_(matfun), dim_(dim), order_(order) {
        set_all_params();
        obj_ = baobzi_init(matfun_wrapper, dim, order, center, half_length, tol);
    }
    ~baobzi() { obj_ = baobzi_free(obj_); }
    double eval(const double *x) {
        set_other_params();
        return baobzi_eval(obj_, x);
    };
    void test() { mexPrintf("Calling test\n"); };

    void set_all_params() {
        sprintf(baobzi_fit_fcn, funname_.c_str());
        set_other_params();
    }

    void set_other_params() {
        baobzi_dim = dim_;
        baobzi_order = order_;
    }

  private:
    baobzi_t obj_;
    int dim_;
    int order_;
    std::string funname_;
};

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    // Get the command string
    char cmd[64];
    if (nrhs < 1 || mxGetString(prhs[0], cmd, sizeof(cmd)))
        mexErrMsgTxt("First input should be a command string less than 64 characters long.");

    // New
    if (!strcmp("new", cmd)) {
        // Check parameters
        mexPrintf("%d\n", nrhs);
        if (nlhs != 1)
            mexErrMsgTxt("New: One output expected.");
        if (nrhs != 2)
            mexErrMsgTxt("New: One input expected.");

        size_t buflen = mxGetN(prhs[1]) + 1;
        int status = mxGetString(prhs[1], baobzi_fit_fcn, (mwSize)buflen);
        mexPrintf("%s\n", baobzi_fit_fcn);

        if (status != 0)
            mexErrMsgIdAndTxt("MATLAB:mexfeval:mxGetString", "Failed to copy function string into allocated memory.");

        int dim = 2;
        int order = 6;
        double center[2] = {0.0, 0.0};
        double half_length[2] = {1.0, 1.0};
        double tol = 1E-8;

        // Return a handle to a new C++ instance
        plhs[0] = convertPtr2Mat<baobzi>(new baobzi(baobzi_fit_fcn, dim, order, center, half_length, tol));
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

    // Got here, so command not recognized
    mexErrMsgTxt("Command not recognized.");
}
