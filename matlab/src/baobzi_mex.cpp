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
        obj_ = baobzi_init(matfun_wrapper, dim_, order_, center, half_length, tol);
    }
    baobzi(const char *matfun, const char *infile) : funname_(matfun) {
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
        if (nlhs != 1)
            mexErrMsgTxt("New: One output expected.");
        if (nrhs != 7)
            mexErrMsgTxt("New: Seven inputs expected.");

        size_t buflen = mxGetN(prhs[1]) + 1;
        int status = mxGetString(prhs[1], baobzi_fit_fcn, (mwSize)buflen);

        if (status != 0)
            mexErrMsgIdAndTxt("MATLAB:mexfeval:mxGetString", "Failed to copy function string into allocated memory.");

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

        size_t buflen = mxGetN(prhs[1]) + 1;
        int status = mxGetString(prhs[1], baobzi_fit_fcn, (mwSize)buflen);

        if (status != 0)
            mexErrMsgIdAndTxt("MATLAB:mexfeval:mxGetString", "Failed to copy function string into allocated memory.");

        buflen = mxGetN(prhs[2]) + 1;
        std::string infile(buflen, '\0');
        status = mxGetString(prhs[2], &infile[0], (mwSize)buflen);

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

        size_t buflen = mxGetN(prhs[2]) + 1;
        std::string outfile(buflen, '\0');
        int status = mxGetString(prhs[2], &outfile[0], (mwSize)buflen);

        baobzi_instance->save(outfile);
        return;
    }

    // Got here, so command not recognized
    mexErrMsgTxt("Command not recognized.");
}
