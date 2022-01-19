#include "baobzi.h"
#include "class_handle.hpp"
#include "mex.h"

#include <list>
#include <mutex>
#include <set>
#include <string>

constexpr int n_slots = 10;
static std::set<int> occupied_slots;
static std::mutex mtx;

struct wrapper_data_t {
    std::string fit_fcn;
    uint16_t dim;
    uint16_t order;
};
static wrapper_data_t wrapper_data_arr[n_slots];

int acquire_slot() {
    mtx.lock();
    for (int i = 0; i < n_slots; ++i) {
        if (!occupied_slots.count(i)) {
            occupied_slots.insert(i);
            mtx.unlock();
            return i;
        }
    }
    mtx.unlock();
    mexErrMsgTxt("Unable to aqcuire free baobzi function slot.");
    return -1;
}

void free_slot(int i) {
    mtx.lock();
    occupied_slots.erase(i);
    mtx.unlock();
}

template <int INDEX>
double matfun_wrapper_template(const double *x) {
    mxArray *prhs[1];
    mxArray *plhs[1];
    constexpr wrapper_data_t &data = wrapper_data_arr[INDEX];
    plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
    prhs[0] = mxCreateDoubleMatrix(data.dim, 1, mxREAL);

    double *xptr = mxGetPr(prhs[0]);
    for (int i = 0; i < data.dim; ++i)
        xptr[i] = x[i];

    int status = mexCallMATLAB(1, plhs, 1, prhs, data.fit_fcn.c_str());
    if (status != 0) {
        free_slot(INDEX);
        mexErrMsgIdAndTxt("MATLAB:mexfeval:mxCallMATLAB", "Failed to execute MATLAB function.");
    }

    return mxGetPr(plhs[0])[0];
}

double matfun_wrapper_0(const double *x) { return matfun_wrapper_template<0>(x); }
double matfun_wrapper_1(const double *x) { return matfun_wrapper_template<1>(x); }
double matfun_wrapper_2(const double *x) { return matfun_wrapper_template<2>(x); }
double matfun_wrapper_3(const double *x) { return matfun_wrapper_template<3>(x); }
double matfun_wrapper_4(const double *x) { return matfun_wrapper_template<4>(x); }
double matfun_wrapper_5(const double *x) { return matfun_wrapper_template<5>(x); }
double matfun_wrapper_6(const double *x) { return matfun_wrapper_template<6>(x); }
double matfun_wrapper_7(const double *x) { return matfun_wrapper_template<7>(x); }
double matfun_wrapper_8(const double *x) { return matfun_wrapper_template<8>(x); }
double matfun_wrapper_9(const double *x) { return matfun_wrapper_template<9>(x); }

double (*matfun_wrapper_arr[n_slots])(const double *) = {
    matfun_wrapper_0, matfun_wrapper_1, matfun_wrapper_2, matfun_wrapper_3, matfun_wrapper_4,
    matfun_wrapper_5, matfun_wrapper_6, matfun_wrapper_7, matfun_wrapper_8, matfun_wrapper_9};

class baobzi {
  public:
    baobzi(const std::string &matfun, int dim, int order, const double *center, const double *half_length,
           const double tol)
        : dim_(dim), order_(order) {

        slot_ = acquire_slot();

        wrapper_data_arr[slot_].fit_fcn = matfun;
        wrapper_data_arr[slot_].dim = dim;
        wrapper_data_arr[slot_].order = order;

        // to prevent matlab from segfaulting if your function doesn't exist, try to call it once
        matfun_wrapper_arr[slot_](center);
        obj_ = baobzi_init(matfun_wrapper_arr[slot_], dim_, order_, center, half_length, tol);
    }
    baobzi(const char *infile) {
        baobzi_header_t header = baobzi_read_header_from_file(infile);
        dim_ = header.dim;
        order_ = header.order;
        slot_ = acquire_slot();

        wrapper_data_arr[slot_].dim = dim_;
        wrapper_data_arr[slot_].order = order_;

        obj_ = baobzi_restore(infile);
    }

    ~baobzi() {
        free_slot(slot_);
        if (obj_)
            obj_ = baobzi_free(obj_);
    }
    double eval(const double *x) { return baobzi_eval(obj_, x); };
    void save(const std::string &fname) { baobzi_save(obj_, fname.c_str()); }

  private:
    int dim_;
    int order_;
    int slot_;
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
        if (nrhs != 7)
            mexErrMsgTxt("New: Seven inputs expected.");

        std::string baobzi_fit_fcn = to_string(prhs[1]);

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
