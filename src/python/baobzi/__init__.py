from ctypes import CDLL, CFUNCTYPE, POINTER, c_double, c_void_p, c_int, Structure, pointer
import os
import sys
import numpy as np

def _find_baobzi_lib():
    pkg_dir = os.path.dirname(__file__)
    if sys.platform == "linux":
        libname = "libbaobzi.so"
    elif sys.platform == "darwin":
        libname = "libbaobzi.dylib"
    elif sys.platform == "win32":
        libname = "libbaobzi.dll"
    else:
        raise RuntimeError("Unsupported platform")
    # Look in package directory
    candidate = os.path.join(pkg_dir, libname)
    if os.path.exists(candidate):
        return candidate
    # Optionally, look in lib/ or lib64/ subdirs
    for subdir in ["lib", "lib64"]:
        candidate = os.path.join(pkg_dir, "..", "..", subdir, libname)
        if os.path.exists(candidate):
            return candidate
    raise OSError("Unable to find 'libbaobzi' in package data.")

baobzi_path = _find_baobzi_lib()
libbaobzi = CDLL(baobzi_path)

INPUT_FUNC = CFUNCTYPE(None, POINTER(c_double), POINTER(c_double), c_void_p)

def _make_callback(pyfunc, m, n):
    def callback(x_ptr, y_ptr, _):
        # Convert input pointer to numpy array
        x = np.ctypeslib.as_array(x_ptr, shape=(m,))
        # Call the Python function
        y = pyfunc(x)
        # Write result to output pointer
        y_out = np.ctypeslib.as_array(y_ptr, shape=(n,))
        y_out[:] = y
    return INPUT_FUNC(callback)


class BAOBZI_STRUCT(Structure):
    _fields_ = (("obj", c_void_p),)


class BAOBZI_INPUT_STRUCT(Structure):
    _fields_ = [("func", INPUT_FUNC),
                ("data", c_void_p),
                ("input_dim", c_int),
                ("output_dim", c_int),
                ("degree", c_int),
                ("tol", c_double),
                ("minimum_leaf_fraction", c_double),
                ("split_multi_eval", c_int),
                ("min_depth", c_int),
                ("max_depth", c_int),
                ("tol_type", c_int),
                ("n_samples_per_dim", c_int),
                ]

baobzi_t = POINTER(BAOBZI_STRUCT)

baobzi_init = libbaobzi.baobzi_init
baobzi_init.restype = baobzi_t
baobzi_init.argtypes = [
    POINTER(BAOBZI_INPUT_STRUCT),
    POINTER(c_double),
    POINTER(c_double)
]

baobzi_eval_multi = libbaobzi.baobzi_eval_multi
baobzi_eval_multi.restype = c_void_p
baobzi_eval_multi.argtypes = [c_void_p, POINTER(c_double), POINTER(c_double), c_int]

baobzi_stats = libbaobzi.baobzi_stats
baobzi_stats.restype = baobzi_t
baobzi_stats.argtypes = [c_void_p]

baobzi_free = libbaobzi.baobzi_free
baobzi_free.restype = baobzi_t
baobzi_free.argtypes = [c_void_p]


class Baobzi:
    def __init__(self,
                 fin=None,
                 input_dim=None,
                 degree=None,
                 center=None,
                 half_length=None,
                 tol=None,
                 output_dim=1,
                 minimum_leaf_fraction=0.0,
                 split_multi_eval=1,
                 min_depth=0,
                 max_depth=50,
                 tol_type=0,
                 n_samples_per_dim=None):
        self.ptr = None
        if fin:
            if not (input_dim and degree and center.size and half_length.size and tol):
                print(
                    "Baobzi: supply dim, order, center, half_length, and tol for init"
                )
            self.input_dim = input_dim
            self.degree = degree
            self.output_dim = output_dim
            n_samples_per_dim = n_samples_per_dim if n_samples_per_dim else degree

            func = _make_callback(fin, input_dim, output_dim)

            inputdata = BAOBZI_INPUT_STRUCT(func, None, input_dim, output_dim, degree, tol,
                                            minimum_leaf_fraction, split_multi_eval, min_depth, max_depth,
                                            tol_type,
                                            n_samples_per_dim,
                                            )

            self.ptr = baobzi_init(pointer(inputdata),
                                   center.ctypes.data_as(POINTER(c_double)),
                                   half_length.ctypes.data_as(POINTER(c_double)))
        else:
            print(
                "Baobzi requires a 'fin' argument"
            )

        if self.ptr[0] is None or self.ptr[0].obj is None:
            raise RuntimeError("Unable to create baobzi object")


    def __del__(self):
        baobzi_free(self.ptr)

    def __call__(self, x):
        xarr = np.array(x, dtype=np.float64)
        n_points = xarr.size // self.input_dim
        res = np.empty(self.output_dim * n_points, dtype=np.float64)
        baobzi_eval_multi(self.ptr, xarr.ctypes.data_as(POINTER(c_double)), res.ctypes.data_as(POINTER(c_double)), n_points)
        return res

    def stats(self):
        baobzi_stats(self.ptr)
