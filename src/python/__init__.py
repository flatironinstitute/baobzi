from ctypes import CDLL, CFUNCTYPE, POINTER, c_double, c_void_p, c_uint16, c_int, c_char_p, Structure, pointer
from ctypes.util import find_library
import numpy as np

baobzi_path = find_library('baobzi')

if not baobzi_path:
    import os
    from sys import platform
    libroot = os.path.sep.join(os.path.realpath(__file__).split(os.path.sep)[:-5])
    if platform == 'linux':
        extension = ".so"
    elif platform == 'darwin':
        extension = ".dylib"
    elif platform == 'win32':
        extension = '.dll'
    else:
        raise RuntimeError("Invalid operating platform found.")
    lib = os.path.join(libroot, "lib", "libbaobzi" + extension)
    lib64 = os.path.join(libroot, "lib64", "libbaobzi" + extension)
    if os.path.exists(lib):
        baobzi_path = lib
    elif os.path.exists(lib64):
        baobzi_path = lib64

if not baobzi_path:
    raise OSError("Unable to find 'libbaobzi'. Add path to its containing directory to your LD_LIBRARY_PATH variable.")

libbaobzi = CDLL(baobzi_path)

INPUT_FUNC = CFUNCTYPE(c_double, POINTER(c_double))


class BAOBZI_STRUCT(Structure):
    _fields_ = (("obj", c_void_p), ("dim", c_int), ("order", c_int),
                ("f_", INPUT_FUNC), ("eval", c_void_p), ("eval_multi", c_void_p),
                ("stats", c_void_p), ("save", c_void_p), ("free", c_void_p),)


class BAOBZI_INPUT_STRUCT(Structure):
    _fields_ = [("func", INPUT_FUNC), ("data", c_void_p), ("dim", c_int),
                ("order", c_int), ("tol", c_double), ("minimum_leaf_fraction", c_double),
                ("split_multi_eval", c_int),("max_depth", c_int)]

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

baobzi_save = libbaobzi.baobzi_save
baobzi_save.restype = c_void_p
baobzi_save.argtypes = [c_void_p, c_char_p]

baobzi_restore = libbaobzi.baobzi_restore
baobzi_restore.restype = baobzi_t
baobzi_restore.argtypes = [c_char_p]


class Baobzi:
    def __init__(self,
                 fin=None,
                 dim=None,
                 order=None,
                 center=None,
                 half_length=None,
                 tol=None,
                 minimum_leaf_fraction=0.0,
                 split_multi_eval=1,
                 max_depth=50,
                 filename=None):
        self.ptr = None
        if filename:
            bfilename = bytes(filename, 'utf-8')
            self.ptr = baobzi_restore(bfilename)
            self.dim = self.ptr[0].dim
            self.order = self.ptr[0].order
        elif fin:
            if not (dim and order and center.size and half_length.size and tol):
                print(
                    "Baobzi: supply dim, order, center, half_length, and tol for init"
                )
            self.dim = dim
            self.order = order
            inputdata = BAOBZI_INPUT_STRUCT(INPUT_FUNC(fin), None, dim, order, tol, minimum_leaf_fraction, split_multi_eval, max_depth)

            self.ptr = baobzi_init(pointer(inputdata),
                                   center.ctypes.data_as(POINTER(c_double)),
                                   half_length.ctypes.data_as(POINTER(c_double)))
        else:
            print(
                "Baobzi requires either a 'filename' argument or a 'fin' argument"
            )

        if self.ptr[0] is None or self.ptr[0].obj is None:
            raise RuntimeError("Unable to create baobzi object")


    def __del__(self):
        baobzi_free(self.ptr)

    def __call__(self, x):
        xarr = np.array(x, dtype=np.float64)
        res = np.empty(xarr.size // self.dim, dtype=np.float64)
        baobzi_eval_multi(self.ptr, xarr.ctypes.data_as(POINTER(c_double)), res.ctypes.data_as(POINTER(c_double)), res.size)
        return res

    def save(self, filename):
        bfilename = bytes(filename, 'utf-8')
        baobzi_save(self.ptr, bfilename)

    def stats(self):
        baobzi_stats(self.ptr)
