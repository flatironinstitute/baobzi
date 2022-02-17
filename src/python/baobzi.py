from ctypes import CDLL, CFUNCTYPE, POINTER, c_double, c_void_p, c_uint16, c_int, c_char_p, Structure, pointer
import numpy as np

libbaobzi = CDLL("libbaobzi.so")
INPUT_FUNC = CFUNCTYPE(c_double, POINTER(c_double))


class BAOBZI_STRUCT(Structure):
    _fields_ = [("obj", c_void_p), ("dim", c_int), ("order", c_int),
                ("f_", INPUT_FUNC), ("eval", c_void_p), ("save", c_void_p),
                ("free", c_void_p)]


class BAOBZI_INPUT_STRUCT(Structure):
    _fields_ = [("func", INPUT_FUNC), ("data", c_void_p), ("dim", c_int),
                ("order", c_int), ("tol", c_double)]


baobzi_t = POINTER(BAOBZI_STRUCT)

baobzi_init = libbaobzi.baobzi_init
baobzi_init.restype = baobzi_t
baobzi_init.argtypes = [
    POINTER(BAOBZI_INPUT_STRUCT),
    POINTER(c_double),
    POINTER(c_double)
]

baobzi_eval = libbaobzi.baobzi_eval
baobzi_eval.restype = c_double
baobzi_eval.argtypes = [c_void_p, POINTER(c_double)]

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
            inputdata = BAOBZI_INPUT_STRUCT(INPUT_FUNC(fin), None, dim, order, tol)

            self.ptr = baobzi_init(pointer(inputdata),
                                   center.ctypes.data_as(POINTER(c_double)),
                                   half_length.ctypes.data_as(POINTER(c_double)))
        else:
            print(
                "Baobzi requires either a 'filename' argument or a 'fin' argument"
            )

    def __del__(self):
        if self.ptr:
            baobzi_free(self.ptr)

    def __call__(self, x):
        if x.size == self.dim:
            return baobzi_eval(self.ptr, (c_double * self.dim)(*x))
        else:
            res = np.empty(x.size // self.dim, dtype=np.float64)
            baobzi_eval_multi(self.ptr, x.ctypes.data_as(POINTER(c_double)), res.ctypes.data_as(POINTER(c_double)), res.size)
            return res

    def save(self, filename):
        bfilename = bytes(filename, 'utf-8')
        baobzi_save(self.ptr, bfilename)

    def stats(self):
        baobzi_stats(self.ptr)
