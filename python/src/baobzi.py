from ctypes import CDLL, CFUNCTYPE, POINTER, c_double, c_void_p, c_uint16, c_int, c_char_p, Structure

libbaobzi = CDLL("libbaobzi.so")
INPUT_FUNC = CFUNCTYPE(c_double, POINTER(c_double))


class BAOBZI_STRUCT(Structure):
    _fields_ = [("obj", c_void_p), ("dim", c_int), ("order", c_int),
                ("f_", INPUT_FUNC), ("eval", c_void_p), ("save", c_void_p),
                ("free", c_void_p)]


baobzi_t = POINTER(BAOBZI_STRUCT)

baobzi_init = libbaobzi.baobzi_init
baobzi_init.restype = baobzi_t
baobzi_init.argtypes = [
    INPUT_FUNC, c_uint16, c_uint16,
    POINTER(c_double),
    POINTER(c_double), c_double
]

baobzi_eval = libbaobzi.baobzi_eval
baobzi_eval.restype = c_double
baobzi_eval.argtypes = [c_void_p, POINTER(c_double)]

baobzi_free = libbaobzi.baobzi_free
baobzi_free.restype = baobzi_t
baobzi_free.argtypes = [c_void_p]

baobzi_save = libbaobzi.baobzi_save
baobzi_save.restype = c_void_p
baobzi_save.argtypes = [c_void_p, c_char_p]

baobzi_restore = libbaobzi.baobzi_restore
baobzi_restore.restype = baobzi_t
baobzi_restore.argtypes = [INPUT_FUNC, c_char_p]


class Baobzi:
    def __init__(self,
                 fin,
                 dim=None,
                 order=None,
                 center=None,
                 half_length=None,
                 tol=None,
                 filename=None):
        self.f = INPUT_FUNC(fin)
        if filename:
            self.ptr = baobzi_restore(self.f, filename)
            self.dim = self.ptr[0].dim
            self.order = self.ptr[0].order
        else:
            self.dim = dim
            self.order = order
            self.ptr = baobzi_init(self.f, self.dim, self.order,
                                   (c_double * dim)(*center),
                                   (c_double * dim)(*half_length), tol)

    def __del__(self):
        baobzi_free(self.ptr)

    def __call__(self, x):
        return baobzi_eval(self.ptr, (c_double * self.dim)(*x))

    def save(self, filename):
        baobzi_save(self.ptr, filename)
