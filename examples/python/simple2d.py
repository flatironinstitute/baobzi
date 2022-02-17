from baobzi import Baobzi
import numpy as np
import time

def py_test_func(x):
    return x[0] * x[1]

center = np.array([0.0, 0.0])
hl = np.array([1.0, 1.0])
point = np.array([0.25, 0.25])

# create, save, and delete baobzi object
test = Baobzi(py_test_func, 2, 6, center, hl, 1E-8)
test.stats()
test.save('test.baobzi')
print(test(point))
del test

# reload baobzi object from disk
test2 = Baobzi(filename='test.baobzi')
print(test2(point))

# time baobzi object megaevals/s
points = 2.0 * (np.array(np.random.uniform(size=20000000)) - 0.5)
st = time.time()
res = test2(points)
dt = time.time() - st
print(points.size / dt / 1E6)
