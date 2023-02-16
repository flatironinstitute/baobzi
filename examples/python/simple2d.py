from baobzi import Baobzi
import numpy as np
import time

def py_test_func(x):
    return x[0]**3 * x[1]**2 + 0.1

center = np.array([1.0, 1.0])
hl = np.array([1.0, 1.0])
point = np.array([0.25, 0.25])

# create, save, and delete baobzi object
test = Baobzi(py_test_func, 2, 6, center, hl, 1E-8)
test.stats()
test.save('test.baobzi')
print("Error in chosen point: ", 1-test(point)/py_test_func(point))
del test

# reload baobzi object from disk
test2 = Baobzi(filename='test.baobzi')

# time baobzi object megaevals/s
n_points = int(1e7)
points = 2.0 * (np.array(np.random.uniform(size=2*n_points))).reshape(n_points, 2)

st = time.time()
res = test2(points)
dt = time.time() - st
print("baobzi MEvals/s:", n_points / dt / 1E6)

st = time.time()
res_numpy = py_test_func(points.T)
dt = time.time() - st
print("numpy MEvals/s:", n_points / dt / 1E6)

print("Max relative error:", np.max(np.abs(1.0-res/res_numpy)))
