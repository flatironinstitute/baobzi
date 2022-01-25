from baobzi import Baobzi

def py_test_func(x):
    return x[0] * x[1]

center = [0.0, 0.0]
hl = [1.0, 1.0]
point = [0.25, 0.25]

test = Baobzi(py_test_func, 2, 6, center, hl, 1E-8)
test.save('test.baobzi')
print(test(point))
del test

test2 = Baobzi(filename='test.baobzi')
print(test2(point))
del test2
