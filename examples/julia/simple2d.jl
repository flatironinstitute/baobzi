import Baobzi

function testfunc(xp::Ptr{Float64})::Cdouble
    x = unsafe_load(xp, 1)
    y = unsafe_load(xp, 2)
    return x * y
end

center = [0.0, 0.0]
hl = [1.0, 1.0]
test_point = [0.25, 0.25]
dim = 2
order = 6
tol = 1E-8
output_file = "simple2d.baobzi"
split_multi_eval = 0

func_approx = Baobzi.init(testfunc, dim, order, center, hl, tol, 0)
Baobzi.stats(func_approx)
println(Baobzi.eval(func_approx, test_point) - testfunc(pointer(test_point)))

Baobzi.save(func_approx, output_file)
Baobzi.free(func_approx)

func_approx = Baobzi.restore(output_file)
println(Baobzi.eval(func_approx, test_point) - testfunc(pointer(test_point)))

points = 2.0 * (rand(Float64, 2000000)) .- 1.0
Baobzi.eval_multi(func_approx, points)
