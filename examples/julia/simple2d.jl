import baobzi

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

func_approx = baobzi.init(testfunc, dim, order, center, hl, tol)
baobzi.stats(func_approx)
println(baobzi.eval(func_approx, test_point) - testfunc(pointer(test_point)))

baobzi.save(func_approx, output_file)
baobzi.free(func_approx)

func_approx = baobzi.restore(output_file)
println(baobzi.eval(func_approx, test_point) - testfunc(pointer(test_point)))

points = 2.0 * (rand(Float64, 2000000)) .- 1.0
baobzi.eval_multi(func_approx, points)
