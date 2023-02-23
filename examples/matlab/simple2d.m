dim = 2;
order = 6;
center = [0.0, 0.0];
hl = [1.0, 1.0];
tol = 1E-8;
minimum_leaf_fraction = 0.0;
split_multi_eval = 1;
max_depth = 50;

func_approx = baobzi('new', 'testfun', dim, order, center, hl, tol, minimum_leaf_fraction, split_multi_eval, max_depth);
display(func_approx.eval([0.25, 0.25]))
func_approx.save('simple2d.baobzi');
clear func_approx
func_approx = baobzi('restore', 'simple2d.baobzi');
func_approx.stats()
display(func_approx.eval([0.25, 0.25]))

x = rand([1000000, 2]);
tic; func_approx.eval(x); toc
