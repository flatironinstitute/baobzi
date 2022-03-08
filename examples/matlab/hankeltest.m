dim = 2;
order = 8;
center = [0.0, 1.0];
hl = [3.0, 0.99999];
tol = 1E-10;
minimum_leaf_fraction = 0.0;
split_multi_eval = 1;

hreal = baobzi('new', 'hankelreal', dim, order, center, hl, tol, minimum_leaf_fraction, split_multi_eval);
himag = baobzi('new', 'hankelimag', dim, order, center, hl, tol, minimum_leaf_fraction, split_multi_eval);

% Get a interweaving array to evaluate on the grid
x = (center(1) - hl(1)):0.001:(center(1) + hl(1)); x(end) = [];
y = (center(2) - hl(2)):0.001:(center(2) + hl(2)); y(end) = [];
[X, Y] = meshgrid(x, y);
Z = [X(:)'; Y(:)']; Z = Z(:);

% level the playing field by not distributing matlab's computation
Nthr = maxNumCompThreads(1);
fprintf("Timing matlab's hankel function (unthreaded)\n");
tic; real(besselh(0, X + 1* i * Y)); dt = toc;
fprintf("Elapsed time: %gs\n", dt);
fprintf("Mevals/s: %g\n\n", length(Z) / 2 / dt / 1E6);

% reset number of threads to not be a jerk
maxNumCompThreads(Nthr);
fprintf("Timing matlab's hankel function (threaded)\n");
tic; real(besselh(0, X + 1* i * Y)); dt = toc;
fprintf("Elapsed time: %gs\n", dt);
fprintf("Mevals/s: %g\n\n", length(Z) / 2 / dt / 1E6);

fprintf("Timing baobzi's hankel function (unthreaded)\n");
tic; hreal.eval(Z); himag.eval(Z); dt = toc;
fprintf("Elapsed time: %gs\n", dt);
fprintf("Mevals/s: %g\n\n", length(Z) / 2 / dt / 1E6);
