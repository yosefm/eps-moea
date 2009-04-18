% Test the archive acceptance function:

% Sample fitness made of a population along two curves:
% 1. x*y = 1
% 2. x*y = 2

fit = zeros(20,2);
fit(1:2:20,1) = 1:10;
fit(2:2:20,1) = 1.1:10.1;
fit(1:2:20,2) = 1./fit(1:2:20,1);
fit(2:2:20,2) = 2./fit(2:2:20,1);

archive = pareto_front(fit);

% For the epsilon-dominance tests:
fine_grid = [0.01; 0.01];
med_grid = [0.1; 0.1];
coarse_grid = [0.5; 0.5];
grids = [fine_grid, med_grid, coarse_grid];

% Test 1: (0,0) is better than everyone, it should be the only one in the new
% archive.
contend_fit = [0, 0];
repl = pop_accept(fit, contend_fit);

for grid = grids
	[new_arch, accepted] = archive_accept(archive, fit, contend_fit, repl, grid');
	assert(accepted);
	assert(new_arch(1) == 1 && all(new_arch(2:end) == 0));
end

% Test 2: (2.1, 0.51) is dominated by the archive, the archive must return 
% unchanged.
contend_fit = [2.1, 0.51];
repl = pop_accept(fit, contend_fit);

for grid = grids
	[new_arch, accepted] = archive_accept(archive, fit, contend_fit, repl, grid');
	assert(~accepted);
	assert(all(new_arch == archive));
end

% Test 3: (1.5, 0.4) dominates the second element of the front, so it is 
% replaced in the archive.
contend_fit = [1.5, 0.4];
repl = pop_accept(fit, contend_fit);

for grid = grids(:,1:2)
	[new_arch, accepted] = archive_accept(archive, fit, contend_fit, repl, grid');
	assert(accepted);
	% because the kicked-out solution was at the front, the representation is 
	% unchanged.
	assert(all(new_arch == archive));
end

% on the coarse grid solutions are lost:
[new_arch, accepted] = archive_accept(archive, fit, contend_fit, repl, coarse_grid');
assert(accepted);
assert(all(new_arch == [1, 0, 1, zeros(1,17)]'))
