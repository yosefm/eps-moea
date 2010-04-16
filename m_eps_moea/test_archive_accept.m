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
    eps_fit = fit - rem(fit, grid(:,ones(20, 1))');
    eps_cont = contend_fit - rem(contend_fit, grid');
	[new_arch, accepted] = archive_accept(archive, fit, eps_fit, ...
        contend_fit, eps_cont);
	assert(accepted);
	assert(all(new_arch == 0));
end

% Test 2: (2.1, 0.51) is dominated by the archive, the archive must return 
% unchanged.
contend_fit = [2.1, 0.51];
repl = pop_accept(fit, contend_fit);

for grid = grids
    eps_fit = fit - rem(fit, grid(:,ones(20, 1))');
    eps_cont = contend_fit - rem(contend_fit, grid');
	[new_arch, accepted] = archive_accept(archive, fit, eps_fit, ...
        contend_fit, eps_cont);
	assert(~accepted);
	assert(all(new_arch == archive));
end

% Test 3: (1.5, 0.4) dominates the second element of the front, so it is 
% replaced in the archive.
contend_fit = [1.5, 0.4];
repl = pop_accept(fit, contend_fit);

for grid = grids(:,1:2)
    eps_fit = fit - rem(fit, grid(:,ones(20, 1))');
    eps_cont = contend_fit - rem(contend_fit, grid');
	[new_arch, accepted] = archive_accept(archive, fit, eps_fit, ...
        contend_fit, eps_cont);

    assert(accepted);
    % Only the place of the second item of the front (index 3 in the
    % population) is changed.
    assert(new_arch(3) == 0);
    new_arch(3) = 1;
    assert(all(new_arch == archive));
end

% on the coarse grid solutions are lost:
eps_fit = fit - rem(fit, coarse_grid(:,ones(20, 1))');
eps_cont = contend_fit - rem(contend_fit, coarse_grid');
[new_arch, accepted] = archive_accept(archive, fit, eps_fit, ...
        contend_fit, eps_cont);
assert(accepted);
assert(all(new_arch == [1, zeros(1,19)]'))
