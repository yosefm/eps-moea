% shake out the pop_accept() function:

fit = zeros(10, 2);
fit(:,1) = 1:10;
fit(:,2) = 1./fit(:,1);

% Test 1: (0,0) fitness dominates everyone, so the first subject will be 
% replaced:
contend_fit = [0, 0];
repl = pop_accept(fit, contend_fit);
assert(repl ~= 0);

% Test 2: (100,100) is dominated by all, it is thrown away.
contend_fit = [100, 100];
repl = pop_accept(fit, contend_fit);
assert(repl == 0);

% Test 3: (1.5, 0.4) dominates the second element, so it is replaced.
contend_fit = [1.5, 0.4];
repl = pop_accept(fit, contend_fit);
assert(repl == 2);
