% exercise preto_front.m

% sample fitness that is in itself a pereto front:
fit1 = zeros(10, 2);
fit1(:,1) = 1:10;
fit1(:,2) = 1./fit1(:,1);
pf = pareto_front(fit1);
assert(all(pf))

fit2 = fit1 + 0.1;
fit2(:,2) = 2./fit2(:,1);
pf = pareto_front(fit2);
assert(all(pf))
pf = pareto_front([fit1; fit2]);
assert(all(pf(1:10)) & ~any(pf(11:end)))
