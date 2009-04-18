% Test function 1 from Tzitzler E. et al, Comparison of Multiobjective 
% Evolutionary Algorithms: Empirical Results, Evolutionary computation, 2000,
% vol 8(2), pp. 173-195.
%
% Arguments: 
% contenders - those whose fitness is to be evaluated
% pop_size - the size of the population in which the contenders live.

function fitness = tau1(contenders)
	fit1 = contenders(:,1);
	g = 1 + 9/(size(contenders, 2) - 1) .* sum(contenders(:,2:end), 2);
	h = 1 - sqrt(fit1 ./ g);
	fitness = [fit1, g.*h];
end
