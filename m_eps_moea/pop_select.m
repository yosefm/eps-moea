% Randomly mates two subjects from a population and selects the dominating
% subject. If no dominance exists, select the second one (this is arbitrary).
% 
% Arguments: 
% fitness - a p by t array, for p individuals and t target 
%	functions, where fitness(p,t) is the value of function t 
%	for individual p.
% archive - a boolean vector stating which of the population is in the
%   archive.
% Returns:
% sel - the index of the subject selected for breeding.

function sel = pop_select(fitness, archive)
    legible = find(~archive);
	% select two indexes randomly
	compete = legible(fix(rand(1, 2) * (size(fitness, 1) - sum(archive)) + 1));
	if all(fitness(compete(1),:) <= fitness(compete(2),:)) && ...
		any(fitness(compete(1),:) < fitness(compete(2),:))
		sel = compete(1);
	else
		sel = compete(2);
	end
end
