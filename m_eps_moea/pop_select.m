% Randomly compares two subjects from a population and selects the dominating
% subject. If no dominance exists, select one of them randomly.
% 
% Arguments: 
% fitness - a p by t array, for p individuals and t target 
%	functions, where fitness(p,t) is the value of function t 
%	for individual p.
% archive - a boolean vector stating which of the population is in the
%   archive.
%
% Returns:
% sel - the index of the subject selected for breeding.

function sel = pop_select(fitness, archive)
	% select two indexes randomly
	compete = fix(rand(1, 2) * size(fitness, 1)) + 1;
	if all(fitness(compete(1),:) <= fitness(compete(2),:)) && ...
		any(fitness(compete(1),:) < fitness(compete(2),:))
		sel = compete(1);
    elseif all(fitness(compete(2),:) <= fitness(compete(1),:)) && ...
		any(fitness(compete(2),:) < fitness(compete(1),:))
		sel = compete(2);
    else
        sel = compete((rand() <= 0.5) + 1);
	end
end
