% Given multiple fitness values of a population, find which
% individuals are in the pareto-front(i.e. the non-dominated
% subset of the population.
% Assumes the problem is in canonical form, i.e. the target 
% functions are all to be minimized.
%
% Arguments:
% fitness - a p by t array, for p individuals and t target 
%	functions, where fitness(p,t) is the value of function t 
%	for individual p.
%
% Returns:
% pf - a boolean vector of length p, saying which individual is in the front.

function pf = pareto_front(fitness)
	num_subjects = size(fitness, 1);
	% pre-allocate result space, mourning MATLAB's lack of proper 'empty()'.
	pf = zeros(num_subjects, 1); % and no 1D arrays. God help me.
	
	for subject = 1:num_subjects
		others = fitness([1:subject-1, subject+1:end],:);
		% MATLAB has no broadcasting, using cludgy index tricks:
		rep_subj_fit = fitness(subject*ones(num_subjects-1,1),:);
		non_dominated = any(rep_subj_fit < others, 2) | all(rep_subj_fit <= others, 2);
		pf(subject) = all(non_dominated);
	end
end
