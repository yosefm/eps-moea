% calculates normalized euclidean distance between each two elements in the population.
%
% Arguments: 
% population - the points in the decision space for which to find distances,
%	normed to the ranges allowed for each variable.
% compare_to - optional element of the population that we compare only to him.
%
% Returns: a square matrix of distance from individual i to individual j, or a 
%	vector of distances to one element if compare_to is given.

function dist = distances(population, compare_to)
	num_subjects = size(population, 1);
	if nargin == 2
		dist = sum((population - compare_to(ones(num_subjects, 1),:)).^2, 2);
		return
	end

	res_size = [num_subjects, num_subjects];
	dist = zeros(res_size);
	% All this setup just because Matlab has no ndenumerate like Python.
	indices = 1:(num_subjects^2);
	
	[ii, jj] = ind2sub(res_size, indices);
	for elem = indices
		dist(ii(elem), jj(elem)) = sum((population(ii(elem),:) - population(jj(elem),:)).^2);
	end
end
