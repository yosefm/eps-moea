% Update the pareto-front in the archive using epsilon-dominance.
%
% Arguments:
% archive - a p by 1 boolean array stating which individual is in the archive.
% fitness - a p by t array whose element (i,j) is the value of the jth fitness
%	function of t functions, for the ith individual of p individuals, after the 
%	last iteration of the population acceptance method.
% contend_fit - the fitness of the contender.
% contend_idx - the index in the population of the contender to replace some of 
%	the archive.
% grid - a vector of length t specifying the edge-length of the n-dim grid in 
%	each of the fitness space axes.
%
% Returns:
% archive - the modified archive.
% accepted - a boolean value stating whether the contender was accepted.

function [archive, accepted] = archive_accept(archive, fitness, contend_fit, contend_idx, grid)
	% We have more than one level of selection from the population,
	% so its convenient to work with indices rather than boolean arrays.
	archive_idxs = find(archive);
    
	% Calculate the grid-fitness of the population:
	archive_size = sum(archive);
    
    dists = rem(fitness(archive_idxs,:), grid(ones(archive_size, 1),:));
    grid_fit = fitness(archive_idxs,:) - dists;

	% and of the contender:
    cont_dist = rem(contend_fit, grid);
    cont_dist_square = sum(cont_dist.^2);
	grid_cont = contend_fit - cont_dist;
	
	% Now check each relevant grid-point in turn.
	for vertex = unique(grid_fit, 'rows')'
		vertex = vertex';
		if all(grid_cont == vertex)
			% This hypercube may have non-dominated members.
			in_grid = find(all(grid_fit == vertex(ones(archive_size, 1),:), 2));
			
			rep_cond_fit = contend_fit(ones(length(in_grid), 1),:);
			top_dogs = all(fitness(archive_idxs(in_grid),:) <= rep_cond_fit, 2) & ...
				any(fitness(archive_idxs(in_grid),:) < rep_cond_fit, 2);
			if any(top_dogs)
				accepted = 0;
                archive(contend_idx) = 0;
				return
			end

			underdogs = all(rep_cond_fit <= fitness(archive_idxs(in_grid),:), 2) & ...
				any(rep_cond_fit < fitness(archive_idxs(in_grid),:), 2);
			archive(archive_idxs(in_grid(underdogs))) = 0;
            if all(underdogs), continue, end
            
            dist_squares = sum(dists(in_grid(~underdogs)).^2, 2);
            remaining = dist_squares < cont_dist_square;
            
            if any(remaining)
                [mn, argmin] = min(dist_squares);
                archive(archive_idxs(in_grid(argmin))) = 1;
                archive(contend_idx) = 0;
                accepted = 0;
                return
            else
                archive(archive_idxs(in_grid(~underdogs))) = 0;
			end

		elseif all(grid_cont >= vertex) && any(grid_cont > vertex)
			% The contender is dominated. The function can end, because
			% it's impossible that others in the archive may still be 
			% dominated.
			archive(contend_idx) = 0; % make sure no leftovers from last iteration.
			accepted = 0;
			return
		
		elseif all(grid_cont <= vertex) && any(grid_cont < vertex)
			% This hypercube is dominated by the new solution, exclude all its
			% members from the archive:
			in_grid = all(grid_fit == vertex(ones(archive_size, 1),:), 2);
			archive(archive_idxs(in_grid)) = 0;
		end
	end %% for vertex = ...
    
    accepted = 1;
    archive(contend_idx) = 1;
end