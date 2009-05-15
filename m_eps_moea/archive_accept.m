% Update the pareto-front in the archive using epsilon-dominance.
%
% Arguments:
% archive - a p by 1 boolean array stating which individual is in the archive.
% fitness - a p by t array whose element (i,j) is the value of the jth fitness
%	function of t functions, for the ith individual of p individuals, after the 
%	last iteration of the population acceptance method.
% contend_fit - the fitness of the contender.
% 
% Returns:
% archive - the modified archive.
% accepted - a boolean value stating whether the contender was accepted.

function [archive, accepted] = archive_accept(archive, fitness, grid_fit,...
    contend_fit, grid_cont)

    % We have more than one level of selection from the population,
    % so its convenient to work with indices rather than boolean arrays.
    archive_idxs = find(archive);
    archive_size = sum(archive);
    arch_grid_fit = grid_fit(archive_idxs,:);
    arch_fit = fitness(archive_idxs,:);
    
    % We look for opportunities to reject the contender, along the way removing 
    % any dominated member of the archive:
    for vertex = unique(arch_grid_fit, 'rows')'
        vertex = vertex';
        
        high = any(grid_cont > vertex);
        low = any(grid_cont < vertex);
        if ~(high || low)
            % This hypercube may have non-dominated members.
            in_grid = find(all(arch_grid_fit == vertex(ones(archive_size, 1),:), 2));
            rep_cont_fit = contend_fit(ones(length(in_grid), 1),:);
            
            lower_fit = any(arch_fit(in_grid,:) < rep_cont_fit);
            higher_fit = any(arch_fit(in_grid,:) > rep_cont_fit);
            if any(lower_fit & ~higher_fit)
                accepted = 0;
                return
            end

            archive(archive_idxs(in_grid)) = 0;
            
            underdogs = higher_fit & ~lower_fit;
            if all(underdogs), break, end
            
            % Of the remaining solutions, the closest to the grid is taken:                
            dist_squares = sum((arch_fit(in_grid(~underdogs),:) - ...
                arch_grid_fit(in_grid(~underdogs),:)).^2, 2);
            remaining = dist_squares <= sum((contend_fit - grid_cont).^2);
            
            if any(remaining)
                archive(archive_idxs(in_grid)) = 1;
                accepted = 0;
                return
            end
        
        elseif high && ~low
            % The contender is dominated. The function can end, because
            % it's impossible that others in the archive may still be dominated.
            accepted = 0;
            return
            
        elseif low && ~high
            % This hypercube is dominated by the new solution, exclude all its
            % members from the archive:
            in_grid = find(all(arch_grid_fit == vertex(ones(archive_size, 1),:), 2));
            archive(archive_idxs(in_grid)) = 0;
        end
    end
    
    accepted = 1;
    return
end