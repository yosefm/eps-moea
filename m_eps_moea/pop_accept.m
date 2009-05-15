% Replace some poor underdog in the population that is pareto-dominated by an 
% up-and-coming contender - if any such underdog exists. Otherwise cast away
% the contender.
%
% Arguments:
% fitness - a p by t matrix for p subjects and t target functions. fitness(i,j)
%	is the value of target function j for subject i.
% contend_fit - a 1 by t vector with the contender's fitness.
%
% Returns:
% repl - the index of the replaced subject, or 0 if no replacement occurred

function repl = pop_accept(fitness, contend_fit)
	% again, due to MATLAB's lack of broadcasting:
	rep_cond_fit = contend_fit(ones(size(fitness, 1), 1),:);
	underdogs = all(rep_cond_fit <= fitness, 2) & ...
		any(rep_cond_fit < fitness, 2);
	
	% If no underdog is found, bail out:
	if any(underdogs)
    	% Replace the first underdog with the contender:
        repl = find(underdogs, 1, 'first');
        return
    end
    
    top_dogs = all(rep_cond_fit >= fitness, 2) & ...
		any(rep_cond_fit > fitness, 2);
    
    if any(top_dogs)
        repl = 0;
        return 
    end
    
    repl = fix(rand()*size(fitness, 1)) + 1;
end
