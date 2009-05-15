% Run an optimization using the epsilon-moea algorithm.
%
% Arguments: 
% creature - a description of the problem as returned from construct creature.
% pop_size - number of individuals to use for the simulation. The more you use,
%	The better is the pareto front, but it's slower and has diminishing 
%	returns.
% conv_gens - after this many generations with no change in the archive 
%	population, the iteration stops.
% nun_gens - maximum total number of generations, with or without convergence.
% objectives - a function that given a population array returns the fitness 
%	array.
% grid - the size of the hypercubes in the epsilon-dominance tests.
%
% Returns:
% population - the population matrix after the latest iteration.
% fitness - the fitness matrix after the latest iteration.
% parent_archive - the final archive.

function [population, fitness, parent_archive] = eps_moea_optimize(creature, pop_size, ...
	conv_gens, num_gens, objectives, grid)
	
	population = gen_population(creature, pop_size);
	archive_stagnation = 0;

	% Initial fitness:
	fitness = objectives(population);
	num_objects = size(fitness, 2);
	parent_archive = pareto_front(fitness);
	grid_fit = fitness - rem(fitness, grid(ones(pop_size, 1),:));
    
	while (archive_stagnation < conv_gens) && (num_gens > 0)
		% generate new solution:
		mama = pop_select(fitness, parent_archive);
		papa = archive_select(parent_archive);
		offsprings = breed(creature, population(mama,:), population(papa,:));
		
        for offs = 1:size(offsprings, 2)
            % accept the new solution to the population and archive:
            contend_fit = objectives(offsprings{offs});
            grid_cont = contend_fit - rem(contend_fit, grid);
            
            [parent_archive, accepted] = archive_accept(...
                parent_archive, fitness, grid_fit, contend_fit, grid_cont);
            repl = pop_accept(fitness, contend_fit);
            
            if accepted
                archive_stagnation = 0;
            end
            if (repl == 0)
                num_gens = num_gens - 1;
                continue
            end

            % prepare next iteration:
            population(repl,:) = offsprings{offs};
            fitness(repl,:) = contend_fit;
            parent_archive(repl) = accepted;
            grid_fit(repl,:) = grid_cont;
        end
        
        if ~accepted
            archive_stagnation = archive_stagnation + 1;
        end
		num_gens = num_gens - 1;
	end
end
