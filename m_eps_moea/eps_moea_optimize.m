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
	
	while (archive_stagnation < conv_gens) && (num_gens > 0)
		% generate new solution:
		mama = pop_select(fitness, parent_archive);
		papa = archive_select(parent_archive);
		offspring = breed(creature, population(mama,:), population(papa,:));
		
		% accept the new solution to the population and archive:
		contend_fit = objectives(offspring);
		repl = pop_accept(fitness, contend_fit);
		if (repl == 0)
			archive_stagnation = archive_stagnation + 1;
			num_gens = num_gens - 1;
			continue
		end
		
		population(repl,:) = offspring;
		[new_archive, accepted] = archive_accept(parent_archive, ...
            fitness, contend_fit, repl, grid);
		
		% prepare next iteration:
		fitness(repl,:) = contend_fit;
		if accepted
			archive_stagnation = 0;
		else
			archive_stagnation = archive_stagnation + 1;
		end
		parent_archive = new_archive;
		num_gens = num_gens - 1;
	end
end
