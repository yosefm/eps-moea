% Run an optimization using the epsilon-moea algorithm.
%
% Arguments: 
% creature - a description of the problem as returned from construct creature.
% pop_size - number of individuals to use for the simulation. The more you use,
%	The better is the pareto front, but it's slower and has diminishing 
%	returns.
% niche_size - the maximum distance between asubjects that generates fitness 
%	punishment (distance in the normed genotype space).
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
	niche_size, conv_gens, num_gens, objectives, grid)
	
	population = gen_population(creature, pop_size);
    niche_size = niche_size^2;
    
	archive_stagnation = 0;

	% Normed values of the genotype, for niching purposes.
	ranges = creature.up_bnds - creature.low_bnds;
	ranges = ranges(ones(pop_size, 1),:);
	normed_pop = (population - creature.low_bnds(ones(pop_size, 1),:))./ranges;
	dist = distances(normed_pop);
	dist(dist == 0) = inf;

	% Initial fitness:
	fitness = objectives(population);
	num_objects = size(fitness, 2);
	parent_archive = pareto_front(fitness);
	
	% Inflict the cholera on crowded populations:
	crowded = dist < niche_size;
	cholera = sum(crowded.*dist/niche_size, 2); % linear niching.
	niched_fit = fitness - cholera(:,ones(num_objects,1));

	while (archive_stagnation < conv_gens) && (num_gens > 0)
		% generate new solution:
		mama = pop_select(niched_fit);
		papa = archive_select(parent_archive);
		offspring = breed(creature, population(mama,:), population(papa,:));
		
		% Offspring niching:
		dist_of = distances(normed_pop, offspring);
		
		% accept the new solution to the population and archive:
		contend_fit = objectives(offspring);
		repl = pop_accept(fitness, contend_fit);
		if (repl == 0)
			% Matlab has no in-place increment. Primitive or what?
			archive_stagnation = archive_stagnation + 1;
			num_gens = num_gens - 1;
			continue
		end
		
		population(repl,:) = offspring;
		normed_pop(repl,:) = (population(repl,:) - creature.low_bnds)./ranges(1,:);
		
		dist(:,repl) = dist_of;
		dist(repl,:) = dist_of;
		dist(repl, repl) = inf;
		
		% Offspring niching:
		% The offspring is punished for being close to existing solutions, but 
		% the existing ones are only judged relative to other existing solutions.
		crowded = dist_of < niche_size;
		cholera = sum(crowded .* dist_of / niche_size); % linear niching.
		niched_cont = contend_fit - cholera;
		
		[new_archive, accepted] = archive_accept(parent_archive, niched_fit, niched_cont, repl, grid);
		
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
