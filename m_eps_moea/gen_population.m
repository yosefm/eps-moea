% Creates a matrix representing the population in a genetic algorithm.
% Arguments: creature - a record describing the rules of creating chromosomes,
%		see <construct_creature()> for the details.
%	num_subjects - the number of creatures in the population.

function pop = gen_population(creature, num_subjects)
	pop = rand(num_subjects, chromosome_length(creature));
	scale = creature.up_bnds - creature.low_bnds;
	index_repl = ones(num_subjects, 1);
	pop = pop .* scale(index_repl,:) + creature.low_bnds(index_repl,:);
end
