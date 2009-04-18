% creates a new genome for a subject by recombination of parent genes, and
% possibly mutation of the result, depending on the creature's mutation
% resistance.
%
% Arguments: 
% creature - a record representing the genome rules, as returned from 
%	construct_creature().
% mama, papa - two breeders from the population, each a vector of genes.
% 
% Returns: 
% offspring - a length-g array with the genotype of the offspring after 
%	recombination and mutation.

function offspring = breed(creature, mama, papa)
	% recombination place:
	recomb = fix(rand() * chromosome_length(creature) + 1);
	offspring = [mama(1:recomb), papa(recomb+1:end)];
	
	% Possibly mutate:
	if rand() > creature.p_mute; return; end
	gene = fix(rand() * chromosome_length(creature) + 1);
	% mutation is done by assigning a new value to the gene in the allowed range.
	offspring(gene) = rand() * (creature.up_bnds(gene) - creature.low_bnds(gene)) + ...
		creature.low_bnds(gene);
end
