% A method of the creature class, giving the cromosome length.
% Arguments: creature - a record representing the creature.

function clen = chromosome_length(creature)
	clen = length(creature.low_bnds);
end