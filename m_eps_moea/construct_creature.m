% Constructor for the creature class. 
%
% Arguments: 
% low_bnds - low-bounds for each of the real-valued genes in the
%	creature's chromosome.
% up_bnds - upper-bounds for the same. Must have the same length as 
%	low_bnds
% mutation_chance - the probability of a mutation happening.
% Returns: a record representing the creature.

function creature = construct_creature(low_bnds, up_bnds, mutation_chance)
	creature.up_bnds = up_bnds;
	creature.low_bnds = low_bnds;
	creature.p_mute = mutation_chance;
	
	% Validation:
	if length(creature.up_bnds) ~= length(creature.low_bnds)
		error('Upper bounds number not equal to lower-bounds number.');
	end
	if any(creature.up_bnds <= creature.low_bnds)
		error('Upper bounds must be > lower bounds.');
	end
end
