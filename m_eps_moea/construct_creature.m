% Constructor for the creature class. 
%
% Arguments: 
% low_bnds - low-bounds for each of the real-valued genes in the
%	creature's chromosome.
% up_bnds - upper-bounds for the same. Must have the same length as 
%	low_bnds
% mutation_chance - the probability of a mutation happening.
% et_m - strength of mutation (see ref [1])
% recomb_chance - chance of any recombination happening
% et_c - strength of recombination distribution function [1].
%
% Returns: a record representing the creature.
%
% Reference:
% [1] Kalyanmoy Deb, An efficient constraint handling method for genetic 
% algorithms, 31 May 2000

function creature = construct_creature(low_bnds, up_bnds, ...
    mutation_chance, et_m, recomb_chance, et_c)

    creature.up_bnds = up_bnds;
	creature.low_bnds = low_bnds;
	creature.p_mute = mutation_chance;
	creature.ranges = up_bnds - low_bnds;
    
    % Matlab's way of default arguments. So primitive.
    if nargin < 6
        et_c = 15;
    end
    if nargin < 5
        recomb_chance = 1;
    end
    if nargin < 4
        et_m = 20;
    end
    creature.p_recomb = recomb_chance;
    creature.et_c = et_c;
    creature.et_m = et_m;
    
	% Validation:
	if length(creature.up_bnds) ~= length(creature.low_bnds)
		error('Upper bounds number not equal to lower-bounds number.');
	end
	if any(creature.up_bnds <= creature.low_bnds)
		error('Upper bounds must be > lower bounds.');
	end
end
