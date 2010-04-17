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
% offsprings - a cell array with two length-g arrays with the genotype of 
%   the offspring after recombination and mutation.

function offsprings = breed(creature, mama, papa)
	offsprings = sbx(creature, mama, papa);
	
    % Mutate:
    for offs = 1:size(offsprings, 2)
        which_genes = rand(1, chromosome_length(creature)) < creature.p_mute;
        if ~any(which_genes), continue, end

        % Polinomial mutation, see ref. [1]
        delta = 1 - min(offsprings{offs}(which_genes) - creature.low_bnds(which_genes), ...
            creature.up_bnds(which_genes) - offsprings{offs}(which_genes)) ...
            ./ creature.ranges(which_genes);

        mute_bases = rand(1, sum(which_genes));
        perturb = zeros(size(mute_bases));
        close = mute_bases <= 0.5;

        perturb(close) = (2*mute_bases(close) + ...
            (1 - 2*mute_bases(close)).*delta(close).^(creature.et_m + 1))...
            .^(1./(creature.et_m + 1)) - 1;
        perturb(~close) = 1 - (2*(1 - mute_bases(~close)) + ...
            2*(mute_bases(~close) - 0.5).*delta(~close).^(creature.et_m + 1))...
            .^(1./(creature.et_m + 1));
        offsprings{offs}(which_genes) = offsprings{offs}(which_genes) + ...
            perturb.*creature.ranges(which_genes);
        
        
        underflow = offsprings{offs} < creature.low_bnds;
        offsprings{offs}(underflow) = creature.low_bnds(underflow);
        overflow = offsprings{offs} > creature.up_bnds;
        offsprings{offs}(overflow) = creature.up_bnds(overflow);
    end
end
