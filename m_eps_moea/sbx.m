% Create an offspring using simulated binary crossover.
%
% Arguments: 
% mama, papa - two breeders from the population, each a vector of genes.
%
% Returns: 
% offsprings - a 2-tuple of length-g arrays, each with the genotype of an 
%   offspring after recombination and mutation.
    
function offsprings = sbx(creature, mama, papa)
    offsprings = {mama, papa};

    if rand() <= creature.p_recomb
        recombed = (rand(1, chromosome_length(creature)) <= 0.5) & (abs(mama - papa) > 1e-14);
        num_recombed = sum(recombed);
        cum_dist = rand(2, num_recombed);

        parents = [mama(recombed); papa(recombed)];
        p_spread = abs(mama(recombed) - papa(recombed));
        correction = 2 - (1 + [
            min(parents) - creature.low_bnds(recombed); 
            creature.up_bnds(recombed) - max(parents)] ./ ...
            p_spread([1 1],:)).^-(creature.et_c + 1);

        spread = zeros(2, num_recombed);
        cont_bool = cum_dist <= 1./correction;
        [contract_i, contract_j] = find(cont_bool);
        [expand_i, expand_j] = find(~cont_bool);
        spread(contract_i, contract_j) = ...
            (correction(contract_i, contract_j).*...
            cum_dist(contract_i, contract_j)).^(1./(creature.et_c + 1));
        spread(expand_i, expand_j) = ...
            (2 - correction(expand_i, expand_j).*...
            cum_dist(expand_i, expand_j)).^(-1./(creature.et_c + 1));

        pm = rand(2, num_recombed);
        pm(pm > 0.5) = 1;
        pm(pm <= 0.5) = -1;
        padd = mama(recombed) + papa(recombed);
        recombed_traits = 0.5.*(padd([1 1],:) + pm.*spread.*p_spread([1 1], :));

        % Damn MATLAB for the next two lines:
        lower = creature.low_bnds([1 1],recombed);
        upper = creature.up_bnds([1 1],recombed);
        [underflow_i, underflow_j] = find(recombed_traits < lower);
        recombed_traits(underflow_i, underflow_j) = ...
            lower(underflow_i, underflow_j);
        [overflow_i, overflow_j] = find(recombed_traits > upper);
        recombed_traits(overflow_i, overflow_j) = ...
            upper(overflow_i, overflow_j);
        
        % Merge into not recombined traits from one of the parents:
        which_recombed = rand(1, num_recombed) <= 0.5;
        offs_recomb = recombed_traits(1,:);
        offs_recomb(which_recombed) = recombed_traits(2,which_recombed);
        offsprings{1}(recombed) = offs_recomb;
        offs_recomb = recombed_traits(2,:);
        offs_recomb(which_recombed) = recombed_traits(1,which_recombed);
        offsprings{2}(recombed) = offs_recomb;

    end
end
