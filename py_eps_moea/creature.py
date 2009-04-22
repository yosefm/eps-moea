# This class represents a creature (problem) genotype (decision space)
#
# References:
# [1] Kalyanmoy Deb, An efficient constraint handling method for genetic 
# algorithms, 31 May 2000

from numpy import random
import numpy as N

class Creature(object):
    def __init__(self, low_bnds, up_bnds, mutation_chance, mutation_prm=20):
        if len(low_bnds) != len(up_bnds):
            raise ValueError("Upper bounds number not equal to lower-bounds number.")
        if (low_bnds > up_bnds).any():
            raise ValueError("Upper bounds must be > lower bounds.")
        
        self._up_bnds = up_bnds
        self._low_bnds = low_bnds
        self._ranges = up_bnds - low_bnds
        self._p_mute = mutation_chance
        self._et_m = mutation_prm
    
    def chromosome_len(self):
        return len(self._low_bnds)
    
    def normalize(self, values):
        """Turn a bunch of gene values into values in range [0,1) normalized by the 
        allowed range for each gene.
        
        Arguments:
        values - a c by p array, for chromosome-length c and some p.
        """
        return (values - self._low_bnds)/self._ranges
        
    def denormalize(self, values):
        """Turn a bunch of values in the [0,1) range to be between the upper and 
        lower bounds of the creature's genotype.
        
        Arguments:
        vaues - a c by p array, for chromosome-length c and some p.
        """
        return values*self._ranges + self._low_bnds
    
    def denormalize_one(self, value, gene):
        """Turn one value in the range [0,1) to a value in the range allowed by 
        one selected gene.
        
        Arguments:
        value - the value to denormalize, in [0,1) range.
        gene - the index into the chromosome where the value is to be found.
        """
        return value*self._ranges[gene] + self._low_bnds[gene]
        
    def gen_population(self, num_subjects):
        pop = random.random_sample((num_subjects, self.chromosome_len()))
        return self.denormalize(pop)
    
    def breed(self, mama, papa):
        """Creates a new genome for a subject by recombination of parent genes, and
        possibly mutation of the result, depending on the creature's mutation 
        resistance.
        
        Arguments: 
        mama, papa - two breeders from the population, each a vector of genes.
        
        Returns: 
        offspring - a length-g array with the genotype of the offspring after 
            recombination and mutation.
        """
        # Recombination place, using one-point crossover:
        recomb = N.random.random_integers(0, self.chromosome_len() - 1)
        offspring = N.r_[mama[:recomb], papa[recomb:]]
        
        # Possibly mutate:
        which_genes = N.random.rand(self.chromosome_len()) < self._p_mute
        if not any(which_genes):
            return offspring
        
        # Polinomial mutation, see ref. [1]
        delta = 1 - N.array((offspring[which_genes] - self._low_bnds[which_genes], \
            self._up_bnds[which_genes] - offspring[which_genes])).min(axis=0) \
            / self._ranges[which_genes]
            
        mute_bases = N.random.random_sample(which_genes.sum())
        perturb = N.empty_like(mute_bases)
        close = mute_bases <= 0.5
        
        perturb[close] = (2*mute_bases[close] + \
            (1 - 2*mute_bases[close])*delta[close]**(self._et_m + 1))**(1./(self._et_m + 1)) - 1
        perturb[~close] = 1 - (2*(1 - mute_bases[~close]) + \
            2*(mute_bases[~close] - 0.5)*delta[~close]**(self._et_m + 1))**(1./(self._et_m + 1))
        offspring[which_genes] += perturb*self._ranges[which_genes]
        
        return offspring
