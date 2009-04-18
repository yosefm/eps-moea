# This class represents a creature (problem) genotype (decision space)
from numpy import random
import numpy as N

class Creature(object):
    def __init__(self, low_bnds, up_bnds, mutation_chance):
        if len(low_bnds) != len(up_bnds):
            raise ValueError("Upper bounds number not equal to lower-bounds number.")
        if (low_bnds > up_bnds).any():
            raise ValueError("Upper bounds must be > lower bounds.")
        
        self._up_bnds = up_bnds
        self._low_bnds = low_bnds
        self._ranges = up_bnds - low_bnds
        self._p_mute = mutation_chance
    
    def chromosome_len(self):
        return len(self._low_bnds)
    
    def normalize(self, values):
        """Turn a bunch of gene values into values in range [0,1) normalized by the 
        allowed range for each gene.
        
        Arguments:
        vaues - a c by p array, for chromosome-length c and some p.
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
        # Recombination place:
        recomb = N.random.random_integers(0, self.chromosome_len() - 1)
        offspring = N.r_[mama[:recomb], papa[recomb:]]
        
        # Possibly mutate:
        if N.random.rand() > self._p_mute: 
            return offspring
        
        gene = N.random.random_integers(0, self.chromosome_len() - 1)
        # mutation is done by assigning a new value to the gene in the allowed range.
        offspring[gene] = self.denormalize_one(N.random.rand(), gene)
        return offspring
    
