# This class represents a creature (problem) genotype (decision space)
#
# References:
# [1] Kalyanmoy Deb, An efficient constraint handling method for genetic 
# algorithms, 31 May 2000

from numpy import random
import numpy as N

class Creature(object):
    def __init__(self, low_bnds, up_bnds, mutation_chance, mutation_prm=20, \
        recomb_chance=1.0, recomb_prm=15.):
        if len(low_bnds) != len(up_bnds):
            raise ValueError("Upper bounds number not equal to lower-bounds number.")
        if (low_bnds > up_bnds).any():
            raise ValueError("Upper bounds must be > lower bounds.")
        
        self._up_bnds = up_bnds
        self._low_bnds = low_bnds
        self._ranges = up_bnds - low_bnds
        
        self._p_mute = mutation_chance
        self._et_m = mutation_prm
        
        self._p_recomb = recomb_chance
        self._et_c = recomb_prm
    
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
    
    def sbx(self, mama, papa):
        """Create an offspring using simulated binary crossover.
        
        Arguments: 
        mama, papa - two breeders from the population, each a vector of genes.
        
        Returns: 
        offsprings - a 2-tuple of length-g arrays, each with the genotype of an 
            offspring after recombination and mutation.
        """
        offsprings = [mama.copy(), papa.copy()]
        
        if random.rand() <= self._p_recomb:
            recombed = (random.rand(self.chromosome_len()) <= 0.5) & (abs(mama - papa) > 1e-14)
            num_recombed = recombed.sum()
            cum_dist = random.rand(2, num_recombed)
            
            parents = N.vstack((mama[recombed], papa[recombed]))
            p_spread = abs(mama[recombed] - papa[recombed])
            correction = 2 - (1 + N.vstack((parents.min(axis=0) - self._low_bnds[recombed], 
                self._up_bnds[recombed] - parents.max(axis=0)))/ \
                p_spread)**-(self._et_c + 1.)
            
            spread = N.empty((2, num_recombed))
            cont_bool = cum_dist <= 1./correction
            contracting = N.where(cont_bool)
            expanding = N.where(~cont_bool)
            spread[contracting] = \
                (correction[contracting]*cum_dist[contracting])**(1./(self._et_c + 1))
            spread[expanding] = \
                (2 - correction[expanding]*cum_dist[expanding])**(-1./(self._et_c + 1))
            
            pm = random.random_integers(0, 1, num_recombed)
            pm[pm == 0] = -1
            recombed_traits = 0.5*(mama[recombed] + papa[recombed] + pm*spread*p_spread)
            
            underflow = N.where(recombed_traits < self._low_bnds[recombed])
            recombed_traits[underflow] = N.tile(self._low_bnds[recombed],(2,1))[underflow]
            overflow = N.where(recombed_traits > self._up_bnds[recombed])
            recombed_traits[overflow] = N.tile(self._up_bnds[recombed],(2,1))[overflow]
            
            # Merge into not recombined traits from one of the parents:
            which_recombed = random.random_integers(0, 1, num_recombed)
            offsprings[0][recombed] = N.where(which_recombed, \
                recombed_traits[0], recombed_traits[1])
            offsprings[1][recombed] = N.where(which_recombed, \
                recombed_traits[1], recombed_traits[0])
            
        return offsprings
        
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
        offsprings = self.sbx(mama, papa)
        
        # Possibly mutate:
        for offspring in offsprings:
            which_genes = N.random.rand(self.chromosome_len()) <= self._p_mute
            if not any(which_genes):
                continue
        
            # Polinomial mutation, see ref. [1]
            delta = 1 - N.array((offspring[which_genes] - self._low_bnds[which_genes], \
                self._up_bnds[which_genes] - offspring[which_genes])) \
                / self._ranges[which_genes]
                
            mute_bases = N.random.random_sample(which_genes.sum())
            perturb = N.empty_like(mute_bases)
            close = mute_bases <= 0.5
            close_ind = N.where(close)
            far_ind = N.where(~close)
            
            perturb[close] = (2*mute_bases[close] + \
                (1 - 2*mute_bases[close])*delta[0,close_ind]**(self._et_m + 1))**(1./(self._et_m + 1)) - 1
            perturb[~close] = 1 - (2*(1 - mute_bases[~close]) + \
                2*(mute_bases[~close] - 0.5)*delta[1,far_ind]**(self._et_m + 1))**(1./(self._et_m + 1))
            offspring[which_genes] += perturb*self._ranges[which_genes]
            
            underflow = offspring < self._low_bnds
            offspring[underflow] = self._low_bnds[underflow]
            overflow = offspring > self._up_bnds
            offspring[overflow] = self._up_bnds[overflow]
            
        return offsprings
