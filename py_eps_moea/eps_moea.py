import numpy as N
from scipy import linalg as LA

def distances(population, compare_to=None):
    """calculates normalized euclidean distance between each two elements in the 
    population.
    
    Arguments: 
    population - the points in the decision space for which to find distances,
        normed to the ranges allowed for each variable.
    compare_to - optional element of the population that we compare only to him.
    
    Returns: a square matrix of distance from individual i to individual j, or a 
        vector of distances to one element if compare_to is given.
    """
    if compare_to is not None:
        diff = (compare_to - population)**2
        return diff.sum(axis=1)
    
    # Otherwise, compare everyone to everyone else:
    dist_func = N.vectorize(lambda ii, jj: ((population[ii] - population[jj])**2).sum())
    return N.fromfunction(dist_func, (population.shape[0], population.shape[0]))

def pareto_front(fitness):
    """Given multiple fitness values of a population, find which individuals are in 
    the pareto-front(i.e. the non-dominated subset of the population.
    Assumes the problem is in canonical form, i.e. the target functions are all to 
    be minimized.
    
    Arguments:
    fitness - a p by t array, where fitness(p,t) is the value of function t for 
        individual p.
    
    Returns:
    pf - a boolean vector of length p, saying which individual is in the front.
    """
    pf = N.empty(fitness.shape[0], dtype=N.bool)
    
    for subject in xrange(fitness.shape[0]):
        others = fitness[N.r_[:subject, (subject+1):fitness.shape[0]]]
        non_dominated = (fitness[subject] < others).any(axis=1) | \
            (fitness[subject] <= others).all(axis=1)
        pf[subject] = non_dominated.all()
        
    return pf

def pop_select(fitness):
    """Randomly mates two subjects from a population and selects the dominating
    subject. If no dominance exists, select the second one (this is arbitrary).
    
    Arguments: 
    fitness - a p by t array, where fitness(p,t) is the value of function t for 
        individual p.
    
    Returns: 
    the index of the subject selected for breeding.
    """
    compete = N.random.random_integers(0, fitness.shape[0] - 1, size=2)
    if (fitness[compete[0]] <= fitness[compete[1]]).all() and \
        (fitness[compete[0]] < fitness[compete[1]]).any():
        return compete[0]
    else:
        return compete[1]

def archive_select(archive_marker):
    """Selects for breeding an individual from the archive. Currently selects 
    randomly.
    
    Arguments:
    archive_marker - a length-p boolean vector stating which of the population
        is in the archive.
    
    Returns:
    the index in the population of the selected archive member.
    """
    archive_size = archive_marker.sum()
    return N.where(archive_marker)[0][N.random.random_integers(0, archive_size - 1)]

def pop_accept(fitness, contend_fit):
    """Find some poor underdog in the population that is pareto-dominated by an 
    up-and-coming contender - if any such underdog exists. Otherwise cast away
    the contender.
    
    Arguments:
    fitness - a p by t matrix for p subjects and t target functions. fitness(i,j)
        is the value of target function j for subject i.
    contend_fit - a 1 by t vector with the contender's fitness.
    
    Returns:
    repl - the index of the replaced subject, or None if no replacement occurred
    """
    underdogs = (contend_fit <= fitness).all(axis=1) & (contend_fit < fitness).any(axis=1)
    if not underdogs.any():
        return None
    return N.where(underdogs)[0][0]

def archive_accept(archive, fitness, contend_fit, contend_idx, grid):
    """Update the pareto-front in the archive using epsilon-dominance.
    
    Arguments:
    archive - a p by 1 boolean array stating which individual is in the archive.
        The archive is modified to the new archive IN-PLACE!
    fitness - a p by t array whose element (i,j) is the value of the jth fitness
        function of t functions, for the ith individual of p individuals, after the 
        last iteration of the population acceptance method.
    contend_fit - the fitness of the contender.
    contend_idx - the index in the population of the contender to replace some of 
        the archive.
    grid - a vector of length t specifying the edge-length of the n-dim grid in 
        each of the fitness space axes.
    
    Returns:
    accepted - a boolean value stating whether the contender was accepted.
    """
    # We have more than one level of selection from the population,
    # so its convenient to work with indices rather than boolean arrays.
    archive_idxs = N.where(archive)[0]
    accepted = 1
    
    # Calculate the grid-fitness of the population:
    archive_size = sum(archive)
    grid_fit = fitness[archive] - N.fmod(fitness[archive], grid)
    # and of the contender:
    grid_cont = contend_fit - N.fmod(contend_fit, grid)
     
    # Instead of Matlab's unique(..., 'rows') we use python's sets.
    for vertex in map(N.array, set(tuple(row) for row in grid_fit)):
        if (grid_cont == vertex).all():
            # This hypercube may have non-dominated members.
            in_grid = N.where((grid_fit == vertex).all(axis=1))[0];
            
            top_dogs = (fitness[archive_idxs[in_grid]] <= contend_fit).all(axis=1) and \
                (fitness[archive_idxs[in_grid]] < contend_fit).any(axis=1);
            if top_dogs.any():
                return 0
            
            underdogs = (contend_fit <= fitness[archive_idxs[in_grid]]).all(axis=1) and \
                (rep_cond_fit < fitness[archive_idxs[in_grid]]).any(axis=1);
            
            archive[archive_idxs[in_grid[underdogs]]] = False;
        
        elif (grid_cont >= vertex).all() and (grid_cont > vertex).any():
            # The contender is dominated. The function can end, because
            # it's impossible that others in the archive may still be dominated.
            archive[contend_idx] = False; # make sure no leftovers from last iteration.
            return 0
            
        elif (grid_cont <= vertex).all() and (grid_cont < vertex).any():
            # This hypercube is dominated by the new solution, exclude all its
            # members from the archive:
            in_grid = N.where((grid_fit == vertex).all(axis=1))[0][0]
            archive[archive_idxs[in_grid]] = False
            archive[contend_idx] = True
        
    return accepted
        
def eps_moea_optimize(creature, pop_size, niche_size, conv_gens, num_gens, objectives, grid):
    """Run an optimization using the epsilon-moea algorithm.
    
    Arguments: 
    creature - a description of the problem in a Creature object.
    pop_size - number of individuals to use for the simulation. The more you use,
        The better is the pareto front, but it's slower and has diminishing returns.
    niche_size - the maximum distance between asubjects that generates fitness 
        punishment (distance in the normed genotype space).
    conv_gens - after this many iterations with no change in the archive population, 
        the iteration stops.
    nun_gens - maximum total number of iterations, with or without convergence.
    objectives - a function that given a population array returns the fitness array.
    grid - the size of the hypercubes in the epsilon-dominance tests.
    
    Returns:
    population - the population array after the latest iteration.
    fitness - the fitness array after the latest iteration.
    archive - the final archive.
    """
    population = creature.gen_population(pop_size)
    niche_size = niche_size**2
    archive_stagnation = 0
    normed_pop = creature.normalize(population)
    dist = distances(normed_pop)
    dist[dist == 0] = N.inf

    # Initial fitness:
    fitness = objectives(population)
    archive = pareto_front(fitness)
    
    # Inflict the cholera on crowded populations:
    crowded = dist < niche_size;
    cholera = (crowded*dist/niche_size).sum(axis=1) # quadratic niching.
    niched_fit = fitness - cholera[:,None];
    
    while (archive_stagnation < conv_gens) and (num_gens > 0):
        # Generate new solution:
        mama = pop_select(niched_fit);
        papa = archive_select(archive);
        offspring = creature.breed(population[mama], population[papa])
        
        # Offspring niching:
        dist_of = distances(normed_pop, offspring)
        
        # Accept the new solution to the population and archive:
        contend_fit = objectives(offspring)
        repl = pop_accept(fitness, contend_fit)
        if repl is None:
            archive_stagnation += 1
            num_gens -= 1
            continue
        
        population[repl] = offspring
        normed_pop[repl] = creature.normalize(offspring)
        dist[:,repl] = dist_of
        dist[repl] = dist_of
        dist[repl, repl] = N.inf
        
        # Offspring niching:
        # The offspring is punished for being close to existing solutions, but 
        # the existing ones are only judged relative to other existing solutions.
        crowded = dist_of < niche_size
        cholera = sum(crowded * dist_of / niche_size) # linear niching.
        niched_cont = contend_fit - cholera
        
        accepted = archive_accept(archive, niched_fit, niched_cont, repl, grid)
        
        # Prepare next iteration:
        fitness[repl] = contend_fit
        if accepted:
            archive_stagnation = 0
        else:
            archive_stagnation += 1
        num_gens -= 1
    
    return population, fitness, archive

if __name__ == "__main__":
    # A little test, with the first test function.
    import test_functions
    from creature import Creature
    
    grid = N.r_[0.01, 0.01]
    cr = Creature(N.zeros(30), N.ones(30), 0.1)
    population, fitness, archive = eps_moea_optimize(cr, 100, 0.43, 600, 25000, \
        test_functions.tau1, grid)
    
    import pylab as P
    archive = N.where(archive)[0]
    P.plot(fitness[archive,0], fitness[archive,1], 'o')
    P.title('Archive fitness')
    
    P.figure()
    P.plot(fitness[:,0], fitness[:,1], 'go')
    P.title('Population fitness')
    
    P.show()
