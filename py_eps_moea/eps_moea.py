import numpy as N

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

def pop_select(fitness, archive):
    """Randomly mates two subjects from a population and selects the dominating
    subject. If no dominance exists, select the second one (this is arbitrary).
    
    Arguments: 
    fitness - a p by t array, where fitness(p,t) is the value of function t for 
        individual p.
    archive - the current archive population gets immunity, this variable is a 
        boolean vector saying which of the population is in the archive.
    
    Returns: 
    the index of the subject selected for breeding.
    """
    compete = N.random.random_integers(0, fitness.shape[0] - 1, size=2)
    
    if (fitness[compete[0]] <= fitness[compete[1]]).all() and \
        (fitness[compete[0]] < fitness[compete[1]]).any():
        return compete[0]
    if (fitness[compete[1]] <= fitness[compete[0]]).all() and \
        (fitness[compete[1]] < fitness[compete[0]]).any():
        return compete[1]
    if N.random.rand() < 0.5:
        return compete[0]
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
    return N.where(archive_marker)[0][\
        N.random.random_integers(0, archive_marker.sum() - 1)]

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
    if underdogs.any():
        return N.where(underdogs)[0][0]
        
    # Is he dominated?
    top_dogs = (contend_fit >= fitness).all(axis=1) & (contend_fit > fitness).any(axis=1)
    if top_dogs.any():
        return None
    
    # Non-domination: select at random.
    return N.random.random_integers(0, fitness.shape[0] - 1)
        
def archive_accept(archive, fitness, grid_fit, contend_fit, grid_cont):
    """Update the pareto-front in the archive using epsilon-dominance.
    
    Arguments:
    archive - a p by 1 boolean array stating which individual is in the archive.
        The archive is modified to the new archive IN-PLACE!
    fitness - a p by t array whose element (i,j) is the value of the jth fitness
        function of t functions, for the ith individual of p individuals, after the 
        last iteration of the population acceptance method.
    grid_fit - the epsilon-fitness value of the population (p by t)
    contend_fit - the fitness of the contender.
    grid_cont - a vector of length t with the epsilon-fitness of the contender
    
    Returns:
    accepted - a boolean value stating whether the contender was accepted.
    """
    # We have more than one level of selection from the population,
    # so its convenient to work with indices rather than boolean arrays.
    archive_idxs = N.where(archive)[0]
    arch_grid_fit = grid_fit[archive_idxs]
    arch_grid_fit.flags.writeable = False # For fast use of set objects
    arch_fit = fitness[archive_idxs]
    
    # We look for opportunities to reject the contender, along the way removing 
    # any dominated member of the archive:
    for vertex in map(N.frombuffer, set(row.data for row in arch_grid_fit)):
        high = (grid_cont > vertex).any()
        low = (grid_cont < vertex).any()
        if not(high or low):
            # This hypercube may have non-dominated members.
            in_grid = N.where((arch_grid_fit == vertex).all(axis=1))[0]
            
            lower_fit = (arch_fit[in_grid] < contend_fit).any(axis=1)
            higher_fit = (arch_fit[in_grid] > contend_fit).any(axis=1)
            if (lower_fit & ~higher_fit).any():
                return False

            archive[archive_idxs[in_grid]] = False
            
            underdogs = higher_fit & ~lower_fit
            if underdogs.all():
                break
            
            # Of the remaining solutions, the closest to the grid is taken:                
            dist_squares = ((arch_fit[in_grid[~underdogs]] - \
                arch_grid_fit[in_grid[~underdogs]])**2).sum(axis=1)
            remaining = dist_squares < ((contend_fit - grid_cont)**2).sum()
            
            if remaining.any():
                archive[archive_idxs[in_grid]] = True
                return False
        
        elif high and not low:
            # The contender is dominated. The function can end, because
            # it's impossible that others in the archive may still be dominated.
            return False
            
        elif low and not high:
            # This hypercube is dominated by the new solution, exclude all its
            # members from the archive:
            in_grid = N.where((arch_grid_fit == vertex).all(axis=1))[0]
            archive[archive_idxs[in_grid]] = False
    
    return True
        
def eps_moea_optimize(creature, pop_size, conv_gens, num_gens, objectives, grid):
    """Run an optimization using the epsilon-moea algorithm.
    
    Arguments: 
    creature - a description of the problem in a Creature object.
    pop_size - number of individuals to use for the simulation. The more you use,
        The better is the pareto front, but it's slower and has diminishing returns.
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
    archive_stagnation = 0

    # Initial fitness:
    fitness = objectives(population)
    archive = pareto_front(fitness)
    grid_fit = fitness - N.fmod(fitness, grid)

    while (archive_stagnation < conv_gens) and (num_gens > 0):
        # Generate new solution:
        mama = pop_select(fitness, archive);
        papa = archive_select(archive);
        offsprings = creature.breed(population[mama], population[papa])
        for offspring in offsprings:
            contend_fit = objectives(offspring)
            grid_cont = contend_fit - N.fmod(contend_fit, grid)
        
            # Accept the new solution to the population and archive:
            accepted = archive_accept(archive, fitness, grid_fit, contend_fit, grid_cont)
            repl = pop_accept(fitness, contend_fit)
            
            if accepted:
                archive_stagnation = 0
            if repl is None:
                continue
            
            # Prepare next iteration:
            population[repl] = offspring
            fitness[repl] = contend_fit
            archive[repl] = accepted
            grid_fit[repl] = grid_cont
        
        if not accepted: # any of the offsprings
            archive_stagnation += 1
        num_gens -= 1
    
    return population, fitness, archive

if __name__ == "__main__":
    # A little test, with the first test function.
    import test_functions
    from creature import Creature
    
    # Select which test to run from the commandline, default TZD1
    import sys
    testfun = test_functions.tau1
    if len(sys.argv) == 2:
        try:
            testfun = test_functions.tzd[int(sys.argv[1])]
        except IndexError:
            pass # Keep the default
    
    import time
    t = time.time()
    grid = N.r_[0.0075, 0.0075]
    cr = Creature(N.zeros(30), N.ones(30), 0.033)
    population, fitness, archive = eps_moea_optimize(cr, 100, 600, 20000, \
        testfun, grid)
    print time.time() - t
    
    import pylab as P
    archive = N.where(archive)[0]
    P.plot(fitness[archive,0], fitness[archive,1], 'o')
    P.title('Archive fitness')
    
    P.figure()
    P.plot(fitness[:,0], fitness[:,1], 'go')
    P.title('Population fitness')
    
    P.show()
