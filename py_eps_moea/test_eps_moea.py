# Test suite for the epsilon-moea algorithm implementation.

import unittest
from numpy import testing
from numpy import array, vstack, hstack, empty, zeros, r_, c_, bool
from eps_moea import *

class TestEpsMOEA(unittest.TestCase):
    def test_pareto_front(self):
        """Pareto front is found correctly"""
        fit1 = empty((10, 2))
        fit1[:,0] = r_[1:11]
        fit1[:,1] = 1./fit1[:,0]
        self.failUnless(pareto_front(fit1).all())

        fit2 = fit1 + 0.1
        fit2[:,1] = 2./fit2[:,0]
        self.failUnless(pareto_front(fit2).all())
        
        pf = pareto_front(vstack((fit1, fit2)));
        self.failUnless(pf[:10].all() and not pf[10:].any())

    def test_pop_accept(self):
        """Acceptence to the population is as expected"""
        fit = empty((10,2))
        fit[:,0] = r_[1:11]
        fit[:,1] = 1./fit[:,0]

        # Test 1: (0,0) fitness dominates everyone, so the first subject will be 
        # replaced:
        contend_fit = r_[0, 0]
        repl = pop_accept(fit, contend_fit)
        self.failIf(repl is None)
        
        # Test 2: (100,100) is dominated by all, it is thrown away.
        contend_fit = r_[100, 100]
        repl = pop_accept(fit, contend_fit)
        self.failUnless(repl is None)

        # Test 3: (1.5, 0.4) dominates the second element, so it is replaced.
        contend_fit = r_[1.5, 0.4]
        repl = pop_accept(fit, contend_fit)
        self.failUnless(repl == 1)
        
    def test_archive_accept(self):
        """Acceptance to the archive works"""
        # Sample fitness made of a population along two curves:
        # 1. x*y = 1
        # 2. x*y = 2
        
        fit = empty((20,2))
        fit[::2,0] = r_[1:11]
        fit[1::2,0] = r_[1.1:11.1]
        fit[::2,1] = 1./fit[::2,0]
        fit[1::2,1] = 2./fit[1::2,0]

        archive = pareto_front(fit)
        
        #For the epsilon-dominance tests:
        fine_grid = r_[0.01, 0.01]
        med_grid = r_[0.1, 0.1]
        coarse_grid = r_[0.5, 0.5]
        grids = c_[fine_grid, med_grid, coarse_grid]

        # Test 1: (0,0) is better than everyone, it should be the only one in the 
        # new archive.
        contend_fit = r_[0, 0]
        repl = pop_accept(fit, contend_fit)

        for grid in grids.T:
            new_arch = archive.copy()
            eps_fit = fit - N.fmod(fit, grid)
            eps_cont = contend_fit - N.fmod(contend_fit, grid)
            accepted = archive_accept(new_arch, fit, eps_fit, contend_fit, eps_cont)
            
            self.failUnless(accepted)
            # Since in this test the contender is not in the population,
            # the entire archive is emptied.
            self.failIf(new_arch[1:].any())

        # Test 2: (2.1, 0.51) is dominated by the archive, the archive must return 
        # unchanged.
        contend_fit = r_[2.1, 0.51]
        repl = pop_accept(fit, contend_fit)

        for grid in grids.T:
            new_arch = archive.copy()
            eps_fit = fit - N.fmod(fit, grid)
            eps_cont = contend_fit - N.fmod(contend_fit, grid)
            accepted = archive_accept(new_arch, fit, eps_fit, contend_fit, eps_cont)
            self.failIf(accepted);
            testing.assert_array_equal(new_arch, archive)

        # Test 3: (1.5, 0.4) dominates the second element of the front, so it is 
        # replaced in the archive.
        contend_fit = r_[1.5, 0.4]
        repl = pop_accept(fit, contend_fit)

        for grid in grids[:,:2]:
            new_arch = archive.copy()
            eps_fit = fit - N.fmod(fit, grid)
            eps_cont = contend_fit - N.fmod(contend_fit, grid)
            accepted = archive_accept(new_arch, fit, eps_fit, contend_fit, eps_cont)
            
            self.failUnless(accepted)
            # Only the place of the second item of the front (index 2 in the 
            # population) is changed.
            self.failIf(new_arch[2])
            new_arch[2] = True
            testing.assert_array_equal(new_arch, archive)

        # on the coarse grid solutions are lost:
        eps_fit = fit - N.fmod(fit, coarse_grid)
        eps_cont = contend_fit - N.fmod(contend_fit, coarse_grid)
        accepted = archive_accept(archive, fit, eps_fit, contend_fit, eps_cont)
        self.failUnless(accepted)
        testing.assert_array_equal(archive, array(([True] + [False]*19)))
