Linearized Solution
===================

The original code solving the linearized 1+1 relativistic multifluid problem, as presented in the paper
"The nonlinear development of the relativistic two-stream instability", Hawke, Comer, Andersson (arXiV:1303.4070),
was an incredibly messy mix of Matlab scripts (parts of which were auto-generated using Maple).

This is a re-implementation using python. The IPython notebooks present examples of how the solution can be used to 

1. test the stability of the system,
2. compute the spacetime solution,
3. compute the time-frequency solution, and adjust for the effect of the numerical differencing in time evolution codes.

These notebooks can take a while to run - at first it will be more effective to view the static versions via
http://nbviewer.ipython.org/github/IanHawke/MultiFluid1d/blob/master/python/Figure1.ipynb
http://nbviewer.ipython.org/github/IanHawke/MultiFluid1d/blob/master/python/InstabilityAnimation.ipynb
http://nbviewer.ipython.org/github/IanHawke/MultiFluid1d/blob/master/python/InstabilityAnimation.ipynb

At the moment these scripts really are a direct re-implementation, and need

1. the auto-generated code to be produced from sympy codegen, rather than via 3 levels of indirection from Maple;
2. speeding up (considerably);
3. Putting in a proper pythonic form.
