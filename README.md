# vqls-bachelor-thesis
This repository contains the source code for a bachelor thesis at the [HISKP](https://www.hiskp.uni-bonn.de/index.php?id=29&L=1) of the [University of Bonn](https://www.uni-bonn.de/en). Topic was the implementation of the [VQLS-algorithm as presented by Bravo-Prieto et al.](https://arxiv.org/abs/1909.05820) in the [Qiskit](https://qiskit.org) quantum simulator. 

## Simulations
In the folder __Simulation__ an implementation for the simulator is stored. It consists of three files:
- `vqls.py` provides all functions implementing the algorithm. 
- `GlobalParameters.py` provides the class `GlobalParameters` which implements the data structure necessary to define the physical problem.
- `execute.py` is a short example on how to use the code for a simulation. 
For a simple demonstration of the code a user has to clone the repository and run `execute.py` on a system with `Python 3` and `Qiskit` installed. If another physical problem shall be simulated the code has to be adapted at two places:
1. In `GlobalParameters.py` the decomposition has to be inserted directly into the classes code, e.g. by altering the corresponding lines in `__init__.
2. In `execute.py`, when initialising `params` as an object of the class `GlobalParameters`, different parameters, e.g. other coefficients or another simulation backend can be selected. This is also the right place if one wants to continue working with the simulation's results.

