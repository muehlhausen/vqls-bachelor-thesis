# vqls-bachelor-thesis
This repository contains the source code for a bachelor thesis at the [HISKP](https://www.hiskp.uni-bonn.de/index.php?id=29&L=1) of the [University of Bonn](https://www.uni-bonn.de/en). Topic was the implementation of the [VQLS-algorithm as presented by Bravo-Prieto et al.](https://arxiv.org/abs/1909.05820) in the [Qiskit](https://qiskit.org) quantum SDK. 

Author of this thesis and implementation: Alexander C. Mühlhausen

Full documentation is provided within the thesis: `thesis-acm-vqls.pdf`.

## Simulations
In the folder __Simulation__ an implementation for the simulator is stored. It consists of three files:
- `GlobalParameters.py` provides the class `GlobalParameters` which implements the data structure necessary to define the physical problem. It also provides an instance of the class, called `params`, that is imported by all other files. This defines the actual problem the code works on.
- `vqls.py` provides all functions that implement the algorithm. 
- `execute.py` is a short example on how to use the code for a simulation. 

For a simple demonstration of the code a user has to clone the repository and run `execute.py` on a system with `Python 3.9` and `Qiskit 0.26.2` installed. If another physical problem shall be simulated the code has to be altered in two places:
1. In `GlobalParameters.py` an alternative decomposition has to be inserted directly into the classes code, e.g. by altering the corresponding lines in `__init__.` Furthermore, all other parameters (e.g. the coefficients) have to be altered in the initialization of the `params` object (that is imported by the other files).
2. In `vqls.py` the function `_U_primitive` is used to define the vector |x_0>. Thus, if one wants to simulate antoher physical problem, this code has to be altered as well.

## Execution on a quantum computer

In the folder __IBMQ__ an implementation for the execution on a real quantum computer, the 5 qubit `ibmq_manila` that is equipped with an [IBM Quantum](https://quantum-computing.ibm.com) `Falcon r5.11L` processor, is provided. The file structure is similar to the one for simulations, cf. to the thesis for a detailed comparison. The most important differences are:
1. A token to access IBM Quantum's systems has to be generated and inserted into the initialization of `params` in `GlobalParameters_IBMQ.py`.
2. When altering the decomposition in `GlobalParameters_IBMQ.py`, one also has to manually alter the code for `_U_Id` or similar functions in `vqls_IBMQ.py`.
