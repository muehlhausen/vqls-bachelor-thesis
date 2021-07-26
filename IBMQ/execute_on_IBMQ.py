# all libraries used by some part of the VQLS-implementation
from qiskit import (
    QuantumCircuit, QuantumRegister, ClassicalRegister, IBMQ,
    Aer, execute, transpile, assemble
    )
from qiskit.circuit import Gate, Instruction
from qiskit.quantum_info.operators import Operator
from qiskit.extensions import ZGate, YGate, XGate, IGate

from scipy.optimize import (
                    minimize, basinhopping, differential_evolution,
                    shgo, dual_annealing
                    )

import random

import numpy as np
import cmath

from typing import List, Set, Dict, Tuple, Optional, Union


# import the params object of the GlobalParameters class
from GlobalParameters_IBMQ import params

# import the vqls algorithm and corresponding code
from vqls_IBMQ import (
    generate_ansatz,
    hadamard_test,
    calculate_beta,
    calculate_delta,
    calculate_local_cost_function,
    minimize_local_cost_function,
    postCorrection,
    _format_alpha,
    _calculate_expectationValue_HadamardTest,
    _U_YZ,
    _U_adjoint_YZ
    )

# The user input for the VQLS-algorithm.
# The decomposition for $A$ has to be manually
# inserted into the code of
# the class GlobalParameters.

print(
    "This program will execute the VQLS-algorithm on the real IBMQ-Manila quantum computer"
    + "with 4 qubits, 5 layers in the Ansatz and YZ (first Y. then Z) gate acting"
    + " on the third qubit (qubit_2).\n"
    + "|x_0> is defined by Hadamard gates acting on qubits 0, 2, 3."
)

# Executing the VQLS-algorithm
alpha_min = minimize_local_cost_function(params.method_minimization)

# Create a circuit for the vqls-result
qr_min = QuantumRegister(params.n_qubits)
circ_min = QuantumCircuit(qr_min)

# generate $V(\vec{alpha})$ and copy $A$
ansatz = generate_ansatz(alpha_min).to_gate()

# apply $V(\vec{alpha})$ and $A$ to the circuit
# this results in the state $\ket{b}$
circ_min.append(ansatz, qr_min)
circ_min.y(2)
circ_min.z(2)

# apply post correction to fix for sign errors and a "mirroring"
# of the result
circ_min = postCorrection(circ_min)

circ_ref = _U_YZ()

# simulate the circuit
backend = Aer.get_backend(
          'statevector_simulator')
t_circ = transpile(circ_min, backend)
qobj = assemble(t_circ)
job = backend.run(qobj)
result = job.result()


print(
    "This is the result of the execution on a real quantum computer.\n"
    + "Reminder: 4 qubits and a YZ-gate on the third qubit."
    + "|x_0> was defined by Hadamard gates acting on qubits 0, 2 and 3.\n"
    + "The return value of the minimizer (alpha_min):\n"
    + str(alpha_min)
    + "\nWith the result a circuit is build an simulated. The resulting statevector is,"
    + " with V(alpha_min) and A and the post correction applied:\n"
    + str(result.get_statevector())
)

t_circ = transpile(circ_ref, backend)
qobj = assemble(t_circ)
job = backend.run(qobj)
result = job.result()

print(
    "And this is the statevector for the reference circuit: A |x_0>\n"
    + str(result.get_statevector())
)
