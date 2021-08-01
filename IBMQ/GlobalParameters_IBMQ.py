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


"""

######### Global Parameters

Those parameters are required by almost all functions but do not need to be
changed during runtime. As one might want to change them between program
executions it is sensible to have them accesible and in one place.

"""


class GlobalParameters:
    def _init_alpha_randomly(self) -> List[float]:
        """
        This function initalizes alpha randomly.
        """
        # how many parameters are required
        n_alpha_i = self.n_qubits + self.n_layers * (self.n_qubits-1) * 2

        alpha_rand = [float(random.randint(0, 6283))/1000 for _ in range(n_alpha_i)]

        return alpha_rand

    def _decompose_asGate(self, gate: str) -> List[Gate]:
        """
        This function returns a list of Gate-Objects of length n_qubits.
        The first item is the Gate applied to the first qubit, the last item
        is the Gate applied to the last qubit.
        """

        dec_asGate: List[Gate] = []
        for i in range(self.n_qubits):
            temp = QuantumCircuit(self.n_qubits)
            if gate == "Z":
                temp.z(i)
            elif gate == "Y":
                temp.y(i)
            elif gate == "X":
                temp.x(i)
            elif gate == "Id":
                temp.x(i)
                temp.x(i)
            elif gate == "ZY":
                temp.z(i)
                temp.y(i)
            elif gate == "YZ":
                temp.y(i)
                temp.z(i)
            temp.to_gate()
            dec_asGate.append(temp)

        return dec_asGate

    def _decompositions(self):
        """
        This helper function is used to prepare some lists with standard decompositions
        that might be used to set together A.
        Those lists are stored in the class.
        """
        # ZY Gates
        self.decomposition_asGate_YZ = self._decompose_asGate("YZ")
        self.decomposition_adjoint_YZ = self._decompose_asGate("ZY")

        # Identity Gates
        self.decomposition_asGate_Id = self._decompose_asGate("Id")
        self.decomposition_adjoint_Id = self.decomposition_asGate_Id

    def __init__(self, n_qubits: int, n_layers: int, coefficients: List[complex],
                 IBMQ_TOKEN: str, COBYLA_maxiter: int = 150, IBMQ_shots: int = 4096,
                 IBMQ_backend: str = 'statevector_simulator',
                 method_minimization: str = 'COBYLA',
                 ):
        self.n_qubits = n_qubits

        # alpha has form: [[0th sublayer], [1st sublayer], [2nd sublayer], ... ]
        self.n_layers = n_layers
        self.n_sublayers = 1 + 2 * self.n_layers

        self.coefficients = coefficients
        self.coefficients_conjugate = [np.conjugate(coeff) for coeff
                                       in self.coefficients]

        self.COBYLA_maxiter = COBYLA_maxiter
        self.method_minimization = method_minimization

        """
        Preparing IBMQ
        """

        self.IBMQ_load_account = IBMQ.enable_account(IBMQ_TOKEN)
        self.IBMQ_account = IBMQ.save_account(IBMQ_TOKEN)
        self.IBMQ_provider = IBMQ.load_account()
        self.IBMQ_backend = self.IBMQ_provider.backend.ibmq_manila
        self.IBMQ_shots = IBMQ_shots

        # initialize the parameters for the Ansatz randomly
        self.alpha_0 = self._init_alpha_randomly()

        self._decompositions()

        """
        Only things below this line require user interaction in the classes code
        (normally).
        """

        """
        The following lines present an example on how to define A and its decomposition.
        """
        self.decomposition_asGate = self.decomposition_asGate_Id
        self.decomposition_adjoint = self.decomposition_adjoint_Id


# This instance of the class is imported and used by all other modules.
# Use it to model your physcal problem.
# To change A or its decomposition you have to modify the code of the class itself.

params = GlobalParameters(
            n_qubits=4,
            n_layers=2,
            coefficients=[complex(0, 0), complex(0, 0),
                          complex(1, 0), complex(0, 0)],
            IBMQ_TOKEN="Insert your token here!",
            COBYLA_maxiter=100,
            IBMQ_shots=8192
            )
