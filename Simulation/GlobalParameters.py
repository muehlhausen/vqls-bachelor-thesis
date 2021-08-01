from qiskit import (
    QuantumCircuit, QuantumRegister, ClassicalRegister,
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
executions it is very helpful to have them accesible and in one place.

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

    def _decompose_asOperator(self, gate: str) -> List[Operator]:
        """
        This function returns a list of Operator-Objects of length n_qubits.
        The first item is the Operator representing a Z gate acting on the
        first qubit. The last item is the Operator representing a Z gate
        acting on the last qubit.
        """

        if gate == "Z":
            Op = Operator(ZGate())
        elif gate == "X":
            Op = Operator(XGate())
        elif gate == "Y":
            Op = Operator(YGate())
        elif gate == "Id":
            Op = Operator(IGate())

        Id = Operator(IGate())
        dec_asOperator: List[Operator] = []

        # Operator acting on the last qubit
        temp1 = Op.copy()
        for _ in range(self.n_qubits - 1):
            temp1 = temp1.tensor(Id)
        dec_asOperator.insert(0, temp1)

        # all other qubits
        for x in range(self.n_qubits - 1):
            i = x+1
            temp = Id
            for j in range(i-1):
                temp = temp.tensor(Id)
            temp = temp.tensor(Op)
            for j in range(self.n_qubits - i - 1):
                temp = temp.tensor(Id)
            dec_asOperator.insert(0, temp)

        return dec_asOperator

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
            temp.to_gate()
            dec_asGate.append(temp)

        return dec_asGate

    def _decompositions(self):
        """
        This helper function is used to prepare some lists with standard decompositions
        that might be used to set together A.
        Those lists are stored in the class.
        """
        # Z Gates
        self.decomposition_asGate_Z = self._decompose_asGate("Z")
        self.decomposition_asOperator_Z = self._decompose_asOperator("Z")
        # [A_0(+), A_1(+), ...]
        self.decomposition_adjoint_asOperator_Z = [op.adjoint() for op in
                                                   self.decomposition_asOperator_Z]
        # Y Gates
        self.decomposition_asGate_Y = self._decompose_asGate("Y")
        self.decomposition_asOperator_Y = self._decompose_asOperator("Y")
        # [A_0(+), A_1(+), ...]
        self.decomposition_adjoint_asOperator_Y = [operator.adjoint() for operator in
                                                   self.decomposition_asOperator_Y]
        # X Gates
        self.decomposition_asGate_X = self._decompose_asGate("X")
        self.decomposition_asOperator_X = self._decompose_asOperator("X")
        # [A_0(+), A_1(+), ...]
        self.decomposition_adjoint_asOperator_X = [operator.adjoint() for operator in
                                                   self.decomposition_asOperator_X]

        # Identity Gates
        self.decomposition_asGate_Id = self._decompose_asGate("Id")
        self.decomposition_asOperator_Id = self._decompose_asOperator("Id")
        # [A_0(+), A_1(+), ...]
        self.decomposition_adjoint_asOperator_Id = [operator.adjoint() for operator in
                                                    self.decomposition_asOperator_Id]

    def __init__(self, n_qubits: int, n_layers: int, coefficients: List[complex],
                 COBYLA_maxiter: int = 150, qiskit_simulation_shots: int = 10**3,
                 qiskit_simulation_backend: str = 'qasm_simulator',
                 method_minimization: str = 'COBYLA'):
        self.n_qubits = n_qubits

        # alpha has form: [[0th sublayer], [1st sublayer], [2nd sublayer], ... ]
        self.n_layers = n_layers
        self.n_sublayers = 1 + 2 * self.n_layers

        self.coefficients = coefficients
        self.coefficients_conjugate = [np.conjugate(coeff) for coeff
                                       in self.coefficients]

        self.COBYLA_maxiter = COBYLA_maxiter
        self.qiskit_simulation_backend = qiskit_simulation_backend
        self.qiskit_simulation_shots = qiskit_simulation_shots
        self.method_minimization = method_minimization

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
        self.decomposition_asOperator = self.decomposition_adjoint_asOperator_Id
        self.decomposition_adjoint_asOperator = self.decomposition_adjoint_asOperator_Id

        # the Operator A
        self.A = self.coefficients[0] * self.decomposition_asOperator[0]
        for coeff, op in zip(self.coefficients[1:], self.decomposition_asOperator[1:]):
            self.A += coeff * op


# This instance of the class is imported and used by all other modules.
# Use it to model your physcal problem.
# To change A or its decomposition you have to modify the code of the class itself.

params = GlobalParameters(
            n_qubits=4,
            n_layers=4,
            coefficients=[complex(0, 0), complex(1, 0),
                          complex(0, 0), complex(0, 0)],
            COBYLA_maxiter=400,
            qiskit_simulation_shots=8192,
            qiskit_simulation_backend="qasm_simulator"
            )
