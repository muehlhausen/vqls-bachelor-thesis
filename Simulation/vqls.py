"""
################################################################################

Import all required libraries.

"""

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

# import the params object of the GlobalParameters class
from GlobalParameters import params


"""
This program provides an implementation of the Variational Quantum Linear Solver
as presented by Bravo-Prieto et al. It is implemented for the QISKIT Quantum
Simulator. This version of the program uses the local cost function.

Author: Alexander Cornelius Muehlhausen

#########################################################################

######### Core functions

The core functions of the VQLS algorithm.

"""


def generate_ansatz(alpha: List[float]) -> QuantumCircuit:
    """
    This function returns a circuit that implements V(alpha).
    """

    qr_ansatz = QuantumRegister(params.n_qubits)
    circ_ansatz = QuantumCircuit(qr_ansatz)

    if not any(isinstance(item, list) for item in alpha):
        # this will reformat the list alpha to the required format (if needed)
        # this is necessary as the minimizer returns a list without sublists
        alpha = _format_alpha(alpha)

    # 0th sublayer
    for qubit in range(0, params.n_qubits):
        circ_ansatz.ry(alpha[0][qubit], qr_ansatz[qubit])

    if params.n_qubits % 2 == 0:
        # all other sublayers
        for sublayer in range(1, 2 * params.n_layers, 2):
            # first sublayer of the current layer
            # controlled Z-Gate pairs
            for qubit_a, qubit_b in zip(qr_ansatz[::2], qr_ansatz[1::2]):
                circ_ansatz.cz(qubit_a, qubit_b)
            for rotation_param, qubit in zip(alpha[sublayer], qr_ansatz):
                circ_ansatz.ry(rotation_param, qubit)

            # second sublayer of the current layer
            for qubit_a, qubit_b in zip(qr_ansatz[1::2], qr_ansatz[2::2]):
                circ_ansatz.cz(qubit_a, qubit_b)
            # and Ry to each qubit except the first and last one
            for rotation_param, qubit in zip(alpha[sublayer+1], qr_ansatz[1::]):
                circ_ansatz.ry(rotation_param, qubit)

    else:
        # all other sublayers
        for sublayer in range(1, 2 * params.n_layers, 2):
            # first sublayer of the current layer
            for qubit_a, qubit_b in zip(qr_ansatz[:params.n_qubits-1:2],
                                        qr_ansatz[1:params.n_qubits-1:2]):
                circ_ansatz.cz(qubit_a, qubit_b)
            for rotation_param, qubit in zip(alpha[sublayer],
                                             qr_ansatz[:params.n_qubits-1:]):
                circ_ansatz.ry(rotation_param, qubit)

            # second sublayer
            for qubit_a, qubit_b in zip(qr_ansatz[1::2], qr_ansatz[2::2]):
                circ_ansatz.cz(qubit_a, qubit_b)
            for rotation_param, qubit in zip(alpha[sublayer+1], qr_ansatz[1::]):
                circ_ansatz.ry(rotation_param, qubit)

    return circ_ansatz


def hadamard_test(
    *, ansatz: Union[Gate, Operator] = None, first: Union[Gate, Operator]
    = None, first_uncontrolled: Union[Gate, Operator] = None, j: int
    = None, second_uncontrolled: Union[Gate, Operator] = None,
    second: Union[Gate, Operator] = None, im=None
     ):
    """
    This function returns a circuit with the Hadamard Test implemented.
    """

    # prep of QuantumCircuit
    qr_had = QuantumRegister(params.n_qubits)
    circ_had = QuantumCircuit(qr_had)
    ancilla = QuantumRegister(1, name="ancilla")
    cr_had = ClassicalRegister(1)

    circ_had.add_register(ancilla)
    circ_had.add_register(cr_had)

    qubits_designation_control = [i for i in range(params.n_qubits)]
    qubits_designation_control.insert(0, params.n_qubits)

    def append_ifExists(obj: Union[Gate, Operator], control=False):
        if isinstance(obj, (Gate, Operator)):
            _obj = obj.copy()
            if isinstance(_obj, Operator):
                _obj = _obj.to_instruction()
            if control is True:
                _obj = _obj.control(1)
                circ_had.append(_obj, qubits_designation_control)
            else:
                circ_had.append(_obj, qr_had)

    # act on the ancilla
    circ_had.h(ancilla)

    # if Im(<>) shall be calculated
    if im is not None:
        circ_had.sdg(ancilla)

    append_ifExists(ansatz)

    append_ifExists(first, True)

    append_ifExists(first_uncontrolled)

    if j is not None:
        circ_had.cz(params.n_qubits, qr_had[j])

    append_ifExists(second_uncontrolled)

    append_ifExists(second, True)

    # last operation on the ancilla & measurement
    circ_had.h(ancilla)

    circ_had.measure(ancilla, cr_had)

    return circ_had


def calculate_beta(alpha: List[float]) -> List[List[complex]]:
    # preparation of the result list
    beta = [[complex(0, 0) for _ in params.coefficients] for _
            in params.coefficients_conjugate]
    # generate Ansatz outside the loops for better performance
    V = generate_ansatz(alpha).to_gate()

    for gate_l, (l, coeff_l) in zip(params.decomposition_asGate,
                                    enumerate(params.coefficients)):
        if coeff_l == 0:
            continue
        for gate_m_adj, (m, coeff_m) in zip(params.decomposition_adjoint_asOperator,
                                            enumerate(params.coefficients_conjugate)):
            if coeff_m == 0:
                continue
            # circuit for Re ( <0| V(alpha)(+) A_m(+) A_l V(alpha) |0>)
            circ_had_re = hadamard_test(ansatz=V, first=gate_l,
                                        second=gate_m_adj)
            # circuit for Im ( <0| V(alpha)(+) A_m(+) A_l V(alpha) |0>)
            circ_had_im = hadamard_test(ansatz=V, first=gate_l,
                                        second=gate_m_adj, im=1)

            # calculate Re and Im
            expV_had_re = _calculate_expectationValue_HadamardTest(circ_had_re)
            expV_had_im = _calculate_expectationValue_HadamardTest(circ_had_im)
            # piece together <> from Re and Im
            expV_had = complex(expV_had_re, expV_had_im)
            beta[l][m] = expV_had

    return beta


def calculate_delta(alpha: List[float], beta: List[List[complex]]) -> List[List[complex]]:
    delta = [[complex(0, 0) for _ in params.coefficients] for _
             in params.coefficients_conjugate]

    V = generate_ansatz(alpha).to_gate()

    U = _U_primitive().to_gate()

    U_dagger = U.copy()
    U_dagger = Operator(U_dagger)
    U_dagger = U_dagger.adjoint()

    for gate_l, (l, coeff_l) in zip(params.decomposition_asGate,
                                    enumerate(params.coefficients)):
        if coeff_l == 0:
            continue
        for gate_m_adj, (m, coeff_m) in zip(
                        params.decomposition_adjoint_asOperator,
                        enumerate(params.coefficients_conjugate)
                         ):
            if coeff_m == 0:
                continue
            temp = beta[l][m]
            # 1/n_qubits sum_j <0| V(+) A_m(+) U (Z_j * 1_{j-bar}) U(+) A_l V |0>
            for j in range(params.n_qubits):
                # Re(<0| V(+) A_m(+) U (Z_j * 1_{j-bar}) U(+) A_l V |0>)
                circ_had_re = hadamard_test(
                              ansatz=V, first=gate_l,
                              first_uncontrolled=U_dagger, j=j,
                              second_uncontrolled=U,
                              second=gate_m_adj)

                # Im(<0| V(+) A_m(+) U (Z_j * 1_{j-bar}) U(+) A_l V |0>)
                circ_had_im = hadamard_test(
                              ansatz=V, first=gate_l,
                              first_uncontrolled=U_dagger, j=j,
                              second_uncontrolled=U,
                              second=gate_m_adj, im=1)

                # calculate Re and Im
                expV_had_re = _calculate_expectationValue_HadamardTest(circ_had_re)
                expV_had_im = _calculate_expectationValue_HadamardTest(circ_had_im)
                # piece together <> from Re and Im
                expV_had = complex(expV_had_re, expV_had_im)

                temp += 1/params.n_qubits * expV_had

            delta[l][m] = temp

    return delta


def calculate_local_cost_function(alpha: List[float]) -> Union[complex, float]:
    """
    returns
    cost = <x| H_local |x> / <psi|psi>

    with
    <x| H_local |x> = Sum_l_m c_l c_m(*) delta
    <psi|psi> = Sum_l_m c_l c_m(*) beta
    """

    beta = calculate_beta(alpha)
    delta = calculate_delta(alpha, beta)
    xHx = 0
    psipsi = 0

    for l, coeff_l in enumerate(params.coefficients):
        for m, coeff_m_conj in enumerate(params.coefficients_conjugate):
            xHx += coeff_l * coeff_m_conj * delta[l][m]
            psipsi += coeff_l * coeff_m_conj * beta[l][m]

    # cost = xHx / psipsi
    cost = abs(xHx/psipsi)

    print(alpha)
    print("local cost function " + str(cost))

    return cost


def minimize_local_cost_function(method: str) -> List[float]:
    """
    minimizes the local cost function. returns the alpha for which the
    approximation
    A V(alpha_out) |0> approx |b>
    is optimal.

    It provides several methods to minimize the cost function:
    Locally using scipy.optimize.minimize
    Basinhopping
    Differential evolution
    SHGO
    dual annealing
    """
    # accepted methods for the standard minimization
    minimization_methods = ["Nelder-Mead", "Powell", "CG", "BFGS",
                            "Newton-CG", "L-BFGS-B", "TNC", "COBYLA",
                            "SLSQP", "trust-constr", "dogleg",
                            "trust-ncg", "trust-exact", "trust-krylov"]


    # use minimize and the methods to find a local minimum
    if method in minimization_methods:
        min = minimize(calculate_local_cost_function, x0=params.alpha_0,
                       method=method,
                       options={'maxiter': params.COBYLA_maxiter})

    # use basinhopping to find the global minimum
    # use method_minimization as set in the section constants for the local
    # minimization
    elif method == "basinhopping":
        min = basinhopping(calculate_local_cost_function, x0=params.alpha_0,
                           niter=params.COBYLA_maxiter, minimizer_kwargs =
                           {'method': method_minimization})

    # global minimum via differnetial evolution
    elif method == "dif_evol":
        min = differential_evolution(calculate_local_cost_function, bounds =
                                     [(0, 2 * np.pi) for i in
                                     range(len(params.alpha_0))])

    # global minimum via shgo
    elif method == "shgo":
        min = shgo(calculate_local_cost_function, bounds = [(0, 2 * np.pi) \
                   for i in range(len(params.alpha_0))])

    # global minimum via dual annealing
    elif method == "dual_annealing":
        min = dual_annealing(calculate_local_cost_function, bounds = \
                             [(0, 2 * np.pi) for i in range(len(params.alpha_0))])

    else:
        print("no valid method given")
        return 0

    print(min)
    alpha_out = min['x']
    print(alpha_out)

    return alpha_out


"""
Fix the result.
"""


def postCorrection(qc: QuantumCircuit) -> QuantumCircuit:
    """
    This function is used to apply post correction to the circuit generated by
    applying V(alpha) and A as the result is not identical to |b>. The result of
    a circuit built using the functions presented above returns the reverse of b
    with random sign errors.
    """

    for i in range(params.n_qubits):
        qc.x(i)
        qc.z(i)
    return qc


"""
################################################################################

######### Helper functions

Functions that do not add to the logic of the algorithm but instead implement
often used features or need to be changed sometimes.

"""


def _format_alpha(alpha_unformated: List[float]) -> List[List[float]]:
    """
    This function formats a list to be in the correct form for the function that
    builds V(alpha). This means it will format the list to be of the form:
    [[0th sublayer], [1st sublayer], [2nd sublayer], ...]

    So for e.g. 4 qubits it will return (- stands for some random value):
    [[-,-,-,-],[-,-,-,-],[-,-],[-,-,-,-],[-,-]]

    and for e.g. 3 qubits:
    [[-,-],[-,-],[-,-],[-,-],[-,-]]
    """

    alpha_formated = []

    if any(isinstance(item, list) for item in alpha_unformated):
        return alpha_unformated

    else:
        if (params.n_qubits % 2) == 0:

            start = 0
            end = params.n_qubits
            alpha_formated.append(alpha_unformated[start:params.n_qubits])

            for _ in range(params.n_layers):
                start = end
                end = start + params.n_qubits
                alpha_formated.append(alpha_unformated[start:end])

                start = end
                end = start + params.n_qubits - 2
                alpha_formated.append(alpha_unformated[start:end])

        else:
            start = 0
            end = params.n_qubits
            alpha_formated.append(alpha_unformated[start:end])

            for _ in range(params.n_layers):
                start = end
                end = start + params.n_qubits-1
                alpha_formated.append(alpha_unformated[start:end])

                start = end
                end = start + params.n_qubits - 1
                alpha_formated.append(alpha_unformated[start:end])

        return alpha_formated


def _calculate_expectationValue_HadamardTest(circ_had: QuantumCircuit) -> float:
    """
    Will return the expectation value for a given circuit for a Hadamard test.
    Supports different backends.
    """

    if params.qiskit_simulation_backend == 'statevector_simulator':

        backend = Aer.get_backend('statevector_simulator')
        t_circ = transpile(circ_had, backend)
        qobj = assemble(t_circ)

        p_0 = 0
        p_1 = 0

        for _ in range(params.qiskit_simulation_shots):
            job = backend.run(qobj)
            result = job.result()
            counts = result.get_counts(circ_had)
            p_0 += counts.get('0', 0)
            p_1 += counts.get('1', 0)

        return ((p_0 - p_1) / params.qiskit_simulation_shots)

    if params.qiskit_simulation_backend == 'qasm_simulator':

        backend = Aer.get_backend('qasm_simulator')
        job = execute(circ_had, backend, shots=params.qiskit_simulation_shots)
        result = job.result()
        counts = result.get_counts(circ_had)

        number = counts.get('1', 0) + counts.get('0', 0)
        p_0 = counts.get('0', 0)/number
        p_1 = counts.get('1', 0)/number

        return (p_0 - p_1)


def _U_primitive() -> QuantumCircuit:
    """
    This function generates a circuit that resembles a U gate that fulfills:
    U |0> = |b> .

    This is primitive because this circuit calculates / defines |b> as:
    |b> = U |0> = A * H(all) * |0>

    H(all) stands for: Hadamard Gate applied to all qubits.
    """

    qr_U_primitive = QuantumRegister(params.n_qubits)
    circ_U_primitive = QuantumCircuit(qr_U_primitive)
    """
    for i in range(n_qubits):
        circ_U_primitive.h(i)
    """

    circ_U_primitive.h(0)
    #circ_U_primitive.h(1)
    #circ_U_primitive.h(2)
    circ_U_primitive.h(3)

    A_copy = params.A.copy()
    if isinstance(params.A, Operator):
        A_copy = A_copy.to_instruction()
    circ_U_primitive.append(A_copy, qr_U_primitive)

    return circ_U_primitive
