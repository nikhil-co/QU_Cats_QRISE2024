"""
This script acts as a host for Q# to implement Quantum Approximate Optimization Algorithm (QAOA).
"""

# Importing libraries
# General imports
import time
from numpy import pi
import numpy as np
import matplotlib.pyplot as plt

# Libraries for Model Formulation
from docplex.mp.model import Model
from scipy.optimize import minimize

# Importing Q#
import qsharp

# Setting Q# root folder path
# This should be set to your Q# project root
# Use "/" backslash instead of "\" as path separator.
qsharp.init(project_root = 'C:/Users/londh/qc/QRISE_QRE/qaoa')

# Qiskit Imports
from qiskit_optimization.converters import QuadraticProgramToQubo
from qiskit_optimization.translators import from_docplex_mp

# Defining helper functions.
def find_most_common_solutions(input_dict, n):
    """
    Sorts the keys of the input dictionary in descending order based on their values and returns the first n keys.

    Parameters:
        input_dict (dict): A dictionary containing the keys and their corresponding values.
        n (int): The number of most common solutions to return.

    Returns:
        list: A list of the n most common keys sorted in descending order based on their values.
    """
    sorted_keys = sorted(input_dict, key=input_dict.get, reverse=True)

    return sorted_keys[:n]


# Building the Model
# We are using docplex to build the model and calculate the quadratics and linear terms.
def build_qubo(arr: list):
    """
    Function to build a QUBO (Quadratic Unconstrained Binary Optimization) model 
    for the Number Partitioning Problem (NPP).

    Args:
        arr (list): a list of integers representing the array from which the QUBO is built.
    Returns:
        a tuple containing the quadratic coefficients, linear coefficients, and the QUBO model.
    """
    # Length of the array - Length of the binary vector x
    n = len(arr)
    # Sum of the array - c
    c = sum(arr)

    # Building the model and its QUBO formulation.
    model = Model()
    x = model.binary_var_list(n)

    # Cost Function for Number Partitioning Problem (NPP)
    Q = (c - 2*sum(arr[i]*x[i] for i in range(n)))**2
    model.minimize(Q)   
    problem = from_docplex_mp(model)

    # QUBO formulation
    converter = QuadraticProgramToQubo()
    qubo = converter.convert(problem)

    # Quadratic and Linear Coefficients
    quadratics = qubo.objective.quadratic.coefficients
    linears = qubo.objective.linear.coefficients

    return quadratics, linears, qubo

def arr_to_str(a):
    """
    Converts a given array to a string representation.

    Args:
        a (list): The array to be converted.

    Returns:
        str: The string representation of the array.

    Example:
        >>> arr_to_str([1, 2, 3])
        '[1,2,3]'
    """
    string =''
    for i in a:
        string += str(i) + ","
    return '[' + string[:-1] + ']'

def integer_to_counts(n: int, result: list) -> dict:
    """
    Convert integers to counts and return a dictionary of counts.
    
    Args:
        n (int): the width of the binary representation
        result (list): a list of integers representing the results
    
    Returns:
        counts (dict):  - a dictionary containing the counts of each binary representation
    Example:
        >>> integer_to_counts(2, [2, 1, 2])
        {'00': 0, '01': 1, '10': 2, '11': 0}
    """
    counts = {}
    for i in range(2**n):
        counts[np.binary_repr(i,width=n)] = 0
    for integer in result:
        counts[np.binary_repr(integer,width=n)] += 1
    return counts

# Global variables for function calls and tracking cost.
func_call = 0
theta = []
cost = []

# Callback function
def callback_func(x):
    theta.append(x)
    return None

# QAOA Function for the Number Partitioning Problem
def qaoa_NPP(arr,layers:int):
    """
    Function implementing the QAOA algorithm for the Number Partitioning Problem.

    Args:
        arr (list): a list of integers.
        layers (int): the number of layers in the QAOA circuit.
    Returns:
        counts (dict): a dictionary containing the counts of each bitstring.

    """
    # Building the QUBO for NPP
    quadratics, linears, qubo = build_qubo(arr)
    num_qubits = len(arr)

    # Preparing the coefficients as flattened numpy arrays
    quadratics = quadratics.toarray().flatten()
    linears = linears.toarray()

    # Objective function we want to optimize
    def expectation_value(theta):
        # Global variables update
        global func_call 
        func_call = func_call + 1

        # Splitting the theta array into gammas and betas
        middle = int(len(theta)/2)
        gammas = theta[:middle]
        betas = theta[middle:]

        # Input string for the Q# function.
        input_str = f"{num_qubits},{layers},{arr_to_str(gammas)},{arr_to_str(betas)},{arr_to_str(quadratics)},{arr_to_str(linears[0])}"

        # Calling the Q# function
        int_results = qsharp.run(f"qaoa.circuit({input_str})",shots=100)

        # Converting the results to counts
        counts = integer_to_counts(num_qubits, int_results)

        # Calculating the expectation value for the best solution
        best_sol = max(counts, key=counts.get)
        exp =  qubo.objective.evaluate(np.array(list(best_sol), dtype='int'))
        
        # Appending to cost array
        cost.append(exp)
        
        # Printing the cost
        print(f'Function call: {func_call} - Cost: {exp}')

        return exp

    # Initial guess for gammas and betas
    init_gamma = np.array([pi / 3] * layers)
    init_beta = np.array([pi / 4] * layers)
    initial_guess = np.concatenate((init_gamma, init_beta))

    # Minimization of the objective function. We are using the COBYLA optimizer.
    start_time = time.time()
    res = minimize(expectation_value, initial_guess, method='COBYLA',callback=callback_func)
    end_time = time.time()
    elapsed_time = end_time - start_time

    print(f'\nElapsed time for QAOA: {elapsed_time} seconds')

    # Preparing the final gammas and betas
    middle = int(len(res.x)/2)
    prime_gammas = res.x[:middle]
    prime_betas = res.x[middle:]

    # Printing the optimal parameters
    print(f'Gammas: {prime_gammas}\nBetas: {prime_betas}\nFinal Cost: {res.fun}\n')

    # Preparing the input string with optimal gammas and betas 
    input_str = f"{num_qubits},{layers},{arr_to_str(prime_gammas)},{arr_to_str(prime_betas)},{arr_to_str(quadratics)},{arr_to_str(linears[0])}"

    # Calling the Q# function with optimal gammas and betas
    results = qsharp.run(f"qaoa.circuit({input_str})",shots=100)
    counts = integer_to_counts(num_qubits, results)
    
    return counts

if __name__ == "__main__":
    # Defining a test array and layers.
    test_array = [2,1,4]
    layers = 4

    # Running QAOA on for Number Partitioning.
    counts = qaoa_NPP(test_array,layers)

    # Plotting the optimal output state.
    plt.figure(figsize=(15, 5))
    plt.bar(range(len(counts)), list(counts.values()), align='center', color='red')
    plt.xticks(range(len(counts)), list(counts.keys()), rotation=90)
    plt.title("QAOA Output State with Optimal Parameters")
    plt.xlabel("Bit strings")
    plt.ylabel("Counts")
    plt.grid()
    plt.show()

    # Plotting Cost vs iterations.
    plt.figure(figsize=(15, 5))
    plt.plot(range(len(cost)),cost,color='g',ls='--',marker='o',lw=2)
    plt.xticks(range(1,len(cost)+1,5))
    plt.title('Cost vs. Iterations')
    plt.xlabel('Iterations')
    plt.ylabel('Cost')
    plt.grid()
    plt.show()

    # Printing Solutions Sets.
    best_sol = find_most_common_solutions(counts,3)
    print(f'\nTop 3 solutions for the array {test_array} and {layers} layers: \n{best_sol}')

    # Calculating S and S_A
    S = []
    S_A = []
    for ind,bit in enumerate(best_sol[0]):
        if bit == '1':
            S.append(ind)
        else:
            S_A.append(ind)

    # Calculating Sum(S) and Sum(S/A)
    sum_S = sum(np.array(test_array)[S])
    sum_S_A = sum(np.array(test_array)[S_A])

    # Printing the best optimal partition.
    print(f'\nBest Partition:\nS: {np.array(test_array)[S]}\nSum(S) = {sum_S}\n\nS/A: {np.array(test_array)[S_A]}\nSum(S/A) = {sum_S_A}')
