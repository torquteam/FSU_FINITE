import ctypes
import numpy as np
import scipy.optimize as sp
import pickle

# Load in c finite library
lib = ctypes.CDLL('./cpylibrary.so')

# Define the argument types and return type for hartree_method
lib.hartree_method.argtypes = [
    ctypes.POINTER(ctypes.c_double),  # fin_couplings (array of 16 doubles)
    ctypes.c_int,                     # A (int)
    ctypes.c_int,                     # Z (int)
    ctypes.c_int,                     # iterations (int)
    ctypes.c_int,                     # gridsize (int)
    ctypes.c_int,                     # meson_iterations (int)
    ctypes.POINTER(ctypes.c_double),  # Observables (array of 7 doubles)
    ctypes.c_double,                  # convergence_help (double)
    ctypes.c_bool,                    # print_densities (bool)
    ctypes.c_bool                     # print_meson_fields (bool)
]
lib.hartree_method.restype = ctypes.c_int  # Return type (int)

# Define the argument types and return type for get_parameters
lib.get_parameters.argtypes = [
    ctypes.c_double, ctypes.c_double, ctypes.c_double, 
    ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double, 
    ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double,
    ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double),
    ctypes.c_bool, ctypes.c_int, ctypes.c_bool
]
lib.get_parameters.restype = ctypes.c_int

# Here I need to wrap the c function to use with lm algo
# Define the Python wrapper function
def call_hartree(fin_couplings, A, Z):
    # Unchanged variables
    print_densities = False 
    print_meson_fields = False
    meson_iterations = 3    
    gridsize = 401      
    iterations = 20  

    # Prepare the input and output arrays
    fin_couplings = np.array(fin_couplings, dtype=np.double)
    Observables = np.zeros(7, dtype=np.float64)
    
    count = 0
    exit_code = -1
    # Call the C function
    while (exit_code != 0):
        exit_code = lib.hartree_method(
            fin_couplings.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            A, Z, iterations, gridsize, meson_iterations,
            Observables.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            pow(1.1,count), print_densities, print_meson_fields
        )
        count = count + 1
        if (count > 7):
            break
    
    if (exit_code != 0):
        results = [0.0, 0.0, 0.0]
        return results

    results = []
    results.append(Observables[0])
    results.append(Observables[3])
    results.append(Observables[5])
    return results

# Here I need to wrap the c function to use with lm algo
# Define the Python wrapper function
# BA, p0, Jtilde, mstar, K, L, Ksym, zeta, xi, lambda_s, fw, fp, masses[4], fin_couplings[16], bool flag, int gd_sol_type, bool delta_coupling)
#bulks = [ms,BA,p0,mstar/m,K,J,L,zeta]
def bulks_to_params(bulks):
    ms = bulks[0]
    BA = bulks[1]
    p0 = bulks[2]
    mstar = bulks[3]
    K = bulks[4]
    J = bulks[5]
    L = bulks[6]
    Ksym = bulks[7]
    zeta = bulks[8]

    # Unchanged variables
    xi = 0.0
    lambda_s = 0.0
    fw = 0.0   
    fp = 0.0  
    masses = [ms,782.5,763.0,980.0]
    delta_coupling = True

    # Prepare the input and output arrays
    masses = np.array(masses, dtype=np.double)
    bulks = np.array(bulks, dtype=np.double)
    fin_couplings = np.zeros(16, dtype=np.double)
    
    # Call the C function
    lib.get_parameters(
        BA, p0, J, mstar*939.0, K, L, Ksym, zeta, xi, lambda_s, fw, fp, 
        masses.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        fin_couplings.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        True, 1, delta_coupling
    )
    return fin_couplings

# Set the Nuclei
A = [16,40,48,68,90,100,116,132,144,208]
Z = [8, 20,20,28,40,50, 50, 50, 62, 82]

# Import exp data
exp_data = np.loadtxt("dat_files/exp_data.txt")

# Set Initial start point
bulks = [500.0,-16.3,0.153,0.57,220.0,32.5,50.0,300.0,0.001]

# function to compute residuals
def residuals(bulks_arr, A, Z, exp_data):
    residuals = []
    couplings = bulks_to_params(bulks_arr)
    for i in range(len(exp_data)):
        y_model = call_hartree(couplings,A[i],Z[i])
        res = (y_model[0] - exp_data[i,0])/exp_data[i,1]
        residuals.append(res)
        if (exp_data[i,2] != -1):
            res = (y_model[1] - exp_data[i,2])/exp_data[i,3]
            residuals.append(res)
        if (exp_data[i,4] != -1):
            res = (y_model[2] - exp_data[i,4])/exp_data[i,5]
            residuals.append(res)
    print(bulks_arr)
    return np.array(residuals)

result = sp.least_squares(residuals,x0=bulks,method='lm',args=(A,Z,exp_data),diff_step=1e-4)

# Save the result to a file
with open('optimization_result_dino2.pkl', 'wb') as f:
    pickle.dump(result, f)
    
# Sample single hartree call
#fin_couplings = [93.5074, 151.6839, 200.5562, 0.0, 5.20326, -0.021739, 0.0007, 0.0, 0.047471, 0.0, 0.0, 0.0, 492.73, 782.5, 763.0, 980.0]
#observs = call_hartree(fin_couplings,A[0],Z[0],1.2)
#print(observs)

#RBM
#502.849  -16.295  0.1525  0.5941  248.354  33.443  61.030  0.0011395

#LM
#502.634  -16.302  0.1495  0.5996  260.543  33.386  64.295  0.0011256