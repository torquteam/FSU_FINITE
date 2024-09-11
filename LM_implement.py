import ctypes
import numpy as np
import scipy.optimize as sp
import pickle
import math
import scipy.integrate as si
import numpy.polynomial.polynomial as poly
import sys

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
    ctypes.c_bool,                    # print_meson_fields (bool)
    ctypes.c_double                   # lgmr perturbation
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
def call_hartree(fin_couplings, A, Z, lgmr):
    # Unchanged variables
    print_densities = True 
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
            pow(1.1,count), print_densities, print_meson_fields,lgmr
        )
        count = count + 1
        if (count > 7):
            break
    
    if (exit_code != 0):
        results = [0.0, 0.0, 0.0]
        print("failed: ",fin_couplings)
        sys.exit()
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
    #zeta = bulks[8]
    zeta = bulks[7]

    # Unchanged variables
    Ksym = 60.0
    xi = 0.0
    lambda_s = 0.0
    fw = 0.0   
    fp = 0.0  
    masses = [ms,782.5,763.0,980.0]
    delta_coupling = False

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

def r2dens(A,Z,couplings):
    lgmr = [-0.0001,0.0001]
    intg = []
    hbar2 = (197.32698)**2
    Es = 13.2695

    dens = np.loadtxt(f"densities{A},{Z}.txt")
    r = dens[:,0]
    rho = dens[:,3] + dens[:,4]
    intg0 = si.simps(rho*r**4,x=r)
    for i in range(len(lgmr)):
        call_hartree(couplings,A,Z,lgmr[i])
        dens = np.loadtxt(f"densities{A},{Z}.txt")
        rho = dens[:,3] + dens[:,4]
        intg.append(si.simps(rho*r**4,x=r))
    
    lgmr.insert(1,0.0)
    intg.insert(1,intg0)
    lam = np.array(lgmr)*Es**3
    intg = np.array(intg)
    p = poly.Polynomial.fit(lam, intg, deg=2)
    p_deriv = p.deriv()
    der = p_deriv(lam[1])

    M1 = 8*math.pi*hbar2/939*pow(hbar2,-2)*intg0
    Mn1 = -2*math.pi*der/hbar2
    return math.sqrt(M1/Mn1)


# Set the Nuclei
A = [16,40,48,68,90,100,116,132,144,208]
Z = [8, 20,20,28,40,50, 50, 50, 62, 82]

# Import exp data
exp_data = np.loadtxt("dat_files/exp_data.txt")

# Set Initial start point
bulks = [500.0,-16.3,0.153,0.57,250.0,32.5,70.0,0.001] # FSU Models
#bulks = [500.0,-16.3,0.153,0.57,220.0,32.5,50.0,300.0,0.001] # DINO models

# function to compute residuals
def residuals(bulks_arr, A, Z, exp_data):
    residuals = []
    couplings = bulks_to_params(bulks_arr)
    print(couplings)
    for i in range(len(exp_data)):
        y_model = call_hartree(couplings,A[i],Z[i],0.0)
        res = (y_model[0] - exp_data[i,0])/exp_data[i,1]
        residuals.append(res)
        if (exp_data[i,2] != -1):
            res = (y_model[1] - exp_data[i,2])/exp_data[i,3]
            residuals.append(res)
        if (exp_data[i,4] != -1):
            res = (y_model[2] - exp_data[i,4])/(exp_data[i,5]/2) # cut uncertainty of form factors in half
            residuals.append(res)
        if (exp_data[i,6] != -1):
            GMR = r2dens(A[i],Z[i],couplings)
            res = (GMR - exp_data[i,6])/exp_data[i,7]
            residuals.append(res)
            print(GMR)
        
    return np.array(residuals)

result = sp.least_squares(residuals,x0=bulks,method='lm',args=(A,Z,exp_data),diff_step=1e-4)
# Save the result to a file
with open('optimization_result_FSU.pkl', 'wb') as f:
    pickle.dump(result, f)


#RBM
#502.849  -16.295  0.1525  0.5941  248.354  33.443  61.030  0.0011395

#LM
#502.634  -16.302  0.1495  0.5996  260.543  33.386  64.295  0.0011256
#couplings = bulks_to_params([ 5.06172335e+02, -1.61767181e+01,  1.48289332e-01,  5.54851251e-01, 2.51656487e+02,  3.32408658e+01,  7.17122557e+01,  1.02460375e-03])
#fin_couplings = [112.1996, 204.5469, 138.4701, 0.0, 1.4203, 0.023762, 0.06, 0.0, 0.030, 0.0, 0.0, 0.0, 491.500, 782.5, 763.0, 980.0]
#fin_couplings = [113.5072, 184.9651, 109.2054, 0.0, 3.5886, -0.015702, 0.001024, 0.0, 0.019583, 0.0, 0.0, 0.0, 506.1723, 782.5, 763.0, 980.0]
#observs = call_hartree(fin_couplings,A[9],Z[9],0.0)
#GMR = r2dens(A[9],Z[9],fin_couplings)
#print(GMR)

#0.0003575970862018177 -2.7699275428449803e-06

#208  82
#BA: -7.925
#Neutron Radius: 5.6961
#Proton Radius: 5.46605
#Charge Radius: 5.51505
#Weak Radius: 5.76953
#Rn - Rp: 0.230046
#Rwk - Rch: 0.254484
#Fch - Fwk: 0.0336714
#Rch: 5.53022

#[ 1.13507253e+02  1.84965110e+02  1.09205464e+02  0.00000000e+00
 # 3.58862159e+00 -1.57029606e-02  1.02460375e-03  0.00000000e+00
 # 1.95839797e-02  0.00000000e+00  0.00000000e+00  0.00000000e+00
 # 5.06172335e+02  7.82500000e+02  7.63000000e+02  9.80000000e+02]