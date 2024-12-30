import ctypes
import numpy as np
import scipy.optimize as sp
import pickle
import math
import scipy.integrate as si
import numpy.polynomial.polynomial as poly
import sys
from pyhmc import hmc

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
    ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double,
    ctypes.c_double, ctypes.c_double, ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double),
    ctypes.c_bool, ctypes.c_int, ctypes.c_bool
]
lib.get_parameters.restype = ctypes.c_int

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
    
    count = 2
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

# Define the Python wrapper function
# BA, p0, Jtilde, mstar, K, L, Ksym, zeta, xi, lambda_s, fw, fp, masses[4], fin_couplings[16], bool flag, int gd_sol_type, bool delta_coupling)
#bulks = [ms,BA,p0,mstar/m,K,J,L,zeta]
def bulks_to_params(bulks):
    ms = bulks[0]*500.0
    BA = bulks[1]*(-16.3)
    p0 = bulks[2]*0.150
    mstar = bulks[3]*0.6
    K = bulks[4]*250.0
    J = bulks[5]*32.0
    L = bulks[6]*80.0
    zeta = bulks[7]*0.01
    #Ksym = bulks[8]*100.0
    #Gh2 = bulks[9]*2.0
    #fp = bulks[8]*(-100.0)
    #bIV = bulks[7]*0.5
    #Gt2 = bulks[9]*0.5
    #xi = bulks[9]*0.5

    
    # Unchanged variables
    Ksym = 15
    xi = 0.0
    lambda_s = 0.0
    fw = 0.0
    fp = 0.0
    Gt2 = 0.0
    #L = 60.0
    Gh2 = 0.05
    bIV = 0.0    
    masses = [ms,782.5,763.0,980.0]
    delta_coupling = False

    # Prepare the input and output arrays
    masses = np.array(masses, dtype=np.double)
    bulks = np.array(bulks, dtype=np.double)
    fin_couplings = np.zeros(19, dtype=np.double)
    
    # Call the C function
    lib.get_parameters(
        BA, p0, J, mstar*939.0, K, L, Ksym, zeta, xi, lambda_s, fw, fp, 
        Gt2, Gh2, bIV, masses.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        fin_couplings.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        True, 1, delta_coupling
    )
    return fin_couplings

# Algorithm to compute the GMR
def r2dens(A,Z,couplings):
    lgmr = [0.001,0.003,0.005]
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
    
    zero_loc = 0
    lgmr.insert(zero_loc,0.0)
    intg.insert(zero_loc,intg0)
    lam = np.array(lgmr)*Es**3
    intg = np.array(intg)
    p = poly.Polynomial.fit(lam, intg, deg=1)
    p_deriv = p.deriv()
    der = p_deriv(lam[zero_loc])

    M1 = 8*math.pi*hbar2/939*pow(hbar2,-2)*intg0
    Mn1 = -2*math.pi*der/hbar2
    return math.sqrt(M1/Mn1)

# function to compute residuals
def residuals(bulks_arr, A, Z, exp_data):
    residuals = []
    couplings = bulks_to_params(bulks_arr)
    print(bulks_arr)
    for i in range(len(exp_data)):
        y_model = call_hartree(couplings,A[i],Z[i],0.0)
        res = (y_model[0] - exp_data[i,0])/exp_data[i,1]
        residuals.append(res)
        if (exp_data[i,2] != -1):
            res = (y_model[1] - exp_data[i,2])/exp_data[i,3]
            residuals.append(res)
        if (exp_data[i,4] != -1):
            res = (y_model[2] - exp_data[i,4])/(exp_data[i,5]/4) # cut uncertainty of form factors by four
            residuals.append(res)
        if (exp_data[i,6] != -1):
            GMR = r2dens(A[i],Z[i],couplings)
            res = (GMR - exp_data[i,6])/exp_data[i,7]
            residuals.append(res)
            print(GMR)
        
    return np.array(residuals)

def logprob(bulks_arr, A, Z, exp_data):
    res = residuals(bulks_arr,A,Z,exp_data)
    logp = 0.5*np.sum(res**2)
    return logp

# Set the Nuclei
A = [16,40,48,68,90,100,116,132,144,208]
Z = [8, 20,20,28,40,50, 50, 50, 62, 82]

# Import exp data
exp_data = np.loadtxt("dat_files/exp_data.txt")

# Set Initial start point
conv = [500.0, -16.3, 0.150, 0.6, 250.0, 32.0, 80.0, 0.01]
bulks = [1.0 , 1.0  , 1.0  , 1.0, 0.88  , 1.0 , 1.0 , 1.0]

# Run Calibration and Save Results
# result = sp.least_squares(residuals,x0=bulks,method='lm',args=(A,Z,exp_data),diff_step=1e-4)
# with open('Isovector_tensor.pkl', 'wb') as f:
#    pickle.dump(result, f)

# Unpack results
with open('Isovector_tensor.pkl', 'rb') as f:
    result = pickle.load(f)
print(result.cost)
bulks = np.array(result.x)*np.array(conv)
print(bulks)

J = np.array(result.jac)
hess = np.matmul(np.transpose(J),J)
cov = np.linalg.inv(hess)
scale = np.diag(conv)
cov = np.matmul(np.matmul(scale,cov),scale)
var = np.diagonal(cov)
std = np.sqrt(var)
print("LSQ std: ", std)
invcov = np.linalg.inv(cov)
for row in invcov:
    print("  ".join(map(str, row)))

# Single Hartree Runs
#couplings = [110.349, 187.695, 192.927, 0.0, 3.26, -0.003551, 0.0235, 0.0, 0.043377, 0.0, 0.0, 0.0, 0.0, 4.0, 0.1, 496.939, 782.5, 763.0, 980.0] #FSUGarnet
couplings = [108.094300, 183.789300, 267.652249, 0.0, 3.0029, -0.000533, 0.0256, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 497.479, 782.5, 763.0, 980.0]
#couplings = bulks_to_params(result.x)
nuclei = 2
#observs = call_hartree(couplings,A[nuclei],Z[nuclei],0.0)