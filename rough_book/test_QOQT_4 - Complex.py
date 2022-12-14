import numpy as np
import matplotlib.pyplot as plt
from qutip import *
import math
from decimal import Decimal
import cmath
from scipy.optimize import least_squares
from scipy.optimize import minimize

def poisson_like(alpha,n):
    return (float(math.exp(-((alpha)**2)/2))*(alpha**n)*(1/float(math.sqrt(Decimal(math.factorial(n))))))

def truncate_hilbert(alpha_expected):
    alpha_expected=abs(alpha_expected)
    if alpha_expected<=2:
        N=math.ceil(3*(((int(alpha_expected))**2)+math.ceil(((int(alpha_expected))**2)/2)))
    if alpha_expected<=4 and alpha_expected>2:
        N=math.ceil(2*(((int(alpha_expected))**2)+math.ceil(((int(alpha_expected))**2)/2)))
    if alpha_expected>4:
        N=math.ceil(1.5*(((int(alpha_expected))**2)+math.ceil(((int(alpha_expected))**2)/2)))
    if alpha_expected==0:
        N=50
        print("Entered value is 0, truncating Hilbert Space to default value")
    print('Truncating Hilbert Space to N=%d dimensions'%(N))
    return N

def N_error_val(N,alpha):
    psi_target=coherent(N+1,alpha)
    #print(psi_target)
    psi_final=destroy(N+1)*psi_target
    #print(psi_final)
    psi_final=(1/alpha)*psi_final
    target_list=[]
    final_list=[]
    error=0
    for i in range(N):
        target_list.append(abs(psi_target[i][0][0]))
        final_list.append(abs(psi_final[i][0][0]))
    #print(target_list)
    #print(final_list)
    for i in range(N):
        error=error+(target_list[i]-final_list[i])**2
    #print(error)
    return error

def SSE_minimizer(N,alpha,error_bound):
    count=0
    while N_error_val(N,alpha)>=error_bound:
        N=N+2
        #print(N)
        count=count+1
    error_calculated=N_error_val(N,alpha)
    error_calculated="{:e}".format(error_calculated)
    error_bound_string="{:e}".format(error_bound)
    if count==0:
        print("Calculated error %s is less than given error bound %s \nHence no changes made to the truncated Hilbert Space"%(error_calculated,error_bound_string))
        print("Truncated Hilbert Space has %d dimensions"%(N))
    if count!=0:
        print("Calculated error was greater than given error bound %s \nAdding dimensions to truncated Hilbert Space...."%(error_bound_string))
        print("Added %d dimensions and the calculated error is: %s which is now lesser than the given error bound %s"%(2*count,error_calculated,error_bound_string))
        print("Truncated Hilbert Space now has %d dimensions"%(N))
    return N

def a(q): #a=lambda q: destroy(q) #Annhilation Operator
    return destroy(q)

def a_dag(q): #a_dag=lambda q: create(q) #Creation Operator
    return create(q)

def I(q): #I=lambda q: qeye(q) #Identity Operator
    return qeye(q)

def vac(q): #vac=lambda q: fock(q,0) #Vacuum state
    return fock(q,0)

#NOTE: There is a bug in tensor product of qutip. Therefore, additional operations are required to define projectors properly for further calculations
def projector(q,w):
    proj=tensor(fock(q,w),fock(q,w).dag())
    proj.dims=[[int(q)],[int(q)]]
    return proj

def D(q,w): #D=lambda q,w: displace(q,w) #Displacement operator displaces the state with complex amplitude "w"
    return displace(q,w)

def SNAP(q,theta,w): #SNAP=lambda q,theta,w: (((cmath.exp(complex(0,1)*(theta*(math.pi)))-1))*projector(q,w))+I(q) #SNAP gate with Berry phase: theta, and conditioned on "w"th Fock state through the projector
    return (((cmath.exp(complex(0,1)*(theta))-1))*projector(q,w))+I(q) #SNAP gate with Berry phase: theta, and conditioned on "w"th Fock state through the projector

def density_matrix(state):
    density_operator=tensor(state,state.dag())
    density_operator.dims=[state.dims[0],state.dims[0]]
    return density_operator

def f_tol():
    ans=input("")
    if ans=="N" or ans=='n' or ans=='no' or ans=='No' or ans=='NO':
        print("f_tol tolerance rate set to default values 1e-08")
        return 1e-8
    elif ans=="Y"or ans=='y' or ans=='yes' or ans=='Yes' or ans=='YES':
        print("Enter the tolerance value for the Cost Function [f_tol] [eg. 1.256e-06,2.56e-09...etc]:")
        tol=float(input(""))
        return tol
    else:
        print("Enter a valid input 'Y' or 'N'")
        return f_tol()

def x_tol():
    ans=input("")
    if ans=="N" or ans=='n' or ans=='no' or ans=='No' or ans=='NO':
        print("x_tol tolerance rate set to default values 1e-08")
        return 1e-8
    elif ans=="Y"or ans=='y' or ans=='yes' or ans=='Yes' or ans=='YES':
        print("Enter the tolerance value for the change in independent variables [x_tol] [eg. 1.256e-06,2.56e-09...etc]:")
        tol=float(input(""))
        return tol
    else:
        print("Enter a valid input 'Y' or 'N'")
        return x_tol()


def fidelity_cf(N,psi_target,x): #1-(tr(sqrt(sqrt(rho)*sigma*sqrt(rho))))^2
    #infidelity_cf([alpha_1.real,alpha_1.imag,theta_SNAP,n_SNAP,alpha_2.real,alpha_2.imag],N,psi_target)
    if len(x)==4:
        alpha_1=float(x[0])
        theta_SNAP=float(x[1])
        n_SNAP=int(x[2])
        alpha_2=float(x[3])
    else:
        alpha_1=complex(x[0],x[1])
        theta_SNAP=float(x[2])
        n_SNAP=int(x[3])
        alpha_2=complex(x[4],x[5])
    psi_initial=vac(N)
    #D(N,alpha) SNAP(N,theta_SNAP,n_SNAP)
    psi_final=D(N,alpha_2)*SNAP(N,theta_SNAP,n_SNAP)*D(N,alpha_1)*psi_initial
    rho=density_matrix(psi_target) #Defining 'rho' to be density matric of the pure target state
    sigma=density_matrix(psi_final) #Defining 'sigma' to be density matric of the psi final state
    sqrt_rho=rho.sqrtm()
    rho_sigma=sqrt_rho*sigma*sqrt_rho
    sqrt_rho_sigma=rho_sigma.sqrtm()
    fidelity=(sqrt_rho_sigma.tr())**2
    #print("Fidelity:%f"%(fidelity))
    return fidelity

def guess_enumeration(n_target,alpha_1,alpha_2,n_SNAP):
    if n_target-1>0: 
        n_SNAP_enumeration=[n_SNAP,n_target-1,n_target,n_target+1]
    if n_target-1<=0:
        n_SNAP_enumeration=[n_SNAP,0,n_target,n_target+1]
    n_SNAP_enumeration=list(set(n_SNAP_enumeration)) #To remove duplicates
    alpha_enumeration=[[alpha_1,alpha_2],[-1*(alpha_1),alpha_2],[alpha_1,-1*(alpha_2)]]
    return n_SNAP_enumeration,alpha_enumeration

def plot_optimal(psi_initial,res):
    if len(res.x)==4:
        alpha_1=float(res.x[0])
        theta_SNAP=float(res.x[1])
        n_SNAP=int(res.x[2])
        alpha_2=float(res.x[3])
        #alpha_1=complex(res.x[0],res.x[1])
        #theta_SNAP=res.x[2]
        #n_SNAP=int(res.x[3])
        #alpha_2=complex(res.x[4],res.x[5])
        optimized_psi_final=D(N,alpha_2)*SNAP(N,theta_SNAP,n_SNAP)*D(N,alpha_1)*psi_initial
        #Plot Ideal Fock State
        plot_wigner(psi_target,colorbar=True, projection='3d')
        plt.title("Wigner function of the ideal Fock State |%d> in 3D"%(n_target))
        plot_wigner(psi_target,colorbar=True, projection='2d')
        plt.title("Wigner function of the ideal Fock State |%d> in 2D"%(n_target))
        #Plot Optimized Final State
        plot_wigner(optimized_psi_final,colorbar=True, projection='3d')
        plt.title("Wigner function of the OPTIMIZED final state with fidelity= %f for alpha 1=%f+%fj, SNAP n=%d, SNAP theta=%f, alpha 2=%f+%fj in 3D"%(1-(res['fun'][0]),alpha_1,0,n_SNAP,theta_SNAP,alpha_2,0))
        plot_wigner(optimized_psi_final,colorbar=True, projection='2d')
        plt.title("Wigner function of the OPTIMIZED final state with fidelity= %f for alpha 1=%f+%fj, SNAP n=%d, SNAP theta=%f, alpha 2=%f+%fj in 2D"%(1-(res['fun'][0]),alpha_1,0,n_SNAP,theta_SNAP,alpha_2,0))
        plt.show()
    else:
        alpha_1=float(res.x[0])
        theta_SNAP=res.x[1]
        n_SNAP=int(res.x[2])
        alpha_2=float(res.x[3])
        #alpha_1=complex(res.x[0],res.x[1])
        #theta_SNAP=res.x[2]
        #n_SNAP=int(res.x[3])
        #alpha_2=complex(res.x[4],res.x[5])
        optimized_psi_final=D(N,alpha_2)*SNAP(N,theta_SNAP,n_SNAP)*D(N,alpha_1)*psi_initial
        #Plot Ideal Fock State
        plot_wigner(psi_target,colorbar=True, projection='3d')
        plt.title("Wigner function of the ideal Fock State |%d> in 3D"%(n_target))
        plot_wigner(psi_target,colorbar=True, projection='2d')
        plt.title("Wigner function of the ideal Fock State |%d> in 2D"%(n_target))
        #Plot Optimized Final State
        plot_wigner(optimized_psi_final,colorbar=True, projection='3d')
        plt.title("Wigner function of the OPTIMIZED final state with fidelity= %f for alpha 1=%f+%fj, SNAP n=%d, SNAP theta=%f, alpha 2=%f+%fj in 3D"%(1-(res['fun'][0]),alpha_1.real,alpha_1.imag,n_SNAP,theta_SNAP,alpha_2.real,alpha_2.imag))
        plot_wigner(optimized_psi_final,colorbar=True, projection='2d')
        plt.title("Wigner function of the OPTIMIZED final state with fidelity= %f for alpha 1=%f+%fj, SNAP n=%d, SNAP theta=%f, alpha 2=%f+%fj in 2D"%(1-(res['fun'][0]),alpha_1.real,alpha_1.imag,n_SNAP,theta_SNAP,alpha_2.real,alpha_2.imag))
        plt.show()

#3&4. Define final state for |1> with (alpha1, n, theta, alpha2). Write expression for Fidelity
print("We now define the Final state (with initial state to be vacuum state) as:")
print("-> psi_initial=|0> (vacuum state)")
print("-> psi_final=D(alpha2)*SNAP(n,theta)*D(alpha1)*|0>")
print("Please enter Sum of Squared Error Bound on the trucation parameter [eg. 1.256e-06,2.56e-09...etc]:")
error_bound=float(input(""))
print("Please enter the Coherent state mean amplitude field 'alpha1' to define D(alpha1):")
alpha_1=complex(input(""))
print("Please enter the Coherent state mean amplitude field 'alpha2' to define D(alpha2):")
alpha_2=complex(input(""))
alpha_1_original=alpha_1
alpha_2_original=alpha_2 
if abs(alpha_1)>=abs(alpha_2):
    print("Using 'alpha1' as a base to truncate the Hilbert Space...")
    N=truncate_hilbert(alpha_1)
    N=SSE_minimizer(N,alpha_1,error_bound)
if abs(alpha_1)<abs(alpha_2):
    print("Using 'alpha2' as a base to truncate the Hilbert Space...")
    N=truncate_hilbert(alpha_2)
    N=SSE_minimizer(N,alpha_2,error_bound)
print("Please enter 'Fock state n' and 'Berry Phase (theta factor to be multiplied by pi)' [eg. entering 2 --> defines 2*pi] to define the SNAP gate:")
n_SNAP=int(input('|n>: '))
theta_SNAP=float(input('theta: '))
theta_SNAP=theta_SNAP*(math.pi)
#Complementing the alpha1 and aplha2 for the gate sequence to work
#Plotting UNOPTIMIZED initial gate sequence
psi_initial=vac(N)
#D(N,alpha) SNAP(N,theta_SNAP,n_SNAP)
psi_final=D(N,alpha_2)*SNAP(N,theta_SNAP,n_SNAP)*D(N,alpha_1)*psi_initial
#plot_wigner(psi_initial,colorbar=True, projection='3d')
#plt.title("Wigner function of the initial state")
#plot_wigner(psi_final,colorbar=True, projection='3d')
#plt.title("Wigner function of the final state after the gate sequence with the input parameters in 3D")
#plot_wigner(psi_final,colorbar=True, projection='2d')
#plt.title("Wigner function of the final state after the gate sequence with the input parameters in 2D")
#plt.show()

print('\n')
print('--------------------------------------------------------------')
print('\n')

#Defining Infidelity as Cost function of (alpha1,(n,theta),alpha2) and its minimization
print("Please enter 'n' to create the Target Fock state to minimize the Cost function and to approximate the optimal parameters:")
n_target=int(input("|n_target>: "))
psi_target=fock(N,n_target) 
print("Do you want to set a custom tolerance rate for the cost function [f_tol] (Y/N)? [DEFAULT: 1e-08]")
cf_tol=f_tol()
print(cf_tol)
cf_tol_string="{:e}".format(cf_tol)
print("Tolerance rate set to %s"%(cf_tol_string))
print("Do you want to set a custom tolerance rate for the change in independent variables [x_tol] (Y/N)? [DEFAULT: 1e-08]")
cx_tol=x_tol()
print(cx_tol)
cx_tol_string="{:e}".format(cx_tol)
print("Tolerance rate set to %s"%(cx_tol_string))

#Defining infidelity_cf
def infidelity_cf(x):
    for i in x:
        #print(i)
        if math.isnan(i)==True: #To avoid spurious differentiation errors
            return 1
    return 1-fidelity_cf(N,psi_target,x)

implicit_fidelity=0.95
implicit_infidelity=0.05

#Minimizing Cost Function with Least Squares: guess_enumeration(n_target,alpha_1,alpha_2,n_SNAP)
n_SNAP_enumeration,alpha_enumeration=guess_enumeration(n_target,alpha_1,alpha_2,n_SNAP)
#n_SNAP_enumeration=[n_SNAP,n_target-1,n_target,n_target+1]
#alpha_enumeration=[[alpha_1,alpha_2],[-1*(alpha_1),-1*(alpha_2)],[-1*(alpha_1),alpha_2],[alpha_1,-1*(alpha_2)]]
result_dict={} #For archiving all the outputs
result_list=[] #Only stores qualified outputs
qualified_list=[] #Failed optimization
min_fun=2 #should keep it greater than 2
print("Optimizing your input...[function evaluation limit [max_nfev]:50]\nPlease do not poweroff your PC...")
#Conditioned on whether we have complex alphas
if alpha_1_original.imag==0 and alpha_2_original.imag==0:
    for i in n_SNAP_enumeration:
        j_counter=0
        for j in alpha_enumeration:
            x0=np.array([j[0].real,float(theta_SNAP),int(i),j[1].real])
            j_counter=j_counter+1
            result=least_squares(infidelity_cf,x0,ftol=cf_tol,xtol=cx_tol,gtol=None,max_nfev=25)
            result_dict[("%d,%d"%(i,j_counter))]=[x0,result] #Archive
            if result['fun'][0]<implicit_infidelity:
                result_list.append([x0,result]) #Store only relevant outputs
            if result['fun'][0]==min_fun:
                qualified_list.append([x0,result]) #Store equal failed outputs
            if result['fun'][0]<min_fun:
                qualified_list=[0] #Discard existing elements
                qualified_list[0]=[x0,result] #Store best failed output
                min_fun=result['fun'][0]
else:
    for i in n_SNAP_enumeration:
        j_counter=0
        for j in alpha_enumeration:
            x0=np.array([j[0].real,j[0].imag,float(theta_SNAP),int(i),j[1].real,j[1].imag])
            j_counter=j_counter+1
            result=least_squares(infidelity_cf,x0,ftol=cf_tol,xtol=cx_tol,gtol=None,max_nfev=25)
            result_dict[("%d,%d"%(i,j_counter))]=[x0,result] #Archive
            if result['fun'][0]<implicit_infidelity:
                result_list.append([x0,result]) #Store only relevant outputs
            if result['fun'][0]==min_fun:
                qualified_list.append([x0,result]) #Store equal failed outputs
            if result['fun'][0]<min_fun:
                qualified_list=[0] #Discard existing elements
                qualified_list[0]=[x0,result] #Store best failed output
                min_fun=result['fun'][0]

#result=minimize(infidelity_cf,x0,method='BFGS',tol=cf_tol)
#result=least_squares(infidelity_cf,x0,ftol=cf_tol,xtol=cx_tol,gtol=None,max_nfev=1000,verbose=2)
#print(result_dict)
x_initial=np.array([alpha_1_original,theta_SNAP,n_SNAP,alpha_2_original])
print("Initial guess parameters are:")
print(x_initial)
if len(result_list)!=0:
    print("The Optimized Outputs reaching the implicit state fidelity threshold of 0.95 are:")
    for i in range(len(result_list)):
        print(result_list[i][1])
        print("\n")
        plot_optimal(psi_initial,result_list[i][1]) #plot_optimal(psi_initial,res)
else:
    print("Proper Optimization NOT ACHIEVED. Output states fidelities less than the implicit state fidelity threshold of 0.95...")
    print("Please revise your initial guess of the gate sequence parameters or change the tolerances [or change the number of function evaluations max_nfev in the code]")
    print("Displaying the outputs for the highest obtained fidelities:")
    for i in range(len(qualified_list)):
        print(qualified_list[i][1])
        print("\n")
        plot_optimal(psi_initial,qualified_list[i][1]) #plot_optimal(psi_initial,res)