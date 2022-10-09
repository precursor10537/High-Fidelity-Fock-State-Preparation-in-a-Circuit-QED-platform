import numpy as np
import matplotlib.pyplot as plt
from qutip import *
import math
from decimal import Decimal
import cmath
from scipy.optimize import least_squares

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
    proj.dims=[[q],[q]]
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
    if ans=="N":
        #print("Tolerance rate set to default values 1e-08")
        return 1e-8
    if ans=="Y":
        print("Enter the tolerance value for the Cost Function [eg. 1.256e-06,2.56e-09...etc]:")
        tol=float(input(""))
        return tol
    if ans!="Y" or ans!="N":
        print("Enter a valid input 'Y' or 'N'")
        return f_tol()

def fidelity_cf(N,psi_target,x): #1-(tr(sqrt(sqrt(rho)*sigma*sqrt(rho))))^2
    #infidelity_cf([alpha_1,theta_SNAP,n_SNAP,alpha_2],N,psi_target)
    alpha_1=x[0]
    theta_SNAP=x[1]
    n_SNAP=int(x[2])
    alpha_2=x[3]
    psi_initial=vac(N)
    #D(N,alpha) SNAP(N,theta_SNAP,n_SNAP)
    psi_final=D(N,alpha_2)*SNAP(N,theta_SNAP,n_SNAP)*D(N,alpha_1)*psi_initial
    rho=density_matrix(psi_target) #Defining 'rho' to be density matric of the pure target state
    sigma=density_matrix(psi_final) #Defining 'sigma' to be density matric of the psi final state
    sqrt_rho=rho.sqrtm()
    rho_sigma=sqrt_rho*sigma*sqrt_rho
    sqrt_rho_sigma=rho_sigma.sqrtm()
    fidelity=(sqrt_rho_sigma.tr())**2
    return fidelity

########################################################################################################################################################################

#1. Define creation and annihilation operators for the cavity, the identity, the vacuum state, the projection operator onto a given Fock state. You will have to pick a
#truncation for the Hilbert space; what number is appropriate? How can you validate your choice?

#[Truncate] Dimensions of Hilbert Space
print('The expected upper bound of Coherent state mean amplitude "aplha" (or) number of Fock states "n" is used in the script for an INITIAL Truncation of the Hilbert Space as:') 
print('z is the upper bound Fock states (or) Coherent state mean amplitude [IMPLICIT ASSUMPTION: We would be using only real values of "aplha"]\n -> N = 3*[(z)^2+((z)^2)/2], z < 2\n -> N = 2*[(z)^2+((z)^2)/2], 2 < z < 4\n -> N = 1.5*[(z)^2+((z)^2)/2], z > 4')
#n_expected=float(input(''))
#N=truncate_hilbert(n_expected)

#Basis of Truncation
print("\n")
print('INITIAL Truncation of Hilbert Space was formulated to optimize the computational complexity based on the following: \n(i) Wigner function formulation for Fock, coherent, and squeezed states\n(ii) Trends seen in various Poisson like distributions of the Coherent state amplitudes as shown in following plot:')
alpha=22
n=200
coherent_amplitudes=[]
counter=-1
mean=[]
for i in range(0,alpha,2):
    coherent_amplitudes.append([])
    counter=counter+1
    mean.append(i)
    for j in range(n):
        coherent_amplitudes[counter].append(poisson_like(i,j))

x=np.linspace(0,200,200)
#print(x)
#print(coherent_amplitudes)
fig,ax=plt.subplots()
for i in range(len(coherent_amplitudes)):
    ax.plot(x,coherent_amplitudes[i],label="alpha=%d"%(mean[i]))
    ax.set_xlabel('Number of Fock states |n>')
    ax.set_ylabel('Poisson Distribution like Coherent state amplitudes')
    ax.set_title('Coherent State amplitude distribution trends')
    plt.legend()
plt.show()

print("We further refine our INITIAL Truncation of Hilbert Space in the computations with the custom input Sum of Squared error bounds throughout the script....")

print('\n')
print('--------------------------------------------------------------')
print('\n')

#Defining Operators as functions of "N":
print("Operators have been defined as function of the computed Hilbert Space dimension 'N'")
#a=lambda q: destroy(q) #Annhilation Operator
#a_dag=lambda q: create(q) #Creation Operator
#I=lambda q: qeye(q) #Identity Operator
#vac=lambda q: fock(q,0) #Vacuum state
#projector=lambda q,w: tensor(fock(q,w),fock(q,w).dag()) #Projection Operator for given a given "w"th Fock state 
#NOTE: There are bugs in outer/tensor product and hence cannot write projector in lambda function

#Print the defined operators
#print(a(N))
#print(a_dag(N))
#print(I(N))
#print(vac(N))
#for i in range(N):
#    print(projector(N,i))    

print('\n')
print('--------------------------------------------------------------')
print('\n')


#Validate your truncation parameter's value
print("Please enter Sum of Squared Error Bound on the trucation parameter [eg. 1.256e-06,2.56e-09...etc]:")
error_bound=float(input(""))
print("Please enter Coherent state mean amplitude 'alpha' to validate our truncation choice:")
alpha=float(input(""))
N=truncate_hilbert(alpha)
N=SSE_minimizer(N,alpha,error_bound)

print('\n')
print('--------------------------------------------------------------')
print('\n')

print('The following is a [static] example to validate our truncation choice [when "alpha" has imaginary part]: ')
print('(i) Wigner function of coherent state |7> in lower dimensional Hilbert space of dims=60 with truncation errors')
print('(ii) Wigner function of coherent state |7> in our calculated "N=111" truncated Hilbert space')
#Coherent state |7> in lower dimensional Hilbert space with truncation errors
plot_wigner(coherent(60,7),alpha_max=15,colorbar=True, projection='3d')
plt.title("Wigner function of coherent state |7> in Hilbert Space dims=60 in 3D")
plot_wigner(coherent(60,7),alpha_max=15,colorbar=True, projection='2d')
plt.title("Wigner function of coherent state |7> in Hilbert Space dims=60 in 2D")
#Coherent state |7> in "N=1.5*[(aplha)^2+((alpha)^2)/2]" dimensional Hilbert space with truncation errors
plot_wigner(coherent(111,7),alpha_max=15,colorbar=True, projection='3d')
plt.title("Wigner function of coherent state |7> in Hilbert Space dims N=111 calculated with our approximation in 3D")
plot_wigner(coherent(111,7),alpha_max=15,colorbar=True, projection='2d')
plt.title("Wigner function of coherent state |7> in Hilbert Space dims N=111 calculated with our approximation in 2D")
plt.show()   

print('\n')
print('--------------------------------------------------------------')
print('\n')

##########################################################################################################################################################################

#2 Define the SNAP and displacement operators. Test their action by applying them to suitable initial states and plotting the Wigner function before and after the gate.
#D=lambda q,w: displace(q,w) #Displacement operator displaces the state with complex amplitude "w"
#SNAP=lambda q,theta,w: (((cmath.exp(complex(0,1)*(theta*(math.pi)))-1))*projector(q,w))+I(q) #SNAP gate with Berry phase: theta, and conditioned on "w"th Fock state through the projector
#NOTE: M+1 = adding Identity Operator "I" to Matrix "M" in qutip

#Testing Displacement operator
print("Displacement Operator Test: \n-> psi_initial=|0> (vacuum state) \n-> D(alpha)*psi_initial=|aplha> (coherent state)")
print("Please enter 'alpha' to displace the vaccum state:")
alpha=float(input(''))
print("Please enter Sum of Squared Error Bound on the trucation parameter [eg. 1.256e-06,2.56e-09...etc]:")
error_bound=float(input(""))
N=truncate_hilbert(alpha)
N=SSE_minimizer(N,alpha,error_bound)
#D(N,alpha)
psi_initial=vac(N)
psi_final=D(N,alpha)*psi_initial
plot_wigner(psi_initial,colorbar=True, projection='3d')
plt.title("Wigner function of the initial state")
plot_wigner(psi_final,colorbar=True, projection='3d')
plt.title("Wigner function of the final state after applying Displacement Operator in 3D")
plot_wigner(psi_final,colorbar=True, projection='2d')
plt.title("Wigner function of the final state after applying Displacement Operator in 2D")
plt.show()

print('\n')
print('--------------------------------------------------------------')
print('\n')

#Testing SNAP operator
print("SNAP Operator Test on arbitrary Coherent state:")
print("Test 1:\n-> psi_initial=|alpha> (arbitrary coherent state) \n-> SNAP(|n>,theta)*psi_initial=|n'> (phase displaced fock state)")
print("Please enter 'alpha' to create an arbitrary initial coherent state:")
alpha=float(input(''))
print("Please enter Sum of Squared Error Bound on the trucation parameter [eg. 1.256e-06,2.56e-09...etc]:")
error_bound=float(input(""))
N=truncate_hilbert(alpha)
N=SSE_minimizer(N,alpha,error_bound)
psi_initial=coherent(N,alpha)
print("Please enter 'Fock state n' and \n'Berry Phase (theta factor to be multiplied by pi)' [eg. entering 2 --> defines 2*pi] \nto define the SNAP gate:")
n_SNAP=int(input('|n>: '))
theta_SNAP=float(input('theta: '))
#SNAP(N,theta_SNAP,n_SNAP)
#print(psi_initial)
#print(SNAP(N,theta_SNAP,n_SNAP))
psi_final=SNAP(N,theta_SNAP,n_SNAP)*psi_initial
plot_wigner(psi_initial,colorbar=True, projection='3d')
plt.title("Wigner function of the initial state")
plot_wigner(psi_final,colorbar=True, projection='3d')
plt.title("Wigner function of the final state after applying SNAP Operator in 3D")
plot_wigner(psi_final,colorbar=True, projection='2d')
plt.title("Wigner function of the final state after applying SNAP Operator in 2D")
plt.show()

print('\n')
print('--------------------------------------------------------------')
print('\n')

#print("Test 2:\n-> psi_initial=|n> (arbitrary fock state) \n-> SNAP(|n>,theta)*psi_initial=|n'> (phase displaced fock state)")
#print("Please enter 'n' to create an arbitrary initial fock state:")
#alpha=int(input(''))
#N=truncate_hilbert(alpha)
#psi_initial=fock(N,alpha)
#print("Please enter 'Fock state n' and 'Berry Phase (theta factor to be multiplied by pi)' [eg. entering 2 --> defines 2*pi] to define the SNAP gate:")
#n_SNAP=int(input('|n>: '))
#theta_SNAP=float(input('theta: '))
#SNAP(N,theta_SNAP,n_SNAP)
#psi_final=SNAP(N,theta_SNAP,n_SNAP)*psi_initial
#plot_wigner(psi_initial,colorbar=True, projection='3d')
#plt.title("Wigner function of the initial state")
#plot_wigner(psi_final,colorbar=True, projection='3d')
#plt.title("Wigner function of the final state after applying SNAP Operator in 3D")
#plot_wigner(psi_final,colorbar=True, projection='2d')
#plt.title("Wigner function of the final state after applying SNAP Operator in 2D")
#plt.show()

#print('\n')
#print('--------------------------------------------------------------')
#print('\n')

##########################################################################################################################################################################
#3&4. Define final state for |1> with (alpha1, n, theta, alpha2). Write expression for Fidelity
print("We now try to create a Final (custom) Fock state with our gate sequence and initial state = vacuum state as:")
print("-> psi_initial=|0> (vacuum state)")
print("-> psi_final=D(alpha2)*SNAP(n,theta)*D(alpha1)*|0>")
print("Please enter Sum of Squared Error Bound on the trucation parameter [eg. 1.256e-06,2.56e-09...etc]:")
error_bound=float(input(""))
print("Please enter the Coherent state mean amplitude field 'alpha1' to define D(alpha1):")
alpha_1=float(input(""))
print("Please enter the Coherent state mean amplitude field 'alpha2' to define D(alpha2):")
alpha_2=float(input("")) 
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
plot_wigner(psi_initial,colorbar=True, projection='3d')
plt.title("Wigner function of the initial state")
plot_wigner(psi_final,colorbar=True, projection='3d')
plt.title("Wigner function of the UNOPTIMIZED final state after the gate sequence with the UNOPTIMIZED input parameters in 3D")
plot_wigner(psi_final,colorbar=True, projection='2d')
plt.title("Wigner function of the UNOPTIMIZED final state after the gate sequence with the UNOPTIMIZED input parameters in 2D")
plt.show()

print('\n')
print('--------------------------------------------------------------')
print('\n')

#Defining Infidelity as Cost function of (alpha1,(n,theta),alpha2) and its minimization
print("Please enter 'n' to create the Target Fock state to minimize the Cost function and to approximate the optimal parameters:")
n_target=int(input("|n_target>: "))
psi_target=fock(N,n_target) 
print("Do you want to set a custom tolerance rate for the cost function (Y/N)? [DEFAULT: 1e-08]")
cf_tol=f_tol()
cf_tol_string="{:e}".format(cf_tol)
print("Tolerance rate set to %s"%(cf_tol_string))

#Defining infidelity_cf
def infidelity_cf(x):
    return 1-fidelity_cf(N,psi_target,x)

#Minimizing Cost Function with Least Squares
if alpha_1<0 and alpha_2<0:
    alpha_2=(-1)*alpha_2
if alpha_1>0 and alpha_2>0:
    alpha_2=(-1)*alpha_2
x0=np.array([alpha_1,theta_SNAP,n_SNAP,alpha_2])
result=least_squares(infidelity_cf,x0,ftol=cf_tol)
print(result)

#Plotting the Optimal Parameter Wigner Functions
print(type(result))
alpha_1=result.x[0]
theta_SNAP=result.x[1]
n_SNAP=int(result.x[2])
alpha_2=result.x[3]
optimized_psi_final=D(N,alpha_2)*SNAP(N,theta_SNAP,n_SNAP)*D(N,alpha_1)*psi_initial
plot_wigner(optimized_psi_final,colorbar=True, projection='3d')
plt.title("Wigner function of the OPTIMIZED final state after the gate sequence with the OPTIMAL input parameters in 3D")
plot_wigner(optimized_psi_final,colorbar=True, projection='2d')
plt.title("Wigner function of the OPTIMIZED final state after the gate sequence with the OPTIMAL input parameters in 2D")
plot_wigner(psi_target,colorbar=True, projection='3d')
plt.title("Wigner function of the ideal Fock State |%d> in 3D"%(n_target))
plot_wigner(psi_target,colorbar=True, projection='2d')
plt.title("Wigner function of the ideal Fock State |%d> in 2D"%(n_target))
plt.show()