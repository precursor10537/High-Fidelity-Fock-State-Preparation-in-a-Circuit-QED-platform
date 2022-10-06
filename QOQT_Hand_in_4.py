import numpy as np
import matplotlib.pyplot as plt
from qutip import *
import math
from decimal import Decimal
import cmath

def poisson_like(alpha,n):
    return (float(math.exp(-((alpha)**2)/2))*(alpha**n)*(1/float(math.sqrt(Decimal(math.factorial(n))))))

def truncate_hilbert(alpha_expected):
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

#NOTE: There is a bug in tensor product of qutip. Therefore, additional operations are required to define projectors properly for further calculations
def projector(q,w):
    proj=tensor(fock(q,w),fock(q,w).dag())
    proj.dims=[[q],[q]]
    return proj

########################################################################################################################################################################

#1. Define creation and annihilation operators for the cavity, the identity, the vacuum state, the projection operator onto a given Fock state. You will have to pick a
#truncation for the Hilbert space; what number is appropriate? How can you validate your choice?

#[Truncate] Dimensions of Hilbert Space
print('The expected upper bound of Coherent state mean amplitude "aplha" (or) number of Fock states "n" is used in the script to truncate the Hilbert Space as:') 
print('z is the upper bound Coherent state mean amplitude (or) Fock states \n -> N = 3*[(z)^2+((z)^2)/2], z < 2\n -> N = 2*[(z)^2+((z)^2)/2], 2 < z < 4\n -> N = 1.5*[(z)^2+((z)^2)/2], z > 4')
#n_expected=float(input(''))
#N=truncate_hilbert(n_expected)

#Basis of Truncation
print('Truncation of Hilbert Space was formulated to optimize the computational complexity based on the following: \n(i) Wigner function formulation for Fock, coherent, and squeezed states\n(ii) Trends seen in various Poisson like distributions of the Coherent state amplitudes as shown in following plot:')
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

print('\n')
print('--------------------------------------------------------------')
print('\n')

#Defining Operators as functions of "N":
print("Operators have been defined as function of the computed Hilbert Space dimension 'N'")
a=lambda q: destroy(q) #Annhilation Operator
a_dag=lambda q: create(q) #Creation Operator
I=lambda q: qeye(q) #Identity Operator
vac=lambda q: fock(q,0) #Vacuum state
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

print('The following is a [static] example to validate our truncation choice:')
print('(i) Wigner function of coherent state |7> in lower dimensional Hilbert space of dims=60 with truncation errors')
print('(ii) Wigner function of coherent state |7> in our calculated "N=111" truncated Hilbert space')
#Validate your truncation parameter's value
#Coherent state |7> in lower dimensional Hilbert space with truncation errors
plot_wigner(coherent(60,7),colorbar=True, projection='3d')
plt.title("Wigner function of coherent state |7> in Hilbert Space dims=60 in 3D")
plot_wigner(coherent(60,7),colorbar=True, projection='2d')
plt.title("Wigner function of coherent state |7> in Hilbert Space dims=60 in 2D")
#Coherent state |7> in "N=1.5*[(aplha)^2+((alpha)^2)/2]" dimensional Hilbert space with truncation errors
plot_wigner(coherent(111,7),colorbar=True, projection='3d')
plt.title("Wigner function of coherent state |7> in Hilbert Space dims N=111 calculated with our approximation in 3D")
plot_wigner(coherent(111,7),colorbar=True, projection='2d')
plt.title("Wigner function of coherent state |7> in Hilbert Space dims N=111 calculated with our approximation in 2D")
plt.show()   

print('\n')
print('--------------------------------------------------------------')
print('\n')

##########################################################################################################################################################################

#2 Define the SNAP and displacement operators. Test their action by applying them to suitable initial states and plotting the Wigner function before and after the gate.
D=lambda q,w: displace(q,w) #Displacement operator displaces the state with complex amplitude "w"
SNAP=lambda q,theta,w: (((cmath.exp(complex(0,1)*(theta*(math.pi)))-1))*projector(q,w))+I(q) #SNAP gate with Berry phase: theta, and conditioned on "w"th Fock state through the projector
#NOTE: M+1 = adding Identity Operator "I" to Matrix "M" in qutip

#Testing Displacement operator
print("Displacement Operator Test: \n-> psi_initial=|0> (vacuum state) \n-> D(alpha)*psi_initial=|aplha> (coherent state)")
print("Please enter 'alpha' to displace the vaccum state:")
alpha=float(input(''))
N=truncate_hilbert(alpha)
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
print("SNAP Operator Tests (1 - Coherent state and 2 - Fock state):")
print("Test 1:\n-> psi_initial=|alpha> (arbitrary coherent state) \n-> SNAP(|n>,theta)*psi_initial=|n'> (phase displaced fock state)")
print("Please enter 'alpha' to create an arbitrary initial coherent state:")
alpha=float(input(''))
N=truncate_hilbert(alpha)
psi_initial=coherent(N,alpha)
print("Please enter 'Fock state n' and \n'Berry Phase (theta factor to be multiplied by pi)' [eg. entering 2 --> defines 2*pi] \nto define the SNAP gate:")
n_SNAP=int(input('|n>: '))
theta_SNAP=float(input('theta: '))
#SNAP(N,theta_SNAP,n_SNAP)
print(psi_initial)
print(SNAP(N,theta_SNAP,n_SNAP))
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

print("Test 2:\n-> psi_initial=|n> (arbitrary fock state) \n-> SNAP(|n>,theta)*psi_initial=|n'> (phase displaced fock state)")
print("Please enter 'n' to create an arbitrary initial fock state:")
alpha=int(input(''))
N=truncate_hilbert(alpha)
psi_initial=fock(N,alpha)
print("Please enter 'Fock state n' and 'Berry Phase (theta factor to be multiplied by pi)' [eg. entering 2 --> defines 2*pi] to define the SNAP gate:")
n_SNAP=int(input('|n>: '))
theta_SNAP=float(input('theta: '))
#SNAP(N,theta_SNAP,n_SNAP)
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

##########################################################################################################################################################################
#3. Define final state for |1> with (alpha1, n, theta, alpha2). Write expression for Fidelity
print("We now define the Final state (with initial state to be vacuum state) as:")
print("-> psi_initial=|0> (vacuum state)")
print("-> psi_final=D(alpha2)*SNAP(n,theta)*D(alpha1)*|0>")
print("Please enter the Coherent state mean amplitude field 'alpha1' to define D(alpha1):")
alpha_1=float(input(""))
print("Please enter the Coherent state mean amplitude field 'alpha2' to define D(alpha2):")
alpha_2=float(input("")) 
if alpha_1>=alpha_2:
    print("Using 'alpha1' as a base to truncate the Hilbert Space...")
    N=truncate_hilbert(alpha_1)
if alpha_1<alpha_2:
    print("Using 'alpha1' as a base to truncate the Hilbert Space...")
    N=truncate_hilbert(alpha_2)
print("Please enter 'Fock state n' and 'Berry Phase (theta factor to be multiplied by pi)' [eg. entering 2 --> defines 2*pi] to define the SNAP gate:")
n_SNAP=int(input('|n>: '))
theta_SNAP=float(input('theta: '))

psi_initial=vac(N)
#D(N,alpha) SNAP(N,theta_SNAP,n_SNAP)
psi_final=D(N,alpha_2)*SNAP(N,theta_SNAP,n_SNAP)*D(N,alpha_1)*psi_initial

plot_wigner(psi_initial,colorbar=True, projection='3d')
plt.title("Wigner function of the initial state")
plot_wigner(psi_final,colorbar=True, projection='3d')
plt.title("Wigner function of the final state after the gate sequence in 3D")
plot_wigner(psi_final,colorbar=True, projection='2d')
plt.title("Wigner function of the final state after the gate sequence in 2D")
plt.show()

print('\n')
print('--------------------------------------------------------------')
print('\n')

#Fidelity: (tr(sqrt(sqrt(rho)*sigma*sqrt(rho))))^2
print("Please enter 'n' to create the Target Fock state:")
n_target=input("|n_target>: ")
psi_target=fock(N,n_target)
#Defining rho to be density matric of the pure target state

