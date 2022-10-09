import numpy as np
import matplotlib.pyplot as plt
import math
from decimal import Decimal
import cmath
from qutip import *

inp=complex(input(""))
print(inp.real)
print(inp.imag)

def projector(q,w):
    proj=tensor(fock(q,w),fock(q,w).dag())
    proj.dims=[[q],[q]]
    return proj

def f_tol():
    ans=input("")
    if ans=="N":
        print("Tolerance rate set to default values 1e-08")
        return 1e-8
    if ans=="Y":
        print("Enter the tolerance value for the Cost Function [eg. 1.256e-06,2.56e-09...etc]:")
        tol=float(input(""))
        return tol
    if ans!="Y" or ans!="N":
        print("Enter a valid input 'Y' or 'N'")
        f_tol()

cf_tol=f_tol()
cf_tol_string="{:e}".format(cf_tol)
print("Tolerance rate set to %s"%(cf_tol_string))

psi_target=coherent(5,1)
print(psi_target.dims[0])
a=destroy(5)
psi_final=a*coherent(5,1)
psi_list=[]
for i in range(5-1):
    print(psi_target[i])
    psi_list.append(abs(psi_target[i][0][0]))

print(psi_list)
print(abs(psi_list[0]))
print(psi_target)
print(psi_final)


print(projector(2,0))
print(projector(2,0)+qeye(2))
print("----------------------------------")


print(math.cos((math.pi)/4))
#print(sigmax()+2)
#print(fock(5,2))
#print(fock(5,2).dag)
z=tensor(fock(2,0),fock(2,0).dag())+1
print(coherent(2,2))
#print(sigmax())
#print(sigmax()*coherent(2,2))
print(z)
#print(z.dims)
z.dims=[[2],[2]]
#print(z.dims)
print(z)
print("----------------")
print(z*coherent(2,2))
#print(complex(1,2)*sigmax())
print(cmath.exp(complex(0,1)*(0.25*(math.pi))))

def poisson_like(alpha,n):
    return (float(math.exp(-((alpha)**2)/2))*(alpha**n)*(1/float(math.sqrt(Decimal(math.factorial(n))))))

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
    plt.legend()
plt.show()
