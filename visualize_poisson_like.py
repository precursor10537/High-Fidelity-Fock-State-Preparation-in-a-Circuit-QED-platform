import numpy as np
import matplotlib.pyplot as plt
import math
from decimal import Decimal
import cmath
from qutip import *

def projector(q,w):
    proj=tensor(fock(q,w),fock(q,w).dag())
    proj.dims=[[q],[q]]
    return proj

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
