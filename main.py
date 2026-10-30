import numpy as np
import matplotlib.pyplot as plt

def S(t):
    return sum(
        j*np.cos((j+1)*t+j)
        for j in range(1, 6)
    )

def dS(t):
    return sum(
        -j*(j+1)*np.sin((j+1)*t+j)
        for j in range(1, 6)
    )

def Shubert(x):
    return S(x[0])*S(x[1])

def dShubert(x):
    return np.array([
        dS(x[0])*S(x[1]),
        S(x[0])*dS(x[1])
    ])
x0 = np.array([-1, 2])
F=Shubert
dF=dShubert
sol, iterx1, itery1,n, Nf, Ng=GrippoMethod(dF, F, x0, Niter=10000, tol=1e-4, M=3)
print("The solution is:", sol)
print("The number of iterations is:", n)
x1, iterx2, itery2,n1, Nf, Ng=Hager_Zhang(x0,F,dF,tol=1e-4, Niter=10000,nuk=0.85,Q0=1)
print("The solution is:", x1)
print("The number of iterations is:", n1)
x2, n2, iterx3, itery3, Nf, Ng=Modified_Hager_Zhang(x0,F,dF,alpha=1,tol=1e-10, Niter=10000,nuk=0.85,Q0=1)
print("The solution is:", x2)
print("The number of iterations is:", n2)
plt.figure()
m=35
plt.plot(iterx1[:m], itery1[:m], '--',c='m',label='Algorithm 1-max')
plt.plot(iterx2[:m], itery2[:m], '--',c='b',label='Algorithm 2-Convex')
plt.plot(iterx3[:m], itery3[:m], '--',c='k',label='Algorithm 3-Mod')
plt.xlabel('Iterations: k')
plt.ylabel('F(xk)')
plt.title('F(xk) vs Iterations k')
plt.legend()
plt.grid(True)
plt.show()
x0 = np.array([0.2, -0.5])
F=Shubert
dF=dShubert
sol, iterx11, itery11,n, Nf, Ng=GrippoMethod(dF, F, x0, Niter=10000, tol=1e-4, M=3)
print("The solution is:", sol)
print("The number of iterations is:", n)
x1, iterx21, itery21,n1, Nf, Ng=Hager_Zhang(x0,F,dF,tol=1e-4, Niter=10000,nuk=0.85,Q0=1)
print("The solution is:", x1)
print("The number of iterations is:", n1)
x2, n2, iterx31, itery31, Nf, Ng=Modified_Hager_Zhang(x0,F,dF,alpha=1,tol=1e-10, Niter=10000,nuk=0.85,Q0=1)
print("The solution is:", x2)
print("The number of iterations is:", n2)
plt.figure()
m=35
m = 40

fig, axes = plt.subplots(1, 2, figsize=(12, 5))

# First initial point
axes[0].plot(iterx1[:m], itery1[:m], 'o--', c='m', label='Algorithm 1')
axes[0].plot(iterx2[:m], itery2[:m], 'o--', c='b', label='Algorithm 2')
axes[0].plot(iterx3[:m], itery3[:m], 'o--', c='k', label='Algorithm 3')

axes[0].set_title('Initial point 1: [0.2, -0.5]')
axes[0].set_xlabel('$k$')
axes[0].set_ylabel('$F(x_k)$')
axes[0].legend()
axes[0].grid(True)

# Second initial point
axes[1].plot(iterx11[:m], itery11[:m], 'o--', c='m', label='Algorithm 1')
axes[1].plot(iterx21[:m], itery21[:m], 'o--', c='b', label='Algorithm 2')
axes[1].plot(iterx31[:m], itery31[:m], 'o--', c='k', label='Algorithm 3')

axes[1].set_title('Initial point 2: x0=[-1, 2]')
axes[1].set_xlabel('$k$')
axes[1].set_ylabel('$F(x_k)$')
axes[1].legend()
axes[1].grid(True)

plt.tight_layout()
plt.show()
import numpy as np

def Rastrigin(x):
    return 20 + sum(
        x[i]**2 - 10*np.cos(2*np.pi*x[i])
        for i in range(2)
    )
def dRastrigin(x):
    return np.array([
        2*x[0] + 20*np.pi*np.sin(2*np.pi*x[0]),
        2*x[1] + 20*np.pi*np.sin(2*np.pi*x[1])
    ])
x0 = np.array([-0.5,0.9])
F=Rastrigin
dF=dRastrigin
sol, iterx1, itery1,n, Nf, Ng=GrippoMethod(dF, F, x0, Niter=10000, tol=1e-4, M=3)
print("The solution is:", sol)
print("The number of iterations is:", n)
x1, iterx2, itery2,n1, Nf, Ng=Hager_Zhang(x0,F,dF,tol=1e-4, Niter=10000,nuk=0.85,Q0=1)
print("The solution is:", x1)
print("The number of iterations is:", n1)
x2, n2, iterx3, itery3, Nf, Ng=Modified_Hager_Zhang(x0,F,dF,alpha=1,tol=1e-10, Niter=10000,nuk=0.85,Q0=1)
print("The solution is:", x2)
print("The number of iterations is:", n2)
x0 = np.array([1.6,-3.2])
F=Rastrigin
dF=dRastrigin
sol, iterx11, itery11,n, Nf, Ng=GrippoMethod(dF, F, x0, Niter=10000, tol=1e-4, M=3)
print("The solution is:", sol)
print("The number of iterations is:", n)
x1, iterx21, itery21,n1, Nf, Ng=Hager_Zhang(x0,F,dF,tol=1e-4, Niter=10000,nuk=0.85,Q0=1)
print("The solution is:", x1)
print("The number of iterations is:", n1)
x2, n2, iterx31, itery31, Nf, Ng=Modified_Hager_Zhang(x0,F,dF,alpha=1,tol=1e-10, Niter=10000,nuk=0.85,Q0=1)
print("The solution is:", x2)
print("The number of iterations is:", n2)
m = 30

fig, axes = plt.subplots(1, 2, figsize=(12, 5))

# First initial point
axes[0].plot(iterx1[:m], itery1[:m], 'o--', c='m', label='Algorithm 1')
axes[0].plot(iterx2[:m], itery2[:m], 'o--', c='b', label='Algorithm 2')
axes[0].plot(iterx3[:m], itery3[:m], 'o--', c='k', label='Algorithm 3')

axes[0].set_title('Initial point 1: [-0.5, 0.9]')
axes[0].set_xlabel('$k$')
axes[0].set_ylabel('$F(x_k)$')
axes[0].legend()
axes[0].grid(True)

# Second initial point
axes[1].plot(iterx11[:m], itery11[:m], 'o--', c='m', label='Algorithm 1')
axes[1].plot(iterx21[:m], itery21[:m], 'o--', c='b', label='Algorithm 2')
axes[1].plot(iterx31[:m], itery31[:m], 'o--', c='k', label='Algorithm 3')

axes[1].set_title('Initial point 2: x0=[1.6, -3.2]')
axes[1].set_xlabel('$k$')
axes[1].set_ylabel('$F(x_k)$')
axes[1].legend()
axes[1].grid(True)

plt.tight_layout()
plt.show()

