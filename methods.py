import numpy as np
def paszhanghager(d1,C0,F, dF, x, taun, c1=10**-4, Niter=5000):
 Fx, dFx = C0, dF(x)
# h=-dF(x)
 h=d1
 alpha=0.5
 i=0
 for n in range(Niter):
   i=i+1
   if F(x + taun*h) <= Fx + c1*taun*np.dot(dFx,h):
     return taun, i
   else:
    taun=alpha*taun
 print("Erreur, l’algorithme zhang-hager n'a pas convergé après", Niter, "itérations")
 def Modified_Hager_Zhang(x0,F,dF,alpha=1,tol=1e-4, Niter=1000,nuk=0.55,Q0=1):
   C0=F(x0)
   x1=x0
   tau=1
   iterx1=np.array([])
   itery1= np.array([])
   Ng=0
   Nf=0
   for n in range(Niter):
    gradFxn = dF(x1)
    er=np.linalg.norm(gradFxn)
    #print(er)
    if  er < tol:
     return x1,n,iterx1,itery1, Ng, Nf
    else:
     #H=Hess(x0)
     #d1= -np.linalg.solve(H, gradFxn)
     iterx1=np.append(iterx1,n)
     itery1=np.append(itery1,F(x0))
     lamda=1/alpha
     d1=-lamda*dF(x0)/tau
     Ng=Ng+1
     gr=dF(x0)
     [tau,i]=paszhanghager(d1,C0,F, dF,x0,1.05)
     Nf=Nf+i
     x1= x0 +tau*d1
     #print(x1)
     Q=nuk*Q0+1
     C=(nuk*Q0*C0+F(x1))/Q
     Nf=Nf+1
     gr1=dF(x1)
     Ng=Ng+1
     y=gr1-gr
     ss=np.dot(gr,gr)
     s1=np.dot(gr,y)
     sk=x1-x0
     #alphak=np.dot(sk,sk)/np.dot(y,sk)
     alpha=-(s1)/(lamda*ss)
     Q0=Q
     C0=C
     x0=x1
   print("Erreur, l’algorithme n’a pas convergé après ", Niter, "itérations")
def paszhanghager(d1,C0,F, dF, x, taun, c1=10**-4, Niter=5000):
 Fx, dFx = C0, dF(x)
# h=-dF(x)
 h=d1
 alpha=0.5
 i=0
 for n in range(Niter):
   i=i+1
   if F(x + taun*h) <= Fx + c1*taun*np.dot(dFx,h):
     return taun, i
   else:
    taun=alpha*taun
 print("Erreur, l’algorithme zhang-hager n'a pas convergé après", Niter, "itérations")
 def Hager_Zhang(x0,F,dF,tol=1e-4, Niter=1000,nuk=0.85,Q0=1):
   C0=F(x0)
   x1=x0
   iterx1=np.array([])
   itery1= np.array([])
   Ng=0
   Nf=0
   for n in range(Niter):
    gradFxn = dF(x1)
    Ng=Ng+1
    er=np.linalg.norm(gradFxn)
    #print(er)
    if  er < tol:
     return x1, iterx1, itery1,n, Nf, Ng
    else:
     d1=-dF(x0)
     [tau, i]=paszhanghager(d1,C0,F, dF,x0,1)
     Nf=Nf+i
     x1= x0 +tau*d1
     iterx1=np.append(iterx1,n)
     itery1=np.append(itery1,F(x0))
     #print(x1)
     Q=nuk*Q0+1
     C=(nuk*Q0*C0+F(x1))/Q
     Nf=Nf+1
     Q0=Q
     C0=C
     x0=x1
   print("Erreur, l’algorithme n’a pas convergé après ", Niter, "itérations")

def mafoncmax(v,n,M):
   k=min(n,M)
   if (k==n):
     return max(v[0:n+1])
     #(vnew)
   else:
      return max(v[n-M:n+1])


#Definition du pas optimal non monotone
def pasWolfeMax(fmax,F, dF, x, taun, c1=10**-4, Niter=5000):
 Fx, dFx = fmax, dF(x)
 h=-dF(x)
 alpha=0.5
 i=0
 xplus=x + taun*h
 for n in range(Niter):
   i=i+1
   if F(xplus) <= Fx + c1*taun*np.dot(dFx,h):
     return [taun,xplus,i]
   else:
    taun=alpha*taun
    xplus=x + taun*h
 print("Erreur, l’algorithme n'a pas convergé après", Niter, "itérations")
#Steepest descent method with backtracking line search non monotone

 def GrippoMethod(dF, F, x0, tol=1e-2, M=10, Niter=5000):
   xn1=x0
   iterx1=np.array([])
   itery1= np.array([])
   fn=np.empty(0)
   fmax=1
   Nf=0
   Ng=0
   for n in range(Niter):
    gradFxn1 = dF(xn1)
    Ng=Ng+1
    er1=np.linalg.norm(gradFxn1)
    if  er1 < tol:
     return xn1, iterx1, itery1,n, Nf, Ng
    else:
     fn=np.append(fn,F(xn1))
     Nf=Nf+1
     fmax=mafoncmax(fn,n,M)
     iterx1=np.append(iterx1,n)
     itery1=np.append(itery1,F(xn1))
     [tau1,xplus,i]=pasWolfeMax(fmax,F, dF, xn1,1)
     Nf=Nf+i
     xn1=xplus
     #xn1= xn1 -tau1*gradFxn1
    #itern1=np.append(itern1,k)
   # iterer1=np.append(iterer1,np.linalg.norm(xn1-xtar))
   print("Erreur, l’algorithme n’a pas convergé après ", Niter, "itérations")
