from random import random,seed,expovariate
from pylab import array
from numpy import linspace,zeros
from matplotlib.pyplot import *
import sympy as sym
from scipy.integrate import odeint
from scipy.linalg import solve_lyapunov
from numpy.linalg import inv
from numpy import cumsum

#simulation d'une chaine de markov continue

class Error(Exception):
    """Base class for eceptions in this module."""
    pass

class NullRate(Error):
    def __init__(self,msg):
        self.msg = msg

def simu_cmtc(taux,L,X,N,time,file="",fix=-1):
    n=len(X)
    nb_trans=len(L)
    t=0
    
    if fix!=-1:
        seed(fix)

    if file!="":
        data=open(file,"w")
    else:
        data=[[0],[X]]
        
    while t<time:
        L_poids=[taux(i,X) for i in range(nb_trans)]
        S=sum(L_poids)
        if S<=1e-14:
            print('System stalled (total rate = 0)')
            t = time
            data[1].append(X)
        else:
            a=random()*S
            l=0
            while a > L_poids[l]:
                a -= L_poids[l]
                l += 1
        
            X = X+(1./N)*L[l]

            t+=expovariate(N*S)

            if file=="":
                data[0].append(t)
                data[1].append(X[:])
            else:
                data.write(str(t)+" ")
                for i in range(n):
                    data.write(str(X[i])+" ")
                    data.write('\n')
    
    if file=="":
        data[1] = array(data[1])
        return(data)
    else:
        data.close()
        return(0)

#récupère les données dans un fichier et renvoit un tableau

def file_to_array(file,n):
    fichier=open(file,"r")
    texte=fichier.read()
    
    tmp=texte.split( )
    tmp=[tmp[i].split(' ',1) for i in range(len(tmp))]
    data=[[],[]]
    for i in range(int(len(tmp)/(n+1))):
        data[0].append(float(tmp[i*(n+1)][0]))
        data[1].append([float(x[0]) for x in tmp[i*(n+1)+1:i*(n+1)+n+1]])
    data[1]=array(data[1])
    return(data)


#Calcule le point fixe a l'aide de la fonction odeint

def fixed_point(taux,liste_transitions,n):
    number_transitions = len(liste_transitions)    
    def drift(x):
        return (sum([liste_transitions[i]*taux(i,x) for i in range(number_transitions)],0))
    
    X0 = zeros(n)
    X0[0] = 1
    t = linspace(0,1000,1000)
    x = odeint( lambda x,t : drift(x), X0, t)
    #figure()
    #plot(x)
    #show()
    return(x[-1,:])


#calcule le coefficient C_i théorique

def theorique(taux,liste_transitions,n,dim):
    number_transitions = len(liste_transitions)  
    X_0 = fixed_point(taux,liste_transitions,n)
    
    Var=array([sym.symbols('x_{}'.format(i)) for i in range(n)])
    
    #print(len(X_0))
    f_x=array([0 for i in range(n)])
    for i in range(number_transitions):
        f_x = f_x + liste_transitions[i]*taux(i,Var)

    if dim==n-1:
        for i in range(n):
            f_x[i]=f_x[i].subs(Var[-1],1-sum(array([Var[i] 
                                                for i in range(n-1)])))

    A=array([[sym.lambdify(Var ,sym.diff(f_x[i],Var[j]))(*[X_0[k] 
                for k in range(n)]) 
              for j in range(dim)] 
             for i in range(dim)])

    B=array([[[sym.lambdify(Var,sym.diff(f_x[j],Var[k],Var[l]))(*[X_0[i] 
                                                                  for i in range(n)]) 
               for l in range(dim)] 
              for k in range(dim)] 
             for j in range(dim)])

    Q=array([[0. for i in range(dim)] for j in range(dim)])

    for l in range(number_transitions):
        Q += array([[liste_transitions[l][p]*liste_transitions[l][m]*taux(l,X_0) 
                     for m in range(dim)] 
                    for p in range(dim)])


    #print('A=',A)
    #print('Q=',Q)
    #print('B=',B)
    
    W = solve_lyapunov(A,Q)
    #print('W=',W)
    #print('sumW=',sum(W,0))

    A_inv=inv(A)
    
    C=[ 0.5*sum([A_inv[i][j]*sum(array([[B[j][k_1][k_2]*W[k_1][k_2] 
                                         for k_2 in range(dim)] 
                                        for k_1 in range(dim)])) 
                 for j in range(dim)]) 
       for i in range(dim)]
    for i in range(len(C)):
        C[i]=sum(C[i])
    return(C)

def print_simu(taux,liste_transition,X,N,time,fix=-1):
    data=simu_cmtc(taux,liste_transition,X,N,time,"",fix)
    K=len(X)
    for i in range(K):
        p=plot(data[0],data[1][:,i],label=str(i),color='C{}'.format(i))
    def drift(x):
        return (sum([liste_transition[i]*taux(i,x) for i in range(len(liste_transition))],0))
    
    t = linspace(0,time,1000)

    x = odeint( lambda x,t : drift(x), X, t)
    for i in range(K):
        plot(t,x[:,i],color='C{}'.format(i),linestyle='--')
    legend()
    show()
