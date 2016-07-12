# -*- coding: utf-8 -*-
"""
Created on Tue Jul 12 15:32:26 2016

@author: rahul_reddy
"""
import numpy as np
import numpy.matlib as matlib
import random
from math import sin, sqrt
import matplotlib.pyplot as plt
''' Problem Definition'''

 
def MOP2(x):


    n=len(x)
    
    z1=1-np.exp(-np.sum(np.power((x-1/sqrt(n)),2)))                                             #ckeck exp for single
    
    z2=1-np.exp(-np.sum((x+1/sqrt(n))**2))
    
    z=np.array([z1,z2])
    return z
def MOP4(x):

    a=0.8
    
    b=3
    
    z1=np.sum(-10*exp(-0.2*sqrt(np.power(x[0:-1],2)+np.power(x[1:],2))))
    
    z2=np.sum(np.power(np.absolute(x),a)+5*(np.power((np.sin(x)),b)))
    
    z=np.array([z1,z2])
    return z

def RouletteWheelSelection(P):

    r=random.random()
    
    C=np.cumsum(P)
    
    i=(np.array([np.nonzero(r<=C)]))[0][0][0]
    
    return i


def CostFunction(x):
    n=len(x)
    f1=4*(x[0]**2)+4*(x[1]**2)
 
    
    den = (n-1)*np.sum(x[1:])
   
    
    
    
    f2 = (x[0]-5)**2+(x[1]-5)**2
    z=np.array([f1,f2])
    
    
    return z
class mty_grid:
    pass
def CreateGrid(pop,nGrid,alpha):

    c=[]
    
    for p in pop:
        c.append(p.Cost)
    c=np.array([c])
    cmin=np.amin(c, axis=1)
    cmax=np.amax(c, axis=1)
    
    dc=cmax-cmin
    cmin=cmin-alpha*dc
    cmax=cmax+alpha*dc
    
    nObj=list(np.shape(c))[2]
    
    
    GridLB=np.zeros((nObj, 2,nGrid+1))
    GridUB=np.zeros((nObj, 2,nGrid+1))                                                                  # change range(1.+1) to 0,+0
    inf_arr =np.array([ float('inf') for i in range(nGrid+1) ])
    
    for j in range(nObj):
        
        cj=np.linspace(cmin[0][j],cmax[0][j],num = nGrid+1)
        
        GridLB[j]=[-(inf_arr),cj]
        GridUB[j]=[cj,inf_arr]
        
    
    return GridLB,GridUB
def FindGridIndex(particle,GridLB,GridUB):

    nObj=len(particle.Cost)
    
    nGrid=len(GridLB[0])
    
    particle.GridSubIndex=np.zeros((1,nObj))
    
    
    for j in range(0,nObj):
        
        
        particle.GridSubIndex[0][j]= (np.array([np.nonzero(particle.Cost[j]<=GridUB[j][0])])[0][0][:])[0]
       

            
                 
    
    
    particle.GridIndex=particle.GridSubIndex[0][0]
    
    for j in range(0,nObj):
        particle.GridIndex=particle.GridIndex-1
        particle.GridIndex=nGrid*particle.GridIndex
        particle.GridIndex=particle.GridIndex+particle.GridSubIndex[0][j]
   
    return particle
    

def Mutate(x,pm,VarMin,VarMax):

    nVar=len(x)
    j=np.random.randint(0, high=nVar)

    dx=pm*(VarMax-VarMin)
    
    lb=x[j]-dx
    if lb<VarMin:
        lb=VarMin
    
    
    ub = x[j]+dx
    if ub>VarMax:
        ub=VarMax
    
    
    xnew = x
    xnew[j]=np.random.uniform(lb,ub)

    return xnew


def Dominates(x,y):
    

    


    x=x.Cost
    
    
   
    y=y.Cost                                                                      #ckeck again
    

    b=np.all(x<=y) and np.any(x<y)
    return b
def SelectLeader(rep,beta):

    # Grid Index of All Repository Members
    GI=np.array([re.GridIndex for re in rep])

    # Occupied Cells
    
    OC = np.unique(GI)
    
    # Number of Particles in Occupied Cells
    N=np.zeros(list(np.shape(OC)))
    
    
    for k in range(len(OC)):
        ind = np.array([np.nonzero(GI==OC[k])])[0][0]
        indx = [i for i in ind]
       
        N[k]=len(indx)
        
   
        
    # Selection Probabilities
    P=np.exp(-beta*N)
    P=P/np.sum(P)
   
    # Selected Cell Index
    sci=RouletteWheelSelection(P)
   
    # Selected Cell
    sc=OC[sci]
   
    # Selected Cell Members
    SCM=np.array([np.nonzero(GI==sc)])[0][0][:]
    
    
    # Selected Member Index
    if len(SCM) == 0:
        smi = 0
    else:
        smi=np.random.randint(0,high = len(SCM))
    
    # Selected Member
    sm=SCM[smi]
    
    # Leader
    leader=rep[sm]

    return leader
def DeleteOneRepMemebr(rep,gamma):
    repd = [re for re in rep]
    #Grid Index of All Repository Members
    GI=np.array([re.GridIndex for re in repd])
    
    # Occupied Cells
    OC=np.unique(GI)
    
    # Number of Particles in Occupied Cells
    N=np.zeros(list(np.shape(OC)))
    for k in range(len(OC)):
        N[k]=len(np.array([np.nonzero(GI==OC[k])])[0][0][:])
    
    # Selection Probabilities
    P=np.exp(gamma*N)
    
    # Selected Cell Index
    sci=RouletteWheelSelection(P)
    
    # Selected Cell
    sc=OC[sci]
    
    # Selected Cell Members
    SCM=np.array([np.nonzero(GI==OC[k])])[0][0][:]
    
    # Selected Member Index
    smi=np.random.randint(1,high = len(SCM))
    # Selected Member
    sm=SCM[smi]
    
    # Delete Selected Member
    rep[sm]=np.array([])

    return rep  
def DetermineDomination(particles):
    nPop = len (particles)
    for p in particles:
        p.IsDominated=False
    for i in range(nPop-1):
        for j in range(i,nPop):
            if Dominates(particles[i],particles[j]):
               particles[j].IsDominated=True
            if Dominates(particles[j],particles[i]):
               particles[i].IsDominated=True
    return particles

def neighbor(solution):
     a =random. random()
     if random.random() > 0.5:
 
         
         if a > 0.5:
             solution[0] = solution[0]-0.01
         else:
             solution[0] = solution[0]+0.01
 
     else:
         if a > 0.5:
             solution[1] = solution[1]-0.01
         else:
             solution[1] = solution[1]+0.01
 
     return solution
 
 
 
def acceptance_probability(old_cost, new_cost, T):
	ap = 0
	return ap
def anneal(particle,T):
    particle.Cost = CostFunction(particle.Position)
    given_posiion = particle.Position
    
    T_min = 0.00001
    alpha = 0.9
    newSol = Particle()
    while T > T_min:
        i = 1
        while i <= 100:
            new_position = neighbor(particle.Position)
            newSol.Position = new_position
            newSol.Cost = CostFunction(newSol.Position)
            ap = 0.5*T
            if Dominates(newSol,particle):
                
                particle.Position=newSol.Position
                particle.Cost=newSol.Cost
            elif Dominates(particle,newSol):
                pass
            else:
                if random.random()> ap:
                    particle.Position=newSol.Position
                    particle.Cost=newSol.Cost
                                        
                
                    
            '''   
            if particle.Position[0] >VarMax or particle.Position[0]<VarMin or particle.Position[1] >VarMax or particle.Position[1]<VarMin:
                print particle.Position,
                particle.Position = given_posiion
                particle.Cost = CostFunction(particle.Position)
                print 'shiffffffffffffffffffffffffffffffffffffffffft'
            '''
            i += 1
        T = T*alpha
    return particle
nVar=2    #Number of Decision Variables
VarSize=np.array([1,nVar])  # Size of Decision Variables Matrix
VarMin=0          # Lower Bound of Variables
VarMax=100        # Upper Bound of Variables
MaxIt=40  # Maximum Number of Iterations
nPop  = 100 # Population Size
nRep=100         # Repository Size
w=0.5              # Inertia Weight
wdamp=0.98        # Intertia Weight Damping Rate
c1=0.5               # Personal Learning Coefficient
c2=1               # Global Learning Coefficient
nGrid=7           # Number of Grids per Dimension
alpha=0.1          # Inflation Rate
beta=2             # Leader Selection Pressure
gamma=2            # Deletion Selection Pressure
mu=0.1             # Mutation Rate
init_factor = 100

class Particle:
    pass
''' Initialization '''
particles = []
Bestparticles =[]
for i in range(nPop):
    p = Particle()
    p.Position=np.array([])
    
    p.Velocity=np.array([])
    p.t = 0
    p.Cost = np.array([])
    p.IsDominated =np.array([])
    p.GridIndex = np.array([])
    p.GridSubIndex =np.array([])
    particles.append(p)
for i in range(nPop):
    p = Particle()
    p.Position=np.array([])
    p.Cost = np.array([])
    Bestparticles.append(p)
for p in particles:
    p.Position=np.array([init_factor*random.random() for _ in range(nVar)] )
    
    p.Velocity=np.array([0 for _ in range(nVar)])
    
    p.Cost=CostFunction(p.Position)
    
    
for i in range(nPop):    # Update Personal Best
    Bestparticles[i].Position=particles[i].Position
    Bestparticles[i].Cost=particles[i].Cost
x,y = [],[]
for r in particles:
    x.append(r.Cost[0])
    y.append(r.Cost[1])
plt.plot(x,y, 'ro')
plt.axis([0, 100, 0, 100])
plt.show()

particles=DetermineDomination(particles)
pre_rep = []
for p in particles:
    pre_rep.append(not p.IsDominated)

particle_array = np.array(particles)
rep=particle_array[np.array(pre_rep)]

x,y = [],[]
for r in rep:
    
    x.append(r.Cost[0])
    
    y.append(r.Cost[1])
plt.plot(x,y, 'ro')
plt.axis([0, 50, 0, 50])
plt.show()    
GridLB,GridUB=CreateGrid(rep,nGrid,alpha)
for i in range(len(rep)):
    rep[i]=FindGridIndex(rep[i],GridLB,GridUB)
for r in rep:
    print r.Cost,r.GridIndex

print SelectLeader(rep,beta).Cost

# MOPSO Main Loop

for it in range(1,MaxIt+1):

    
    for i in range(nPop):
        leader=SelectLeader(rep,beta)
        
        particles[i].Velocity = w*particles[i].Velocity  \
            +np.multiply(c1*np.random.rand(1,nVar),(Bestparticles[i].Position-particles[i].Position)) \
            +np.multiply(c2*np.random.rand(1,nVar),(leader.Position-particles[i].Position))
        
        particles[i].Position = particles[i].Position + particles[i].Velocity[0]
        
        
        particles[i].Position = np.maximum(particles[i].Position, VarMin)
        particles[i].Position = np.minimum(particles[i].Position, VarMax)
        
        particles[i].Cost = CostFunction(particles[i].Position)

                    
        
        
        if Dominates(particles[i],Bestparticles[i]):
            
            a = particles[i].Position
            Bestparticles[i].Position=a
            b = particles[i].Cost
            Bestparticles[i].Cost=b
            
        elif Dominates(Bestparticles[i],particles[i]):
            pass
            
        else:
            if random.random()<0.5:
                
                a = particles[i].Position
                
                Bestparticles[i].Position=a
                b = particles[i].Cost
                Bestparticles[i].Cost=b
        
     
    particles = DetermineDomination(particles)    
    dom_particles = np.array([ not p.IsDominated for p in particles])
    # Add Non-Dominated Particles to REPOSITORY
    parray = np.array(particles)
    
    rep=np.array(parray[dom_particles]) #ok
    
    
    # Determine Domination of New Resository Members
    rep=DetermineDomination(rep)
    
    # Keep only Non-Dminated Memebrs in the Repository
    pre_rep1 = []
    
    for p in rep:
        pre_rep1.append(not p.IsDominated)
    pre_rep1    = np.array(pre_rep1) 
    
 
        
    rep=rep[pre_rep1]
    
    
       
    # Update Grid
    GridLB,GridUB=CreateGrid(rep,nGrid,alpha)

    # Update Grid Indices
    for i in range(0,len(rep)):
        rep[i]=FindGridIndex(rep[i],GridLB,GridUB)
    
    '''
    # Check if Repository is Full
    if len(rep)>nRep:
        
        Extra=len(rep)-nRep
        for e in range(1,Extra+1):
            rep=DeleteOneRepMemebr(rep,gamma)
       
    '''    
    
    
    x,y = [],[]
    for r in particles:
        x.append(r.Cost[0])
        y.append(r.Cost[1])
    plt.plot(x,y, 'ro')
    plt.axis([0, 120, 0, 50])
    plt.show()
            
               
   
 
    
    # Show Iteration Information
    print 'Iteration ',str(it),': Number of Rep Members = ',str(len(rep)),leader.Cost,leader.Position
    
    # Damping Inertia Weight
    w= w*wdamp
temperature_arr = []
for p in particles:
    p.t = sqrt(p.Velocity[0][0]**2+p.Velocity[0][1]**2)
    
    temperature_arr.append(p.t)
for fillslot in range(len(temperature_arr)-1,0,-1):
    positionOfMax=0
    for location in range(1,fillslot+1):
        if temperature_arr[location]>temperature_arr[positionOfMax]:
            positionOfMax = location

    temp = temperature_arr[fillslot]
    temp1 = particles[fillslot]
    temperature_arr[fillslot] = temperature_arr[positionOfMax]
    particles[fillslot] = particles[positionOfMax]
    temperature_arr[positionOfMax] = temp
    particles[positionOfMax] = temp1

t1,t2,t3 = 0,0,0
for p in particles[0:(len(particles))/3]:
    t1 = t1 + (p.t)**2
for p in particles[(len(particles))/3:2*(len(particles))/3]:
    t2 = t2 + (p.t)**2
for p in particles[2*(len(particles))/3:]:
    t3 = t3 + (p.t)**2
print t1,t2,t3
x,y = [],[]
for r in rep:
    x.append(r.Cost[0])
    y.append(r.Cost[1])
plt.plot(x,y, 'ro')
plt.axis([0, 120, 0, 50])
plt.show()
for r in particles[0:(len(particles))/3]:
    r = anneal(r,t1)
for r in particles[(len(particles))/3:2*(len(particles))/3]:
    r = anneal(r,t2)
for r in particles[2*(len(particles))/3:]:
    r = anneal(r,t3)
a,b,c,d,e,f =[],[],[],[],[],[]
for r in particles[0:(len(particles))/3]:
    a.append(r.Cost[0])
    b.append(r.Cost[1])
for r in particles[(len(particles))/3:2*(len(particles))/3]:
    c.append(r.Cost[0])
    d.append(r.Cost[1])
for r in particles[2*(len(particles))/3:]:
    e.append(r.Cost[0])
    f.append(r.Cost[1])
plt.plot(a,b, 'ro',c,d,'bs',e,f,'g^')
plt.axis([0, 120, 0, 50])
plt.show()
particles=DetermineDomination(particles)
    
    # Keep only Non-Dminated Memebrs in the Repository
pre_rep2 = []
    
for p in particles:
    pre_rep2.append(not p.IsDominated)
pre_rep2    = np.array(pre_rep2) 
    
 
particles = np.array(particles)     
rep=particles[pre_rep2]
print 'annealed'
x,y = [],[]
for r in rep:
    x.append(r.Cost[0])
    y.append(r.Cost[1])
plt.plot(x,y, 'ro')
plt.axis([0, 120, 0, 50])
plt.show()

