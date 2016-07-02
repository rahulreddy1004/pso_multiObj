import numpy as np
import random
from math import sin, sqrt
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

    r=random()
    
    C=np.cumsum(P)
    
    i=(np.array([np.nonzero(r<=c)[0][0][:]]))[0] 
    return i


def CostFunction(x):
    n=len(x)
    f1=x[0]
    g=  1+9/(n-1)*np.sum(x[1:end])
    h=1-sqrt(f1/g)
    f2=g*h
    z=np.array([f1,f2])
    return z
def CreateGrid(pop,nGrid,alpha):

    c=[]
    for p in pop:
        c.append(p.Cost)
    
    cmin=np.amax(c, axis=1)
    cmax=np.amin(c, axis=1)
    
    dc=cmax-cmin
    cmin=cmin-alpha*dc
    cmax=cmax+alpha*dc
    
    nObj=list(np.shape(c))[0]
    
    empty_grid.LB=[]
    empty_grid.UB=[]
    Grid=numpy.matlib.repmat(empty_grid,nObj,1)
                                                                                         # change range(1.+1) to 0,+0
    for j in range(nObj):
        
        cj=np.linspace(cmin[j],cmax[j],nGrid+1)
        Grid[j].LB=np.array([-inf,cj])
        Grid[j].UB=np.array([cj,+inf])
        
    return Grid
def FindGridIndex(particle,Grid):

    nObj=len(particle.Cost)
    
    nGrid=len(Grid[0].LB)
    
    particle.GridSubIndex=np.zeros(1,nObj)
    
    for j in range(0,nObj):
        
        particle.GridSubIndex[j]= (numpy.array([numpy.nonzero(particle.Cost[j]<Grid[j].UB)])[0][0][:])[0]
    

            
                 
    

    particle.GridIndex=particle.GridSubIndex[0]
    for j in range(1,nObj):
        particle.GridIndex=particle.GridIndex-1
        particle.GridIndex=nGrid*particle.GridIndex
        particle.GridIndex=particle.GridIndex+particle.GridSubIndex[j]
    return particle
    

def Mutate(x,pm,VarMin,VarMax):

    nVar=len(x)
    j=np.array([random.random() for _ in range(nVar)] )

    dx=pm*(VarMax-VarMin)
    
    lb=x[j]-dx
    if lb<VarMin:
        lb=VarMin
    
    
    ub = x[j]+dx
    if ub>VarMax:
        ub=VarMax
    
    
    xnew = x
    xnew[j]=unifrnd(lb,ub)

    return xnew


def Dominates(x,y):
    

    if type(x) is struct:


        x=x.Cost
    
    
    if type(y) is struct:
        y=y.Cost
    

    b=np.all(x<=y) and np.any(x<y)
    return b
def SelectLeader(rep,beta):

    # Grid Index of All Repository Members
    GI=np.array([rep.GridIndex])
    
    # Occupied Cells
    OC=np.unique(GI)
    
    # Number of Particles in Occupied Cells
    N=np.zeros(list(np.shape(OC)))
    
    
    for k in range(len(OC)):
        N(k)=len(numpy.array([numpy.nonzero(GI==OC(k))])[0][0][:])
    # Selection Probabilities
    P=np.exp(-beta*N)
    P=P/np.sum(P)
    
    # Selected Cell Index
    sci=RouletteWheelSelection(P)
    
    # Selected Cell
    sc=OC[sci]
    
    # Selected Cell Members
    SCM=numpy.array([numpy.nonzero(GI==sc)])[0][0][:]
    
    
    # Selected Member Index
    smi=np.random.randini(1,high = len(SCM))
    
    # Selected Member
    sm=SCM[smi]

    # Leader
    leader=rep[sm]
    
def DeleteOneRepMemebr(rep,gamma):

    #Grid Index of All Repository Members
    GI=np.array([rep.GridIndex])
    
    # Occupied Cells
    OC=numpy.unique(GI)
    
    # Number of Particles in Occupied Cells
    N=np.zeros(list(np.shape(OC)))
    for k in range(len(OC)):
        N(k)=len(numpy.array([numpy.nonzero(GI==OC(k))])[0][0][:])
    
    # Selection Probabilities
    P=np.exp(gamma*N)
    
    # Selected Cell Index
    sci=RouletteWheelSelection(P)
    
    # Selected Cell
    sc=OC[sci]
    
    # Selected Cell Members
    SCM=numpy.array([numpy.nonzero(GI==OC(k))])[0][0][:]
    
    # Selected Member Index
    smi=np.random.randini(1,high = len(SCM))
    # Selected Member
    sm=SCM[smi]
    
    # Delete Selected Member
    rep[sm]=np.array([])

    return rep  
def DetermineDomination(particles):
    nPop = len (particles)
    for p in particles:
        p.IsDominated=false
    for i in range(nPop-1):
        for j in range(i,nPop):
            if Dominates(particles[i],particles[j]):
               particles[j].IsDominated=true
            if Dominates(particles[j],particles[i]):
               particles[i].IsDominated=true
    return particles


nVar=5             #Number of Decision Variables
VarSize=np.array([1,nVar])  # Size of Decision Variables Matrix
VarMin=0          # Lower Bound of Variables
VarMax=1          # Upper Bound of Variables
MaxIt=200           # Maximum Number of Iterations
nPop=200           # Population Size
nRep=100           # Repository Size
w=0.5              # Inertia Weight
wdamp=0.99         # Intertia Weight Damping Rate
c1=1               # Personal Learning Coefficient
c2=2               # Global Learning Coefficient
nGrid=7           # Number of Grids per Dimension
alpha=0.1          # Inflation Rate
beta=2             # Leader Selection Pressure
gamma=2            # Deletion Selection Pressure
mu=0.1             # Mutation Rate
class Particle:
    pass
''' Initialization '''
particles = []
for i in range(pop_size):
    p = Particle()
    p.Position=np.array([])
    p.Best.Position=np.array([])
    p.Velocity=np.array([])
    p.Best.Cost = np.array([])
    p.Cost = np.array([])
    p.IsDominated =np.array([])
    p.GridIndex = np.array([])
    p.GridSubIndex =np.array([])
    particles.append(p)
for p in particles:
    p.Position=np.array([random.random() for _ in range(nVar)] )
    
    p.Velocity=np.array([0 for _ in range(nVar)])
    
    p.Cost=CostFunction(p.Position)
    
    
    # Update Personal Best
    p[i].Best.Position=p[i].Position;
    p[i].Best.Cost=p[i].Cost;
# Determine Domination
particles=DetermineDomination(particles)
rep=pop(~[pop.IsDominated]);
Grid=CreateGrid(rep,nGrid,alpha);
for i in range(len(rep)):
    rep(i)=FindGridIndex(rep(i),Grid);



# MOPSO Main Loop

for it in range(1,MaxIt+1):

    
    for i in range(nPop):
        
        leader=SelectLeader(rep,beta)
        
        particles[i].Velocity = w*particles[i].Velocity  \
            +np.multiply(c1*rand(VarSize),(pop(i).Best.Position-pop(i).Position)) \
            +np.multiply(c2*rand(VarSize),(leader.Position-pop(i).Position))
        
        particles[i].Position = particles[i].Position + particles[i].Velocity
        
        particles[i].Position = np.maximum(particles[i].Position, VarMin)
        particles[i].Position = np.minimum(particles[i].Position, VarMax)
        
        particles[i].Cost = CostFunction(particles[i].Position)

        NewSol = Particle()
        # Apply Mutation
        pm=(1-(it-1)/(MaxIt-1))**(1/mu)
        if rand<pm:

            NewSol.Position=Mutate(particles[i].Position,pm,VarMin,VarMax);
            NewSol.Cost=CostFunction(NewSol.Position)
            if Dominates(NewSol,particles[i]):
                particles[i].Position=NewSol.Position
                particles[i].Cost=NewSol.Cost

            elif Dominates(particles[i],NewSol):
                pass

            else:

                if rand<0.5:
                    particles[i].Position=NewSol.Position
                    pop[i].Cost=NewSol.Cost
                
            
        
        
        if Dominates(particles[i],particles[i].Best):
            particles[i].Best.Position=particles[i].Position
            particles[i].Best.Cost=particles[i].Cost
            
        elif Dominates(pop(i).Best,pop(i)):
            pass
            
        else:
            if rand<0.5:
                particles[i].Best.Position=particles[i].Position
                particles[i].Best.Cost=particles[i].Cost
            
        
        
    
    
    # Add Non-Dominated Particles to REPOSITORY
    rep=np.array([rep,pop(~[pop.IsDominated])]) #ok
    
    # Determine Domination of New Resository Members
    rep=DetermineDomination(rep)
    
    # Keep only Non-Dminated Memebrs in the Repository
    rep=rep(~[rep.IsDominated])
    
    # Update Grid
    Grid=CreateGrid(rep,nGrid,alpha)

    # Update Grid Indices
    for i in range(1,num+1):
        rep[i]=FindGridIndex(rep[i],Grid)
    
    
    # Check if Repository is Full
    if len(rep)>nRep:
        
        Extra=len(rep)-nRep
        for e in range(1,Extra+1):
            rep=DeleteOneRepMemebr(rep,gamma)
        
        

    
    # Plot Costs
    
    
    # Show Iteration Information
    print 'Iteration ',num2str(it),': Number of Rep Members = ',num2str(numel(rep))
    
    # Damping Inertia Weight
    w= w*wdamp
    


