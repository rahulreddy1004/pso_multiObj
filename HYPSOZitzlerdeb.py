import numpy as np
import numpy.matlib as matlib
import random
from math import sin, sqrt
import matplotlib.pyplot as plt


def RouletteWheelSelection(P):
    r=random.random()
    C=np.cumsum(P)
    i=(np.array([np.nonzero(r<=C)]))[0][0][0]
    return i

           '''define cost function here'''
def CostFunction(x):
    res1 = -0.2*sqrt((x[0]*x[0]) + (x[1]*x[1]))
    res2 = -0.2*sqrt((x[1]*x[1]) + (x[2]*x[2]))
    f1 = -10.0*( 2.71828**(res1) + 2.71828**(res2))
    f2 = 0.0
    for i in range(3):
        f2 += pow(abs(x[i]),0.8) + 5.0*sin(pow(x[i],3.0))
    z=np.array([f1,f2])
    return z

class mty_grid:
    pass
     '''Creating the hypercubes'''

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
    GridUB=np.zeros((nObj, 2,nGrid+1))                                                                  
    inf_arr =np.array([ float('inf') for i in range(nGrid+1) ])
    for j in range(nObj):
        cj=np.linspace(cmin[0][j],cmax[0][j],num = nGrid+1)
        GridLB[j]=[-(inf_arr),cj]
        GridUB[j]=[cj,inf_arr]
    return GridLB,GridUB

     '''index the hypercubes according to densities'''
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
'''  Domination Rule '''

def Dominates(x,y):
    x=x.Cost   
    y=y.Cost                                                                      
    b=np.all(x<=y) and np.any(x<y)
    return b
''' Select Leader'''
def SelectLeader(rep,beta):
    # Grid Indx
    GI=np.array([re.GridIndex for re in rep])
    # Cells that are filled
    OC = np.unique(GI)
    # Number of Particles in Filled Cells
    N=np.zeros(list(np.shape(OC)))
    for k in range(len(OC)):
        ind = np.array([np.nonzero(GI==OC[k])])[0][0]
        indx = [i for i in ind]
        N[k]=len(indx)
    # Probabilities for selections
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
    sm=SCM[smi]
    leader=rep[sm]
    return leader
       ''' Determine domination the given array of partilcles'''

	
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
        '''to find the nearest neighbour'''

	
def neighbor(solution,T):
    nVar= len(x)
    j=np.random.randint(0, high=nVar)
    temp = x[j]
    if random.random()>0.5:
        x[j] = x[j] + 0.0003*T
    else:
        x[j] = x[j] - 0.0003*T
    if x[j] < VarMin or x[j] > VarMax:
        x[j] = temp
    return x

def acceptance_probability(old_cost, new_cost, T):
    ap = 0.05*T
    return ap

''' main function of simmulated annealing'''
def anneal(particle,T):
    particle.Cost = CostFunction(particle.Position)
    given_posiion = particle.Position
    T_min = 0.00001
    alpha = 0.9
    newSol = Particle()
    while T > T_min:
        i = 0
        while i <= 3000:
            new_position = neighbor(particle.Position,T)
            newSol.Position = new_position
            newSol.Cost = CostFunction(newSol.Position)
            ap = 0
            if newSol.Position[0] >VarMax or newSol.Position[0]<VarMin:
                i += 1
                continue
            if Dominates(newSol,particle):
                particle.Position=newSol.Position
                particle.Cost=newSol.Cost
                print "kkkklk"
            elif Dominates(particle,newSol):
                pass
            else:
                if random.random()< ap:
                    particle.Position=newSol.Position
                    particle.Cost=newSol.Cost
                    print "llllllk"
            '''    
            if particle.Position[0] >VarMax or particle.Position[0]<VarMin or particle.Position[1] >VarMax or particle.Position[1]<VarMin:
                print particle.Position,
                particle.Position = given_posiion
                particle.Cost = CostFunction(particle.Position)
                print 'shiffffffffffffffffffffffffffffffffffffffffft'
            '''
            i += 1
        T = T*alpha

	
nVar=3   #Number of Decision Variables
VarSize=np.array([1,nVar])  # Size of Decision Variables Matrix
VarMin=-5         # Lower Bound of Variables
VarMax=5    # Upper Bound of Variables
MaxIt= 200 # Maximum Number of Iterations
nPop  =1000 # Population Size
nRep=350         # Repository Size
w=0.9           # Inertia Weight
wdamp=0.98        # Intertia Weight Damping Rate
c1=1           # Personal Learning Coefficient
c2=1           # Global Learning Coefficient
nGrid=7           # Number of Grids per Dimension
alpha=0.1          # Inflation Rate
beta=2             # Leader Selection Pressure
gamma=2            # Deletion Selection Pressure
mu=0.1             # Mutation Rate
init_factor = 1
class Particle:
    pass
def initialise(p):
    p.Position=np.array([random.uniform(-5,5) for _ in range(nVar)] )
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
    initialise(p)
    p.Velocity=np.array([0 for _ in range(nVar)])
    p.Cost=CostFunction(p.Position)
for i in range(nPop):    # Update Personal Best
    Bestparticles[i].Position=particles[i].Position
    Bestparticles[i].Cost=particles[i].Cost
'''
x,y = [],[]
for r in particles:
    x.append(r.Cost[0])
    y.append(r.Cost[1])
m,n = [],[]
plt.plot(x,y, 'ro')
plt.axis([0, 140, 0,50])
plt.show()
'''
particles=DetermineDomination(particles)


pre_rep = []
for p in particles:
    pre_rep.append(not p.IsDominated)
particle_array = np.array(particles)
rep=particle_array[np.array(pre_rep)]
GridLB,GridUB=CreateGrid(rep,nGrid,alpha)


for i in range(len(rep)):
    rep[i]=FindGridIndex(rep[i],GridLB,GridUB)
for r in rep:
    print r.Cost,r.GridIndex
print SelectLeader(rep,beta).Cost


RepX = []
for i in rep:
    RepX.append(i.Position)
RepF = []
for i in rep:
    RepX.append(i.Cost)


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
    # Show Iteration Information
    print 'Iteration ',str(it),': Number of Rep Members = ',str(len(rep)),leader.Cost,leader.Position
    # Damping Inertia Weight
    w= 0.9 -(0.5) * (float(it)/MaxIt)


l,m = [],[]
for r in rep:
    l.append(r.Cost[0])
    m.append(r.Cost[1])
plt.plot(l,m, 'ro')
plt.axis([-20,-14,-12,2])
plt.show()                                               #plot and print MOPSO output
print 'rep length = ',len(rep)

for i in range (len(l)):
    print l[i] ,m[i]


dist_arr=[]

for p in particles:
    a = p.Cost[0]
    c = p.Cost[1]
    dis = []
    for r in rep:
        b = r.Cost[0]
        d = r.Cost[1]
        dis.append(sqrt((a - b)**2+(c - d)**2))     
    min(dis)
    dist_arr.append(min(dis))    
i=0
temperature_arr = []
for p in particles:
    p.t = dist_arr[i] 
    i=i+1
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
print particles[:10]

t1,t2,t3 = 0,0,0

for p in particles[0:(len(particles))/3]:
    t1 = t1 + (p.t)
for p in particles[(len(particles))/3:2*(len(particles))/3]:
    t2 = t2 + (p.t)
for p in particles[2*(len(particles))/3:]:
    t3 = t3 + (p.t)
print "before",t1,t2,t3

def normalize(t1,t2,t3): 
    if t3 > 15:
        t1 = t1/5 + 0.01
        t2 = t2/5 + 0.01
        t3 = t3/5 + 0.01
    if t3 < 1:
        t1 = t1 + 0.01
        t2 = (t2 + 0.01)*2
        t3 = (t3 + 0.01)*5           
    if t3 > 15:
        t1,t2,t3 = normalize(t1,t2,t3)
    return t1,t2,t3    

t1,t2,t3 = normalize(t1,t2,t3)


x,y = [],[]
for r in rep:
    x.append(r.Cost[0])
    y.append(r.Cost[1])
pa = 0
for r in particles[0:(len(particles))/3]:
    r = anneal(r,t1)             
    pa = pa+ 1
    print pa
for r in particles[(len(particles))/3:2*(len(particles))/3]:
    r = anneal(r,t2)
    pa = pa+ 1
    print pa
for r in particles[2*(len(particles))/3:]:
    r = anneal(r,t3)
    pa = pa+ 1
    print pa


a,b,c,d,e,f =[],[],[],[],[],[]
for r in particles[0:(len(particles))/3]:
    a.append(r.Cost[0])
    b.append(r.Cost[1])
for r in particles[(len(particles))/3:2*(len(particles))/3]:                      #plot temperature range particles output                   
    c.append(r.Cost[0])
    d.append(r.Cost[1])
for r in particles[2*(len(particles))/3:]:
    e.append(r.Cost[0])
    f.append(r.Cost[1])
plt.plot(a,b, 'ro',c,d,'bs',e,f,'g^')
plt.axis([-20,-14,-12,2])
plt.show()
particles=DetermineDomination(particles)
    # Keep only Non-Dminated Memebrs in the Repository
pre_rep2 = []

for p in particles:
    pre_rep2.append(not p.IsDominated)
pre_rep2    = np.array(pre_rep2) 
particles = np.array(particles)     
rep=particles[pre_rep2]
print 'annealed','sima',sima,'simb',simb

x,y = [],[]
for r in rep:
    x.append(r.Cost[0])
    y.append(r.Cost[1])
plt.plot(x,y, 'bs',l,m,'ro')
plt.axis([-20,-14,-12,2])
plt.show()
print 'rep length = ',len(rep)

for i in range (len(x)):
    print x[i] ,y[i] 
