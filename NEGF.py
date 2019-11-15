from rdkit import Chem
import numpy as np 
from functools import reduce
import matplotlib.pyplot as plt
import sys

# input SMILES and linker element
SMILES_in = input("SMILES: ")
linker = input("linker element: ")

m = Chem.MolFromSmiles(SMILES_in)

# electronegativity, w/ respect to carbon
EN = {'C':0,'N':-2.7391,'O':-5.3315,'S':0.1539}
# ionization energy w/ respect to carbon
IE = {'C':0,'N':-1.63692,'O':-1.17888,'S':0.45014}
# FIND LINKERS
# make list of possible linker atoms
linker_list = []
i = 0

atoms=[]
for atom in m.GetAtoms():
    atoms.append(atom.GetSymbol())
    if atom.GetSymbol() == linker:
        linker_list.append(i)
    i += 1

if len(linker_list)<2:
    print("ERROR: did not find linker")
    sys.exit()

# find atoms with the longest graph distance
DistMatrix = Chem.GetDistanceMatrix(m)
dist_max = 0
for i in range(len(linker_list)):
    j = linker_list[i]
    for k in linker_list[i+1:]:
        if DistMatrix[j][k] > dist_max:
            l = j
            r = k
            dist_max = DistMatrix[j][k]

# replace linkers with carbon?
DoReplace = input("replace linker with Carbon? (y/n) : ")
if DoReplace=='y':
    atoms[l]='C'
    atoms[r]='C'

# make adjacency matrix
useBO = input("use bond order? (y/n) ")
useBO = useBO=='y'
if(not useBO):
    print('not using bond order')

AdjMatrix = np.array(Chem.GetAdjacencyMatrix(m,useBO=useBO))
# compute transmission function
# input parameters
# eps = float(input("on-site energy: "))
Gam = float(input("lead coupling: "))
if(useBO):
    tau1 = float(input("single bond coupling: "))
    tau2 = float(input("double bond coupling: "))
    tau3 = float(input("triple bond coupling: "))
    AdjMatrix = tau1*(AdjMatrix==1)+tau2*(AdjMatrix==2)+tau3*(AdjMatrix==3)    
else:
    tau = float(input("atomic coupling: "))
    AdjMatrix = tau*AdjMatrix

N = AdjMatrix.shape[0]
# make Hamiltonian
# H = AdjMatrix+eps*np.identity(N)
H = AdjMatrix
OnsiteEnergy = input("Onsite energy type? (EN/IE) : ")
if OnsiteEnergy=='EN' or OnsiteEnergy=='IE':
    for i,atom in enumerate(atoms):
        if OnsiteEnergy=='EN' :
            H[i][i]=EN[atom]
        elif OnsiteEnergy=='IE' :
            H[i][i]=IE[atom]

# make Coupling matrices
GamL=np.zeros((N,N))
GamR=np.zeros((N,N))
GamL[l][l]=Gam
GamR[r][r]=Gam

# Green's function
# for each energy, calculate the Green's function matrix
maxE = float(input("maximum energy: "))
energy = np.linspace(-maxE,maxE,num=1000)
G_list = [np.linalg.inv(E*np.identity(N)-(H-0.5j*GamL-0.5j*GamR)) for E in energy]

# Calculate Transmission
T = [np.trace(reduce(np.matmul,[GamL,G,GamR,np.conjugate(G)])) for G in G_list]
T = list(map(np.absolute,T))

# write to file
name = input('molecule name: ')
file_name = name+'.txt'
f = open(file_name,'w')
for i in range(len(energy)):
    f.write("\t".join(map(str,(energy[i],T[i]))))
    f.write("\n")
f.close()

# plot transmission function
plt.plot(energy,T)
plt.yscale('log')
plt.xlabel('Energy')
plt.ylabel('Transmission')
plt.ylim(1e-5,1)
plt.savefig(name+'.png')