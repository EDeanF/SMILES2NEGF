from rdkit import Chem
import numpy as np 
from functools import reduce
import matplotlib.pyplot as plt

# input SMILES and linker element
SMILES_in = input("SMILES: ")
linker = input("linker element: ")

m = Chem.MolFromSmiles(SMILES_in)

# FIND LINKERS
# make list of possible linker atoms
linker_list = []
i = 0
for atom in m.GetAtoms():
    if atom.GetSymbol() == linker:
        linker_list.append(i)
    i += 1

if len(linker_list)<2:
    print("ERROR: did not find linker")
    quit()
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

# make adjacency matrix
AdjMatrix = np.array(Chem.GetAdjacencyMatrix(m))#,useBO=True))

# compute transmission function
# input parameters
eps = float(input("on-site energy: "))
tau = float(input("atomic coupling: "))
Gam = float(input("lead coupling: "))

N = AdjMatrix.shape[0]
# make Hamiltonian
H = tau*AdjMatrix+eps*np.identity(N)
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
plt.savefig(name+'.png')