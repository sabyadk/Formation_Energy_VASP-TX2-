import numpy as np
import scipy
import matplotlib.pyplot as plt
import cmath

def area(lattice_vec,direction):

   print direction
   t=[0,1,2]
   for i in range(len(t)):
     if t[i] not in direction:
        X=t[i]
     
   V=np.dot(lattice_vec[X,:],(np.cross(lattice_vec[direction[0],:],lattice_vec[direction[1],:])))
   return(V/np.linalg.norm(lattice_vec[X,:]))




def mat_calc(atom_pos,natoms,materials):
  near_neigh=[]
  mat=[]
  l=0
  prev=0

  for i in range(len(atom_pos[:,0])):
    B=[]
    if i>=natoms[l]+prev:
      prev=prev+natoms[l]
      l=l+1

    if(i<(natoms[l]+prev)):
      mat.append(materials[l])

  return(mat)




def extract_POSCAR(file1,direction):
  natoms2=0
  lattice_vec=np.zeros((3,3),dtype=float)



  with open(file1) as f:
    j=0
    for line in f:

      if j==1:
        print(line.split()) 
        a=float(line.split()[0])

      elif((j>1)&(j<5)):
         lattice_vec[j-2,0]=float(line.split()[0])
         lattice_vec[j-2,1]=float(line.split()[1])
         lattice_vec[j-2,2]=float(line.split()[2])
      elif j==5:
          materials=np.array(line.split())
      elif j==6:
          natoms=np.array(line.split()).astype('int')

          for i in range(len(natoms)):
              natoms2=natoms[i]+natoms2
          atom_pos=np.zeros((natoms2,3),dtype=float)
      elif((j>7)&(j<8+natoms2)):
        atom_pos[j-8,0]=float(line.split()[0])
        atom_pos[j-8,1]=float(line.split()[1])
        atom_pos[j-8,2]=float(line.split()[2])

      j=j+1
  mat_typ=int(len(materials))
  atom_pos2=np.zeros((len(atom_pos[:,0]),len(atom_pos[0,:])),dtype=float)

  for j in range(len(atom_pos[:,0])):


        atom_pos2[j,:]=atom_pos[j,0]*lattice_vec[0,:]*a+atom_pos[j,1]*lattice_vec[1,:]*a+atom_pos[j,2]*lattice_vec[2,:]*a



  plt.show()
  mat=mat_calc(atom_pos,natoms,materials)
#  near_neigh=[]
#  for i in range(len(min_bond)-1):
 #    near_neigh.append(N_near_neigh_calc(atom_pos,atom_pos2,lattice_vec,a,min_bond[i],min_bond[i+1]))



#
  #near_neigh=np.array(near_neigh)
#  print near_neigh[2,57]
  A=area(lattice_vec,direction)
  return(atom_pos2,mat,materials,natoms,A)


