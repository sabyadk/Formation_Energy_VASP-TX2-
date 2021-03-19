import numpy as np
import scipy
import matplotlib.pyplot as plt



def open_OUTCAR(natoms,arg='TOTEN',file_read='OUTCAR'):
 
    
    if arg=='TOTEN':
    
        TOTEN=read_TOTEN(file_read)
        
        print('total energy:'+str(TOTEN))
        return TOTEN
    if arg=='magnetization':
        
        mag=mag_read(file_read,natoms)
       
        return mag

def read_TOTEN(file_read):

   with open(file_read) as f:

       for line in reversed(f.readlines()):

           if line.rsplit('=')[0]=='  free  energy   TOTEN  ':
             print(line.rsplit('='))
             try:
              toten=float(line.rsplit('=')[1].split(' ')[-2])
             except ValueError:
              toten=float(line.rsplit('=')[1].split('  ')[-2]) 
             return toten

def mag_read(file_read,natoms):
     
 
  with open(file_read) as f:
  
         j=0
         mag=[]
         for line in f:
           #print line.strip() 
           if(line.strip()=='magnetization (x)'): 
             j=j+1
           if j>0:
             j=j+1
#           if j==5:
#             print 'magnetization (x)'
           if((j>5)&(j<(natoms+6))):
 
              mag.append(float(line.strip().split('  ')[-1]))
 #             print(line.strip())
         mag=np.array(mag)
         
         if len(mag)==0:
             print('no magnetization')
         return mag



