import numpy as np
import scipy 
import matplotlib.pyplot as plt
import argparse
from read_OUTCAR import read_TOTEN
from extract_position import extract_POSCAR

def substance(T):


 metal=['Fe','Sc','V','Ti','Cr','Mn','Ni','Co','Cu']
  
 tran_metal=['W','Mo']
  
 chalcogen=['Se','S','Te']

 if T in metal:
   
    return 'M'
  
 elif T in tran_metal:
 
    return 'T'

 elif T in chalcogen:

    return 'X' 


def read_files(strings,direction):

     l=0
     L=len(strings)
 
     
     for f in strings:

            if l==0:

               atom_pos1,mat1,materials1,natoms1,A=extract_POSCAR(f,direction)

            if l==1:
               print f
               Et=read_TOTEN(f)

       
            l=l+1

     return(materials1,natoms1,Et,A)  


def plot_line(E,mu_X,A,font=16):
  from matplotlib import colors, ticker, cm
  from matplotlib.mlab import bivariate_normal
  import matplotlib as mpl
  mpl.rcParams['axes.linewidth'] = 1.5
  from matplotlib import rc

  rc('text',usetex=True)
  rc('font',family='serif')



  from matplotlib import rcParams
  rcParams.update({'figure.autolayout': True})



  legend=[]
  colors=['b','g','orange','c','m','y','k']
  #l1=list(np.zeros((len(E),1)))
  E=np.array(E)
 
  for i in range(len(E[:,0])):
    E[i,:]=E[i,:]/A[i]

  fig=plt.figure()
  for i in range(len(E[:,0])):
  
    plt.plot(mu_X,E[i,:],linewidth=2.0,color=colors[i])
    legend.append('Structure :'+str(i+1))
  plt.legend(np.array(legend[:]),fontsize=font)
  plt.plot(mu_X,np.zeros((len(mu_X),1),dtype=float),linestyle='--',color='black')
  plt.axis('tight') 
  plt.xlabel('$\mu_{\\rm S}$(eV)',fontsize=font)
  plt.ylabel('$E_f$(eV/\\rm{\AA})',fontsize=font)
  plt.xticks(fontsize=font)
  plt.yticks(fontsize=font)
   
  plt.show()


def calculate_formation_TX2(Et,matt,nt,Eb,matb,nb,mu_M,mu_X,U):


      nnt=sum(nt)
      nnb=sum(nb)
      mult=int(nnt)/nnb
      Ed=np.zeros((len(mu_X),1),dtype=float)
      nt_M=0
      nt_T=0
      nt_C=0
      nb_M=0
      nb_T=0
      nb_C=0

      for i in range(len(matt)):

        if substance(matt[i])=='M':

           nt_M=nt[i]

        elif substance(matt[i])=='T':

           nt_T=nt[i]
        
        elif substance(matt[i])=='X':
 
           nt_X=nt[i] 

      for i in range(len(matb)):

        if substance(matb[i])=='M':

           nb_M=nb[i]

        elif substance(matb[i])=='T':

           nb_T=nb[i]

        elif substance(matb[i])=='X':

           nb_X=nb[i]

 
     
      del_X=-(nb_X*mult-nt_X)
      del_T=-(nb_T*mult-nt_T)
      del_M=-(nb_M*mult-nt_M)

#      print([del_X,del_T,del_M])
      print('Energy difference:'+str(Et-mult*Eb+U*del_T))

      for i in range(len(mu_X)):

         Ed[i]=(Et-mult*Eb-del_T*(Eb)-(del_X-2*del_T)*mu_X[i]-del_M*mu_M-U*del_M)/np.absolute(del_M)
#         Ed[i]=Et-mult*Eb-del_M*mu_X[i]-U*del_M

 #        print(-(del_X-2*del_T)*mu_X[i]-del_M*mu_M+U*del_T)        
      return Ed 



def main():

     from globvar import MAT,mus
     parser = argparse.ArgumentParser()

     parser.add_argument('strings', type=str, nargs='+')
     parser.add_argument("--mu", "-d", type=float, required=True)
     parser.add_argument("--nfile", "-dt",type=int,   required=True)
     parser.add_argument("-typ", "-t",type=str,   required=True)
     parser.add_argument("-U", "-U",type=float,   required=False)
     parser.add_argument("--font", "-font",type=int,   required=False)
     parser.add_argument('--direction', nargs='*')
     parser.add_argument("--mat",type=str, required=True)     
     args = parser.parse_args()
#     print args
     mat=args.mat
     nfiles=args.nfile
     natoms=[]
     materials=[]
     A=[]
     E=np.zeros((nfiles,1),dtype=float)


     if args.direction!= True:
       direction=np.array([0,1])
     else:
       direction=args.direction


     print('lattice plain :'+str(direction[0])+','+str(direction[1]))
 
     for i in range(nfiles):

       materials2,natoms2,E[i],A1=read_files(args.strings[i*2:i*2+2],direction)
       materials.append(materials2)
       natoms.append(natoms2)
       A.append(A1)

     A=np.array(A)
     materials=np.array(materials)

#  Bulk material attributes

     Eb=E[-1]
     nb=natoms[-1]
     matb=materials[-1]     
     for i in range(len(MAT)):

       if(mat==MAT[i]):

        mu_M=mus[i]
        print('Material chosen : '+str(MAT[i]))
        print('Chemical potential of '+str(MAT[i])+' : '+str(mu_M)+' eV') 
        break

 #    mu_X=np.linspace(-5.69,-3.29,10)     #(-21.83+10.45)/2<mu_X(Selenium)<-3.29(energy of a bare Se atom from A flow)
     mu_X=np.linspace(-5.59,-3.29,10)

     Et=[]
     
     for i in range(nfiles-1):

       if args.typ=='TX2':


           E2=calculate_formation_TX2(E[i],materials[i],natoms[i],Eb,matb,nb,mu_M,mu_X,args.U)         
     
       Et.append(E2)
     if args.font!=None:
       plot_line(Et,mu_X,A)
     else: 
       plot_line(Et,mu_X,A,args.font)

if __name__=='__main__':


# Parse bulk structure in the end else this code will calculate bullshitt! 
# Its a bad way of coding but who cares!!!!!!!!
# Default material is WSe_2

     main()







