import os
import sys
import argparse
import logging
import numpy as np
import subprocess

from lsdmap.rw import reader


class SelectionStep(object):
    """
    SelectionStep()

    A class used to perform the selection step of Extended DM-d-MD 
    """

    def create_arg_parser(self):

        parser = argparse.ArgumentParser(description="select n new configurations uniformily along the first and second DCs..")
        parser.add_argument('npoints', metavar='n', type=int, help='number of configurations to be selected')

        # required options
        parser.add_argument("-s",
           type=str,
           dest="kernel",
           required=True,
           help="File containing the whole kernel from lsdmap: .kernel")

        parser.add_argument("-c",
           type=str,
           dest="inputfile",
           required=True,
           help="File containing  old structures .gro")

        parser.add_argument("-w",
           type=str,
           dest="wfile",
           required=True,
           help="File containing the old weights .w")

        parser.add_argument("-o",
           type=str,
           dest="outputfile",
           required=True,
           help="where save new structures .gro")

        parser.add_argument("-nc",
           type=str,
           dest="ncfile",
           required=False,
           help='File containing a single column with the number of copies to be made for each configuration (output, opt.): nc')
        
        parser.add_argument("-r",
           type=str,
           dest="startfile",
           required=True,
           help="for mdanalysis reference structure .gro")

        return parser

   
    def run(self):

        parser = self.create_arg_parser()
        args = parser.parse_args()

        kernel = np.loadtxt(args.kernel)
        size=kernel.shape[0]
        
        max_v=np.zeros((size,size))
        max_ind=np.zeros((size,size))

        for n in range(1,10):
          #print "n", n
          for i1 in range(size):
           for i2 in range(size): 
            if max_v[i1,i2]<kernel[i1,i2]:
             max_v[i1,i2]=kernel[i1,i2]
             max_ind[i1,i2]=n
          new=np.dot(kernel,kernel)
          kernel=np.copy(new)
        
        #np.save('max_ind.npy',max_ind)

        test=np.zeros((size,size))
        test[max_ind<4]=1
        for i in range(size):
        	test[i,i]=0
        ncopies=np.full((size),1)


        for i in range(1000):
          array=np.dot(test,ncopies)
          #print ncopies.sum(),ncopies.min(),ncopies.max(), array[np.argmax(array)], array[np.argmin(array)]
          choice_max=np.random.choice(np.where(test[np.argmax(array),:]==1)[0],100)
          n=100
          for i2 in range(100):
            choice=choice_max[i2]
            if ncopies[choice]>0:
              ncopies[choice]=ncopies[choice]-1
            else:
              n=n-1
          choice_min=np.random.choice(np.where(test[np.argmin(array),:]==1)[0],100)
          for i3 in range(n):
              ncopies[choice_min[i3]]=ncopies[choice_min[i3]]+1        


        np.savetxt(args.ncfile,ncopies,fmt='%0i')
       
        #convert inputfile to xtc so later can write fast
        p2=subprocess.call("aprun -n 1 -N 1 -d 1 trjconv_mpi -f "+str(args.inputfile)+" -o tmp.xtc 1>/dev/null 2>/dev/null",shell=True)
        #read in structures
        import MDAnalysis
        u = MDAnalysis.Universe(startfile,'tmp.xtc')
      
        #continue to calculate new weights, which structures
        weights=np.loadtxt(wfile)

        test=np.zeros((size,size))
        test[max_ind<2]=1
        for i in range(size):
                test[i,i]=0
        #first distribute weights from killed replicas
        for i in range(size):
          if int(ncopies[i])==0:
            choice1=np.where(test[i,:]==1)[0]
            choice2=np.random.choice(np.where(ncopies[choice1[:]]>0)[0],100)
            #choice1[choice2[:]]
            for j in range(choice1[choice2[:]].shape[0]):
              j2=choice1[choice2[j]]
              weights[j2]=weights[j2]+weights[i]/100.0
            weights[i]=0
        #then distribute weight along duplicated weights, generate structure file
        #and save new structures 
        with MDAnalysis.coordinates.GRO.GROWriter(outputfile, u.atoms.n_atoms) as w:
         weights_new=np.zeros(int(ncopies.sum()))
         old_index=np.zeros(int(ncopies.sum()))
         k=0
         for i in range(size):
           if int(ncopies[i])>0:
             for j in range(int(ncopies[i])):
               weights_new[k]=weights[i]/ncopies[i]
               old_index[k]=i
               k=k+1
               w.write(u, frame=i)
 

        #save new weights
        np.savetxt(wfile,weights_new,fmt='%0.15f')
 

if __name__ == "__main__":
    SelectionStep().run()
