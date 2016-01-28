import os
import sys
import argparse
import logging
import numpy as np

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

        parser.add_argument("-o",
           type=str,
           dest="ncfile",
           required=False,
           help='File containing a single column with the number of copies to be made for each configuration (output, opt.): nc')

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
          for i in range(100):
            choice=choice_max[i]
            if ncopies[choice]>0:
              ncopies[choice]=ncopies[choice]-1
            else:
              n=n-1
          choice_min=np.random.choice(np.where(test[np.argmin(array),:]==1)[0],100)
          for i in range(n):
              ncopies[choice_min[i]]=ncopies[choice_min[i]]+1        


       np.savetxt(args.ncfile,ncopies,fmt='%0i')


if __name__ == "__main__":
    SelectionStep().run()
