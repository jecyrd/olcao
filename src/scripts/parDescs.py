#!/usr/bin/env python3
#
# This module handles calculations of the blacs parameters and the 
# local array parameters for the parallel distribution.
#
# Author: James Currie <jecyrd@mail.umkc.edu>
#

import math as m

def numroc(n,nb,iproc,isrcproc,nprocs):
    ''' 
    n - The number of rows/columns in the distributed matrix
    nb - Block size, size of the blocks the distributed matrix is split into
    iproc - The coordinate of the process whose local array row or column
            is to be determined.
    isrcproc - The coordinate of the process that possesses the first row
               or column of the distributed matrix.
    nprocs - The total number of processes over which the matrix is 
             distributed.
    '''

    # Figure proc's distance from source process
    mydist = (nprocs+iproc-isrcproc) % nprocs

    # Figure the total number of whole NB blocks N is split up into
    nblocks = n // nb

    # Figure the minimum number of rows/cols a process can have
    num = (nblocks // nprocs) * nb

    # See if there are any extra blocks
    extrablks = nblocks % nprocs

    # If I have an extra block
    if (mydist < extrablks):
        num += nb
    # If I have last block, it may be a partial block
    elif (mydist == extrablks):
        num += n%nb

    return num

def localCoordToGlobal(l,mb,Pr,pr,x):
    '''
    (a,b) are the coordinates with the local array.
    (l,m) are the coordinates of a block within the local array .
    (x,y) are the coordinates within a local block in the local array.
    (I,J) are the coordinates in the global array.
    (Pr,Pc) are the dimensions of the process grid.
    (pr,pc) are the coordinates of a process in the process grid.

    Per the netlib documentation the global coordinates can be caluclated by:
    I = (l*Pr + pr)+mb+x
    J = (m*Pc + pc)+nb+y
    '''

    return (l*Pr + pr) + mb + x

def localToGlobal(l,m,x,y,desc):
    '''
    See notation notes above in localCoordToGlobal
    '''

    I = localCoordToGlobal(l, desc.mb, desc.prows, desc.myprow, x)
    J = localCoordToGlobal(m, desc.nb, desc.pcols, desc.mypcol, y)

    return I,J

def calcProcGridLoc(I,J,mb,nb,Pr,Pc):
    pr = ((I-1) // nb) % Pr
    pc = ((J-1) // nb) % Pc

    return pr,pc

class procDescriptor(object):
    def __init__(self,mpirank,mpisize,I,J):
        self.mpirank = mpirank
        self.mpisize = mpisize
        self.I = I
        self.J = J

        self.prows = 0
        self.pcols = 0
        self.myprow = 0
        self.mypcol = 0

        self.mb = 0
        self.nb = 0

        self.nrblocks = 0
        self.ncblocks = 0
        self.extraRows = 0
        self.extraCols = 0

        self.ilocal = 0
        self.jlocal = 0

        self.setup()
    
    def calcProcGrid(self):
        if (self.mpisize > 1):
            self.prows = int(m.sqrt(self.mpisize))
            self.pcols = int(self.mpisize // self.prows)
        else:
            self.prows = 1
            self.pcols = 1

    def calcBlockingFactors(self):
        self.mb = int(self.I // self.prows)
        self.nb = int(self.J // self.pcols)
    
    def calcMyCoords(self):
        # There's a way to calculate this and not loop over all blocks.
        #cnt = 0
        #for i in range(self.prows):
        #    for j in range(self.pcols):
        #        if (cnt == self.mpirank):
        #            self.myprow = i
        #            self.mypcol = j
        #        cnt+=1
        self.myprow = self.mpirank // self.pcols
        self.mypcol = self.mpirank % self.pcols

    def calcOtherMetaData(self):
        self.ilocal = numroc(self.I, self.mb, self.myprow, 0, self.prows)
        self.jlocal = numroc(self.J, self.nb, self.mypcol, 0, self.pcols)

        extraRow = (self.I // self.mb) % self.prows
        extraCol = (self.J // self.nb) % self.pcols

        self.nrblocks = self.I // self.mb // self.prows

        if ( (extraRow > 0) and (extraRow == self.myprow) ):
            self.extraRows = self.I % self.mb
            self.nrblocks += 1
        else:
            self.extraRows = 0
        if ( (extraCol > 0) and (extraCol == self.mypcol) ):
            self.extraCols = self.J % self.nb
            self.ncblocks += 1
        else:
            self.extraCols = 0

    def setup(self):
        self.calcProcGrid()
        self.calcBlockingFactors()
        self.calcMyCoords()

        self.calcOtherMetaData()
