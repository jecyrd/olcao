#!/usr/bin/env python3
#
# This module replicates work done in the fortran code to make information 
# available about system to the routines that decide the best parallel
# layout of the interaction matrices.
#
# Author: James Currie <jecyrd@mail.umkc.edu>
#

import re
import math as m
import sys
import copy

class AtomSiteType():
    def __init__(self):
        self.atomTypeAssn = 0
        self.cartPos = [0.0]*3
        self.cumulCoreStates = 0
        self.cumulValeStates = 0

def getAtomicSitesImplicitInfo(atomSites, atomTypes):

    numAtomSites = len(atomSites)

    # Initialize the core and valence dimensions defined by the number of
    # core or valence states for each atomic type multiplied by the number
    # of atoms for that type.
    valeDim = 0
    coreDim = 0
    maxDim  = 0

    # Loop over the number of atoms in the system. At each interation:
    # 1) Increment the coreDim and valeDim by the number of states for the
    #    current atom.
    # 2) Store that value with the atom's data structure.
    for i in range(numAtomSites):

        # Get the current type of this atom
        currentType = atomSites[i].atomTypeAssn -1

        # Store the current summation value with the current atom's data.
        # Note that the cumulative sum is off by one atom so that atom
        # 1 has a sum of zero.
        atomSites[i].cumulValeStates = valeDim
        atomSites[i].cumulCoreStates = coreDim

        # Cumulatively add the valence and core dimensions of each atom.
        valeDim += atomTypes[currentType].numValeStates
        coreDim += atomTypes[currentType].numCoreStates

    maxdim = max(valeDim,coreDim)

    return coreDim, valeDim


def readAtomicSites(sDatFile):
    f = open(sDatFile,'r')
    structdat = [line.strip() for line in f.readlines()]
    f.close()

    # Read the number of atomic sites
    line=0
    while not ('NUM_ATOM_SITES' in structdat[line]):
        line+=1
    line+=1
    numAtomSites = int(structdat[line])

    atomSites = [AtomSiteType() for i in range(numAtomSites)]

    line +=2
    for i in range(numAtomSites):
        splist = structdat[line].split()
        counter = int(splist[0])
        atomSites[i].atomTypeAssn = int(splist[1])
        atomSites[i].cartPos[0] = float(splist[2])
        atomSites[i].cartPos[1] = float(splist[3])
        atomSites[i].cartPos[2] = float(splist[4])
        atomName = splist[5]

        line+=1
        if (counter != i+1):
            print('Potential site list is out of order at i = ',i)
            sys.exit()

    return atomSites

if __name__=="__main__":
    pass
