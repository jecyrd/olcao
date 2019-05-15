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

class AtomTypeType():
    def __init__(self):
        maxOrbitals = 4
        self.typeLabel = ''
        # Identity information
        self.elementName = ''
        self.elementID = 0
        self.speciesID = 0
        self.typeID = 0
        # Basis Wave function information
        self.numCoreRadialFns = 0
        self.numValeRadialFns = 0
        self.coreRadialFns = [[]]
        self.valeRadialFns = [[]]
        # Alpha information
        self.numOrbAlphas = []
        self.alphas = []
        # Radial function information
        self.maxCoreQN_l=0
        self.maxValeQN_l=0
        self.coreQN_nList = []
        self.valeQN_nList = []
        self.coreQN_lList = []
        self.valeQN_lList = []
        self.coreQN_2jList = []
        self.valeQN_2jList = []
        self.numQN_lCoreRadialFns = [0]*maxOrbitals
        self.numQN_lValeRadialFns = [0]*maxOrbitals
        # State information
        self.numCoreStates = 0
        self.numValeStates = 0
        self.numTotalStates = 0

def getAtomicTypeImplicitInfo(atomTypes):
    minAtomicAlpha      = 100
    maxNumAtomAlphas    = 0
    maxNumCoreRadialFns = 0
    maxNumValeRadialFns = 0 
    maxNumCoreStates    = 0
    maxNumValeStates    = 0
    maxNumStates        = 0
    numElements         = 0

    numAtomTypes = len(atomTypes)

    for i in range(numAtomTypes):
        # Obtain the element name for each atomic type
        atomTypes[i].elementName = re.findall('\b[A-Za-z]{1,2}',\
                atomTypes[i].typeLabel)

        # Compare the smallest alpha for this type to the smallest alpha
        # seen yet
        minAtomicAlpha = min(minAtomicAlpha,atomTypes[i].alphas[0])

        # Compare the number of alphas for this type to the largest
        # number of alphas seen for any atomic type yet.
        if (atomTypes[i].numOrbAlphas[0] > maxNumAtomAlphas):
            maxNumAtomAlphas = atomTypes[i].numOrbAlphas[0]

        # Track the number of s,p,d,f radial functions for the core and
        # valence of this type
        atomTypes[i].maxCorenQN_l, \
        atomTypes[i].numCoreStates, \
        atomTypes[i].numCoreRadialFns = \
            countStatesAndFns(atomTypes[i].numQN_lCoreRadialFns,\
                    atomTypes[i].numCoreRadialFns,\
                    atomTypes[i].coreQN_lList)
        atomTypes[i].maxValenQN_l, \
        atomTypes[i].numValeStates, \
        atomTypes[i].numValeRadialFns = \
            countStatesAndFns(atomTypes[i].numQN_lValeRadialFns,\
                    atomTypes[i].numValeRadialFns,\
                    atomTypes[i].valeQN_lList)

        # Sum the results of the counters for the number of valence and
        # core states into one value
        atomTypes[i].numTotalStates = atomTypes[i].numValeStates + \
                atomTypes[i].numCoreStates

        # Compare results for this atomic type to the maximum values yet
        # seen.
        maxNumCoreStates = max(maxNumCoreStates,atomTypes[i].numCoreStates)
        maxNumValeStates = max(maxNumValeStates,atomTypes[i].numValeStates)
        maxNumStates = max(maxNumStates,atomTypes[i].numTotalStates)

        # Compare the number of core radial functions for this type to the
        # maximum value yet seen.
        if (atomTypes[i].numCoreRadialFns > maxNumCoreRadialFns):
            maxNumCoreRadialFns = atomTypes[i].numCoreRadialFns

        # Compare the number of valence radial functions for this type to
        # the maximum value yet seen
        if (atomTypes[i].numValeRadialFns > maxNumValeRadialFns):
            maxNumValeRadialFns = atomTypes[i].numValeRadialFns

        # Compare the element ID number for this type the the highest yet
        # seen.
        if (atomTypes[i].elementID > numElements):
            numElements = atomTypes[i].elementID

    return atomTypes
 
def countStatesAndFns(numQN_lRadialFns, numRadialFns, QN_lList):

    numQN_lRadialFns[:] = [0]*len(numQN_lRadialFns)

    maxQN_l = -1

    numStates = 0

    # Loop over the radial functions for this given type (core or valence).
    for i in range(numRadialFns):
        # Increment the valence s,p,d,f function counter
        numQN_lRadialFns[QN_lList[i]] = numQN_lRadialFns[QN_lList[i]] + 1

        # Increment the total state counter
        numStates = numStates + 2*QN_lList[i]+1

        # Determine if this is hte highest angular momentum value yet seen
        # for this atomic type (core or valence).
        if (QN_lList[i] > maxQN_l):
            maxQN_l = QN_lList[i]

    return maxQN_l, numStates, numRadialFns
        

def readAtomicTypes(oDatFile, basisCode):
    maxOrbitals = 4
    f = open(oDatFile,'r')
    olcaodat = [line.strip() for line in f.readlines()]
    f.close()

    # Find the
    line=0
    while (not ('NUM_ATOM_TYPES' in olcaodat[line])):
        line+=1
    num_atom_types = int(olcaodat[line+1])

    atomTypes = []
    line += 2
    for i in range(num_atom_types):
        atomTypes.append(AtomTypeType())

        # Read down until we find the next "ATOM_TYPE_ID__SEQUENTIAL_NUMBER"
        while not ('ATOM_TYPE_ID__SEQUENTIAL_NUMBER' in olcaodat[line]):
            line+=1
        line +=1
        # Line after ATOM_TYPE_ID__SEQUENTIAL_NUMBER 
        typeIDs = olcaodat[line].split()
        atomTypes[i].elementID = int(typeIDs[0])
        atomTypes[i].speciesID = int(typeIDs[1])
        atomTypes[i].typeID    = int(typeIDs[2])
        # The last value of typeIDs equals loop index i

        # Read the type label for the atom type and adjust the spacing
        line += 2
        atomTypes[i].typeLabel = olcaodat[line]

        # Read the number of gaussian terms for the s,p,d,f radial basis
        # functions. They are indexed as 'alphas' from exp(-alpha*r^2).
        line += 2
        atomTypes[i].numOrbAlphas = list(map(int,olcaodat[line].split()))

        # Check that hte number of alphas is monotonically decreasing
        for j in range(maxOrbitals-1):
            if (atomTypes[i].numOrbAlphas[j] < \
                    atomTypes[i].numOrbAlphas[j+1]):
                print('Num of alphas not monotonically decreasing.')

        # Read the alphas for this atomic type.
        line += 2
        nlines = int(m.ceil(atomTypes[i].numOrbAlphas[0]/4.0))
        for j in range(nlines):
            atomTypes[i].alphas+=list(map(float,olcaodat[line].split()))
            line+=1
        
        # Read in the number of core radial basis functions for this atomic
        # type.
        line += 1
        tempIntArray = list(map(int,olcaodat[line].split()))
        atomTypes[i].numCoreRadialFns = tempIntArray[basisCode]

        # Allocate space
        atomTypes[i].coreQN_nList = [0]*atomTypes[i].numCoreRadialFns
        atomTypes[i].coreQN_lList = [0]*atomTypes[i].numCoreRadialFns
        atomTypes[i].coreQN_2jList = [0]*atomTypes[i].numCoreRadialFns

        atomTypes[i].coreRadialFns = [[0]*atomTypes[i].numOrbAlphas[0]]*\
                atomTypes[i].numCoreRadialFns


        # Read in the core radial basis functions. We pass this off to a
        # subroutine because the procedure for reading in the functions
        # for the core and valence is the same. They just record the data
        # to different data structions
        line += 2
        readRadialFns(olcaodat, line, tempIntArray[2], basisCode,\
                atomTypes[i].numOrbAlphas,\
                atomTypes[i].coreQN_nList,\
                atomTypes[i].coreQN_lList,\
                atomTypes[i].coreQN_2jList,\
                atomTypes[i].coreRadialFns)
        
        # Read forward to NUM_VALE_RADIAL_FNS
        while not ('NUM_VALE_RADIAL_FNS' in olcaodat[line]):
            line+=1
        line +=1
        # Read in the num of valence radial basis functions for this atomic
        # type.
        tempIntArray = list(map(int,olcaodat[line].split()))
        atomTypes[i].numValeRadialFns = tempIntArray[basisCode]
        
        # Allocate space
        atomTypes[i].valeQN_nList = [0]*atomTypes[i].numValeRadialFns
        atomTypes[i].valeQN_lList = [0]*atomTypes[i].numValeRadialFns
        atomTypes[i].valeQN_2jList = [0]*atomTypes[i].numValeRadialFns

        atomTypes[i].valeRadialFns = [[0]*atomTypes[i].numOrbAlphas[0]]*\
                atomTypes[i].numValeRadialFns
        
        # Read in the vale radial basis functions. 
        line += 2
        readRadialFns(olcaodat, line, tempIntArray[2], basisCode, \
                atomTypes[i].numOrbAlphas,\
                atomTypes[i].valeQN_nList,\
                atomTypes[i].valeQN_lList,\
                atomTypes[i].valeQN_2jList,\
                atomTypes[i].valeRadialFns)

    return atomTypes


def readRadialFns(olcaodat, line, numRadialFns, basisCode, numOrbAlphas,\
        QN_nList, QN_lList, QN_2jList, radialFns):
    #readunit, writeunit, tempintarray(3), 
    #atomtypes[i].numOrbAlphas, atomTypes[i].coreQN_nList, 
    #atomTypes[i].coreQN_lList, atomTypes[i].coreQN_2jList
    #atomTypes[i].coreRadialFns)

    i=-1
    for j in range(numRadialFns):
        tempIntArray = [0]*5

        tempIntArray[:2] = list(map(int,olcaodat[line].split()))[0:]
        if (tempIntArray[1] <= basisCode+1):
            i+=1
            line+=1
            tempIntArray = list(map(int,olcaodat[line].split()))
            QN_nList[i] = tempIntArray[0]
            QN_lList[i] = tempIntArray[1]
            QN_2jList[i] = tempIntArray[2]

            #line+=1
            nlines = int(m.ceil(numOrbAlphas[QN_lList[i]]/4.0))
            c = 0
            for k in range(nlines):
                line+=1
                for l in list(map(float,olcaodat[line].split())):
                    if (c<numOrbAlphas[QN_lList[i]]):
                        radialFns[i][c] = l
                    c+=1
            line+=1
        else:
            # Should be at the line with two integers representing component
            # and basis. If we move 2 lines fowardward we'll be at the line
            # where the radial functions begin
            line+=2
            # Now we need to move ahead the numbers of lines with radial
            # functions on them
            line += int(m.ceil(numOrbAlphas[QN_lList[i]]/4.0))
            # Now we move ahead one more line to be at the start of the
            # next block
            #line +=1


if __name__=="__main__":
    pass
