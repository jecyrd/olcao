#!/usr/bin/env python3
#
# This module calculates the the distribution information for each if the 
# desired processors during the Fortran execution
#
# BIG NOTE: Remember that as is in the SCALAPACK documentation at:
# https://www.netlib.org/scalapack/slug/node76.html Global Matrix and 
# Local Matrix Indices start with 1. However, the process grid, and 
# local block coordinate, start with indices start with 0.

import os
import sys
import h5py
import math as m
import numpy as np
import collections

import atomicTypes as aType
import atomicSites as aSite
from parDescs import procDescriptor

def createAtomicStructs(fort4, fort5, verbose=False):
    atomTypes = aType.readAtomicTypes(fort5,1)

    atomTypes = aType.getAtomicTypeImplicitInfo(atomTypes)

    atomSites = aSite.readAtomicSites(fort4)

    coreDim, valeDim = aSite.getAtomicSitesImplicitInfo(atomSites, atomTypes)

    if verbose:
        print('types')
        for atype in atomTypes:
            print(atype.numCoreStates,atype.numValeStates,atype.numTotalStates)
        print('total', coreDim, valeDim)

        print('sites')
        print(len(atomSites))
        for i in range(len(atomSites)):
            print(atomSites[i].cumulCoreStates, atomSites[i].cumulValeStates)


    return atomTypes, atomSites, coreDim, valeDim

def createDescriptors(coreDim,valeDim,nprocs=2):
    ccDescs = [procDescriptor(i,nprocs,coreDim,coreDim) for i in range(nprocs)]
    cvDescs = [procDescriptor(i,nprocs,coreDim,valeDim) for i in range(nprocs)]
    vvDescs = [procDescriptor(i,nprocs,valeDim,valeDim) for i in range(nprocs)]

    return ccDescs,cvDescs,vvDescs

def rank2grid(rank,pgrid=(2,2)):
    prow = int(m.floor(rank/pgrid[1]))
    pcol = int(rank%pgrid[1])

    return (prow,pcol)

def grid2rank(myrow,mycol,pgrid):
    return (myrow*pgrid[1])+mycol

def calcProc(i,j,mb,nb, pgrid):
    pr = m.floor((i-1)/mb) % pgrid[0]
    pc = m.floor((j-1)/nb) % pgrid[1]
    return (pr, pc)

def atomPairs(atomSites, atomTypes, descs, dim1, dim2):
    atomDict = { i:[] for i in range(descs[0].mpisize) }
    mb = descs[0].mb
    nb = descs[0].nb
    pgrid = (descs[0].prows,descs[0].pcols)
    for i in range(len(atomSites)):
        for j in range(len(atomSites)):
            if dim1 is 'c':
                glor = atomSites[i].cumulCoreStates+1
                ghir = glor+atomTypes[atomSites[i].atomTypeAssn-1].numCoreStates
                ghir -= 1
            elif dim1 is 'v':
                glor = atomSites[i].cumulValeStates+1
                ghir = glor+atomTypes[atomSites[i].atomTypeAssn-1].numValeStates
                ghir -= 1

            if dim2 is 'c':
                gloc = atomSites[j].cumulCoreStates+1
                ghic = gloc+atomTypes[atomSites[j].atomTypeAssn-1].numCoreStates
                ghic -= 1
            if dim2 is 'v':
                gloc = atomSites[j].cumulValeStates+1
                ghic = gloc+atomTypes[atomSites[j].atomTypeAssn-1].numValeStates
                ghic -= 1
            
            lopr,lopc = calcProc(glor,gloc,mb,nb,pgrid)
            hipr,hipc = calcProc(ghir,ghic,mb,nb,pgrid)
            
            if (lopr <= hipr):
                rows = list(range(lopr,hipr+1))
            else:
                rows = [hipr]+list(range(lopr,pgrid[0]))

            if (lopc <= hipc):
                cols = list(range(lopc,hipc+1))
            else:
                cols = [hipc] + list(range(lopc,pgrid[1]))


            for k in rows:
                for l in cols:
                    rank = grid2rank(k,l,pgrid)
                    atomDict[rank].append((i,j))
    return atomDict

def getAtomPairs(atomSites, atomTypes, ccDescs, cvDescs, vvDescs,verbose=False):
    ccAtomDict = atomPairs(atomSites,atomTypes,ccDescs,'c','c')
    cvAtomDict = atomPairs(atomSites,atomTypes,cvDescs,'c','v')
    vvAtomDict = atomPairs(atomSites,atomTypes,vvDescs,'v','v')

    if verbose:
        for key in sorted(vvAtomDict.keys()):
            print(f'Proc {key} has pairs: {sorted(vvAtomDict[key])}')

    atomDict = { i:[] for i in range(len(ccDescs)) }
    for i in range(ccDescs[0].mpisize):
        atomDict[i] = set(ccAtomDict[i]) | \
                      set(cvAtomDict[i]) | \
                      set(vvAtomDict[i])

    return atomDict

def mostSquareGridAllProcs(nprocs):
    if (m.sqrt(nprocs)%1 != 0):
        x = m.floor(m.sqrt(nprocs))
        while (nprocs%x!=0 and x>=1):
            x-=1
        y = m.floor(nprocs/x)
    else:
        x = int(m.sqrt(nprocs))
        y = x
    return (x,y)

def mostSquareGridLessProcs(nprocs):
    x = int(m.sqrt(nprocs))
    y = int(m.floor(nprocs/x))
    return (x,y)

def checkGrid(nprocs, verbose=False):
    # This number is arbitrary, need performance analysis on the effect.
    maxGridDiff = 2 
    nnew = nprocs
    agrid = mostSquareGridAllProcs(nprocs)
    if verbose:
        print(f'The most square grid using all available procs is: {agrid}')

    if (abs(agrid[0]-agrid[1])>maxGridDiff):
        lgrid = mostSquareGridLessProcs(nprocs)
    
        nnew = lgrid[0]*lgrid[1]
        if verbose:
            print(f'Instead using a grid of {lgrid} with {nnew} proccesses')

    return nnew

def createParallelDistribution(nprocs=2,verbose=False):
    fort4 = './inputs/structure.dat'
    fort5 = './inputs/olcao.dat'

    types,sites,coreDim,valeDim = \
            createAtomicStructs(fort4,fort5,verbose)

    ccDescs,cvDescs,vvDescs = createDescriptors(coreDim,valeDim,nprocs)

    if verbose:
        for vv in vvDescs:
            ostr = f"rank:{vv.mpirank}, pgrid: ({vv.prows},{vv.pcols}),"\
                    f" gridloc: ({vv.myprow},{vv.mypcol})"+\
                    f" bFactor: ({vv.mb},{vv.nb}))"
            print(ostr)

    atomPairs = getAtomPairs(sites, types, ccDescs, cvDescs, vvDescs)

    if verbose:
        for key in atomPairs.keys():
            print(f'Rank {key} has pairs: {sorted(atomPairs[key])}')

    pgrid = (vvDescs[0].prows,vvDescs[1].pcols)
    createProcInputFile(atomPairs,pgrid)

def createProcInputFile(atomPairs,pgrid=(2,2)):
    hf = h5py.File('polcaoDistribution.hdf5')

    hf.create_dataset('nprocs',data=len(atomPairs.keys()))
    hf.create_dataset('pgrid',data=np.array(pgrid))
    gp = hf.create_group('ProcessAtomList')
    glen = gp.create_group('AtomLen')
    galist = gp.create_group('AtomList')
    for key in sorted(atomPairs.keys()):
        #print(key,len(atomPairs[key]),type(atomPairs[key]))
        ds = np.array(sorted(atomPairs[key]))
        ds = np.add(ds,1)
        glen.create_dataset(f'{key}',data=np.array(len(atomPairs[key])))
        galist.create_dataset(f'{key}',data=ds) #,\
                #compression="gzip",compression_opts=9)
       
    hf.close()

if __name__=="__main__":
    usage = f"usage: {sys.argv[0]} <nprocs>"

    fort4 = './inputs/structure.dat'
    fort5 = './inputs/olcao.dat'

    try:
        nprocs = int(sys.argv[1])
        f=open(fort4,'r')
        del f
        f=open(fort5,'r')
        del f

    except ValueError:
        print("First arguement not a valid number of processes.")
        print(usage)
        sys.exit(1)
    except IndexError:
        print("Number of processes not specified.")
        print(usage)
        sys.exit(1)
    except FileNotFoundError:
        print(f"Could not find {fort4}, or {fort5} for reading")
        sys.exit(1)

    nprocs = checkGrid(nprocs)
    createParallelDistribution(nprocs)
