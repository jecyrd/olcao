#!/usr/bin/python

import math as m


def gridCalc(numprocs):
  prow = int(m.sqrt(numprocs))
  #pcol = int(m.sqrt(numprocs))

  #if (numprocs-(prow*pcol) >= prow):
  #  prow+=1
  pcol = numprocs/prow

  return prow,pcol


for x in xrange(2,10):
  nprow, npcol = gridCalc(x)
  print "Numprocs: ", x, 'Extra:',(x)-nprow*npcol, "nprow: ", nprow, "npcol: ", npcol

#  print "\t\tStuff: ", x/npcol, (x-1)/npcol
