#!/usr/bin/python

import math as m


def gridCalc(numprocs):
  npcol = 2*int(m.sqrt(numprocs))
  for x in range(npcol,0,-1):
    if (numprocs % x == 0):
      npcol = x
      break
  
  nprow = numprocs / npcol
  if (nprow > npcol):
    x = npcol
    npcol = nprow
    nprow = x

  return nprow, npcol


for x in xrange(1023):
  #nprow, npcol = gridCalc(x+2)
  #print "Numprocs: ", x+2, "nprow: ", nprow, "npcol: ", npcol
