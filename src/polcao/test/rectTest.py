#!/user/bin/python




def checkOverlap(lo1, hi1, lo2, hi2):
  if ( (lo1[1]<=hi2[1]) and (hi1[1]>=lo2[1]) and \
       (-lo1[0]>=-hi2[0]) and (-hi1[0]<=-lo2[0]) ):
    return 1

  return 0

def getOverlap(lo1, hi1, lo2, hi2):
  loOvlp = [0,0]
  hiOvlp = [0,0]

  loOvlp[0] = max(lo1[0], lo2[0])
  loOvlp[1] = max(lo1[1], lo2[1])
  hiOvlp[0] = min(hi1[0], hi2[0])
  hiOvlp[1] = min(hi1[1], hi2[1])

  return loOvlp, hiOvlp

def inIntersection(i, j, loOvlp, hiOvlp):
  if ( ((i>=loOvlp[0]) and (i<=hiOvlp[0])) and\
       ((j>=loOvlp[1]) and (j<=hiOvlp[1])) ):
       return 1

  return 0

def main():
  lo1 = [1,1]
  hi1 = [4,4]
  
  lo2 = [3,3]
  hi2 = [6,6]
  
  lo3 = [4,2]
  hi3 = [6,3]

  indices = []
  for i in xrange(lo1[0], hi1[0]+1):
    for j in xrange(lo1[1], hi1[1]+1):
      indices.append([i,j])

  indices.append([0,0,0])

  if (checkOverlap(lo1,hi1,lo2,hi2)):
    loOvlp, hiOvlp = getOverlap(lo1,hi1,lo2,hi2)
    print loOvlp,hiOvlp

    for i in xrange(lo2[0], hi2[0]+1):
      for j in xrange(lo2[1], hi2[1]+1):
        if (not inIntersection(i,j, loOvlp, hiOvlp)):
          indices.append([i,j])
  indices.append([0,0,0])

  if (checkOverlap(lo1,hi1,lo3,hi3)):
    print 'true'
    loOvlp1, hiOvlp1 = getOverlap(lo1,hi1,lo3,hi3)

  if (checkOverlap(lo2,hi2,lo3,hi3)):
    print 'true'
    loOvlp2, hiOvlp2 = getOverlap(lo2,hi2,lo3,hi3)
    
    for i in xrange(lo3[0], hi3[0]+1):
      for j in xrange(lo3[1], hi3[1]+1):
        if ((not inIntersection(i,j, loOvlp1, hiOvlp1)) and
            (not inIntersection(i,j, loOvlp2, hiOvlp2)) ):
          indices.append([i,j])
  
  for x in xrange(len(indices)):  
    print indices[x]

if __name__=="__main__":
  main()
