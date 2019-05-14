#!/usr/bin/env python3

import math as m

for i in range(1,101):
    x = m.sqrt(i)%1
    if (x != 0):
        x = m.floor(m.sqrt(i))
        while (i%x!=0 and x>=1):
            x-=1
        y = m.floor(i/x)
    else:
        x = int(m.sqrt(i))
        y=x
    print(f'{i}, {x}, {y}')
        

