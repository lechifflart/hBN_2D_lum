import numpy as np


iz = 0.000000
# example : grid is 9x9x1
grid_size = 9
int_half = grid_size // 2
xx = np.array([x for x in range(-int_half,int_half+1)])/9.
with open('db_grid.txt','w') as f :
    for ix in xx:
        for iy in xx:
            f.write(f" {ix:.6f} {iy:.6f} {iz:.6f} \n")
