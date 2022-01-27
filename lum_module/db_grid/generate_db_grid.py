import numpy as np

xx = np.linspace(-0.5,0.5,9)
yy = np.linspace(-0.5,0.5,9)
iz = 0.000000
with open('db_grid.txt','w') as f :
    for ix in xx:
        for iy in yy:
            f.write(f" {ix:.6f} {iy:.6f} {iz:.6f} \n")
