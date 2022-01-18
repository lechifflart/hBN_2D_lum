import numpy as np

def ypp_output_parser(file):
    with open(file) as f:
        data = np.genfromtxt(fname=file,comments='#')
