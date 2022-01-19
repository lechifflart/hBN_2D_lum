import numpy as np

def matdyn_freq_parser(file):
    # get a the q point and corresponding frequencies from QE matdyn, in cm-1
    q_tilde = []
    frequencies = []
    with open(file,'r') as f:
        lines=f.readlines()
    for i in range(1,len(lines),2):
        q_tilde.append(np.array([float(x) for x in lines[i].split()]))
        frequencies.append(np.array([float(x) for x in lines[i+1].split()]))
    return(q_tilde,frequencies)

#print(matdyn_freq_parser('matdyn_dbg.freq'))
