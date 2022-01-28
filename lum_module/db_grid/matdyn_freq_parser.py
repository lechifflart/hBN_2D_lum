import numpy as np
#
#   TODO : case where number of frequencies is not a multiple of 6
#               |_ then the length of the lists are not the same and creates trouble when using .flatten() or ravel()
#

def matdyn_freq_parser(file):
    # get a the q point and corresponding frequencies from QE matdyn, in cm-1
    with open(file,'r') as f:
        lines=f.readlines()

    # extract nb of phonons modes and nb of k pts
    nb_modes_k = [int(s) for s in lines[0].replace(',',' ').split() if s.isdigit()]

    # loop for every point. Frequencies are 6 by line
    q_tilde = np.empty([nb_modes_k[1],3])
    frequencies = np.empty([nb_modes_k[1],nb_modes_k[0]])

    step = -(nb_modes_k[0]//-6)  # ceiling division ; frequencies are at most 6 by line

    for i in range(1,len(lines)-step,step+1):
        q_index = i//(step+1)
        q_tilde[q_index]=np.array([float(x) for x in lines[i].split()])
        frequencies[q_index] = np.array([float(x) for x in np.concatenate(np.array([line.split() for line in lines[i+1:i+step+1]],dtype=object))])

    # check if all points were parsed
    if len(q_tilde)==nb_modes_k[1]:
        return(q_tilde,frequencies)
    else :
        raise ValueError('[ERROR] wrong number of parsed kpoints in matdyn output')

#matdyn_freq_parser('matdyn_dbg.freq')
#matdyn_freq_parser('test_other_matdyn_output')
#print(matdyn_freq_parser('another_matdyn_output_test'))
