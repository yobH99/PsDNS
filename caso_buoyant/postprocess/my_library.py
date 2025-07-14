import numpy as np 
import matplotlib.pyplot as plt
import subprocess
import os
import scienceplots


def read_file_to_store_index(file):
    '''
    Details relevant to the code:
    In the output files of PSDNS, the data are saved as functions of the form phi(z, t).
    At each time t = time, the profile phi(z, time) is saved as a column.
    There are therefore two headers for each time dump that start with the character #:
        #t = 0
        #z, epsxx, ... (in general, all the names of the variables being saved).
    Purpose of the code:
    Returns:
    nt : Number of time steps saved.
    names_profiles : Names of the quantities being saved.
    dim : Should be the vertical dimension (nz).
    time_vec : Physical times at which data have been saved.
    '''
    time_vec = []
    names_profiles = []
    dim = 0

    with open(file, 'r') as file:
        for line in file:
            line = line.strip()

            # Check for the time variable header saved in the file
            if '=' in line:
                dim = 0
                t = line.split()[-1]
                time_vec.append(float(t))

            # Check for the names of variables saved in the file
            # Redundant but consistent since the output format is fixed
            elif '#' in line and '=' not in line:
                names_profiles = line.split()[1:]

            # Count rows for the vertical dimension nz
            elif '#' not in line and '=' not in line:
                data_row = list(map(float, line.split()))
                if len(data_row) == 0:
                    continue
                else:
                    dim += 1

    # Number of time steps saved
    nt = len(time_vec)

    return nt, names_profiles, dim, time_vec

def generate_data_structure(file,nt,dim,names_profiles):
    
    nvar = len(names_profiles)
    data = np.zeros((nt,nvar,dim))
    i = -1
    k = 0 

    with open(file, 'r') as file:
        for line in file:
            if '#' not in line and '=' not in line:
                data_row = list(map(float, line.split()))
                if(len(data_row)==0):
                    continue
                else:
                    for j in range(nvar):
                        data[i,j,k] = data_row[j]
                    k = k + 1 
            elif '=' in line:
                i = i + 1 
                k = 0 
            
    data = np.array(data)   
    return data


def read_file_to_store_index_pt2():
    '''
    Details relevant to the code:
    In the output files of PSDNS, the data are saved as functions of the form phi(z, t).
    At each time t = time, the profile phi(z, time) is saved as a column.
    There are therefore two headers for each time dump that start with the character #:
        #t = 0
        #z, epsxx, ... (in general, all the names of the variables being saved).
    Purpose of the code:
    Returns:
    nt : Number of time steps saved.
    names_profiles : Names of the quantities being saved.
    dim : Should be the vertical dimension (nz).
    time_vec : Physical times at which data have been saved.
    '''
    time_vec = []
    names_profiles = []
    dim = 0

    with open(file, 'r') as file:
        for line in file:
            line = line.strip()

            # Check for the time variable header saved in the file
            if '=' in line:
                dim = 0
                t = line.split()[-1]
                time_vec.append(float(t))

            # Check for the names of variables saved in the file
            # Redundant but consistent since the output format is fixed
            elif '#' in line and '=' not in line:
                names_profiles = line.split()[1:]

            # Count rows for the vertical dimension nz
            elif '#' not in line and '=' not in line:
                data_row = list(map(float, line.split()))
                if len(data_row) == 0:
                    continue
                else:
                    dim += 1

    # Number of time steps saved
    nt = len(time_vec)

    return nt, names_profiles, dim, time_vec

