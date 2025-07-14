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

def latex_format_name(names_profiles,target_old_name,target_new_name):
    #length of the common name
    name_vec  = []
    nl = len(target_old_name)
    for name in names_profiles:
        if(name.startswith(target_old_name)):
            tmp = "\\" +  target_new_name + "_{" + name[nl:] + "}"
            name_vec.append(tmp)
        else:
            name_vec.append(name)
    return name_vec

def save_plot_data(data, name_old, name_new, nt, DeltaNt, dt, figr_folder,param_plot):

    nvariables = data.shape[1]
    zc = data[0, 0, :]

    for nvar in range(1, nvariables):

        name_old_tmp = name_old[nvar]
        name_new_tmp = name_new[nvar]

        fig, ax = plt.subplots()
        fig_name = name_old_tmp + ".jpg"

        for p in range(0, nt, DeltaNt):
            time = p * dt
            with plt.style.context(['science', 'ieee']):
                ax.plot(zc, data[p, nvar, :],label=time)
        
        ax.autoscale(tight=True)
        ax.set(**param_plot)
        ax.legend()
        fig.savefig(figr_folder + fig_name, dpi=300)
        plt.close(fig)


data_folder = "../data/"
figr_folder = "./figures/"
dt = 1
DeltaNt = 1

#clean the folder where data will be saved
path_fgr = os.path.abspath(figr_folder) 
subprocess.run(['rm', '-r',path_fgr])
subprocess.run(['mkdir',path_fgr])


filename_vec          = ["spectra.dat"]
target_old_vec_name   = ["none"]
target_new_vec_name   = ["none"]


#-------------------------------------------------------------------------------------------------#
#MANIPULATE SPECTRA
file            = data_folder + filename_vec[0]
target_old_name = target_old_vec_name[0]
target_new_name = target_new_vec_name[0] 

nt, name_old, nz,time_vec = read_file_to_store_index(file)
data                      = generate_data_structure(file,nt,nz,name_old)
name_new                  = latex_format_name(name_old,target_old_name,target_new_name)
save_plot_data(data, name_old, name_new, nt,DeltaNt,dt,figr_folder,param_plot)

#keep verbosity for clarity
k_vec  = data[:,0,:]
E_u    = data[:,1,:]; E_ku   = data[:,2,:]
E_v    = data[:,3,:]; E_kv   = data[:,4,:]
E_w    = data[:,5,:]; E_kw   = data[:,6,:]
E_c    = data[:,7,:]; E_kc   = data[:,8,:]

lambda_u = np.sum(E_ku,axis=1) / np.sum(E_u,axis=1)
lambda_v = np.sum(E_kv,axis=1) / np.sum(E_v,axis=1)
lambda_w = np.sum(E_kw,axis=1) / np.sum(E_w,axis=1)
lambda_c = np.sum(E_kc,axis=1) / np.sum(E_c,axis=1)


nvariables = data.shape[1]
    zc = data[0, 0, :]

    for nvar in range(1, nvariables):

        name_old_tmp = name_old[nvar]
        name_new_tmp = name_new[nvar]

        fig, ax = plt.subplots()
        fig_name = name_old_tmp + ".jpg"

        for p in range(0, nt, DeltaNt):
            time = p * dt
            with plt.style.context(['science', 'ieee']):
                ax.plot(zc, data[p, nvar, :],label=time)
        
        ax.autoscale(tight=True)
        param_plot = dict(ylabel=r"$||\boldsymbol{k}||$", xlabel=fr"${{{target_new_name}}}$")
        ax.set(**param_plot)
        ax.legend()
        fig.savefig(figr_folder + fig_name, dpi=300)
        plt.close(fig)




#filename = data_folder + "pressure.dat"
#nt, name_old_pressure, nz, time_vec = read_file_to_store_index(filename)
#data_pressure                       = generate_data_structure(filename,nt,nz,name_old_pressure)
#name_new_pressure                   = latex_format_name(name_old_pressure,"pres","pressure")
#save_plot_data(data_pressure, name_old_pressure, name_new_pressure, nt, DeltaNt, dt, figr_folder)
             






