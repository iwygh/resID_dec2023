#%% imports block
import pandas as pd
import numpy as np
import os
#%% delete .out files
# First, delete the .out files cluttering up the folder from the previous run on Puma.
thisdir = os.getcwd()
files_in_dir = os.listdir(thisdir)
for this_file in files_in_dir:
    if this_file.endswith('.out'):
        os.remove(os.path.join(thisdir,this_file))
#%% read clone classifications, keep only objects with ZERO Resonant clones
cloned_objects_file = 'cloned_objects.csv'
df = pd.read_csv(cloned_objects_file)
des_list = df['cloned_objects'].tolist()
Nobj = df.shape[0]
outfile = 'clone_counts.csv'
resonant_count_list = []
classical_count_list = []
scattering_count_list = []
detached_count_list = []
other_count_list = []
for iobj in range(Nobj):
    des = des_list[iobj]
    classes_file = 'class_lists_' + des + '.csv'
    df = pd.read_csv(classes_file)
    categories = df['category'].tolist()
    Nclones = df.shape[0]
    resonant = 0
    classical = 0
    scattering = 0
    detached = 0
    other = 0
    for iclone in range(Nclones):
        print(iobj,Nobj,iclone,Nclones)
        category = categories[iclone]
        category = category.lstrip()
        category = category.rstrip()
        if category == 'Resonant':
            resonant = resonant + 1
        elif category == 'Scattering':
            scattering = scattering + 1
        elif category == 'Classical':
            classical = classical + 1
        elif category == 'Detached':
            detached = detached + 1
        else:
            other = other + 1
    resonant_count_list.append(resonant)
    classical_count_list.append(classical)
    scattering_count_list.append(scattering)
    detached_count_list.append(detached)
    other_count_list.append(other)
print(np.max(np.array(resonant_count_list)))
print(np.max(np.array(classical_count_list)))
print(np.max(np.array(scattering_count_list)))
print(np.max(np.array(detached_count_list)))
print(np.max(np.array(other_count_list)))
dictionary = {'packed_designation':des_list,\
              'resonant_count':resonant_count_list,\
              'classical_count':classical_count_list,\
              'scattering_count':scattering_count_list,\
              'detached_count':detached_count_list,\
              'other_count':other_count_list}
df2 = pd.DataFrame.from_dict(dictionary)
df2.to_csv(outfile,index=None)
