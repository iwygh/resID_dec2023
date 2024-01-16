#%%
# all objects
# asteroids
# all objects
# exclude comet fragments
# e<=1
# a range 25:200
# asteroid basic
# select all
# full precision
# This file finds objects in the SBDB and MPC databases that might be Plutinos
# and prepares for them to be classified on the Puma cluster. After running
# this file, simply upload the entire folder to Puma and run the command
# 'sbatch slurm_classify_clones.slurm'.
import aa_utilities as ut
import pandas as pd
import numpy as np
# This is the date we downloaded the SBDB and MPC databases. It is the epoch
# that will be used to compute population-level statistics and make population-
# level figures after we have identified the Plutino sample.
year = 2023
month = 12
day = 2
hour = 0
minute = 0
second = 0
Nclones = 300
Njobs = 1000
# Save settings for epoch and number of clones.
dictionary = {'year':[year],'month':[month],'day':[day],'hour':[hour],'minute':[minute],\
              'second':[second],'Nclones':[Nclones]}
df_out = pd.DataFrame(dictionary)
df_out.to_csv('epoch_settings.csv',index=False)
# Make a string of the date, for use in filenames.
datestr = ut.datestr([year,month,day,hour,minute,second],'short')
# We will need the Julian date to look up orbital elements at the epoch
# chosen for making population-level figures and statistics.
JD = ut.JDfun(year,month,day,hour,minute,second)
# Retrieve all the information we really need from the MPCORB file.
print('starting ut.dat2csv')
ut.dat2csv(datestr)
print('done with ut.dat2csv')
# Reduce entire SBDB catalog to a smaller population as described in README.docx.
# Each object that remains will be checked to see if it is Resonant or not.
print('starting ut.sbdb_reduce_34_150')
Nstart,Nobj = ut.sbdb_reduce_34_150(JD,datestr)
print('starting number of objects',Nstart)
print('objects in amin,amax range',Nobj)
# Generate clones as described in README.docx.
Nobj_cloned = ut.generate_clones(JD,datestr,Nclones)
print('number of objects cloned',Nobj_cloned)
# Figure out how many jobs to run.
obj_per_instance = int(np.ceil(Nobj_cloned/Njobs))
Njobs_real = int(np.ceil(Nobj_cloned/obj_per_instance))
print('nominal number of jobs',Njobs)
print('objects per job',obj_per_instance)
print('actual number of jobs',Njobs_real)
# # Generate as many instances of classify_clones_#.py as there are jobs to run.
template_file = '01a_classify_clones_template.py'
with open(template_file,'r') as file:
    template_data = file.read()
for i in range(Njobs_real):
    ct = i + 1
    ctstr = str(ct)
    outfile = 'classify_clones_' + ctstr + '.py'
    filedata = template_data
    filedata = filedata.replace('THIS_INSTANCE = 1','THIS_INSTANCE = '+str(int(i+1)))
    filedata = filedata.replace('Njobs = 300','Njobs = '+str(Njobs_real))
    with open(outfile,'w') as file:
        file.write(filedata)
    start_obj = (ct-1) * obj_per_instance
    stop_obj = ct * obj_per_instance
    if stop_obj > Nobj_cloned:
        stop_obj = Nobj_cloned
    print(i+1,Njobs_real,start_obj,stop_obj,Nobj_cloned)
# Generate a batch submission file for the cluster to classify the clones with.
# On the cluster, run "sbatch slurm_classify_clones.slurm".
template_file = '01b_slurm_classify_clones_template.slurm'
with open(template_file,'r') as file:
    template_data = file.read()
filedata = template_data
filedata = filedata.replace('NJOBS',str(Njobs_real))
outfile = 'slurm_classify_clones.slurm'
with open(outfile,'w') as file:
    file.write(filedata)
