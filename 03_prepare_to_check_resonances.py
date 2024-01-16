import pandas as pd
import numpy as np
import aa_utilities as ut
#%%
clone_counts_file = 'clone_counts.csv'
df = pd.read_csv(clone_counts_file)
des_list = df['packed_designation'].tolist()
resonant_list = df['resonant_count'].tolist()
resonant_array = np.array(resonant_list)
resonant_max = np.max(resonant_array)
half_max = resonant_max/2
Nobj = len(des_list)
Njobs = 1000
obj_per_instance = int(np.ceil(Nobj/Njobs))
Njobs_real = int(np.ceil(Nobj/obj_per_instance))
epoch_file = 'epoch_settings.csv'
df2 = pd.read_csv(epoch_file)
year = int(df2['year'][0])
month = int(df2['month'][0])
day = int(df2['day'][0])
hour = int(df2['hour'][0])
minute = int(df2['minute'][0])
second = int(df2['second'][0])
Nclones = int(df2['Nclones'][0])
datestr = ut.datestr([year,month,day,hour,minute,second],'short')
JD = ut.JDfun(year,month,day,hour,minute,second)
template_file = '03a_check_resonances_template.py'
with open(template_file,'r') as file:
    template_data = file.read()
for i in range(Njobs_real):
    ct = i + 1
    ctstr = str(ct)
    outfile = 'check_resonances_' + ctstr + '.py'
    filedata = template_data
    filedata = filedata.replace('THIS_INSTANCE = 1','THIS_INSTANCE = '+str(int(i+1)))
    filedata = filedata.replace('Njobs = 300','Njobs = '+str(Njobs_real))
    filedata = filedata.replace('Nclones = 300','Nclones = '+str(Nclones))
    with open(outfile,'w') as file:
        file.write(filedata)
    start_obj = (ct-1) * obj_per_instance
    stop_obj = ct * obj_per_instance
    if stop_obj > Nobj:
        stop_obj = Nobj
    print(i+1,Njobs_real,start_obj,stop_obj,Nobj)
# Generate a batch submission file for the cluster to classify the clones with.
# On the cluster, run "sbatch slurm_classify_clones.slurm".
template_file = '03b_slurm_check_resonances_template.slurm'
with open(template_file,'r') as file:
    template_data = file.read()
filedata = template_data
filedata = filedata.replace('NJOBS',str(Njobs_real))
outfile = 'slurm_check_resonances.slurm'
with open(outfile,'w') as file:
    file.write(filedata)
