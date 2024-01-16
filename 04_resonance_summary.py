import pandas as pd
cloned_objects_file = 'cloned_objects.csv'
df = pd.read_csv(cloned_objects_file)
des_list = df['cloned_objects'].tolist()
Nobj = df.shape[0]
which_obj_list = []
which_resonance_list = []
which_conditions_checked_list = []
count_not_resonant = 0
count_other = 0
count_32 = 0
count_21 = 0
count_52 = 0
for iobj in range(Nobj):
    des = des_list[iobj]
    infile = 'resonance_summary_' + des + '.csv'
    df2 = pd.read_csv(infile)
    which_obj_list.append(df2['which_obj'][0])
    which_resonance_list.append(df2['which_resonance'][0])
    foo = df2['which_resonance'][0]
    if foo == 'not_resonant':
        count_not_resonant = count_not_resonant + 1
    elif foo == 'other':
        count_other = count_other + 1
    elif str(foo) == '32':
        count_32 = count_32 + 1
    elif str(foo) == '21':
        count_21 = count_21 + 1
    elif str(foo) == '52':
        count_52 = count_52 + 1
    else:
        print('what the heck is happening here')
    which_conditions_checked_list.append(df2['which_conditions_checked'][0])
dictionary = {'which_des':des_list,\
              'which_obj':which_obj_list,\
              'which_resonance':which_resonance_list,\
              'which_conditions_checked':which_conditions_checked_list}
outfile = 'resonance_summary_000_all.csv'
df3 = pd.DataFrame.from_dict(dictionary)
df3.to_csv(outfile,index=False)
print(Nobj)
print(count_not_resonant,count_other,count_32,count_21,count_52)
print(count_not_resonant+count_other+count_32+count_21+count_52)
