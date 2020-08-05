import pandas as pd
import os

path = '/home/amukhopadhyay/nucleic_acids/test_results'

all_files = []
for i in os.listdir(path):
    if os.path.isfile(os.path.join(path,i)) and 'angles' in i:
        all_files.append(i)
print(all_files)
dfs = list()
for filename in all_files:
    print(filename)
    df = pd.read_csv(os.path.join(path, filename))
    dfs.append(df)
frame = pd.concat(dfs, axis=0, ignore_index=True)
print(frame.shape)
frame.to_csv('merged_test_angles.csv', index=False)
