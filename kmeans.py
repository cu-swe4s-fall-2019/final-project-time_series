import time
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# preprocessing data for clust
raw_counts= pd.read_csv('raw_counts.txt',sep='\t')
# remove columns
raw_counts = raw_counts.drop(['AFFY_ID','ACCNUM'], axis=1)
cols = raw_counts.columns[1:]
for i in cols:
    raw_counts[i] = pd.to_numeric(raw_counts[i], errors='raise')
raw_counts['SYMBOL'] = raw_counts['SYMBOL'].replace([' ',',',';'], '_', regex=True)
raw_counts.to_csv('raw_count2.txt',columns =['SYMBOL', '0.15H.IFN_1', '0.5H.IFN_1', '0.75H.IFN_1',
       '1H.IFN_1', '1.25H.IFN_1', '1.5H.IFN_1', '2.25H.IFN_1', '2.75H.IFN_1',
       '3H.IFN_1', '3.5H.IFN_1', '5H.IFN_1', '5.5H.IFN_1', '6H.IFN_1',
       '6.5H.IFN_1', '7H.IFN_1', '8H.IFN_1', '9H.IFN_1', '10H.IFN_1',
       '11H.IFN_1', '12H.IFN_1', '13H.IFN_1', '14H.IFN_1', '15H.IFN_1'], index = False, sep='\t')
