#!/usr/bin/env python3

import pandas as pd 

table=pd.read_csv("/srv/nfs/ngs-stockage/NGS_Virologie/HadrienR/PIPELINE_NGS/COUNT_MAPPING/S20107_S17.txt",sep=";",names=['COUNT','MAPPING'])  
table = table.loc[table['MAPPING'] != "*"]
print(table)
SubType = table['MAPPING'].values[0]
print(SubType)

