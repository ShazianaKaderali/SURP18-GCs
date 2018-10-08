from gaia_tools import query
import matplotlib
matplotlib.use("AGG")
import matplotlib.pyplot as plt
import numpy as np

import pandas as pd
import numpy as np

GlobClust_Log = pd.read_csv('~/GlobClust_l-b-r_h-r_c.csv')
#GlobClust_Log = GlobClust_Log.set_index("ID", drop = True) #--> List by ID for index
GlobClust_Log_short = GlobClust_Log.loc[0:10,"ID":"B"]
#GlobClust_Log_short = GlobClust_Log.loc["NGC 104":"NGC 362","L":"B"]
GlobClust_Log_l = GlobClust_Log_short.loc[:,"L"]
GlobClust_Log_b = GlobClust_Log_short.loc[:,"B"]
GlobClust_Log_noID = GlobClust_Log_short.loc[:,"L":]
GlobClust_Log_ID = GlobClust_Log_short.loc[:,"ID"]

for i in range(0,10): #Should be better
    globquery="""select ra, dec, l, b, pmra, pmdec from gaiadr2.gaia_source where abs(l-%f)<0.1 and\
    abs(b-(%f))<0.1""" %(GlobClust_Log_l[i],GlobClust_Log_b[i])
    
    out=query.query(globquery,local=True,timeit=True)
    
