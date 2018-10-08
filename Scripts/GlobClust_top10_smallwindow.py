from gaia_tools import query
import matplotlib
matplotlib.use("AGG")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from astropy.table import Table
import pyexcel

GlobClust_Log = pd.read_csv('~/GlobClust_l-b-r_h.csv')
#GlobClust_Log = GlobClust_Log.set_index("ID", drop = True) #--> List by ID for index
GlobClust_Log_short = GlobClust_Log.loc[0:10,"ID":"r_h"]
#GlobClust_Log_short = GlobClust_Log.loc["NGC 104":"NGC 362","L":"B"]
GlobClust_Log_l = GlobClust_Log_short.loc[:,"L"]
GlobClust_Log_b = GlobClust_Log_short.loc[:,"B"]
GlobClust_Log_noID = GlobClust_Log_short.loc[:,"L":]
GlobClust_Log_ID = GlobClust_Log_short.loc[:,"ID"]
GlobClust_Log_r_h = GlobClust_Log_short.loc[:,'r_h']

for i in range(0,10): #Should be better
    print(GlobClust_Log_ID[i])
    print(GlobClust_Log_l[i])
    print(GlobClust_Log_b[i])
    
    globquery="""select ra, dec, l, b, pmra, pmdec from gaiadr2.gaia_source where abs(l-%f)<%f and\
    abs(b-(%f))<%f""" %(GlobClust_Log_l[i],GlobClust_Log_r_h[i]/60*3,GlobClust_Log_b[i],GlobClust_Log_r_h[i]/60*3)
    
    out=query.query(globquery,local=True,timeit=True)
    
    df_t=pd.DataFrame(np.array(out))
    df_t.to_csv('out'+ GlobClust_Log_ID[i] +'full_r_h-small.csv', index=True, header=True)
    
    fig=plt.figure(figsize=(16,12))
    plt.scatter(x=out['l'], y=out['b'], c=out['pmra'], cmap='inferno')#, vmax=10, vmin=-10)
    plt.xlabel('l (deg)', fontsize=15)
    plt.ylabel('b (deg)', fontsize=15)
    plt.title('Small Window (3*r_h):' + 'b (deg) vs l (deg) of ' + GlobClust_Log_ID[i], fontsize=15) #GlobClust_Log_ID[i]
    cb=plt.colorbar()
    cb.ax.tick_params(labelsize=15)
    cb.set_label(r'$\mathrm{pmra}$',fontsize=15)
    plt.gca().tick_params(labelsize=15)
    #plt.gca().set_ylim([-20.4,-20.55])
    #plt.gca().set_xlim([307.25,307.45])
    plt.show()
    fig.savefig('GC_'+ GlobClust_Log_ID[i] +'full_pmra_r_h-small.pdf')
    
    fig=plt.figure(figsize=(16,12))
    plt.scatter(x=out['l'], y=out['b'], c=out['pmdec'], cmap='inferno')#, vmax=10, vmin=-10)
    plt.xlabel('l (deg)', fontsize=15)
    plt.ylabel('b (deg)', fontsize=15)
    plt.title('Small Window (3*r_h):' + 'b (deg) vs l (deg) of ' + GlobClust_Log_ID[i], fontsize=15) #GlobClust_Log_ID[i]
    cb=plt.colorbar()
    cb.ax.tick_params(labelsize=15)
    cb.set_label(r'$\mathrm{pmdec}$',fontsize=15)
    plt.gca().tick_params(labelsize=15)
    #plt.gca().set_ylim([-20.4,-20.55])
    #plt.gca().set_xlim([307.25,307.45])
    plt.show()
    fig.savefig('GC_'+ GlobClust_Log_ID[i] +'full_pmdec_r_h-small.pdf')
