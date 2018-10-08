from gaia_tools import query
import matplotlib
matplotlib.use("AGG")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from astropy.table import Table
import pyexcel
from astropy import units as u
from astropy.coordinates import SkyCoord

GlobClust_Log_ra_dec = pd.read_csv('~/GlobClust_ra-dec.csv')
GlobClust_Log_ra_dec = GlobClust_Log_ra_dec.loc[0:7,"ID":"DEC"]
GlobClust_Log_RA = GlobClust_Log_ra_dec.loc[:,"RA"]
GlobClust_Log_DEC = GlobClust_Log_ra_dec.loc[:,"DEC"]
GlobClust_Log_ID = GlobClust_Log_ra_dec.loc[:,"ID"]
GlobClust_Log_r_t = GlobClust_Log_ra_dec.loc[:,"r_t"]

coord_ra_dec = SkyCoord(GlobClust_Log_RA, GlobClust_Log_DEC, frame='icrs')

for i in range(0,8): #Should be better
    print(GlobClust_Log_ID[i])
    print(GlobClust_Log_RA[i])
    print(GlobClust_Log_DEC[i])

    globquery="""select ra, dec, l, b, pmra, pmdec, phot_g_mean_mag, phot_rp_mean_mag, phot_bp_mean_mag from gaiadr2.gaia_source where abs(ra-%f)<%f and abs(dec-(%f))<%f""" %(coord_ra_dec.ra.deg[i],GlobClust_Log_r_t[i]/60*4,coord_ra_dec.dec.deg[i],GlobClust_Log_r_t[i]/60*4)

    out=query.query(globquery,local=True,timeit=True)

    df_t=pd.DataFrame(np.array(out))
    df_t.to_csv('out_'+ GlobClust_Log_ID[i] +'_full_4-rt_ra_dec.csv', index=True, header=True)
    
    fig=plt.figure(figsize=(16,12))
    plt.hexbin(x=out['ra'], y=out['dec'], cmap='plasma', bins='log')
    plt.xlabel('ra (deg)', fontsize=15)
    plt.ylabel('dec (deg)', fontsize=15)
    plt.title('Window (4*r_t):' +'ra (deg) vs dec (deg) of ' + GlobClust_Log_ID[i], fontsize=15)
    cb=plt.colorbar()
    cb.ax.tick_params(labelsize=15)
    cb.set_label(r'$\mathrm{pmra}$',fontsize=15)
    plt.gca().tick_params(labelsize=15)
    plt.show()
    fig.savefig('GC_'+ GlobClust_Log_ID[i] +'_full_pmra_4-rt_ra_dec.pdf')

    fig=plt.figure(figsize=(16,12))
    plt.scatter(x=out['ra'], y=out['dec'], c=out['pmdec'], cmap='inferno')
    plt.xlabel('ra (deg)', fontsize=15)
    plt.ylabel('dec (deg)', fontsize=15)
    plt.title('Window (4*r_t):' + 'ra (deg) vs dec (deg) of ' + GlobClust_Log_ID[i], fontsize=15)
    cb=plt.colorbar()
    cb.ax.tick_params(labelsize=15)
    cb.set_label(r'$\mathrm{pmdec}$',fontsize=15)
    plt.gca().tick_params(labelsize=15)
    plt.show()
    fig.savefig('GC_'+ GlobClust_Log_ID[i] +'_full_pmdec_4-rt_ra_dec.pdf')
