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

import sklearn
from sklearn import cluster
from astropy import units as u
from astropy.coordinates import SkyCoord
from matplotlib.patches import Circle
#from matplotlib import patches

#%pylab inline

import pywt
from galpy.util import bovy_plot

# csv-pd set up for DBSCANN
GlobClust_Log_ra_dec = pd.read_csv('~/GlobClust_ra-dec.csv')
GlobClust_Log_r_t = GlobClust_Log_ra_dec.loc[0:0,"r_t"]
GlobClust_Log_radec = GlobClust_Log_ra_dec.loc[0:0,"RA":"DEC"]
GlobClust_Log_ID = GlobClust_Log_ra_dec.loc[0:0,"ID"]
GlobClust_Log_r_c = GlobClust_Log_ra_dec.loc[0:0,"r_c"]
GlobClust_Log_r_t = GlobClust_Log_ra_dec.loc[0:0,"r_t"]
#GlobClust_Log_dec = GlobClust_Log_ra_dec.loc[:,"dec"]
#GlobClust_Log_RA = GlobClust_Log_ra_dec.loc[0:7,"RA"]
#GlobClust_Log_DEC = GlobClust_Log_ra_dec.loc[0:7,"DEC"]
#print(GlobClust_Log_r_t)

GlobClust_Log_2 = pd.read_csv('~/out_NGC 288_full_12-rt_ra_dec.csv')
GlobClust_Log_pmra = GlobClust_Log_2.loc[:,"pmra"]
GlobClust_Log_pmdec = GlobClust_Log_2.loc[:,"pmdec"]
GlobClust_Log_ra = GlobClust_Log_2.loc[:,"ra"]
GlobClust_Log_dec = GlobClust_Log_2.loc[:,"dec"]

GlobClust_Log_g = GlobClust_Log_2.loc[:,"phot_g_mean_mag"]
GlobClust_Log_rp = GlobClust_Log_2.loc[:,"phot_rp_mean_mag"]
GlobClust_Log_bp = GlobClust_Log_2.loc[:,"phot_bp_mean_mag"]
GlobClust_Log_rpbp = GlobClust_Log_bp - GlobClust_Log_rp
GlobClust_Log_ra_out = GlobClust_Log_2.loc[:,"ra"]
GlobClust_Log_dec_out = GlobClust_Log_2.loc[:,"dec"]
GlobClust_Log_pmdec_out = GlobClust_Log_2.loc[:,"pmdec"]
GlobClust_Log_pmra_out = GlobClust_Log_2.loc[:,"pmra"]

#GlobClust_Log_pmra = GlobClust_Log_pmra.as_matrix()
#GlobClust_Log_pmdec = GlobClust_Log_pmdec.as_matrix()
#GlobClust_Log_ra = GlobClust_Log_ra.as_matrix()
#GlobClust_Log_dec = GlobClust_Log_dec.as_matrix()

GlobClust_Log_parallax_out = GlobClust_Log_2.loc[:,"parallax"]
GlobClust_Log_parallaxerror_out = GlobClust_Log_2.loc[:,"parallax_error"]
GlobClust_Log_parallaxe_out = GlobClust_Log_2.loc[:,"parallax":"parallax_error"]

# DBSCANN - parallax before scan
GlobClust_Log_parallaxover_out = GlobClust_Log_parallaxerror_out/GlobClust_Log_parallax_out

#parallax_indx = (GC_parallaxover>0.20) | ((1/GC_parallax>5.)*(GC_parallaxover<0.20))
parallax_indx=((1/GlobClust_Log_parallax_out<5.)*(GlobClust_Log_parallaxover_out<0.10))
parallax_indx=[not i for i in parallax_indx]

coord_ra_dec_deg=[]
coord_ra_dec = SkyCoord(GlobClust_Log_radec.loc[:,"RA"], GlobClust_Log_radec.loc[:,"DEC"], frame='icrs')

Clust_Edge = GlobClust_Log_r_c[0]
print(Clust_Edge)
d_pc=8900

#M=m-5(log(d)-1)
GlobClust_Log_M = GlobClust_Log_g-5*((np.log10(d_pc)) - 1)

indx = (np.sqrt((np.fabs(GlobClust_Log_ra_out-coord_ra_dec.ra.deg)**2+np.fabs(GlobClust_Log_dec_out-coord_ra_dec.dec.deg)**2))<Clust_Edge/60)

# DBSCAN
DB_Params =[]
DB_Params=np.transpose(DB_Params)

DF={'M': GlobClust_Log_M[parallax_indx]/5, 'rpbp': GlobClust_Log_rpbp[parallax_indx]/2.5, 
    'pmra': GlobClust_Log_pmra[parallax_indx]/20, 'pmdec': GlobClust_Log_pmdec[parallax_indx]/20}

DB_Params = pd.DataFrame(data=DF)
DB_Params=DB_Params.dropna()

print(np.shape(DB_Params))

db=sklearn.cluster.DBSCAN(eps=0.031, min_samples=18, metric='euclidean', metric_params=None, algorithm='auto', 
                      leaf_size=30, p=None, n_jobs=1).fit(DB_Params) #original eps=0.1625

print(db.labels_)
print(np.shape(db.labels_))

fig=plt.figure(figsize=(18,12))
plt.scatter(GlobClust_Log_rpbp/2.5, GlobClust_Log_M/5, c=GlobClust_Log_rpbp/2.5, cmap='RdBu_r')
#plt.scatter(db,DB_Params)
#plt.show()

core_samples_mask = np.zeros_like(db.labels_, dtype=bool)
core_samples_mask[db.core_sample_indices_] = True
labels = db.labels_

# Number of clusters in labels, ignoring noise if present.
n_clusters_ = len(set(labels)) - (1 if -1 in labels else 0)

print('Estimated number of Subgroups: %d' % n_clusters_)

x_param3=[]
y_param3=[]

x_param4=[]
y_param4=[]

# Black removed and is used for noise instead.
unique_labels = set(labels)
colors = [plt.cm.Spectral(each)
          for each in np.linspace(0, 1, len(unique_labels))]
for k, col in zip(unique_labels, colors):
    if k == -1:
        # Black used for noise.
        col = [0, 0, 0, 1]

    class_member_mask = (labels == k)
    
    #NOISE
    #xy = DB_Params[class_member_mask & ~core_samples_mask]
    #plt.plot(xy.loc[:, "rpbp"], xy.loc[:, "M"], 'o', markerfacecolor=tuple(col),
    #         markeredgecolor='k', markersize=6)
    
    #CLEANED DATA
    xy = DB_Params[class_member_mask & core_samples_mask]
    
    plt.plot(xy.loc[:, "rpbp"], xy.loc[:, "M"], 'o', markerfacecolor=tuple(col),
             markeredgecolor='k', markersize=14)
    
    x_param3.append(xy.loc[:, "rpbp"]*2.5)
    y_param3.append(xy.loc[:, "M"]*5)

    #x_param4.append(xy.loc[:, "pmra"]*0.35)
    #y_param4.append(xy.loc[:, "pmdec"]*0.3)
    
    #x_param5.append(xy.loc[:, "ra"]*0.35)
    #y_param5.append(xy.loc[:, "dec"]*0.3)

plt.gca().invert_yaxis()
plt.title('Estimated number of Subgroups: %d' % n_clusters_)
#cb=plt.colorbar()
#cb.ax.tick_params(labelsize=15)
plt.xlabel('bp-rp', fontsize=15)
plt.ylabel('M', fontsize=15)

plt.legend()
plt.show()
fig.savefig('GC_DBSCAN_'+ GlobClust_Log_ID[0] +'_CMD_raw.png')
plt.close()

x_out3 = pd.DataFrame(x_param3)#, index=None)#, mangle_dupe_cols=True)
#x_out3.reset_index()
x_out3 = x_out3.transpose() 
x_out3 = x_out3.loc[:,~x_out3.columns.duplicated()]
#x_out3 = x_out3["rpbp"].as_matrix()
x_out3

y_out3 = pd.DataFrame(y_param3)#, index=None)#, mangle_dupe_cols=True)
#y_out3.reset_index()
y_out3 = y_out3.transpose() 
y_out3 = y_out3.loc[:,~y_out3.columns.duplicated()]
#y_out3 = y_out3["M"].as_matrix()
print(np.shape(y_out3))

GC_ra = pd.DataFrame({'ra':GlobClust_Log_ra})
GC_dec = pd.DataFrame({'dec':GlobClust_Log_dec})
GC_pmra = pd.DataFrame({'pmra':GlobClust_Log_pmra})
GC_pmdec = pd.DataFrame({'pmdec':GlobClust_Log_pmdec})

ra_DB = GC_ra.reindex(x_out3.index)
dec_DB = GC_dec.reindex(y_out3.index)
pmra_DB = GC_pmra.reindex(x_out3.index)
pmdec_DB = GC_pmdec.reindex(y_out3.index)
M_DB = GlobClust_Log_M.reindex(y_out3.index)
rpbp_DB = GlobClust_Log_rpbp.reindex(y_out3.index)

pmra_DB_norm = (pmra_DB-pmra_DB.mean())/(pmra_DB.mean())
pmdec_DB_norm = (pmdec_DB-pmdec_DB.mean())/(pmdec_DB.mean())

# DISTANCE LOOP FROM NEW CMD OF NGC 288 USING r_c AS CONSTRAINT FOR OVERPLOTTED CLUSTER STARS ON CMD (MATCHING)
# NEEDS BE RUN AFTER PREVIOUS CELL!!!
import time

t0 = time.time()
Cluster = np.zeros(len(GlobClust_Log_M), dtype='bool')
M_i = GlobClust_Log_M
M_k = M_DB
rpbp_i = GlobClust_Log_rpbp
rpbp_k = rpbp_DB

M_i = M_i.as_matrix()
M_k = M_k.as_matrix()
rpbp_i = rpbp_i.as_matrix()
rpbp_k = rpbp_k.as_matrix()

#M_k[1]

#sqrt(np.fabs(M_i[1]-M_k[1])**2+(rpbp_i[1]-rpbp_k[1])**2)
#rpbp_i

#dist=[]
#dist=(fabs(rpbp_i[1]-rpbp_k[1])**2)
#dist
#dist = sqrt(np.fabs((M_i[0]-M_k[0])**2+(rpbp_i[0]-rpbp_k[0])**2))

M_final=[]
rpbp_final=[]

print('Marching Start.....')

for i in range(0, len(M_i)): #,21399):#len(M_i)):
    #Cluster[i]=False
    for k in range(0, len(M_k)): #,3302):#len(M_k)):
#        print([i], [k])
#        if M_k[k] == 243:
#            print([i], [k], ': 10% done')
#        elif M_k[k] == 243*2:
#            print([i], [k], ': 20% done')
#        elif M_k[k] == 243*3:
#            print([i], [k], ': 30% done')
#        elif M_k[k] == 243*4:
#            print([i], [k], ': 40% done')
#        elif M_k[k] == 243*5:
#            print([i], [k], ': 50% done')
#        elif M_k[k] == 243*6:
#            print([i], [k], ': 60% done')
#        elif M_k[k] == 243*7:
#            print([i], [k], ': 70% done')
#        elif M_k[k] == 243*8:
#            print([i], [k], ': 80% done')
#        elif M_k[k] == 243*9:
#            print([i], [k], ': 90% done')
#        elif M_k[k] == 243*10:
#            print([i], [k], ': 100% done')
        #print(M_k[k])
        #print(rpbp_k[k])
        #print(M_i[i])
        #print(rpbp_i[i])
        #print(M_i[i][~np.isnan(M_i[i])])
        #print(rpbp_i[i][~np.isnan(rpbp_i[i])])
        #if rpbp_i[i][~np.isnan(rpbp_i[i])] & rpbp_k[k][~np.isnan(rpbp_k[k])]:
        dist = np.sqrt(np.fabs((M_i[i]-M_k[k])**2+(rpbp_i[i]-rpbp_k[k])**2))
        #print(dist)
        
        if dist < 0.1/2:
            Cluster[i]=True
            #Cluster[i]
            #np.append(M_final,M_k[k])
            #np.append(rpbp_final,rpbp_k[k])
            
            M_final.append(M_i[i])
            rpbp_final.append(rpbp_i[i])
            
            break


t1 = time.time()
total = t1-t0
print(total)
print(Cluster)
print(np.shape(Cluster))
print(np.sum(Cluster))

# MUST BE RUN AFTER PREVIOUS CELLS

#Clust_indx = (sqrt(np.fabs((M_i-M_k)**2+(rpbp_i-rpbp_k)**2))<0.1)

rpbp_ic=rpbp_i[Cluster] #np.multiply(rpbp_i,Cluster)
M_ic=M_i[Cluster] #multiply(M_i,Cluster)

#for i in range(0, len(M_ic)):
#    if rpbp_ic[i]!=0. and M_ic[i]!=0.:
#        rpbp_ic[i]=rpbp_ic[i]
#        M_ic[i]=M_ic[i]



fig=plt.figure(figsize=(16,12))
plt.scatter(x=GlobClust_Log_rpbp, y=GlobClust_Log_M, c=GlobClust_Log_rpbp, cmap='RdBu_r', marker='.', label='Raw Data')
plt.scatter(x=rpbp_final, y=M_final, c='purple', alpha=0.7, marker='.', label = 'After matching')
plt.scatter(x=rpbp_ic, y=M_ic, c='g', alpha=0.2, marker='.', label='After matching')
plt.xlabel('bp-rp', fontsize=15)
plt.ylabel('M', fontsize=15)
#plt.gca().invert_yaxis()
#plt.title('CMD of ' + GC_IDs[i], fontsize=15)
#cb=plt.colorbar()
#cb.ax.tick_params(labelsize=15)
#cb.set_label(r'$\mathrm{pmdec}$',fontsize=15)
plt.gca().tick_params(labelsize=15)
plt.gca().invert_yaxis()
plt.legend()
plt.show()
fig.savefig('GC_DBSCAN_'+ GlobClust_Log_ID[0] +'_CMD_matched.png')
plt.close()

#GC_StarsPlot = {'bprp indx': GlobClust_Log_rpbp[indx], 'M indx': GlobClust_Log_M[indx], 'ra': GlobClust_Log_ra_out,\
#                    'dec': GlobClust_Log_dec_out}
#GC_StarscutPlot = pd.DataFrame(data=GC_StarsPlot)
#GC_StarscutPlot.to_csv('out_GCStars_'+ 'NGC 288' +'_cutplot2_r_c.csv', index=True, header=True)

# PLOT RA-DEC FROM MATCHING
# RUN AFTER PREVIOUS CELLS

ra_ic=GlobClust_Log_ra_out[Cluster] #np.multiply(rpbp_i,Cluster)
dec_ic=GlobClust_Log_dec_out[Cluster] #multiply(M_i,Cluster)

fig=plt.figure(figsize=(16,12))
plt.scatter(x=GlobClust_Log_ra_out, y=GlobClust_Log_dec_out, c='c', alpha=0.3, marker='.', label = 'Raw Data')
#plt.scatter(x=rpbp_final2, y=M_final2, c='purple', alpha=0.7)
plt.scatter(x=ra_ic, y=dec_ic, c='purple', marker='.', label='After Matching')#, alpha=0.2)
plt.xlabel('ra', fontsize=15)
plt.ylabel('dec', fontsize=15)
#plt.gca().invert_yaxis()
#plt.title('CMD of ' + GC_IDs[i], fontsize=15)
#cb=plt.colorbar()
#cb.ax.tick_params(labelsize=15)
#cb.set_label(r'$\mathrm{pmdec}$',fontsize=15)
plt.gca().tick_params(labelsize=15)
plt.legend()
plt.show()
fig.savefig('GC_DBSCAN_'+ GlobClust_Log_ID[0] +'_scatterplot.png')
plt.close()
