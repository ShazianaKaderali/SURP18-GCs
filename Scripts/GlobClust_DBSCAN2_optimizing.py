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
from matplotlib.patches import Circle
#from matplotlib import patches
#%pylab inline
import pywt
from galpy.util import bovy_plot

import PyPDF2
from PyPDF2 import PdfFileReader, PdfFileWriter

from scipy import stats as st
from matplotlib import colors, ticker, cm

import csv

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

GlobClust_Log_2 = pd.read_csv('~/out_NGC 288_full_15rt_ra_dec.csv')
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

GlobClust_Log_pmra_error = GlobClust_Log_2.loc[:,"pmra_error"]
GlobClust_Log_pmdec_error = GlobClust_Log_2.loc[:,"pmdec_error"]
GlobClust_Log_ra_error = GlobClust_Log_2.loc[:,"ra_error"]
GlobClust_Log_dec_error = GlobClust_Log_2.loc[:,"dec_error"]
GlobClust_Log_vpu = GlobClust_Log_2.loc[:,"visibility_periods_used"]
GlobClust_Log_aen = GlobClust_Log_2.loc[:,"astrometric_excess_noise"]

#indx = (np.sqrt((np.fabs(GlobClust_Log_ra_out-coord_ra_dec.ra.deg)**2+np.fabs(GlobClust_Log_dec_out-coord_ra_dec.dec.deg)**2))<Clust_Edge/60)

# DBSCANN - parallax before scan
GlobClust_Log_parallaxover_out = GlobClust_Log_parallaxerror_out/GlobClust_Log_parallax_out

#parallax_indx = (GC_parallaxover>0.20) | ((1/GC_parallax>5.)*(GC_parallaxover<0.20)), /3600000
parallax_indx=((1/(GlobClust_Log_parallax_out)<5.)*(GlobClust_Log_parallaxover_out<0.10))
parallax_indx=[not i for i in parallax_indx]

coord_ra_dec_deg=[]
coord_ra_dec = SkyCoord(GlobClust_Log_radec.loc[:,"RA"], GlobClust_Log_radec.loc[:,"DEC"], frame='icrs')

Clust_Edge = GlobClust_Log_r_c[0]
print(Clust_Edge)
Tidal_Edge = GlobClust_Log_r_t[0]
print(Tidal_Edge)
d_pc=8900

#M=m-5(log(d)-1)
GlobClust_Log_M = GlobClust_Log_g-5*((np.log10(d_pc)) - 1)

indx = (np.sqrt((np.fabs(GlobClust_Log_ra_out-coord_ra_dec.ra.deg)**2+np.fabs(GlobClust_Log_dec_out-coord_ra_dec.dec.deg)**2))<Clust_Edge/60)

r_t = GlobClust_Log_r_t.as_matrix()
r_t = r_t[0]

# Creating a routine that appends files to the output file
def append_pdf(input,output):
    [output.addPage(input.getPage(page_num)) for page_num in range(input.numPages)]

# Creating an object where pdf pages are appended to
output = PdfFileWriter()

data = np.loadtxt('n288_gaia.dat')
np.shape(data)

data_pd = pd.DataFrame({'id': data[:,0], 'x': data[:,1], 'y': data[:,2], 'z': data[:,3],
                       'vx': data[:,4], 'vy': data[:,5], 'vy': data[:,6], 
                       'ra': data[:,7], 'dec':data[:,8], 'pmra': data[:,9], 'pmdec': data[:,10]})

ra_n = data_pd.loc[:,"ra"]
#print(ra_n)
dec_n = data_pd.loc[:,"dec"]
#print(dec_n)
pmra_n = data_pd.loc[:,"pmra"]
#print(pmra_n)
pmdec_n = data_pd.loc[:,"pmdec"]
#print(pmdec_n)

core_indx = (np.sqrt((np.fabs(ra_n.as_matrix()-coord_ra_dec.ra.deg-0.85)**2+np.fabs(dec_n.as_matrix()-coord_ra_dec.dec.deg-0.1)**2)) <= GlobClust_Log_r_c[0]/60)

tidal_indx = (np.sqrt((np.fabs(ra_n.as_matrix()-coord_ra_dec.ra.deg-0.85)**2+\
                       np.fabs(dec_n.as_matrix()-coord_ra_dec.dec.deg-0.1)**2)) <= GlobClust_Log_r_t[0]/60)

tidal_indx_15 = (np.sqrt((np.fabs(ra_n.as_matrix()-coord_ra_dec.ra.deg-0.85)**2+\
                       np.fabs(dec_n.as_matrix()-coord_ra_dec.dec.deg-0.1)**2)) <= GlobClust_Log_r_t[0]/60*15)


pmra_n_norm = pmra_n-pmra_n[core_indx].mean()
pmdec_n_norm = pmdec_n-pmdec_n[core_indx].mean()

tidal_indx2 = (np.sqrt((np.fabs(ra_n.as_matrix()-coord_ra_dec.ra.deg-0.85)**2+\
                       np.fabs(dec_n.as_matrix()-coord_ra_dec.dec.deg-0.1)**2)) <= GlobClust_Log_r_t[0]/60*15) *\
    (np.sqrt((np.fabs(ra_n.as_matrix()-coord_ra_dec.ra.deg-0.85)**2+\
                       np.fabs(dec_n.as_matrix()-coord_ra_dec.dec.deg-0.1)**2)) > GlobClust_Log_r_t[0]/60)

#Optimization
epsilon = np.arange(0.01, 0.052, 0.001)
min_samps = np.arange(10, 52, 1)

pm_indx=(GlobClust_Log_pmra_error[parallax_indx] < 1) * (GlobClust_Log_pmdec_error[parallax_indx] < 1) *\
        (GlobClust_Log_vpu[parallax_indx] > 8) * (GlobClust_Log_aen[parallax_indx] == 0) *\
        (GlobClust_Log_pmra[parallax_indx] > 2.5) * (GlobClust_Log_pmra[parallax_indx] < 6) *\
        (GlobClust_Log_pmdec[parallax_indx] > -7.5) * (GlobClust_Log_pmdec[parallax_indx] < -3.5) #*\
        #(GlobClust_Log_M[parallax_indx] < 4.5)

tidal_indx_sq = (np.fabs(ra_n-coord_ra_dec.ra.deg-0.85) < r_t/60*15 ) *\
                (np.fabs(dec_n-coord_ra_dec.dec.deg-0.1) < r_t/60*15 )
    
tidal_indx_sq2 = (np.fabs(ra_n-coord_ra_dec.ra.deg-0.85)<r_t/60*15) * (np.fabs(dec_n-coord_ra_dec.dec.deg-0.1)<r_t/60*15) *\
                 (np.sqrt((np.fabs(ra_n-coord_ra_dec.ra.deg-0.85))**2+(np.fabs(dec_n-coord_ra_dec.dec.deg-0.1))**2) > r_t/60)

for ii in range(0,len(epsilon)):
    for j in range(0,len(min_samps)):
        print([ii], [j])
        # DBSCANN - parallax before scan
        GlobClust_Log_parallaxover_out = GlobClust_Log_parallaxerror_out/GlobClust_Log_parallax_out
        
        #parallax_indx = (GC_parallaxover>0.20) | ((1/GC_parallax>5.)*(GC_parallaxover<0.20))
        parallax_indx2=((1/GlobClust_Log_parallax_out<5.)*(GlobClust_Log_parallaxover_out<0.10))
        parallax_indx2=[not i for i in parallax_indx2]

        coord_ra_dec_deg=[]
        coord_ra_dec = SkyCoord(GlobClust_Log_radec.loc[:,"RA"], GlobClust_Log_radec.loc[:,"DEC"], frame='icrs')
        print(coord_ra_dec)

        DB_Params =[]
        DB_Params=np.transpose(DB_Params)

        DF={'M': GlobClust_Log_M[parallax_indx][pm_indx]/7, 'rpbp': GlobClust_Log_rpbp[parallax_indx][pm_indx], 
            'pmra': GlobClust_Log_pmra[parallax_indx][pm_indx]/50, 'pmdec': GlobClust_Log_pmdec[parallax_indx][pm_indx]/50}


        DB_Params = pd.DataFrame(data=DF)
        DB_Params=DB_Params.dropna()

        print(DB_Params.shape)

        db=sklearn.cluster.DBSCAN(eps=epsilon[ii], min_samples=min_samps[j], metric='euclidean', metric_params=None, algorithm='auto',
                      leaf_size=30, p=None, n_jobs=1).fit(DB_Params) #original eps=0.1625

        #print(db.labels_)
        print((db.labels_).shape)

        fig=plt.figure(figsize=(18,12))
        plt.hexbin(x=GlobClust_Log_rpbp[parallax_indx], y=GlobClust_Log_M[parallax_indx]/7, gridsize=500,bins='log',cmap='gist_yarg', rasterized=True)#, extent=[-0.5,0.5,  0,2.5])
        
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

            plt.plot(xy.loc[:, "rpbp"], xy.loc[:, "M"], '.', markerfacecolor=tuple(col), rasterized=True)
                     #markeredgecolor='k', markersize=14,  )

            x_param3.append(xy.loc[:, "rpbp"])
            y_param3.append(xy.loc[:, "M"]*7)


        plt.gca().invert_yaxis()
        plt.title('Estimated number of Subgroups: %d, [%.3f] [%d]' % (n_clusters_, epsilon[ii], min_samps[j]))
        #cb=plt.colorbar()
        #cb.ax.tick_params(labelsize=15)
        plt.xlabel('bp-rp', fontsize=15)
        plt.ylabel('M', fontsize=15)

        #plt.legend()
        #plt.show()
        fig.savefig('GC_DBSCAN_'+ GlobClust_Log_ID.loc[0] +'_CMD_15rt_'+ str(epsilon[ii])  + '_' +  str(min_samps[j]) +'.pdf')
        plt.close()


        x_out3 = pd.DataFrame(x_param3)#, index=None)#, mangle_dupe_cols=True)
        x_out3.reset_index()
        x_out3 = x_out3.transpose()
        x_out3 = x_out3.loc[:,~x_out3.columns.duplicated()]
        
        y_out3 = pd.DataFrame(y_param3)#, index=None)#, mangle_dupe_cols=True)
        y_out3.reset_index()
        y_out3 = y_out3.transpose()
        y_out3 = y_out3.loc[:,~y_out3.columns.duplicated()]
        #y_out3 = y_out3["M"].as_matrix()
        #print(shape(y_out3))

        GC_ra = pd.DataFrame({'ra':GlobClust_Log_ra})
        GC_dec = pd.DataFrame({'dec':GlobClust_Log_dec})
        GC_pmra = pd.DataFrame({'pmra':GlobClust_Log_pmra})
        GC_pmdec = pd.DataFrame({'pmdec':GlobClust_Log_pmdec})
        GC_pmra_error = pd.DataFrame({'pmra_error':GlobClust_Log_pmra_error})
        GC_pmdec_error = pd.DataFrame({'pmdec_error':GlobClust_Log_pmdec_error})

        ra_DB = GC_ra.reindex(x_out3.index)
        dec_DB = GC_dec.reindex(y_out3.index)
        pmra_DB = GC_pmra.reindex(x_out3.index)
        pmdec_DB = GC_pmdec.reindex(y_out3.index)
        pmra_DB_error = GC_pmra_error.reindex(x_out3.index)
        pmdec_DB_error = GC_pmdec_error.reindex(y_out3.index)
        M_DB = GlobClust_Log_M.reindex(y_out3.index)
        rpbp_DB = GlobClust_Log_rpbp.reindex(y_out3.index)
        
        pmra_DB_norm = (pmra_DB-pmra_DB.mean())#/(pmra_DB.mean())
        pmdec_DB_norm = (pmdec_DB-pmdec_DB.mean())#/(pmdec_DB.mean())

        fig = plt.figure(figsize=(16,12))
        ax = plt.gca()
        #plt.scatter(ra_DB-ra_DB.mean(), dec_DB-dec_DB.mean(), c='r', alpha=0.2, label='original from DBSCANN' )
        plt.scatter(ra_DB-coord_ra_dec.ra.deg, dec_DB-coord_ra_dec.dec.deg, c=pmdec_DB_norm, marker='o', s=100,
                cmap='inferno_r', vmin=-1, vmax=1, label='eps= %.3f, minsamp= %d, error < 1' %(epsilon[ii], min_samps[j]),  rasterized=True)#, label='with parallax criteria after DBSCANN')
        circle = plt.Circle((0, 0), GlobClust_Log_r_t[0]/60, color='k', fill=False, lw=5)
        #circle2 = plt.Circle((0, 0), 4*GlobClust_Log_r_t[0]/60, color='b', fill=False, lw=5)
        ax.add_artist(circle)
        #ax.add_artist(circle2)
        plt.title('[%d]: eps=%.3f; [%d]: min samps=%f' % (ii, epsilon[ii], j, min_samps[j]))
        #plt.plot(x, y, label='y = %.2f x + %.2f' %(A, B))
        #('eps='+ epsilon[i] +', min samps='+ min_samps[i]) #'out_GCStars_'+ GC_IDs[i] +'_cutplot_r_c.csv'
        cb=plt.colorbar()
        cb.ax.tick_params(labelsize=15)
        cb.set_label(r'$\mathrm{pmdec}$',fontsize=15)
        plt.xlabel('ra', fontsize=15)
        plt.ylabel('dec', fontsize=15)
        #plt.legend()
        #plt.show()
        fig.savefig('GC_DBSCAN_'+ GlobClust_Log_ID.loc[0] +'_radec-pmdec_15rt_'+ str(epsilon[ii])  + '_' +  str(min_samps[j]) +'.pdf')
        plt.close()

        fig = plt.figure(figsize=(16,12))
        plt.scatter(pmra_DB, pmdec_DB, marker='.', rasterized=True)
        plt.title('pmdec vs pmra, [%d] [%d], [%.3f] [%d]' %(ii,j, epsilon[ii], min_samps[j]))
        plt.xlabel('pmra', fontsize=15)
        plt.ylabel('pmdec', fontsize=15)
        fig.savefig('GC_DBSCAN_'+ GlobClust_Log_ID.loc[0] +'_pms_15rt_'+ str(epsilon[ii]) + '_' +  str(min_samps[j]) +'.pdf')
        plt.close()
        
        clust_indx = (np.sqrt((np.fabs(ra_DB.as_matrix()-coord_ra_dec.ra.deg)**2+np.fabs(dec_DB.as_matrix()-coord_ra_dec.dec.deg)**2))
              >Tidal_Edge/60)
        
        quad_I_n = ((ra_n[tidal_indx_sq2]-0.85) > coord_ra_dec.ra.deg[0]) *\
                   ((dec_n[tidal_indx_sq2]-0.1) > coord_ra_dec.dec.deg[0])
    
        quad_II_n = ((ra_n[tidal_indx_sq2]-0.85) < coord_ra_dec.ra.deg[0]) *\
                    ((dec_n[tidal_indx_sq2]-0.1) > coord_ra_dec.dec.deg[0])
    
        quad_III_n = ((ra_n[tidal_indx_sq2]-0.85) < coord_ra_dec.ra.deg[0]) *\
                     ((dec_n[tidal_indx_sq2]-0.1) < coord_ra_dec.dec.deg[0])
    
        quad_IV_n = ((ra_n[tidal_indx_sq2]-0.85) > coord_ra_dec.ra.deg[0]) *\
                    ((dec_n[tidal_indx_sq2]-0.1) < coord_ra_dec.dec.deg[0])
    

        quad_I_DB = ((ra_DB.as_matrix()[clust_indx]) > coord_ra_dec.ra.deg[0]) *\
                    ((dec_DB.as_matrix()[clust_indx]) > coord_ra_dec.dec.deg[0])
    
        quad_II_DB = ((ra_DB.as_matrix()[clust_indx]) < coord_ra_dec.ra.deg[0]) *\
                     ((dec_DB.as_matrix()[clust_indx]) > coord_ra_dec.dec.deg[0])
    
        quad_III_DB = ((ra_DB.as_matrix()[clust_indx]) < coord_ra_dec.ra.deg[0]) *\
                      ((dec_DB.as_matrix()[clust_indx]) < coord_ra_dec.dec.deg[0])
    
        quad_IV_DB = ((ra_DB.as_matrix()[clust_indx]) > coord_ra_dec.ra.deg[0]) *\
                     ((dec_DB.as_matrix()[clust_indx]) < coord_ra_dec.dec.deg[0])
        
        
        D, QI_pmdec_ks_p = st.stats.ks_2samp(pmdec_n_norm[tidal_indx_sq2][quad_I_n], np.asarray(pmdec_DB_norm[clust_indx][quad_I_DB])[:,0])
        D, QII_pmdec_ks_p = st.stats.ks_2samp(pmdec_n_norm[tidal_indx_sq2][quad_II_n], np.asarray(pmdec_DB_norm[clust_indx][quad_II_DB])[:,0])
        D, QIII_pmdec_ks_p = st.stats.ks_2samp(pmdec_n_norm[tidal_indx_sq2][quad_III_n], np.asarray(pmdec_DB_norm[clust_indx][quad_III_DB])[:,0])
        D, QIV_pmdec_ks_p = st.stats.ks_2samp(pmdec_n_norm[tidal_indx_sq2][quad_IV_n], np.asarray(pmdec_DB_norm[clust_indx][quad_IV_DB])[:,0])
        
        D, QI_pmra_ks_p = st.stats.ks_2samp(pmra_n_norm[tidal_indx_sq2][quad_I_n], np.asarray(pmra_DB_norm[clust_indx][quad_I_DB])[:,0])
        D, QII_pmra_ks_p = st.stats.ks_2samp(pmra_n_norm[tidal_indx_sq2][quad_II_n], np.asarray(pmra_DB_norm[clust_indx][quad_II_DB])[:,0])
        D, QIII_pmra_ks_p = st.stats.ks_2samp(pmra_n_norm[tidal_indx_sq2][quad_III_n], np.asarray(pmra_DB_norm[clust_indx][quad_III_DB])[:,0])
        D, QIV_pmra_ks_p = st.stats.ks_2samp(pmra_n_norm[tidal_indx_sq2][quad_IV_n], np.asarray(pmra_DB_norm[clust_indx][quad_IV_DB])[:,0])
                
        D, QI_dec_ks_p = st.stats.ks_2samp(dec_n[tidal_indx_sq2][quad_I_n], np.asarray(dec_DB[clust_indx][quad_I_DB])[:,0])
        D, QII_dec_ks_p = st.stats.ks_2samp(dec_n[tidal_indx_sq2][quad_II_n], np.asarray(dec_DB[clust_indx][quad_II_DB])[:,0])
        D, QIII_dec_ks_p = st.stats.ks_2samp(dec_n[tidal_indx_sq2][quad_III_n], np.asarray(dec_DB[clust_indx][quad_III_DB])[:,0])
        D, QIV_dec_ks_p = st.stats.ks_2samp(dec_n[tidal_indx_sq2][quad_IV_n], np.asarray(dec_DB[clust_indx][quad_IV_DB])[:,0])

        D, QI_ra_ks_p = st.stats.ks_2samp(ra_n[tidal_indx_sq2][quad_I_n], np.asarray(ra_DB[clust_indx][quad_I_DB])[:,0])
        D, QII_ra_ks_p = st.stats.ks_2samp(ra_n[tidal_indx_sq2][quad_II_n], np.asarray(ra_DB[clust_indx][quad_II_DB])[:,0])
        D, QIII_ra_ks_p = st.stats.ks_2samp(ra_n[tidal_indx_sq2][quad_III_n], np.asarray(ra_DB[clust_indx][quad_III_DB])[:,0])
        D, QIV_ra_ks_p = st.stats.ks_2samp(ra_n[tidal_indx_sq2][quad_IV_n], np.asarray(ra_DB[clust_indx][quad_IV_DB])[:,0])
        
        row = [str(ii+1),\
               str(QI_pmdec_ks_p), str(QII_pmdec_ks_p), str(QIII_pmdec_ks_p), str(QIV_pmdec_ks_p),\
               str(QI_dec_ks_p), str(QII_dec_ks_p), str(QIII_dec_ks_p), str(QIV_dec_ks_p),\
               str(QI_pmra_ks_p), str(QII_pmra_ks_p), str(QIII_pmra_ks_p), str(QIV_pmra_ks_p),\
               str(QI_ra_ks_p), str(QII_ra_ks_p), str(QIII_ra_ks_p), str(QIV_ra_ks_p)]
        
        with open('ks_test_ps.csv', 'w') as csvFile:
            writer = csv.writer(csvFile)
            writer.writerow(row)
            
        csvFile.close()

        str_name_cmd='GC_DBSCAN_'+ GlobClust_Log_ID.loc[0] +'_CMD_15rt_'+ str(epsilon[ii]) + '_' + str(min_samps[j]) + '.pdf'
        str_name_radec='GC_DBSCAN_'+ GlobClust_Log_ID.loc[0] +'_radec-pmdec_15rt_'+ str(epsilon[ii]) + '_' + str(min_samps[j]) +'.pdf'
        str_name_pms="GC_DBSCAN_"+ GlobClust_Log_ID.loc[0] +"_pms_15rt_"+ str(epsilon[ii]) + '_' + str(min_samps[j]) +".pdf"

        # Appending PDFs of GC's of Interest
        CMD_stream = open(str_name_cmd, "rb")
        radec_stream = open(str_name_radec, "rb")
        pms_stream = open(str_name_pms, "rb")
        
        CMD_input = PdfFileReader(CMD_stream)
        radec_input = PdfFileReader(radec_stream)
        pms_input = PdfFileReader(pms_stream)
        
        #append_pdf(PdfFileReader(open(str_name, "rb")), output)
        #append_pdf(PdfFileReader(open("GC_DBSCAN_"+ GlobClust_Log_ID.loc[0] +"_CMD_15rt_"+ str(ii)+ '_' +str(j) + ".pdf","rb")),output)
        #append_pdf(PdfFileReader(open("GC_DBSCAN_"+ GlobClust_Log_ID.loc[0] +"_radec-pmdec_15rt_"+ str(ii)+ '_' +str(j) +".pdf","rb")),output)
        #append_pdf(PdfFileReader(open("GC_DBSCAN_"+ GlobClust_Log_ID.loc[0] +"_pms_15rt_"+ str(ii)+ '_' +str(j) +".pdf","rb")),output)
        #append_pdf(PdfFileReader(open("GC_DBSCAN_"+ GlobClust_Log_ID.loc[0] +"_density_15rt_"+ str(ii)+ '_' +str(j) +".pdf","rb")),output)
        
        append_pdf(CMD_input, output)
        append_pdf(radec_input, output)
        append_pdf(pms_input, output)

        # Writing all the collected pages to a file
        output.write(open("GC_NGC288_DBSCANs2.pdf","wb"))
        
        CMD_stream.close()
        radec_stream.close()
        pms_stream.close()
        
        
        #str_name='GC_DBSCAN_'+ GlobClust_Log_ID.loc[0] +'_CMD_12rt_'+ str(ii) + str(j) + '.pdf'

        # Appending PDFs of GC's of Interest
        #append_pdf(PdfFileReader(open(str_name, "rb")), output)
        #append_pdf(PdfFileReader(open("GC_DBSCAN_"+ GlobClust_Log_ID.loc[0] +"_CMD_12rt_"+ str(ii) + str(j) +".pdf","rb")),output)
        #'GC_DBSCAN_'+ GlobClust_Log_ID.loc[0] +'_CMD_12rt_'+ str(ii) + str(j) +
        #append_pdf(PdfFileReader(open("GC_DBSCAN_"+ GlobClust_Log_ID.loc[0] +"_CMD_12rt_"+ str(ii) + str(j) + ".pdf","rb")),output)
        #append_pdf(PdfFileReader(open("GC_DBSCAN_"+ GlobClust_Log_ID.loc[0] +"_radec-pmdec_12rt_"+ str(ii) + str(j) +".pdf","rb")),output)
        #append_pdf(PdfFileReader(open("GC_DBSCAN_"+ GlobClust_Log_ID.loc[0] +"_pms_12rt_"+ str(ii) + str(j) +".pdf","rb")),output)
        #append_pdf(PdfFileReader(open("GC_DBSCAN_"+ GlobClust_Log_ID.loc[0] +"_density_12rt_"+ str(ii) + str(j) +".pdf","rb")),output)


        # Writing all the collected pages to a file
        #output.write(open("GC_NGC288_DBSCANs.pdf","wb"))
