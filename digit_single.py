# -*- coding: utf-8 -*-
"""
Created on Tue Jul 28 12:42:15 2020
@author: KFH
This script is for comparing bmd in multiple digits, using .npy data files 
created from processing of microCt images

NOTE: this depends on variables and data created by mdbands.py
@author: KFH
"""
from scipy.stats import norm
import scipy.stats as stats
import numpy as np
import matplotlib.pyplot as plt

path="your_ct_img_path"      #this is the path to your microCT image files
imgdir="your_ct_img_folder"
outpath="your_file_outpath"       #this is the path to where your output files go

first=imgdir    
curr_name=imgdir.rsplit(' ',1)
firstname=curr_name[0]

"""
Plotting
"""
s_1=8
a_1=.8
a_2=.4

fig=plt.figure(figsize=(7.0,3.5), constrained_layout=True)
Xa=np.load(outpath+first+"X.npy")
Ya=np.load(outpath+first+"Y.npy")
Za=np.load(outpath+first+"Z.npy")
Da=np.load(outpath+first+"D"+".npy")
resol=3   
roii=resol
fdiva=np.size(Da,1)
GARa=np.load(outpath+first+"GAR.npy")
VVa=(GARa[:,3]-57.7)*.0177 + 0.75    # corrected 4/28/20

#violin plot
ax = fig.add_subplot(121)
title1=firstname
labels=(title1,)
vio=ax.violinplot([VVa],showmeans=False,showmedians=True)
ax.set_ylabel('Mineral Density (g/cm^3)',fontsize=9)
ax.set_ylim(0.5, 4.0)
#ax.set_aspect(1.0)
plt.show()

for pc in vio['bodies']:
    pc.set_facecolor('skyblue')
    pc.set_edgecolor('black')
    pc.set_alpha(1)

def set_axis_style(ax, labels):
    ax.get_xaxis().set_tick_params(direction='out')
    ax.xaxis.set_ticks_position('bottom')
    ax.set_xticks(np.arange(1, len(labels) + 1))
    ax.set_xticklabels(labels)
    ax.set_xlim(0.25, len(labels) + 0.75)
    
    
set_axis_style(ax,labels)

#now the histogram
ax = fig.add_subplot(122)
hist1=ax.hist(VVa,bins=32,density=True, facecolor='dimgrey', edgecolor="black", linewidth=0.5, alpha=.75)
#add bestfit line
(mu_a, sigma_a)=norm.fit(VVa)
bfy_a=norm.pdf(hist1[1], mu_a, sigma_a)
younghind=VVa
#standard error of the mean
sem_yh=stats.sem(younghind)
#standard deviation
stdev_yh=np.std(younghind, axis=0)
#average
avg_yh=younghind.mean()
bf1=ax.plot(hist1[1],bfy_a,'k--', linewidth=1.5,label='Normal Fit')   #old non bw
title=firstname
ax.legend((hist1[2][0],bf1[0],), 
          (title,'Normal Fit\nMean=%.3f \nSD=%.3f' % (avg_yh,stdev_yh),
          ),fontsize=9)

ax.set_xlabel('Mineral Density (g/cm^3)',fontsize=9)
ax.set_xlim(0.6, 3.4)
ax.set_ylabel('Probability',fontsize=9)
#ax.set_ylim(0.0, 3.0)
plt.show()


plt.savefig(outpath+first + '_d_single', bbox_inches='tight',dpi=500)
print("script complete")
