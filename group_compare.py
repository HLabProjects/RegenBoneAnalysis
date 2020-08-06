# -*- coding: utf-8 -*-
"""
Created on Wed Jul 29 13:19:57 2020
@author: KFH
Script for comparing multiple digitsin pooled fashion, using .npy data files
and plotting different looks at at mineral density in mouse bone regeneration images
"""

from scipy.stats import norm
from scipy import stats
import re
import numpy as np
import matplotlib.pyplot as plt
import os
import math

path="your_ct_img_path"      #this is the path to your microCT image files
imgdir="your_ct_img_folder"
outpath="your_file_outpath"       #this is the path to where your output files go

"""
section to automatically find the subfolders needed, using regex
"""
search_key="naming_part"
title1=('Your_Group_Folder_1')
single1_path=(path+title1+'\\')      
single1_dirlist=[x for x in os.listdir(single1_path) if re.search(search_key,x)]  
group1=single1_dirlist

title2=('Your_Group_Folder_2')
single2_path=(path+title2+'\\')      
single2_dirlist=[x for x in os.listdir(single2_path) if re.search(search_key,x)]  
group2=single2_dirlist

#these are the groups chosen
num1=group1
num2=group2

#what resolution?
resol=3    ######


def average_datacall(file_list_name):
   for i in range(0,len(file_list_name)):
      if i==0:
            first=file_list_name[0]
            second=file_list_name[1]
            GAR1=np.load(outpath+first+"GAR.npy")
            VV1=(GAR1[:,3]-57.7)*.0177 + 0.75
            GAR2=np.load(outpath+second+"GAR.npy")
            VV2=(GAR2[:,3]-57.7)*.0177 + 0.75
            VVtot=np.concatenate([VV1,VV2])
      if i>1:
            current=file_list_name[i]
            GAR1=np.load(outpath+current+"GAR.npy")
            VV1=(GAR1[:,3]-57.7)*.0177 + 0.75
            VVtot=np.concatenate([VVtot,VV1])
      
   return VVtot

digi_stats_list2=[]
def average_datacall2(file_list_name):
   digi_stats_list=[]
   for i in range(0,len(file_list_name)):
      first=file_list_name[i]
      GAR1=np.load(outpath+first+"GAR.npy")
      VV1=(GAR1[:,3]-57.7)*.0177 + 0.75
      g_one=VV1
      #standard deviation
      stdev_yh=np.std(g_one, axis=0)
      #population variance
      yh_var1=np.var(g_one, axis=0)
      yh_var2=g_one.var(ddof=1)
      #square root of the count
      sqrt_cnt_yh=math.sqrt(len(g_one))
      #standard error of the mean
      sem_yh=stats.sem(g_one)
      #average
      avg_yh=g_one.mean()
      #count
      yh_size=g_one.size
      digi_stats_list.append([first,avg_yh,yh_size,sqrt_cnt_yh,stdev_yh,sem_yh,yh_var1])
   return digi_stats_list

#Calculations
lowerlim=0.65
upperlim=3.8
bin_number=32
bin_array=np.arange(lowerlim,upperlim,(upperlim-lowerlim)/(bin_number+1))

g_one=average_datacall(num1)
g_two=average_datacall(num2)
stat_list=average_datacall2(num2)
np.save(outpath+title1+"pool",num1)
np.save(outpath+title2+"pool",num2)


#add bestfit line
#run statistics

#standard deviation
stdev_yh=np.std(g_one, axis=0)
stdev_oh=np.std(g_two, axis=0)

#population variance
yh_var1=np.var(g_one, axis=0)
oh_var1=np.var(g_two, axis=0)
yh_var2=g_one.var(ddof=1)
oh_var2=g_two.var(ddof=1)

#square root of the count
sqrt_cnt_yh=math.sqrt(len(g_one))
sqrt_cnt_oh=math.sqrt(len(g_two))

#standard error of the mean
sem_yh=stats.sem(g_one)
sem_oh=stats.sem(g_two)

#average
avg_yh=g_one.mean()
avg_oh=g_two.mean()

#count
yh_size=g_one.size
oh_size=g_two.size

"""
Plotting
"""
fig, axes=plt.subplots(1,3, figsize=(10.5,3.5), constrained_layout=True)

hist1=axes[0].hist(g_one,bins=32,density=True, facecolor='dimgrey', edgecolor="black", linewidth=0.5, alpha=.75)
hist2=axes[0].hist(g_two,bins=32,density=True, facecolor='lightgrey', edgecolor="black", linewidth=0.5, alpha=.5)
(mu_a, sigma_a)=norm.fit(g_two)
(mu_b, sigma_b)=norm.fit(g_one)
bfy_a=norm.pdf(hist1[1], mu_a, sigma_a)
bfy_b=norm.pdf(hist2[1], mu_b, sigma_b)
bf1=axes[0].plot(hist2[1],bfy_a,'k--', linewidth=1.5,label='Best Fit 1')
bf2=axes[0].plot(hist1[1],bfy_b,color='dimgray',ls='-.', linewidth=1.5,label='Best Fit 2')
title1a=title1
title2a=title2
axes[0].legend((hist1[2][0],bf1[0],hist2[2][0],bf2[0],), 
          (title1a,'Normal Fit\nMean=%.3f \nSD=%.5f' % (avg_yh,stdev_yh),title2a,
           'Normal Fit\nMean=%.3f \nSD=%.5f' % (avg_oh,stdev_oh),),fontsize=9)

axes[0].set_xlabel('Mineral Density (g/cm^3)',fontsize=9)
axes[0].set_xlim(0.8, 3.8)
axes[0].set_ylabel('Probability',fontsize=9)
#plt.tight_layout()
#plt.show()


#violin and box
labels=(title1, title2)
vio=axes[1].violinplot([g_two,g_one],showmeans=False,showmedians=True)
axes[1].set_ylabel('Mineral Density (g/cm^3)',fontsize=9)

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
    #ax.set_xlabel('Sample name')
    
red_square = dict(markerfacecolor='r', marker='s')
bp=axes[2].boxplot([g_two,g_one], showmeans=False, meanprops=red_square, showfliers=True, patch_artist=True)
axes[2].set_title('box plot')
axes[2].set_ylabel('Mineral Density (g/cm^3)',fontsize=9)
for box in bp['boxes']:
   # change outline color
    box.set( color='darkblue', linewidth=2)
    # change fill color
    box.set( facecolor = 'skyblue' )
            
## change color and linewidth of the caps
for cap in bp['caps']:
    cap.set(color='black', linewidth=2)

## change color and linewidth of the medians
for median in bp['medians']:
    median.set(color='red', linewidth=3)
               

#the following is to annotate figs with stat sig
def stars(p):
   if p < 0.0001:
       return "****"
   elif (p < 0.001):
       return "***"
   elif (p < 0.01):
       return "**"
   elif (p < 0.05):
       return "*"
   else:
       return "-"
    
def t_test_call(groupa,groupb):
   t, p = stats.ttest_ind(groupa, groupb, equal_var=False)
   return t,p

def stat_sig_graph(group_uno, group_duo,ggset,axis):
   y_max = np.max(np.concatenate((group_uno, group_duo)))
   y_min = np.min(np.concatenate((group_uno, group_duo)))
   ttestr=t_test_call(group_uno, group_duo)
   x_ss=ggset
   buff=y_max*1.2
   if ttestr[1] < 0.0001:
      pcapt="< 0.0001"
   elif ttestr[1] < 0.001:
      pcapt="< 0.001"
   elif ttestr[1] > 0.001:
      pcapt=("= "+str(ttestr[1]))
   axes[axis].set_ylim(0.5,buff)
   axes[axis].annotate("", xy=(x_ss[0], y_max), xycoords='data',
           xytext=(x_ss[1], y_max), textcoords='data',
           arrowprops=dict(arrowstyle="-", ec='#aaaaaa',lw=2, 
                           connectionstyle="bar,fraction=0.2"))
   axes[axis].text((x_ss[0]+x_ss[1])/2, y_max + abs(y_max - y_min)*0.13, stars(ttestr[1]),
           horizontalalignment='center',
           verticalalignment='center', size=14)
   axes[axis].text((x_ss[0]+x_ss[1])/2, y_max + abs(y_max - y_min)*0.03, "p "+pcapt,
           horizontalalignment='center',
           verticalalignment='center', size=9)
   
ggset=[1,2]
group1=g_two
group2=g_one
stat_sig_graph(group1, group2,ggset,1)
stat_sig_graph(group1, group2,ggset,2)
set_axis_style(axes[1],labels)
set_axis_style(axes[2],labels)
#plt.tight_layout()
title=('CompMD '+title1+' '+title2)
plt.savefig(outpath+title, dpi=500)

print("plotting complete")
print("script complete")




