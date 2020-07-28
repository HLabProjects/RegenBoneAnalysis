# -*- coding: utf-8 -*-
"""
Created on Tue Jul 28 13:01:22 2020
@author: KFH
Script to calculate and construct the negative void volume of the mouse digits 
using microCT images 
Note: depends on data from mdbands.py

"""
import numpy as np
from scipy import ndimage
import cv2
import os
import datetime
import glob
import math
from skimage import morphology

ipath="your_ct_img_path"      #this is the path to your microCT image files
imgdir="your_ct_img_folder"
outpath="your_file_outpath"       #this is the path to where your output files go

bmpmatches=glob.glob(ipath+'/*[0-9].tif')
numbmp=len(bmpmatches)
sample_img=bmpmatches[numbmp//2]

#check if files exist
chk=os.path.exists(sample_img)
if chk==True:      #if there then proceed to import
    print('confirmed')
    #print(image+' submitted...')
    print(datetime.datetime.now())
    img = cv2.imread(sample_img,0)
    
if chk==False:
    print('no file found')

"""
-Processing to account for microCT stack dimensions - make largest x-y dimension into h
-Setting size of computable roi to iterate through image stacks
   -initilize arrays for data in x,y, and z
"""
hh, ww = img.shape[:2]
if hh > ww:
   h=hh
elif ww > hh:
   h=ww

resol=3    #reconstruction resolution
div=numbmp/resol

roii=resol
zsized=numbmp/roii
xsized=hh/roii
ysized=ww/roii
if zsized > xsized and zsized > ysized:
   fdiv=int(zsized)
elif xsized > zsized and xsized > ysized:
   fdiv=int(xsized)
elif ysized > xsized and ysized > zsized:
   fdiv=int(ysized)

hhindex=math.floor(hh/roii)
wwindex=math.floor(ww/roii)

"""
Voxel Calibration
"""
pxc=3.9  #from microCT scanner datafile, this will change with different machines
h2=hh*3.9
w2=ww*3.9
roii_area=(resol**(2))*pxc
roii_vol=(resol**(3))*pxc


"""
Main processing
-see flowchart in published manuscript for additional information
"""

M=np.zeros([fdiv,fdiv,fdiv])
Ma=M>=25
Mc=M>=25
Mint=M>=25
ctr=0
reg_tot=0
closed_tot=0
reg_vol_tot=0
closed_vol_tot=0
internal_vol_tot=0
internal_tot=0
hhindex=math.floor(hh/roii)
wwindex=math.floor(ww/roii)
threshold=55  

a_thres=27000/(resol**2)   #relate area closing parameters to resolution

for k in range(0,round(div)-1):
    image_current=bmpmatches[k*resol]
    img = cv2.imread(image_current,0)
    for j in range(0,wwindex):
        for i in range(0,hhindex):   
            roi=img[i*roii:(i+1)*roii, j*roii:(j+1)*roii]
            hold=np.average(roi, axis=1)
            value=np.average(hold)
            thresh=threshold  
            if value > thresh:
                M[k,i,j]=value
                
    mask1=M[k,:,:]>=55
    Ma[k,:,:]=mask1
    #core operations: 
    opened_mask = ndimage.binary_opening(mask1, iterations=1)
    closed_mask = ndimage.binary_closing(opened_mask, iterations=11)
    clean=morphology.remove_small_holes(closed_mask, area_threshold=a_thres)
    chull=morphology.convex_hull_image(mask1)
    chull_diff=np.logical_xor(chull,clean)
    cdiff_close=ndimage.binary_dilation(chull_diff, iterations=4)
    interior_p1=np.logical_xor(clean,mask1)
    extra1=np.logical_and(interior_p1,mask1)
    interior_p2=np.logical_xor(interior_p1,extra1)
    extra2=np.logical_and(interior_p2,cdiff_close)
    interior=np.logical_xor(interior_p2,extra2)
    Mint[k,:,:]=interior
    #now volumes sections
    reg_sum=np.sum(mask1)
    reg_tot=reg_tot+reg_sum
    closed_sum=np.sum(clean)
    closed_tot=closed_tot+closed_sum
    internal_sum=np.sum(interior)
    internal_tot=internal_tot+internal_sum
                   
  

vfract = internal_tot/(internal_tot+reg_tot)   #calculate void fraction

#reporting
print('regular volume is: ',reg_tot, ' voxels')
print('internal volume is: ',internal_tot, ' voxels')
print('total volume is: ',closed_tot, ' voxels')

#output
ocheck=1
if ocheck==1:
   np.save(outpath+imgdir+"IVTot_r.npy",internal_tot)
   np.save(outpath+imgdir+"CVTot_r.npy",closed_tot)
   np.save(outpath+imgdir+"RVTot_r.npy",reg_tot)
   np.save(outpath+imgdir+"M_r.npy",M)
   np.save(outpath+imgdir+"Ma_r.npy",Ma)
   np.save(outpath+imgdir+"Mint_r.npy",Mint)
   np.save(outpath+imgdir+"vfract.npy",vfract)
   


print("script complete")
