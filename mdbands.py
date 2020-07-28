# -*- coding: utf-8 -*-
"""
Created on Mon Jul 27 23:21:22 2020
@author: KFH
script to import greyscale images from CT scans of mouse bones and calculate
bone mineral density in 3d space for a given pixel/voxel/increment size
"""

from scipy.stats import norm
import cv2
import numpy as np
import os
import math
import time
import datetime
import glob
from tqdm import tqdm
from shapely.geometry import Point


path="your_ct_img_path"      #this is the path to your microCT image files
imgdir="your_ct_img_folder"
outpath="your_file_outpath"       #this is the path to where your output files go


bmpmatches=glob.glob(path+'/*[0-9].tif')   #grabs list of all the ctimages in the folder
#bmpmatches=glob.glob(path+'/*[0-9].bmp')
numbmp=len(bmpmatches)
sample_img=bmpmatches[numbmp//2]

chk=os.path.exists(sample_img)  
if chk==True:      #if sample image exists then proceed to import images in folder
    print('confirmed')
    img = cv2.imread(sample_img,0)

if chk==False:
    print('no file found')

hh, ww = img.shape[:2]
print('img directory is: ',imgdir)

"""
Voxel Calibration
"""
pxc=3.9   #from microCT scanner datafile, this will change with different machines
h2=hh*3.9
w2=ww*3.9

"""
-Processing to account for microCT stack dimensions - make largest x-y dimension into h
-Setting size of computable roi to iterate through image stacks
   -initilize arrays for data in x,y, and z
"""
if hh > ww:
   h=hh
elif ww > hh:
   h=ww
   
resol=3   #reconstruction resolution: number of voxel lengths per reconstruction length
div=numbmp/resol
imgdir=imgdir+str(resol)

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


"""create the roii maxtrix 
-pad out numbers or cut at edges to get to the roii x div
-proceding in left to right from top of image
-img goes by (height, width)
-apply thresholding for non mineralized (bone or cartileged or question)
"""
D=np.zeros([fdiv,fdiv,fdiv])
D1=np.zeros([fdiv,fdiv,fdiv])
D2=np.zeros([fdiv,fdiv,fdiv])
D3=np.zeros([fdiv,fdiv,fdiv])
D4=np.zeros([fdiv,fdiv,fdiv])
G=[]
G1=[]
G2=[]
G3=[]
G4=[]
ctr=0
hhindex=math.floor(hh/roii)
wwindex=math.floor(ww/roii)

for k in tqdm(range(0,round(div)-1),position=0,leave=True):
    image_current=bmpmatches[k*resol]
    img = cv2.imread(image_current,0)
    for j in range(0,wwindex):
        for i in range(0,hhindex):    #this is hh
            roi=img[i*roii:(i+1)*roii, j*roii:(j+1)*roii]
            #print(i,j)
            hold=np.average(roi, axis=1)
            value=np.average(hold)
            thresh=55  #This is about 0.7 g/cm3 for value of 55
            if value > thresh:
                band1=[thresh,90]  
                band2=[91,112] 
                band3=[113,134] 
                band4=[135,255] 
                D[k,i,j]=value
                G.append([i,j,k,value])
                if value <= band1[1]:
                    D1[k,i,j]=value
                    G1.append([i,j,k,value])
                if value >= band2[0] and value <= band2[1]:
                    D2[k,i,j]=value
                    G2.append([i,j,k,value])
                if value >= band3[0] and value <= band3[1]:
                    D3[k,i,j]=value
                    G3.append([i,j,k,value])
                if value >= band4[0]: 
                    D4[k,i,j]=value
                    G4.append([i,j,k,value])
              

"""
Reporting
"""

print('img height(px): ',hh,' img width: ',ww)
print('img height(mc): ',h2,' img width: ',w2)

#Here we turn the list of certified bone points (G) into an array
GAR=np.array(G)
X=GAR[:,0]
Y=GAR[:,1]
Z=GAR[:,2]
VV=(GAR[:,3]-57.7)*.0177 + 0.75    #calibration for greyscale to density values
Vold=GAR[:,3]

#to check for no larger mineral values
if len(G3) == 0:
    G3.append([0,0,0,0])
if len(G4) == 0:
    G4.append([0,0,0,0])

GA1=np.array(G1)
GA2=np.array(G2)
GA3=np.array(G3)
GA4=np.array(G4)

VV1=(GA1[:,3]-57.7)*.0177 + 0.75
VV2=(GA2[:,3]-57.7)*.0177 + 0.75
VV3=(GA3[:,3]-57.7)*.0177 + 0.75
VV4=(GA4[:,3]-57.7)*.0177 + 0.75

b1c=[np.min(VV1),np.max(VV1)]    
b2c=[np.min(VV2),np.max(VV2)] 
b3c=[np.min(VV3),np.max(VV3)] 
b4c=[np.min(VV4),np.max(VV4)] 



#we modify the D array
#remember, for D[k,i,j]=value format
D[:,0]=D[:,0]*resol
D[:,1]=D[:,1]*roii
D[:,2]=D[:,2]*roii

#some stats work for histogram work
(mu, sigma)=norm.fit(VV)

sss=time.time()
def centroid_points(arr):
    length = arr.shape[0]
    sum_x = np.sum(arr[:, 0])
    sum_y = np.sum(arr[:, 1])
    sum_z = np.sum(arr[:, 2])
    return sum_x/length, sum_y/length, sum_z/length

cent1=centroid_points(GAR)


#now to calculate distance of each point to centroid.
def getcendist(GARX):
    dcen=np.zeros([len(GARX)])
    for i in tqdm(range(0,len(GARX)),position=0, leave=True):
        pointi=Point(GARX[i,0]*roii,GARX[i,1]*roii,GARX[i,2]*resol)
        poi_dist=pointi.distance(Point(cent1))
        dcen[i]=poi_dist
    return dcen
 
GAR1=np.array(G1)
GAR2=np.array(G2)
GAR3=np.array(G3)
GAR4=np.array(G4)

dcen_total=getcendist(GAR)
dcen1=getcendist(GAR1)
dcen2=getcendist(GAR2)
dcen3=getcendist(GAR3)
dcen4=getcendist(GAR4)

#to scale all the GAR correctly...


#print("centroid distance time -- %s seconds ---" % (time.time() - c1time))

#print("Starting reporting processes")
#reporting parameters and writing to a .txt file
f= open(outpath+imgdir+".txt","w+")
f.write('This is the report file for processing of bone microCT images\n')
f.write('Using HLab Projects code, '+ str(datetime.datetime.now()))
f.write('\n')
f.write('files from: '+imgdir+'\n')
f.write('\n')
f.write('The folowing are parameters used in visualization and calculations \n')
f.write('img height(px): '+str(hh)+'    img width: '+str(ww))
f.write('\n')
f.write('img height(mc): '+str(h2)+'    img width: '+str(w2))
f.write('\n')
f.write('Length per pixel/resolution (mc): '+str(pxc))
f.write('\n')
f.write('Resolution factor is: '+str(resol))
f.write('\n')
f.write('size of 3d scatter array: '+str(fdiv))
f.write('\n')
f.write('Band 1 is between '+str(b1c[0])+' g/cm^3 and '+str(b1c[1])+' g/cm^3 \n')
f.write('Band 2 is between '+str(round(b2c[0],2))+' g/cm^3 and '+str(round(b2c[1],2))+' g/cm^3 \n')
f.write('Band 3 is between '+str(b3c[0])+' g/cm^3 and '+str(b3c[1])+' g/cm^3 \n')
f.write('Band 4 is between '+str(b4c[0])+' g/cm^3 and '+str(b4c[1])+' g/cm^3 \n')
f.close()

#section on creating values to save for later calculation, using folder names
np.save(outpath+imgdir+"fdiv",fdiv)
np.save(outpath+imgdir+"X",X)
np.save(outpath+imgdir+"Y",Y)
np.save(outpath+imgdir+"Z",Z)
np.save(outpath+imgdir+"D",D)
np.save(outpath+imgdir+"GAR",GAR)
np.save(outpath+imgdir+"GAR1",GAR1)
np.save(outpath+imgdir+"GAR2",GAR2)
np.save(outpath+imgdir+"GAR3",GAR3)
np.save(outpath+imgdir+"GAR4",GAR4)
cent11=np.array(cent1)
np.save(outpath+imgdir+"cent11",cent11)
np.save(outpath+imgdir+"dcen_total",dcen_total)
np.save(outpath+imgdir+"dcen1",dcen1)
np.save(outpath+imgdir+"dcen2",dcen2)
np.save(outpath+imgdir+"dcen3",dcen3)
np.save(outpath+imgdir+"dcen4",dcen4)


print("calculations complete")
