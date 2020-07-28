# -*- coding: utf-8 -*-
"""
Created on Tue Jul 28 13:21:10 2020
@author: KFH
-script to skeletonize interior void volume.
*Note: depends on data from previous scripts
"""
import numpy as np
from skimage.morphology import skeletonize, skeletonize_3d
from skan import skeleton_to_csgraph
from skan import Skeleton, summarize

ipath="your_ct_img_path"      #this is the path to your microCT image files
imgdir="your_ct_img_folder"
outpath="your_file_outpath"       #this is the path to where your output files go

"""
Processing
"""
Mint=np.load(outpath+imgdir+"Mint_r"+".npy")
Mint_uint8=Mint.astype('uint8')
skel3=skeletonize_3d(Mint_uint8)   
"""
Note:
#Mint is the 3d binary output of the internal void shape
#skel stands for skeleton
#needs a 0 or 1 binary array
"""

#section using skan package, see reference in manuscript for more information
pixel_graph, coordinates, degrees = skeleton_to_csgraph(skel3)
branch_data= summarize(Skeleton(skel3))
branch_data.head()

tot_blength=np.sum(branch_data["branch-distance"])

#output
np.save(outpath+imgdir+"3Dslen",tot_blength)
np.save(outpath+imgdir+"skel3_",skel3)

print("script complete")