# RegenBoneAnalysis

This is the readme file of the Regenerated Bone Analysis repository.

Contained in this repository are scripts developed for improved analysis of regenerated bone, created and utilized in submitted manuscript:  
Hoffseth,K.F., Simkin, J., Busse, E., Stewart, K., Watt, J., Hargrove, A., Sammarco, M. "A new approach to analyzing regenerated bone quality in the mouse digit amputation model using semi-automatic processing of microCT data" BONE. Submitted August 2020.

The manuscript details a new approach to analyze the regenerated bone tissue, based on semi-automatic methods combining image processing, solid modeling, and numerical calculations to analyze bone tissue at a more granular level using µCT image data from a mouse digit model of bone regeneration.  These scripts reflect the use of those methds, excepting select modeling done through MeshMixer.  

Each script is structured to be run as a stand-alone script.  The file "mdbands.py", should be run first, as all other scripts depend on the data it generates, stored in numpy format.  The user must choose their input and out directories for correct import of microCT images, and export of data, respectively.  

Note:  The scripts assume the P3 bone has been segmented properly prior to import into "mdbands.py". 

Main Dependencies:
  -python > 3.7
  -numpy > 1.18.1
  -opencv > 4.1
  -matplotlib > 3.1.3
  -scikit-image > 0.16.2
  -skan > ?
  -scipy > 1.4.1

Feel free to utilize the code, we simply ask you cite the work as follows:
Hoffseth,K.F., Simkin, J., Busse, E., Stewart, K., Watt, J., Hargrove, A., Sammarco, M. "A new approach to analyzing regenerated bone quality in the mouse digit amputation model using semi-automatic processing of microCT data" BONE. Submitted August 2020.

Script Descriptions:

Under construction
