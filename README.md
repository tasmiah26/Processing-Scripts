# Fetal MRI Processing Scripts

This repository contains python scripts that automate the fetal MRI processing pipeline. Detailed information about the pipeline and its steps can be found here: [Manual](https://docs.google.com/document/d/1HlpgPguOVPi5-OvLSErXkho5lGzMlBmZahVWLOd30-M/edit?usp=sharing) 


The pipeline is divided into 2 parts: 

1) The first script performs the steps from masking the ROI (brain) to automatic segmentation of the inner/outer coritcal plates
2) The second script is for surface and volume extraction. It should be used after performing any manual corrections needed for masking, reconstruction, and segmentation

The scripts can only be accessed and used on an FNNDSC machine after setting up the environment. Both scripts are executable so they can be used by calling the path and using the necessary flags. 

## **Script 1- auto_segmentation**

Current version is  v2.0

The path to the file is:
``/neuro/labs/grantlab/research/MRI_processing/tasmiah/script_1/auto_segmentation_v2.0.py``

### **Usage**
This script can be used for initial processing of raw MRI scans using the --all flag.

``` 
python3 /neuro/labs/grantlab/research/MRI_processing/tasmiah/script_1/auto_segmentation_v2.0.py --input_fol ${file} --all 
```

`--input_fol` is a required argument and should be given the path to the folder containing the MRI scans, shown as the variable "file" above. 

Flag options (replace “`--all`”) for each step to run individually or to start at a particular step and continue are in the table below.  The `--help` or `-h` flag will display information about each flag:

Flag         | Description
------------ | -------------
--masking | creates masks and moves them into a folder labeled "masks" and puts the masked regions inside a folder labeled "brain" [(source code)](https://github.com/sofia-urosa/brain-masking). Also creates a  "verify" folder to store .png files of the masked ROIs 
--remask or --from_remask| use after manual mask correction, to create a new brain folder with corrected regions
--NUC or --from_NUC | performs non-uniformity correction and puts them inside a folder "nuc" [(source code)](https://github.com/FNNDSC/pl-ANTs_N4BiasFieldCorrection) 
--QA or --from_QA | creates a folder "Best_Images_Crop" with the highest quality images and the quality evaluation is exported to quality_assessment.csv [(source code)](https://github.com/FNNDSC/pl-fetal-brain-assessment)
--recon or --from_recon |  performs 3 reconstructions using the top 3 targets with the highest quality scores
--alignment or --from_alignment | reorients the reconstructed images inside the temp_recon_#/alignment_temp folder. Files necessary for segmentation and surface extraction steps are created and stored inside temp_recon_#/ (recon* files)
--segment | performs automatic segmentation **

**automatic segmentation has another flag option that can be used if the default resulted poorly: ```--segmentation_WO_att```. This will overwrite the recon_to31_nuc_deep_agg.nii file, so be sure to save the previous segmentation by changing the name


## **Script 2- surface_processing**

Current version is  v3.3

The path to the file is:
`/neuro/labs/grantlab/research/MRI_processing/tasmiah/script_2/surface_processing_v3.2.py
`
### **Usage**
``` 
python3 /neuro/labs/grantlab/research/MRI_processing/tasmiah/script_2/surface_processing_v3.2.py --input_fol ${file} --all 
```
Flag options for each step to run individually or to start at a particular step and continue are shown in the table below. The `--help `or `-h` flag will display the information about each flag:

Flag         | Description
------------ | -------------
--extract | extracts surfaces from segmentation_to31_final.nii and transforms surfaces to mni and native size [(source code)](https://github.com/FNNDSC/pl-fetal-surface-extract)
--registration or --from_registration   | performs surface registration to templates (29w, 31w, and adult)
--resample | resamples original surface to template, and transforms the resampled surfaces to 31w template and native space
--surface_measures | calculates surface area, sulcal depth, and mean curvature. Whole brain measure will be saved as Area_Depth_aMC.rsl.s5.txt (outputs in the following order: left surface area, right surface area, left sulcal depth, right sulcal depth, left absolute mean curvature, and right absolute mean curvature), inside "surfaces" folder. Also creates a verification s5.png file that shows the extracted surface in various heat maps for visualization [(source code)](https://github.com/FNNDSC/pl-surfigures)
--volume_measures | measures tissue volumes and saves them in Volume_measures.txt (output order: left inner volume, right inner volume, left CP volume, and right CP volume), found in "recon_segmentation" folder 
--gyrification_index | calculates left/right/whole gyrification indices and stores them in  GI_info_final.txt, inside "surfaces" folder
