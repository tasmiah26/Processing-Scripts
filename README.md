# Processing-Scripts

This repository contains python scripts that automate the fetal MRI processing pipeline. Detailed information about the pipeline and its steps can be found here: [Manual](https://docs.google.com/document/d/1HlpgPguOVPi5-OvLSErXkho5lGzMlBmZahVWLOd30-M/edit?usp=sharing) 


The pipeline is divided into 2 parts: 

1) The first script performs the steps from masking the ROI to automatic segmentation of the inner/outer coritcal plates
2) The second script is for surface and volume extraction. It should be used after performing any manual corrections needed for masking, reconstruction, and segmentation

Both scripts are executable so they can be used by calling the path and using the necessary flags. 

## **Script 1- auto_segmentation**

Current version is  v2.0

The path to the file is:
/neuro/labs/grantlab/research/MRI_processing/tasmiah/script_1/auto_segmentation_v2.0.py

### **Usage**
This script can be used for initial processing of raw MRI scans using the --all flag.

``` 
python3 /neuro/labs/grantlab/research/MRI_processing/tasmiah/script_1/auto_segmentation_v2.0.py --input_fol ${file} --all 
```

--input_fol is a required argument and should be given the path to the folder containing the MRI scans, shown as the variable ${file} above. 

Flag options (replace “--all”) for each step to run individually or to start at a particular step and continue are in the table below.  The --help or -h flag will display information about each flag:

Flag         | Description
------------ | -------------
--masking | creates masks and moves them into a folder labeled "masks" and puts the masked regions inside a folder labeled "brain". Also creates a  "verify" folder to store .png files of the masked ROIs 
--remask or --from_remask| use after manual mask correction, to create a new brain folder with corrected regions
--NUC or --from_NUC | performs non-uniformity correction and puts them inside a folder "nuc"
--QA or --from_QA | creates a folder "Best_Images_Crop" with the highest quality images and quality evaluation is exported to quality_assessment.csv 
--recon or --from_recon |  performs 3 reconstructions
--alignment or --from_alignment | reorients the reconstructed images inside the alignment_temp/ folder for each temp_recon_#/ directory. Files needed for segmentation and surface extraction steps are created and stored inside temp_recon_#/ (recon* files). 
--segment | performs automatic segmentation



## **Script 2- surface_processing**

Current version is  v3.2


The path to the file is:
/neuro/labs/grantlab/research/MRI_processing/tasmiah/script_2/surface_processing_v3.2.py

### **Usage**
``` 
python3 /neuro/labs/grantlab/research/MRI_processing/tasmiah/script_2/surface_processing_v3.2.py --input_fol ${file} --all 
```
Flag options (replace “--all”) for each step to run individually or to start at a particular step and continue. The --help or -h flag will display information about each flag:

Flag         | Description
------------ | -------------
--extract | extracts surfaces from segmentation_to31_final.nii. and transforms surfaces to mni and native size
--registration or --from_registration   | performs surface registration to template (29w, 31w, and adult)
--resample | resample original surface to template and transform the resampled surfaces to 31w template and native space.
--surface_measures | calculates surface area, sulcal depth, and mean curvature. Whole brain measures (l/r area, depth, and absolute mean curvature) will be saved as Area_Depth_aMC.rsl.s5.txt, inside "surfaces" folder
--volume_measures | measures tissue volumes and saves them in Volume_measures.txt, found in "recon_segmentation" folder 
--gyrification_index | calculates left/right/whole gyrification indices and stores them in  GI_info_final.txt, inside "surfaces" folder
