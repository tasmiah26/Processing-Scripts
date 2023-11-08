#!/usr/bin/env python3

#v3.4 - "--all_young" flag added to process fetuses younger than 28.5 GA

#surface extraction after segmentation_to31_final.nii created 
'''
TO DO:
- use Tafoya's flag option for extraction WO subsample
'''
import os
import sys
import math
import argparse
import numpy as np

parser = argparse.ArgumentParser()

parser.add_argument('--input_fol',
    nargs='?',			
    required=True,
    help='path to the data directory')

parser.add_argument('--extraction', 
    dest='extrac', 
	action='store_true',
    help='creates minc files, performs surface extract singularity, and transforms input objects with the corresponding .xfm transform')

parser.add_argument('--all_young', 
    dest='all_young', 
	action='store_true',
    help='creates minc files, performs surface extract singularity, and transforms input objects with the corresponding .xfm transform')

parser.add_argument('--from_registration', 
    dest='reg__', 
	action='store_true',
    help='function to register scan files to 28-31 GW templates')

parser.add_argument('--registration', 
    dest='surf_reg', 
	action='store_true',
    help='function to register scan files to 28-31 GW templates')

parser.add_argument('--resample', 
    dest='resampling',
	action='store_true',
    help='sphere resampling and creating objects with transfometrics for surface extraction')

parser.add_argument('--surface_measures', 
    dest='surf_meas',
	action='store_true',
    help='calculates sulcal depth, surface area, mean curvature')

parser.add_argument('--volume_measures', 
    dest='vol_meas',
	action='store_true',
    help='volume measures stored in txt file')

parser.add_argument('--gyrification_index', 
    dest='GI',
	action='store_true',
    help='calculates and stores the gyrification index values in txt file')

parser.add_argument('--all', 
    dest='all',
	action='store_true',
    help='does all steps')


#surface extraction
def extraction_woSubsample():		
	os.system('nii2mnc '+input_fol+'/recon_segmentation/segmentation_to31_final.nii '+input_fol+'/temp/segmentation_to31_final.mnc;')
	os.system('mincmath -clobber -eq -const 160  '+input_fol+'/temp/segmentation_to31_final.mnc '+input_fol+'/temp/inner_left.mnc;'
				'mincmath -clobber -eq -const 161  '+input_fol+'/temp/segmentation_to31_final.mnc '+input_fol+'/temp/inner_right.mnc;')
		
	os.system('singularity exec docker://fnndsc/pl-fetal-cp-surface-extract extract_cp --adapt_object_mesh 0,100,0,0 --mincmorph-iterations 1 '+input_fol+'/temp/ '+input_fol+'/temp/;')
	
	os.system('adapt_object_mesh '+input_fol+'/temp/inner_left._81920.obj '+input_fol+'/surfaces/lh.smoothwm.to31.obj 0 100 0 0;')
	os.system('adapt_object_mesh '+input_fol+'/temp/inner_right._81920.obj '+input_fol+'/surfaces/rh.smoothwm.to31.obj 0 100 0 0;')

	os.system('/neuro/labs/grantlab/research/HyukJin_MRI/code/obj2asc '+input_fol+'/surfaces/lh.smoothwm.to31.obj '+input_fol+'/surfaces/lh.smoothwm.to31.asc;')
	os.system('/neuro/labs/grantlab/research/HyukJin_MRI/code/obj2asc '+input_fol+'/surfaces/rh.smoothwm.to31.obj '+input_fol+'/surfaces/rh.smoothwm.to31.asc;')

	os.system('transform_objects '+input_fol+'/surfaces/lh.smoothwm.to31.obj '+input_fol+'/recon_segmentation//recon_native.xfm '+input_fol+'/surfaces/lh.smoothwm.native.obj;')
	os.system('transform_objects '+input_fol+'/surfaces/rh.smoothwm.to31.obj '+input_fol+'/recon_segmentation//recon_native.xfm '+input_fol+'/surfaces/rh.smoothwm.native.obj;')
	os.system('transform_objects '+input_fol+'/surfaces/lh.smoothwm.to31.obj /neuro/labs/grantlab/research/HyukJin_MRI/Fetal_template/xfm/template-31toMNI.xfm '+input_fol+'/surfaces/lh.smoothwm.mni.obj;')
	os.system('transform_objects '+input_fol+'/surfaces/rh.smoothwm.to31.obj /neuro/labs/grantlab/research/HyukJin_MRI/Fetal_template/xfm/template-31toMNI.xfm '+input_fol+'/surfaces/rh.smoothwm.mni.obj;')

def extrac():		
	os.system('nii2mnc '+input_fol+'/recon_segmentation/segmentation_to31_final.nii '+input_fol+'/temp/segmentation_to31_final.mnc;')
	os.system('mincmath -clobber -eq -const 160  '+input_fol+'/temp/segmentation_to31_final.mnc '+input_fol+'/temp/inner_left.mnc;'
				'mincmath -clobber -eq -const 161  '+input_fol+'/temp/segmentation_to31_final.mnc '+input_fol+'/temp/inner_right.mnc;')
		
	os.system('singularity exec docker://fnndsc/pl-fetal-cp-surface-extract extract_cp --adapt_object_mesh 0,100,0,0 --mincmorph-iterations 1 --subsample '+input_fol+'/temp/ '+input_fol+'/temp/;')
	
	os.system('adapt_object_mesh '+input_fol+'/temp/inner_left._81920.obj '+input_fol+'/surfaces/lh.smoothwm.to31.obj 0 100 0 0;')
	os.system('adapt_object_mesh '+input_fol+'/temp/inner_right._81920.obj '+input_fol+'/surfaces/rh.smoothwm.to31.obj 0 100 0 0;')

	os.system('/neuro/labs/grantlab/research/HyukJin_MRI/code/obj2asc '+input_fol+'/surfaces/lh.smoothwm.to31.obj '+input_fol+'/surfaces/lh.smoothwm.to31.asc;')
	os.system('/neuro/labs/grantlab/research/HyukJin_MRI/code/obj2asc '+input_fol+'/surfaces/rh.smoothwm.to31.obj '+input_fol+'/surfaces/rh.smoothwm.to31.asc;')

	os.system('transform_objects '+input_fol+'/surfaces/lh.smoothwm.to31.obj '+input_fol+'/recon_segmentation//recon_native.xfm '+input_fol+'/surfaces/lh.smoothwm.native.obj;')
	os.system('transform_objects '+input_fol+'/surfaces/rh.smoothwm.to31.obj '+input_fol+'/recon_segmentation//recon_native.xfm '+input_fol+'/surfaces/rh.smoothwm.native.obj;')
	os.system('transform_objects '+input_fol+'/surfaces/lh.smoothwm.to31.obj /neuro/labs/grantlab/research/HyukJin_MRI/Fetal_template/xfm/template-31toMNI.xfm '+input_fol+'/surfaces/lh.smoothwm.mni.obj;')
	os.system('transform_objects '+input_fol+'/surfaces/rh.smoothwm.to31.obj /neuro/labs/grantlab/research/HyukJin_MRI/Fetal_template/xfm/template-31toMNI.xfm '+input_fol+'/surfaces/rh.smoothwm.mni.obj;')

#surface registration 	
def surf_reg():		
	#template_num=str(sys.arv[2])
	for i in templ_num:
		os.system('bestsurfreg.pl -clobber -min_control_mesh 80 -max_control_mesh 81920 -blur_coef 1.25 -neighbourhood_radius 2.8 -maximum_blur 1.9 /neuro/users/mri.team/fetal_mri/Surface_template/template-'+i+'/bh.smoothwm.mni.obj ./'+input_fol+'/surfaces/lh.smoothwm.mni.obj  ./'+input_fol+'/surfaces/template-'+i+'/lh.smoothwm.toT'+i+'.sm;')
		os.system('bestsurfreg.pl -clobber -min_control_mesh 80 -max_control_mesh 81920 -blur_coef 1.25 -neighbourhood_radius 2.8 -maximum_blur 1.9 /neuro/users/mri.team/fetal_mri/Surface_template/template-'+i+'/bh.smoothwm.mni.obj ./'+input_fol+'/surfaces/rh.smoothwm.mni.obj  ./'+input_fol+'/surfaces//template-'+i+'/rh.smoothwm.toT'+i+'.sm;')

#surface resample
def resampling():
	for i in templ_num:
		os.system('sphere_resample_obj -clobber ./'+input_fol+'/surfaces/lh.smoothwm.mni.obj ./'+input_fol+'/surfaces/template-'+i+'/lh.smoothwm.toT'+i+'.sm ./'+input_fol+'/surfaces/template-'+i+'/lh.smoothwm.mni.rsl'+i+'.obj;')
		os.system('sphere_resample_obj -clobber ./'+input_fol+'/surfaces/rh.smoothwm.mni.obj ./'+input_fol+'/surfaces/template-'+i+'/rh.smoothwm.toT'+i+'.sm ./'+input_fol+'/surfaces/template-'+i+'/rh.smoothwm.mni.rsl'+i+'.obj;')
		
		if i == "29":
			os.system('cp ./'+input_fol+'/surfaces/template-29/lh.smoothwm.mni.rsl29.obj ./'+input_fol+'/surfaces/template-29/lh.smoothwm.mni.rsl.obj;')
			os.system('cp ./'+input_fol+'/surfaces/template-29/rh.smoothwm.mni.rsl29.obj ./'+input_fol+'/surfaces/template-29/rh.smoothwm.mni.rsl.obj;')
		else:
			os.system('sphere_resample_obj -clobber ./'+input_fol+'/surfaces/template-'+i+'/lh.smoothwm.mni.rsl'+i+'.obj /neuro/users/mri.team/fetal_mri/Surface_template/template-'+i+'/bh.smoothwm.T29.sm ./'+input_fol+'/surfaces/template-'+i+'/lh.smoothwm.mni.rsl.obj;'
			'sphere_resample_obj -clobber ./'+input_fol+'/surfaces/template-'+i+'/rh.smoothwm.mni.rsl'+i+'.obj /neuro/users/mri.team/fetal_mri/Surface_template/template-'+i+'/bh.smoothwm.T29.sm ./'+input_fol+'/surfaces/template-'+i+'/rh.smoothwm.mni.rsl.obj;')

		os.system('transform_objects ./'+input_fol+'/surfaces/template-'+i+'/lh.smoothwm.mni.rsl.obj /neuro/users/mri.team/fetal_mri/Surface_template/xfm/template-31toMNI_inv.xfm  ./'+input_fol+'/surfaces/template-'+i+'/lh.smoothwm.to31.rsl.obj;')
		os.system('transform_objects ./'+input_fol+'/surfaces/template-'+i+'/rh.smoothwm.mni.rsl.obj /neuro/users/mri.team/fetal_mri/Surface_template/xfm/template-31toMNI_inv.xfm  ./'+input_fol+'/surfaces/template-'+i+'/rh.smoothwm.to31.rsl.obj;')
		os.system('transform_objects '+input_fol+'/surfaces/template-'+i+'/lh.smoothwm.to31.rsl.obj ./'+input_fol+'/recon_segmentation/recon_native.xfm ./'+input_fol+'/surfaces/template-'+i+'/lh.smoothwm.native.rsl.obj;')
		os.system('transform_objects '+input_fol+'/surfaces/template-'+i+'/rh.smoothwm.to31.rsl.obj ./'+input_fol+'/recon_segmentation/recon_native.xfm ./'+input_fol+'/surfaces/template-'+i+'/rh.smoothwm.native.rsl.obj;')

#surface measures
def surf_meas():
	for i in templ_num:
#sulcal depth 			
		os.system('python3 /neuro/labs/grantlab/research/HyukJin_MRI/code/ADT/ADT_white_vFetal_final_rsl.py '+input_fol+'/surfaces/template-'+i+'/;')
		os.system('depth_potential -smooth 10 '+input_fol+'/surfaces/template-'+i+'/lh.smoothwm.native.rsl.depth '+input_fol+'/surfaces/template-'+i+'/lh.smoothwm.mni.rsl.obj  '+input_fol+'/surfaces/template-'+i+'/lh.smoothwm.native.rsl.depth.s10;')
		os.system('depth_potential -smooth 10 '+input_fol+'/surfaces/template-'+i+'/rh.smoothwm.native.rsl.depth '+input_fol+'/surfaces/template-'+i+'/rh.smoothwm.mni.rsl.obj  '+input_fol+'/surfaces/template-'+i+'/rh.smoothwm.native.rsl.depth.s10;')	
		os.system('depth_potential -smooth 20 '+input_fol+'/surfaces/template-'+i+'/lh.smoothwm.native.rsl.depth '+input_fol+'/surfaces/template-'+i+'/lh.smoothwm.mni.rsl.obj  '+input_fol+'/surfaces/template-'+i+'/lh.smoothwm.native.rsl.depth.s20;')
		os.system('depth_potential -smooth 20 '+input_fol+'/surfaces/template-'+i+'/rh.smoothwm.native.rsl.depth '+input_fol+'/surfaces/template-'+i+'/rh.smoothwm.mni.rsl.obj  '+input_fol+'/surfaces/template-'+i+'/rh.smoothwm.native.rsl.depth.s20;')	

# surface area
		os.system('depth_potential -area_voronoi '+input_fol+'/surfaces/template-'+i+'/lh.smoothwm.native.rsl.obj '+input_fol+'/surfaces/template-'+i+'/lh.smoothwm.native.rsl.area;')
		os.system('depth_potential -area_voronoi '+input_fol+'/surfaces/template-'+i+'/rh.smoothwm.native.rsl.obj '+input_fol+'/surfaces/template-'+i+'/rh.smoothwm.native.rsl.area;')
		os.system('depth_potential -smooth 10 '+input_fol+'/surfaces/template-'+i+'/lh.smoothwm.native.rsl.area '+input_fol+'/surfaces/template-'+i+'/lh.smoothwm.mni.rsl.obj '+input_fol+'/surfaces/template-'+i+'/lh.smoothwm.native.rsl.area.s10;')
		os.system('depth_potential -smooth 10 '+input_fol+'/surfaces/template-'+i+'/rh.smoothwm.native.rsl.area '+input_fol+'/surfaces/template-'+i+'/rh.smoothwm.mni.rsl.obj '+input_fol+'/surfaces/template-'+i+'/rh.smoothwm.native.rsl.area.s10;')
		os.system('depth_potential -smooth 20 '+input_fol+'/surfaces/template-'+i+'/lh.smoothwm.native.rsl.area '+input_fol+'/surfaces/template-'+i+'/lh.smoothwm.mni.rsl.obj '+input_fol+'/surfaces/template-'+i+'/lh.smoothwm.native.rsl.area.s20;')
		os.system('depth_potential -smooth 20 '+input_fol+'/surfaces/template-'+i+'/rh.smoothwm.native.rsl.area '+input_fol+'/surfaces/template-'+i+'/rh.smoothwm.mni.rsl.obj '+input_fol+'/surfaces/template-'+i+'/rh.smoothwm.native.rsl.area.s20;')

# mean curvature
		os.system('depth_potential -mean_curvature '+input_fol+'/surfaces/template-'+i+'/lh.smoothwm.mni.rsl.obj  '+input_fol+'/surfaces/template-'+i+'/lh.smoothwm.mni.rsl.mc;')
		os.system('depth_potential -mean_curvature '+input_fol+'/surfaces/template-'+i+'/rh.smoothwm.mni.rsl.obj  '+input_fol+'/surfaces/template-'+i+'/rh.smoothwm.mni.rsl.mc;')
		os.system('depth_potential -smooth 10 '+input_fol+'/surfaces/template-'+i+'/lh.smoothwm.mni.rsl.mc '+input_fol+'/surfaces/template-'+i+'/lh.smoothwm.mni.rsl.obj '+input_fol+'/surfaces/template-'+i+'/lh.smoothwm.mni.rsl.mc.s10;')
		os.system('depth_potential -smooth 10 '+input_fol+'/surfaces/template-'+i+'/rh.smoothwm.mni.rsl.mc '+input_fol+'/surfaces/template-'+i+'/rh.smoothwm.mni.rsl.obj '+input_fol+'/surfaces/template-'+i+'/rh.smoothwm.mni.rsl.mc.s10;')
		os.system('depth_potential -smooth 20 '+input_fol+'/surfaces/template-'+i+'/lh.smoothwm.mni.rsl.mc '+input_fol+'/surfaces/template-'+i+'/lh.smoothwm.mni.rsl.obj '+input_fol+'/surfaces/template-'+i+'/lh.smoothwm.mni.rsl.mc.s20;')
		os.system('depth_potential -smooth 20 '+input_fol+'/surfaces/template-'+i+'/rh.smoothwm.mni.rsl.mc '+input_fol+'/surfaces/template-'+i+'/rh.smoothwm.mni.rsl.obj '+input_fol+'/surfaces/template-'+i+'/rh.smoothwm.mni.rsl.mc.s20;')

#s5 calculations inside /surfaces folder
	os.system('python3 /neuro/labs/grantlab/research/HyukJin_MRI/code/ADT/ADT_white_vFetal_final.py '+input_fol+'/surfaces/;')
	os.system('depth_potential -smooth 5 '+input_fol+'/surfaces/lh.smoothwm.native.depth '+input_fol+'/surfaces/lh.smoothwm.mni.obj  '+input_fol+'/surfaces/lh.smoothwm.native.depth.s5;')
	os.system('depth_potential -smooth 5 '+input_fol+'/surfaces/rh.smoothwm.native.depth '+input_fol+'/surfaces/rh.smoothwm.mni.obj  '+input_fol+'/surfaces/rh.smoothwm.native.depth.s5;')

	os.system('depth_potential -area_voronoi '+input_fol+'/surfaces/lh.smoothwm.native.obj '+input_fol+'/surfaces/lh.smoothwm.native.area;')
	os.system('depth_potential -area_voronoi '+input_fol+'/surfaces/rh.smoothwm.native.obj '+input_fol+'/surfaces/rh.smoothwm.native.area;')
	os.system('depth_potential -smooth 5 '+input_fol+'/surfaces/lh.smoothwm.native.area '+input_fol+'/surfaces/lh.smoothwm.mni.obj  '+input_fol+'/surfaces/lh.smoothwm.native.area.s5;')
	os.system('depth_potential -smooth 5 '+input_fol+'/surfaces/rh.smoothwm.native.area '+input_fol+'/surfaces/rh.smoothwm.mni.obj  '+input_fol+'/surfaces/rh.smoothwm.native.area.s5;')

	os.system('depth_potential -mean_curvature '+input_fol+'/surfaces/lh.smoothwm.mni.obj  '+input_fol+'/surfaces/lh.smoothwm.mni.mc;')
	os.system('depth_potential -mean_curvature '+input_fol+'/surfaces/rh.smoothwm.mni.obj  '+input_fol+'/surfaces/rh.smoothwm.mni.mc;')
	os.system('depth_potential -smooth 5 '+input_fol+'/surfaces/lh.smoothwm.mni.mc '+input_fol+'/surfaces/lh.smoothwm.mni.obj '+input_fol+'/surfaces/lh.smoothwm.mni.mc.s5;')
	os.system('depth_potential -smooth 5 '+input_fol+'/surfaces/rh.smoothwm.mni.mc '+input_fol+'/surfaces/rh.smoothwm.mni.obj '+input_fol+'/surfaces/rh.smoothwm.mni.mc.s5;')

# whole brain depth, area, and curvature -> txt 
	os.system('echo `python3 /neuro/labs/grantlab/research/HyukJin_MRI/code/vertex_mean_sum.py -sum '+input_fol+'/surfaces/lh.smoothwm.native.area.s5` '
		'`python3 /neuro/labs/grantlab/research/HyukJin_MRI/code/vertex_mean_sum.py -sum '+input_fol+'/surfaces/rh.smoothwm.native.area.s5` '
		'`python3 /neuro/labs/grantlab/research/HyukJin_MRI/code/vertex_mean_sum.py -mean '+input_fol+'/surfaces/lh.smoothwm.native.depth.s5` '
		'`python3 /neuro/labs/grantlab/research/HyukJin_MRI/code/vertex_mean_sum.py -mean '+input_fol+'/surfaces/rh.smoothwm.native.depth.s5` '
		'`python3 /neuro/labs/grantlab/research/HyukJin_MRI/code/vertex_mean_sum.py -ab_mean '+input_fol+'/surfaces/lh.smoothwm.mni.mc.s5` '
		'`python3 /neuro/labs/grantlab/research/HyukJin_MRI/code/vertex_mean_sum.py -ab_mean '+input_fol+'/surfaces/rh.smoothwm.mni.mc.s5` >  '+input_fol+'/surfaces/Area_Depth_aMC.txt;')

# verify surface measures 
	os.system('mkdir '+input_fol+'/surfaces/verify')
	os.system('cp '+input_fol+'/surfaces/*.mni.obj* '+input_fol+'/surfaces/*.s5* '+input_fol+'/surfaces/verify;')
	os.system('apptainer exec docker://fnndsc/pl-surfigures:latest surfigures --range .area.s5:0.0:1.0,.depth.s5:0.0:20.0,.mc.s5:-0.5:0.5 -s .s5 -o s5.png surfaces/verify surfaces/verify')

#volume measures
def vol_meas():			
	seg_final= input_fol+'/recon_segmentation/segmentation_to31_final.nii'
	final_nii= nib.load(seg_final)
	img_data= final_nii.get_fdata()

	vox_inRange_1 = ((img_data > 159) & (img_data < 161)).sum()
	vox_inRange_2 = ((img_data > 160) & (img_data < 162)).sum()
	vox_inRange_3 = ((img_data > 41) & (img_data < 43)).sum()
	vox_inRange_4 = ((img_data > 0) & (img_data < 2)).sum()

	print("num voxels:", vox_inRange_1, vox_inRange_2, vox_inRange_3, vox_inRange_4)

	voxel_sizes = np.abs(final_nii.header.get_zooms())
	print("vox size:", voxel_sizes)
    
	vol1 = vox_inRange_1 * np.prod(voxel_sizes)
	vol2 = vox_inRange_2 * np.prod(voxel_sizes)
	vol3 = vox_inRange_3 * np.prod(voxel_sizes)
	vol4 = vox_inRange_4 * np.prod(voxel_sizes)
	
	xfm_path = input_fol+'/recon_segmentation/recon_native.xfm'
			
	os.system('xfm2param '+input_fol+'/recon_segmentation/recon_native.xfm | grep scale | cut -d" " -f11 > '+input_fol+'/temp/temp2.txt;')
	os.system('xfm2param '+input_fol+'/recon_segmentation/recon_native.xfm | grep scale | cut -d" " -f15 >> '+input_fol+'/temp/temp2.txt;')
	os.system('xfm2param '+input_fol+'/recon_segmentation/recon_native.xfm | grep scale | cut -d" " -f19 >> '+input_fol+'/temp/temp2.txt;')
	
	with open(input_fol+'/temp/temp2.txt') as temp2:
		scale_temp = [0,0,0]
		scale_temp=temp2.readlines()

	scale_av = float(scale_temp[0]) * float(scale_temp[1]) * float(scale_temp[2])
	print(scale_av)

	v1xScale = (vol1 * scale_av)
	v2xScale = (vol2 * scale_av)
	v3xScale = (vol3 * scale_av)
	v4xScale = (vol4 * scale_av)

	with open(input_fol+'/recon_segmentation/Volume_measures.txt', "w") as file:
		file.write(str(v1xScale) + "\n" + str(v2xScale)+ "\n" + str(v3xScale) + "\n" +str(v4xScale)  )

#gyrification index
def GI():
	os.system('mkdir ./'+input_fol+'/temp/hull/;')
	os.system('mri_distance_transform ./'+input_fol+'/recon_segmentation/segmentation_to31_final.nii 160 10 1 ./'+input_fol+'/temp/lh.dist10.nii;'
			'mri_threshold ./'+input_fol+'/temp/lh.dist10.nii 10 ./'+input_fol+'/temp/lh.dist10.mask.nii;'
			'mri_distance_transform ./'+input_fol+'/temp/lh.dist10.mask.nii 0 10 2 ./'+input_fol+'/temp/lh.dist10.mask.dist10.nii;'
			'fslmaths ./'+input_fol+'/temp/lh.dist10.mask.dist10.nii -mul -1 ./'+input_fol+'/temp/lh.dist10.mask.dist10.inv.nii;'
			'gunzip ./'+input_fol+'/temp/lh.dist10.mask.dist10.inv.nii.gz;'
			'mri_threshold ./'+input_fol+'/temp/lh.dist10.mask.dist10.inv.nii 10 ./'+input_fol+'/temp/lh.closing.nii;')
			
	os.system('mri_distance_transform ./'+input_fol+'/recon_segmentation/segmentation_to31_final.nii 161 10 1 ./'+input_fol+'/temp/rh.dist10.nii;'
			'mri_threshold ./'+input_fol+'/temp/rh.dist10.nii 10 ./'+input_fol+'/temp/rh.dist10.mask.nii;'
			'mri_distance_transform ./'+input_fol+'/temp/rh.dist10.mask.nii 0 10 2 ./'+input_fol+'/temp/rh.dist10.mask.dist10.nii;'
			'fslmaths ./'+input_fol+'/temp/rh.dist10.mask.dist10.nii -mul -1 ./'+input_fol+'/temp/rh.dist10.mask.dist10.inv.nii;'
			'gunzip ./'+input_fol+'/temp/rh.dist10.mask.dist10.inv.nii.gz;'
			'mri_threshold ./'+input_fol+'/temp/rh.dist10.mask.dist10.inv.nii 10 ./'+input_fol+'/temp/rh.closing.nii;')

	os.system('nii2mnc ./'+input_fol+'/temp/lh.closing.nii ./'+input_fol+'/temp/lh.closing.mnc;'
				'nii2mnc ./'+input_fol+'/temp/rh.closing.nii ./'+input_fol+'/temp/rh.closing.mnc;')
				
	os.system('mincmath -clobber -gt -const 0 ./'+input_fol+'/temp/lh.closing.mnc ./'+input_fol+'/temp/hull/closing_bin_left.mnc;'
				'mincmath -clobber -gt -const 0 ./'+input_fol+'/temp/rh.closing.mnc ./'+input_fol+'/temp/hull/closing_bin_right.mnc;')

				
	os.system('singularity exec docker://fnndsc/pl-fetal-cp-surface-extract extract_cp --adapt_object_mesh 0,100,0,0 '+input_fol+'/temp/hull/ '+input_fol+'/temp/hull/;')
	
	os.system('cp '+input_fol+'/temp/hull/closing_bin_left._81920.obj '+input_fol+'/temp/lh.hull.to31.obj;'
				'cp '+input_fol+'/temp/hull/closing_bin_right._81920.obj '+input_fol+'/temp/rh.hull.to31.obj;')

	os.system('measure_surface_area ./'+input_fol+'/surfaces/lh.smoothwm.to31.obj > ./'+input_fol+'/temp/GI_inf2.txt;'
			'measure_surface_area ./'+input_fol+'/temp/lh.hull.to31.obj >> ./'+input_fol+'/temp/GI_inf2.txt;'
			'measure_surface_area ./'+input_fol+'/surfaces/rh.smoothwm.to31.obj >> ./'+input_fol+'/temp/GI_inf2.txt;'
			'measure_surface_area ./'+input_fol+'/temp/rh.hull.to31.obj >> ./'+input_fol+'/temp/GI_inf2.txt;')
			
	with open(input_fol+'/temp/GI_inf2.txt') as GI_temp:
		temp=GI_temp.read()
		temp=temp.split('Area: '); 
		
		bh_surf= float(temp[1])+float(temp[3])
		bh_hull= float(temp[2])+float(temp[4])
		lh_GI=float(temp[1])/float(temp[2])
		rh_GI=float(temp[3])/float(temp[4])
		
		GI_temp.close

		bh_GI=str(bh_surf/bh_hull)
		lh_GI=str(lh_GI)
		rh_GI=str(rh_GI)
			
	os.system('echo '+lh_GI+' '+rh_GI+' '+bh_GI+' > ./'+input_fol+'/surfaces/GI_info_final.txt;')
	

def main():
	args = parser.parse_args()
	print(args)

	global input_fol
	input_fol = args.input_fol		## input_fol should be: ./$file

	if not os.path.exists(input_fol+'/surfaces/'):		#/ use os.mkdir()?
		os.system('mkdir '+input_fol+'/surfaces/')		
	if not os.path.exists(input_fol+'/temp/'):
		os.system('mkdir '+input_fol+'/temp/')

	global templ_num
	global i 
	templ_num=("29","31","adult")		##default to template 29
	for i in templ_num:
		os.makedirs(input_fol+'/surfaces/template-'+i+'/', exist_ok=True)
	
	extraction = args.extrac
	from_reg = args.reg__
	registration = args.surf_reg
	resample = args.resampling
	S_meas = args.surf_meas
	V_meas = args.vol_meas
	gi = args.GI
	allSteps = args.all
	all_young = args.all_young


	rnum=("1","2","3")
	#creates recon_native.xfm if not already inside appropriate folder(s) 
	for r in rnum:
		if os.path.exists(input_fol+'/temp_recon_'+r) and not os.path.exists(input_fol+'/temp_recon_'+r+'/recon_native.xfm'):
			os.system('convert_xfm -omat '+input_fol+'/temp_recon_'+r+'/recon_to31_inv.xfm -inverse '+input_fol+'/temp_recon_'+r+'/recon_to31.xfm;')
			os.system('echo `avscale '+input_fol+'/temp_recon_'+r+'/recon_to31_inv.xfm | grep Scales` > '+input_fol+'/temp/temp.txt;')
			scales=open(input_fol+'/temp/temp.txt', encoding='utf-8')
			scales=scales.read()
			os.system('param2xfm -clobber -scales '+scales[16:-1]+' '+input_fol+'/temp_recon_'+r+'/recon_native.xfm;')

		if os.path.exists(input_fol+'/temp_recon_'+r+'/segmentation_to31_final.nii') and not os.path.exists(input_fol+'/recon_segmentation/segmentation_to31_final.nii'):
			##use os.rename instead of mv?
			os.system('mv -T '+input_fol+'/temp_recon_'+r+' '+input_fol+'/recon_segmentation')
			print('segmentation_to31_final.nii file found in temp_recon_'+r+' and folder renamed to "recon_segmentation"')		

	if not os.path.exists(input_fol+'/recon_segmentation/recon_native.xfm'):
		os.system('convert_xfm -omat '+input_fol+'/recon_segmentation/recon_to31_inv.xfm -inverse '+input_fol+'/recon_segmentation/recon_to31.xfm;')
		os.system('echo `avscale '+input_fol+'/recon_segmentation/recon_to31_inv.xfm | grep Scales` > '+input_fol+'/temp/temp.txt;')
		scales=open(input_fol+'/temp/temp.txt', encoding='utf-8')
		scales=scales.read()
		os.system('param2xfm -clobber -scales '+scales[16:-1]+' '+input_fol+'/recon_segmentation/recon_native.xfm;')

	if extraction == True:
		extrac()
		print('surface extraction step complete')
	if from_reg == True:
		surf_reg()
		resampling()
		surf_meas()
		vol_meas()
		GI()
	if registration == True:
		surf_reg()
	if resample == True:
		resampling()
	if S_meas == True:
		surf_meas() 
	if V_meas == True:
		vol_meas()
		
	if gi == True:
		GI()

	if allSteps == True:
		extrac()
		surf_reg()
		resampling()
		surf_meas()
		vol_meas()
		GI() 

if __name__ == '__main__':
    main()
