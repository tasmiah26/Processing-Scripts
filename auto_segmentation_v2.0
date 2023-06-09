#!/usr/bin/env python3

#folder: /neuro/labs/grantlab/research/MRI_processing/tasmiah/script_1/
#v 2.0- updated folder structure 

### add error messages

    
import numpy as np
import math
from numpy import zeros
import nibabel as nib
import matplotlib.pyplot as plt
import os
import csv
import sys
import glob
import shutil
import argparse

parser = argparse.ArgumentParser()

parser.add_argument('--input_fol',
    nargs='?',			
    required=True,
    help='relative path to the data directory ($file variable)')

parser.add_argument('--masking', 
    dest='masks', 
	action='store_true',
    help='creates masks of the raw scans and moves the mask files into a folder')

parser.add_argument('--remask', 
    dest='remask', 
	action='store_true',
    help='creates brain folder and extracts brain region using manually corrected masks')

parser.add_argument('--NUC', 
    dest='NUC', 
	action='store_true',
    help='performs non uniformity corrrection')

parser.add_argument('--QA', 
    dest='QA',
	action='store_true',
    help='creates quality assessment .csv')

parser.add_argument('--recon', 
    dest='recon',
	action='store_true',
    help='performs 3 reconstructions using different targets')  ##make # of recons a flag option; default to 3 

parser.add_argument('--alignment', 
    dest='align',
    # nargs='+',
    # default='1','2','3',
	action='store_true',
    help='aligns the reconstructed images')

parser.add_argument('--segment', 
    dest='auto_seg',
	action='store_true',
    help='automatically segments the reconstructed images')

parser.add_argument('--all', 
    dest='all',
	action='store_true',
    help='does all steps from masking')

parser.add_argument('--from_remasking', 
    dest='remask__',
	action='store_true',
    help='extracts corrected brain region and does following steps')

parser.add_argument('--from_NUC', 
    dest='nuc__',
	action='store_true',
    help='NUC to auto segmentation and resize')

parser.add_argument('--from_QA', 
    dest='qa__',
	action='store_true',
    help='QA to resize')

parser.add_argument('--from_recon', 
    dest='recon__',
	action='store_true',
    help='recon to resize')

parser.add_argument('--from_alignment', 
    dest='align__',
	action='store_true',
    help='alignment to resize')

def verify():
    img_list = np.asarray(sorted(glob.glob(input_fol+'/verify/*_brain.nii*')))
    def auto_crop_image(input_name, output_name, reserve):
        nim = nib.load(input_name)
        image = nim.get_data()
        if np.mean(image) == 0:
            print(input_name,'\t Passed')
            return 0
        # else:
        #     print(input_name, '\t Worked')
        image = np.pad(image, [(50,50),(50,50),(16,16)], 'constant')
        X, Y, Z = image.shape[:3]

        # Detect the bounding box of the foreground
        idx = np.nonzero(image > 0)
        x1, x2 = idx[0].min() - reserve[0,0], idx[0].max() + reserve[0,1] + 1
        y1, y2 = idx[1].min() - reserve[1,0], idx[1].max() + reserve[1,1] + 1
        z1, z2 = idx[2].min() - reserve[2,0], idx[2].max() + reserve[2,1] + 1
        # print('Bounding box')
        # print(input_name+'\t'+str([x2-x1, y2-y1, z2-z1]))
        # return [x2-x1, y2-y1, z2-z1]
        # print('  bottom-left corner = ({},{},{})'.format(x1, y1, z1))
        # print('  top-right corner = ({},{},{})'.format(x2, y2, z2))

        # Crop the image
        image = image[x1:x2, y1:y2, z1:z2]

        # Update the affine matrix
        affine = nim.affine
        affine[:3, 3] = np.dot(affine, np.array([x1, y1, z1, 1]))[:3]
        nim2 = nib.Nifti1Image(image, affine)
        nib.save(nim2, output_name)
        return image

    for i in range(len(img_list)):
        f,axarr = plt.subplots(1,6)#,figsize=(len(in_img_list),9))
        f.patch.set_facecolor('k')
        #imsize = np.zeros([len(in_img_list),3])
        img = auto_crop_image(img_list[i], img_list[i].replace('.nii.gz','').replace('.nii','')+'_crop.nii.gz', np.array([[0,0],[0,0],[0,0]]))
        if isinstance(img,(list, tuple, np.ndarray)) == False:
            continue
        
        img = nib.load(img_list[i].replace('.nii.gz','').replace('.nii','')+'_crop.nii.gz').get_data()
        hdr = nib.load(img_list[i].replace('.nii.gz','').replace('.nii','')+'_crop.nii.gz').header
        axarr[0].imshow(np.rot90(img[:,:,np.int_(img.shape[-1]*0.3)]),cmap='gray')
        axarr[0].axis('off')
        axarr[0].set_title(str(img_list[i]),size=5,color='white')

        axarr[1].imshow(np.rot90(img[:,:,np.int_(img.shape[-1]*0.4)]),cmap='gray')
        axarr[1].axis('off')

        axarr[2].imshow(np.rot90(img[:,:,np.int_(img.shape[-1]*0.5)]),cmap='gray')
        axarr[2].axis('off')

        axarr[3].imshow(np.rot90(img[:,:,np.int_(img.shape[-1]*0.6)]),cmap='gray')
        axarr[3].axis('off')

        axarr[4].imshow(np.rot90(img[:,np.int_(img.shape[-2]*0.5),:]),cmap='gray',aspect=str(hdr['pixdim'][3]/hdr['pixdim'][2]))
        axarr[4].axis('off')

        axarr[5].imshow(np.rot90(img[np.int_(img.shape[0]*0.5),:,:]),cmap='gray',aspect=str(hdr['pixdim'][3]/hdr['pixdim'][1]))
        axarr[5].axis('off')

        f.subplots_adjust(wspace=0, hspace=0)
        plt.savefig(img_list[i].replace('.nii.gz','').replace('.nii','')+'_crop_verify.png', facecolor=f.get_facecolor(), pad_inches=0, dpi=300)
        plt.close()
    return 0

    def _interslice_inorm(img_array):
        for i in range(img_array.shape[-1]):
            if i==0 :
                n_mean = np.mean(img_array[:,:,i+1][img_array[:,:,i+1]>0])
            elif i==img_array.shape[-1]-1:
                n_mean = np.mean(img_array[:,:,i-1][img_array[:,:,i-1]>0])
            else:
                n_mean = np.mean(img_array[:,:,[i-1,i+1]][img_array[:,:,[i-1,i+1]]>0])
            loc = np.where(img_array[:,:,i])
            img_array[:,:,i][loc]=img_array[:,:,i][loc]-np.mean(img_array[:,:,i][loc])+n_mean
        return img_array

    def _3dN4(img_array):
        import SimpleITK as sitk
        sitk_img = sitk.GetImageFromArray(img_array)
        maskImage = sitk.GetImageFromArray((img_array>0).astype(int))
        corrector = sitk.N4BiasFieldCorrectionImageFilter()
        sitk_img = sitk.Cast(sitk_img, sitk.sitkFloat32)
        maskImage = sitk.Cast(maskImage, sitk.sitkUInt8)
        output = corrector.Execute( sitk_img, maskImage )
        img_array = sitk.GetArrayFromImage(output)
        return img_array

    def _2dN4(img_array):
        import SimpleITK as sitk
        result = np.zeros(np.shape(img_array))
        for i in range(result.shape[-1]):
            sitk_img = sitk.GetImageFromArray(img_array[:,:,i])
            maskImage = sitk.GetImageFromArray((img_array[:,:,i]>0).astype(int))
            corrector = sitk.N4BiasFieldCorrectionImageFilter()
            sitk_img = sitk.Cast(sitk_img, sitk.sitkFloat32)
            maskImage = sitk.Cast(maskImage, sitk.sitkUInt8)
            output = corrector.Execute( sitk_img, maskImage )
            result[:,:,i] = sitk.GetArrayFromImage(output)
        return result

#create masks in raw folder, then move mask.nii into /masks and moves masked region into /brain
def masks():
    os.system('singularity run --no-home -B ./'+input_fol+'/raw:/data /neuro/labs/grantlab/research/MRI_processing/sofia.urosa/mask_project/singularity/brain_masking.sif /data;')
    os.system('mv '+input_fol+'/raw/*mask.nii '+input_fol+'/masks;')
        
    img_list= np.asarray(sorted(glob.glob(input_fol+'/masks/*mask.nii')))
    for i in range(len(img_list)):
        vol = nib.load(img_list[i])
        vol_data = vol.get_data()
        if np.max(vol_data)>0.01:
            os.system('mri_mask '+img_list[i].replace('masks/','raw/').replace('_mask.nii','.nii')+' '+img_list[i]+'  '+img_list[i].replace('masks/','brain/').replace('_mask.nii','_brain.nii'))
   ##verify images
    os.system('cp -r ./'+input_fol+'/brain/ ./'+input_fol+'/verify/;')
    verify()

### REMASKING option// creates new brain folder after manual mask correction 
def remask():
    if os.path.exists(input_fol+'/masks') and os.path.exists(input_fol+'/raw'):
        img_list= np.asarray(sorted(glob.glob(input_fol+'/masks/*mask.nii')))
        for i in range(len(img_list)):
            vol = nib.load(img_list[i])
            vol_data = vol.get_data()
            if np.max(vol_data)>0.01:
                os.system('mri_mask '+img_list[i].replace('masks/','raw/').replace('_mask.nii','.nii')+' '+img_list[i]+'  '+img_list[i].replace('masks/','brain/').replace('_mask.nii','_brain.nii'))
    ##verify images
    os.system('cp -r ./'+input_fol+'/brain/ ./'+input_fol+'/verify/;')
    verify()

#NUC
def nuc():
    img_list = np.asarray(sorted(glob.glob(input_fol+'/brain/*.nii')))
    for i in range(len(img_list)):
        os.system('~/arch/Linux64/packages/ANTs/current/bin/N4BiasFieldCorrection -d 3 -o '+img_list[i].replace('/brain/','/nuc/')+' -i '+img_list[i])

#Quality assessment 
def qa():
    os.system('singularity exec docker://fnndsc/pl-fetal-brain-assessment:1.3.0 fetal_brain_assessment ./'+input_fol+'/nuc/ ./'+input_fol+'/;')

#reconstruction
def recon():
    threshold = 0.4
    Best_list = []
    with open(input_fol+'/quality_assessment.csv', 'r') as file:
        reader = csv.DictReader(file)
        for row in reader:
            score = float(row["quality"]) #Reading QC score and cast it to float
            if score >= threshold: 
                Best_list.append(row)

    Best_list_sorted = sorted(Best_list, key=lambda row: row['quality'], reverse=True)
    
    Second_list_sorted = Best_list_sorted[0:3]
    Second_list_sorted = [Second_list_sorted[i] for i in [1, 0, 2]]
    Second_list_sorted.extend(Best_list_sorted[3:])

    Third_list_sorted = Best_list_sorted[0:3]
    Third_list_sorted = [Third_list_sorted[i] for i in [2, 0, 1]]
    Third_list_sorted.extend(Best_list_sorted[3:])
  
    cmd = [
    'reconstruction', 
    input_fol+'/temp_recon_1/recon.nii',
    str(len(Best_list_sorted)),
    [row['filename'] for row in Best_list_sorted],
    [(['id'] * len(Best_list_sorted))],
    '-thickness',
    [row['slice_thickness'] for row in Best_list_sorted],
    '-packages',
    [(['1'] * len(Best_list_sorted))]
    ]
    os.system('/neuro/arch/Linux64/packages/irtk/build-tbb/bin/'+repr(cmd).replace(",","").replace("'","").replace("[","").replace("]",""))

    cmd = [
    'reconstruction', 
    input_fol+'/temp_recon_2/recon.nii',
    str(len(Second_list_sorted)),
    [row['filename'] for row in Second_list_sorted],
    [(['id'] * len(Second_list_sorted))],
    '-thickness',
    [row['slice_thickness'] for row in Second_list_sorted],
    '-packages',
    [(['1'] * len(Second_list_sorted))]
    ]
    os.system('/neuro/arch/Linux64/packages/irtk/build-tbb/bin/'+repr(cmd).replace(",","").replace("'","").replace("[","").replace("]",""))

    cmd = [
    'reconstruction', 
    input_fol+'/temp_recon_3/recon.nii',
    str(len(Third_list_sorted)),
    [row['filename'] for row in Third_list_sorted],
    [(['id'] * len(Third_list_sorted))],
    '-thickness',
    [row['slice_thickness'] for row in Third_list_sorted],
    '-packages',
    [(['1'] * len(Third_list_sorted))]
    ]
    os.system('/neuro/arch/Linux64/packages/irtk/build-tbb/bin/'+repr(cmd).replace(",","").replace("'","").replace("[","").replace("]",""))  

#alignment
def align():
    rnum=("1","2","3")
    for r in rnum:
        os.makedirs(input_fol+'/temp_recon_'+r+'/alignment_temp/', exist_ok = True)

        if not os.path.exists(input_fol+'/temp_recon_'+r+'/alignment_temp/recon.nii'):
            os.system('cp '+input_fol+'/temp_recon_'+r+'/recon.nii '+input_fol+'/temp_recon_'+r+'/alignment_temp/recon.nii')
    
        os.system('cp -R /neuro/labs/grantlab/research/HyukJin_MRI/templates_for_alginment/ ~')

    def alignment_func(temp_recon_a):
        a_temp=("23", "24", "25", "26", "27", "28", "29", "30", "31", "32")
        for t in a_temp:
            os.system(\
                'flirt -in /neuro/labs/grantlab/research/HyukJin_MRI/templates_for_alginment/template-'+t+'/template-'+t+'.nii -ref '+temp_recon_a+'/recon.nii \
                -out '+temp_recon_a+'/Temp-Recon-7dof-'+t+'.nii -omat '+temp_recon_a+'/Temp-Recon-7dof-'+t+'.xfm -dof 7  -searchrx -180 180 -searchry -180 180 -searchrz -180 180;')
            os.system('flirt -in /neuro/labs/grantlab/research/HyukJin_MRI/templates_for_alginment/template-'+t+'/csf-'+t+'.nii -ref '+temp_recon_a+'/recon.nii \
                -out '+temp_recon_a+'/csf-aligned'+t+'.nii -init '+temp_recon_a+'/Temp-Recon-7dof-'+t+'.xfm -applyxfm;')

        recon = nib.load(temp_recon_a+'/recon.nii') # Load reconstruction image 
        size=recon.get_fdata().shape # Get the dimensions of the volume 
        meas = [0,0,0]
        beginning=[0,0,0]
    
        for i in range (0, 3):
            beginning[i]=int(round(size[i]/10.0))
            temp=size[i]-beginning[i]
            meas[i]=temp-beginning[i]
            
            if (i==0):
                dimensions=zeros([size[1],size[2]])
                coronal= dict.fromkeys(range(0, meas[i]),dimensions)
            if (i==1):
                dimensions=zeros([size[0],size[2]])
                sagital= dict.fromkeys(range(0, meas[i]),dimensions)
            if (i==2):
                dimensions=zeros([size[0],size[1]])
                axial= dict.fromkeys(range(0, meas[i]),dimensions)
        
        ccl= meas[0]+meas[1]+meas[2]

        im = {'corrcoef':zeros([1,ccl*10]), 'greatest':zeros([1,ccl]), 'template':zeros([1,ccl])}

        for i in range (0, meas[0]): #Gets the number of slides choosen
            a= beginning[0]+i #Starts in the slide selected as beginning and ends passing the one selected as end
            coronal[i]=np.uint8(np.squeeze(recon.get_fdata()[a,:,:])/4)  #Get the slide of the reconstruction image.

        for i in range (0, meas[1]):    
            a= beginning[1]+i #Starts in the slide selected as beginning and ends passing the one selected as end
            sagital[i]=np.uint8(np.squeeze(recon.get_fdata()[:,a,:])/4)  #Get the slide of the reconstruction image.

        for i in range (0, meas[2]):    
            a= beginning[2]+i #Starts in the slide selected as beginning and ends passing the one selected as end
            axial[i]=np.uint8(np.squeeze(recon.get_fdata()[:,:,a])/4)  #Get the slide of the reconstruction image.
                    
        number=22 #Stablishes the base to number the templates

        def mean2(value):
            mean2value=np.sum(value)/np.size(value)
            return mean2value

        def corr2(R, T):
            R=R-mean2(R)
            T=T-mean2(T)
            corr=((R*T).sum())/(math.sqrt((R*R).sum()*(T*T).sum()))
            return corr
        
        for i in range (0,10): #Says that the process will repeat for the 10 templates
            number=number+1 #The first template will be 23
            volume=''+temp_recon_a+'/csf-aligned%d.nii.gz' %number #Construct the name of the template volume that will be loaded
            volume=nib.load(volume) #Load the template volume
            meascor=(meas[0]*i)
            meassag=(meas[1]*i)+(meas[0]*9)
            measax=(meas[2]*i)+(meas[0]*9)+(meas[1]*9)
            
            for j in range (0,ccl):
                if (j<meas[0]):
                    t=meascor+j #Define the position in which the results will be stored\
                    a=beginning[0]+j; #Select the slide that will be taken
                    slide=volume.get_fdata()[a,:,:] #Loads the slide
                    
                if (j>=meas[0] and j<(meas[0]+meas[1])):
                    t=meassag+j #Define the position in which the results will be stored\
                    a=beginning[0]+j-meas[0]; #Select the slide that will be taken
                    slide=volume.get_fdata()[:,a,:] #Loads the slide
                
                if (j>=(meas[0]+meas[1])):
                    t=measax+j #Define the position in which the results will be stored\
                    a=beginning[0]+j-meas[0]-meas[1]; #Select the slide that will be taken
                    slide=volume.get_fdata()[:,:,a] #Loads the slide    
                    
            #normalize the slide
                slide=np.uint8(256*(slide-slide.min())/(slide.max()-slide.min()))
            
                index=np.nonzero(slide) #get the positions of the non zero values 
                csfi=slide[np.nonzero(slide)] #get the values of the no zero values
            
                temp=csfi.shape #gets the number of non zero values 
                temp=temp[0]
                reconi=[]
                
                if (j<meas[0]):
                    for n in range (0,temp): #get the same indexes of the recon image
                        slide=coronal[j].item(index[0][n],index[1][n])
                        reconi.append(slide)
                                
                if (j>=meas[0] and j<(meas[0]+meas[1])):
                    for n in range (0,temp): #get the same indexes of the recon image
                        jj=j-meas[0]
                        slide=sagital[jj].item(index[0][n],index[1][n])
                        reconi.append(slide)
                        
                if (j>=(meas[0]+meas[1])):
                    for n in range (0,temp): #get the same indexes of the recon image
                        jj=j-meas[0]-meas[1]
                        slide=axial[jj].item(index[0][n],index[1][n])
                        reconi.append(slide)
                        
                im['corrcoef'][0,t]=corr2(reconi, csfi)
            
            a=0
                
        for a in range (0,ccl):
            im['greatest'][0,a]=-5
            
            for i in range (0,10):
                if (a<meas[0]):
                        j=(meas[0]*i)+a
                                        
                if (a>=meas[0] and a<(meas[1]+meas[0]) ):
                    j=(meas[1]*i)+(meas[0]*9)+a
                    
                if (a>=(meas[1]+meas[0])):    
                    j=(meas[2]*i)+(meas[0]*9)+(meas[1]*9)+a            
                    
                if im['corrcoef'][0,j]>im['greatest'][0,a]:
                    im['greatest'][0,a]=im['corrcoef'][0,j]
                    im['template'][0,a]=i+23

        t=im['template']
        t=t[0]
        a=t.tolist()
        j=im['greatest']
        jj=max(j[0])-((max(j[0]))/2.5)

        for i in range (0,ccl):
            if (j[0][i]>jj):
                a.append(im['template'][0,i])

        n=np.histogram(t,bins=[23,24,25,26,27,28,29,30,31,32,33])
        temp=np.argsort(n[0])[::-1]    
        tempi=temp+23
        temp=str(tempi[0])

        os.system('convert_xfm -omat '+temp_recon_a+'/InvAligned-'+temp+'.xfm -inverse '+temp_recon_a+'/Temp-Recon-7dof-'+temp+'.xfm')
        os.system('convert_xfm -omat '+temp_recon_a+'/recon_to31.xfm -concat /neuro/labs/grantlab/research/HyukJin_MRI/templates_for_alginment/template-'+temp+'/template-'+temp+'to31.xfm '+temp_recon_a+'/InvAligned-'+temp+'.xfm')
        
        os.system('flirt -in '+temp_recon_a+'/recon.nii -ref /neuro/labs/grantlab/research/HyukJin_MRI/templates_for_alginment/template-31/template-31.nii \
                    -out '+temp_recon_a+'/recon_to31.nii.gz -init '+temp_recon_a+'/recon_to31.xfm -applyxfm')
    
        os.system('gunzip '+temp_recon_a+'/recon_to31.nii.gz')

    alignment_func(input_fol+'/temp_recon_1/alignment_temp/')
    alignment_func(input_fol+'/temp_recon_2/alignment_temp/')
    alignment_func(input_fol+'/temp_recon_3/alignment_temp/')

   
    for r in rnum:
        os.system('cp '+input_fol+'/temp_recon_'+r+'/alignment_temp/recon_to31.* '+input_fol+'/temp_recon_'+r)

        os.system('~/arch/Linux64/packages/ANTs/current/bin/N4BiasFieldCorrection -d 3 -o '+input_fol+'/temp_recon_'+r+'/recon_to31_nuc.nii -i '+input_fol+'/temp_recon_'+r+'/recon_to31.nii;')

        os.system('convert_xfm -omat '+input_fol+'/temp_recon_'+r+'/recon_to31_inv.xfm -inverse '+input_fol+'/temp_recon_'+r+'/recon_to31.xfm;')
        os.system('echo `avscale '+input_fol+'/temp_recon_'+r+'/recon_to31_inv.xfm | grep Scales` > '+input_fol+'/temp/temp.txt;')
        scales=open(input_fol+'/temp/temp.txt', encoding='utf-8')
        scales=scales.read()
        os.system('param2xfm -clobber -scales '+scales[16:-1]+' '+input_fol+'/temp_recon_'+r+'/recon_native.xfm;')

#auto segmentation
def auto_seg():
    os.system('singularity run --no-home -B '+input_fol+'/temp_recon_1/:/data --nv /neuro/labs/grantlab/research/MRI_processing/sungmin.you/MRI_SIF/fetal_cp_seg_att.sif recon_to31_nuc.nii . 1;')
    os.system('singularity run --no-home -B '+input_fol+'/temp_recon_2/:/data --nv /neuro/labs/grantlab/research/MRI_processing/sungmin.you/MRI_SIF/fetal_cp_seg_att.sif recon_to31_nuc.nii . 1;')
    os.system('singularity run --no-home -B '+input_fol+'/temp_recon_3/:/data --nv /neuro/labs/grantlab/research/MRI_processing/sungmin.you/MRI_SIF/fetal_cp_seg_att.sif recon_to31_nuc.nii . 1;')
                    # add ((optional)) seg_.sif version for attention model


def main():
    args = parser.parse_args()
    print(args)

    global input_fol
    input_fol = args.input_fol		## input_fol should be: ./$file

    #create subfolders inside main data directory 
    foldr=['raw', 'masks', 'brain', 'nuc', 'temp', 'temp_recon_1', 'temp_recon_2', 'temp_recon_3']
    for items in foldr:
	    os.makedirs(input_fol+'/'+items, exist_ok=True)
    
    #move raw scans into 'raw' folder from data dir
    dest=(input_fol+'/raw')		#dest folder for raw images
    if os.path.isdir(input_fol) and os.path.isdir(dest):
        for nii_file in glob.glob(input_fol+'/*nii*'):
            shutil.move(nii_file, dest) 

    masking = args.masks
    remasking = args.remask
    NUC = args.NUC
    QA = args.QA
    reconstruction = args.recon
    alignment = args.align
    seg = args.auto_seg
    from_remask = args.remask__
    fromNUC = args.nuc__
    fromQA = args.qa__
    from_recon = args.recon__
    from_align = args.align__
    allSteps = args.all
    rm = args.rm


    if masking == True:
        masks()
    if remasking == True: 
        remask()
    if NUC == True:
        nuc()
    if QA == True:
        qa()
    if reconstruction == True:
        recon()
    if alignment == True:
        align()
    if seg == True:
        auto_seg()
	
    if from_remask == True:
        remask()
        nuc()
        qa()
        recon()
        align()
        auto_seg()

    if fromNUC == True:
        nuc()
        qa()
        recon()
        align()
        auto_seg()

    if fromQA == True:
        qa()
        recon()
        align()
        auto_seg()

    if from_recon == True:
        recon()
        align()
        auto_seg()

    if from_align == True:
        align()
        auto_seg()

    if allSteps == True:
            masks()
            nuc()
            qa()
            recon()
            align()
            auto_seg()

    if rm == True:
        if os.path.exists(input_fol+'/temp_recon_1'):
            if input("Delete temp_recon files? (y/n)") == "y":
                os.system('rm -r '+input_fol+'/temp_recon_?/')
            exit()
        if input("Delete all other files & folders except for masks/ and raw/? (y/n)") == "y":
            os.system('rm -r '+input_fol+'/Best_Images_crop/ '+input_fol+'/brain/ '+input_fol+'/nuc/ '+input_fol+'/quality_assessment.csv '+input_fol+'/recon/ '+input_fol+'/temp_recon_?/')
            exit()


if __name__ == '__main__':
    main()

