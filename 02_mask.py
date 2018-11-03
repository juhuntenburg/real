import os
import nipype.interfaces.afni as afni
import nipype.interfaces.ants as ants
import nipype.interfaces.fsl as fsl
import scipy.io as sio
import nibabel as nb

dataset = "CL181030fmrssouris3"
scan = "26"

data_dir = "/home/julia/projects/real_data/mouse_visual/%s/"%dataset
in_file = data_dir + "raw/%s/converted/%s_%s.nii.gz"%(scan, dataset, scan)
mean_file = data_dir + "raw/%s/converted/%s_%s_mean.nii.gz"%(scan, dataset, scan)
bias_file = data_dir + "raw/%s/converted/%s_%s_mean_corr.nii.gz"%(scan, dataset, scan)
mask_file = data_dir + "processed/func/%s/func_mask.nii.gz"%scan
mask_mat = data_dir + "processed/func/%s/func_mask.mat"%scan

if not os.path.isdir(data_dir + "processed/func/%s"%scan):
    os.mkdir(data_dir + "processed/func/%s"%scan)

print("Compute mean")
mean = fsl.MeanImage(dimension='T',
                     in_file=in_file,
                     out_file=mean_file)
mean.run()

print("Bias field correction")
biasfield = ants.N4BiasFieldCorrection(dimension=3,
                                       n_iterations=[150,100,50,30],
                                       convergence_threshold=1e-11,
                                       bspline_fitting_distance = 10,
                                       bspline_order = 4,
                                       shrink_factor = 2,
                                       input_image = mean_file,
                                       output_image= bias_file)
biasfield.run()

print("Compute mask")
mask = afni.Automask(outputtype='NIFTI_GZ', in_file=bias_file, out_file=mask_file)
mask.run()

os.remove("%s_%s_mean_corr_masked.nii.gz"%(dataset, scan))

print("Create .mat version")
data=nb.load(mask_file).get_data()
sio.savemat(mask_mat, {'data':data})
