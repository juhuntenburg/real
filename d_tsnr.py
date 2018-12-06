import nibabel as nb
import numpy as np
import scipy.io as sio

dataset = "CL181030fmrssouris2"
scan = "30"
data_dir = "/home/julia/projects/real_data/mouse_visual/%s/"%(dataset)

real_file  = data_dir+"raw/%s/converted/data_real.mat"
mag_file = data_dir + "raw/%s/converted/%s/%s_%s/%s_%s.nii.gz"%(scan, dataset, dataset, scan, dataset, scan)
mask_file = data_dir + "processed/func/%s/func_mask.nii.gz"%scan

mask = np.asarray(nb.load(mask_file).get_data(), dtype='float64')
aff = nb.load(mag_file).affine
hdr = nb.load(mag_file).header

# calculate tsnr real
data = np.nan_to_num(sio.loadmat(real_file%scan)['data'])
tsnr = np.mean(data,axis=3)/np.std(data)
tsnr *= mask
nb.Nifti1Image(data, aff, hdr).to_filename(data_dir+"processed/func/%s/data_real.nii.gz"%(scan))
nb.Nifti1Image(tsnr, aff, hdr).to_filename(data_dir+"processed/func/%s/data_real_tsnr.nii.gz"%(scan))

# calculate tsnr mag
data = np.nan_to_num(nb.load(mag_file).get_data())
tsnr = np.mean(data,axis=3)/np.std(data)
tsnr *= mask
nb.Nifti1Image(data, aff, hdr).to_filename(data_dir+"processed/func/%s/data_mag.nii.gz"%(scan))
nb.Nifti1Image(tsnr, aff, hdr).to_filename(data_dir+"processed/func/%s/data_mag_tsnr.nii.gz"%(scan))
