import nibabel as nb
import numpy as np
import scipy.io as sio

dataset = "CL181102fmrssouris7"
scan = "27"
data_dir = "/home/julia/projects/real_data/mouse_visual/%s/"%(dataset)

hdr_file = data_dir + "raw/%s/converted/%s/%s_%s/%s_%s.nii.gz"%(scan, dataset, dataset, scan, dataset, scan)
mask_file = data_dir + "processed/func/%s/func_mask.nii.gz"%scan

mask = np.asarray(nb.load(mask_file).get_data(), dtype='float64')
aff = nb.load(hdr_file).affine
hdr = nb.load(hdr_file).header

# calculate tsnr
for i in ['real', 'mag']:
    data = np.nan_to_num(sio.loadmat(data_dir+"raw/%s/converted/data_%s.mat"%(scan,i))['data'])
    tsnr = np.mean(data,axis=3)/np.std(data)
    tsnr *= mask
    nb.Nifti1Image(data, aff, hdr).to_filename(data_dir+"processed/func/%s/data_%s.nii.gz"%(scan,i))
    nb.Nifti1Image(tsnr, aff, hdr).to_filename(data_dir+"processed/func/%s/data_%s_tsnr.nii.gz"%(scan,i))
