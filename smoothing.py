from nipype.interfaces.spm import Smooth
from glob import glob


in_files=glob('/home/julia/projects/real_data/mouse_visual/*/processed/func/*/*/func_moco.nii')

smooth = Smooth(in_files=in_files,
                fwhm=[0.145, 0.145, 0.145],
                paths=['/home/julia/software/spm/spm12'])
res = smooth.run()
