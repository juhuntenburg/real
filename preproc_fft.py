from nipype.pipeline.engine import Node, Workflow
import nipype.interfaces.fsl as fsl
import nipype.interfaces.utility as util
import nipype.interfaces.io as nio
import nipype.algorithms.rapidart as ra
from functions.functions import motion_regressors, median, nilearn_denoise, selectindex, fft_analysis

dataset = 'CL181028fmrsmouse1'
scans = ['9','33']
recons = ['real', 'mag']

bandpass = [None, 0.01]
vol_to_remove = 10
motion_norm = 0.6
z_thr = 3
tr = 1
peaks = [6,7,8,13,14,15]


# directories
working_dir = '/home/julia/projects/real_data/working_dir/'
data_dir= '/home/julia/projects/real_data/mouse_visual/%s/'%dataset
out_dir = '/home/julia/projects/real_data/mouse_visual/%s/processed/func/'%dataset

# main workflow
preproc_func = Workflow(name='preproc_func')
preproc_func.base_dir = working_dir
preproc_func.config['execution']['crashdump_dir'] = preproc_func.base_dir + "/crash_files"

# iterate over scans
scan_infosource = Node(util.IdentityInterface(fields=['scan']),
                  name='scan_infosource')
scan_infosource.iterables=[('scan', scans)]

# iterate over recons
recon_infosource = Node(util.IdentityInterface(fields=['recon']),
                  name='recon_infosource')
recon_infosource.iterables=[('recon', recons)]

# select files
templates = {'func' : 'processed/func/{scan}/{recon}/sfunc_moco.nii',
             'par' : 'processed/func/{scan}/{recon}/confounds/data_{recon}_roi_masked.nii.gz.par',
             'mask' : 'processed/func/{scan}/func_mask.nii.gz',
            }
selectfiles = Node(nio.SelectFiles(templates, base_directory=data_dir),
                   name="selectfiles")
preproc_func.connect([(recon_infosource, selectfiles, [('recon', 'recon')]),
                      (scan_infosource, selectfiles, [('scan', 'scan')])])

# compute median
median = Node(util.Function(input_names=['in_files'],
                       output_names=['median_file'],
                       function=median),
              name='median')

preproc_func.connect([(selectfiles, median, [('func', 'in_files')])])

# artefact detection
artefact = Node(ra.ArtifactDetect(save_plot=True,
                                  use_norm=True,
                                  parameter_source='NiPy',
                                  mask_type='file',
                                  norm_threshold=motion_norm,
                                  zintensity_threshold=z_thr,
                                  use_differences=[True,False]),
                 name='artefact')

preproc_func.connect([(selectfiles, artefact, [('func', 'realigned_files'),
                                               ('par', 'realignment_parameters'),
                                               ('mask', 'mask_file')])
                     ])

# calculate motion regressors
motreg = Node(util.Function(input_names=['motion_params', 'order','derivatives'],
                            output_names=['out_files'],
                            function=motion_regressors),
                 name='motion_regressors')
motreg.inputs.order=1
motreg.inputs.derivatives=1
preproc_func.connect([(selectfiles, motreg, [('par','motion_params')])])

# noise regression
regress = Node(util.Function(input_names=['in_file', 'brain_mask',
                                          'motreg_file', 'outlier_file',
                                          'bandpass', 'tr'],
                             output_names=['denoised_img', 'denoised_data',
                                           'confounds'],
                             function=nilearn_denoise), name='regress')
regress.inputs.tr = tr
regress.inputs.bandpass = bandpass

preproc_func.connect([(selectfiles, regress, [('func', 'in_file'),
                                               ('mask', 'brain_mask')]),
                      (motreg, regress, [(('out_files',selectindex,[0]), 'motreg_file')]),
                      (artefact, regress, [('outlier_files', 'outlier_file')])
                      ])

fft = Node(util.Function(input_names=['in_file', 'peaks', 'tr'],
                         output_names=['out_file'],
                         function=fft_analysis), name='fft')
fft.inputs.peaks = peaks
fft.inputs.tr = tr
preproc_func.connect([(regress, fft, [('denoised_img','in_file')])])

def makebase(scan, out_dir):
    return out_dir + scan + '/'

# Sink relevant files
func_sink = Node(nio.DataSink(parameterization=False),name='func_sink')
preproc_func.connect([(recon_infosource, func_sink, [('recon', 'container')]),
                      (scan_infosource, func_sink, [(('scan', makebase, out_dir), 'base_directory')]),
                      (median, func_sink, [('median_file', '@median')]),
                      (artefact, func_sink, [('norm_files', 'fft.@norm_motion'),
                                             ('outlier_files', 'fft.@outlier_files'),
                                             ('plot_files', 'fft.@outlier_plots')]),
                      (regress, func_sink, [('confounds', 'fft.@confounds'),
                                            ('denoised_img', 'fft.@img')]),
                      (fft, func_sink, [('out_file', 'fft.@median')]),
                     ])


preproc_func.run(plugin='MultiProc', plugin_args={'n_procs' : 3})
