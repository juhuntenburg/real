from nipype.pipeline.engine import Node, Workflow
import nipype.interfaces.fsl as fsl
import nipype.interfaces.utility as util
import nipype.interfaces.io as nio
import nipype.interfaces.nipy as nipy
import nipype.interfaces.ants as ants
import nipype.interfaces.afni as afni
import nipype.interfaces.freesurfer as fs
import nipype.algorithms.rapidart as ra
from functions import strip_rois_func, motion_regressors, median

dataset = 'CL181030fmrssouris3'
scans = ['4', '26']
recons = ['real', 'mag']

vol_to_remove = 10
motion_norm = 0.04
z_thr = 2
TR = 1


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
templates = {'func' : 'processed/func/{scan}/data_{recon}.nii.gz',
             'mask' : 'processed/func/{scan}/func_mask.nii.gz',
             'struct' : 'processed/struct/struct_masked.nii.gz'
            }
selectfiles = Node(nio.SelectFiles(templates, base_directory=data_dir),
                   name="selectfiles")
preproc_func.connect([(recon_infosource, selectfiles, [('recon', 'recon')]),
                      (scan_infosource, selectfiles, [('scan', 'scan')])])


# Remove first volumes
remove_vol = Node(util.Function(input_names=['in_file','t_min'],
                                output_names=["out_file"],
                                function=strip_rois_func),
                  name='remove_vol')
remove_vol.inputs.t_min = vol_to_remove
preproc_func.connect([(selectfiles, remove_vol, [('func', 'in_file')])])

# motion correction
moco = Node(nipy.SpaceTimeRealigner(slice_times='asc_alt_2', tr=TR, slice_info=[2,1]),name="moco")
preproc_func.connect([(remove_vol, moco, [('out_file', 'in_file')])])

# compute median
median = Node(util.Function(input_names=['in_files'],
                       output_names=['median_file'],
                       function=median),
              name='median')

preproc_func.connect([(moco, median, [('out_file', 'in_files')])])

# artefact detection
artefact = Node(ra.ArtifactDetect(save_plot=True,
                                  use_norm=True,
                                  parameter_source='NiPy',
                                  mask_type='file',
                                  norm_threshold=motion_norm,
                                  zintensity_threshold=z_thr,
                                  use_differences=[True,False]),
                 name='artefact')

preproc_func.connect([(moco, artefact, [('out_file', 'realigned_files'),
                                        ('par_file', 'realignment_parameters')]),
                      (selectfiles, artefact, [('mask', 'mask_file')]),
                 ])

# calculate motion regressors
motreg = Node(util.Function(input_names=['motion_params', 'order','derivatives'],
                            output_names=['out_files'],
                            function=motion_regressors),
                 name='motion_regressors')
motreg.inputs.order=2
motreg.inputs.derivatives=2
preproc_func.connect([(moco, motreg, [('par_file','motion_params')])])

# nuissance regression
# regress = Node(util.Function(input_names=['in_file', 'brain_mask',
#                                           'motreg_file', 'outlier_file',
#                                           'bandpass', 'tr'],
#                              output_names=['denoised_img', 'denoised_data',
#                                            'confounds'],
#                              function=nilearn_denoise), name='regress')
# regress.inputs.tr = tr
# regress.inputs.bandpass = bandpass
#
# preproc.connect([(moco, regress, [('out_file', 'in_file')]),
#                  (func_mask, regress, [('out_file', 'brain_mask')]),
#                  (motreg, regress, [(('out_files',selectindex,[0]), 'motreg_file')]),
#                  (artefact, regress, [('outlier_files', 'outlier_file')])
#                  ])



coreg = Node(ants.Registration(output_warped_image = out_dir + 'func2struct_mean.nii.gz',
                         output_transform_prefix = out_dir + 'func2struct_',
                         dimension = 3,
                         transforms = ['Rigid', 'SyN'],
                         metric = ['MI', 'CC'],
                         transform_parameters = [(0.1,),(0.1, 3.0, 0)],
                         metric_weight = [1,1],
                         radius_or_number_of_bins = [32,4],
                         sampling_percentage = [0.33, None],
                         sampling_strategy = ['Regular', None],
                         convergence_threshold = [1.e-11, 1.e-6],
                         convergence_window_size = [10,10],
                         smoothing_sigmas = [[0],[0,0]],
                         sigma_units = ['vox', 'vox'],
                         shrink_factors = [[1],[2,1]],
                         use_estimate_learning_rate_once = [False, False],
                         use_histogram_matching = [False, True],
                         number_of_iterations = [[300], [50,10]],
                         collapse_output_transforms = True,
                         winsorize_lower_quantile = 0.05,
                         winsorize_upper_quantile = 0.95,
                         args = '--float',
                         num_threads = 1,
                         initial_moving_transform_com = True,
                         ), name='coreg')

preproc_func.connect([(selectfiles, coreg, [('struct','fixed_image')]),
                    (median, coreg, [('median_file','moving_image')])
                    ])

def makebase(scan, out_dir):
    return out_dir + scan + '/'

# Sink relevant files
func_sink = Node(nio.DataSink(parameterization=False),name='func_sink')
func_sink.inputs.base_directory = out_dir
func_sink.inputs.regexp_substitutions = [('corr_.*_roi','func_moco')]
preproc_func.connect([(recon_infosource, func_sink, [('recon', 'container')]),
                      (scan_infosource, func_sink, [(('scan', makebase, out_dir), 'base_directory')]),
                      (moco, func_sink, [('out_file', '@realigned_file'),
                                         ('par_file', 'confounds.@orig_motion')]),
                      (artefact, func_sink, [('norm_files', 'confounds.@norm_motion'),
                                       ('outlier_files', 'confounds.@outlier_files'),
                                       ('intensity_files', 'confounds.@intensity_files'),
                                       ('statistic_files', 'confounds.@outlier_stats'),
                                       ('plot_files', 'confounds.@outlier_plots')]),
                      # (motreg, func_sink, [('out_files', 'confounds.@motreg')]),
                      # (regress, func_sink, [('denoised_img', '@denoised_img'),
                      #                  ('denoised_data', '@denoised_data'),
                      #                  ('confounds', 'confounds.@confounds')])
                      (coreg, func_sink, [('warped_image', 'coreg.@masked'),
                                         ('forward_transforms', 'coreg.@fwd'),
                                         ('reverse_transforms', 'coreg.@rvs')])
                    ])


preproc_func.run(plugin='MultiProc', plugin_args={'n_procs' : 2})
