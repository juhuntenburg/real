from glob import glob
from nipype.pipeline.engine import Node, Workflow, MapNode
import nipype.interfaces.fsl as fsl
import nipype.interfaces.utility as util
import nipype.interfaces.io as nio
import nipype.interfaces.nipy as nipy
import nipype.interfaces.ants as ants
import nipype.interfaces.afni as afni
import nipype.interfaces.freesurfer as fs
import nipype.algorithms.rapidart as ra
from functions import strip_rois_func, motion_regressors, median, selectindex, nilearn_denoise, weighted_avg, pca_denoising

dataset = 'CL181030fmrssouris3'
struct = ''
recons = ['real', 'mag']

vol_to_remove = 10
motion_norm = 0.04
z_thr = 2
TR = 1


# directories
working_dir = '/home/julia/projects/real_data/working_dir/'
data_dir= '/home/julia/projects/real_data/mouse_visual/'
out_dir = '/home/julia/projects/real_data/mouse_visual/''

# main workflow
preproc = Workflow(name='preproc')
preproc.base_dir = working_dir
preproc.config['execution']['crashdump_dir'] = preproc.base_dir + "/crash_files"

# iterate over sessions
recon_infosource = Node(util.IdentityInterface(fields=['recon']),
                  name='session_infosource')
recon_infosource.iterables=[('recon', recons)]

# select files
templates = {'func' : '{dataset}/recon/data_{recon}.nii.gz',
             'struct' : ''
           }
selectfiles = Node(nio.SelectFiles(templates, base_directory=data_dir),
                   name="selectfiles")
selectfiles.inputs.dataset = dataset
selectfiles.inputs.struct = struct
preproc.connect([(recon_infosource, selectfiles, [('recon', 'recon')])])


############################
# Functional preprocessing #
############################

# Remove first volumes
remove_vol = Node(util.Function(input_names=['in_file','t_min'],
                                output_names=["out_file"],
                                function=strip_rois_func),
                  name='remove_vol')
remove_vol.inputs.t_min = vol_to_remove
preproc.connect([(selectfiles, remove_vol, [('func', 'in_file')])])

# motion correction
moco = Node(nipy.SpaceTimeRealigner(slice_times='asc_alt_2', tr=TR, slice_info=[2,1]),name="moco")
preproc.connect([(remove_vol, moco, [('out_file', 'in_file')])])

# compute median
median = Node(util.Function(input_names=['in_files'],
                       output_names=['median_file'],
                       function=median),
              name='median')

preproc.connect([(moco, median, [('out_file', 'in_files')])])

# bias field correction
biasfield = Node(ants.segmentation.N4BiasFieldCorrection(dimension=3,
                    n_iterations=[150,100,50,30], convergence_threshold=1e-11,
                    bspline_fitting_distance = 10, bspline_order = 4,
                    shrink_factor = 2, output_image='func_median.nii.gz'),name='biasfield')
preproc.connect([(median, biasfield, [('median_file', 'input_image')])])

# compute functional mask
func_mask = Node(afni.Automask(dilate=1, args='-peels 3',
                               outputtype='NIFTI_GZ', out_file='func_mask.nii.gz'),
                 name='func_mask')

preproc.connect([(biasfield, func_mask, [('output_image', 'in_file')])])


# artefact detection
artefact = Node(ra.ArtifactDetect(save_plot=True,
                                  use_norm=True,
                                  parameter_source='NiPy',
                                  mask_type='file',
                                  norm_threshold=motion_norm,
                                  zintensity_threshold=z_thr,
                                  use_differences=[True,False]),
                 name='artefact')

preproc.connect([(moco, artefact, [('out_file', 'realigned_files'),
                                   ('par_file', 'realignment_parameters')]),
                 (func_mask, artefact, [('out_file', 'mask_file')]),
                 ])

# calculate motion regressors
motreg = Node(util.Function(input_names=['motion_params', 'order','derivatives'],
                            output_names=['out_files'],
                            function=motion_regressors),
                 name='motion_regressors')
motreg.inputs.order=2
motreg.inputs.derivatives=2
preproc.connect([(moco, motreg, [('par_file','motion_params')])])

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


############################
# Structural preprocessing #
############################

# Bias field correction
struct_bias = Node(ants.N4BiasFieldCorrection(dimension=3,
                                              n_iterations=[100,100,100,100],
                                              convergence_threshold=0.0,),
                                              name='struct_bias')

preproc.connect([(selectfiles, struct_bias, [('struct', 'input_image')])])


# Skull stripping
skullstrip = Node(afni.SkullStrip(outputtype='NIFTI_GZ',
                                  args='-rat -push_to_edge'),
                  name='skullstrip')
preproc.connect([(struct_bias, skullstrip, [('output_image','in_file')])])

# Binarize mask
struct_mask = Node(fs.Binarize(out_type = 'nii.gz', min=0.1,
                   binary_file='struct_mask.nii.gz'), name='struct_mask')
preproc.connect([(skullstrip, struct_mask, [('out_file','in_file')])])

# Apply mask
apply_struct = Node(fsl.ApplyMask(out_file='struct_masked.nii.gz'), name='apply_struct')

preproc.connect([(struct_mask, apply_struct, [('binary_file','mask_file')]),
                 (struct_bias, apply_struct, [('output_image','in_file')])
                ])


################
# Registration #
################

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
                         num_threads = 3,
                         initial_moving_transform_com = True,
                         )

preproc.connect([(apply_struct, coreg, [('out_file','fixed_image')]),
                 (apply_func, coreg, [('out_file','moving_image')])
                ])


#############
# Save data #
#############

# Sink relevant files
func_sink = Node(nio.DataSink(parameterization=False),name='func_sink')
func_sink.inputs.base_directory = out_dir + dataset
func_sink.inputs.regexp_substitutions = [('corr_.*_roi_denoised', 'func_final'),
                                         ('corr_.*_roi','func_moco')]
preproc.connect([(recon_infosource, func_sink, [('recon', 'container')]),
                 (moco, func_sink, [('out_file', '@realigned_file'),
                                    ('par_file', 'confounds.@orig_motion')]),
                 (func_mask, func_sink, [('out_file', '@mask')]),
                 (biasfield, func_sink, [('output_image', '@median')]),
                 (artefact, func_sink, [('norm_files', 'confounds.@norm_motion'),
                                   ('outlier_files', 'confounds.@outlier_files'),
                                   ('intensity_files', 'confounds.@intensity_files'),
                                   ('statistic_files', 'confounds.@outlier_stats'),
                                   ('plot_files', 'confounds.@outlier_plots')]),
                 (apply_func, func_sink, [('out_file', '@masked')]),
                 # (motreg, func_sink, [('out_files', 'confounds.@motreg')]),
                 # (regress, func_sink, [('denoised_img', '@denoised_img'),
                 #                  ('denoised_data', '@denoised_data'),
                 #                  ('confounds', 'confounds.@confounds')])
                 (coreg, func_sink, [('warped_image', 'coreg.@masked'),
                                     ('forward_transforms', 'coreg.@fwd'),
                                     ('reverse_transforms', 'coreg.@rvs')])
                ])

struct_sink = Node(nio.DataSink(parameterization=False),name='struct_sink')
struct_sink.inputs.base_directory = out_dir + dataset + '/struct'
preproc.connect([(struct_mask, struct_sink, [('binary_file', '@mask')]),
                 (struct_bias, struct_sink, [('output_image', '@corrected')]),
                 (apply_struct, struct_sink, [('out_file', '@masked')])
                 ])


preproc.run(plugin='MultiProc', plugin_args={'n_procs' : 2})
