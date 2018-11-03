from nipype.pipeline.engine import Node, Workflow
import nipype.interfaces.fsl as fsl
import nipype.interfaces.io as nio
import nipype.interfaces.ants as ants
import nipype.interfaces.afni as afni
import nipype.interfaces.freesurfer as fs

dataset = 'CL181030fmrssouris3'
struct = '3'


# directories
working_dir = '/home/julia/projects/real_data/working_dir/'
data_dir= '/home/julia/projects/real_data/mouse_visual/%s/raw/'%dataset
out_dir = '/home/julia/projects/real_data/mouse_visual/%s/processed/struct/'%dataset

# main workflow
preproc_struct = Workflow(name='preproc_struct')
preproc_struct.base_dir = working_dir
preproc_struct.config['execution']['crashdump_dir'] = preproc_struct.base_dir + "/crash_files"


# select files
templates = {'struct' : '{struct}/converted/*{struct}.nii.gz',
           }
selectfiles = Node(nio.SelectFiles(templates, base_directory=data_dir),
                   name="selectfiles")
selectfiles.inputs.struct = struct


# Bias field correction
struct_bias = Node(ants.N4BiasFieldCorrection(output_image='struct_corr.nii.gz',
                                              dimension=3,
                                              n_iterations=[100,100,100,100],
                                              convergence_threshold=0.0,),
                                              name='struct_bias',
                                              )

preproc_struct.connect([(selectfiles, struct_bias, [('struct', 'input_image')])])


# Reslice
reslice = Node(afni.Resample(voxel_size= (0.1,0.1,0.1), resample_mode='Cu',
                             out_file='struct_resliced.nii.gz'),
               name='reslice')
preproc_struct.connect([(selectfiles, reslice, [('struct', 'in_file')])])

# Skull stripping
skullstrip = Node(afni.SkullStrip(outputtype='NIFTI_GZ',
                                  args='-rat -surface_coil'),
                  name='skullstrip')
preproc_struct.connect([(reslice, skullstrip, [('out_file','in_file')])])

# Binarize mask
struct_mask = Node(fs.Binarize(out_type = 'nii.gz', min=0.1, erode=1),
                   name='struct_mask')
preproc_struct.connect([(skullstrip, struct_mask, [('out_file','in_file')])])

#Reslice mask
reslice_mask = Node(afni.Resample(resample_mode='NN', out_file='struct_mask.nii.gz'),
               name='reslice_mask')
preproc_struct.connect([(struct_bias, reslice_mask, [('output_image', 'master')]),
                        (struct_mask, reslice_mask, [('binary_file', 'in_file')])
                        ])

# Apply mask
apply_struct = Node(fsl.ApplyMask(out_file='struct_masked.nii.gz'), name='apply_struct')

preproc_struct.connect([(reslice_mask, apply_struct, [('out_file','mask_file')]),
                        (struct_bias, apply_struct, [('output_image','in_file')])
                        ])

# Sink relevant files
struct_sink = Node(nio.DataSink(parameterization=False),name='struct_sink')
struct_sink.inputs.base_directory = out_dir
preproc_struct.connect([(reslice_mask, struct_sink, [('out_file', '@mask')]),
                        (struct_bias, struct_sink, [('output_image', '@corrected')]),
                        (apply_struct, struct_sink, [('out_file', '@masked')])
                        ])


preproc_struct.run(plugin='MultiProc', plugin_args={'n_procs' : 1})
