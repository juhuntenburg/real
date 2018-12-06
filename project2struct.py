from nipype.pipeline.engine import Node, Workflow
import nipype.interfaces.fsl as fsl
import nipype.interfaces.utility as util
import nipype.interfaces.io as nio
import nipype.interfaces.nipy as nipy
import nipype.interfaces.ants as ants
import nipype.interfaces.afni as afni



dataset = 'CL181030fmrssouris3'
scans = ['4','26']
recons = ['real', 'mag']

# directories
working_dir = '/home/julia/projects/real_data/working_dir/'
data_dir = '/home/julia/projects/real_data/mouse_visual/%s/processed/'%dataset
out_dir ='/home/julia/projects/real_data/mouse_visual/%s/processed/func/'%dataset


# main workflow
func2struct = Workflow(name='func2struct')
func2struct.base_dir = working_dir
func2struct.config['execution']['crashdump_dir'] = func2struct.base_dir + "/crash_files"

# iterate over scans
scan_infosource = Node(util.IdentityInterface(fields=['scan']),
                  name='scan_infosource')
scan_infosource.iterables=[('scan', scans)]

# iterate over recons
recon_infosource = Node(util.IdentityInterface(fields=['recon']),
                  name='recon_infosource')
recon_infosource.iterables=[('recon', recons)]

# select files
templates = {'spm' : 'func/{scan}/{recon}/fft/sfunc_moco_denoised_fft.nii',
             'struct' : 'struct/struct_masked.nii.gz',
             'trans' : 'func/{scan}/mag/func2struct_0GenericAffine.mat',
            }
selectfiles = Node(nio.SelectFiles(templates, base_directory=data_dir),
                   name="selectfiles")
func2struct.connect([(recon_infosource, selectfiles, [('recon', 'recon')]),
                      (scan_infosource, selectfiles, [('scan', 'scan')])])



trans = Node(ants.ApplyTransforms(dimension=3,
                                  invert_transform_flags=[False],
                                  output_image='fft2struct.nii',
                                  interpolation = 'BSpline'),name='trans')


func2struct.connect([(selectfiles, trans, [('spm', 'input_image'),
                                           ('struct', 'reference_image'),
                                           ('trans', 'transforms')])])

def makebase(scan, out_dir):
    return out_dir + scan + '/'

# Sink relevant files
sink = Node(nio.DataSink(parameterization=False),name='sink')
func2struct.connect([(recon_infosource, sink, [('recon', 'container')]),
                      (scan_infosource, sink, [(('scan', makebase, out_dir), 'base_directory')]),
                      (trans, sink, [('output_image', 'fft.@trans')])
                    ])



func2struct.run()
