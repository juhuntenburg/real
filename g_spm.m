%%%%%%%%%%%%%%%%%%%%%% Input %%%%%%%%%%%%%%%%%%%%%

data_dir='/home/julia/projects/real_data/mouse_visual/CL181102fmrssouris7/processed/func/';
scan = {'3', '27'};
recon = {'mag', 'real'};



spm('defaults','fMRI')


%%% Specify model

job_spec.timing.units='scans';
job_spec.timing.RT=1;
job_spec.timing.fmri_t=100;  % Own HRF: Microtime resolution = round((time-bins of hrf)/(hrf duration in seconds) * (TR in seconds)/ 1 scan) = 2000/20 * 1/1 = 100 time-bins/scan
                        % Default SPM HRF: Microtime resolution = number of slices
job_spec.timing.fmri_t0=56; % Own HRF: Microtime onset = round((microtime resolution)/(number of slices)* ((reference slice position in sliceorder)- 1/2))
                      % Default SPM HRF: Microtime onset = reference slice position in sliceorder
job_spec.sess.cond.name='visual';
job_spec.sess.cond.onset=[49,121,193,264,337]-1; % On either scale, the (start of the) first scan will be timepoint zero (not 1!).
job_spec.sess.cond.duration=[24,24,24,24,24];
job_spec.sess.cond.tmod=0;
job_spec.sess.cond.pmod=struct([]);
job_spec.sess.cond.orth=1;
job_spec.sess.multi=cell(1);
job_spec.sess.regress=struct([]);
job_spec.sess.hpf=90;
job_spec.fact=struct([]);
job_spec.bases.hrf.derivs=[0,0]; %job_spec.bases.none=1;
job_spec.volt=1;
job_spec.global='None';
job_spec.mthresh=1.000000000000000e-04;
job_spec.cvi='AR(1)';

%%% Estimate model

job_est.method.Classical=1;
job_est.write_residuals=0;

%%% Contrasts
job_con.consess{1}.tcon.name='1';
job_con.consess{1}.tcon.weights=1;
job_con.consess{1}.tcon.sessrep='none';
job_con.delete=1;

%%% Save maps
job_maps.conspec.titlestr='';
job_maps.conspec.contrasts=Inf;
job_maps.conspec.threshdesc='none';
job_maps.conspec.thresh=0.05;
job_maps.conspec.extent=8;
job_maps.conspec.mask.none=1;
job_maps.units=1;
job_maps.print='ps';
job_maps.export{1}.tspm.basename='005_8';
% job_maps.write.tspm.basename='005_8'; %previous spm version


%%%%%%%%%%%%%%%%%% Computations %%%%%%%%%%%%%%%%%%

for s=1:2
    for r=1:2
        job_spec.dir={[data_dir scan{s} filesep recon{r} filesep 'spm']};
        h1 = spm_vol([data_dir scan{s} filesep recon{r} filesep 'func_moco.nii']);
        h2 = struct2cell(h1);
        job_spec.sess.scans = h2(1,1);
        job_spec.sess.multi_reg={[data_dir scan{s} filesep recon{r} filesep 'confounds' filesep 'all_confounds.txt']};
        job_spec.mask={[data_dir scan{s} filesep recon{r} filesep 'func_mask.nii']};
        spm_run_fmri_spec(job_spec);

        job_est.spmmat={[data_dir scan{s} filesep recon{r} filesep 'spm' filesep 'SPM.mat']};
        spm_run_fmri_est(job_est);

        job_con.spmmat={[data_dir scan{s} filesep recon{r} filesep 'spm' filesep 'SPM.mat']};
        spm_run_con(job_con);

        job_maps.spmmat={[data_dir scan{s} filesep recon{r} filesep 'spm' filesep 'SPM.mat']};
        spm_run_results(job_maps);

    end
end
