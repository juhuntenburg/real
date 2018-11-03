%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% QualityGuidedUnwrap2D implements 2D quality guided path following phase
% unwrapping algorithm.
%
% Inputs: 1. Complex image in .mat double format
%         2. Binary mask (optional)
% Outputs: 1. Unwrapped phase image
%          2. Phase quality map
%
% This code can easily be extended for 3D phase unwrapping.
% Technique adapted from:
% D. C. Ghiglia and M. D. Pritt, Two-Dimensional Phase Unwrapping:
% Theory, Algorithms and Software. New York: Wiley-Interscience, 1998.
%
% Posted by Bruce Spottiswoode on 22 December 2008
% 2010/07/23  Modified by Carey Smith
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%% Consider adding thermal noise removal %%%%%%%%%%%

clear all;
addpath(genpath('/home/julia/workspace/real/'));

dataset='CL181030fmrssouris3';
scan='26';
recon='3';
data_dir=strcat('/home/julia/projects/real_data/mouse_visual/', dataset, filesep, 'raw/', scan, filesep);
mask_mat=strcat('/home/julia/projects/real_data/mouse_visual/', dataset, filesep, 'processed/func/', scan, filesep, 'func_mask.mat');

acq_params=readBrukerParamFile(strcat(data_dir,'acqp'));
method_params=readBrukerParamFile(strcat(data_dir,'method'));
visu_params=readBrukerParamFile(strcat(data_dir, filesep, 'pdata', filesep, recon, filesep, 'visu_pars'));
[data, visu] = readBruker2dseq(strcat(data_dir, filesep, 'pdata', filesep, recon, filesep,'2dseq'), visu_params);
data=squeeze(data);


%%%%% load mask mask
mask_struct=load(mask_mat);
mask=mask_struct.data;

[nx,ny,nz,nt,ncoil]=size(data);

corrected_real_data = zeros([nx,ny,nz,nt]);
% add_mag_data = zeros([nx,ny,nz,nt]);
sos_mag_data = zeros([nx,ny,nz,nt]);

rowref = round(nx./2);
colref = round(ny./2);

parfor z=1:nz %%%% loop over slices

    disp(z)
    im_mask = mask(:,:,z);

    for t=1:nt %%%% loop over timepoints

        UPIm = zeros([nx,ny,ncoil]);
        mag_coils = zeros([nx,ny,ncoil]);
        orig_data = squeeze(data(:,:,z,t,:));

        %%% check if FFT needs to be applied (if the middle point has a
        %%% lower value than the edges, present in all coils so just
        %%% checking for 1)
        if abs(orig_data(rowref, colref, 1)) < 2.*(abs(orig_data(rowref,  ny-2)));
            for c=1:ncoil;
                orig_data(:,:,c) = fftshift(orig_data(:,:,c),2);
            end;
        end;

        for c=1:ncoil;   %%%% loop over coils

            IM=orig_data(:,:,c);
            im_mag   = abs(IM);                  %Magnitude image
            im_phase = angle(IM);                %Phase image

            tmp = im_mag.*double(im_mask);
            mag_max = max(tmp(:));

            im_unwrapped = nan(size(IM));        %Initialze the output unwrapped version of the phase
            adjoin = zeros(size(IM));            %Zero starting matrix for adjoin matrix
            unwrapped_binary = zeros(size(IM));  %Binary image to mark unwrapped pixels

            %% Calculate phase quality map
            im_phase_quality = PhaseDerivativeVariance_r1(im_phase);

            %% Automatically (default) or manually identify starting seed point on a phase quality map
%            minp = im_phase_quality(2:end-1, 2:end-1); minp = min(minp(:));
%            maxp = im_phase_quality(2:end-1, 2:end-1); maxp = max(maxp(:));
% %             %if(0)    % Chose starting point interactively
%            figure, imagesc(im_phase_quality,[minp maxp]), colormap(gray), colorbar, axis square, axis off; title('Phase quality map')
% %             uiwait(msgbox('Select known true phase reference phase point. Black = high quality phase; white = low quality phase.','Phase reference point','modal'));
% %             [xpoint,ypoint] = ginput(1);                %Select starting point for the guided floodfill algorithm
% %             colref = round(xpoint);
% %             rowref = round(ypoint);
% %             close;  % close the figure;
% %             else   % Chose starting point = max. intensity, but avoid an edge pixel
%              [rowrefn,colrefn] = find(im_mag(2:end-1, 2:end-1) >= 0.99*mag_max);
%              rowref = rowrefn(1)+1; % choose the 1st point for a reference (known good value)
%              colref = colrefn(1)+1; % choose the 1st point for a reference (known good value)
% %            end

            %% Unwrap
            im_unwrapped(rowref,colref) = im_phase(rowref,colref);                          %Save the unwrapped values
            unwrapped_binary(rowref,colref,1) = 1;
            if im_mask(rowref-1, colref, 1)==1;  adjoin(rowref-1, colref, 1) = 1; end       %Mark the pixels adjoining the selected point
            if im_mask(rowref+1, colref, 1)==1;  adjoin(rowref+1, colref, 1) = 1; end
            if im_mask(rowref, colref-1, 1)==1;  adjoin(rowref, colref-1, 1) = 1; end
            if im_mask(rowref, colref+1, 1)==1;  adjoin(rowref, colref+1, 1) = 1; end
            im_unwrapped = GuidedFloodFill_r1(im_phase, im_mag, im_unwrapped, unwrapped_binary, im_phase_quality, adjoin, im_mask);    %Unwrap

            UPIm(:,:,c) = im_unwrapped;
            mag_coils(:,:,c) = im_mag;
        end

        corrected_real_data(:,:,z,t) = sqrt(mean(real(orig_data.*(exp(-1i.*UPIm))),3));
        % add_mag_data(:,:,z,t) = sum(abs(mag_coils), 3);
        sos_mag_data(:,:,z,t) = sqrt(sum(abs(mag_coils).^2, 3));
    end
end

data = corrected_real_data;
save(strcat(data_dir,'data_real.mat'), 'data');
data = sos_mag_data;
save(strcat(data_dir,'data_mag.mat'),'data');
