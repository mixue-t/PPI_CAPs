spm_file = 'D:\data\211214_additionalPilots\fmri\data4analyses\sub-01\FFX\FFX_model_block_30s\SPM.mat';
matlabbatch{1}.spm.util.voi.spmmat = cellstr(fullfile(deblank(spm_file)));
matlabbatch{1}.spm.util.voi.adjust = 1;
matlabbatch{1}.spm.util.voi.session = 1;
matlabbatch{1}.spm.util.voi.name = 'Thalamus';
% Auditory3 is in T2 space but will automatically resliced to functional space by SPM
matlabbatch{1}.spm.util.voi.roi{1}.mask.image = {'D:\data\211214_additionalPilots\fmri\data4analyses\sub-01\FFX\FFX_model_block_30s\thalamus.nii,1'}; 
matlabbatch{1}.spm.util.voi.roi{1}.mask.threshold = 0.5;
matlabbatch{1}.spm.util.voi.expression = 'i1';
spm_jobman('run',matlabbatch);