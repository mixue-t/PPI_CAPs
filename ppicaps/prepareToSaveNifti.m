function [voxel_size, voxel_shift] = prepareToSaveNifti(analyze_dir, VOI_name)

VOI = spm_select('FPListRec',analyze_dir,[VOI_name '_1.mat$']); 
VOI = load(VOI);

voxel_size = diag(VOI.xY.spec.mat);
voxel_size = voxel_size(1:end-1)';

voxel_shift = VOI.xY.spec.mat(:,4);
voxel_shift = voxel_shift(1:end-1)';

end
