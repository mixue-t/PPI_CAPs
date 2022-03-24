function [inGM, nonGM] = findGMvoxels(scriptsDir)

mask      = load([scriptsDir '/dependencies/MNIBrainMask.mat']); % Load MNIBrainMask
nonGM     = find(mask.MNIBrainMask(:) < 1); 
inGM  = find(mask.MNIBrainMask(:) == 1);

end