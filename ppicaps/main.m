% -------------------------------------------------------
%
%    main - function to execute PPI-CAPs workflow
%
%    Created by:         Lorena Freitas 
%    Last checked:       28.09.2019
%
%   
% ------------------------------------------------------

function main

% ------------------------------------------------------
% Set up
% ------------------------------------------------------

% Distance used for kmeans 
dist        = 'cosine';

% Subjects in analysis
sujets      = {'Sujet01', 'Sujet10', 'Sujet16', 'Sujet17', 'Sujet19' ...
               'Sujet21', 'Sujet26', 'Sujet27', 'Sujet31', 'Sujet32', ...
               'Sujet09', 'Sujet14', 'Sujet20', 'Sujet23', 'Sujet25' };

ppicapsDir = 'D:\data\211214_additionalPilots\fmri\ppi_caps';  
pilot_name = 'AP';
mod_name = 'FFX_model_block_30s';
mod_path = fullfile('FFX', mod_name);
VOI_name = 'Thalamus';
% TODO load SPM to get dimensions of initial image
analyze_dir = 'D:\data\211214_additionalPilots\fmri\data4analyses\sub-01\FFX\FFX_model_block_30s';
SPM = spm_select('FPListRec',analyze_dir,'SPM.mat$'); load(SPM);
x_dim = SPM.xVol.DIM(1);
y_dim = SPM.xVol.DIM(2);
z_dim = SPM.xVol.DIM(3);
 
% ------------------------------------------------------
% Load data for all subjects
% ------------------------------------------------------

% Initialise variables
supraFrames=[];subjInfo=[];timeInfo=[];taskInfo=[]; seedSignLabels=[]; 

% Loop through subjects' data to make super subject. The following
% information must have been saved for each subject: the suprathreshold
% frames (selected based on seed's most (de)active moments after whole brain
% deconvolution), the indeces of each suprathreshold frame within the 
% original time course (timeInfo), and task labels for each suprathreshold 
% frame indicating which task was being performed at that point in time 
% (taskInfo); 
for i=1:length(sujets)
    
    
     % Load subject's variables, paths, etc
   
    thisSubject = initialize_vars(i, pilot_name);
    analyze_dir = fullfile(thisSubject.dataDir, mod_path);

    
    % Select suprathresholdFrames
    [PPIframes, suprathresholdFramesIx, correspondingTask4Frames] = thresholdFrames(i, pilot_name, 'PCC', 60);
    save(fullfile(analyze_dir, ['PPIframes_sub-' num2str(i, '%02d') '.mat']), 'PPIframes', '-v7.3');
    save(fullfile(analyze_dir, ['suprathresholdFramesIx_sub-' num2str(i, '%02d') '.mat']), 'suprathresholdFramesIx', '-v7.3');
    save(fullfile(analyze_dir, ['correspondingTask4Frames-' num2str(i, '%02d') '.mat']), 'correspondingTask4Frames', '-v7.3');

%     load([thisSubject.dataDir thisSubject.curSubj '.mat']); % loads PPIframes, suprathresholdFramesIx and correspondingTask4Frames for each subject
    
   load(fullfile(analyze_dir, ['PPIframes_sub-' num2str(i, '%02d') '.mat']));
    load(fullfile(analyze_dir, ['suprathresholdFramesIx_sub-' num2str(i, '%02d') '.mat']));
    load(fullfile(analyze_dir, ['correspondingTask4Frames-' num2str(i, '%02d') '.mat']));

    
    supraFrames      = [supraFrames, PPIframes']; % suprathreshold frames from each subject in the shape #voxels x #frames
    subjInfo         = [subjInfo, ones(1,size( PPIframes',2))*i]; % subject indeces corresponding to each supraframe
    timeInfo         = [timeInfo, suprathresholdFramesIx]; % each frame's corresponding index within a subject's timecourse
    taskInfo         = [taskInfo, correspondingTask4Frames]; % task each frame belongs to
    
    seedSignLabels = [seedSignLabels suprathresholdFramesSeedSign];
    
end

% Find voxels in and out of Grey Matter
[voxels, nonVoxels] = findGMvoxels;

% Find matrix rotations for saving nifti files
[voxel_size, voxel_shift] = prepareToSaveNifti(analyze_dir, VOI_name);


% ------------------------------------------------------
% Cluster timepoints, using frames as feature vectors
% ------------------------------------------------------

% Choose number of clusters (PPI-CAPs)
k = 4;

% Run kmeans++ on 
dist = 'cosine';
[PPI_CAPs, PPI_CAPs_Ix, ~, ~, ~]  = kmeanspp(supraFrames(voxels,:), k, dist, 100); 

for i = 1:k
        CAPix{i}              = find(PPI_CAPs == i); % #frames x 1; index of frames for this CAP
        CAPl{i}               = PPI_CAPs == i;  % boolean mask of frames for this CAP
        scf{i}                = CAPl{i}.* PPI_CAPs_Ix; %signed cluster's frames
        cluster{i}            = supraFrames(:,CAPix{i}); % #voxels x #frames
        CAPmean{i}               = (sum(supraFrames(:,scf{i}==1),2) - sum(supraFrames(:,scf{i}==2),2))/sum(CAPl{i} ); % whole-brain centroid (average)
        CAPmean{i}(nonVoxels,:)  = 0;
        
        polarityLabelsTmp= scf{i};
        polarityLabels{i} = polarityLabelsTmp(CAPix{i});
        taskLabels{i} = taskInfo(CAPix{i});
        seedSigns{i}   = seedSignLabels(CAPix{i});
        
        
        % Temporal normalization, as done by Karahanoglu et al., 2015
        if strcmp(dist, 'cosine')
            sd1 = std(cluster{i}, 0, 2);
            CAPmean{i}(sd1(:)>0) =  bsxfun(@rdivide,  CAPmean{i}(sd1(:)>0), sd1(sd1(:)>0));
        end
        
        % Spatial nomalization (voxelwise, for visualisation purposes only)
        sd2 = std(CAPmean{i}(CAPmean{i}~=0));
        CAPmean{i}(CAPmean{i}~=0)       = (CAPmean{i}(CAPmean{i}~=0) - mean(CAPmean{i}(CAPmean{i}~=0)))/sd2;
        
        % Reshape back to 3D
        CAP_vol{i}     = reshape(CAPmean{i} , 53, 63, 52); % save it back in the original 3D shape
        
        
        % Creates a new NIFTI for the CAP
        tmp_cap = make_nii(CAP_vol{i},voxel_size,round(-voxel_shift./voxel_size));
        tmp_cap.hdr.dime.datatype=64;
        tmp_cap.hdr.dime.bitpix=64;
        
        % Saves the cap nifti
        capNiftiFile = fullfile(ppicapsDir,['ppicaps_2Level_cap' num2str(i) '_' date '.nii']);
        save_nii(tmp_cap,capNiftiFile);
      
end

end


% From here on, statistical analysis is performed using the same code as in
% that used in ppicaps_toyExample.m for the statistical analysis part.




% ------------------------------------------------------
% HELPER FUNCTIONS
% ------------------------------------------------------

% Finds indeces of voxels in or out of Grey Matter (GM) 
function [inGM, nonGM] = findGMvoxels

mask      = load([pwd '/dependencies/MNIBrainMask.mat']); % Load MNIBrainMask
nonGM     = find(mask.MNIBrainMask(:) < 1); 
inGM  = find(mask.MNIBrainMask(:) == 1);

end

function [voxel_size, voxel_shift] = prepareToSaveNifti

VOI = spm_select('FPListRec',analyze_dir,[VOI_name '_1.mat$']); 
VOI = load(VOI);

voxel_size = diag(VOI.xY.spec.mat);
voxel_size = voxel_size(1:end-1)';

voxel_shift = VOI.xY.spec.mat(:,4);
voxel_shift = voxel_shift(1:end-1)';

end
