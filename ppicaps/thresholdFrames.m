% -------------------------------------------------------
%
%    thresholdFrames - function to threshold frames
%
%    Created by:         Lorena Freitas
%    Last checked:       28.09.2019
%
%
% ------------------------------------------------------

function [PPIframes, thresholdedIx, task4frames, seedSign] = thresholdFrames(thisSubject, pilot, seed, percent)

% ----------------------
% Variables' set up
% ----------------------

% SPM's seed name
VOI_Name    = ['VOI_' seed '_1'];
i_sess = 1;

% Folder where analysis files were stored
b           = initialize_vars(thisSubject, pilot);
analyze_dir = [b.dataDir 'analysis/'];
analyze_dir = 'D:\data\211214_additionalPilots\fmri\data4analyses\sub-01\FFX\FFX_model_block_30s';

% Total number of frames
nFrames     = 887;
nScans = nFrames;

% Load subject's data (deconvData)
load(fullfile(analyze_dir, [b.curSubj '_deconvolvedData.mat']));
deconvData = de; 
% Reshape matrix into 2D to facilitate computations
subjectData2D      = zscore(reshape(deconvData, [], nFrames), 0, 2)'; % z-score each voxel across frames/time


% ----------------------
% Create seed timecourse
% ----------------------

% Obtain seed
VOI         = spm_select('FPListRec',analyze_dir,[VOI_Name '.mat$']); VOI = load(VOI);
VOIxyzMNI   = VOI.xY.XYZmm; % Dimensions: [3 x #voxels double]
VOIxyzMat   = round(VOI.xY.spec.mat \ [VOIxyzMNI; ones(1, length(VOI.xY.XYZmm))]);
VOIxyzMat   = VOIxyzMat(1:3,:); % XYZ indeces of seed's voxels in "MATLAB matrix" space
clear VOI_Name VOI VOIxyzMNI

% Average the seed signal from all voxels
meanSeedSignal = averageSignal(deconvData, VOIxyzMat, nFrames);
% TODO: verify
seedSign = sign(meanSeedSignal);

% ----------------------
% Create task timecourse
% ----------------------

% Repetition time
TR_MW     = 0.72;

% Onset times for task A (fun scences), in seconds
funcFolder = 'D:\data\211214_additionalPilots\fmri\data4analyses\sub-01\func\model_block_30s';
onsetList = dir(fullfile(funcFolder, '*_events.tsv'));
onsetInfo = readtable(fullfile(funcFolder, onsetList(i_sess).name),'FileType','text','Delimiter','\t','TreatAsEmpty',{'N/A','n/a'});
ind_cond = find(cellfun(@(x) strcmp(x, '4'), string(onsetInfo.trial_type)));
onset_cond4 = onsetInfo.onset(ind_cond);
ind_cond = find(cellfun(@(x) strcmp(x, '1'), string(onsetInfo.trial_type)));
onset_cond1 = onsetInfo.onset(ind_cond);

% Duration of blockes. First row: dur_fun; Second row: dur_sci
dur_MW = [repmat(30, 1, 3); repmat(30, 1, 3)];

% Get task labels for all frames
[~, taskLabels]    = createTaskRegressor(onset_cond4, onset_cond1, dur_MW, TR_MW, nScans);

% Find index of frames to threshold
[thresholdedIx, ~] = thresholdSeed(meanSeedSignal, nFrames, percent);

% Find suprathreshold frames to keep
PPIframes     = subjectData2D(thresholdedIx,:);

% Find task labels for suprathreshold frames
task4frames   = taskLabels(thresholdedIx);

end





% ----------------
% HELPER FUNCTIONS
% ----------------


function meanSignal = averageSignal(X, VOIxyzMat, no_timepoints)

% Initialise variables
nVoxels = size(VOIxyzMat,2);
signal  = zeros(nVoxels, no_timepoints);
%X       = zscore(X);
X       = zscore(X, 0, 4); % z-score each voxel across time !

% Get signal from voxels
for i=1:nVoxels
    x = VOIxyzMat(1,i); y = VOIxyzMat(2,i); z = VOIxyzMat(3,i);
    signal(i,:) = squeeze(X(x,y,z,:))';
end

% Average signal from voxels
meanSignal = mean(signal); % meanSeedSignal = [1 x # nscan]

end


function [taskRegressor, taskLabels] = createTaskRegressor(onset_cond4, onset_cond1, dur_MW, TR_MW, nScans)

scanCond4 = ceil((onset_cond4/TR_MW)); % round up (?)
taskRegressor1    = zeros(1,nScans);

scanCond1 = ceil((onset_cond1/TR_MW)); % round up (?)
taskRegressor2   = zeros(1,nScans);

for i=1:length(scanCond4)
    thisOnset_4    = scanCond4(i);
    thisOnset_1    = scanCond1(i);
%     thisDuration = round(dur_MW(1,i)/TR_MW);
    thisDuration = ceil(30/TR_MW);
    taskRegressor1(thisOnset_4:(thisOnset_4+thisDuration-1)) = 1;
    taskRegressor2(thisOnset_1:(thisOnset_1+thisDuration-1)) = 1;
end

% % Flip zeros and ones
% taskRegressor2 = 1 - taskRegressor1; 
% 
% % After the 178th scan it's just resting state
% taskRegressor2([1,2,175:end]) = 0; 
taskRegressor = taskRegressor1;

taskLabels = taskRegressor;
taskLabels(taskRegressor2==1) = 2;
 % create task regressor with 1 / -1 values
taskRegressor(taskRegressor2==1) = -1;

end



function [thresholdedIx, threshold] = thresholdSeed(meanSeedSignal, nFrames, percent)

nFrames2keep              = floor(nFrames*percent / 100);
[sortedSeedS,sortingIxS]  = sort(meanSeedSignal,'descend');
thresholdedIx             = sortingIxS(1:ceil(nFrames2keep/2));
threshold                 = min(meanSeedSignal(thresholdedIx));

% Now take also moments where the seed is deactive
thresholdedIx            = sort([thresholdedIx sortingIxS(end-ceil(nFrames2keep/2):end)]);

end