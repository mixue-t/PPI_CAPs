function step_2_main_threshold_frames_sl
% -------------------------------------------
%    Updated & modified by: Serafeim Loukas on Sep 10, 2021
%    For PPI-CAPs project with Joana (& Lorena)
% -------------------------------------------
%    step_2_main_threshold_frames_sl - function to execute ONLY the thresholding workflow
%
%    BUG FIX 1: misplaced parenthesis causing vector expansion (line 90)
%

clc;
clear all;

%% General paths setup (only the main path "ppicapsDir" needs modification for this function to work smoothly)
%  Set main path (universal path)
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

% Go to subject folder
addpath(genpath('/home/loukas/spm12/'));
addpath(genpath('C:\Users\tanm\Documents\stealth_game\analysis_git\PPI_CAPs'));
mainPath = [ppicapsDir, '/../data4analyses/'];
cd(mainPath)

S = dir([mainPath '/sub*']);
dirFlags = [S.isdir];
subFolders = S(dirFlags);

%% T2 template
%T2_template_file = [ppicapsDir, 'template/', 'template_T2_fused40w_n20.nii'];

%% Create cell with the subject folder names
n_subjects= size(S,1);
subject_name= {};
for k = 1 : n_subjects-1
    subject_name{k} = subFolders(k).name;
end
 
%% Initialise variables
seed = 'Thalamus';
percent = [5, 50, 60]; % desired percentages of percentage of frames to retain
percent = [50];
%% Load data (if exist) instead of running the loop just below this line
%load([ppicapsDir '/PPIallInfo_5.mat']);

%% Loop through subjects' data to make super subject. The following
% information must have been saved for each subject: the suprathreshold
% frames (selected based on seed's most (de)active moments after whole brain
% deconvolution), the indeces of each suprathreshold frame within the 
% original time course (timeInfo), and task labels for each suprathreshold 
% frame indicating which task was being performed at that point in time 
% (taskInfo);

for perc = 1:length(percent)
    disp(' ');
    disp(['The percentage value is:' ' ' num2str(percent(perc))]);
    
    % Initialise variables, will be cleared for each unique desired percentage value
    PPIframesALL = []; thresholdedIxALL = [];task4framesALL = [];seedSignALL = [];
    groupLabelALL = []; % 1 PTC, 2 PTM, 3 T
    subjectLabelALL = [];
    
    % main subject loop
    for i = 1:n_subjects
        
        % initialize the paths for the specific/current subject
        thissubj = initialize_vars(subject_name{i}, 'AP', false);

        disp(thissubj.curSubj)
        deconvData = deconvolveSubj_parallel(thissubj); % deconvolution has already been performed -- no need to re-run
        
        % TODO what is seedSign (output)?
        % Perform the suprathresholding procedure for this subject and for the desired % of frames to retain
        [PPIframes, thresholdedIx, task4frames, seedSign] = thresholdFrames(thissubj.curSubj, 'AP', seed, percent(perc));

        disp("Estimating the suprathersholded frames for the current subject...");
        
        % Append results from each subject to build group-level suprathresholded matrix of frames
        PPIframesALL = [PPIframesALL PPIframes']; % suprathreshold frames from each subject in the shape #voxels x #frames
        thresholdedIxALL = [thresholdedIxALL thresholdedIx]; % each frame's corresponding index within a subject's timecourse
        task4framesALL = [task4framesALL task4frames]; % task each frame belongs to
        seedSignALL = [seedSignALL seedSign]; % seed signs
        
        % group & subject labels construction
        groupLabelALL = [groupLabelALL ones(1,length(thresholdedIx))];   % PTC
        subjectLabelALL = [subjectLabelALL ones(1,length(thresholdedIx)) * str2num(thissubj.curSubj(5:6))];
% 
%         if thissubj.curSubj(1) == 'T' % FT
%             groupLabelALL = [groupLabelALL ones(1,length(thresholdedIx))*3];
%             subjectLabelALL = [subjectLabelALL ones(1,length(thresholdedIx)) * str2num(thissubj.curSubj(5:6))];
%         else
%             subjectLabelALL = [subjectLabelALL ones(1,length(thresholdedIx)) * str2num(thissubj.curSubj(7:8))];
%             if thissubj.curSubj(3) == 'C'
%                 groupLabelALL = [groupLabelALL ones(1,length(thresholdedIx))];   % PTC
%             else
%                 groupLabelALL = [groupLabelALL ones(1,length(thresholdedIx))*2]; % PTM
%             end
%         end
    end
    disp(['Saving results for percentage value:' ' ' num2str(percent(perc))]);
    save([ppicapsDir 'PPIallInfo_' num2str(percent(perc)) '.mat'], 'PPIframesALL', 'thresholdedIxALL', 'task4framesALL', 'seedSignALL', 'groupLabelALL', 'subjectLabelALL', '-v7.3');
    disp(['Results were saved !']);
    disp(' ');
end
end
