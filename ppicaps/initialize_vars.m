% This function initialises the subject variables for the batch analysis of
% the Movie Watching data.
% _________________________________________________________________________
% Last checked: September 2019

% Lorena Freitas
% $Id: initialize_vars.m 11 2016-12-01 16:26:24F Lorena $


function [b] = initialize_vars(subject, pilot, b_server)

% SPM info
%--------------------------------------------------------------------------
b.spmDir = fileparts(which('spm')); % path to SPM installation


% Directory information
%--------------------------------------------------------------------------
if strcmp(pilot, 'QV')
    b.pilot = '210707_questionnaireValidation';
elseif strcmp(pilot, 'AP')
    b.pilot = '211214_additionalPilots'; 
end

if b_server
    b.root = '/home/users/t/tanm/data/'
else
    b.root = 'D:\data\'
end

dataDir = [b.root, b.pilot, '\fmri\data4analyses\'];


% Subject information
%--------------------------------------------------------------------------
% b.group = group;
b.curSubj = num2str(subject,'%02.f') ; % subject code (e.g., 'sujet01')
b.dataDir = strcat(dataDir, b.curSubj, filesep); % make data directory subject-specific

% Data type to use (deconvolved = 1 on non-deconvolved = 0)
%--------------------------------------------------------------------------
b.deconvolveData = 1;

% Folder for functional files
%--------------------------------------------------------------------------
b.dir_func  = [b.dataDir, 'func' filesep];


end
