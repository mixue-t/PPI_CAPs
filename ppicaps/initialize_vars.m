% This function initialises the subject variables for the batch analysis of
% the Movie Watching data.
% _________________________________________________________________________
% Last checked: September 2019

% Lorena Freitas
% $Id: initialize_vars.m 11 2016-12-01 16:26:24F Lorena $


function [b] = initialize_vars(subject,pilot)

% SPM info
%--------------------------------------------------------------------------
b.spmDir = fileparts(which('spm')); % path to SPM installation


% Directory information
%--------------------------------------------------------------------------
% dataDir = '<GENERAL PATH FOR ALL SUBJECTS DATA>';
if strcmp(pilot, 'QV')
    b.pilot = '210707_questionnaireValidation';
elseif strcmp(pilot, 'AP')
    b.pilot = '211214_additionalPilots'; 
end
dataDir = ['D:\data\', b.pilot, '\fmri\data4analyses'];

% Subject information
%--------------------------------------------------------------------------
% b.group = group;
b.curSubj = ['sub-' num2str(subject,'%02.f')] ; % subject code (e.g., 'sujet01')
b.dataDir = strcat(dataDir, filesep , b.curSubj, filesep); % make data directory subject-specific

% Data type to use (deconvolved = 1 on non-deconvolved = 0)
%--------------------------------------------------------------------------
b.deconvolveData = 1;

% Folder for functional files
%--------------------------------------------------------------------------
% b.dir_fonc = '<PATH TO FUNCTIONAL DATA, e.g., /warped/ >';
b.dir_fonc  = [b.dataDir, 'func'];

end
