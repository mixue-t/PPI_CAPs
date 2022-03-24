% Function based on the relevant section of SPM8's spm_peb_ppi.m 
% function for devonvolution, usually used in the context of PPI analyses.
%
% Last checked: September 2019

function deconvData = deconvolveSubj_parallel(subj)


% Load SPM
%--------------------------------------------------------------------------
% analyze_dir = [subj.dataDir 'GLM_SciFun/'];
analyze_dir= [subj.root subj.pilot '\fmri\data4analyses\' subj.curSubj '\FFX\FFX_model_block_30s'];
SPM = spm_select('FPListRec',analyze_dir,'SPM.mat$'); 
load(SPM);
tempVOI = load([subj.root subj.pilot '\fmri\ppi_caps\VOI_Thalamus_1.mat']);
i_sess = 1; 

xY = tempVOI.xY;


% Setup variables
%--------------------------------------------------------------------------
RT      = SPM.xY.RT;
dt      = SPM.xBF.dt;
NT      = round(RT/dt);
% N       = length(xY(1).u);
N       = SPM.nscan(i_sess);
k       = 1:NT:N*NT;
Sess    = SPM.Sess(i_sess);
VI =  SPM.xY.VY(1:N); % mapped image volumes (.img file info)

% Create basis functions and hrf in scan time and microtime
%--------------------------------------------------------------------------
hrf = spm_hrf(dt); % 1 corresponds to the microtime resolution I want


% Create convolved explanatory {Hxb} variables in scan time
%--------------------------------------------------------------------------
xb  = spm_dctmtx(N*NT + 128,N); %xb = DCT matrix
Hxb = zeros(N,N);
for i = 1:N
    Hx       = conv(xb(:,i),hrf);
    Hxb(:,i) = Hx(k + 128);
end
xb = xb(129:end,:); % TODO why?

% Get confounds (in scan time) and constant term
%--------------------------------------------------------------------------
X0 = xY(1).X0;
M  = size(X0,2);

% Specify covariance components; assume neuronal response is white
% treating confounds as fixed effects
%--------------------------------------------------------------------------
Q = speye(N,N)*N/trace(Hxb'*Hxb);
Q = blkdiag(Q, speye(M,M)*1e6);


% Get whitening matrix (NB: confounds have already been whitened)
% used to give weighted least squares estimates (WLS). 
%--------------------------------------------------------------------------
W = SPM.xX.W(Sess.row,Sess.row);


% Create structure for spm_PEB
%--------------------------------------------------------------------------
clear P
P{1}.X = [W*Hxb X0];        % Design matrix for lowest level
P{1}.C = speye(N,N)/4;      % i.i.d assumptions
P{2}.X = sparse(N + M,1);   % Design matrix for parameters (0's)
P{2}.C = Q;

% Deconvolution
%=========================================================================
tic
% Iterate over voxels and deconvolve
% % counter = 0; percentage = 0;
x= SPM.xY.VY(1).dim(1);
y= SPM.xY.VY(1).dim(2);
z= SPM.xY.VY(1).dim(3);
t = SPM.nscan(i_sess);
% % total = 153594;
% total = x*y*z;
% % t = SPM.nscan;
% deconvData = zeros(x,y,z,t);
% deconvData_x = zeros(x);
% deconvData_y = zeros(y,z,t);


% cd into analyze_dir to be able to calculate beta
cd(analyze_dir);
% matrix of trials and betas (columns) in each trial
SPM.xX.xKXs.X = SPM.xX.xKXs.X(1:N,:);
% u as in X u*diag(ds)*v'
SPM.xX.xKXs.u = SPM.xX.xKXs.u(1:N,:);
% v as in X u*diag(ds)*v'
SPM_xX_K_sess = SPM.xX.K(i_sess);

% avoid broadcast variables
% struct array of beta image handles
SPM_Vbeta = SPM.Vbeta;
% Contrast definitions structure array 
SPM_xCon_xY_Ic = SPM.xCon(xY.Ic);
% low frequency confound: high-pass cutoff (secs)
SPM_xX_xKXs = SPM.xX.xKXs;


deconvData_raw = cell(x,1);
parfor ix = 1:x
    disp (['ix: ', int2str(ix)]) 
    deconvData_y = cell(y,1);
    for iy = 1:y
        disp (['iy: ', int2str(iy)])
        deconvData_z = zeros(z,t);
        for iz = 1:z
            thisVoxel = spm_data_read(VI,'xyz',[ix;iy;iz]);
            thisVoxel = spm_filter(SPM_xX_K_sess,W*thisVoxel);
            
            % Remove null space of contrast
            beta  = spm_data_read(SPM_Vbeta,'xyz',[ix;iy;iz]);
            if ~isnan(mean(beta)) 
                thisVoxel = thisVoxel - spm_FcUtil('Y0',SPM_xCon_xY_Ic,SPM_xX_xKXs,beta);
            end
            % Simple deconvolution
            % -------------------------------------------------------------
            C  = spm_PEB(thisVoxel,P);
            xn = xb*C{2}.E(1:N);
            xn = spm_detrend(xn);
            
            % Save variables (NOTE: xn is in microtime and does not account for
            % slice timing shifts). To convert to BOLD signal convolve with a hrf.
            % Use a microtime to scan time index to convert to scan time: e.g.,
            % k = 1:NT:N*NT; where NT = number of bins per TR = TR/dt or SPM.xBF.T
            % and N = number of scans in the session. Finally account for slice
            % timing effects by shifting the index accordingly.
            %----------------------------------------------------------------------
            deconvVoxel = downsample(xn,16); 
            deconvData_z(iz,:) = deconvVoxel;
%             deconvData(ix, iy, iz,:) = deconvVoxel;
            
%             counter = counter +1;
%             percentage = (counter / total ) * 100;
%            
        end
        deconvData_y{iy,1} = deconvData_z;
    end
    deconvData_raw{ix,1} = deconvData_y;
end
toc

% CD back into working folder to save file
cd(pwd);
save([subj.root subj.pilot '/ppi_caps/' subj.curSubj '/' subj.curSubj '_deconvolvedData_raw.mat'], 'deconvData_raw', '-v7.3');

% reshape 
deconvData = zeros(x,y,z,t);
for ix = 1:x
    for iy = 1:y
    deconvData(ix,iy,:,:) = deconvData_raw{ix,1}{iy,1};
    end
end

save([subj.root subj.pilot '/ppi_caps/' subj.curSubj '/' subj.curSubj '_deconvolvedData.mat'], 'deconvData', '-v7.3');
% save(['sub-01_deconvolvedData_final.mat'], 'deconvData', '-v7.3');


end









