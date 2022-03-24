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

% extract values from a given session
SPM.xX.xKXs.X = SPM.xX.xKXs.X(1:N,:);
SPM.xX.xKXs.u = SPM.xX.xKXs.u(1:N,:);
SPM_xX_K_sess = SPM.xX.K(i_sess); 

% avoid broadcast variables
SPM_Vbeta = SPM.Vbeta;
SPM_xCon_xY_Ic = SPM.xCon(xY.Ic);
SPM_xX_xKXs = SPM.xX.xKXs;

deconvData = cell(x,1);
parfor ix = 1:x
    disp (['ix: ', int2str(ix)]) 
    deconvData_y = cell(10,1);
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
    deconvData{ix,1} = deconvData_y;
end
toc
