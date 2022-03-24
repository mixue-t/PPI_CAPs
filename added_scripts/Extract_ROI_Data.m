function ROI_data = Extract_ROI_Data(ROI, Contrast)

    Y = spm_read_vols(spm_vol(ROI),1);
    indx = find(Y>0);
    [x,y,z] = ind2sub(size(Y),indx);

    XYZ = [x y z]';

    ROI_data = nanmean(spm_get_data(Contrast, XYZ),2)

end