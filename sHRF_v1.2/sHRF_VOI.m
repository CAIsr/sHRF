%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% 
% sHRF: supplementary hemodynamic response function (HRF) modeling for SPM.
%
% sHRF_TC:Estimate HRF from the fMRI time course generated from the SPM
%
% sHRF_VOI:Estimate HRF from the SPM{T} map and ROI definition images
%
% Author: Zuyao Shan
% Create: July 2011
% Last modified: Sep 25 2012
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
clear all
close all
%--Set Up Interface--------------------------------------------------------
Finter = spm_figure('FindWin', 'Interactive');
set(Finter, 'Name', 'sHRF Set Up');
sInput = str2mat ('Estimate HRF from SPMt map and VOI images');    
InputType = spm_input('The time course is defined as', '1', 'm',sInput);
%%
%-----------------------------------------------------------------------
%------User specify output and relavent parameters----------------------
%-----------------------------------------------------------------------
%---The output directory------------------------------------------------
sInput = str2mat ('Specify an output directory');
InputType = spm_input('The output directory is specified', '1', 'm',sInput);
Out_dir = spm_select(1, 'dir', 'Select an output directory');
%---The VOI image-------------------------------------------------------
sInput = str2mat ('Select the VOI image');
InputType = spm_input('', '1', 'm',sInput);
Mask_Img = spm_select(1, 'image', 'Select the VOI defining image');
%---The group activation map--------------------------------------------
sInput = str2mat ('Select the group activation map');
InputType = spm_input('', '1', 'm',sInput);
Common_Img = spm_select(1, 'image', 'Select the activation map');
%---The individual SPMt map---------------------------------------------
sInput = str2mat ('Select the individual SPMt map');
InputType = spm_input('', '1', 'm',sInput);
SPMt_Img = spm_select(1, 'image', 'Select the individual SPMt map');
%---The fMRI volumes----------------------------------------------------
sInput = str2mat ('Select the all the fMRI volumes in time order');
InputType = spm_input('', '1', 'm',sInput);
Img = spm_select(Inf, 'image', 'Select all the fMRI images in time order');
%---Processing parameters-----------------------------------------------
thres = spm_input('Group activation map threshold:', 6);
LPF = spm_input('Low pass filter in seconds:', 7);
tau = spm_input('Reliable group activation percentage:', 8);
%---Experimental designs------------------------------------------------
TR=spm_input('The repetition time of image:', 9);
Volume_start_time = ...
    spm_input('First volume acquisition time:', 10);
V_ns=spm_input('Specify fMRI volumes for baseline', 11);
T = spm_input('The length of experiment:', 12);
Sti_onsets = spm_input('The stimuli onset time in seconds', 13);
num_sti = spm_input('The number of the stimuli per block',14);
%%
%-----------------------------------------------------------------------
%-----User specify stimulus function and construct stimulus function----
%-----------------------------------------------------------------------
sSti = str2mat(...
    'Stimulus of boxcar function',...
    'Stimulus of stick function');
sStiType = spm_input('Select stimulus function...', 14, 'm', sSti);
if (sStiType == 1)
    boxtime = spm_input('The stimulus duration in seconds', 15);
    Sti_tc = sSti_box(Sti_onsets, boxtime, T, num_sti);
else
    Sti_tc = sSti_Stick(Sti_onsets,T,num_sti);
end
%%
%----------------------------------------------------------------------
%----------Define spatial coordinates according to VOI and group SPMt--
%----------------------------------------------------------------------
V_M = spm_vol(Mask_Img);
V_C = spm_vol(Common_Img);
V_S = spm_vol(SPMt_Img);
n=1;
for kk= 1:V_M.dim(3)
    for jj=1:V_M.dim(2)
        for ii=1:V_M.dim(1)
            Cord=[ii jj kk]';
            if (spm_get_data(V_M,Cord)) > 0
                if (spm_get_data(V_C,Cord) - thres) > 0
                    Locations(n,:)=[Cord', spm_get_data(V_M,Cord),...
                        spm_get_data(V_S, Cord)];
                    n=n+1;
                end
            end
        end
    end
end
roi_ID=unique(Locations(:,4));
%%
%----------------------------------------------------------------------
%-----------Set up output environment----------------------------------
%----------------------------------------------------------------------
fout = strcat(Out_dir, 'HRF.txt');
fp = fopen(fout, 'w');
%%
%----------------------------------------------------------------------
%---------Define the local maxima cluster------------------------------
for kk=(1:1:size(roi_ID,1))
    Common_Activated = Locations(Locations(:,4)==roi_ID(kk),:);
    Common_Activated = sortrows(Common_Activated, -5);
    roi_size = floor(size(Common_Activated,1)*tau);
    roi = Common_Activated(1:roi_size,:);
    %---Extract the time course of signal from fMRIs
    V_fMRI=spm_vol(Img);
    m=numel(V_fMRI);
    fMRI=zeros(1,m);
    fMRI_times=zeros(1,m);
    for ii=1:m
        for jj=(1:1:size(roi,1))
            fMRI(ii) = spm_get_data(V_fMRI(ii),...
                (roi(jj,1:3))')/size(roi,1) + fMRI(ii);
        end
        fMRI_times(ii) = Volume_start_time + (ii-1)*TR;
    end
    %----------Lowpass filter--------------
    K.RT=TR;
    K.row=1:m;
    K.HParam = LPF;
    fMRI=spm_filter(K,fMRI');
    fMRI=fMRI';
    %---------Determine th base intensity--------
    Sum_na=0;
    Count_na=0;
    for ii=1:size(V_ns,2)
        Sum_na=Sum_na+fMRI(V_ns(ii));
        Count_na=Count_na+1;
    end
    Base=Sum_na/Count_na;
    %--------fMRI signal changes----------------
    fMRI=(fMRI-Base)*100/Base;
    fMRI_tc=[fMRI_times;fMRI];
    fMRI_tc=fMRI_tc';
    %-------Plot fMRI signal changes-----------
    figure('Units','normalized','Position',[1 0 1/2 1/3]);
    plot(Sti_tc(:,1), Sti_tc(:,2)*1000, 'color', [.8,.8,.8]);
    xlabel('Times');
    ylabel('Signal changes in percentage');
    hold on;
    axis([0 T+10 min(fMRI)-0.2 max(fMRI)+0.2]);
    SName = num2str(roi_ID(kk),'Structure %d');
    Str = strcat(SName, ' fMRI signal changes');
    title(Str);
    plot(fMRI_tc(:,1), fMRI_tc(:,2), '--cs', 'LineWidth', 2,...
        'MarkerFaceColor','c', 'MarkerSize',3);
    %-------Estimate HRF--------------------------
    [hrf, fit, e, param, aic, VM] = Fit_NL666(fMRI_tc,Sti_tc);
    fit_times=(0:0.1:T);
    HRF_fit=[fit_times;fit(1:(T*10+1))];
    plot(HRF_fit(1,:),HRF_fit(2,:), 'm', 'LineWidth', 2);
    hold off;
    figure;
    hold on;
    xlabel('Peri-stimulus times');
    ylabel('Signal changes in percentage');
    Str = strcat(SName, ' estimated HRF');
    title(Str);
    plot(hrf(:,1), hrf(:,2), 'LineWidth', 2);
    hold off;
    %----Output the results---------------------------------
    fprintf(fp, 'Strcture ID: %d \n',roi_ID(kk));
    fprintf(fp, 'Cluster size: %d \n', size(roi,1));
    fprintf(fp, 'HRF VM: %f %f %f %f %f %f %f %f %f\n', VM);
    fprintf(fp, 'Simplfied HRF parameter of height delay width onsets AUC: %f %f %f %f %f \n', param);
    fprintf(fp, 'The RSS and AIC: %f %f\n', e, aic);
    fprintf(fp, '\n');
end
    

