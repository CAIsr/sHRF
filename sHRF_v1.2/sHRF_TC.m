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
sInput = str2mat ('Estimate HRF from SPM cluster time course');    
InputType = spm_input('The time course is defined as', '1', 'm',sInput);
%%
%-----------------------------------------------------------------------
%------User specify output and relavent parameters----------------------
%-----------------------------------------------------------------------
sInput = str2mat ('Select the time couse file');
InputType = spm_input('The time course is defined as', '1', 'm',sInput);
TC_mat = ...
spm_select(1, 'mat', ...
        'Select time course rendered by SPM');
load (TC_mat);
sInput = str2mat ('Specify an output directory');
InputType = spm_input('The output directory is specified', '1', 'm',sInput);
Out_dir = ...
    spm_select(1, 'dir', 'Select an output directory');
TR=spm_input('The repetition time of image:', 4);
Volume_start_time = ...
    spm_input('First volume acquisition time:', 5);
V_ns=spm_input('Specify fMRI volumes for baseline', 6);
%%
%-----------------------------------------------------------------------
%-----User specify stimulus function and construct stimulus function----
%-----------------------------------------------------------------------
sSti = str2mat(...
    'Stimulus of boxcar function',...
    'Stimulus of stick function');
sStiType = spm_input('Select stimulus function...', 7, 'm', sSti);
T=spm_input('The length of experiment:', 8);
Sti_onsets = spm_input('The stimuli onset time in seconds', 9);
num_sti = spm_input('The number of the stimuli per block',10);
if (sStiType == 1)
    boxtime = spm_input('The stimulus duration in seconds', 11);
    Sti_tc = sSti_box(Sti_onsets, boxtime, T, num_sti);
else
    Sti_tc = sSti_Stick(Sti_onsets,T,num_sti);
end
%%
%-----------------------------------------------------------------------
%-----Load fMRI signal change time course using SPM eingenvalues----
%-----------------------------------------------------------------------
m=size(xY.u, 1);
fMRI=xY.u';
% Approximate the first and last eigenvalue with the nearest neigbor
fMRI(1)=fMRI(2);
fMRI(m)=fMRI(m-1);
fMRI_times=(Volume_start_time:TR:TR*(m-1));
%%
%-----------------------------------------------------------------------
%-----Calculate fMRI signal change time course--------------------------
%-----------------------------------------------------------------------
%---------Determine th base intensity-----------------------------------
Sum_na=0;
Count_na=0;
for i=1:size(V_ns,2)
    Sum_na=Sum_na+fMRI(V_ns(i));
    Count_na=Count_na+1;
end
Base=Sum_na/Count_na;
fMRI=(fMRI-Base)*100/Base;
fMRI_tc=[fMRI_times;fMRI];
fMRI_tc=fMRI_tc';
figure('Units','normalized','Position',[1 0 1/2 1/3]);
plot(Sti_tc(:,1), Sti_tc(:,2)*1000, 'color', [.8,.8,.8]);
xlabel('Times');
ylabel('Signal changes in percentage');
hold on;
axis([0 T+10 min(fMRI)-0.2 max(fMRI)+0.2]);
title('fMRI signal changes');
plot(fMRI_tc(:,1), fMRI_tc(:,2), '--cs', 'LineWidth', 2,...
    'MarkerFaceColor','c', 'MarkerSize',3);
%%
%----------------------------------------------------------------------
%----Fitting the HRF---------------------------------------------------
%----------------------------------------------------------------------
[hrf, fit, e, param, aic, VM] = Fit_NL666(fMRI_tc,Sti_tc);
fit_times=(0:0.1:T);
HRF_fit=[fit_times;fit(1:(T*10+1))];
plot(HRF_fit(1,:),HRF_fit(2,:), 'm', 'LineWidth', 2);
hold off;
figure(2);
hold on;
xlabel('Peri-stimulus times');
ylabel('Signal changes in percentage');
title('Estimated HRF');
plot(hrf(:,1), hrf(:,2), 'LineWidth', 2);
hold off;
%%
%----------------------------------------------------------------------
%----Output the results------------------------------------------------
%----------------------------------------------------------------------
[f_token, f_remain] =strtok(TC_mat,'/');
while true
    [f_token, f_remain] =strtok(f_remain,'/');
    if isempty(f_remain), break; end
end
f_token=strtok(f_token,'.');
fout = strcat(Out_dir, f_token, '_HRF.txt');
fp = fopen(fout, 'w');
fprintf(fp, 'HRF VM: %f %f %f %f %f %f %f %f %f\n', VM);
fprintf(fp, '\n');
fprintf(fp, 'The fitting RSS and model AIC: %f %f\n', e, aic);
fprintf(fp,...
    'HRF parameter of height delay width onsets AUC: %f %f %f %f %f',...
    param);
fprintf(fp, '\n');

