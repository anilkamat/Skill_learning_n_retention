%% setups
%################# REGIONAL LEVEL ANALYSIS: retention and transfer #########################
clc;clear all; close all;warning off;
% add paths for the data folder and support files for the toolbox

addpath(genpath('C:\Users\_Kamat_\Desktop\RPI\ResearchWork\Papers_\Effective_Connectivity\Suturing_LC_phase2\EEG_codes_results\Codes')) %Scripts
addpath 'C:\Users\_Kamat_\Desktop\RPI\ResearchWork\Papers_\Effective_Connectivity\GCCA'             % toolbox directory
addpath 'C:\Users\_Kamat_\Desktop\RPI\ResearchWork\Papers_\Effective_Connectivity\GCCA\utilities'   % toolbox utility directory
addpath 'C:\Users\_Kamat_\Desktop\RPI\ResearchWork\Papers_\Effective_Connectivity\bsmart'           % toolbox utility directory
addpath 'C:\Users\_Kamat_\Desktop\RPI\ResearchWork\Papers_\Effective_Connectivity\toolbox_original\mvgc_v1.0\demo' % toolbox directory
connection_GC = {'LPFC-->RPFC','LPFC-->LPMC','LPFC-->RPMC','LPFC-->SMA',...
    'RPFC-->LPFC','RPFC-->LPMC','RPFC-->RPMC','RPFC-->SMA',...
    'LPMC-->LPFC','LPMC-->RPFC','LPMC-->RPMC','LPMC-->SMA',...
    'RPMC-->LPFC','RPMC-->RPFC','RPMC-->LPMC','RPMC-->SMA',...
    'SMA-->LPFC','SMA-->RPFC','SMA-->LPMC','SMA-->RPMC'};
addpath(genpath('C:\Users\_Kamat_\Desktop\RPI\ResearchWork\Papers_\Effective_Connectivity\EEG\Phase2_data\LC_EEG_CSD\Retention_Transfer')) %dataset
load('C:\Users\_Kamat_\Desktop\RPI\ResearchWork\Papers_\Effective_Connectivity\EEG\Phase2_data\HTA\HTA_suturing_LC_retTransfer.mat'); % new HTA timestamps
phases_idx = [1 3 6 11 14];                                                                                              % phases of HTA
LPFC = []; RPFC = []; LSMA = []; RSMA = []; SMA = []; LPMC = []; RPMC = [];Index_temp = [];
Cnames = {'LPFC','RPFC','LPMC','RPMC','SMA'};
%% ############################################# GC Analysis  ########################################################%%%
fg = 0;         % Flag for loading 
fs = 250;       % sampling frequency 
if isfile('C:\Users\_Kamat_\Desktop\RPI\ResearchWork\Papers_\Effective_Connectivity\Suturing_LC_phase2\EEG_codes_results\Results\Ph2_LC_retTransf_eeg_sut_GC.mat')
    load('C:\Users\_Kamat_\Desktop\RPI\ResearchWork\Papers_\Effective_Connectivity\Suturing_LC_phase2\EEG_codes_results\Results\Ph2_LC_retTransf_eeg_sut_GC.mat');
    fg = 1; 
    fprintf('Loading the precomputed analysis, skiping recomputation....');
else
    nsub = 28; %       % total no.of subjects
    nT = 3; %3         % number of trials
    freq = [8:0.05:34];           % Frequency range in Alpha or beta range based on analysis
    r = 1;
    fileList = sort(getAllFiles('C:\Users\_Kamat_\Desktop\RPI\ResearchWork\Papers_\Effective_Connectivity\EEG\Phase2_data\LC_EEG_CSD\Retention_Transfer'));
    nF = numel(fileList);
    cF = 0;  % counts number of subjct matched in the folder
    tic
    %normal_performance = readcell('C:\Users\_Kamat_\Desktop\RPI\ResearchWork\Papers_\Effective_Connectivity\Stressor_arm\Datasets\Normal_Trial_Success-Performance-Phase 2_(1).xlsx');
    nST = 4;     % Number of subtasks
    for st = 1:nST              % loop over subtasks
        h = 1;  % counter of subjects each subtask
        for sub_N = 1:nsub %nD              % loop over subjects
            subject = sprintf('P%02dReten_Transf',sub_N);  %P02Reten_Transf.mat
            for f= 1:nF
                fileName = cell2mat(fileList(f));  % name of the file on the folder
                if strcmp(fileName(1:15),subject)== 1 && strcmp(fileName(16:19),'.mat')== 1 && isfield(HTA_suturing_LC_retTransfer,fileName(1:3))==1
        %            Index = find(contains(normal_performance(:,6),fileName(1:end-4)));
        %             class = cell2mat(normal_performance(Index,7));
                    fprintf('%s::: %s :::%d \n',fileName, subject);
                    cF = cF+1;
                    subj_data = importdata(fileName); % load signals
                    nT = size(HTA_suturing_LC_retTransfer.(fileName(1:3)),1);
                    %fprintf('%s found, #trials %d \n',fileName, nT);
                    data_Reg = get_LS_channelData(subj_data);  % aggregated regional signals
                    for t = 1:nT                % loop over Trials
                        ST = sprintf('ST%d',st);
%                         sst = HTA_suturing_LC_retTransfer.(fileName(1:3))(t,st);     %start of the subtask
%                         est = HTA_suturing_LC_retTransfer.(fileName(1:3))(t,st+1);     %end of the subtask
                        
                        sst = HTA_suturing_LC_retTransfer.(fileName(1:3))(t,phases_idx(st));     %start of the subtask
                        est = HTA_suturing_LC_retTransfer.(fileName(1:3))(t,phases_idx(st+1));     %end of the subtask
                        Lst = est-sst;                                      %Length of the subtask
                        if ~isnan(sst) && ~isnan(est) && Lst >=0
                            %fprintf('\n filename: %s, ST: %d sst:%d, est:%d, Lst:%d \n', fileName, st, sst, est, Lst);
                            Y_subtask = data_Reg(:,sst*250:est*250);          % with sampling frequency of 250Hz (sf of data processed by BU is 250)                    
                            Y_subtask_allCh = subj_data.data(:,sst*250:est*250);          % with sampling frequency of 250Hz (sf of data processed by BU is 250)                    
                            
                            %Write the eeg signals to csv file for
                            %non-linear granger causality analysis.
                            csvfilename2 = sprintf('C:\\Users\\_Kamat_\\Desktop\\RPI\\ResearchWork\\Papers_\\Effective_Connectivity\\Suturing_LC_phase2\\EEG_codes_results\\Results\\RetentionTransf_allChaEEG_Phases\\P%d\\%sTr_%d_alphafreqBP.csv',st,fileName(1:3),t); %sub_N,
                            Y_wind_allChan_alpha = bandpass(Y_subtask_allCh',[8 12],fs)';  % alpha band pass
                            csvwrite(csvfilename2,Y_wind_allChan_alpha)
                            
%                             NL = size(Y_subtask,2);
%                             [bic,aic] = cca_find_model_order(Y_subtask,2,9);
%                             mo = min(bic,aic);                                      % selection of model order
%                             [GW,COH,pp,waut,cons]= cca_pwcausal(Y_subtask,1,NL,mo,250,freq,1); 
%                             %GC_fqmean.ExpFLS.(sub_trial_name).(trial).(W)   = mean(GW(:,:,idx3:idx4),3);
%                             idx1 = find(freq == 8);                                  % mean in the neurophysiology frequency band.
%                             idx2 = find(freq == 12);
%                             GC_temp = mean(GW(:,:,idx1:idx2),3);
%                             sz = numel(GC_temp);
%                             Suturing_GC_retTransf.(ST).sub(h,:)  = [reshape(GC_temp,1,sz)];
                        end
                        h = h+1;  
                    end
                end
            end
        end
    end
    fprintf('%d Number of subjects in data folder \n',cF)
end
toc
%%  If the precomputed file is loaded, then no need to delet the diagonals
if fg == 0
    for j =1:size(fieldnames(Suturing_GC_retTransf))
        ST_temp = sprintf('ST%d',j);        % Subtask          
        temp = all(Suturing_GC_retTransf.(ST_temp).sub,2);
        zero_rows = find(temp ==0);
        Suturing_GC_retTransf.(ST_temp).sub(zero_rows, :) = [];
        Suturing_GC_retTransf.(ST_temp).sub(:,[1 7 13 19 25],:) = [];       % Remove all the diagonals
    end
    save('..\Results\Ph2_LC_retTransf_eeg_sut_GC.mat','Suturing_GC_retTransf')  % saves the connectivity to the local current directory if it doesn't exist
end
%%  ############ Box plot #############
fpath = '..\Results\Box_Plots_GC_retTransfer';
% plot for the subjects
for j =1:size(fieldnames(Suturing_GC_retTransf))
    ST = sprintf('ST%d',j);
    if isfield(Suturing_GC_retTransf.(ST),'sub') && size(Suturing_GC_retTransf.(ST).sub,1)>1 % check if the file exists and is greater than 1 subject
        close all;
        f = figure('visible','off');
        size(Suturing_GC_retTransf.(ST).sub(:,1:end-1))
        fprintf('%s %s',subject,ST)
        boxplot(Suturing_GC_retTransf.(ST).sub(:,1:end),'Labels',connection_GC, 'Colors','b','jitter',0, 'symbol', ''); % 'whisker', inf,  'PlotStyle','Compact'
        set(gca,'FontSize',10,'XTickLabelRotation',45)
        xlabel('Connection')
        ylabel('Granger Causality')
        ylim([-0.01 0.15])
        ttl = sprintf('Suturing LC-retTransfer: Subtask-%d',j);
        title(ttl);
        baseFileName = sprintf('SuturingLC_retenTransfer_ST%d.png',j);
        fullFileName = fullfile(fpath, baseFileName);
        saveas(f,fullFileName)
    end
end
%% Write connectivity  to csv for SVM analysis
for st =1:13  % loop ove subtasks  
    ST = sprintf('ST%d',st);
    fileName_ = sprintf('C:\\Users\\_Kamat_\\Desktop\\RPI\\ResearchWork\\Papers_\\Effective_Connectivity\\Suturing_LC_phase2\\EEG_codes_results\\Results\\Connectivities\\retenTransfer_suturingEEG_GC_ST%d.csv',st);
    if isfield(Suturing_GC_retTransf,ST)
        T = array2table(Suturing_GC_retTransf.(ST).sub);
        T.Properties.VariableNames(1:20) = connection_GC;
        writetable(T,fileName_)
    end
end
fprintf('\n writing to connectivites done!')
%%
function MatFileList = getAllFiles(dirName)
    dirInfo = dir(dirName);
    dirIndex = [dirInfo.isdir];
    fileList = {dirInfo(~dirIndex).name}';
    
    % Filter only files with ".mat" extension
    MatFileList = fileList(endsWith(fileList, '.mat'));
    
    subDirs = {dirInfo(dirIndex).name};
    validIndex = ~ismember(subDirs, {'.', '..'});
    
    for i = find(validIndex)
        nextDir = fullfile(dirName, subDirs{i});
        MatFileList = [MatFileList; getAllFiles(nextDir)];
    end
end
