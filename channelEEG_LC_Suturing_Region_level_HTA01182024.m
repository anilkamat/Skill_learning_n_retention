%% setups
%################# REGIONAL LEVEL ANALYSIS #########################
clc; clear all; close all; warning off;
% add paths for the data folder and support files for the toolbox
addpath(genpath('C:\Users\_Kamat_\Desktop\RPI\ResearchWork\Papers_\Effective_Connectivity\EEG\Phase2_data\LC_EEG_CSD')) %dataset
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
load('C:\Users\_Kamat_\Desktop\RPI\ResearchWork\Papers_\Effective_Connectivity\EEG\Phase2_data\HTA\HTA_suturing_LC.mat'); % new HTA timestamps
phases_idx = [1 3 6 11 14];  
%% ############################################# GC Analysis  ########################################################%%%
fg = 0;                     % Flag for loading 
fs = 250;       % sampling frequency 
if isfile('C:\Users\_Kamat_\Desktop\RPI\ResearchWork\Papers_\Effective_Connectivity\Suturing_LC_phase2\EEG_codes_results\Results\Ph2_LC_eeg_sut_GC.mat')
    load('C:\Users\_Kamat_\Desktop\RPI\ResearchWork\Papers_\Effective_Connectivity\Suturing_LC_phase2\EEG_codes_results\Results\Ph2_LC_eeg_sut_GC.mat');
    fg = 1; 
    fprintf('Loading the precomputed analysis, skiping recomputation....');
else
    nD = 15; %15;      % total no.of days
    freq = [8:0.05:34];           % Frequency range in Alpha or beta range based on analysis
    r = 1;
    %fileList = sort(getAllFiles('C:\Users\_Kamat_\Desktop\RPI\ResearchWork\Papers_\Effective_Connectivity\EEG\Phase2_data\LC_EEG_CSD\LC'));
    fileList = sort(getAllFiles('C:\Users\_Kamat_\Desktop\RPI\ResearchWork\Papers_\Effective_Connectivity\EEG\Phase2_data\LC_EEG_CSD\LC\P15'));
    nF = numel(fileList);
    cF = 0;  % counts number of subjct matched in the folder
    tic
    %for d = 1:nD              % loop over days for all subject
    for d = 13:13              % loop over days for individual day and subject
        h = 1;  % counter of subjects each day
        day = sprintf('Day_%d',d);
        subject = sprintf('D%02d',d);
        for f = 1:nF
            fileName = cell2mat(fileList(f));  % name of the file on the folder
            %if strcmp(fileName(4:6),subject)== 1 &&
            %strcmp(fileName(7:10),'.mat')== 1 &&
            %isfield(HTA_suturing_LC,fileName(1:6))==1 % for all subjects
            %at once
            if strcmp(fileName(4:6),subject)== 1 && strcmp(fileName(1:3),'P15')== 1 && strcmp(fileName(7:10),'.mat')== 1 && isfield(HTA_suturing_LC,fileName(1:6))==1 % for individual level of missing subjects
    %            Index = find(contains(normal_performance(:,6),fileName(1:end-4)));
    %             class = cell2mat(normal_performance(Index,7));
                %fprintf('%s::: %s :::%d \n',fileName, subject);
                cF = cF+1;
                subj_data = importdata(fileName); % load signals
                nT = size(HTA_suturing_LC.(fileName(1:6)),1);
                %fprintf('%s found, #trials %d \n',fileName, nT);
                data_Reg = get_LS_channelData(subj_data);  % aggregated regional signals
                for t = 1:nT                % loop over Trials
                    trial = sprintf('Tr_%d',t);
                    %fprintf('Trial# %d \n',t);
                    W = sprintf('Trial_%d',t);
                    s = 1;      % counter for successful subjects
                    us = 1;     % counter for unsuccessful subjects
                    as = 1;     % counter for all subjects
                    nST = 13;     % Number of subtasks
                    
                    for st = 1:nST              % loop over subtasks
                        ST = sprintf('ST%d',st);        
                        % For subtask level extraction
                        sst = HTA_suturing_LC.(fileName(1:6))(t,st);     %start of the subtask
                        est = HTA_suturing_LC.(fileName(1:6))(t,st+1);     %end of the subtask
                        % For phase level extraction
%                         sst = HTA_suturing_LC.(fileName(1:6))(t,phases_idx(st));     %start of the phase
%                         est = HTA_suturing_LC.(fileName(1:6))(t,phases_idx(st+1));     %end of the phase                      
                        
                        Lst = est-sst;                  %Length of the subtask
                        %fprintf('ST# %d, lst: %d \n',st,Lst);
                        if ~isnan(sst) && ~isnan(est) && Lst >=0
                            %fprintf('\n filename: %s, ST: %d sst:%d, est:%d, Lst:%d \n', fileName, st, sst, est, Lst);
                            Y_subtask = data_Reg(:,sst*250:est*250);          % with sampling frequency of 250Hz (sf of data processed by BU is 250)                    
                            
                            Y_subtask_allCh = subj_data.data(:,sst*250:est*250);          % with sampling frequency of 250Hz (sf of data processed by BU is 250)
                            %Write the subtask level eeg signals to csv file for
                            %non-linear granger causality analysis.
                            csvfilename2 = sprintf('C:\\Users\\_Kamat_\\Desktop\\RPI\\ResearchWork\\Papers_\\Effective_Connectivity\\Suturing_LC_phase2\\EEG_codes_results\\Results\\LC_allChaEEG\\D%d\\ST%d\\%sTr_%d_alphafreqBP.csv',d,st,fileName(1:3),t); %sub_N,
                            %csvfilename2 = sprintf('C:\\Users\\_Kamat_\\Desktop\\RPI\\ResearchWork\\Papers_\\Effective_Connectivity\\Suturing_LC_phase2\\EEG_codes_results\\Results\\Connectivities_LSTMED_multiChanROI_individualLevel\\D%dST%d%sTr_%d_alphafreqBP.csv',d,st,fileName(1:3),t); %sub_N,
                            Y_wind_allChan_alpha = bandpass(Y_subtask_allCh',[8 12],fs)';  % alpha band pass
                            csvwrite(csvfilename2,Y_wind_allChan_alpha)
                            
                            %Write phase level data
%                             csvfilename2 = sprintf('C:\\Users\\_Kamat_\\Desktop\\RPI\\ResearchWork\\Papers_\\Effective_Connectivity\\Suturing_LC_phase2\\EEG_codes_results\\Results\\LC_allChaEEG_Phase\\D%d\\P%d\\%sTr_%d_alphafreqBP.csv',d,st,fileName(1:3),t); %sub_N,
%                             Y_wind_allChan_alpha = bandpass(Y_subtask_allCh',[8 12],fs)';  % alpha band pass
%                             csvwrite(csvfilename2,Y_wind_allChan_alpha)
                            
%                             NL = size(Y_subtask,2);
%                             [bic,aic] = cca_find_model_order(Y_subtask,2,9);
%                             mo = min(bic,aic);                                      % selection of model order
%                             [GW,COH,pp,waut,cons]= cca_pwcausal(Y_subtask,1,NL,mo,250,freq,1);   % in fNIRS the low model order was selected to make the algorithm work
%                             %GC_fqmean.ExpFLS.(sub_trial_name).(trial).(W)   = mean(GW(:,:,idx3:idx4),3);
%                             idx1 = find(freq == 8);                                  % mean in the neurophysiology frequency band.
%                             idx2 = find(freq == 12);
%                             GC_temp = mean(GW(:,:,idx1:idx2),3);
%                             sz = numel(GC_temp);
%                             Suturing_GC.(day).(ST).sub(h,:)  = [reshape(GC_temp,1,sz)];
                        end                    
                    end
                    h = h+1;  
                end
            end
        end
    end
    fprintf('%d Number of subjects in data folder \n',cF)
end
toc
%%  If the precomputed file is loaded, then no need to delet the diagonals
if fg == 0
    for i = 1:size(fieldnames(Suturing_GC),1)
        d = sprintf('Day_%d',i);
        for j =1:size(fieldnames(Suturing_GC.(d)))
            ST_temp = sprintf('ST%d',j);        % Subtask          
            temp = all(Suturing_GC.(d).(ST_temp).sub,2);
            zero_rows = find(temp ==0);
            Suturing_GC.(d).(ST_temp).sub(zero_rows, :) = [];
            Suturing_GC.(d).(ST_temp).sub(:,[1 7 13 19 25],:) = [];       % Remove all the diagonals
        end
    end
    save('C:\Users\_Kamat_\Desktop\RPI\ResearchWork\Papers_\Effective_Connectivity\Suturing_LC_phase2\EEG_codes_results\Results\Ph2_LC_eeg_sut_GC.mat','Suturing_GC')  % saves the connectivity to the local current directory if it doesn't exist
    
    %C:\Users\_Kamat_\Desktop\RPI\ResearchWork\Papers_\Effective_Connectivity\Suturing_LC_phase2\EEG_codes_results\Results
end
%%  ############ Box plot #############
fpath = 'C:\Users\_Kamat_\Desktop\RPI\ResearchWork\Papers_\Effective_Connectivity\Suturing_LC_phase2\EEG_codes_results\Results\Box_Plots_GC_';
% plot for the subjects
for i = 1:size(fieldnames(Suturing_GC),1)
    d = sprintf('Day_%d',i);
    for j =1:size(fieldnames(Suturing_GC.(d)))
        ST = sprintf('ST%d',j);
        if isfield(Suturing_GC.(d).(ST),'sub') && size(Suturing_GC.(d).(ST).sub,1)>1 % check if the file exists and is greater than 1 subject
            close all;
            f = figure('visible','off');
            size(Suturing_GC.(d).(ST).sub(:,1:end-1))
            fprintf('%s %s',d,ST)
            boxplot(Suturing_GC.(d).(ST).sub(:,1:end),'Labels',connection_GC, 'Colors','b','jitter',0, 'symbol', ''); % 'whisker', inf,  'PlotStyle','Compact'
            set(gca,'FontSize',10,'XTickLabelRotation',45)
            xlabel('Connection')
            ylabel('Granger Causality')
            ylim([-0.01 0.1])
            ttl = sprintf('Suturing LC successful: Day-%d, Subtask-%d',i,j);
            title(ttl);
            baseFileName = sprintf('SuturingLC_succ_D%d_T%d.png',i,j);
            fullFileName = fullfile(fpath, baseFileName);
            saveas(f,fullFileName)
        end
    end
end
%% ############# Plot the means of all the trails ##############
% plot the boxplot successful
q = 1; % grand mean of all the trials
fpath = 'C:\Users\_Kamat_\Desktop\RPI\ResearchWork\Papers_\Effective_Connectivity\Suturing_LC_phase2\EEG_codes_results\Results\Mean_connectivity_';
for con= 1:numel(connection_GC)  % loop over connectivity
    close all
    p =1; % counter of total number of trials in entire experiment
    for i=1:size(fieldnames(Suturing_GC),1)
        d = sprintf('Day_%d',i);
        for j=1:size(fieldnames(Suturing_GC.(d)),1)
            t = sprintf('ST%d',j);
            if isfield(Suturing_GC.(d).(t),'sub') && size(Suturing_GC.(d).(t).sub,1) > 1
                M(p) = mean(Suturing_GC.(d).(t).sub(:,con));
                STD(p) = std(Suturing_GC.(d).(t).sub(:,con));
                GM(q) = M(p);
                p = p+1;
                q = q+1;
            end
        end
        mDay(i) = p-1 ;% marker for days
    end
    f = figure('visible','off','position',[10 10 700 400]);
    errorbar(M,STD,"-s","MarkerSize",10,"MarkerEdgeColor","blue","MarkerFaceColor",[0.65 0.85 0.90])
    ylim([-0.05 0.1])
    title('Connectivity: Day1')
    legend(connection_GC(con))
    xlabel('Day')
    ylabel('Connectivity')
    %vline2([mDay])
    %text([mDay(6:6)-1], 1*[1 ], { 'D6'})
    %text([mDay(7:7)-1], 1.1*[1], { 'D7'})
    %text([mDay(8:end)-3], 1*[ 1 1 1 1 1 1 1 1], { 'D8', 'D9', 'D10', 'D11', 'D12', 'D13', 'D14', 'D15'})
    baseFileName = sprintf('Conn%d_D1-D15.png',con);
    fullFileName = fullfile(fpath, baseFileName);
    saveas(f,fullFileName)
end
%% Histogram of all the connectivities
f = figure(); %('visible','off');
histogram(GM,30)
title('Distribution: connectivity')
xlabel('Granger causality')
ylabel('samples')
%% Write connectivity of sucessful and unsuccessful subjects together
GC_alDays = [];
for d = 1:15  % loop over days
    D = sprintf('Day_%d',d);
    for st =1:13  % loop ove trials, 3 trials taken after 3 trials the perdiction of the model has been found to be deteoriate 
        fileName_ = sprintf('C:\\Users\\_Kamat_\\Desktop\\RPI\\ResearchWork\\Papers_\\Effective_Connectivity\\Suturing_LC_phase2\\EEG_codes_results\\Results\\Connectivities_learningDays\\EEG_GC_D%dST%d_T_All.csv',d,st);
        GC_data1 = [];
        ST = sprintf('ST%d',st);
        if isfield(Suturing_GC.(D),ST)
            GC_data1 = [GC_data1 ; Suturing_GC.(D).(ST).sub];
        end
        T = array2table(GC_data1);
        T.Properties.VariableNames(1:20) = connection_GC;
        writetable(T,fileName_)
    end
end
fprintf('writing to connectivites done!')
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
