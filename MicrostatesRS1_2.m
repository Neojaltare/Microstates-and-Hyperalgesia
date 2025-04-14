
%% Microstate Analysis for Resting-State EEG
%
% This script performs microstate analysis on resting-state EEG data 
% collected during two sessions (RS1 and RS2) in a sample of participants.
% The analysis includes:
%
% - Loading and preprocessing EEG datasets using EEGLAB
% - Selecting GFP peaks for microstate segmentation
% - Performing clustering to identify microstate prototypes (k = 3–6)
% - Backfitting the identified prototypes to individual datasets
% - Computing statistics (GEV, duration, occurrence, etc.) for each microstate
% - Exporting the results to an Excel file
% - Optional visualization of prototype maps and segmented time series
%
% Note:
% - The script is tailored to the file structure and naming conventions of the original project
% - Paths must be updated to match the local directory structure
% - This code uses EEGLAB and the Microstate EEGLAB Toolbox
%
% Author: Ketan Jaltare
% Email: ketan.jaltare@kuleuven.be

%% Microstates analysis for all resting state data
clear;clc;

addpath('/Path to EEGLAB/eeglab2023.0')
addpath(genpath('/Matlab Code')) % This folder contains some useful functions (isnumber, etc)

% start EEGLAB to load all dependent paths
eeglab nogui
%% Establish participants to analyse
% Parameters for the analysis


RS = ["RS1","RS2"]; % RS1 or RS2
Sample = "Control"; % Define which subjects to analyse - All? Or just the Control group
Toexclude = [19,41,47,35]; % These are definitely to be excluded
suspect = [21,30,40,46,57]; % These are the participants where the behavioural data is suspect
Excludesus = 0;



HFSdat = readtable('/Path to HFS data/Exp3HFSdata.xlsx');
EEGdir = '/Path to EEG data/';
contents = dir(EEGdir);
dirFlags = [contents.isdir];
subFolders = contents(dirFlags);
subjectFolders = {subFolders(4:end).name};

if Sample == "Control"

    Controlsubs = HFSdat.participant(string(HFSdat.group) == 'C');
    subjectFolders = string(subjectFolders);
    Subs = strings(size(Controlsubs,1),1);
    ind = 1;
    for i = 1:length(subjectFolders)
        subjectno = subjectFolders(i);
        subjectnoidx = isnumber(subjectno{1});
        subjectno = double(string(subjectno{1}(subjectnoidx)));
        if ismember(subjectno,Controlsubs)
            Subs(ind) = string(subjectFolders(i));
            ind = ind + 1;
        end
    end

elseif Sample == "All"

    Subs = string(subjectFolders);

end

%% 3.3 Data selection and aggregation
%% 3.3.1 Loading datasets in EEGLAB

% Initialize table to store data
preprocvars = {'Subject','RS','Condition','Analysis','MS','TotalGEV','GFP','Duration','Occurance','Coverage',...
    'GEV','TP1','TP2','TP3','TP4','TP5'};
d = 1;
rec = "1";
Condition = "string";
Analysis = "string";
MS = "string";
TotalGEV = 1;
GFP = 1;
Duration = 1;
Occurance = 1;
Coverage = 1;
GEV = 1;
tp1 = 1;
tp2 = 1;
tp3 = 1;
tp4 = 1;
tp5 = 1;
MicrostateData  = table(d,rec,Condition,Analysis,MS,TotalGEV,GFP,Duration,Occurance,Coverage,GEV,...
    tp1,tp2,tp3,tp4,tp5,'VariableNames', preprocvars);


ind = 1;
for k = 1:2
    for i=1:length(Subs)
        subjectno = Subs(i);
        subjectnoidx = isnumber(subjectno{1});
        subjectno = double(string(subjectno{1}(subjectnoidx)));

        if ismember(subjectno,Toexclude)
            continue
        end

        if Excludesus == 1
            if ismember(subjectno,suspect)
                continue
            end
        end
        % retrieve a list of all EEG Files in EEGdir
        EEGFile = strcat('Cleaned',Subs(i),RS(k),'.set');
        EEG = pop_loadset('filename',char(EEGFile),'filepath',char(strcat(EEGdir,Subs(i),'/')));
        
        [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, ind );
        MicrostateData.Subject(ind) = subjectno;
        MicrostateData.RS(ind) = RS(k);
        MicrostateData.Condition(ind) = string(HFSdat.group(HFSdat.participant == subjectno));
        MicrostateData.Analysis(ind) = Sample;
        ind = ind+1;
        %eeglab redraw % updates EEGLAB datasets
    end
end
ALLEEG(ind:end) = [];
numsubs = length(ALLEEG);
%% 3.3.2 Select data for microstate analysis

%  % Use this for loop if you want to do it one dataset at a time
%  for i = 1:length(ALLEEG)
%      [ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,'retrieve',i,'study',0);
%      [EEG, ALLEEG] = pop_micro_selectdata( EEG, ALLEEG, 'datatype', 'spontaneous',...
%          'avgref', 1, ...
%          'normalise', 0, ...
%          'MinPeakDist', 10, ...
%          'Npeaks', 3000, ...
%          'GFPthresh', 1, ...
%          'dataset_idx', [] );
%      % store data in a new EEG structure
%      [ALLEEG, EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
%      %  eeglab redraw % updates EEGLAB datasets
%  end

% Use this if you are doing all loaded datasets all together
[EEG, ALLEEG] = pop_micro_selectdata( EEG, ALLEEG, 'datatype', 'spontaneous',...
    'avgref', 1, ...
    'normalise', 0, ...
    'MinPeakDist', 10, ...
    'Npeaks', 1000, ... % This needs to be perhaps confirmed. How many datapoints per participant?
    'GFPthresh', 1, ...
    'dataset_idx', 1:numsubs );



%% 3.4 Microstate segmentation
% select the "GFPpeak" dataset and make it the active set
[ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,'retrieve',numsubs+1,'study',0);
%eeglab redraw
% Perform the microstate segmentation
EEG = pop_micro_segment( EEG, 'algorithm', 'modkmeans', ...
    'sorting', 'Global explained variance', ...
    'Nmicrostates', 3:6, ...
    'verbose', 1, ...
    'normalise', 0, ...
    'Nrepetitions', 50, ... % remember to make this 50 for the actual analysis
    'max_iterations', 2000, ...
    'threshold', 1e-06, ...
    'fitmeas', 'CV',...
    'optimised',1);
[ALLEEG, EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);

%  % select the other "GFPpeak" dataset and make it the active set
%  [ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, 4,'retrieve',6,'study',0);
%  eeglab redraw
%  % Perform the microstate segmentation
%  EEG = pop_micro_segment( EEG, 'algorithm', 'modkmeans', ...
%      'sorting', 'Global explained variance', ...
%      'Nmicrostates', 2:8, ...
%      'verbose', 1, ...
%      'normalise', 0, ...
%      'Nrepetitions', 50, ...
%      'max_iterations', 1000, ...
%      'threshold', 1e-06, ...
%      'fitmeas', 'CV',...
%      'optimised',1);
%  [ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
% [ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,'retrieve',1,'study',0);

%% 3.5 Review and select microstate segmentation
%  %% 3.5.1 Plot microstate prototype topographies
%  figure;MicroPlotTopo( EEG, 'plot_range', [] );
%% 3.5.2 Select active number of microstates
EEG = pop_micro_selectNmicro( EEG);
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);

[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,'retrieve',numsubs+1,'study',0);

EEG = pop_saveset( EEG, 'filename',char(strcat(Sample,"MicroProtoTogether")),'filepath',char(strcat(EEGdir,"Microstates/",Sample,"/Together")));

%% Import microstate prototypes from other dataset to the datasets that should be back-fitted
% note that dataset number 5 is the GFPpeaks dataset with the microstate
% prototypes


for i = 1:numsubs
    fprintf('Importing prototypes and backfitting for dataset %i\n',i)
    [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,'retrieve',i,'study',0);
    EEG = pop_micro_import_proto( EEG, ALLEEG, numsubs+1);

    %% 3.6 Back-fit microstates on EEG
    EEG = pop_micro_fit( EEG, 'polarity', 0 );
    %% 3.7 Temporally smooth microstates labels
    EEG = pop_micro_smooth( EEG, 'label_type', 'backfit', ...
        'smooth_type', 'reject segments', ...
        'minTime', 30, ...
        'polarity', 0 );
    %% 3.9 Calculate microstate statistics
    EEG = pop_micro_stats( EEG, 'label_type', 'backfit', ...
        'polarity', 0 );
    [ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);

    if contains(EEG.filename,"RS1","IgnoreCase",true)
        EEG = pop_saveset( EEG, 'filename',char(strcat("RS1",Sample,"Subject",string(MicrostateData.Subject(i)))),'filepath',char(strcat(EEGdir,"Microstates/",Sample,"/Together")));
    elseif contains(EEG.filename,"RS2","IgnoreCase",true)
        EEG = pop_saveset( EEG, 'filename',char(strcat("RS2",Sample,"Subject",string(MicrostateData.Subject(i)))),'filepath',char(strcat(EEGdir,"Microstates/",Sample,"/Together")));
    end

end
%% 3.8 Illustrating microstate segmentation
% Plotting GFP of active microstates for the first 1500 ms for subject 1.
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,'retrieve',1,'study',0);
figure;MicroPlotSegments( EEG, 'label_type', 'backfit', ...
    'plotsegnos', 'first', 'plot_time', [20 2400], 'plottopos', 1 );
% eeglab redraw

% topoplot(ALLEEG(5).data(:,3000),EEG.chanlocs)

%% Organise the data into tables for each prototype


RS1Microstate1 = MicrostateData;
for i = 1:numsubs
    RS1Microstate1.TotalGEV(i) = ALLEEG(i).microstate.stats.GEVtotal;
    RS1Microstate1.GFP(i) = ALLEEG(i).microstate.stats.Gfp(1);
    RS1Microstate1.Duration(i) = ALLEEG(i).microstate.stats.Duration(1);
    RS1Microstate1.Occurance(i) = ALLEEG(i).microstate.stats.Occurence(1);
    RS1Microstate1.Coverage(i) = ALLEEG(i).microstate.stats.Coverage(1);
    RS1Microstate1.GEV(i) = ALLEEG(i).microstate.stats.GEV(1);
    RS1Microstate1.TP2(i) = ALLEEG(i).microstate.stats.TP(2,1);
    RS1Microstate1.TP3(i) = ALLEEG(i).microstate.stats.TP(3,1);
    RS1Microstate1.TP4(i) = ALLEEG(i).microstate.stats.TP(4,1);
    RS1Microstate1.TP5(i) = ALLEEG(i).microstate.stats.TP(5,1);
    RS1Microstate1.MS(i) = "MS1";

end

RS1Microstate2 = MicrostateData;
for i = 1:numsubs
    RS1Microstate2.TotalGEV(i) = ALLEEG(i).microstate.stats.GEVtotal;
    RS1Microstate2.GFP(i) = ALLEEG(i).microstate.stats.Gfp(2);
    RS1Microstate2.Duration(i) = ALLEEG(i).microstate.stats.Duration(2);
    RS1Microstate2.Occurance(i) = ALLEEG(i).microstate.stats.Occurence(2);
    RS1Microstate2.Coverage(i) = ALLEEG(i).microstate.stats.Coverage(2);
    RS1Microstate2.GEV(i) = ALLEEG(i).microstate.stats.GEV(2);
    RS1Microstate2.TP1(i) = ALLEEG(i).microstate.stats.TP(1,2);
    RS1Microstate2.TP3(i) = ALLEEG(i).microstate.stats.TP(3,2);
    RS1Microstate2.TP4(i) = ALLEEG(i).microstate.stats.TP(4,2);
    RS1Microstate2.TP5(i) = ALLEEG(i).microstate.stats.TP(5,2);
    RS1Microstate2.MS(i)  = "MS2";
end

RS1Microstate3 = MicrostateData;
for i = 1:numsubs
    RS1Microstate3.TotalGEV(i) = ALLEEG(i).microstate.stats.GEVtotal;
    RS1Microstate3.GFP(i) = ALLEEG(i).microstate.stats.Gfp(3);
    RS1Microstate3.Duration(i) = ALLEEG(i).microstate.stats.Duration(3);
    RS1Microstate3.Occurance(i) = ALLEEG(i).microstate.stats.Occurence(3);
    RS1Microstate3.Coverage(i) = ALLEEG(i).microstate.stats.Coverage(3);
    RS1Microstate3.GEV(i) = ALLEEG(i).microstate.stats.GEV(3);
    RS1Microstate3.TP1(i) = ALLEEG(i).microstate.stats.TP(1,3);
    RS1Microstate3.TP2(i) = ALLEEG(i).microstate.stats.TP(2,3);
    RS1Microstate3.TP4(i) = ALLEEG(i).microstate.stats.TP(4,3);
    RS1Microstate3.TP5(i) = ALLEEG(i).microstate.stats.TP(5,3);
    RS1Microstate3.MS(i)  = "MS3";
end

RS1Microstate4 = MicrostateData;
for i = 1:numsubs
    RS1Microstate4.TotalGEV(i) = ALLEEG(i).microstate.stats.GEVtotal;
    RS1Microstate4.GFP(i) = ALLEEG(i).microstate.stats.Gfp(4);
    RS1Microstate4.Duration(i) = ALLEEG(i).microstate.stats.Duration(4);
    RS1Microstate4.Occurance(i) = ALLEEG(i).microstate.stats.Occurence(4);
    RS1Microstate4.Coverage(i) = ALLEEG(i).microstate.stats.Coverage(4);
    RS1Microstate4.GEV(i) = ALLEEG(i).microstate.stats.GEV(4);
    RS1Microstate4.TP1(i) = ALLEEG(i).microstate.stats.TP(1,4);
    RS1Microstate4.TP2(i) = ALLEEG(i).microstate.stats.TP(2,4);
    RS1Microstate4.TP3(i) = ALLEEG(i).microstate.stats.TP(3,4);
    RS1Microstate4.TP5(i) = ALLEEG(i).microstate.stats.TP(5,4);
    RS1Microstate4.MS(i)  = "MS4";
end

RS1Microstate5 = MicrostateData;
for i = 1:numsubs
    RS1Microstate5.TotalGEV(i) = ALLEEG(i).microstate.stats.GEVtotal;
    RS1Microstate5.GFP(i) = ALLEEG(i).microstate.stats.Gfp(5);
    RS1Microstate5.Duration(i) = ALLEEG(i).microstate.stats.Duration(5);
    RS1Microstate5.Occurance(i) = ALLEEG(i).microstate.stats.Occurence(5);
    RS1Microstate5.Coverage(i) = ALLEEG(i).microstate.stats.Coverage(5);
    RS1Microstate5.GEV(i) = ALLEEG(i).microstate.stats.GEV(5);
    RS1Microstate5.TP1(i) = ALLEEG(i).microstate.stats.TP(1,5);
    RS1Microstate5.TP2(i) = ALLEEG(i).microstate.stats.TP(2,5);
    RS1Microstate5.TP3(i) = ALLEEG(i).microstate.stats.TP(3,5);
    RS1Microstate5.TP4(i) = ALLEEG(i).microstate.stats.TP(4,5);
    RS1Microstate5.MS(i)  = "MS5";
end

Microstats = [RS1Microstate1; RS1Microstate2; RS1Microstate3; RS1Microstate4; RS1Microstate5];
writetable(Microstats, strcat(EEGdir,"Microstates/",Sample,"/Together/",'Microstats.xlsx'))



%% Plot the 5 prototypes
clear;clc;

addpath('/Path to EEGLAB/eeglab2023.0')
addpath(genpath('/Matlab Code'))

% start EEGLAB to load all dependent paths
eeglab nogui
Path = '/Path to EEG data/Microstates/Control/Together/';
Name = 'ControlMicroProtoTogether.set';
EEG = pop_loadset('filename',Name,'filepath',Path);

for i = 1:5
figure(i)
topoplot(EEG.microstate.prototypes(:,i),EEG.chanlocs,'electrodes','off'...
    ,'style','map','plotrad',.6)
set(gcf,'color','w');
title(strcat('MS ',num2str(i)))
end


%% Plot an example time course of back fit microstates

Path = '/Path to EEG data/Microstates/Control/Together/';
Name = 'RS1ControlSubject3.set';
EEG = pop_loadset('filename',Name,'filepath',Path);

figure;MicroPlotSegments( EEG, 'label_type', 'backfit', ...
    'plotsegnos', 'first', 'plot_time', [20 2400], 'plottopos', 1 );

