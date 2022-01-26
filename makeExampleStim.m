% This script just goes through one run without actually playing the audio.
% The purpose is to same sound data for debuging. It saves the .mat file as usual,
% and also a .wav file. 

clc;
clear;

% make sure we got access to all the required functions and inputs
initEnv();

% Define the task = 'RhythmFT', 'RhythmBlock', 'Nonmetric'
% Get parameters by providing task name
cfg = getParams('NonmetricFT');

%% Init Experiment
% get time point at the beginning of the script (machine time)
cfg.timing.scriptStartTime = GetSecs();

% Init the experiment
[cfg] = initPTB(cfg);

% wait for trigger from fMRI
cfg.experimentStart = 0000;

%% Start Experiment

% take the runNb corresponding sequence
iSequence = cfg.subject.runNb;

% prep for BIDS saving structures
currSeq = struct();
responseEvents = struct();

% construct sequence
currSeq = makeSequenceGrahn2007(cfg, iSequence);

onset = 0000; 

%% save timing and sequence info
% ===========================================
% log into matlab structure
% ===========================================
cfg.timing.sequenceNb = iSequence;
cfg.timing.sequenceStart = onset;
cfg.timing.experimentStart = cfg.experimentStart;
cfg.data(iSequence).sequenceStart = onset;
cfg.data(iSequence).ptbVolume = PsychPortAudio('Volume', cfg.audio.pahandle);
cfg.data(iSequence).seq = currSeq;

matFile = fullfile(cfg.dir.output, ...
                    strrep(cfg.fileName.events, 'tsv', 'mat'));
                
save(matFile, '-v7.3');


wavFile = fullfile(cfg.dir.output, ...
                    strrep(cfg.fileName.events, 'tsv', 'wav'));
                
audiowrite(wavFile, currSeq(1).outAudio, cfg.audio.fs); 


cleanUp;

