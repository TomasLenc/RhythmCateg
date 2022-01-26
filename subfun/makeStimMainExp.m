function [s,env] = makeStimMainExp(pattern, patternIOIs, cfg, currGridIOI, currF0, taskF0idx, isTask, varargin)
% this function creates pattern cycles according to the grid that was
% provided
% if nCycles = 1, it will create only 1 time repeated pattern

% ------
% INPUT
% ------
%   pattern:        vector of grid representation of the rhtthmic
%                   pattern to create (e.g. [1,1,1,0,1,1,1,0,0,1,...])
%   cfg:
%
% ------
% OUTPUT
% ------
%   s:             audio waveform of the output
%   env:           envelope of the output

% I added as varargin in case you are using this function somewhere else
% than the tapMainExperiment
if nargin<8
    currAmp = 1;
else
    currAmp = varargin{1};
end

if strcmp(cfg.testingDevice,'mri')
    %use rms ratio for target
    currAmp = cfg.isTask.rmsRatio;
end


% get the time of each tone event relative to the start of the pattern 
soundOnsetTimes = cumsum([0, patternIOIs]); 

%% make envelope for the individual sound event

% number of samples for the whole tone duration 
eventSamples   = round(cfg.pattern.eventDur * cfg.fs);

% number of samples for the onset ramp (proportion of gridIOI)
ramponSamples   = round(cfg.pattern.eventRampon * cfg.fs);

% number of samples for the offset ramp (proportion of gridIOI)
rampoffSamples  = round(cfg.pattern.eventRampoff * cfg.fs);

% individual sound event duration defined as proportion of gridIOI
envEvent = ones(1, round(cfg.pattern.gridIOIs * cfg.fs));

% make the linear ramps
envEvent(1:ramponSamples) = envEvent(1:ramponSamples) .* linspace(0,1,ramponSamples);

envEvent(eventSamples-rampoffSamples+1:eventSamples) = ... 
    envEvent(eventSamples-rampoffSamples+1:eventSamples) .* linspace(1,0,rampoffSamples);

envEvent(eventSamples+1:end) = 0; 

%% make time vector for individual sound event 

tEvent = [0:length(envEvent)-1]/cfg.fs; 

%% calculate the rms for normalisation

% generate task sound
if isTask
    
    % indices of sound events in the pattern
    idxTargets = find(pattern);
    nSoundsEvents = length(idxTargets); 

    % find first N non-zero element
    % number of targets (numEvent) in a pattern is defined in getParam.m
    if cfg.isTask.numEvent < nSoundsEvents
        % the second tone (and onwards) in the pattern will be targets
        idxTargets = idxTargets([2:1+cfg.isTask.numEvent]); 
    end
    
    % check the current F0 & find the corresponding cfg.targetSound
    if length(taskF0idx) == 1
        targetSoundIdx = taskF0idx;  
    else
        targetSoundIdx = taskF0idx(idxTargets);
    end
    
    currTargetS = cfg.isTask.targetSounds{targetSoundIdx}; 

end

%% synthesize whole pattern

% construct time vector
t = [0 : (soundOnsetTimes(end)+currGridIOI) * cfg.fs] / cfg.fs; 

% construct envelope for the whole pattern
env = zeros(1,length(t));
s = zeros(1,length(t));

for iSound=1:length(soundOnsetTimes)

    %get the idx for inserting events
    idx = round(soundOnsetTimes(iSound)*cfg.fs);

    %insert sound event envelope
    env(idx+1:idx+length(envEvent)) = envEvent;
    
    % insert piano key when there's target
    if isTask && ismember(iSound,targetSoundIdx) % target - piano key
        % apply envelop to the target
        soundEvent = currTargetS.*envEvent;

    % no target = sine wave
    else 
        if length(currF0)>1
            F0 = currF0(iSound); 
        else
            F0 = currF0; 
        end
        soundEvent = sin(2*pi*F0*tEvent);
        soundEvent = soundEvent.* envEvent;
        % apply the amplitude
        soundEvent = soundEvent.*currAmp;
    end
    
    %insert sound event
    s(idx+1:idx+length(soundEvent)) = soundEvent;
    
end

    
% to visualise 1 pattern
%  figure; plot(t,s);
% ylim([-1.5,1.5])



end






