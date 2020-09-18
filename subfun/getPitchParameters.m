function     cfg = getPitchParameters(cfg)
% this function generates audio sequences to be played in the man
% experiment


% % %
% start the sequence with one B-category segment that will be discarded during analysis
% % %



%% contruct individual sound events (that will make up each pattern)

% Define envelope shape of the individual sound event. 
% All parameters are defined in seconds. 

% total sound duration _/```\_  
cfg.pattern.eventDur             = 0.190; % s
% onset ramp duration  _/     
cfg.pattern.eventRampon          = 0.010; % s
% offset ramp duration       \_ 
cfg.pattern.eventRampoff         = 0.020; % s

% Make sure the total ramp durations are not longer than tone duration. 
if (cfg.pattern.eventRampon+cfg.pattern.eventRampoff) > cfg.pattern.eventDur
    error(sprintf('The summed duration of onset+offset ramps (%g ms) is longer than requensted tone duration (%g ms).',...
                  (cfg.pattern.eventRampon + cfg.pattern.eventRampoff)*1e3, ...
                  cfg.pattern.eventDur * 1e3)); 
end


%% construct pattern (smallest item in sequence)
cfg.pattern.nGridPoints = 12; % length(pat_complex(1).pattern)

% the grid interval can vary across steps or segments (gridIOI selected 
% randomly from a set of possible values for each new step or segment) 
cfg.pattern.minGridIOI 	= 0.190;  % minimum possible grid IOI 
cfg.pattern.maxGridIOI 	= 0.190; % maximum possible grid IOI 
cfg.pattern.nGridIOI 	= 1; 	% number of unique IOI values between the limits
cfg.pattern.gridIOIs 	= linspace((cfg.pattern.minGridIOI),...
                            (cfg.pattern.maxGridIOI),cfg.pattern.nGridIOI); 

cfg.pattern.interPatternInterval = cfg.pattern.nGridPoints * cfg.pattern.maxGridIOI; 

%================================================================
% The gridIOI changes are controlled by the Boolean variables below. 
% change gridIOI for each segment
cfg.pattern.changeGridIOISegm       = 0;           
% change gridIOI for each segment-type (every time A changes to B or the other way around)
cfg.pattern.changeGridIOICategory   = 0;    
% change gridIOI for each step
cfg.pattern.changeGridIOIStep       = 0;         


% Make sure the tone duration is not longer than smallest gridIOI. 
if cfg.pattern.eventDur >  cfg.pattern.minGridIOI
    error(sprintf('Requested tone duration (%g ms) is longer than shortest gridIOI (%g ms).',...
                  cfg.pattern.eventDur * 1e3, ...
                  cfg.pattern.minGridIOI * 1e3)); 
end

%% construct segment

% % how many times the pattern will be repeated/cycle through
% it set by default = 1 in makeStimMainExp.m
% cfg.nCyclesPerPattern = 1;

% how many pattern cycles are within each step of [ABBB]
% how many pattern in each segment A or B.
cfg.pattern.nPatternPerSegment = 4;

% if the gridIOI can vary across pattern cycles, we need to set the time 
% interval between two successive segments to a fixed value (this must be 
% greater or equal to the maximum possible segment duration)
cfg.pattern.interSegmInterval = cfg.pattern.nPatternPerSegment * ...
                                cfg.pattern.interPatternInterval; 

% there can be a pause after all segments for category A are played 
% (i.e. between A and B)
cfg.pattern.delayAfterA = 0; 
% there can be a pause after all segments for category B are played 
% (i.e. between B and A)
cfg.pattern.delayAfterB = 0; 

%% construct step [ABBB]
% how many successive segments are presented for category A
% manw many times segment A will be sequentially repeated
cfg.pattern.nSegmentA = 1; 

% how many successive segments are presented for category B
% manw many times segment B will be sequentially repeated
cfg.pattern.nSegmentB = 3; 

% number of segments for each step
cfg.pattern.nSegmPerStep = cfg.pattern.nSegmentB + cfg.pattern.nSegmentA; %4; 

%calculation duration (s) of a step 
cfg.pattern.interStepInterval = (cfg.pattern.interSegmInterval * ...
                                cfg.pattern.nSegmPerStep) + ...
                                (cfg.pattern.delayAfterA * ...
                                cfg.pattern.nSegmentA) + ... 
                                (cfg.pattern.delayAfterB * ...
                                cfg.pattern.nSegmentB);

%% construct whole sequence
% how many steps are in the whole sequence
% how many repetition of grouped segment [ABBB] in the whole sequence
% if debug for behav exp, cut it short! 
if cfg.debug.do && strcmpi(cfg.testingDevice,'pc')
    cfg.pattern.nStepsPerSequence = 1;
else
    cfg.pattern.nStepsPerSequence = 5;
end


% calculate trial duration (min)
cfg.pattern.SequenceDur = (cfg.pattern.interStepInterval * ...
                           cfg.pattern.nStepsPerSequence); 
fprintf('\n\nsequence duration is: %.1f minutes\n',cfg.pattern.SequenceDur/60);

%================================================================
% The pitch changes are controlled by the Boolean variables below. 
% NOTE: the parameters work together hierarchically, i.e. if you set
% cfg.changePitchCycle = TRUE, then it's obvious that pitch will be changed
% also every segment and step...

% change pitch for each new pattern cycle
cfg.pattern.changePitchCycle 	= 1;
% change pitch for each segment
cfg.pattern.changePitchSegm 	= 0;           
% change pitch for each segment-category (every time A changes to B or the other way around)
cfg.pattern.changePitchCategory = 0;    
% change pitch for each step
cfg.pattern.changePitchStep 	= 0;  

% refuse to pitch-change in categB
cfg.fixedPitchCategB   = 1;

%% construct pitch features of the stimulus 
% the pitch (F0) of the tones making up the patterns can vary 
% (it can be selected randomly from a set of possible values)

%for this design, we want categ A and B having sifferent pitches (A will
%have 5 pitches differing, B will only have 1 pitch)
cfg.pattern.minF0 	= 349.228; % 350 or 349.228 minimum possible F0
cfg.pattern.maxF0 	= 880; % 900 or 880 maximum possible F0
cfg.pattern.nF0 	= 5; % number of unique F0-values between the limits
cfg.pattern.F0s 	= logspace(log10(cfg.pattern.minF0),...
                      log10(cfg.pattern.maxF0),cfg.pattern.nF0); 

cfg.differF0 = 277.183; % this is also logspaced
% calculate required amplitude gain
% butchered this part with adding cfg.differ - change in the future ! ! !
if cfg.equateSoundAmp
    cfg.pattern.F0sAmpGain =  equalizePureTones([cfg.F0s,cfg.pattern.differF0],[], []);
else
    cfg.pattern.F0sAmpGain = ones(1,cfg.pattern.nF0);
end


% use the requested gain of each tone to adjust the base amplitude
cfg.pattern.F0sAmps = cfg.baseAmp * cfg.pattern.F0sAmpGain; 

%% create two sets of patterns
cfg = readPatternText(cfg);
%% generate sequence

% get pattern IDs for all sequences used in the experiment
% this is to make sure each pattern is used equal number of times in the
% whole experiment
cfg.pattern.labelCategA = 'simple'; 
cfg.pattern.labelCategB = 'complex';
% for exp like fmri that we will present 1 sequence per run, we are
% creating full exp design in the first run and saving it for the other
% runs to call .mat file

% create randomized sequence for 9 runs when run =1
% overwrites cfg.pattern.seqDesignFullExp
cfg = makefMRISeqDesign(cfg);

% overwrite the base amp
cfg = normaliseEvent(cfg);
cfg.pattern.F0sAmps = cfg.baseAmp * cfg.pattern.F0sAmpGain * ...
    cfg.isTask.rmsRatio;

% provide an error if it amp above 1 !
% % %

% % %
% can I normalise the target sound to make the max [-1 1]
% to increase the amplitude of whole sounds?
% % %



%% Task Instructions
% fMRI instructions
cfg.fmriTaskInst = ['Fixate to the cross & count the deviant tone\n \n\n'];

% ------------------------------------------------
% instruction showing info about sequence curation 
% ------------------------------------------------

cfg.trialDurInstruction = [sprintf('Trial duration will be: %.1f minutes\n\n',cfg.SequenceDur/60), ...
                            'Set your volume now. \n\n\nThen start the experiment whenever ready...\n\n']; 
                           
% ------------------------------
% sequence-specific instructions
% ------------------------------

% this is general instruction displayed after each sequence
cfg.generalDelayInstruction = ['The %d out of %d is over!\n\n', ...
                            'You can have a break. \n\n',...
                            'Good luck!\n\n']; 


    
end

function cfg = readPatternText(cfg)

% read from txt files
grahnPatSimple = loadIOIRatiosFromTxt(fullfile('stimuli','Grahn2007_simple.txt')); 
grahnPatComplex = loadIOIRatiosFromTxt(fullfile('stimuli','Grahn2007_complex.txt')); 

% get different metrics of the patterns
cfg.pattern.patternSimple = getPatternInfo(grahnPatSimple, 'simple',cfg); 
cfg.pattern.patternComplex = getPatternInfo(grahnPatComplex, 'complex', cfg); 

end

function cfg = normaliseEvent(cfg)

%mini subfunction to normalise sound event according to the target rms

% make the env and sound for 1 event
[s, EventEnv] = makeEvent(cfg);
s = s .*cfg.baseAmp;

% calculate the rms of an event
cfg.isTask.rmsEvent = rms(s);


% rms the target
% find the biggest rms among the target sounds
for i = length(cfg.isTask.targetSounds)
    % apply env
    currTargetS = cfg.isTask.targetSounds{i}.*EventEnv.*cfg.baseAmp;%
    % take rms of all target sounds
    rmsAllTarget(i) = rms(currTargetS); %#ok<AGROW>
end
% use the smallest target rms as reference
cfg.isTask.rmsTarget = max(rmsAllTarget);

%overwrite the baseAmp with normalised value
cfg.isTask.rmsRatio = cfg.isTask.rmsTarget/cfg.isTask.rmsEvent;


end