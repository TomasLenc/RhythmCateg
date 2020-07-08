function     [cfg,expParam] = getMainExpParameters(cfg,expParam)
% this function generates audio sequences to be played in the man
% experiment


% % %
% start the sequence with one B-category segment that will be discarded during analysis
% % %



%% contruct individual sound events (that will make up each pattern)

% Define envelope shape of the individual sound event. 
% All parameters are defined in seconds. 

% total sound duration _/```\_  
cfg.soundDur             = 0.190; % s
% onset ramp duration  _/     
cfg.eventRampon          = 0.010; % s
% offset ramp duration       \_ 
cfg.eventRampoff         = 0.020; % s

% Make sure the total ramp durations are not longer than tone duration. 
if (cfg.eventRampon+cfg.eventRampoff) > cfg.soundDur
    error(sprintf('The summed duration of onset+offset ramps (%g ms) is longer than requensted tone duration (%g ms).',...
                  (cfg.eventRampon+cfg.eventRampoff)*1e3, ...
                  cfg.soundDur*1e3)); 
end

%% construct pattern (smallest item in sequence)
cfg.nGridPoints = 12; % length(pat_complex(1).pattern)

% the grid interval can vary across steps or segments (gridIOI selected 
% randomly from a set of possible values for each new step or segment) 
cfg.minGridIOI 	= 0.190;  % minimum possible grid IOI 
cfg.maxGridIOI 	= 0.190; % maximum possible grid IOI 
cfg.nGridIOI 	= 1; 	% number of unique IOI values between the limits
cfg.gridIOIs 	= linspace((cfg.minGridIOI),(cfg.maxGridIOI),cfg.nGridIOI); 

cfg.interPatternInterval = cfg.nGridPoints * cfg.maxGridIOI; 

%================================================================
% The gridIOI changes are controlled by the Boolean variables below. 
% change gridIOI for each segment
cfg.changeGridIOISegm       = 0;           
% change gridIOI for each segment-type (every time A changes to B or the other way around)
cfg.changeGridIOICategory   = 0;    
% change gridIOI for each step
cfg.changeGridIOIStep       = 0;     


% Make sure the tone duration is not longer than smallest gridIOI. 
if cfg.soundDur >  cfg.minGridIOI
    error(sprintf('Requested tone duration (%g ms) is longer than shortest gridIOI (%g ms).',...
                  cfg.soundDur*1e3, ...
                  cfg.minGridIOI*1e3)); 
end

%% construct segment

% % how many times the pattern will be repeated/cycle through
% it set by default = 1 in makeStimMainExp.m
% cfg.nCyclesPerPattern = 1;

% how many pattern cycles are within each step of [ABBB]
% how many pattern in each segment A or B.
cfg.nPatternPerSegment = 4;


% there can be a pause after all segments for category A are played 
% (i.e. between A and B)
cfg.delayAfterA = 0; 
% there can be a pause after all segments for category B are played 
% (i.e. between B and A)
cfg.delayAfterB = 0; 

% if the gridIOI can vary across pattern cycles, we need to set the time 
% interval between two successive segments to a fixed value (this must be 
% greater or equal to the maximum possible segment duration)
cfg.interSegmInterval = cfg.nPatternPerSegment * cfg.interPatternInterval; 

%durSeg = cfg.interSegmInterval;

%% construct step [ABBB]
% how many successive segments are presented for category A
% manw many times segment A will be sequentially repeated
cfg.nSegmentA = 1; 

% how many successive segments are presented for category B
% manw many times segment B will be sequentially repeated
cfg.nSegmentB = 3; 

% number of segments for each step
cfg.nSegmPerStep = cfg.nSegmentB + cfg.nSegmentA; %4; 

%calculation duration (s) of a step 
cfg.interStepInterval = (cfg.interSegmInterval * cfg.nSegmPerStep) + ...
                        (cfg.delayAfterA * cfg.nSegmentA) + ... 
                        (cfg.delayAfterB * cfg.nSegmentB);

%% construct whole sequence
% how many steps are in the whole sequence
% how many repetition of grouped segment [ABBB] in the whole sequence
if cfg.debug
    cfg.nStepsPerSequence = 1;
else
    cfg.nStepsPerSequence = 5;
end


% calculate trial duration (min)
cfg.SequenceDur = (cfg.interStepInterval * cfg.nStepsPerSequence); 
fprintf('\n\nsequence duration is: %.1f minutes\n',cfg.SequenceDur/60);



%% construct pitch features of the stimulus 
% the pitch (F0) of the tones making up the patterns can vary 
% (it can be selected randomly from a set of possible values)
cfg.minF0 	= 350; % minimum possible F0
cfg.maxF0 	= 900; % maximum possible F0
cfg.nF0 	= 5; % number of unique F0-values between the limits
cfg.F0s 	= logspace(log10(cfg.minF0),log10(cfg.maxF0),cfg.nF0); 

% calculate required amplitude gain
if expParam.equateSoundAmp
    cfg.F0sAmpGain = equalizePureTones(cfg.F0s,[], []);
else
    cfg.F0sAmpGain = ones(1,cfg.nF0);
end

% use the requested gain of each tone to adjust the base amplitude
cfg.F0sAmps = cfg.baseAmp * cfg.F0sAmpGain; 


%================================================================
% The pitch changes are controlled by the Boolean variables below. 
% NOTE: the parameters work together hierarchically, i.e. if you set
% cfg.changePitchCycle = TRUE, then it's obvious that pitch will be changed
% also every segment and step...

% change pitch for each new pattern cycle
cfg.changePitchCycle 	= 1;
% change pitch for each segment
cfg.changePitchSegm 	= 0;           
% change pitch for each segment-category (every time A changes to B or the other way around)
cfg.changePitchCategory = 0;    
% change pitch for each step
cfg.changePitchStep 	= 0;     


%% create two sets of patterns

% read from txt files
grahnPatSimple = loadIOIRatiosFromTxt(fullfile('stimuli','Grahn2007_simple.txt')); 
grahnPatComplex = loadIOIRatiosFromTxt(fullfile('stimuli','Grahn2007_complex.txt')); 

% get different metrics of the patterns
cfg.patternSimple = getPatternInfo(grahnPatSimple, 'simple',cfg); 
cfg.patternComplex = getPatternInfo(grahnPatComplex, 'complex', cfg); 


%% generate sequence

% get pattern IDs for all sequences used in the experiment
% this is to make sure each pattern is used equal number of times in the
% whole experiment

cfg.labelCategA = 'simple'; 
cfg.labelCategB = 'complex'; 
%%%%%%%%%%%%
% ! important, the order of arguments matters ! -> getAllSeq(categA, categB, ...)
%%%%%%%%%%%%
cfg.seqDesignFullExp = getAllSeqDesign(cfg.patternSimple, cfg.patternComplex, cfg, expParam); 


%% generate example audio for volume setting
% added F0s-amplitude because the relative dB set in volume adjustment in
% PychPortAudio will be used in the mainExp
cfg.volumeSettingSound = repmat(makeStimMainExp(ones(1,16), cfg,...
    cfg.gridIOIs(end), cfg.F0s(end), cfg.F0sAmps(end) ), 2,1); 




%% Task Instructions

% fMRI instructions
expParam.fmriTaskInst = ['Fixate to the cross & count the deviant tone\n \n\n'];




% behavioral instructions

loadPathInstr = fullfile('lib','instr','mainExp'); 

% -------------------
% intro instructions
% -------------------
% These need to be saved in separate files, named: 'instrMainExpIntro#'
% The text in each file will be succesively (based on #) displayed on 
% the screen at the begining of the experiment. Every time, the script 
% will wait for a keypress. 

dirInstr = dir(fullfile(loadPathInstr,'instrMainExpIntro*')); 
expParam.introInstruction = cell(1,length(dirInstr)); 
for i=1:length(dirInstr)
    instrFid = fopen(fullfile(loadPathInstr, dirInstr(i).name),'r','n','UTF-8'); 
    while ~feof(instrFid)
        expParam.introInstruction{i} = [expParam.introInstruction{i}, fgets(instrFid)]; 
    end
    fclose(instrFid); 
end

% ------------------------
% general task instructions
% ------------------------
% This is a general summary of the instructions. Participants can toggle
% these on the screen between sequences if they forget, or want to make
% sure they understand their task. 

dirInstr = dir(fullfile(loadPathInstr,'instrMainExpGeneral')); 
expParam.generalInstruction = ''; 
instrFid = fopen(fullfile(loadPathInstr, dirInstr.name),'r','n','UTF-8'); 
while ~feof(instrFid)
    expParam.generalInstruction = [expParam.generalInstruction, fgets(instrFid)]; 
end
fclose(instrFid); 



% ------------------------------------------------
% instruction showing info about sequence curation 
% ------------------------------------------------

expParam.trialDurInstruction = [sprintf('Trial duration will be: %.1f minutes\n\n',cfg.SequenceDur/60), ...
                            'Set your volume now. \n\n\nThen start the experiment whenever ready...\n\n']; 
    

                        
% ------------------------------
% sequence-specific instructions
% ------------------------------

% this is general instruction displayed after each sequence
expParam.generalDelayInstruction = ['The %d out of %d is over!\n\n', ...
                            'You can have a break. \n\n',...
                            'Good luck!\n\n']; 


% For each sequence, there can be additional instructions. 
% Save as text file with name: 'instrMainExpDelay#', 
% where # is the index of the sequence after which the
% instruction should appear. 

dirInstr = dir(fullfile(loadPathInstr,'instrMainExpDelay*')); 
expParam.seqSpecificDelayInstruction = cell(1, expParam.numSequences); 
for i=1:length(dirInstr)
    
    targetSeqi = regexp(dirInstr(i).name, '(?<=instrMainExpDelay)\d*', 'match'); 
    targetSeqi = str2num(targetSeqi{1}); 
    instrFid = fopen(fullfile(loadPathInstr, dirInstr(i).name),'r','n','UTF-8'); 
    while ~feof(instrFid)
        expParam.seqSpecificDelayInstruction{targetSeqi} = [expParam.seqSpecificDelayInstruction{i}, fgets(instrFid)]; 
    end
    fclose(instrFid); 
end



    
end


