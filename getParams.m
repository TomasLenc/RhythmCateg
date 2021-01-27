function cfg = getParams(task)
% NOTE: in order to use behav + fMRI with 1 getParams, we are using
% getParam with 3 aguments so it won't interfere with behav script
% CB edit 0n 09/07/2020


% Initialize the parameters variables
% Initialize the general configuration variables
% =======
% INPUT: 
% =======
%     task:             string specifying the current task to get parameters for
%                       (tapTraining or tapMainExp)
% =======
% OUTPUT: 
% =======
%     cfg:
%     expParameters:



%% Init parameter structures
cfg = struct(); 

%% set the type of your computer
if IsWin
    cfg.stimComp='windows';
elseif ismac
    cfg.stimComp = 'mac';
elseif IsLinux
    cfg.stimComp = 'linux';
end


%logfile folder - behav exp
cfg.dir.output = fullfile(...
    fileparts(mfilename('fullpath')), ...
    'output');

%% Debug mode settings
cfg.debug.do        = false ;  
cfg.debug.transpWin = false ;     % To test the script with trasparent full size screen
cfg.debug.smallWin  = false;
cfg.verbose         = true;        % add here and there some explanations with if verbose is ON. 

    
%% MRI settings

cfg.testingDevice = 'mri';
cfg.eyeTracker.do = false;          % Set to 'true' if you are testing in MRI and want to record ET data

% PsychPortAudio('GetDevices')
%  cfg.audio.devIdx = 5;
cfg.audio.do = true;

%set visual
cfg = setMonitor(cfg);

% set audio
cfg = setAudio(cfg);

% Keyboards
cfg = setKeyboards(cfg);

% MRI settings
cfg = setMRI(cfg);

cfg.pacedByTriggers.do = false;
    
%% general configuration
%for BIDS format: 
cfg.task.name = task;                % should be calling behav or fmri
cfg.subject.askGrpSess = [1 0]; % it won't ask you about group or session
    
% set and load all the subject input to run the experiment
cfg = userInputs(cfg);
cfg = createFilename(cfg);

%% Timing 

% these are for behavioral exp delays
cfg.sequenceDelay = 1; %wait in between sequences? y/n
cfg.pauseSeq = 1; % give a pause of below seconds in between sequences

% define ideal number of sequences to be made
% multiple of 3 is balanced design
cfg.pattern.numSequences = 6;

if cfg.debug.do
     cfg.pattern.numSequences = 2;
end

if strcmpi(cfg.testingDevice,'mri')
    cfg.pattern.numSequences = 9;
    cfg.pattern.numSeq4Run = 1; % for an fMRI run time calculation
    cfg.pattern.extraSeqNum = 3; % extra session for piloting
end



%% fMRI task
% display a fixation cross during the fMRI run


if strcmpi(cfg.testingDevice,'mri')
    
    %Fixation Cross
    % Used Pixels here since it really small and can be adjusted during the experiment
    cfg.fixation.type                   = 'bestFixation'; %
    cfg.fixation.width                  = .2;   % Set the length of the lines (in Pixels) of the fixation cross
    cfg.fixation.lineWidthPix           = 5;    % Set the line width (in Pixels) for our fixation cross
    cfg.fixation.xDisplacement          = 0;    % Manual displacement of the fixation cross
    cfg.fixation.yDisplacement          = 0;    % Manual displacement of the fixation cross
    cfg.fixation.color                  = cfg.color.white;
    
%     cfg.ctrlscreen.idx = min(Screen('Screens'));
%     
%     if cfg.screen.winWidth < 1920
%         cfg.fixation.width = .2;
%     end

    %Task
    cfg.task.instruction = ['If you saw a shiny ! point, don''t mind it. It''s OK! \n\n\n'...
                           'Your task is: Fixate to the cross & count the piano tones\n\n\n'];
    cfg.task.instructionPress = ['Please indicate by pressing button, '...
                                 'how many times you detected piano tones\n\n\n'];
    cfg.task.instructionEnd = ['You are at the VERY LAST BLOCK.' ...
                               ':)\n\n\n Soon we will take you out!'];
    cfg.task.instructionCont = ['This run is over. We will shortly start'...
                                ' the following ! \n\n\n'];
    % set default for no-task 
    cfg.isTask.Idx = 0;
    
    % how many targets within 1 pattern
    cfg.isTask.numEvent = 1;
    
    %logfile columns
    cfg.extraColumns = {'sequenceNum', 'segmentNum', 'segmentOnset', ...
    'stepNum', 'stepOnset', 'patternID', 'segmentCateg', 'F0', 'isTask', ...
    'gridIOI', 'patternAmp', 'minPE4', 'rangePE4', 'minLHL24', ...
    'rangeLHL24', 'LHL24', 'PE4'};

    % response columns
    cfg.responseExtraColumns = {'keyName', 'pressed', 'target'};

end


%% Load needed stimuli files

%load target tones, missing stimuli, ...
cfg = loadTargetTones(cfg);

%% more parameters to get according to the type of experiment
% this part is solely for behavioral exp
% control fMRI script has its getxxx.m instead of in here (getParam.m)
if strcmp(cfg.task.name,'tapTraining')
    
    % get tapping training parameters
    cfg = getTrainingParameters(cfg);
    
elseif strcmp(cfg.task.name,'tapMainExp') || strcmp(cfg.task.name,'RhythmFT')
    
    % get main experiment parameters
    cfg = getMainExpParameters(cfg);
    
elseif strcmp(cfg.task.name,'RhythmBlock')
    % get main experiment parameters
    cfg = getBlockParameters(cfg);
    
elseif strcmp(cfg.task.name,'PitchFT')
    
    cfg = getPitchParameters(cfg);
    
end




%% differentiating response button (subject) from keyboard(experimenter)
% cfg.responseBox would be the testingDevice used by the participant to give his/her response: 
%   like the button box in the scanner or a separate keyboard for a behavioral experiment
%
% cfg.keyboard is the keyboard on which the experimenter will type or press the keys necessary 
%   to start or abort the experiment.
%   The two can be different or the same.

% Using empty vectors should work for linux when to select the "main"
% keyboard. 


% After connecting fMRI response button to the laptop, 
% LOOK what are the experimenters' keyboard & fmri responseKey 
% by using GetKeyboardIndices below.
% copy-paste these this into command window 
% then you can assign device number to the main keyboard or the response
% box. 
% otherwise it's set for PTB to assign. 
[keyboardNumbers, keyboardNames] = GetKeyboardIndices;

% disp(keyboardNumbers);
% disp(keyboardNames);



end

function cfg = setKeyboards(cfg)

cfg.keyboard.escapeKey = 'ESCAPE';
cfg.keyboard.responseKey = {'d', 'a', 'c', 'b'};
cfg.keyboard.keyboard = [];
cfg.keyboard.responseBox = [];


% %behavioral exp keys to check
% cfg.keywait         = KbName({'RETURN'}); % press enter to start bloc
% cfg.keyquit         = KbName('ESCAPE'); % press ESCAPE at response time to quit
% cfg.keyToggleInstr  = KbName({'I'}); % press I to show/remove general instructions from the screen
% cfg.keytap          = KbName('SPACE');
% cfg.keyVolUp        = KbName('UpArrow');
% cfg.keyVolDown      = KbName('DownArrow');
% cfg.keyAudioPlay    = KbName('p');
% cfg.keyAudioStop    = KbName('s');
% cfg.keyInstrBack    = KbName('b');
% cfg.keyInstrNext    = KbName('n');


if strcmpi(cfg.testingDevice, 'mri')
    cfg.keyboard.keyboard = [];
    cfg.keyboard.responseBox = [];
end
end

function cfg = setMRI(cfg)

% BIDS compatible logfile folder
cfg.dir.output = fullfile(...
    fileparts(mfilename('fullpath')),'..', ...
    'output');
    
% letter sent by the trigger to sync stimulation and volume acquisition
cfg.mri.triggerKey = 's';
cfg.mri.triggerNb = 1; % for hyberpand insert 4 here! ! ! 

% json sidecar file for bold data
cfg.mri.repetitionTime = 1.75;
cfg.bids.MRI.Instructions = 'Fixate to the cross & count the piano tones';
cfg.bids.MRI.TaskDescription = [];
cfg.bids.mri.SliceTiming = [0, 0.9051, 0.0603, 0.9655, 0.1206, 1.0258, 0.181, ...
                      1.0862, 0.2413, 1.1465, 0.3017, 1.2069, 0.362, ...
                      1.2672, 0.4224, 1.3275, 0.4827, 1.3879, 0.5431, ...
                      1.4482, 0.6034, 1.5086, 0.6638, 1.5689, 0.7241, ...
                      1.6293, 0.7844, 1.6896, 0.8448, 0, 0.9051, 0.0603, ...
                      0.9655, 0.1206, 1.0258, 0.181, 1.0862, 0.2413, ...
                      1.1465, 0.3017, 1.2069, 0.362, 1.2672, 0.4224, ...
                      1.3275, 0.4827, 1.3879, 0.5431, 1.4482, 0.6034, ...
                      1.5086, 0.6638, 1.5689, 0.7241, 1.6293, 0.7844, ...
                      1.6896, 0.8448];

%Number of seconds before the rhythmic sequence (exp) are presented
cfg.timing.onsetDelay = 2 *cfg.mri.repetitionTime; %5.2s
% Number of seconds after the end of all stimuli before ending the fmri run!
cfg.timing.endDelay = 4 * cfg.mri.repetitionTime; %10.4s


% ending timings for fMRI
cfg.timing.endScreenDelay = 2; %end the screen after thank you screen
% delay for script ending
cfg.timing.endResponseDelay = 10; % wait for participant to response for counts



end


function cfg = setMonitor(cfg)


% Monitor parameters for PTB
cfg.skipSyncTests = 1;

% Text format 
cfg.text.font         = 'Arial'; %'Courier New'
cfg.text.size         = 48; %18


% Monitor parameters for PTB
cfg.color.white = [255 255 255];
cfg.color.black = [0 0 0];
cfg.color.red = [255 0 0];
cfg.color.grey = mean([cfg.color.black; cfg.color.white]);
cfg.color.background = cfg.color.grey;
cfg.text.color = cfg.color.white;
%cfg.color.foreground =  [127 127 127];

% % Monitor parameters
% if strcmpi(cfg.testingDevice, 'mri')
%     cfg.screen.monitorWidth = 69.8;
%     cfg.screen.monitorDistance = 170;
% end

end


function cfg = setAudio(cfg)

%% audio other parameters
% sampling rate
cfg.audio.fs = 44100; 
%cfg.audio.initVolume = 1;
%cfg.audio.requestedLatency = 2;

%  boolean for equating the dB across different tones for behavioral exp
cfg.equateSoundAmp = 1;
% sound levels
cfg.baseAmp = 0.5;
cfg.PTBInitVolume = 0.3;
   
if strcmpi(cfg.testingDevice, 'mri')  
    
    cfg.baseAmp = 0.99; 
    cfg.PTBInitVolume = 1; 

    cfg.equateSoundAmp = 0;
      
end

end

function cfg = loadTargetTones(cfg)

%download missing stimuli (.wav)
cfg.soundpath = fullfile(fileparts(mfilename('fullpath')));
checkSoundFiles(cfg.soundpath);

% piano keys
% read the audio files and insert them into cfg
targetList = dir(fullfile(cfg.soundpath,'stimuli/Piano*.wav'));

for isound = 1:length(targetList)
    [S,cfg.fs] = audioread(fullfile('stimuli',targetList(isound).name));
    cfg.isTask.targetSounds{isound} = S';
end

end


