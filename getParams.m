function cfg = getParams(task,cfg)
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

  %% set the type of your computer
  if IsWin
    cfg.stimComp = 'windows';
  elseif ismac
    cfg.stimComp = 'mac';
  elseif IsLinux
    cfg.stimComp = 'linux';
  end

  % logfile folder - behav exp
  cfg.dir.output = fullfile( ...
                            fileparts(mfilename('fullpath')), '..', ...
                            'output');

  %% Debug mode settings
  cfg.debug.do        = false;
  cfg.debug.transpWin = true;     % To test the script with trasparent full size screen
  cfg.debug.smallWin  = false;
  cfg.verbose         = false;        % add here and there some explanations with if verbose is ON.

  %% MRI settings

  cfg.testingDevice = 'pc';
  cfg.eyeTracker.do = false;          % Set to 'true' if you are testing in MRI and want to record ET data

  cfg.audio.do = true;

  % set visual
  cfg = setMonitor(cfg);

  % set audio
  cfg = setAudio(cfg);
  
  % set audio with device
 if strcmpi(cfg.testingDevice, 'pc') && ~strcmpi(task, 'tapTraining')
     cfg = setAudioExtend(cfg);
 end

  % Keyboards
  cfg = setKeyboards(cfg);

  cfg.pacedByTriggers.do = false;

  % task - set default for no-task
  cfg.isTask.Idx = 0;
  %% general configuration
  % for BIDS format:
  cfg.task.name = task;                % should be calling behav or fmri
  cfg.subject.askGrpSess = [0 0]; % it won't ask you about group or session
  
  % ask user input if it is tapTraining exp/session
  if strcmpi(task, 'tapTraining') 
    cfg = userInputs(cfg);
  end
      
  % set runNb  to 1 for debugging 
  if cfg.debug.do == 1
    cfg.subject.runNb = 1; 
  end
  cfg = createFilename(cfg);

  %% Timing

  % these are for behavioral exp delays
  cfg.timing.startDelay = 2; % wait before starting exp
  cfg.timing.breakDelay = 2; % give a pause of below seconds in between sequences
  cfg.timing.stopDelay = 1; % wait before ending the experiment

  % define ideal number of sequences to be made
  % multiple of 3 is balanced design
  cfg.pattern.numSequences = 6;

  if cfg.debug.do
    cfg.pattern.numSequences = 2;
  end

  if strcmpi(cfg.testingDevice, 'mri')
    cfg.pattern.numSequences = 9;
    cfg.pattern.numSeq4Run = 1; % for an fMRI run time calculation
    cfg.pattern.extraSeqNum = 3; % extra session for piloting
  end

  %% Logfile

  % logfile columns
  cfg.extraColumns = {'sequenceNum', 'segmentNum', 'segmentOnset', ...
                      'stepNum', 'stepOnset', 'patternID', 'segmentCateg', 'F0', 'isTask', ...
                      'gridIOI', 'patternAmp', 'minPE4', 'rangePE4', 'minLHL24', ...
                      'rangeLHL24', 'LHL24', 'PE4', 'triggerValue'};

  %% fMRI task
  % display a fixation cross during the fMRI run

  if strcmpi(cfg.testingDevice, 'mri')

    % Fixation Cross
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

    % Task
    cfg.task.instruction = ['If you saw a shiny ! point, don''t mind it. It''s OK! \n\n\n'...
                            'Your task is: Fixate to the cross & count the piano tones\n\n\n'];
    cfg.task.instructionPress = ['Please indicate by pressing button, '...
                                 'how many times you detected piano tones\n\n\n'];
    cfg.task.instructionEnd = ['You are at the VERY LAST BLOCK.' ...
                               ':)\n\n\n Soon we will take you out!'];
    cfg.task.instructionCont = ['This run is over. We will shortly start'...
                                ' the following ! \n\n\n'];

    % how many targets within 1 pattern
    cfg.isTask.numEvent = 1;

    % response columns
    cfg.responseExtraColumns = {'keyName', 'pressed', 'target'};

  end

  %% behavioral 
  
    % mapping of trigger values onto tasks
    cfg.beh.trigTaskMapping = containers.Map({'RhythmFT','RhythmBlock'}, ...
                                             {1, 2});
  
  %% Load needed stimuli files

  % load target tones, missing stimuli, ...
  cfg = loadTargetTones(cfg);

  %% more parameters to get according to the type of experiment
  
  % control fMRI script has its getxxx.m instead of in here (getParam.m)
  % get main experiment parameters
  if strcmp(cfg.task.name, 'RhythmFT')
    
      cfg = getMainExpParameters(cfg);
    
    % load instructions
    [cfg] = makeBehavInstruction(cfg);
    
  % if main exp is pitchFT or Block:
  elseif strcmp(cfg.task.name, 'RhythmBlock')
    
      cfg = getBlockParameters(cfg);
    
    % load instructions
    [cfg] = makeBehavInstruction(cfg);
    
  elseif strcmp(cfg.task.name, 'PitchFT')
    cfg = getPitchParameters(cfg);
    
    % load instructions
    [cfg] = makeBehavInstruction(cfg);
    
  elseif strcmp(cfg.task.name, 'tapTraining')
      
    % get tapping training parameters
    cfg = getTrainingParameters(cfg);
      
  end


  %% behavioral instructions
 

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

  KbName('UnifyKeyNames') % copmatibility across OS

  cfg.keyboard.escapeKey = 'ESCAPE';
  cfg.keyboard.responseKey = {'d', 'a', 'c', 'b'};
  cfg.keyboard.keyboard = [];
  cfg.keyboard.responseBox = [];

  % behavioral exp keys to check
  cfg.keyboard.wait         = KbName({'RETURN'}); % press enter to start bloc
  cfg.keyboard.quit         = KbName('ESCAPE'); % press ESCAPE at response time to quit
  cfg.keyboard.toggleInstr  = KbName({'I'}); % press I to show/remove general instructions from the screen
  cfg.keyboard.tap          = KbName('SPACE');
  cfg.keyboard.volUp        = KbName('UpArrow');
  cfg.keyboard.volDown      = KbName('DownArrow');
  cfg.keyboard.audioPlay    = KbName('p');
  cfg.keyboard.audioStop    = KbName('s');
  cfg.keyboard.instrBack    = KbName('b');
  cfg.keyboard.instrNext    = KbName('n');

  if strcmpi(cfg.testingDevice, 'mri')
    cfg.keyboard.keyboard = [];
    cfg.keyboard.responseBox = [];
  end
end

function cfg = setMRI(cfg)

  % BIDS compatible logfile folder
  cfg.dir.output = fullfile( ...
                            fileparts(mfilename('fullpath')), '..', ...
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

  % Number of seconds before the rhythmic sequence (exp) are presented
  cfg.timing.onsetDelay = 2 * cfg.mri.repetitionTime; % 5.2s
  % Number of seconds after the end of all stimuli before ending the fmri run!
  cfg.timing.endDelay = 4 * cfg.mri.repetitionTime; % 10.4s

  % ending timings for fMRI
  cfg.timing.endScreenDelay = 2; % end the screen after thank you screen
  % delay for script ending
  cfg.timing.endResponseDelay = 10; % wait for participant to response for counts

end

function cfg = setMonitor(cfg)

  % Monitor parameters for PTB
  cfg.skipSyncTests = 1;

  % Text format
  cfg.text.font = 'Arial';
  cfg.text.size = 30;

  % Monitor parameters for PTB
  cfg.color.white = [255 255 255];
  cfg.color.black = [0 0 0];
  cfg.color.red = [255 0 0];
  cfg.color.grey = mean([cfg.color.black; cfg.color.white]);
  cfg.color.background = cfg.color.grey;
  cfg.text.color = cfg.color.white;
  % cfg.color.foreground =  [127 127 127];

  % Monitor parameters
  if strcmpi(cfg.testingDevice, 'mri')
    % save the distance - not really important for our exp
    cfg.screen.monitorWidth = 69.8;
    cfg.screen.monitorDistance = 170;

    % adapt text size for MRI monitor
    cfg.text.size = 48;
  end

end

function cfg = setAudio(cfg)

  %% audio other parameters
  % sampling rate
  cfg.audio.fs = 44100;
  % cfg.audio.initVolume = 1;
  % cfg.audio.requestedLatency = 2;
  cfg.audio.channels = [2];
  
  cfg.audio.requestTimeOffset = 1;
  cfg.audio.requestSampleOffset = cfg.audio.fs * cfg.audio.requestTimeOffset;
  
  cfg.audio.pushTime = 0.01;
  cfg.audio.pushSample = cfg.audio.fs * cfg.audio.pushTime;
  

  % boolean for equating the dB across different tones for behavioral exp
  cfg.equateSoundAmp = 1;
  
  % sound levels
  cfg.audio.baseAmp = 0.5;
  cfg.audio.initVolume = 0.1; % CAREFUL, safety first with in-ears

  if strcmpi(cfg.testingDevice, 'mri')

    cfg.audio.baseAmp = 0.99;
    cfg.audio.initVolume = 1;

    cfg.equateSoundAmp = 0;

  end

end

function cfg = setAudioExtend(cfg)
  
  cfg.audio.useDevice = true; 
  cfg.audio.deviceName = 'Fireface'; 
  
  % open only 8 channels in 8 out (max 18 on RME)
  cfg.audio.channels = [8,8];
  
  % 1: playback, 2: capture, 3: simult playback+capture
  cfg.audio.playbackMode = 3;

  % 3: most drastic setting
  cfg.audio.requestedLatency = 3;

  % downsampling frequency (to log tap force data)
  cfg.audio.fsDs = 200;
  
  % mapping of trigger values onto audio output channels
  cfg.audio.trigChanMapping = containers.Map({1, 2, 3}, ...
                                             {[3], [4], [3,4]});
  
  % each small buffer push small duration only (e.g. 0.100 s)
  cfg.audio.pushDur  = 0.200;
  
  % first push will be longer (e.g. 5 s)
  cfg.audio.initPushDur = 5;
  
  % if we're doing capture, we need to initialize buffer with enough space to
  % store enough data between pushes (when we also pull everything and
  % reset the acquisition buffer) 
  cfg.audio.tapBuffDur = 30; % 30s is fine


end

function cfg = loadTargetTones(cfg)

  % download missing stimuli (.wav)
  cfg.soundpath = fullfile(fileparts(mfilename('fullpath')));
  checkSoundFiles(cfg.soundpath);

  % piano keys
  % read the audio files and insert them into cfg
  targetList = dir(fullfile(cfg.soundpath, 'stimuli/Piano*.wav'));

  for isound = 1:length(targetList)
    [S, cfg.fs] = audioread(fullfile('stimuli', targetList(isound).name));
    cfg.isTask.targetSounds{isound} = S';
  end

end
