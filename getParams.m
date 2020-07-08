function [cfg,expParam] = getParams(task)
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
cfg = struct; 
expParam = struct;

%% set the type of your computer
if IsWin
    cfg.stimComp='windows';
elseif ismac
    cfg.stimComp = 'mac';
elseif IsLinux
    cfg.stimComp = 'linux';
end

%% Debug mode settings
cfg.debug               = 1 ;  % To test the script with trasparent full size screen 
expParam.verbose        = true; % add here and there some explanations with if verbose is ON. 

%% MRI settings
cfg.device        = 'scanner';       % 'PC': does not care about trigger(for behav) - otherwise use 'Scanner'
cfg.triggerKey    = 's';        % Set the letter sent by the trigger to sync stimulation and volume acquisition
cfg.numTriggers   = 4;          % first #Triggers will be dummy scans
cfg.eyeTracker    = false;      % Set to 'true' if you are testing in MRI and want to record ET data

%% general configuration
expParam.fmriTask = true; % the task is behav exp or fMRI? 

% it'll only look for space press -
% later on change with the responseBox indices/numbers! ! !
expParam.responseKey = {'space'};

%it should be calling behav or fmri - important for BIDS format.
expParam.task = task;

%it won't ask you about group or session
expParam.askGrpSess = [0 0];

%esc key for both behav and fmri exp
cfg.keyquit         = KbName('ESCAPE'); % press ESCAPE at response time to quit


%% monitor
% Monitor parameters - fMRI - CHANGE with fMRI parameters
cfg.monitorWidth  	  = 42;  % Monitor Width in cm
cfg.screenDistance    = 134; % Distance from the screen in cm
cfg.diameterAperture  = 8;   % Diameter/length of side of aperture in Visual angles

% Monitor parameters for PTB
cfg.white  = [255 255 255];
cfg.black  = [ 0   0   0 ];
cfg.gray   = mean([cfg.black; cfg.white]);
cfg.backgroundColor  = cfg.gray;
cfg.textColor        = cfg.white;
cfg.textFont         = 'Arial'; %'Courier New'
cfg.textSize         = 30; %18
%cfg.textStyle        = 1;


    
%% sound levels
% assuming that participant will do the task with headphones
cfg.baseAmp = 0.5; 

% i think this cannot be smaller than cfg.Amp! ! !
cfg.PTBInitVolume = 0.3; 


if strcmpi(cfg.device, 'scanner')
    
    %  boolean for equating the dB across different tones for behavioral exp
    expParam.equateSoundAmp = 0;
    
    % BIDS compatible logfile folder
    expParam.outputDir = fullfile(...
        fileparts(mfilename('fullpath')),'..', ...
        'output');
else
    
    %  boolean for equating the dB across different tones for behavioral exp
    expParam.equateSoundAmp = 1;
    
    % BIDS non-compatible logfile folder
    expParam.outputDir = fullfile(...
        fileparts(mfilename('fullpath')), ...
        'output');
    
end

%% audio other parameters
% sampling rate
cfg.fs = 44100; 
% channels is 1 for mono sound or 2 for stereo sound
cfg.audio.channels = 2;

%% download missing stimuli (.wav)
checkSoundFiles();

%% Timing 
% % %
% convert waitSecs according to the TR
% TR = 2.28s
% waitsecs(3*TR);
expParam.onsetDelay = 3; %Number of seconds before the rhythmic sequence (exp) are presented
expParam.endDelay = 0; % Number of seconds after the end of all stimuli before ending the fmri run! 
% % %

% ending timings for fMRI
expParam.endResponseDelay = 13; % wait for participant to response for counts
expParam.endScreenDelay = 2; %end the screen after thank you screen

% these are for behavioral exp delays
expParam.sequenceDelay = 1; %wait in between sequences? y/n
expParam.pauseSeq = 1; % give a pause of below seconds in between sequences


% define ideal number of sequences to be made
if cfg.debug
    expParam.numSequences = 2; % multiples of 3
else
    expParam.numSequences = 6;
end


%% fMRI task
% if fmriTask == true, it'll display a fixation cross during the fMRI run

if expParam.fmriTask
    
    % Used Pixels here since it really small and can be adjusted during the experiment
    expParam.fixCrossDimPix               = 10;   % Set the length of the lines (in Pixels) of the fixation cross
    expParam.lineWidthPix                 = 4;    % Set the line width (in Pixels) for our fixation cross
    expParam.xDisplacementFixCross        = 0;    % Manual displacement of the fixation cross
    expParam.yDisplacementFixCross        = 0;    % Manual displacement of the fixation cross
    expParam.fixationCrossColor           = cfg.white;
    
    %calculate the location coord for cross
    cfg.xCoords = [-expParam.fixCrossDimPix expParam.fixCrossDimPix 0 0] ...
        + expParam.xDisplacementFixCross;
    cfg.yCoords = [0 0 -expParam.fixCrossDimPix expParam.fixCrossDimPix] ...
        + expParam.yDisplacementFixCross;
    cfg.allCoords = [cfg.xCoords; cfg.yCoords];
    
end

%% more parameters to get according to the type of experiment
if strcmp(expParam.task,'tapTraining')
    
    % get tapping training parameters
    [cfg,expParam] = getTrainingParameters(cfg,expParam);
    
elseif strcmp(expParam.task,'tapMainExp')
    
    % get main experiment parameters
    [cfg,expParam] = getMainExpParameters(cfg,expParam);
    
end




%% differentiating response button (subject) from keyboard(experimenter)
% cfg.responseBox would be the device used by the participant to give his/her response: 
%   like the button box in the scanner or a separate keyboard for a behavioral experiment
%
% cfg.keyboard is the keyboard on which the experimenter will type or press the keys necessary 
%   to start or abort the experiment.
%   The two can be different or the same.

% Using empty vectors should work for linux when to select the "main"
%   keyboard. You might have to try some other values for MacOS or Windows

% % %

% After connecting fMRI response button to the laptop, 
% LOOK what are the experimenters' keyboard & fmri responseKey 
% by using GetKeyboardIndices below.
cfg.keyboard = [];
cfg.responseBox = [];

% % %


[cfg.keyboardNumbers, cfg.keyboardNames] = GetKeyboardIndices;
cfg.keyboardNumbers
cfg.keyboardNames 


switch lower(cfg.device)
    
    
    % this part might need to be adapted because the "default" device
    % number might be different for different OS or set up
    
    case 'pc'
        
        cfg.keyboard = [];
        cfg.responseBox = [];
        
        %behavioral exp keys to check
        cfg.keywait         = KbName({'RETURN'}); % press enter to start bloc
        cfg.keyToggleInstr  = KbName({'I'}); % press I to show/remove general instructions from the screen
        cfg.keytap          = KbName('SPACE');
        cfg.keyVolUp        = KbName('UpArrow');
        cfg.keyVolDown      = KbName('DownArrow');
        cfg.keyAudioPlay    = KbName('p');
        cfg.keyAudioStop    = KbName('s');
        cfg.keyInstrBack    = KbName('b');
        cfg.keyInstrNext    = KbName('n');
        
        if ismac
            cfg.keyboard = [];
            cfg.responseBox = [];
        end
        
    case 'scanner'
    % do nothing because in line 177-178 should have been assigned to keys
    % check
        
    otherwise
        
        cfg.keyboard = [];
        cfg.responseBox = [];
        
end


%%
if cfg.debug
    fprintf('\n\n\n\n')
    fprintf('######################################## \n')
    fprintf('##  DEBUG MODE, NOT THE ACTUAL EXP CODE  ## \n')
    fprintf('######################################## \n\n')    
end



