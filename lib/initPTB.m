function cfg = initPTB(cfg)


% check for octave:
if IsOctave
    checkOctave()
end


% check for OpenGL compatibility, abort otherwise:
AssertOpenGL;


%%
PsychJavaTrouble;

% Make sure keyboard mapping is the same on all supported operating systems
% Apple MacOS/X, MS-Windows and GNU/Linux:
KbName('UnifyKeyNames');


% Don't echo keypresses to Matlab window
% enable listening & additional keypress will be suppressed in command window
% use CTRL+C to reenable keyboard input when necessary
ListenChar(-1);



%% mouse
% Hide the mouse cursor:
HideCursor;
% Type ShowCursor if it's lost ;-)

%% init Visual

cfg.screen             = max(Screen('Screens'));

% Open a fullscreen, onscreen window with gray background. Enable 32bpc
% floating point framebuffer via imaging pipeline on it.
PsychImaging('PrepareConfiguration');

% init PTB with different options in concordance to the Debug Parameters
if cfg.debug

    % set to one because we don not care about time
    Screen('Preference', 'SkipSyncTests', 2);
    Screen('Preference', 'Verbosity', 0);
    Screen('Preferences', 'SuppressAllWarnings', 2);

    PsychDebugWindowConfiguration
    [cfg.win, cfg.winRect] = PsychImaging('OpenWindow', cfg.screen, cfg.backgroundColor);
        
else
    % we do not need high accuracy for the screen atm
    % Screen('Preference', 'SkipSyncTests', 1);
    Screen('Preference','SkipSyncTests', 0);
    [cfg.win, cfg.winRect] = PsychImaging('OpenWindow', cfg.screen, cfg.backgroundColor);

end



[cfg.winWidth,cfg.winHeight] = Screen('WindowSize',cfg.win);
%[cfg.winWidth,cfg.winHeight] = WindowSize(cfg.win);

% Get the Center of the Screen
cfg.center = [cfg.winRect(3), cfg.winRect(4)]/2;

% % % CONSIDER CHANGES ACCORDUNG TO THE FMRI MONITOR

% % %

% Enable alpha-blending, set it to a blend equation useable for linear
% superposition with alpha-weighted source.
Screen('BlendFunction',cfg.win,GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);

%% text
Screen('TextFont',cfg.win,cfg.textFont);
Screen('TextSize',cfg.win,cfg.textSize);


%% Timing
% Set priority for script execution to realtime priority:
Priority(MaxPriority(cfg.win));

%% init Audio


InitializePsychSound(1);

if any(strcmp(cfg.stimComp,{'mac','linux'}))
    
    
    % pahandle = PsychPortAudio('Open' [, deviceid][, mode][, reqlatencyclass][, freq] ...
    %       [, channels][, buffersize][, suggestedLatency][, selectchannels][, specialFlags=0]);
    % Try to get the lowest latency that is possible under the constraint of reliable playback
    cfg.pahandle = PsychPortAudio('Open', [], [], 3, cfg.fs, cfg.audio.channels);
    
    
else
    
    % get audio device list
    audio_dev       = PsychPortAudio('GetDevices');
    
    % find output device using WASAPI deiver
    idx             = find([audio_dev.NrInputChannels] == 0 & ...
        [audio_dev.NrOutputChannels] == 2 & ...
        ~cellfun(@isempty, regexp({audio_dev.HostAudioAPIName},'WASAPI')));
    
    % save device ID
    cfg.audio.i     = audio_dev(idx).DeviceIndex;
    
    % get device's sampling rate
    cfg.fs          = audio_dev(idx).DefaultSampleRate;
    
    % the latency is not important - but consistent latency is! Let's try with WASAPI driver.
    cfg.pahandle    = PsychPortAudio('Open', cfg.audio.i, 1, 3, cfg.fs, cfg.audio.channels);
    % cfg.pahandle = PsychPortAudio('Open', [], [], 0, cfg.fs, cfg.audio.channels);
    
end

% set initial PTB volume for safety (participants can adjust this manually
% at the begining of the experiment)
PsychPortAudio('Volume', cfg.pahandle, cfg.PTBInitVolume);

cfg.audio.pushsize  = cfg.fs*0.010; %! push N ms only
cfg.requestoffsettime = 1; % offset 1 sec
cfg.reqsampleoffset = cfg.requestoffsettime*cfg.fs; %


% playing parameters

% sound repetition
cfg.PTBrepet = 1;

% Start immediately (0 = immediately)
cfg.PTBstartCue = 0;

% Should we wait for the device to really start (1 = yes)
cfg.PTBwaitForDevice = 1;


end

