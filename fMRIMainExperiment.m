
%%
% Clear all the previous stuff
if ~ismac
    close all;
    clear Screen;
else
    clc;
    clear;
end

% make sure we got access to all the required functions and inputs
initEnv();

% Define the task = 'RhythmFT', 'RhythmBlock', 'Nonmetric'
% Get parameters by providing task name
cfg = getParams('Nonmetric');



%% Experiment

% Safety loop: close the screen if code crashes
try
    
    %% Init Experiment
    % get time point at the beginning of the script (machine time)
    cfg.timing.scriptStartTime = GetSecs();

    % Init the experiment
    [cfg] = initPTB(cfg);

    % create  logfile with extra columns to save - BIDS
    logFile.extraColumns = cfg.extraColumns;
    [logFile]  = saveEventsFile('init', cfg, logFile); %dummy initialise

    % set the real length of columns
    logFile(1).extraColumns.LHL24.length = 12;
    logFile(1).extraColumns.PE4.length = 12;

    % actual inititalization
    logFile = saveEventsFile('open', cfg, logFile);

    
    % Show instructions for fMRI task & wait for space press
    standByScreen(cfg);

    % prepare the KbQueue to collect responses
    getResponse('init', cfg.keyboard.responseBox, cfg);
    getResponse('start', cfg.keyboard.responseBox);

    % wait for trigger from fMRI
    cfg.experimentStart = waitForTrigger(cfg);
    
    %% Start Experiment
    % show fixation cross + get timestamp
    %cfg = getFixationCross(cfg);
    drawFixation(cfg);
    Screen('Flip', cfg.screen.win);

    % take the runNb corresponding sequence
    iSequence = cfg.subject.runNb;

    % prep for BIDS saving structures
    currSeq = struct();
    responseEvents = struct();

    % construct sequence
    currSeq = makeSequence(cfg, iSequence);
    
    
    %% play sequences
    % fill the buffer % start sound presentation
    PsychPortAudio('FillBuffer', cfg.audio.pahandle, ...
        [currSeq.outAudio;currSeq.outAudio]);
    
    % wait for baseline delays and then start the audio
    onset = PsychPortAudio('Start', cfg.audio.pahandle, [], ...
        cfg.experimentStart + cfg.timing.onsetDelay,1);
    
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
    
    % ===========================================
    % stimulus save for BIDS
    % ===========================================
    target = collectAndSaveEvents(cfg, logFile, currSeq,iSequence, onset);
   
    %% Wait for audio and delays to catch up
    % stay here till audio stops & check esc key press
    waitAndCheckEsc(cfg);

    %%
    % record exp ending time
    cfg.timing.fMRIDuration = GetSecs - cfg.experimentStart;

    %% Check button presses
    % instructions to press
    displayInstr(cfg.task.instructionPress, cfg);

    % wait for participant to press button
    WaitSecs(cfg.timing.endResponseDelay);

    % save response & target
    cfg.target = target;
    responseEvents = collectAndSave(cfg, ...
        logFile, cfg.experimentStart);
    

    %% wrapping up
    % last screen
    if cfg.debug.do || cfg.subject.runNb == cfg.pattern.numSequences
        displayInstr(cfg.task.instructionEnd, cfg);
    else
        displayInstr(cfg.task.instructionCont, cfg);
    end

    % wait for ending the screen/exp
    WaitSecs(cfg.timing.endScreenDelay);

    % record script ending time
    cfg.timing.scriptDuration = GetSecs - cfg.experimentStart;

    % clear the buffer & stop key checks
    getResponse('stop', cfg.keyboard.responseBox);
    getResponse('release', cfg.keyboard.responseBox);


    %% save
    % Close the logfiles (tsv)   - BIDS
    saveEventsFile('close', cfg, logFile);

    % save the whole workspace
    matFile = fullfile(cfg.dir.output, ...
        strrep(cfg.fileName.events, 'tsv', 'mat'));
    if IsOctave
        save(matFile, '-mat7-binary');
    else
        save(matFile, '-v7.3');
    end

    createJson(cfg, 'func');
    
    % clean the workspace
    cleanUp;

catch

    % save everything into .mat file
    matFile = fullfile(cfg.dir.output, ...
        strrep(cfg.fileName.events, 'tsv', 'mat'));
    if IsOctave
        save(matFile, '-mat7-binary');
    else
        save(matFile, '-v7.3');
    end


    % clean the workspace
    cleanUp;
    psychrethrow(psychlasterror);
end
