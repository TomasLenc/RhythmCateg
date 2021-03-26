function [whereIsData] = tapMainExperiment(task, subjectNb, runNb)
% TapMainExperiment script which runs the exp and present auditory sequence
% and records the tapping


% make sure we got access to all the required functions and inputs
initEnv();

cfg.subject.subjectGrp = '';
cfg.subject.subjectNb = subjectNb;
cfg.subject.sessionNb = 1;
cfg.subject.runNb = runNb;
    
 
% Define the task = 'RhythmFT', 'RhythmBlock'
% Get task specific parameters by providing task name
cfg = getParams(task,cfg);

%% Experiment

% Safety loop: close the screen if code crashes
try

  % get time point at the beginning of the experiment (machine time)
  cfg.experimentStartTime = GetSecs();

  % Init the experiment
  [cfg] = initPTB(cfg);

  % instructions
  if cfg.subject.runNb == 1 
    % show instructions and do initial volume setting
    cfg = setVolume(cfg);

    % more instructions
    displayInstr(cfg.trialDurInstruction, cfg, 'setVolume');
    
  else
      fbkToDisp = sprintf(cfg.generalDelayInstruction, ...
                          (cfg.subject.runNb-1), ...
                          cfg.pattern.numSequences);
      
      displayInstr(fbkToDisp, cfg, ...
                   'setVolumeToggleGeneralInstr', ...
                   'generalInstrTxt', cfg.generalInstruction);
  end
  
  % change screen to "GET READY" instruction
  displayInstr('GET READY', cfg);
  
  % if there's wait time,..wait
  WaitSecs(cfg.timing.startDelay);

  % get a list of sequence indices to run (from user defined runNb to the
  % end) 
 % seqToRun = [cfg.subject.runNb:cfg.pattern.numSequences]; 
  
  %% play sequences in a loop
  for iSequence = cfg.subject.runNb %seqToRun
    
    % set run to the current iSequence
    cfg.subject.runNb = iSequence; 
    cfg.iSequence = iSequence; 

    cfg = createFilename(cfg);

    logFile = struct();

    % create events with extra columns to save - BIDS
    logFile(1).extraColumns = cfg.extraColumns;
    
    % dummy initialise
    logFile = saveEventsFile('init', cfg, logFile); 
    
    % set the real length of columns
    logFile(1).extraColumns.LHL24.length = 12;
    logFile(1).extraColumns.PE4.length = 12;
    
    % actual inititalization
    logFile = saveEventsFile('open', cfg, logFile);    
    
    % construct sequence
    currSeq = makeSequence(cfg, iSequence);

    % change screen to "TAP" instruction
    displayInstr('TAP', cfg, 'instrAndQuitOption');


    %% present stimulus, record tapping
                               
    trigVal = cfg.beh.trigTaskMapping(task);  
    
    trigChan = cfg.audio.trigChanMapping(trigVal); 
    
    % play sound
    [tapData, sequenceOnset, trialTerminated] = playSound(currSeq(1).outAudio, ...
                                                        cfg.audio.fs, ...
                                                        cfg.audio.pahandle, ...
                                                        cfg.audio.channels(1), ...
                                                        cfg.audio.channels(2), ...
                                                        trigChan, ...
                                                        cfg.audio.initPushDur, ...
                                                        cfg.audio.pushDur, ...
                                                        cfg.keyboard.quit); 
       
    % wait for playback end                                                
    PsychPortAudio('Stop', cfg.audio.pahandle, 1);  
    
    % show quick message to the participant 
    displayInstr('... saving data, please wait ...', cfg)
    
    %% log data 
                                                                        
    % ===========================================
    % BIDS
    % ===========================================            
    collectAndSaveEvents(cfg, logFile, currSeq, iSequence, sequenceOnset, task);
    
    % ---- SAVE tapping as stim file ---- 
    
    % downsample tapping
    [P,Q] = rat(cfg.audio.fsDs/cfg.audio.fs); 
    tapDataDs = num2cell(resample(tapData(1,:), P, Q)); 
    
    % get stimulus envelope and downsample 
    env = abs(hilbert(tapData(2,:))); 
    soundDataDs = num2cell(resample(env, P, Q)); 
            
    stimFile = []; 
    
    % set columns for data
    stimFile(1).extraColumns = {'tapForce', 'sound'};

    stimFile(1).StartTime = 0; 
    stimFile(1).SamplingFrequency = cfg.audio.fsDs; 

    % init stim file
    stimFile = saveEventsFile('init_stim', cfg, stimFile); 

    % output directory (created automatically)
    stimFile = saveEventsFile('open', cfg, stimFile); 
    
    % prepre data
    [stimFile(1:length(tapDataDs),1).tapForce] = tapDataDs{:}; 
    [stimFile(1:length(tapDataDs),1).sound] = soundDataDs{:}; 

    % write
    stimFile = saveEventsFile('save', cfg, stimFile);
    
    % Close the logfiles (tsv)   - BIDS
    saveEventsFile('close', cfg, logFile);
  
    % ===========================================
    % audiofile
    % ===========================================
    fileNameTap = [cfg.fileName.base, ...
                   cfg.fileName.suffix.run, ...
                   '_recording-tapping_physio', ...
                   '_date-',  cfg.fileName.date, ...
                   '.wav']; 

    % save tapping with original sampling rate          
    audiowrite(fullfile(cfg.dir.outputSubject, fileNameTap),...
               tapData', ...
               cfg.audio.fs); % 'BitsPerSample' -> only 8 bit resolution to save faster...

    %% 
    % show pause screen in between sequences
    showPauseScreen(cfg);

  end

  % bids files were saved in
  whereIsData = fullfile( cfg.dir.outputSubject, ...
            cfg.fileName.modality, ...
            logFile.filename);

  % save the whole workspace
  currSeq = []; 
  env = []; 
  tapData = []; 
  tapDataDs = []; 
  soundDataDs = []; 
  
  matFile = fullfile(cfg.dir.output, ...
                     strrep(cfg.fileName.events, 'tsv', 'mat'));
  if IsOctave
    save(matFile, '-mat7-binary');
  else
    save(matFile, '-v7.3');
  end

  % clean the workspace
  cleanUp();

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
  cleanUp();
  psychrethrow(psychlasterror);
end
