

% Clear all the previous stuff
if ~ismac
    close all;
    clear Screen;
else 
    clc; clear;
end

% make sure we got access to all the required functions and inputs
addpath(genpath(fullfile(pwd, 'lib')))


% Get parameters

% % %
% provide here to specifiy the sequence length and/or call getParamsMainExp
% in here with the arguments for rhythmic sequence, control1 and control2
% designs
% % %
[cfg,expParam] = getParams('tapMainExp');

% set and load all the subject input to run the experiment
expParam = userInputs(cfg,expParam);
expParam = createFilename(cfg,expParam);



% get time point at the beginning of the script (machine time)
expParam.scriptStartTime = GetSecs();

%% Experiment

% Safety loop: close the screen if code crashes
try
    % Init the experiment
    [cfg] = initPTB(cfg);

    
    
    % Prepare for the output logfiles - BIDS
    % saving 2 arrays long-form
    % open events logfile
    logFile  = saveEventsFile('open', expParam,[],'sequenceNum',...
        'segmentNum','segmentOnset','stepNum','stepOnset','patternID',...
        'category','F0','gridIOI','patternAmp','PE4_01', 'PE4_02',...
        'PE4_03', 'PE4_04', 'PE4_05', 'PE4_06','PE4_07', 'PE4_08',...
        'PE4_09', 'PE4_10', 'PE4_11', 'PE4_12','minPE4','rangePE4',...
        'LHL24_01', 'LHL24_02', 'LHL24_03', 'LHL24_04', 'LHL24_05',...
        'LHL24_06','LHL24_07', 'LHL24_08', 'LHL24_09', 'LHL24_10',...
        'LHL24_11', 'LHL24_12', 'minLHL24','rangeLHL24');
    
    % open stimulation logfile
    countFile  = saveEventsFile('open_stim', expParam,[],'target',...
        'key_name','pressed');


    % Show instructions for fMRI task
    if expParam.fmriTask
        displayInstr(expParam.fmriTaskInst,cfg);
    end
    
    % wait for space key to be pressed by the experimenter
    % to make the script more verbose
    pressSpace4me
    
    % prepare the KbQueue to collect responses
    getResponse('init', cfg, expParam);
    getResponse('start',cfg,expParam);
    
    % wait for trigger from fMRI
    wait4Trigger(cfg);
    
    % show fixation cross 
    if expParam.fmriTask
        drawFixationCross(cfg,expParam, expParam.fixationCrossColor);
        Screen('Flip',cfg.win);
    end
    

    % wait for dummy fMRI scans
    % and collect the timestamp
    expParam.experimentStart = GetSecs;
    
    WaitSecs(expParam.onsetDelay);
    
    
    
    %% play sequences
    for seqi = 1:expParam.numSequences

        currSeq = struct();
        responseEvents = struct();
        

        % construct sequence
        currSeq = makeSequence(cfg,seqi);

        
        % ===========================================
        % stimulus save for BIDS
        % ===========================================
        % we save sequence by sequence so we clear this variable every loop
        currSeq(1).fileID = logFile.fileID;
        
        % adding columns in currSeq for BIDS format
        for iPattern = 1:length(currSeq)
            currSeq(iPattern,1).trial_type  = 'dummy';
            currSeq(iPattern,1).duration    = 0;
            currSeq(iPattern,1).sequenceNum = seqi;            
        end
        
        saveEventsFile('save', expParam, currSeq,'sequenceNum',...
        'segmentNum','segmentOnset','stepNum','stepOnset','patternID',...
        'segmCateg','F0','gridIOI','patternAmp','PE4','minPE4',...
        'rangePE4','LHL24','minLHL24','rangeLHL24');


        %% present stimulus, accidential button press during the sequence

        % response save for BIDS (set up)
        responseEvents(1).fileID = logFile.fileID;            

        % fill the buffer
        PsychPortAudio('FillBuffer', cfg.pahandle, [currSeq.outAudio;currSeq.outAudio]);

        % start playing
        currSeqStartTime = PsychPortAudio('Start', cfg.pahandle, cfg.PTBrepet,...
            cfg.PTBstartCue, cfg.PTBwaitForDevice);

        
        % keep collecting tapping until sound stops (log as you go)
        expParam.seqi = seqi;
        expParam.currSeqStartTime = currSeqStartTime;
        
        % record response in case accidential press
        [tapOnsets, responseEvents] = mb_getResponse(cfg, ...
            expParam, ...
            responseEvents, ...
            currSeq);
        
        
        % response save for BIDS (write)
        if isfield(responseEvents,'onset')

            saveEventsFile('save', expParam, responseEvents,'sequenceNum',...
                'segmentNum','segmentOnset','stepNum','stepOnset','patternID',...
                'segmCateg','F0','gridIOI','patternAmp','PE4','minPE4',...
                'rangePE4','LHL24','minLHL24','rangeLHL24');

        end
        
        % ===========================================
        % log everything into matlab structure
        % ===========================================

        % save (machine) onset time for the current sequence
        expParam.data(seqi).currSeqStartTime = currSeqStartTime;

        % save PTB volume
        expParam.data(seqi).ptbVolume = PsychPortAudio('Volume',cfg.pahandle);

        % save current sequence information (without the audio, which can
        % be easily resynthesized)
        currSeq(1).outAudio = [];
        expParam.data(seqi).seq = currSeq;

        % save all the taps for this sequence
        expParam.data(seqi).taps = tapOnsets;


    end % sequence loop


    % wait while fMRI is ongoing
    WaitSecs(expParam.endDelay);

    % % %
    % ask for the button press tot times at the end if fmri run & give visual feedback?
    % % %
    displayInstr('Please indicate by pressing button, how many times you detected pitch changes\n\n\n',cfg);
    

    % flush the previous button presses
    %getResponse('flush', cfg, expParam);
    
    % wait 3 seconds for participant to press button
    WaitSecs(13);
    
    % write down buffered responses
    countEvents = getResponse('check', cfg, expParam);

    if ~isempty(countEvents(1).onset)
        
        countEvents(1).fileID = countFile.fileID;
        
        for iResp = 1:size(countEvents,1)
            countEvents(iResp).onset = countEvents(iResp).onset - expParam.experimentStart;
                countEvents(iResp).target = 8; % assign the correct target number

        end
        
        % saving in the tsv file
        saveEventsFile('save', expParam, countEvents,...
            'target','key_name','pressed');
    end
    
    
    

    
    % % make a if loop for the finaly run: 
    if expParam.runNb == 666 %change this with the known final run#
        displayInstr('DONE. \n\n\nTHANK YOU FOR PARTICIPATING :)\n\n\n Soon we will take you out!',cfg);
    else
        displayInstr('This run is over. We will shortly start the following!',cfg);
    end
    

    % wait 2 seconds for ending the screen/exp
    WaitSecs(2);
    
    % Close the logfiles (tsv)   - BIDS
    saveEventsFile('close', expParam, logFile);
    saveEventsFile('close', expParam, countFile);

    
    
    % save the whole workspace 
    matFile = fullfile(expParam.outputDir, strrep(expParam.fileName.events,'tsv', 'mat'));
    if IsOctave
        save(matFile, '-mat7-binary');
    else
        save(matFile, '-v7.3');
    end
    
    
    % clean the workspace
    cleanUp(cfg);



catch

    % save everything into .mat file
    matFile = fullfile(expParam.outputDir, strrep(expParam.fileName.events,'tsv', 'mat'));
    if IsOctave
        save(matFile, '-mat7-binary');
    else
        save(matFile, '-v7.3');
    end
    
    % Close the logfiles - BIDS
    saveEventsFile('close', expParam, logFile);
    saveEventsFile('close', expParam, countFile);



    % clean the workspace
    cleanUp(cfg);

    psychrethrow(psychlasterror);
end
