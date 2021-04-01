

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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


% Define the task = 'RhythmFT', 'RhythmBlock'
% Get task specific parameters by providing task name
cfg = struct();
cfg = getParams('tapTraining',cfg);

% datalogging structure
datalog = [];


% get time point at the beginning of the experiment (machine time)
datalog.experimentStartTime = GetSecs();


try
  [cfg] = initPTB(cfg);

  % Prepare for the output logfiles
  datalog = saveOutput(cfg, 'open');

  % show instructions and do initial volume setting
  currInstrPage = 1;
  nInstrPages = length(cfg.taskInstruction);
  while 1
    % display instructions and wait for action
    subAction = displayInstr(cfg.taskInstruction{currInstrPage}, cfg, 'setVolumePrevNext', ...
                             'currInstrPage', currInstrPage, ...
                             'nInstrPages', nInstrPages);
    % go one instruction page forward or backward (depending on subject's action)
    if strcmp(subAction, 'oneInstrPageForward')
      currInstrPage = min(currInstrPage + 1, length(cfg.taskInstruction));
    elseif strcmp(subAction, 'oneInstrPageBack')
      currInstrPage = max(currInstrPage - 1, 1);
    elseif strcmp(subAction, 'done')
      break
    end
  end

  % simultenaous feedback
  fbkOnScreen = false;

  % index (counter) of current pattern that is used in the stimulus
  % sequence
  currPatterni = 4;

  %% loop over patterns (atm, n pattern = $)
  while 1

    if ischar(cfg.patterns{currPatterni})
      currPatternStr = cfg.patterns{currPatterni};
    else
      currPatternStr = num2str(cfg.patterns{currPatterni}, '%d');
    end

    % once the pattern is repeated 4 times, the script looks back in
    % this time window, this is an index (counter) of analysis windows
    % for the current sequence
    currWini = 1;

    % dB changes in the loop - tapping cue (decreases over time if the
    % error rate is low(er))
    % which dB sound will be played choose the index in the list
    cueDBleveli = 1;

    % number of grip points between two tapping cue sounds
    % convert the grip interval to target inter-tap-interval in seconds.
    % every 800ms there'll be a cue sound to tap along
    cuePeriodTime = cfg.cuePeriodGrid(currPatterni) * cfg.gridIOI(currPatterni);

    % total duration of the analysis window (4 x pattern) in seconds
    winDur = cfg.winDur(currPatterni);

    % time from the start of the analysis window (in secs) taken to analyse the tapping
    taketapDur = winDur - 0.200;

    % to be used later on to calculate the min number that participant
    % has to tap -
    maxPossibleNtaps = floor(taketapDur / cuePeriodTime);

    % minimum required number of taps in an analysis window
    % (70% of max possible taps in the window)
    minNtaps = floor(maxPossibleNtaps * cfg.minNtapsProp);

    % counters for calculating the tapping accuracy later on
    performStatus = 0; % if good +1, if bad -1
    taps = [];
    istap = false;

    % counter for audio-track index (used only for tracks, but we can
    % pass it to makeStimTrain all the time)
    soundIdx = 0;

    %% allocate datalog variables
    winIdxs            = [];
    cueDBs             = [];
    winStartTimes      = [];
    feedbacks          = {};
    instructions       = {};

    %% make stimuli
    % get audio for the first step/window (4 x pattern)
    [seq] = makeStimTrain(cfg, currPatterni, cueDBleveli, currWini, soundIdx);

    % update audio index (only relevant for audio tracks)
    soundIdx = seq.idxEnd;

    %% display sequence-specific instructions

    % show instructions before the sequence
    currInstrPage = 1;
    nInstrPages = length(cfg.beforeSeqInstruction{currPatterni});
    while 1
      % display instructions and wait for action
      subAction = displayInstr(cfg.beforeSeqInstruction{currPatterni}{currInstrPage}, cfg, 'setVolumePrevNext', ...
                               'currInstrPage', currInstrPage, ...
                               'nInstrPages', nInstrPages);
      % go one instruction page forward or backward (depending on subject's action)
      if strcmp(subAction, 'oneInstrPageForward')
        currInstrPage = min(currInstrPage + 1, length(cfg.beforeSeqInstruction{currPatterni}));
      elseif strcmp(subAction, 'oneInstrPageBack')
        currInstrPage = max(currInstrPage - 1, 1);
      elseif strcmp(subAction, 'done')
        break
      end
    end

    % during sequence (part 1)
    if currWini <= cfg.nWinNoCue(currPatterni)
      % if part 1 of the sequence
      % look for specific instructions and task
      instr2disp = cfg.duringSeqInstruction_part1{currPatterni};
      currTask = cfg.taskPart1{currPatterni};
      % if not available, use default
      if isempty(instr2disp)
        instr2disp = 'LISTEN';
      end
      if isempty(currTask)
        currTask = 'listen';
      end
    else
      % if part 2 of the sequence
      % look for specific instructions and task
      instr2disp = cfg.duringSeqInstruction_part2{currPatterni};
      currTask = cfg.taskPart2{currPatterni};
      % if not available, use default
      if isempty(instr2disp)
        instr2disp = 'TAP';
      end
      if isempty(currTask)
        currTask = 'tap';
      end
    end

    % display the instruction on the screen
    displayInstr(instr2disp, cfg, 'instrAndQuitOption');

    %% fill the buffer
    % first, fill the buffer with 60s silence

    % to get a buffer longer than what you are pushing
    % allocate the buffer

    % if case some stays too long in the while loop, we will need this
    % buffer to allocate
    PsychPortAudio('FillBuffer', cfg.audio.pahandle, zeros(2, 60 * cfg.fs));

    %% start playback

    % start playback (note: set repetitions=0, otherwise it will not allow you to seamlessly push more data into the buffer once the sound is playing)
    % starts to play whats in the buffer and play on whatever is in on
    % a seamlessly in the loop
    currSeqStartTime = PsychPortAudio('Start', cfg.audio.pahandle, 0, [], 1);
    % startTime = PsychPortAudio('Start', pahandle [, repetitions=1] [, when=0] [, waitForStart=0] [, stopTime=inf] [, resume=0]);

    % 1 sound input into 1 channels also works
    audio2push = [seq.s; seq.s];

    % silence is going, then we will upload to the buffer audio sequence after the
    % 1s of silent has started
    [underflow] = PsychPortAudio('FillBuffer', cfg.audio.pahandle, audio2push, 1, cfg.audio.requestSampleOffset);

    % and update start time (by offset)
    % start time = actual time of audio seq presented
    currSeqStartTime = currSeqStartTime + cfg.audio.requestTimeOffset;
    currWinStartTime = currSeqStartTime;
    nSamplesAudio2push = 0;
    idx2push = 1;

    %% loop over pattern windows (in which we may change dB levels atm)
    while 1

      %% tapping while loop
      while GetSecs < (currWinStartTime + taketapDur)

        % collect tapping
        [~, tapOnset, keyCode] = KbCheck(-1);

        % terminate if quit-button pressed
        if find(keyCode) == cfg.keyboard.quit
          error('Experiment terminated by user...');
        end

        % if they did not press delete, it looks for any response
        % button and saves the time
        if ~istap && any(keyCode)
          taps = [taps, tapOnset - currSeqStartTime];
          istap = true;

          % -------------------- log ----------------------------
          % now we have some time before they tap again so let's
          % write to the log file
          fprintf(datalog.fidTapTrainer, '%s\t%d\t%s\t%s\t%f\t%f\t%f\t%d\t%f\t%f\n', ...
                  cfg.subject.subjectNb, ...               % subject id
                  currPatterni, ...                   % pattern
                  currPatternStr, ...                 % name of the current pattern/track
                  currTask, ...                      % instruction
                  currSeqStartTime, ...                % machine time of sequence audio start
                  cuePeriodTime, ...                   % cue (i.e. metronome) period (N of grid-points)
                  seq.cueDB, ...                       % cue (i.e. metronome) level in dB (SNR)
                  currWini, ...                       % index (count) of this analysis window (for this sequence)
                  currWinStartTime - currSeqStartTime, ...  % analysis window start time wrt sequence start
                  taps(end));                         % tap onset time relative to sequence start time
          % -----------------------------------------------------
        end

        % it counts as tap if reponse buttons were released
        % initially
        if istap && ~any(keyCode)
          istap = false;
        end

        % if there is any audio waiting to be pushed, push it to the buffer!
        if nSamplesAudio2push
          if idx2push + cfg.audio.pushSample > nSamplesAudio2push
            pushdata = audio2push(:, idx2push:end);
            nSamplesAudio2push = 0;
          else
            pushdata = audio2push(:, idx2push:idx2push + cfg.audio.pushSample - 1);
            idx2push = idx2push + cfg.audio.pushSample;
          end
          [curunderflow, ~, ~] = PsychPortAudio('FillBuffer', cfg.audio.pahandle, pushdata, 1);
        end

        % if there is overdue feedback on the screen, remove it
        if fbkOnScreen
          if (GetSecs - fbkOnScreenTime) > cfg.fbkOnScreenMaxtime
            fbkOnScreen = false;
            if currWini <= cfg.nWinNoCue(currPatterni)
              % if part 1 of the sequence
              % look for specific instructions and task
              instr2disp = cfg.duringSeqInstruction_part1{currPatterni};
              currTask = cfg.taskPart1{currPatterni};
              % if not available, use default
              if isempty(instr2disp)
                instr2disp = 'LISTEN';
              end
              if isempty(currTask)
                currTask = 'listen';
              end
            else
              % if part 2 of the sequence
              % look for specific instructions and task
              instr2disp = cfg.duringSeqInstruction_part2{currPatterni};
              currTask = cfg.taskPart2{currPatterni};
              % if not available, use default
              if isempty(instr2disp)
                instr2disp = 'TAP';
              end
              if isempty(currTask)
                currTask = 'tap';
              end
            end
            % display the instruction on the screen
            displayInstr(instr2disp, cfg, 'instrAndQuitOption');
          end
        end

      end % end of tapping while loop

      %%

      % let's look at the last window to analyse the tapping
      % evaluate tapping and decide on the next step/window parameters
      % look at which window you are in (curr_step_
      % taps you will consider : step_dur
      minTaketapTime = (currWini - 1) * winDur;
      maxTaketapTime = (currWini) * winDur;

      % we will only take the taps which are relevant - in our
      % current window
      % test test
      % curr_taps = [0:0.8:6.2] + randn(1,8)*0.1
      currTaps = taps(taps > minTaketapTime & taps < maxTaketapTime) - minTaketapTime;

      % making tapping vector
      disp(currTaps');
      currTapsN = length(currTaps);

      % round to the closest beat position
      % there are target positions set by the cue. But they can skip
      % a tap and still be good in synch. round the mto the closest
      % beat positions. that will give you the target position. and
      % you calculate the how far your target + your current taps
      targetTapTimes = round(currTaps / cuePeriodTime) * cuePeriodTime;

      % calculate asynchronies
      tapAsynch = currTaps - targetTapTimes;
      % first normalise it by the interval. because your variability
      % is proportional to the legnth of the interval (equalise across
      % different time interval/tempi)
      % then std it
      tapCvAsynch = std(tapAsynch / cuePeriodTime);

      %% update tapping performance
      % (wrt cvASY threshold and n-taps)
      % tap_perform_status = -1, 0, 1, 2, 3, ...(correctness)

      if currWini <= cfg.nWinNoCue(currPatterni)
        % this was a listening-only window, don't evalueate
        % performance
        currPerform = 'NA';
        performStatus = 0;

      elseif (tapCvAsynch < cfg.tapCvAsynchThr) && currTapsN >= minNtaps
        % good performance, one up!
        currPerform = 'good';
        performStatus = max(0, performStatus); % if negative make 0
        performStatus = performStatus + 1;
      else
        % bad performance, one down...
        currPerform = 'bad';
        performStatus = min(0, performStatus); % if positive make 0
        performStatus = performStatus - 1;
      end

      %% update datalog variables

      % index (count) of the current analysis window
      winIdxs            = [winIdxs, currWini];
      % cue dB level (SNR)
      cueDBs             = [cueDBs, cfg.cueDB(cueDBleveli)];
      % start time of the analysis window (wrt sequence start time)
      winStartTimes      = [winStartTimes, currWinStartTime - currSeqStartTime];
      % feedback for the current window (good/bad)
      feedbacks          = [feedbacks, currPerform];
      % instructions for the current window (listen/tap)
      instructions       = [instructions, currTask];

      %% update next window parameters
      % staircase here to adapt

      % if we've used all samples in the audio track
      % of we have reached time out
      if isfield(seq, 'AUDIO_END') | (GetSecs - currSeqStartTime) > cfg.timeOut(currPatterni)

        % stop the audio
        PsychPortAudio('Stop', cfg.audio.pahandle, 1);

        % end the loop over pattern windows (we will continue with
        % the following pattern after participant has a break)
        break

        % this was a listening-only window, don't evalueate performance
      elseif currWini <= cfg.nWinNoCue(currPatterni)

        if currWini + 1 <= cfg.nWinNoCue(currPatterni)
          % if part 1 of the sequence
          % look for specific instructions and task
          instr2disp = cfg.duringSeqInstruction_part1{currPatterni};
          currTask = cfg.taskPart1{currPatterni};
          % if not available, use default
          if isempty(instr2disp)
            instr2disp = 'LISTEN';
          end
          if isempty(currTask)
            currTask = 'listen';
          end
        else
          % if part 2 of the sequence
          % look for specific instructions and task
          instr2disp = cfg.duringSeqInstruction_part2{currPatterni};
          currTask = cfg.taskPart2{currPatterni};
          % if not available, use default
          if isempty(instr2disp)
            instr2disp = 'TAP';
          end
          if isempty(currTask)
            currTask = 'tap';
          end
        end
        % display the instruction on the screen
        displayInstr(instr2disp, cfg, 'instrAndQuitOption');

        % if this window was the last dB level, and the last-dB-level
        % counter is equal to the goal number
      elseif (cueDBleveli == cfg.nCueDB) && (performStatus == cfg.nWinUp_lastLevel)

        % stop the audio
        PsychPortAudio('Stop', cfg.audio.pahandle, 1);

        % end the loop over pattern windows (we will continue with
        % the following pattern after participant has a break)
        break

        % if we are not yet in the last level, and we have enough good
        % successive windows, we need to move one db-level up
      elseif (cueDBleveli ~= cfg.nCueDB) && (performStatus == cfg.nWinUp)

        % reset the performance counter to start next level from 0
        performStatus = 0;

        % increase the dB level one step up
        cueDBleveli = cueDBleveli + 1;

        % Give positive feedback.
        errStr = getErrorStr(tapCvAsynch, cfg);
        fbktxt = sprintf([sprintf('level %d out of %d.\n\n', cueDBleveli, cfg.nCueDB), ...
                          'your performance: \n', ...
                          errStr, ...
                          '\n\nWell done!']);
        displayInstr(instr2disp, cfg, 'fbktxt', fbktxt, 'instrAndQuitOption');
        fbkOnScreenTime = GetSecs;
        fbkOnScreen = 1;

        % disregarding which level you are in, if the last N successive steps
        % were bad (N == n_steps_down) -> decrease level
      elseif performStatus <= -cfg.nWinDown

        % reset the performance counter to start the next (decreased) level from 0
        performStatus = 0;

        % decrease the dB level one step down (don't change if it
        % is already at the lowest possible level)
        cueDBleveli = max(cueDBleveli - 1, 1);

        % Give negative feedback.
        errStr = getErrorStr(tapCvAsynch, cfg);
        fbktxt = [sprintf('level %d out of %d.\n\n', cueDBleveli, cfg.nCueDB), ...
                  'your performance: \n', ...
                  errStr, ...
                  '\n\nKeep trying :)'];
        displayInstr(instr2disp, cfg, 'fbktxt', fbktxt, 'instrAndQuitOption');
        fbkOnScreenTime = GetSecs;
        fbkOnScreen = 1;

        % otherwise just give feedback and continue...
      else
        errStr = getErrorStr(tapCvAsynch, cfg);
        fbktxt = [sprintf('level %d out of %d.\n\n', cueDBleveli, cfg.nCueDB), ...
                  'your performance: \n', ...
                  errStr, ...
                  sprintf('\n')];
        displayInstr(instr2disp, cfg, 'fbktxt', fbktxt, 'instrAndQuitOption');
        fbkOnScreenTime = GetSecs;
        fbkOnScreen = 1;
      end

      fprintf('current metronome-SNR level = %d\n', cueDBleveli);
      % ---- end of next window parameters update ----

      %% create audio with the new dB level for the next step/window

      % Update current time!
      % This will be used as the start time of the following window
      currWinStartTime = currWinStartTime + winDur;

      % update step counter (because we're already constructing for
      % the next window)
      currWini = currWini + 1;

      % get a new audio sequence for the next window
      [seq] = makeStimTrain(cfg, currPatterni, cueDBleveli, currWini, soundIdx);

      % update audio index (only relevant for audio tracks, dummy for patterns)
      soundIdx = seq.idxEnd;

      % check if we're pushing tha last samples in the audio track
      if isfield(seq, 'AUDIO_END')
        % if yes, update the window length so we end the  tapping while loop on time
        taketapDur = length(seq.s) / seq.fs;
      end

      % update the push variable
      audio2push = [seq.s; seq.s];
      nSamplesAudio2push = size(audio2push, 2);
      idx2push = 1;

    end % end of tapping loop

    % ====================== update datalog =============================

    % save (machine) onset time for the current sequence
    datalog.data(currPatterni).currSeqStartTime = currSeqStartTime;

    % save current sequence information (without the audio, which can
    % be easily resynthesized)
    datalog.data(currPatterni).seq      = seq;
    datalog.data(currPatterni).seq.s    = [];

    % save all the taps for this sequence
    datalog.data(currPatterni).taps = taps;

    % save PTB volume
    datalog.data(currPatterni).ptbVolume = PsychPortAudio('Volume', cfg.audio.pahandle);

    % save other window-level variables
    datalog.data(currPatterni).wini             = winIdxs;
    datalog.data(currPatterni).cueDBs           = cueDBs;
    datalog.data(currPatterni).winStartTimes    = winStartTimes;
    datalog.data(currPatterni).feedbacks        = feedbacks;
    datalog.data(currPatterni).instructions     = instructions;

    % ========================= instructions ===============================

    % show instructions after the sequence
    currInstrPage = 1;
    nInstrPages = length(cfg.afterSeqInstruction{currPatterni});
    while 1
      % display instructions and wait for action
      subAction = displayInstr(cfg.afterSeqInstruction{currPatterni}{currInstrPage}, cfg, 'setVolumePrevNext', ...
                               'currInstrPage', currInstrPage, ...
                               'nInstrPages', nInstrPages);
      % go one instruction page forward or backward (depending on subject's action)
      if strcmp(subAction, 'oneInstrPageForward')
        currInstrPage = min(currInstrPage + 1, length(cfg.afterSeqInstruction{currPatterni}));
      elseif strcmp(subAction, 'oneInstrPageBack')
        currInstrPage = max(currInstrPage - 1, 1);
      elseif strcmp(subAction, 'done')
        break
      end
    end

    if currPatterni == cfg.nPatterns
      % end of experient
      displayInstr('The training is over now. \n\n\nTHANK YOU!\n\n\nPlease continue to the Main Experiment.', cfg);
      % wait 3 seconds and end the experiment
      WaitSecs(3);
      break
    end

    % ========================= update counter ===============================

    % we will move on to the next pattern in the following loop iteration
    currPatterni = currPatterni + 1;

  end % end of loop over patterns

  saveOutput(cfg, 'savemat');
  saveOutput(cfg, 'close');

  cleanUp();
catch

  cleanUp();

  saveOutput(cfg, 'savemat');
  saveOutput(cfg, 'close');

  psychrethrow(psychlasterror);
end

