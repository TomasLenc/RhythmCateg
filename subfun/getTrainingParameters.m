function cfg = getTrainingParameters(cfg)

  % % % % % % % % % % % % % % % % % % % % % % % % % % %
  %
  % check if some parameters could be inserted  into
  % getParams.m script instead
  %
  % % % % % % % % % % % % % % % % % % % % % % % % % % %

  % there are 3 types is "sequence" that can be used (specified in cfg.patterns):
  %
  %     1) pattern:               1 cycle specified as a vector of ones (sounds) and zeros (silences)
  %
  %     2) audio track:           provide name of the audio file as string (must be saved in "stimuli"
  %                               folder in a .wav or .mp3 format)
  %
  %     3) nonrepeating sequence: use string 'GrahnComplex' to use complex patterns,
  %                               or 'GrahnSimple' to use simple patterns
  %
  % patterns/ audio tracks will be tried from the first to the last
  cfg.patterns = {'Flabaire-Alpha2(120BPM).mp3', ...
                  'YCreate-IDontWantToBe(1.364Hz-cut3min).mp3', ...
                  'Bugz_4_Hugz_Dub.mp3', ...
                  'GrahnComplex', ...
                  'GrahnComplex'};

  % number of patterns
  cfg.nPatterns   = length(cfg.patterns);

  % find which items are audio tracks
  cfg.isTrackIdx = find(cellfun(@(x)strcmp(class(x), 'char'), cfg.patterns));

  % grahn-like sequences are not audio tracks!
  cfg.isTrackIdx(find(~cellfun(@isempty, regexpi(cfg.patterns(cfg.isTrackIdx), 'grahn')))) = [];

  nTracks = length(cfg.isTrackIdx);

  %% tapping training parameters

  % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % ---------
  % IMPORTANT:    Below are pattern/track-specific parameters!
  % ---------     Make sure you change all of them accordingly, after
  %               changing the track or pattern order above!!!
  %

  % number of pattern cycles in each step/window of pattern: how many cycles
  % of the pattern will be repeated before tapping performance is evaluated
  % !!! for audio tracks, there are no cycles defined so  instead,
  % this is the number of beats in a window (beat duration = cuePeriodGrid *
  % gridIOI) !!!
  cfg.nCyclesPerWin = [8, 8, 8, 4, 4];

  % tapping cue sounds (metronome)
  cfg.cuePeriodGrid = [4, 4, 4, 4, 3]; % each pattern needs a metronome period assigned (units: N-grid-ticks)

  % time-interval of one grid-tick (IOI between events)
  % it's not necessarily the duration of the sound.
  % This needs to be set separately for each pattern (or track)
  cfg.gridIOI = [0.125, 1 / (1.364) / 4, 0.125, 0.190, 0.190];

  % max number of seconds the subject can spend on one pattern/track. If the
  % tapping is still not good after this time, the sequence will terminate
  % and continue to the next pattern/track.
  % Important: this defines duration of listen-only tracks! (otherwise
  % they would never end).
  cfg.timeOut = [180, 180, 120, 180, 120];

  % Each sequence (either pattern or audio track) has two parts.
  % ------ PART 1 ------
  % No metronome cue is presented.
  % Tapping is recorded, but not evalueated, and no feedback is
  % provided to the participant.

  % The duration of part 1 is defined by the parameter nWinNoCue, as the
  % number of windows from the begining of each pattern/track, where no cue
  % is presented
  cfg.nWinNoCue = [4, 6, 4, 4, 0];

  % instruction for logging only for part 1 (if empty, default is
  % used: 'listen')
  cfg.taskPart1 = {'tap_alone', ...
                   'tap_alone', ...
                   'tap_alone', ...
                   'tap_alone', ...
                   'tap_with_bass'};

  % ------ PART 2 ------
  % Metronome cue is presented at the dB level defined by participants prior
  % tapping performance. Feedback about tapping is provided at the end of
  % each analysis window.
  % The duraion of part 2 is till the end of the current sequence.

  % instruction for logging only for part 2 (if empty, default is
  % used: 'tap')
  cfg.taskPart2 = {'tap_with_claps', ...
                   'tap_with_claps', ...
                   'tap_with_claps', ...
                   'tap_with_bass', ...
                   'tap_with_bass'};

  % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % decreasing the DB along with the high accuracy of tapping
  cfg.cueDB = [0, -15, -Inf]; % [0, -14, -25, -Inf] SNRs between rhythm and metronome audio (to be used across levels)

  % to calculate how many difficulty levels there (in dB)
  % or number of SNR-levels
  cfg.nCueDB = length(cfg.cueDB);

  % duration (in seconds) of the each step/window of pattern
  % (note: the window duration for audio tracks will be changed later in this
  % script, the vale assigned here is not valid)
  cfg.winDur = cfg.nCyclesPerWin .* cellfun(@length, cfg.patterns) .* cfg.gridIOI;

  % threshold for coefficient of variation (cv) of tap-beat asynchronies (defines good/bad tapping performance in each step)
  % the taps to be registered as correct "it" should be below the threshold
  % within the 4 cycle of pattern representation the error is calculated and
  % this is the shift of tapping variation should be in the range of 10%
  cfg.tapCvAsynchThr = 0.160;

  % minimum N taps for the step/cycle/window to be valid (units: proportion of max. possible N taps considering the beat period)
  % if they tapped less than 70% of the maximum possible number of taps
  % it would be a bad trial below this point (70%)
  cfg.minNtapsProp = 0.7;

  %% staircase parameters
  % how many trials/windows do you need to go through staircase procedure
  % N successive steps that need to be "good tapping" to move one SNR level up
  cfg.nWinUp = 2;

  % N successive steps that need to be "bad tapping" to move one SNR level down
  cfg.nWinDown = 1;

  % N successive steps/windows that need to be "good tapping" for the final level to finish
  % this is in the  last level (level cfg.nCueDB)
  cfg.nWinUp_lastLevel = 3;

  % duration (secs) for which real-time feedback will be displayed on screen during tapping
  cfg.fbkOnScreenMaxtime = Inf;

  %% load wav files to make sounds

  % load audio samples
  soundPattern    = audioread(fullfile('.', 'stimuli', 'tone440Hz_10-50ramp.wav')); % rimshot_015
  soundPattern    = 1 / 3 * soundPattern; % set amplitude to 1/3 to prevent clipping after adding pattern+metronome

  soundBeat       = audioread(fullfile('.', 'stimuli', 'Kick8.wav'));
  soundBeat       = mean(soundBeat, 2); % average L and R channels
  soundBeat       = 1 / 3 * soundBeat; % set amplitude to 1/3 to prevent clipping after adding pattern+metronome

  soundGrid       = audioread(fullfile('.', 'stimuli', 'Perc5_cut.wav'));
  soundGrid       = mean(soundGrid, 2); % average L and R channels
  soundGrid       = 1 / 3 * soundGrid; % set amplitude to 1/3 to prevent clipping after adding pattern+metronome

  soundBeatTrack  = audioread(fullfile('.', 'stimuli', 'clap_005.wav'));
  soundBeatTrack  = mean(soundBeatTrack, 2); % average L and R channels
  soundBeatTrack  = 1 / 3 * soundBeatTrack; % set amplitude to 1/3 to prevent clipping after adding pattern+metronome

  % equalize RMS
  rmsPat          = rms(soundPattern);
  rmsBeat         = rms(soundBeat);
  rmsGrid         = rms(soundGrid);
  maxAllowedRms   = min([rmsPat, rmsBeat, rmsGrid]);

  cfg.soundPattern    = soundPattern / rmsPat * maxAllowedRms;
  cfg.soundBeat       = soundBeat / rmsBeat * maxAllowedRms;
  cfg.soundGrid       = soundGrid / rmsGrid * maxAllowedRms;
  cfg.soundBeatTrack  = soundBeatTrack / rmsGrid * maxAllowedRms;

  %% load audio tracks

  % prepare cell for audio tracks
  cfg.soundTracks = cell(1, length(cfg.patterns));
  % prepare cell for metronome tracks
  cfg.soundTrackBeat = cell(1, length(cfg.patterns));

  for tracki = 1:nTracks
    % load audio track
    sTrack = audioread(fullfile('.', 'stimuli', cfg.patterns{cfg.isTrackIdx(tracki)}));
    % make it mono and transpose to row vector
    sTrack = mean(sTrack, 2)';
    % set RMS
    sTrack = sTrack / rms(sTrack) * maxAllowedRms;
    % assign it to config
    cfg.soundTracks{cfg.isTrackIdx(tracki)} = sTrack;

    % generate beat sequence as audio track! (we don't use grid here)
    gridIOI     = cfg.gridIOI(cfg.isTrackIdx(tracki));
    beatPeriod  = cfg.cuePeriodGrid(cfg.isTrackIdx(tracki)) * gridIOI;
    seqDur      = round(length(sTrack) / cfg.fs);
    nBeatsInSeq = floor(seqDur / beatPeriod);
    beatTimes   = beatPeriod * [0:nBeatsInSeq - 1];
    seqBeat     = zeros(1, length(sTrack));

    sBeatIdx = round(beatTimes * cfg.fs);
    for i = 1:length(sBeatIdx)
      seqBeat(sBeatIdx(i) + 1:sBeatIdx(i) + length(cfg.soundBeatTrack)) = cfg.soundBeatTrack;
    end

    cfg.soundTrackBeat{cfg.isTrackIdx(tracki)} = seqBeat;

    % get window duration in seconds
    cfg.winDur(cfg.isTrackIdx(tracki)) = cfg.nCyclesPerWin(cfg.isTrackIdx(tracki)) * beatPeriod;

  end

  %% Grahn(2007) patterns

  cfg.pattern.fs = cfg.fs;

  % read from txt files
  grahnPatSimple = loadIOIRatiosFromTxt(fullfile('stimuli', 'Grahn2007_simple.txt'));
  grahnPatComplex = loadIOIRatiosFromTxt(fullfile('stimuli', 'Grahn2007_complex.txt'));

  % get different metrics of the patterns
  cfg.pattern.patternSimple = getPatternInfo(grahnPatSimple, 'simple', cfg);
  cfg.pattern.patternComplex = getPatternInfo(grahnPatComplex, 'complex', cfg);

  % the grid interval can vary across steps or segments (gridIOI selected
  % randomly from a set of possible values for each new step or segment)
  cfg.pattern.gridIOI = 0.190;

  % Define envelope shape of the individual sound event.
  % All parameters are defined in seconds.

  % total sound duration _/```\_
  cfg.pattern.eventDur             = 0.190; % s
  % onset ramp duration  _/
  cfg.pattern.eventRampon          = 0.010; % s
  % offset ramp duration       \_
  cfg.pattern.eventRampoff         = 0.020; % s

  % Make sure the total ramp durations are not longer than tone duration.
  if (cfg.pattern.eventRampon + cfg.pattern.eventRampoff) > cfg.pattern.eventDur
    error(sprintf('The summed duration of onset+offset ramps (%g ms) is longer than requensted tone duration (%g ms).', ...
                  (cfg.pattern.eventRampon + cfg.pattern.eventRampoff) * 1e3, ...
                  cfg.pattern.eventDur * 1e3));
  end
  % Make sure the tone duration is not longer than smallest gridIOI.
  if cfg.pattern.eventDur >  cfg.pattern.gridIOI
    error(sprintf('Requested tone duration (%g ms) is longer than shortest gridIOI (%g ms).', ...
                  cfg.pattern.eventDur * 1e3, ...
                  cfg.pattern.gridIOI * 1e3));
  end

  % construct pattern (smallest item in sequence)
  cfg.pattern.nGridPoints = 12; % length(pat_complex(1).pattern)

  cfg.pattern.interPatternInterval = cfg.pattern.nGridPoints * cfg.pattern.gridIOI;

  % construct pitch features of the stimulus
  % the pitch (F0) of the tones making up the patterns can vary
  % (it can be selected randomly from a set of possible values)
  cfg.pattern.minF0   = 350; % minimum possible F0
  cfg.pattern.maxF0   = 900; % maximum possible F0
  cfg.pattern.nF0       = 5; % number of unique F0-values between the limits
  cfg.pattern.F0s       = logspace(log10(cfg.pattern.minF0), log10(cfg.pattern.maxF0), cfg.pattern.nF0);

  % use the requested gain of each tone to adjust the base amplitude
  cfg.pattern.F0sAmpGain = equalizePureTones(cfg.pattern.F0s, [], []);
  cfg.pattern.F0sAmps = 1 / sqrt(2) * maxAllowedRms * cfg.pattern.F0sAmpGain;

  %% generate example stimulus/sequence for only volume setting

  volTestSound = makeStimTrain(cfg, 1, 1, 0, 0);
  % make sequence for 2 channels
  cfg.volumeSettingSound = repmat(volTestSound.s, 2, 1);

  %% Instructions
  % !!! NOTE: use UTF-8 encoding, otherwise there will be problem with quotation marks

  loadPathInstr = fullfile('subfun', 'instr', 'tapTrainer');

  % -----------------------------------
  % general task instructions and intro
  % -----------------------------------
  % These need to be saved in separate files, named: 'instrTrainingIntro#'
  % The text in each file will be succesively displayed on the screen at the
  % begining of the experiment. Every time, the script will wait for
  % keypress.

  dirInstr = dir(fullfile(loadPathInstr, 'instrTrainingIntro*'));
  cfg.taskInstruction = cell(1, length(dirInstr));
  for i = 1:length(dirInstr)
    instrFid = fopen(fullfile(loadPathInstr, dirInstr(i).name), 'r', 'n', 'UTF-8');
    while ~feof(instrFid)
      cfg.taskInstruction{i} = [cfg.taskInstruction{i}, fgets(instrFid)];
    end
    fclose(instrFid);
  end

  % ------------------------------
  % sequence-specific instructions
  % ------------------------------
  % for each pattern/track (each sequence), there can be specific
  % instructions, which explains some important concepts that should be learned by the
  % participants.
  cfg.beforeSeqInstruction = cell(1, length(cfg.patterns));
  cfg.afterSeqInstruction = cell(1, length(cfg.patterns));

  % also, for each sequence, we can specify unique instructions to be
  % displayed during sound presentation
  cfg.duringSeqInstruction_part1 = cell(1, length(cfg.patterns));
  cfg.duringSeqInstruction_part2 = cell(1, length(cfg.patterns));

  % get filenames for the instruction textfiles of different categories (sort
  % them by filenames using natural sort)
  dirInstrBefore = dir(fullfile(loadPathInstr, sprintf('instrTrainingBeforeSeq*')));
  [~, sortNatIdx] = sortNatural({dirInstrBefore.name});
  dirInstrBefore = dirInstrBefore(sortNatIdx);

  dirInstrDuring1 = dir(fullfile(loadPathInstr, sprintf('instrTrainingDuringSeq*_part1')));
  [~, sortNatIdx] = sortNatural({dirInstrDuring1.name});
  dirInstrDuring1 = dirInstrDuring1(sortNatIdx);

  dirInstrDuring2 = dir(fullfile(loadPathInstr, sprintf('instrTrainingDuringSeq*_part2')));
  [~, sortNatIdx] = sortNatural({dirInstrDuring2.name});
  dirInstrDuring2 = dirInstrDuring2(sortNatIdx);

  dirInstrAfter = dir(fullfile(loadPathInstr, sprintf('instrTrainingAfterSeq*')));
  [~, sortNatIdx] = sortNatural({dirInstrAfter.name});
  dirInstrAfter = dirInstrAfter(sortNatIdx);

  % loop to load the instructions for each pattern/track
  for pati = 1:length(cfg.patterns)

    % --- BEFORE each sequence ---

    % find relevant text files with instructions
    instrIdx = find(~cellfun(@isempty, regexp({dirInstrBefore.name}, ...
                                              sprintf('instrTrainingBeforeSeq%d(($)|(-\\d*))', pati))));
    beforeSeqInstructionTmp = cell(1, length(instrIdx));
    if ~isempty(instrIdx)
      % if you can find text file(s), load it/them
      for instri = 1:length(instrIdx)
        instrFid = fopen(fullfile(loadPathInstr, dirInstrBefore(instrIdx(instri)).name), 'r', 'n', 'UTF-8');
        tmptxt = [];
        while ~feof(instrFid)
          tmptxt = [tmptxt, fgets(instrFid)];
        end
        fclose(instrFid);
        % assign each instruction text file to temporary cell
        beforeSeqInstructionTmp{instri} = tmptxt;
      end
    else
      % if not, just write empty text
      cfg.beforeSeqInstruction{pati} = '';
      warning(sprintf('no instructions found before pattern %d', pati));
    end
    % assignn the resulting cell to the expParam structure
    cfg.beforeSeqInstruction{pati} = beforeSeqInstructionTmp;

    % --- AFTER each sequence ---

    % find relevant text files with instructions
    instrIdx = find(~cellfun(@isempty, regexp({dirInstrAfter.name}, ...
                                              sprintf('instrTrainingAfterSeq%d(($)|(-\\d*))', pati))));
    afterSeqInstructionTmp = cell(1, length(instrIdx));
    if ~isempty(instrIdx)
      % if you can find text file(s), load it/them
      for instri = 1:length(instrIdx)
        instrFid = fopen(fullfile(loadPathInstr, dirInstrAfter(instrIdx(instri)).name), 'r', 'n', 'UTF-8');
        tmptxt = [];
        while ~feof(instrFid)
          tmptxt = [tmptxt, fgets(instrFid)];
        end
        fclose(instrFid);
        % assign each instruction text file to temporary cell
        afterSeqInstructionTmp{instri} = tmptxt;
      end
    else
      % if not, just write empty text
      cfg.afterSeqInstruction{pati} = '';
      warning(sprintf('no instructions found after pattern %d', pati));
    end
    % assignn the resulting cell to the expParam structure
    cfg.afterSeqInstruction{pati} = afterSeqInstructionTmp;

    % --- DURING each sequence (part 1) ---

    if exist(fullfile(loadPathInstr, sprintf('instrTrainingDuringSeq%d_part1', pati)))
      % if you can find a text file, load it
      instrFid = fopen(fullfile(loadPathInstr, sprintf('instrTrainingDuringSeq%d_part1', pati)), 'r', 'n', 'UTF-8');
      tmptxt = [];
      while ~feof(instrFid)
        tmptxt = [tmptxt, fgets(instrFid)];
      end
      fclose(instrFid);
      cfg.duringSeqInstruction_part1{pati} = tmptxt;
    else
      % if not, just write empty text
      cfg.duringSeqInstruction_part1{pati} = '';
      warning(sprintf('no instructions found during pattern %d (part 1)', pati));
    end

    % --- DURING each sequence (part 2) ---

    if exist(fullfile(loadPathInstr, sprintf('instrTrainingDuringSeq%d_part2', pati)))
      % if you can find a text file, load it
      instrFid = fopen(fullfile(loadPathInstr, sprintf('instrTrainingDuringSeq%d_part2', pati)), 'r', 'n', 'UTF-8');
      tmptxt = [];
      while ~feof(instrFid)
        tmptxt = [tmptxt, fgets(instrFid)];
      end
      fclose(instrFid);
      cfg.duringSeqInstruction_part2{pati} = tmptxt;
    else
      % if not, just write empty text
      cfg.duringSeqInstruction_part2{pati} = '';
      warning(sprintf('no instructions found during pattern %d (part 2)', pati));
    end

  end

end
