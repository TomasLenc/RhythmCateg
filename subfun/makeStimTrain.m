function [seq] = makeStimTrain(cfg, currPatterni, cueDBleveli, currWini, soundIdx)
  %
  % Input
  % -----
  % cfg : struct
  %     configuration structure
  % currPatterni: int
  %     index of the currently played pattern (or track)
  % cueDBleveli: int
  %     index of the currently used cue dB level
  % currWini: int
  %     index of the current analysis window (for the current sequence)
  %     if currWini <= cfg.nWinNoCue, the cue will be silenced
  % audioIdx: int
  %     optional, current index position direclty after the audio that's been already
  %     pushed to buffer and played
  %
  % Output
  % ------
  % seq : struct
  %     structure with the generated audio and information
  %
  % ============================================================================================

  seq                 = [];
  seq.fs              = cfg.fs;

  if currWini <= cfg.nWinNoCue(currPatterni)
    seq.cueDB = -Inf;
  else
    seq.cueDB = cfg.cueDB(cueDBleveli);
  end

  %% this is an audio track that has been loaded
  if ismember(currPatterni, cfg.isTrackIdx)

    % check that we have the audioIdx in the input
    if nargin < 4
      error('you need to supply audioIdx when calling makeStimTrain with audio track!');
    end

    seq.idxEnd = min(soundIdx + round(cfg.winDur(currPatterni) * cfg.fs), ...
                     length(cfg.soundTracks{currPatterni}));

    % scale beat track rms
    soundBeat = cfg.soundTrackBeat{currPatterni}(soundIdx + 1:seq.idxEnd);
    rmsBeat = rms(soundBeat);
    soundBeat = soundBeat / rmsBeat * rmsBeat * 10^(seq.cueDB / 20);

    % add sound track + beat track
    seq.s = cfg.soundTracks{currPatterni}(soundIdx + 1:seq.idxEnd) + soundBeat;

    % check if we are at the end of audio track, if yes, send a flag
    if seq.idxEnd == length(cfg.soundTracks{currPatterni})
      seq.AUDIO_END = 1;
    end

    %% this is a sequence made of Grahn2007 patterns
  elseif strcmpi(class(cfg.patterns{currPatterni}), 'char')

    nSamplesPat     = round(cfg.pattern.interPatternInterval * cfg.fs);
    seq.dur         = cfg.pattern.interPatternInterval * cfg.nCyclesPerWin(currPatterni);
    seq.nSamples    = round(seq.dur * cfg.fs);
    seq.nCycles     = cfg.nCyclesPerWin(currPatterni);
    seq.idxEnd      = 0; % dummy, only used for audiotracks
    seq.patternIDs  = cell(1, seq.nCycles);
    seq.pattern     = nan(1, seq.nCycles * cfg.pattern.nGridPoints);
    seq.cue         = repmat([1, zeros(1, cfg.cuePeriodGrid(currPatterni) - 1)], ...
                             1, floor(length(seq.pattern) / cfg.cuePeriodGrid(currPatterni)));

    seqPattern = nan(1, seq.nSamples);

    % get the F0 for each pattern (without direct repetition)
    F0s2use             = zeros(1, cfg.nCyclesPerWin(currPatterni));
    amps2use            = zeros(1, cfg.nCyclesPerWin(currPatterni));
    availF0idx          = 1:length(cfg.pattern.F0s);
    currF0idx           = Inf;
    for pati = 1:cfg.nCyclesPerWin(currPatterni)
      availF0idx      = 1:length(cfg.pattern.F0s);
      availF0idx(availF0idx == currF0idx) = [];
      currF0idx       = randsample(availF0idx, 1);
      F0s2use(pati)   = cfg.pattern.F0s(currF0idx);
      amps2use(pati)  = cfg.pattern.F0sAmps(currF0idx);

    end
    seq.F0s = F0s2use;

    % get pattern IDs (without direct repetition)
    if ~isempty(regexpi(cfg.patterns{currPatterni}, 'simple'))

      availPatIdx = 1:length(cfg.pattern.patternSimple);
      currPatIdx = Inf;
      idx = 0;
      for pati = 1:cfg.nCyclesPerWin(currPatterni)
        % indices of patterns are available to choose from (to avoid
        % succesive repetition)
        availPatIdx = 1:length(cfg.pattern.patternSimple);
        availPatIdx(availPatIdx == currPatIdx) = [];
        % choose random pattern
        currPatIdx = randsample(availPatIdx, 1);
        currPattern = cfg.pattern.patternSimple(currPatIdx).pattern;
        seq.patternIDs{pati} = cfg.pattern.patternSimple(currPatIdx).ID;
        seq.pattern((pati - 1) * cfg.pattern.nGridPoints + 1:pati * cfg.pattern.nGridPoints) = currPattern;
        % generate audio
        seqPattern(idx + 1:idx + nSamplesPat) = amps2use(pati) * makeStimMainExp(currPattern, ...
                                                                                 cfg, ...
                                                                                 cfg.pattern.gridIOI, ...
                                                                                 F0s2use(pati));
        % update audio index
        idx = idx + nSamplesPat;
      end

    elseif ~isempty(regexpi(cfg.patterns{currPatterni}, 'complex'))
      availPatIdx = 1:length(cfg.pattern.patternComplex);
      currPatIdx = Inf;
      idx = 0;
      for pati = 1:cfg.nCyclesPerWin(currPatterni)
        % indices of patterns are available to choose from (to avoid
        % succesive repetition)
        availPatIdx = 1:length(cfg.pattern.patternComplex);
        availPatIdx(availPatIdx == currPatIdx) = [];
        % choose random pattern
        currPatIdx = randsample(availPatIdx, 1);
        currPattern = cfg.pattern.patternComplex(currPatIdx).pattern;
        seq.patternIDs{pati} = cfg.pattern.patternComplex(currPatIdx).ID;
        seq.pattern((pati - 1) * cfg.pattern.nGridPoints + 1:pati * cfg.pattern.nGridPoints) = currPattern;
        % generate audio
        seqPattern(idx + 1:idx + nSamplesPat) = amps2use(pati) * makeStimMainExp(currPattern, ...
                                                                                 cfg, ...
                                                                                 cfg.pattern.gridIOI, ...
                                                                                 F0s2use(pati));
        % update audio index
        idx = idx + nSamplesPat;
      end

    end

    % set the requested cue dB
    rmsBeat = rms(cfg.soundBeat);
    cfg.soundBeat = cfg.soundBeat / rmsBeat * rmsBeat * 10^(seq.cueDB / 20);

    rmsGrid = rms(cfg.soundGrid);
    cfg.soundGrid = cfg.soundGrid / rmsGrid * rmsGrid * 10^(seq.cueDB / 20);

    % further attenuate grid sound (fixed attenuation)
    cfg.soundGrid = cfg.soundGrid * 1 / 4;

    % generate metrononme sequence
    seqBeat = zeros(1, seq.nSamples);
    sBeatIdx = round((find(seq.cue) - 1) * cfg.gridIOI(currPatterni) * seq.fs);
    for i = 1:length(sBeatIdx)
      seqBeat(sBeatIdx(i) + 1:sBeatIdx(i) + length(cfg.soundBeat)) = cfg.soundBeat;
    end

    % generate grid sequence
    seqGrid = zeros(1, seq.nSamples);
    sGridIdx = round((find(ones(1, length(seq.pattern))) - 1) * cfg.gridIOI(currPatterni) * seq.fs);
    for i = 1:length(sGridIdx)
      seqGrid(sGridIdx(i) + 1:sGridIdx(i) + length(cfg.soundGrid)) = cfg.soundGrid;
    end

    % add them together
    seq.s = seqPattern + seqBeat + seqGrid;

    %% this is a regular pattern that needs to be synthesized
  elseif strcmpi(class(cfg.patterns{currPatterni}), 'double')

    % add parameters to output structure
    seq.pattern         = cfg.patterns{currPatterni};
    seq.cue             = repmat([1, zeros(1, cfg.cuePeriodGrid(currPatterni) - 1)], ...
                                 1, floor(length(seq.pattern) / cfg.cuePeriodGrid(currPatterni)));
    seq.nCycles         = cfg.nCyclesPerWin(currPatterni);
    seq.dur             = length(seq.pattern) * seq.nCycles * cfg.gridIOI(currPatterni);
    seq.nSamples        = round(seq.dur * seq.fs);

    seq.idxEnd          = 0; % dummy, only used for audiotracks

    % set the requested cue dB
    rmsBeat = rms(cfg.soundBeat);
    cfg.soundBeat = cfg.soundBeat / rmsBeat * rmsBeat * 10^(seq.cueDB / 20);

    rmsGrid = rms(cfg.soundGrid);
    cfg.soundGrid = cfg.soundGrid / rmsGrid * rmsGrid * 10^(seq.cueDB / 20);

    % further attenuate grid sound (fixed attenuation)
    cfg.soundGrid = cfg.soundGrid * 1 / 4;

    % generate pattern sequence
    seqPattern = zeros(1, seq.nSamples);
    sPatIdx = round((find(repmat(seq.pattern, 1, seq.nCycles)) - 1) * cfg.gridIOI(currPatterni) * seq.fs);
    for i = 1:length(sPatIdx)
      seqPattern(sPatIdx(i) + 1:sPatIdx(i) + length(cfg.soundPattern)) = cfg.soundPattern;
    end

    % generate metrononme sequence
    seqBeat = zeros(1, seq.nSamples);
    sBeatIdx = round((find(repmat(seq.cue, 1, seq.nCycles)) - 1) * cfg.gridIOI(currPatterni) * seq.fs);
    for i = 1:length(sBeatIdx)
      seqBeat(sBeatIdx(i) + 1:sBeatIdx(i) + length(cfg.soundBeat)) = cfg.soundBeat;
    end

    % generate grid sequence
    seqGrid = zeros(1, seq.nSamples);
    sGridIdx = round((find(ones(1, seq.nCycles * length(seq.pattern))) - 1) * cfg.gridIOI(currPatterni) * seq.fs);
    for i = 1:length(sGridIdx)
      seqGrid(sGridIdx(i) + 1:sGridIdx(i) + length(cfg.soundGrid)) = cfg.soundGrid;
    end

    % add them together
    seq.s = seqPattern + seqBeat + seqGrid;

  end

  % check sound amplitude for clipping
  if max(abs(seq.s)) > 1
    warning('sound amplitude larger than 1...normalizing');
    seq.s = seq.s ./ max(abs(seq.s));
  end
