function target = collectAndSaveEvents(cfg, logFile, currentSequence, SeqNum, audioOnset, task)

  % write into logfile
  currentSequence(1).isStim = false; 
  currentSequence(1).fileID = logFile(1).fileID;
  currentSequence(1).extraColumns = logFile(1).extraColumns;

  % adding columns in currSeq for BIDS format
  for iPattern = 1:numel(currentSequence)

    % adding compulsory BIDS structures
    currentSequence(iPattern, 1).trial_type = currentSequence(iPattern, 1).segmentCateg;        
    
    currentSequence(iPattern, 1).duration = cfg.pattern.interPatternInterval;
    
    % adding other interest
    currentSequence(iPattern, 1).sequenceNum = SeqNum;
    target(iPattern, 1) = currentSequence(iPattern, 1).isTask;

  end

  % add dummy event that marks the PTB onset time of the sequence
  % this can be use to troubleshoot trigger order in biosemi recording to
  % identify individual trials
  nEvents = size(currentSequence,1); 
  currentSequence(nEvents+1, 1).trial_type = 'sequence_onset'; 
  currentSequence(nEvents+1, 1).onset = audioOnset; 
  currentSequence(nEvents+1, 1).duration = cfg.pattern.SequenceDur; 
  currentSequence(nEvents+1, 1).triggerValue = cfg.beh.trigTaskMapping(task); 
  
  saveEventsFile('save', cfg, currentSequence); 

end
