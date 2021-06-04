function seq = makeSequence(cfg,seqi,varargin)
% This function constructs a stimulus sequence.
% by using makeStimMainExp.m script

% it's also depending on getAllSeqDesign.m output with the given 
% run number==seqi (fmri) or expParam.numSequences == seqi (behav)

% ------
% INPUT
% ------
%     cfg:          structure with confuguration info
%     seqi:         sequence iterator index (which sequence in the
%                   experiment this is?)
% 
% ------
% OUTPUT
% ------
%     seq:          structure with stimulus sequence and info about it
%
%


if cfg.debug.do
    seqi = 1;
end
    
%% allocate variables to log

% main output structure (we'll put everything else into it)
seq = struct(); 

% vector of F0 (pitch) values for each pattern in the sequence
seq.F0 = zeros(1, cfg.pattern.nPatternPerSegment * ...
    cfg.pattern.nSegmPerStep * cfg.pattern.nStepsPerSequence); 

% vector of gridIOI values for each pattern in the sequence
seq.gridIOI = zeros(1, cfg.pattern.nPatternPerSegment * ...
    cfg.pattern.nSegmPerStep * cfg.pattern.nStepsPerSequence); 

% segment-category (A or B) for each pattern in the sequence 
seq.segmentCateg = cell(1, cfg.pattern.nPatternPerSegment * ...
    cfg.pattern.nSegmPerStep * cfg.pattern.nStepsPerSequence); 

% onset time of each pattern
% i'm still conflicting how to make this onset both BIDS compatible &
% explicit that it's PATTERN ONSET WE ARE RECORDING
seq.onset = nan(1, cfg.pattern.nPatternPerSegment * ...
    cfg.pattern.nSegmPerStep * cfg.pattern.nStepsPerSequence); 

% put together all the patterns from both categories, we will pick from
% this using the unique ID of each pattern (we know which IDs we want from
% the output of getAllSeq function. 
patterns2choose = [cfg.pattern.patternA,cfg.pattern.patternB]; 


% each pattern will have its own ID (integer number; patterns with the same
% ID number from category A vs. B will have the same inter-onset intervals,
% just rearranged in time)

% TL: I know the same thing is saved twice (i.e. the category information is, saved 
% in the "category" colum, but is also in the rhythm ID string. But I'd perhaps keep this
% redundancy for the sake of safety :) Because this way I can directly check if there is
% any discrepancy between what the script thinks is category A and the
% actual patterns that are used for that category. So this info is primarily for visual
% checking that the script is doing what it's meant to. If I need the
% rhythm ID number at any point I can always get it out of the string with regexp. 

seq.patternID = cell(1, cfg.pattern.nPatternPerSegment * ...
    cfg.pattern.nSegmPerStep * cfg.pattern.nStepsPerSequence); 


% audio waveform of the sequence
seq.outAudio = zeros(1,round(cfg.pattern.SequenceDur*cfg.fs)); 

% % audio envelop  of the sequence
% seq.outEnvelop = zeros(1,round(cfg.pattern.SequenceDur*cfg.fs)); 

% carries pitches length 
numPitch = cfg.pattern.nF0;

%% initialize counters 

% currently chosen F0 index (indexing value in cfg.pattern.F0s, initialize to 1)
currF0idx = 1; 

% currently chosen gridIOI index (indexing value in cfg.pattern.gridIOIs, initialize to 1)
currGridIOIidx = 1; 

% current time point in the sequence as we go through the loops
currTimePoint = 0; 

% pitch counter 
cPitch = 1;

% loop over steps
for stepi=1:cfg.pattern.nStepsPerSequence
    
    % take the timestamp for logging the current step time
    stepOnset  = currTimePoint;
    
    % loop over segments in 1 step
    for segmi=1:cfg.pattern.nSegmPerStep
        
        % take the  timestamp for logging the current segment time
        segmentOnset  = currTimePoint;
        
        % find which patterns are part of the current segment 
        nPatBeforeThis = (stepi-1)*cfg.pattern.nSegmPerStep*cfg.pattern.nPatternPerSegment ...
                         + (segmi-1)*cfg.pattern.nPatternPerSegment; 
                     
        patIdx = [nPatBeforeThis+1 : nPatBeforeThis+cfg.pattern.nPatternPerSegment]; 
        
        % loop over patterns in 1 segment                
        for pati=1:cfg.pattern.nPatternPerSegment
                        
            % find the pattern ID from the seqDesignFullExp (output of
            % getAllSeq function)
            currPatID = cfg.pattern.seqDesignFullExp{seqi,stepi,segmi,pati}; 
            currPatIdx = find(strcmp(currPatID,{patterns2choose.ID}));

            % find if this pattern categ, simple or complex
            currPatternCateg = regexp(patterns2choose(currPatIdx(1)).ID, '\D*(?=\d.*)', 'match'); 
            currPatternCateg = currPatternCateg{1}; 
            currCategLabel = currPatternCateg;
            
            %find segment info, if it's A/B
            currSegmentLabel = cfg.pattern.seqDesignSegment{seqi,stepi,segmi,pati}; 
            
            % --------------------------------------------------
            % ----- determine if gridIOI needs to be changed ---
            % --------------------------------------------------
            CHANGE_IOI = 0; 
            
            % gridIOI change requested every segment (and this is the first
            % cycle in the segment)
            if cfg.pattern.changeGridIOISegm && pati == 1
                CHANGE_IOI = 1; 
                
            % gridIOI change requested every category (and this is the first
            % cycle in the segment, and a category just changed from A->B
            % or B->A)
            elseif cfg.pattern.changeGridIOICategory && pati == 1 && ...
                    ( segmi == 1 || segmi == cfg.pattern.nSegmentA+1 )
                CHANGE_IOI = 1; 
                
            % gridIOI change requested every step 
            elseif cfg.pattern.changeGridIOIStep && segmi==1
                CHANGE_IOI = 1; 
                
            end
                        
            % if change of gridIOI requested, randomly choose a new gridIOI
            % do this only if there is more than 1 gridIOI to choose from
            if CHANGE_IOI && length(cfg.pattern.gridIOIs)>1
                % get gridIOI to choose from 
                gridIOI2ChooseIdx = 1:length(cfg.pattern.gridIOIs); 
                % remove gridIOI used in the previous iteration (to prevent
                % repetition in the sequence) 
                gridIOI2ChooseIdx(gridIOI2ChooseIdx==currGridIOIidx) = []; 
                % randomly select new IOI
                currGridIOIidx = randsample(gridIOI2ChooseIdx); 
            end
            
            currGridIOI = cfg.pattern.gridIOIs(currGridIOIidx); 
            
            
            % --------------------------------------------------
            % ----- determine if pitch needs to be changed -----
            % --------------------------------------------------
            CHANGE_PITCH = 0; 
            
            % pitch change requested in every pattern cycle
            if cfg.pattern.changePitchCycle
                CHANGE_PITCH = 1; 
                
            % pitch change requested in every segment (and this is the first
            % pattern cycle in the segment)
            elseif cfg.pattern.changePitchSegm && pati==1
                CHANGE_PITCH = 1; 
                
            % pitch change requested every category (and this is the first
            % pattern cycle in the segment, and a category just changed from A->B
            % or B->A)
            elseif cfg.pattern.changePitchCategory && pati==1 && ...
                    ( segmi==1 || segmi==cfg.pattern.nSegmentA+1 )
                CHANGE_PITCH = 1; 
                
            % pitch change requested only in every step 
            elseif cfg.pattern.changePitchStep && segmi==1
                CHANGE_PITCH = 1; 
                
            end
            
            if isfield(cfg.pattern,'fixedPitchCategB')
                
                % only categB is with fixed pitch
                if cfg.pattern.fixedPitchCategB && strcmpi(currSegmentLabel,'B')
                    %assign to the 5th pitch to all categB patterns
                    currF0 = cfg.pattern.differF0;
                    currAmp = cfg.pattern.F0sAmps(end);
                    currF0idx = cfg.pattern.nF0 + 1; % 5th pitch !!!!!!!!
                    
                elseif cfg.pattern.fixedPitchCategB && strcmpi(currSegmentLabel,'A')
                    
                    currF0idx = squeeze(cfg.pattern.seqDesignToneF0(seqi,stepi,segmi,pati,:)); 

                    currF0idxMask = currF0idx > 0; 
                    
                    % assign the index to current F0 index
                    currF0idx(currF0idx==0) = []; 
                    
                    %assign the randomly chosen ones to current pitch
                    currF0 = zeros(size(currF0idxMask)); 
                    currF0(currF0idxMask) = cfg.pattern.F0s(currF0idx);
                    
                    currAmp = zeros(size(currF0idxMask)); 
                    currAmp(currF0idxMask) = cfg.pattern.F0sAmps(currF0idx);
                    
                    currF0idx = squeeze(cfg.pattern.seqDesignToneF0(seqi,stepi,segmi,pati,:)); 

                end
            end 

            if ~isfield(cfg.pattern,'fixedPitchCategB') || ~ cfg.pattern.fixedPitchCategB
                % if fixedPitchCategB is not defined
                if CHANGE_PITCH && length(cfg.pattern.F0s)>1
                    
                    % counterbalance the pitches across patterns
                    if mod(pati, numPitch) == 1
                        
                        %reset the counter
                        cPitch = 1;
                        
                        % shuffle the F0 array & get one F0
                        arrayPitchIdx = Shuffle(1:numPitch);
                        pitch2ChooseIdx = arrayPitchIdx(cPitch);
                        
                        % prevent repetition of pitch in sequential
                        % patterns
                        while pitch2ChooseIdx == currF0idx
                            
                            % shuffle the F0 array & get one F0
                            arrayPitchIdx = Shuffle(1:numPitch);
                            pitch2ChooseIdx = arrayPitchIdx(cPitch);
                            
                        end
                        
                    else
                        % increase pitch counter
                        cPitch = cPitch+1;
                        % get the following F0
                        pitch2ChooseIdx = arrayPitchIdx(cPitch);
                    end
                    
                    % assign the index to current F0 index
                    currF0idx = pitch2ChooseIdx; 
                    
                    %assign the randomly chosen ones to current pitch
                    currF0 = cfg.pattern.F0s(currF0idx);
                    currAmp = cfg.pattern.F0sAmps(currF0idx);

                else
                    currF0 = cfg.pattern.F0s(currF0idx);
                    currAmp = cfg.pattern.F0sAmps(currF0idx);
                    
                end
            end
        
            % -----------------------------------
            % -----    save pattern data   -----
            % -----------------------------------
            
            % get the pattern
            currPattern = patterns2choose(currPatIdx).pattern;
                                                            
            % First, check for fmri task exists
            if isfield(cfg.pattern,'taskIdxMatrix')
                cfg.isTask.Idx = cfg.pattern.taskIdxMatrix(seqi,stepi,segmi,pati);
                
                % the current F0s index is used for finding the
                % taskSound
                cfg.isTask.F0Idx = currF0idx;
            else 
                cfg.isTask.Idx = 0;
            end
            
            seq(patIdx(pati),1).patternID   = currPatID;
            seq(patIdx(pati),1).segmentCateg   = currCategLabel;
            seq(patIdx(pati),1).segmentLabel   = currSegmentLabel;
            seq(patIdx(pati),1).segmentNum  = segmi;
            seq(patIdx(pati),1).segmentOnset = segmentOnset;
            seq(patIdx(pati),1).stepNum     = stepi;
            seq(patIdx(pati),1).stepOnset   = stepOnset;
            seq(patIdx(pati),1).taskF0idx   = currF0idx;
            seq(patIdx(pati),1).isTask      = cfg.isTask.Idx;
            
            seq(patIdx(pati),1).pattern     = currPattern; 
            seq(patIdx(pati),1).nSounds     = length(find(currPattern)); 
            seq(patIdx(pati),1).F0          = currF0(1);
            seq(patIdx(pati),1).gridIOI     = currGridIOI;
            seq(patIdx(pati),1).patternAmp  = currAmp(1);

            % get pattern info e.g. PE and LHL
            seq(patIdx(pati),1).PE4        = patterns2choose(currPatIdx).PE4;
            seq(patIdx(pati),1).minPE4     = patterns2choose(currPatIdx).minPE4;
            seq(patIdx(pati),1).rangePE4   = patterns2choose(currPatIdx).rangePE4;
            seq(patIdx(pati),1).LHL24      = patterns2choose(currPatIdx).LHL24;
            seq(patIdx(pati),1).minLHL24   = patterns2choose(currPatIdx).minLHL24;
            seq(patIdx(pati),1).rangeLHL24 = patterns2choose(currPatIdx).rangeLHL24;
                        
        end % pattern loop 
        
        

        % ------------------------------------------------------------------
        % ----- get inter-onset interval representation for each pattern ---
        % ------------------------------------------------------------------
        
        currSegmPatterns = [seq(patIdx).pattern]; 
%         nSoundsPerPattern = [seq(patIdx).nSounds]; 

        % convert grid representation inter-onset interval representation
        currSegmIOIs = diff(find(currSegmPatterns)); 

        % find what different inter-onset intervals are in the sequence (1,2,3,4)
        IOIClasses = unique(currSegmIOIs); 

        
        % boool flag, do the nonmetric scambling 
        DO_NONMETRIC = 0;
        if isfield(cfg.pattern,'doNonMetric') && strcmp(currCategLabel,'complex')
                DO_NONMETRIC = 1;
                ioiScrambleRatio = cfg.pattern.ioiScrambleRatio;
        end
                
        if DO_NONMETRIC

            % allocate new nonmetric IOI representation 
            currSegmIOIsNonmetric = currSegmIOIs; 

            % loop over IOI categories (i.e. 1,2,3,4 in Grahn's rhythms)
            for IOIClass=IOIClasses

                % ignore intervals with gridIOI value, because we can't shorten them (the
                % successive tones would overlap because we've set tone duration to be
                % pretty much gridIOI)
                if IOIClass==1
                    continue
                end

                % which IOIs in the pattern have the current value 
                idxIOIClass = find(currSegmIOIs==IOIClass); 

                if mod(length(idxIOIClass),2)==0
                    % if even, shorten half and lenghten half 
                    idxToLenghten = randsample(idxIOIClass, length(idxIOIClass)/2); 
                    idxToShorten = setdiff(idxIOIClass,idxToLenghten); 
                else
                    % if odd, keep one intact and do the same 
                    idxToKeep = randsample(idxIOIClass,1); 
                    idxIOIClass(idxIOIClass==idxToKeep) = []; 
                    % use round, incase length == 1
                    idxToLenghten = randsample(idxIOIClass, round(length(idxIOIClass)/2)); 
                    idxToShorten = setdiff(idxIOIClass,idxToLenghten); 
                end
                currSegmIOIsNonmetric(idxToLenghten) = currSegmIOIs(idxToLenghten) * 1+ioiScrambleRatio; 
                currSegmIOIsNonmetric(idxToShorten) = currSegmIOIs(idxToShorten) * 1-ioiScrambleRatio; 
            end

            currSegmIOIs = currSegmIOIsNonmetric; 
            
            % re-write the segment category if it is nonmetric                               
            [seq(patIdx,1).segmentCateg]   = deal('nonmetric');
        end
        



        % -----------------------------------------------------------------------
        % ----- synthesize audio for the whole segment and assign to sequence ---
        % -----------------------------------------------------------------------
        c = 0; 

        segmAudio = zeros(1, round(cfg.pattern.interSegmInterval*cfg.fs)); 

        for pati=1:length(patIdx)

            % get sound event IOIs and convert to seconds 
            if pati<length(patIdx)
                currPatIOIs = currSegmIOIs(c+1:c+seq(patIdx(pati),1).nSounds) * seq(patIdx(pati),1).gridIOI; 
            else
                currPatIOIs = currSegmIOIs(c+1:c+seq(patIdx(pati),1).nSounds-1) * seq(patIdx(pati),1).gridIOI;        
            end

            % from IOIs, get onset time of each sound event wrt current segment onset 
            if isfield(seq,'IOIs')
                patternOnsetTimeWrtSegm = sum( [seq(patIdx(patIdx<patIdx(pati)),1).IOIs] ); 
            else
                patternOnsetTimeWrtSegm = 0; 
            end
            
            %add into structure
            seq(patIdx(pati),1).soundOnsetTimesWrtSegm = patternOnsetTimeWrtSegm + ...
                                                 cumsum([0, currPatIOIs(1:end-1)]);  


            % if this is the last pattern in the segment, we need to add one IOI for completeness                               
            if pati==length(patIdx)
                % time from last sound event in the current segment to next segment onset
                ioiToNextSegm = cfg.pattern.interSegmInterval - ...
                                seq(patIdx(end),1).soundOnsetTimesWrtSegm(end); 

                currPatIOIs(end+1) = ioiToNextSegm; 
            end

            % assign the resutlting IOIs for this pattern 
            seq(patIdx(pati),1).IOIs = currPatIOIs; 

            % update sound-event counter for this segment 
            c = c + seq(patIdx(pati),1).nSounds; 

            % make audio for the pattern 
            [patternAudio,~] = makeStimMainExp(seq(patIdx(pati),1).pattern, ...
                                               currPatIOIs(1:end-1), ...
                                               cfg, ...
                                               seq(patIdx(pati),1).gridIOI, ...
                                               seq(patIdx(pati),1).F0,...                                      
                                               seq(patIdx(pati),1).taskF0idx, ...
                                               seq(patIdx(pati),1).isTask, ... 
                                               seq(patIdx(pati),1).patternAmp);

            % update pattern onset time wrt whole sequence                                
            seq(patIdx(pati),1).onset = currTimePoint+seq(patIdx(pati),1).soundOnsetTimesWrtSegm(1);

            % get current audio index in the sequence, and append the audio
            currAudioIdx = round( seq(patIdx(pati),1).soundOnsetTimesWrtSegm(1) * cfg.fs); 

            % assign pattern audio to segment audio
            segmAudio(currAudioIdx+1:currAudioIdx+length(patternAudio)) = patternAudio; 

        end

        % get index for the current segment onset wrt sequence
        currAudioIdx = round(currTimePoint*cfg.fs); 

        % we only put the audio data in the first structure in the
        % array of structures to save memory...
        seq(1).outAudio(currAudioIdx+1:currAudioIdx+length(segmAudio)) = segmAudio; 

        % update current time point in the sequence 
        currTimePoint = currTimePoint + cfg.pattern.interSegmInterval; 

        % add delay after each category (if applicable)
        % by shifting the current time point with delay
        if strcmpi(currCategLabel,cfg.pattern.labelCategA)
            currTimePoint = currTimePoint + cfg.pattern.delayAfterA;         
        elseif strcmpi(currCategLabel,cfg.pattern.labelCategB)
            currTimePoint = currTimePoint + cfg.pattern.delayAfterB;         
        end

        
    end % segment loop
    

end % step loop




