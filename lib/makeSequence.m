function seq = makeSequence(cfg,seqi,varargin)
% This function constructs a stimulus sequence.
% by using makeStimMainExp.m script

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

%% allocate variables to log

% main output structure (we'll put everything else into it)
seq = struct(); 

% vector of F0 (pitch) values for each pattern in the sequence
seq.F0 = zeros(1, cfg.nPatternPerSegment * cfg.nSegmPerStep * cfg.nStepsPerSequence); 

% vector of gridIOI values for each pattern in the sequence
seq.gridIOI = zeros(1, cfg.nPatternPerSegment * cfg.nSegmPerStep * cfg.nStepsPerSequence); 

% segment-category (A or B) for each pattern in the sequence 
seq.segmCateg = cell(1, cfg.nPatternPerSegment * cfg.nSegmPerStep * cfg.nStepsPerSequence); 

% onset time of each pattern
% i'm still conflicting how to make this onset both BIDS compatible &
% explicit that it's PATTERN ONSET WE ARE RECORDING
seq.onset = nan(1, cfg.nPatternPerSegment * cfg.nSegmPerStep * cfg.nStepsPerSequence); 

% put together all the patterns from both categories, we will pick from
% this using the unique ID of each pattern (we know which IDs we want from
% the output of getAllSeq function. 
patterns2choose = [cfg.patternSimple,cfg.patternComplex]; 


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

seq.patternID = cell(1, cfg.nPatternPerSegment * cfg.nSegmPerStep * cfg.nStepsPerSequence); 

% % THIS IS UNUSED
% % cell array, each element is a grid representation of the chosen pattern
% % (successively as the sequence unfolds)
% seq.outPatterns = cell(1, cfg.nPatternPerSegment * cfg.nSegmPerStep * cfg.nStepsPerSequence);
% % % 


% audio waveform of the sequence
seq.outAudio = zeros(1,round(cfg.SequenceDur*cfg.fs)); 

% % audio envelop  of the sequence
% seq.outEnvelop = zeros(1,round(cfg.SequenceDur*cfg.fs)); 

%% initialize counters 

% currently chosen F0 index (indexing value in cfg.F0s, initialize to 1)
currF0idx = 1; 

% currently chosen gridIOI index (indexing value in cfg.gridIOIs, initialize to 1)
currGridIOIidx = 1; 

% pattern counter (over the whole sequence)
cPat = 1; 

% current time point in the sequence as we go through the loops
currTimePoint = 0; 




%% loop over steps
for stepi=1:cfg.nStepsPerSequence
    
    
    % take the timestamp for logging the current step time
    stepOnset  = currTimePoint;
    
    
    %% loop over segments in 1 sequence
    % to make 
    for segmi=1:cfg.nSegmPerStep
        
        % take the  timestamp for logging the current segment time
        segmentOnset  = currTimePoint;
        
        % Determine which segment category this is (A or B), 
        % the first 'cfg.nSegmentA' segments will be category A, 
        % the rest will be B
        if ismember(segmi, [1:cfg.nSegmentA])            
            currCategLabel = cfg.labelCategA; 
        else
            currCategLabel = cfg.labelCategB; 
        end
        
        
        
        %% loop over pattern cycles in 1 segment
        % to create a segment
        for pati=1:cfg.nPatternPerSegment
    
            
            % --------------------------------------------------
            % ----- determine if gridIOI needs to be changed ---
            % --------------------------------------------------
            CHANGE_IOI = 0; 
            
            % gridIOI change requested every segment (and this is the first
            % cycle in the segment)
            if cfg.changeGridIOISegm && pati == 1
                CHANGE_IOI = 1; 
                
            % gridIOI change requested every category (and this is the first
            % cycle in the segment, and a category just changed from A->B
            % or B->A)
            elseif cfg.changeGridIOICategory && pati == 1 && ...
                    ( segmi == 1 || segmi == cfg.nSegmentA+1 )
                CHANGE_IOI = 1; 
                
            % gridIOI change requested every step 
            elseif cfg.changeGridIOIStep && segmi==1
                CHANGE_IOI = 1; 
                
            end
                        
            % if change of gridIOI requested, randomly choose a new gridIOI
            % do this only if there is more than 1 gridIOI to choose from
            if CHANGE_IOI && length(cfg.gridIOIs)>1
                % get gridIOI to choose from 
                gridIOI2ChooseIdx = 1:length(cfg.gridIOIs); 
                % remove gridIOI used in the previous iteration (to prevent
                % repetition in the sequence) 
                gridIOI2ChooseIdx(gridIOI2ChooseIdx==currGridIOIidx) = []; 
                % randomly select new IOI
                currGridIOIidx = randsample(gridIOI2ChooseIdx); 
            end
            
            currGridIOI = cfg.gridIOIs(currGridIOIidx); 
            
            
            % --------------------------------------------------
            % ----- determine if pitch needs to be changed -----
            % --------------------------------------------------
            CHANGE_PITCH = 0; 
            
            % pitch change requested in every pattern cycle
            if cfg.changePitchCycle
                CHANGE_PITCH = 1; 
                
            % pitch change requested in every segment (and this is the first
            % pattern cycle in the segment)
            elseif cfg.changePitchSegm && pati==1
                CHANGE_PITCH = 1; 
                
            % pitch change requested every category (and this is the first
            % pattern cycle in the segment, and a category just changed from A->B
            % or B->A)
            elseif cfg.changePitchCategory && pati==1 && ...
                    ( segmi==1 || segmi==cfg.nSegmentA+1 )
                CHANGE_PITCH = 1; 
                
            % pitch change requested only in every step 
            elseif cfg.changePitchStep && segmi==1
                CHANGE_PITCH = 1; 
                
            end
            
            % if change of pitch requested, randomly choose a new pitch
            % do this only if there is more than 1 F0 to choose from
            if CHANGE_PITCH && length(cfg.F0s)>1
                % get F0s to choose from 
                pitch2ChooseIdx = 1:length(cfg.F0s); 
                % remove F0 used in the previous iteration (to prevent
                % repetition in the sequence) 
                pitch2ChooseIdx(pitch2ChooseIdx==currF0idx) = []; 
                % randomly select new F0 idx
                currF0idx = randsample(pitch2ChooseIdx,1); 
            end
            
            currF0 = cfg.F0s(currF0idx); 
            currAmp = cfg.F0sAmps(currF0idx);
            
            
            % --------------------------------------------------
            % ----------------- make the audio -----------------
            % --------------------------------------------------
            
            % find the pattern ID from the seqDesignFullExp (output of
            % getAllSeq function)
            currPatID = cfg.seqDesignFullExp{seqi,stepi,segmi,pati}; 
            currPatIdx = find(strcmp(currPatID,{patterns2choose.ID})); 
            
            % do a quick check that the assigment of category labels is
            % consistent, if not, give a warning
            currPatternCateg = regexp(patterns2choose(currPatIdx).ID, '\D*(?=\d.*)', 'match'); 
            currPatternCateg = currPatternCateg{1}; 
            if ~strcmpi(currPatternCateg,currCategLabel)
                warning('mismatching category labels during sequence construction...'); 
            end
            
            % get the pattern
            currPattern = patterns2choose(currPatIdx).pattern;

            % make audio 
            [patternAudio,~] = makeStimMainExp(currPattern, ...
                                               cfg, currGridIOI, ...
                                               currF0,currAmp); 
            

            % get current audio index in the sequence, and append the audio
            currAudioIdx = round(currTimePoint*cfg.fs); 
            
           
            % we only put the audio data in the first structure in the
            % array of structures to save memory...
            seq(1).outAudio(currAudioIdx +1:currAudioIdx +...
                length(patternAudio)) = patternAudio; 
            
            seq(cPat,1).patternID   = currPatID;
            seq(cPat,1).segmCateg   = currCategLabel;
            seq(cPat,1).onset       = currTimePoint;
            seq(cPat,1).segmentNum  = segmi;
            seq(cPat,1).segmentOnset = segmentOnset;
            seq(cPat,1).stepNum    = stepi;
            seq(cPat,1).stepOnset  = stepOnset;
            
            seq(cPat,1).pattern     = currPattern; 
            seq(cPat,1).F0          = currF0;
            seq(cPat,1).gridIOI     = currGridIOI;
            seq(cPat,1).patternAmp  = currAmp;

            % get pattern info e.g. PE and LHL
            seq(cPat,1).PE4        = patterns2choose(currPatIdx).PE4;
            seq(cPat,1).minPE4     = patterns2choose(currPatIdx).minPE4;
            seq(cPat,1).rangePE4   = patterns2choose(currPatIdx).rangePE4;
            seq(cPat,1).LHL24      = patterns2choose(currPatIdx).LHL24;
            seq(cPat,1).minLHL24   = patterns2choose(currPatIdx).minLHL24;
            seq(cPat,1).rangeLHL24 = patterns2choose(currPatIdx).rangeLHL24;
            
            % --------------------------------------------------
            % update current time point
            currTimePoint = currTimePoint + cfg.interPatternInterval;         

            % increase pattern counter
            cPat = cPat+1; 
            
                       
            
        end % segment loop 
        
        
        % add delay after each category (if applicable)
        % by shifting the current time point with delay
        if strcmpi(currCategLabel,cfg.labelCategA)
            currTimePoint = currTimePoint + cfg.delayAfterA;         
        elseif strcmpi(currCategLabel,cfg.labelCategB)
            currTimePoint = currTimePoint + cfg.delayAfterB;         
        end
        

    end
    

end



