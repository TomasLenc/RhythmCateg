function [tapOnsets, responseEvents] = mb_getResponse(cfg, expParam, responseEvents, currSeq)

% allocate vector of tap times
tapOnsets = [];

% counter to count how many taps were made
cTap = 1; 

% boolean helper variable used to determine if the button was just
% pressed (and not held down from previous loop iteration)
istap = false;


currSeqStartTime = expParam.currSeqStartTime;


% stay in the loop until the sequence ends
while GetSecs < (currSeqStartTime+cfg.SequenceDur)

        % check if key is pressed
        [~, tapTime, keyCode] = KbCheck(cfg.keyboard);

        % terminate if quit-button pressed
        if find(keyCode)==cfg.keyquit
            error('Experiment terminated by user...');
        end

        % check if tap and save time (it counts as tap if
        % reponse buttons were released initially)
        if ~istap && any(keyCode)
            % tap onset relative to sequence start time
            tapTime = tapTime-currSeqStartTime; 
            
            % append to the vector of tap times
            tapOnsets = [tapOnsets,tapTime];
            
            % set flag
            istap = true;

            
            % ------ bids -------
            % Find latest pattern for which current tap time is larger 
            % than it's onset. 
            currPatIdx = max( find(tapTime > [currSeq.onset]) ); 

            responseEvents(cTap,1).sequenceNum     = expParam.seqi;
            responseEvents(cTap,1).onset           = tapTime;
            responseEvents(cTap,1).duration        = 0;
            responseEvents(cTap,1).trial_type      = 'response';
            responseEvents(cTap,1).patternID       = currSeq(currPatIdx).patternID;
            responseEvents(cTap,1).segmCateg       = currSeq(currPatIdx).segmCateg;
            responseEvents(cTap,1).segmentNum      = currSeq(currPatIdx).segmentNum;
            responseEvents(cTap,1).segmentOnset    = currSeq(currPatIdx).segmentOnset;
            responseEvents(cTap,1).stepNum         = currSeq(currPatIdx).stepNum;
            responseEvents(cTap,1).stepOnset       = currSeq(currPatIdx).stepOnset;
            
            responseEvents(cTap,1).F0              = currSeq(currPatIdx).F0;
            responseEvents(cTap,1).gridIOI         = currSeq(currPatIdx).gridIOI;
            responseEvents(cTap,1).patternAmp      = currSeq(currPatIdx).patternAmp;
            
            % get pattern info
            responseEvents(cTap,1).PE4        = currSeq(currPatIdx).PE4;
            responseEvents(cTap,1).minPE4     = currSeq(currPatIdx).minPE4;
            responseEvents(cTap,1).rangePE4   = currSeq(currPatIdx).rangePE4;
            responseEvents(cTap,1).LHL24      = currSeq(currPatIdx).LHL24;
            responseEvents(cTap,1).minLHL24   = currSeq(currPatIdx).minLHL24;
            responseEvents(cTap,1).rangeLHL24 = currSeq(currPatIdx).rangeLHL24;
            
            
            % increase tap counter
            cTap = cTap+1; 
        end
        
        if istap && ~any(keyCode)
            istap = false;
        end

end
