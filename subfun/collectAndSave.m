function collectAndSave(responseEvents, cfg, logFile, experimentStart)

    target = cfg.target;
%     responseEvents = getResponse('check', cfg.keyboard.responseBox, cfg);

    if isfield(responseEvents(1), 'onset') && ~isempty(responseEvents(1).onset)

        for iResp = 1:size(responseEvents, 1)
            responseEvents(iResp).onset = ...
                responseEvents(iResp).onset - experimentStart;
%             
%             responseEvents(iResp).sequenceNum = 'n/a';
%             responseEvents(iResp).segmentNum = 'n/a';
%             responseEvents(iResp).segmentOnset = 'n/a';
%             responseEvents(iResp).stepNum = 'n/a';
%             responseEvents(iResp).stepOnset = 'n/a';
%             responseEvents(iResp).patternID = 'n/a';
%             responseEvents(iResp).segmentCateg = 'n/a';
%             responseEvents(iResp).F0 = 'n/a';
%             responseEvents(iResp).isTask = 'n/a';
%             responseEvents(iResp).gridIOI = 'n/a';
%             responseEvents(iResp).patternAmp = 'n/a';
%             responseEvents(iResp).minPE4 = 'n/a';
%             responseEvents(iResp).rangePE4 = 'n/a';
%             responseEvents(iResp).minLHL24 = 'n/a';
%             responseEvents(iResp).rangeLHL24 = 'n/a';
%             responseEvents(iResp).LHL24 = 'n/a';
%             responseEvents(iResp).PE4 = 'n/a';
            responseEvents(iResp).target = sum(target);

        end

        responseEvents(1).fileID = logFile.fileID;
        responseEvents(1).extraColumns = logFile.extraColumns;
        
        responseEvents(1).isStim = logFile.isStim;
        
        saveEventsFile('save', cfg, responseEvents);

    end
end
