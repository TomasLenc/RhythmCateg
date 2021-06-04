function collectAndSave(responseEvents, cfg, logFile, experimentStart)

    target = cfg.target;
    

    if isfield(responseEvents(1), 'onset') && ~isempty(responseEvents(1).onset)

        for iResp = 1:size(responseEvents, 1)
            responseEvents(iResp).onset = ...
                responseEvents(iResp).onset - experimentStart;
            responseEvents(iResp).target = target;

        end

        responseEvents(1).isStim = logFile.isStim;
        responseEvents(1).fileID = logFile.fileID;
        responseEvents(1).extraColumns = logFile.extraColumns;
                
        saveEventsFile('save', cfg, responseEvents);

    end
end
