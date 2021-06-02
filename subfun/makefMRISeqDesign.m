function cfg = makefMRISeqDesign(cfg)

% this function creates counterbalanced audio sequences for fMRI in
% by using subfunction as:
% makeSequecen.m
% getAllSeqDesign.m
% makeStimMainExp.m

% It is embedded into getParams.m script but since it's
% depending on the expParam.runNb parameter, one should be causious on
% changing the function's location.(e.g. after the script gets runNb)

%%%%%%%%%%%%
% ! important, the order of arguments matters ! -> getAllSeq(categA, categB, ...)
% getAllSeqDesign(cfg.patternSimple, cfg.patternComplex, cfg, expParam);
%%%%%%%%%%%%

%%%%%%%%%%%%
% keep in mind:
% DesginFullExp (runNum, stepNum,segmentNum,patternNum)
%%%%%%%%%%%%


% if debug, put back run =1 so in the main script sequence =runNb ==1
runNb = cfg.subject.runNb;

if cfg.debug.do
    runNb = 1;
end

% path to save design matrix
saveFileName = [cfg.fileName.base,'_SeqDesign'];
saveFile = fullfile(fileparts(mfilename('fullpath')),'../../output',saveFileName);


%% Get counterbalanced sequences according to the total fMRI RUNs


if runNb == 1
    
    % get the design
    [DesignCateg, DesignSegment, DesignToneF0] = getAllSeqDesign(...
        cfg.pattern.patternA,...
        cfg.pattern.patternB, cfg);
    
    % get design pseudorandomized A/B if it's BlockDesign
    if strcmp(cfg.task.name,'RhythmBlock')
        [DesignCateg, DesignSegment] = addRandomCategOrder(cfg, ...
            DesignCateg, ...
            DesignSegment);
    end
       
    % add Task according to SegmentLabels
    cfg = addRandomizedTask(cfg,DesignSegment,cfg.pattern.numSequences);
    
    %save the Design
    
    %check if .mat file exists & give an error for overwrite
    if exist([saveFile,'.mat'],'file')
        reply = input('Do you want to re-run design matrix by runNb = 1? y/n :','s');
        if strcmp(reply,'y')
            disp('Okay! I''m overwriting your design matrix, runNb ==1\n');
        else
            error('Stopping the overwrite, you clearly did press the wrong runNb!');
        end
    end
    
    save(saveFile,'DesignCateg','DesignSegment','DesignToneF0','cfg');
    cfg.pattern.seqDesignFullExp = DesignCateg;
    cfg.pattern.seqDesignSegment = DesignSegment;
    cfg.pattern.seqDesignToneF0 = DesignToneF0;
    
else
    
    design = load(saveFile);
    cfg.pattern.seqDesignFullExp = design.DesignCateg;
    cfg.pattern.taskIdxMatrix = design.cfg.pattern.taskIdxMatrix; 
    cfg.pattern.seqDesignSegment = design.DesignSegment;
    cfg.pattern.seqDesignToneF0 = design.DesignToneF0;

    
end


%% Adding extra sessions design matrix
if runNb > cfg.pattern.numSequences && mod(runNb,3)==1
    
    %create design matrix
    [extraSeqDesign,extraSeqSegment,extraSeqToneF0] = getAllSeqDesign(cfg.pattern.patternA,...
        cfg.pattern.patternB, cfg, cfg.pattern.extraSeqNum);
    
    %create task matrix
    extracfg = addRandomizedTask(cfg,extraSeqSegment,cfg.pattern.extraSeqNum);
    
    % add and assign new design with task
    cfg.pattern.seqDesignFullExp = [cfg.pattern.seqDesignFullExp; ...
        extraSeqDesign];
    cfg.pattern.taskIdxMatrix = [cfg.pattern.taskIdxMatrix; ...
        extracfg.pattern.taskIdxMatrix];
    
    cfg.pattern.seqDesignSegment = [cfg.pattern.seqDesignSegment; ...
        extraSeqSegment]; 
    
    cfg.pattern.seqDesignToneF0 = [cfg.pattern.seqDesignToneF0; ...
        extraSeqToneF0];
    
    
    %save to be called later in runs
    DesignCateg = cfg.pattern.seqDesignFullExp;
    DesignSegment = cfg.pattern.seqDesignSegment;
    DesignToneF0 = cfg.pattern.seqDesignToneF0;
    save(saveFile,'DesignCateg','DesignSegment',...
                                'cfg','extracfg','extraSeqSegment',...
                                'DesignToneF0','extraSeqToneF0');
    
    fprintf('new sequence design and task added! Wohoo!\n\n');
end


end

function cfg = addRandomizedTask(cfg,Design,numSequence)

%create an empty cell to store the task==1s and 0s
taskMatrix =  zeros(...
    numSequence, ...
    cfg.pattern.nStepsPerSequence,...
    cfg.pattern.nSegmPerStep, ...
    cfg.pattern.nPatternPerSegment);
    
% find the categA and categB from SegmentLabel
idxCategA = contains(Design(:),cfg.pattern.labelSegmentA);
idxCategB = contains(Design(:),cfg.pattern.labelSegmentB);

%count the number of patterns categA and categB
categANum = sum(idxCategA);
categBNum = sum(idxCategB);


% take the 10%
% of patterns
cfg.pattern.categANumTarget = round(categANum*0.1);
cfg.pattern.categBNumTarget = round(categBNum*0.1);


%create zero array
categA = zeros(categANum,1);
categB = zeros(categBNum,1);


%assign 1s to indicate the targets
categA(1:cfg.pattern.categANumTarget) = 1;
categB(1:cfg.pattern.categBNumTarget) = 1;


% make shuffled matrix with 1s and 0s
taskMatrix = shuffleMe(taskMatrix, idxCategA, idxCategB, categA, categB);


% control point 1
% for all the beginning on runs == beginning of
for iRun=1:length(taskMatrix)
    
    while taskMatrix(iRun) == 1
        
        taskMatrix = shuffleMe(taskMatrix, idxCategA, idxCategB, categA, categB);
        
    end
    if taskMatrix(iRun)
        sprintf('There''s a target in the first pattern!');
    end
    
end

% control point 2
% 2. omit repetition 1 1 : sequential targets
% find the repeting 1s
taskIndices = find(taskMatrix == 1);

% shuffle the matrix till no more sequential repeating 1s
while taskIndices(diff(taskIndices)==1)
    
    taskMatrix = shuffleMe(taskMatrix, idxCategA, idxCategB, categA, categB);
    taskIndices = find(taskMatrix == 1);
    
end
    

% control point 3 not yet implemented:
% if sum([taskIdxMatrix(irun,:,:,:)]) > 12 shuffle again

% lastly, assign it into cfg to call later
cfg.pattern.taskIdxMatrix = taskMatrix;

end

function [DesignCateg, DesignSegment] = addRandomCategOrder(cfg, DesignCateg, DesignSegment)

numRun = cfg.pattern.numSequences;

% if subject nb is odd choose 4 runs to change order else 5 runs
% numRunToChange = round(numRun/2);
if mod(cfg.subject.subjectNb,2)
    numRunToChange = 4;
else
    numRunToChange = 5;
end


%get shuffled
idxRuns = Shuffle(1:numRun);
idxRuns = idxRuns(1:numRunToChange);

for iRun= 1:numRunToChange
    
    currRun = idxRuns(iRun);
    
    % change Design and Segment
    % dims: [seq x step x segm x pat]
    % 1x5x2x4
    currRunDesign = DesignCateg(currRun,:,:,:);
    % flip segments - dim:3
    flipCurrRunDesign = flip(currRunDesign, 3);
    DesignCateg(currRun,:,:,:) = flipCurrRunDesign;
    
    % do the same for SegmentLabel
    currRunSegment = DesignSegment(currRun,:,:,:);
    flipCurrRunSegment = flip(currRunSegment, 3);
    DesignSegment(currRun,:,:,:) = flipCurrRunSegment;
end
end


function taskMatrix = shuffleMe(taskMatrix, idxCategA, idxCategB, categA, categB)

%and shuffle the order or target across seq (runs), steps, segments, ...
idxCategATarget = Shuffle(categA);
idxCategBTarget = Shuffle(categB);

% put conditions together in order
taskMatrix(idxCategA)= idxCategATarget;
taskMatrix(idxCategB)= idxCategBTarget;

end








