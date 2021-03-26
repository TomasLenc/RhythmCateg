% let's run behavioral experiment


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

subject = 77;

tasks = {'RhythmFT','RhythmBlock'};

repetitionNb = 6;
taskLabel = [1 2];

runStartFrom = 1;
runEndIn = repetitionNb;

% for eve run numbered subjects, flip the order
if mod(subject,2) == 0
   taskLabel = [2 1];
end


% run the experiments
count = 1;
for irun = runStartFrom:runEndIn
    runNb = irun;
    
    for iTask = 1: length(taskLabel)
        
        task = tasks{taskLabel(iTask)};
        
      whereIsData = tapMainExperiment(task, subject, runNb);
       
      expOrderOfSubject(count).taskName = task;
      expOrderOfSubject(count).LogFileLocation = whereIsData;
       count = count + 1;
    end
        
end


% save
save(fullfile(fileparts(whereIsData),['sub-0',num2str(subject),'expOrder.mat']),'expOrderOfSubject');