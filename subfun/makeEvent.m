function [s, envEvent] = makeEvent(cfg,varargin)

% mini function to only make the ramp-on and off shape

% INPUT cfg and F0
% it has to be provided:
% cfg.evenRampon duration
% cfg.Rampoff duration
% cfg.eventDur event duration 
% cfg.fs sampling frequency

% OUTPUT envEvent : event envelope shape

if nargin<2
    currF0 = cfg.pattern.F0s(1);
else
    currF0 =varargin{1};
end

% number of samples for the whole tone duration 
eventSamples   = round(cfg.pattern.eventDur * cfg.fs);

% number of samples for the onset ramp (proportion of gridIOI)
ramponSamples   = round(cfg.pattern.eventRampon * cfg.fs);

% number of samples for the offset ramp (proportion of gridIOI)
rampoffSamples  = round(cfg.pattern.eventRampoff * cfg.fs);

% individual sound event duration defined as proportion of gridIOI
envEvent = ones(1, round(cfg.pattern.gridIOIs * cfg.fs));

% make the linear ramps
envEvent(1:ramponSamples) = envEvent(1:ramponSamples) .* linspace(0,1,ramponSamples);

envEvent(eventSamples-rampoffSamples+1:eventSamples) = ... 
    envEvent(eventSamples-rampoffSamples+1:eventSamples) .* linspace(1,0,rampoffSamples);

envEvent(eventSamples+1:end) = 0; 

t = [0 : round(cfg.pattern.gridIOIs * cfg.fs)-1]/cfg.fs;

% create carrier
s = sin(2*pi*currF0*t);
% apply envelope to the carrier
s = s.* envEvent;
end

% % to visualise 1 pattern
% figure; plot(t,s);

