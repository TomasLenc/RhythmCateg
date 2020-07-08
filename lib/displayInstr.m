function [subAction] = displayInstr(instrTxt,cfg,varargin)
% varargin:   
% 
%     'setVolume' displays intrTxt and offers possibility for the
%     participant to set the volume 
% 
%     'setVolumeAndGeneralInstrOption' displays intrTxt and offers possibility for the
%     participant to set the volume, or toggle general instructions for the
%     expeiment
%         - !!! in varargin, provide also: 
%                 generalInstrTxt :   string, 
%                                     text of the genral task instruction to display
%                                     if toggled
% 
%     'waitForKeypress' displays intrTxt and waits for keypress to continue
% 
%     'instrAndQuitOption' displays intrTxt and info about experiment
%     termination (quit option)
% 
%     if no varargin provided, the function only displays intrTxt on the screen and returns


% This optional variable can be used to inform the script about particpant's action.
% E.g. if participants wants to move one instruction-page back, we can let the script
% know through this output variable. 
subAction = []; 


Screen('TextFont',cfg.win,cfg.textFont);
Screen('TextSize',cfg.win,cfg.textSize);

%% volume setting  
if any(strcmpi(varargin,'setVolume'))
    
    % put volume-setting audio into buffer
    PsychPortAudio('FillBuffer',cfg.pahandle,cfg.volumeSettingSound); 
    
    % set allowed keys  
    allowedKeys = [cfg.keyAudioPlay, cfg.keyAudioStop, ...
                   cfg.keyVolDown, cfg.keyVolUp, ...
                   cfg.keywait, cfg.keyquit]; 
               
    % wait for participatnt's response          
    while 1
        
        % display instructions in the center of cfg.screen 
        DrawFormattedText(cfg.win,instrTxt,'center','center',cfg.white); 

        % display continue option on the bottom of the screen
        DrawFormattedText(cfg.win,'Press [ENTER] to continue...','center',cfg.winHeight*0.9,cfg.white); 

        % display quit option in the cfg.screen corner
        txt = sprintf([sprintf(' Press [%s] to play test sound\n',KbName(cfg.keyAudioPlay)), ...
               sprintf(' Press [%s] to stop test sound\n',KbName(cfg.keyAudioStop)), ...
               sprintf(' Press [%s] to increase volume\n',KbName(cfg.keyVolUp)), ...
               sprintf(' Press [%s] to decrease volume\n',KbName(cfg.keyVolDown)), ...
               sprintf('\n Volume = %d%%\n',round(PsychPortAudio('Volume',cfg.pahandle)*100)), ...
               ]); 
        tbx     = Screen('TextBounds', cfg.win, txt);
        width   = tbx(3);
        height  = tbx(4);
        r = [0 0 width height + Screen('TextSize', cfg.win)];
        r = AlignRect(r,cfg.winRect,RectLeft,RectTop);
        DrawFormattedText(cfg.win, txt, r(RectLeft), r(RectBottom), cfg.white);

        
        % display quit option 
        txt = sprintf('press [%s] to quit the whole experiment  ',KbName(cfg.keyquit)); 
        tbx     = Screen('TextBounds', cfg.win, txt);
        width   = tbx(3);
        height  = tbx(4);
        r = [0 0 width height + Screen('TextSize', cfg.win)];
        r = AlignRect(r,cfg.winRect,RectRight,RectTop);
        DrawFormattedText(cfg.win, txt, r(RectLeft), r(RectBottom), cfg.white);
        
        
        Screen('Flip', cfg.win);      

        % wait for keypress
        keyCodePressed = waitForKeyKbCheck(allowedKeys); 
        
        if ismember(keyCodePressed, cfg.keyAudioPlay)
            PsychPortAudio('Stop',cfg.pahandle); 
            PsychPortAudio('Start',cfg.pahandle); 
            
        elseif ismember(keyCodePressed, cfg.keyAudioStop)
            PsychPortAudio('Stop',cfg.pahandle); 

        elseif ismember(keyCodePressed, cfg.keyVolUp)         
            oldVolume = PsychPortAudio('Volume',cfg.pahandle); 
            newVolume = min(oldVolume+0.05, 1); 
            PsychPortAudio('Volume',cfg.pahandle,newVolume);
            
        elseif ismember(keyCodePressed, cfg.keyVolDown)
            oldVolume = PsychPortAudio('Volume',cfg.pahandle); 
            newVolume = max(oldVolume-0.05, 0); 
            PsychPortAudio('Volume',cfg.pahandle,newVolume);
 
        elseif ismember(keyCodePressed, cfg.keywait)
            PsychPortAudio('Stop',cfg.pahandle); 
            break
            
        elseif ismember(keyCodePressed, cfg.keyquit)
            PsychPortAudio('Stop',cfg.pahandle); 
            error('experiment terminated by user'); 
            
        end
        
    end
    
    
    
    
 
%% volume setting and option to go previous/next in the instruction pages
elseif any(strcmpi(varargin,'setVolumePrevNext'))
%     
%      !!! IMPORTANT !!! -> specify in varargin: 
%         currInstrPage :     int
%                             which instruction page are we currently at? 
%         nInstrPages :       int
%                             how many instruction pages are there before
%                             continuing? 
    
    % put volume-setting audio into buffer
    PsychPortAudio('FillBuffer',cfg.pahandle,cfg.volumeSettingSound); 
    
    % set allowed keys  
    allowedKeys = [cfg.keyAudioPlay, cfg.keyAudioStop, ...
                   cfg.keyVolDown, cfg.keyVolUp, ...
                   cfg.keyInstrBack, cfg.keyInstrNext, ...
                   cfg.keywait, cfg.keyquit]; 
   
    % get current instruction page index from varargin
    if any(strcmpi(varargin,'currInstrPage'))
        currInstrPage = varargin{find(strcmpi(varargin,'currInstrPage'))+1}; 
    else
        warning('instruction page index not specified'); 
        currInstrPage = [0]; 
    end
    % get total number of instruction pages from varargin
    if any(strcmpi(varargin,'nInstrPages'))
        nInstrPages = varargin{find(strcmpi(varargin,'nInstrPages'))+1}; 
    else
        warning('number of instruction pages not specified'); 
        nInstrPages = [0]; 
    end
               
    % wait for participatnt's response          
    while 1
        
        % display instructions in the center of cfg.screen 
        DrawFormattedText(cfg.win,instrTxt,'center','center',cfg.white); 

        % display continue option on the bottom of the screen
        % ! only if we're at the last instruction page ! 
        if currInstrPage == nInstrPages
            DrawFormattedText(cfg.win,'Press [ENTER] to continue...','center',cfg.winHeight*0.9,cfg.white); 
        end
        
        % display quit option in the cfg.screen corner
        txt = sprintf([sprintf(' Press [%s] to play test sound\n',KbName(cfg.keyAudioPlay)), ...
               sprintf(' Press [%s] to stop test sound\n',KbName(cfg.keyAudioStop)), ...
               sprintf(' Press [%s] to increase volume\n',KbName(cfg.keyVolUp)), ...
               sprintf(' Press [%s] to decrease volume\n',KbName(cfg.keyVolDown)), ...
               sprintf('\n Volume = %d%%\n',round(PsychPortAudio('Volume',cfg.pahandle)*100)), ...
               ]); 
        tbx     = Screen('TextBounds', cfg.win, txt);
        width   = tbx(3);
        height  = tbx(4);
        r = [0 0 width height + Screen('TextSize', cfg.win)];
        r = AlignRect(r,cfg.winRect,RectLeft,RectTop);
        DrawFormattedText(cfg.win, txt, r(RectLeft), r(RectBottom), cfg.white);
        
        % display quit option 
        txt = sprintf('press [%s] to quit the whole experiment  ',KbName(cfg.keyquit)); 
        tbx     = Screen('TextBounds', cfg.win, txt);
        width   = tbx(3);
        height  = tbx(4); 
        r = [0 0 width height + Screen('TextSize', cfg.win)];
        r = AlignRect(r,cfg.winRect,RectRight,RectTop);
        DrawFormattedText(cfg.win, txt, r(RectLeft), r(RectBottom), cfg.white);
        
        % display previous option for instruction pages
        txt = sprintf( [sprintf('page %d/%d \n', currInstrPage, nInstrPages), ...
                        sprintf('press [%s] to go forward\n', KbName(cfg.keyInstrNext)), ...
                        sprintf('press [%s] to go back\n', KbName(cfg.keyInstrBack)), ...
                        ]); 
        nLines  = 3; 
        tbx     = Screen('TextBounds', cfg.win, txt);
        width   = tbx(3);
        height  = tbx(4);
        r = [0 0 width height + nLines*Screen('TextSize', cfg.win)];
        r = AlignRect(r,cfg.winRect,RectRight,RectBottom);
        DrawFormattedText(cfg.win, txt, r(RectLeft), r(RectTop), cfg.white);
        
        
        
        Screen('Flip', cfg.win);      

        % wait for keypress
        keyCodePressed = waitForKeyKbCheck(allowedKeys); 
        
        if ismember(keyCodePressed, cfg.keyAudioPlay)
            PsychPortAudio('Stop',cfg.pahandle); 
            PsychPortAudio('Start',cfg.pahandle); 
            
        elseif ismember(keyCodePressed, cfg.keyAudioStop)
            PsychPortAudio('Stop',cfg.pahandle); 

        elseif ismember(keyCodePressed, cfg.keyVolUp)         
            oldVolume = PsychPortAudio('Volume',cfg.pahandle); 
            newVolume = min(oldVolume+0.05, 1); 
            PsychPortAudio('Volume',cfg.pahandle,newVolume);
            
        elseif ismember(keyCodePressed, cfg.keyVolDown)
            oldVolume = PsychPortAudio('Volume',cfg.pahandle); 
            newVolume = max(oldVolume-0.05, 0); 
            PsychPortAudio('Volume',cfg.pahandle,newVolume);
 
        elseif ismember(keyCodePressed, cfg.keywait) && currInstrPage == nInstrPages
            PsychPortAudio('Stop',cfg.pahandle); 
            subAction = 'done'; 
            break
            
        elseif ismember(keyCodePressed, cfg.keyInstrNext)
            PsychPortAudio('Stop',cfg.pahandle); 
            subAction = 'oneInstrPageForward'; 
            break

        elseif ismember(keyCodePressed, cfg.keyInstrBack)
            subAction = 'oneInstrPageBack'; 
            break

        elseif ismember(keyCodePressed, cfg.keyquit)
            PsychPortAudio('Stop',cfg.pahandle); 
            error('experiment terminated by user'); 
            
        end
        
    end
 
    
    
    
       
    
%% display instructions, give option to 1) set volume, 2) display general instruction  
elseif any(strcmpi(varargin,'setVolumeToggleGeneralInstr'))
    
    % look for general intruction text in varargin
    if any(strcmpi(varargin,'generalInstrTxt'))
        generalInstrTxt = varargin{find(strcmpi(varargin,'generalInstrTxt'))+1}; 
    else
        error('You need to call displayInstr function with ''generalInstrTxt'' in the varargin. Otherwise I don''t know what general instructions to display...'); 
    end
    
    % put volume-setting audio into buffer
    PsychPortAudio('FillBuffer',cfg.pahandle,cfg.volumeSettingSound); 
    
    % set allowed keys  
    allowedKeys = [cfg.keyAudioPlay, cfg.keyAudioStop, ...
                   cfg.keyVolDown, cfg.keyVolUp, ...
                   cfg.keywait, cfg.keyToggleInstr, cfg.keyquit]; 
      
   % boolean used to toggle general instruction display 
   generalInstrOnScreen = false; 
    
    % wait for participatnt's response          
    while 1
        
        if generalInstrOnScreen
        % we are displaying general instructions for the experiment
            
            % display instructions in the center of cfg.win 
            DrawFormattedText(cfg.win,generalInstrTxt,'center','center',cfg.white); 

            % display continue option on the bottom of the screen
            DrawFormattedText(cfg.win,'Press [I] to toggle instructions, or [ENTER] to continue...','center',cfg.winHeight*0.9,cfg.white); 
        
            
        else
        % we are displaying normal setVolume display 
            
            % display instructions in the center of cfg.screen 
            DrawFormattedText(cfg.win,instrTxt,'center','center',cfg.white); 

            % display continue option on the bottom of the screen
            DrawFormattedText(cfg.win,'Press [I] to toggle instructions, or [ENTER] to continue...','center',cfg.winHeight*0.9,cfg.white); 

            % display quit option in the cfg.screen corner
            txt = sprintf([sprintf(' Press [%s] to play test sound\n',KbName(cfg.keyAudioPlay)), ...
                   sprintf(' Press [%s] to stop test sound\n',KbName(cfg.keyAudioStop)), ...
                   sprintf(' Press [%s] to increase volume\n',KbName(cfg.keyVolUp)), ...
                   sprintf(' Press [%s] to decrease volume\n',KbName(cfg.keyVolDown)), ...
                   sprintf('\n Volume = %d%%\n',round(PsychPortAudio('Volume',cfg.pahandle)*100)), ...
                   ]); 
            tbx     = Screen('TextBounds', cfg.win, txt);
            width   = tbx(3);
            height  = tbx(4);
            r = [0 0 width height + Screen('TextSize', cfg.win)];
            r = AlignRect(r,cfg.winRect,RectLeft,RectTop);
            DrawFormattedText(cfg.win, txt, r(RectLeft), r(RectBottom), cfg.white);

            % display quit option 
            txt = sprintf('press [%s] to quit the whole experiment  ',KbName(cfg.keyquit)); 
            tbx     = Screen('TextBounds', cfg.win, txt);
            width   = tbx(3);
            height  = tbx(4);
            r = [0 0 width height + Screen('TextSize', cfg.win)];
            r = AlignRect(r,cfg.winRect,RectRight,RectTop);
            DrawFormattedText(cfg.win, txt, r(RectLeft), r(RectBottom), cfg.white);
            
                        
        end
        
        Screen('Flip', cfg.win);      

        % wait for keypress
        keyCodePressed = waitForKeyKbCheck(allowedKeys); 
        
        if ismember(keyCodePressed, cfg.keyAudioPlay)
            PsychPortAudio('Stop',cfg.pahandle); 
            PsychPortAudio('Start',cfg.pahandle); 
            
        elseif ismember(keyCodePressed, cfg.keyAudioStop)
            PsychPortAudio('Stop',cfg.pahandle); 

        elseif ismember(keyCodePressed, cfg.keyVolUp)         
            oldVolume = PsychPortAudio('Volume',cfg.pahandle); 
            newVolume = min(oldVolume+0.05, 1); 
            PsychPortAudio('Volume',cfg.pahandle,newVolume);
            
        elseif ismember(keyCodePressed, cfg.keyVolDown)
            oldVolume = PsychPortAudio('Volume',cfg.pahandle); 
            newVolume = max(oldVolume-0.05, 0); 
            PsychPortAudio('Volume',cfg.pahandle,newVolume);
 
        elseif ismember(keyCodePressed, cfg.keyToggleInstr)

            % update the bool
            generalInstrOnScreen = ~generalInstrOnScreen; 
                    
        elseif ismember(keyCodePressed, cfg.keywait)
            PsychPortAudio('Stop',cfg.pahandle); 
            break
            
        elseif ismember(keyCodePressed, cfg.keyquit)
            PsychPortAudio('Stop',cfg.pahandle); 
            error('experiment terminated by user'); 
            
        end
        
        
    end
            
    
    
    
    
    
    
    
    
    
    
    
    
%% display instructions and wait for keypress 
elseif any(strcmpi(varargin,'waitForKeypress'))
    
    % display instructions in the center of cfg.screen 
    DrawFormattedText(cfg.win,instrTxt,'center','center',cfg.white); 

    % display continue option on the bottom of the screen
    DrawFormattedText(cfg.win,'Press [ENTER] to continue...','center',cfg.winHeight*0.9,cfg.white); 

    % display quit option 
    txt = sprintf('press [%s] to quit the whole experiment  ',KbName(cfg.keyquit)); 
    tbx     = Screen('TextBounds', cfg.win, txt);
    width   = tbx(3);
    height  = tbx(4);
    r = [0 0 width height + Screen('TextSize', cfg.win)];
    r = AlignRect(r,cfg.winRect,RectRight,RectTop);
    DrawFormattedText(cfg.win, txt, r(RectLeft), r(RectBottom), cfg.white);

    
    Screen('Flip', cfg.win);      
    
    
    % wait for participant's keypress to continue 
    keyCodePressed = waitForKeyKbCheck([cfg.keywait,cfg.keyquit]);  
    
    if ismember(keyCodePressed, cfg.keywait)
        PsychPortAudio('Stop',cfg.pahandle); 
        return
        
    elseif ismember(keyCodePressed, cfg.keyquit)
        PsychPortAudio('Stop',cfg.pahandle); 
        error('experiment terminated by user'); 

    end


    
    
%% display instructions and quit option but don't wait
elseif any(strcmpi(varargin,'instrAndQuitOption'))
%     if there is additional feedback you'd like to display, you can specify it in varargin
%         fbktxt:     string
%                     text of additional feecback that will appear on the bottom of
%                     the screen
    
    % look for additional feedback, if requested display it 
    if any(strcmpi(varargin,'fbktxt'))
        fbktxt = varargin{find(strcmpi(varargin,'fbktxt'))+1}; 
        % switcht to mono font
        Screen('TextFont', cfg.win, 'Consolas');
        DrawFormattedText(cfg.win, fbktxt, 'center', cfg.winHeight*0.75, cfg.white);          
        % switch font back to default
        Screen('TextFont',cfg.win,cfg.textFont);
    end
    
    % display instructions in the center of cfg.screen 
    DrawFormattedText(cfg.win,instrTxt,'center','center',cfg.white); 

    % display small-font quit option 
    Screen('TextSize',cfg.win,cfg.textSize*0.7);
    txt = sprintf('(in case of emergency, press [%s] to terminate the experiment)  ',KbName(cfg.keyquit)); 
    tbx     = Screen('TextBounds', cfg.win, txt,[],[],[],[]);
    width   = tbx(3);
    height  = tbx(4);
    r = [0 0 width height + Screen('TextSize', cfg.win)];
    r = AlignRect(r,cfg.winRect,RectRight,RectTop);
    DrawFormattedText(cfg.win, txt, r(RectLeft), r(RectBottom), cfg.white);    
    Screen('TextSize',cfg.win,cfg.textSize);
    
    Screen('Flip',cfg.win); 
    
    
%% display instructions 
else
    
    % display instructions in the center of cfg.screen 
    DrawFormattedText(cfg.win,instrTxt,'center','center',cfg.white); 

    Screen('Flip',cfg.win); 
    
end





