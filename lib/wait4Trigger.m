function wait4Trigger(cfg)

if strcmp(cfg.device,'scanner')
    
    fprintf('Waiting for trigger \n');
    
    % display instructions in the center of cfg.screen
    displayInstr('Waiting For Trigger...',cfg);
    
    
    triggerCounter=0;
    
    while triggerCounter < cfg.numTriggers
        
        [keyIsDown, ~, keyCode, ~] = KbCheck(-1);
        
        if strcmp(KbName(keyCode),cfg.triggerKey)
            
            triggerCounter = triggerCounter+1 ;
            
            fprintf('Trigger %s \n', num2str(triggerCounter));
            
            %display the trigger count
            displayInstr(['Experiment starting in ',num2str(4-triggerCounter),'...'],cfg);
            
            while keyIsDown
                [keyIsDown, ~, ~, ~] = KbCheck(-1);
            end
            
        end
    end
end
end
