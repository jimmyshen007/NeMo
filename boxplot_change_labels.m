fontSize = 18;    
    tickLabelStr = {'Label alpha','Label beta','Label chi','Label delta',...
                        'Label epsilon','Label fish','Label gamma','Label hallo','Label ingo'}
    % generate data
    final_res = 10*randn(300,9)+10;
    
    % group boxes
    width = 0.5;
    sh = 0.1; 
    
    pos = [1+sh 2-sh 3+sh 4-sh 5+sh 6-sh 7+sh 8-sh 9];
    wid = width * ones(1,length(pos));
    
    % boxplot
    figure
    boxplot(final_res, 'notch', 'on', ...
            'positions', pos,...
            'widths', wid) 
    
    % label, change fontsize
    % y-axis
    set(gca, 'FontSize', fontSize)
    ylim([-.5 40])
    ylabel('Error [%]', 'FontSize', fontSize)
    
    %x-labels
    text_h = findobj(gca, 'Type', 'text');
    rotation = 45;
    
    for cnt = 1:length(text_h)
        set(text_h(cnt),    'FontSize', fontSize,...
                            'Rotation', rotation, ...
                            'String', tickLabelStr{length(tickLabelStr)-cnt+1}, ...
                            'HorizontalAlignment', 'right')
    end
    
    % 'VerticalAlignment', 'cap', ...
    
    % smaller box for axes, in order to un-hide the labels
    squeeze = 0.2;
    left = 0.02;
    right = 1;
    bottom = squeeze;
    top = 1-squeeze;
    set(gca, 'OuterPosition', [left bottom right top])
        
    % remove outliers
    hout = findobj(gca,'tag','Outliers');
    for out_cnt = 1 : length(hout)
        set(hout(out_cnt), 'Visible', 'off')
    end