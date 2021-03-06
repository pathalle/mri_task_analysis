function [ onsets, half, correct, astim, stimon_mri, feedon_mri ] = get_eread_onsets_fbl( logfile )
    fileID = fopen(logfile);
    content = textscan(fileID,'%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s ','Delimiter','\t');
    
    %congruent = content{5}(3:end);
    
    condition = content{7}(2:end);
    condition(char(condition) == '2') = {'0'};
    
    stimon_mri = str2num(char(content{12}(2:end)))./1000;
    feedon_mri = str2num(char(content{14}(2:end)))./1000;
    rt         = str2num(char(content{6}(2:end)))./1000;
    astim      = content{4}(2:end);
    
    half = [];
    
    s = containers.Map
    s('1')=0;
    s('2')=0;
    s('3')=0;
    s('4')=0;
    s('5')=0;
    for i=1:length(astim)
        aud = astim(i);
        aud = cell2mat(aud);
        s(aud) = s(aud)+1;
        if s(aud) < 5
            half(i) = 1;
        else
            half(i) = 2;
        end
    end
    half = transpose(half);
    half = string(half);
    
    onsets = {};
    
    cond_names = {'First half', 'Second half','Feedback positive','Feedback negative'};
    
    for i = 1:length(cond_names)
        switch i
            case 1
                % find stimuli presented <= 4 times
                index = find(strcmp(half,'1'));
                onsets{1} = stimon_mri(index) ;
            case 2
                % find stimuli presented > 4 times
                index = find(strcmp(half,'2'));
                onsets{2} = stimon_mri(index) ;
            case 3
                % find incorrect trials
                index = find(strcmp(condition,'0'));
                onsets{3} = feedon_mri(index) ;
            case 4
                % find correct trials
                index = find(strcmp(condition,'1'));
                onsets{4} = feedon_mri(index) ;    
        end   
    end   
    
    % rt congruent
    index = find(strcmp(condition,'1'));
    rt_pos_fb = rt(index);
    
    % rt incongruent
    index = find(strcmp(condition,'0'));
    rt_neg_fb = rt(index);
    
    % convert congruent
    condition = str2num(char(condition));
    
    % convert correct / result
    correct = str2num(char(content{7}(2:end)));
    correct(correct == 2) = 0;
    
    % convert Vstim
    astim = str2num(char(astim));
    
    onsets{5} =rt_pos_fb;
    onsets{6} =rt_neg_fb ; 
    
    %onsets = vec2mat(str2num(str2mat(onsets{:})),5)'./1000;
    fclose(fileID);
end