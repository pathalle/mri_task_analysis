function [ onsets, congruent, correct, vstim, stimon_mri, feedon_mri ] = get_eread_onsets_4( logfile )
    fileID = fopen(logfile);
    content = textscan(fileID,'%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s ','Delimiter','\t');
    
    %congruent = content{5}(3:end);
    
    condition = content{7}(2:end);
    
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
    
    onsets = [];
    cond_names = {'Feedback positive','Feedback negative'};
    
    for i = 1:length(cond_names)
        switch i
            case 1
                index = find(strcmp(congruent,'1'));
                onsets = [ onsets stimon_mri(index) ];
            case 2
                index = find(strcmp(congruent,'0'));
                onsets = [ onsets stimon_mri(index) ];
            case 3
                index = find(strcmp(congruent,'1'));
                onsets = [ onsets feedon_mri(index) ];
            case 4
                index = find(strcmp(congruent,'0'));
                onsets = [ onsets feedon_mri(index) ];    
        end   
    end   
    
    % rt congruent
    index = find(strcmp(congruent,'1'));
    rt_congruent = rt(index);
    
    % rt incongruent
    index = find(strcmp(congruent,'0'));
    rt_incongruent = rt(index);
    
    % convert congruent
    congruent = str2num(char(congruent));
    
    % convert correct / result
    correct = str2num(char(content{8}(2:end)));
    correct(correct == 2) = 0;
    
    % convert Vstim
    vstim = str2num(char(content{3}(2:end)));
    
    onsets = [ onsets rt_congruent rt_incongruent ]; 
    
    %onsets = vec2mat(str2num(str2mat(onsets{:})),5)'./1000;
    fclose(fileID);
end