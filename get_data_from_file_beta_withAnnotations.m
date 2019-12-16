%The function obtains a sequence as a clean char array from the input file
%(which should be .gb or .ape) Moreover, the function can give a template
%(annt) for a new (potentially assembled) ape/gb-file, or only the features
%(features_only) (which is useful for inserts) 

function [file_path, annt, features_only] = get_data_from_file(file_path)
for jj = 1:length(file_path) 
    file = fileread(file_path{jj}); 
    seq = file((strfind(file, 'ORIGIN')+6):end); %retrieve the sequence part
    %Clean seq from trash characters and spaces:
    trash = {'0','1','2','3','4','5','6','7','8','9', '/'};
    for ii = 1:length(trash)
        seq = erase(seq, trash{ii});
    end
    seq(isspace(seq)) = [];
    file_path{jj} = seq; %rewrite the initial array of pathways to files 
    %with sequences
    
    
    
    
    %Get needy annotations and features   
    annt = {'LOCUS', {}; [], 'bp ds-DNA'; 'circular', upper(date); ...
        'FEATURES', 'Location/Qualifiers'};
    misc = strfind(file, '..'); %.. is specific for two numbers separation
    features_only = cell(length(misc)*2,2);
    punctuation = {'_','-','–','.','/', ',','(', ')'}; %to erase these symbols
    nn = 1;
    for ii = 1:(length(misc)- numel(strfind(file, 'join(')))
       %Extract order # of nucleotides to change them in future (e.g. after
       %insertion)
       n1 = file((misc(ii)-10):(misc(ii)-1)); 
       if isempty(regexpi(n1, 'join')) == 0 %if the feature determined by joined bp
           n1(isspace(n1)) = []; n1(isletter(n1)) = []; 
           for j = 1:length(punctuation)
               n1 = erase(n1, punctuation{j});
           end
           n2 = file((misc(ii)+2):(misc(ii)+6));
           n2 = n2(1:(regexpi(n2, ',')-1));
           n2(isspace(n2)) = [];
           for j = 1:length(punctuation)
               n2 = erase(n2, punctuation{j});
           end
           n1 = str2double(n1); n2 = str2double(n2);
           features_only{nn, 2} = (min([n1 n2])):(max([n1 n2])); 
           feat_name = file((misc(ii)-20):(misc(ii)-6));
           feat_name = erase(feat_name, {'join'});
           feat_name(isspace(feat_name)) = [];
           features_only{nn, 1} = feat_name;
           features_only{nn, 1} = ['      ', features_only{nn, 1}]; %needed for correct 
       %ape/gb-file syntax (otherwise ApE does not inread the file
       misc(ii+1) = [];
       features_only(end-1:end, :) = []; %delete unneedy cells in this case
       else
           n1(isspace(n1)) = []; n1(isletter(n1)) = []; 
           for j = 1:length(punctuation)
               n1 = erase(n1, punctuation{j});
           end
           n2 = file((misc(ii)+2):(misc(ii)+7)); 
           n2(isspace(n2)) = []; n2(isletter(n2)) = [];
           for j = 1:length(punctuation)
               n2 = erase(n2, punctuation{j});
           end
           features_only{nn, 2} = str2double(n1):str2double(n2);
           feat_name = file((misc(ii)-20):(misc(ii)-6));
           feat_name(isspace(feat_name)) = [];
           features_only{nn, 1} = feat_name;
           features_only{nn, 1} = ['      ', features_only{nn, 1}]; %needed for correct 
       %ape/gb-file syntax (otherwise ApE does not inread the file)
       end
       nn = nn + 2; %record each second line
    end
    misc = strfind(file, '/locus_tag'); %record features info in analogous way
    nn = 2;
    for ii = 1:(length(misc))
        if ii == length(misc)
            features_only{nn, 1} = extractBetween(file, misc(ii), ...
                regexpi(file, 'ORIGIN')-1);
        features_only{nn, 1} = ['                     ', char(features_only{nn, 1})];
        else
        features_only{nn, 1} = extractBetween(file, misc(ii), misc(ii+1)-51);
        features_only{nn, 1} = ['                     ', char(features_only{nn, 1})];
        nn = nn + 2;
        end
    end
    %Put all the rest fields needed + space for the future assembled seq
    annt = [annt; features_only; {'ORIGIN', {}; {}, {}; {'//'}, {}}];
    
end
end