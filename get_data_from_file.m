%The function obtains a sequence as a clean char array from the input file
%(which should be .gb or .ape) Moreover, the function can give a template
%(annt) for a new (potentially assembled) ape/gb-file, or only the features
%(features_only) (which is useful for inserts) 

function [file_path] = get_data_from_file(file_path)
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
    
end
end



