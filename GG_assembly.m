%% Put your data here (files supported: .gb and .ape)

%1. Choose TypeIIs restriction enzyme: 1 - BsaI, 2 - BbsI, 3 - Esp3I
enzyme = 1; 
%2.  Put pathway* to the destination vector
vector = {'ptypeiis.gb'};
%3. Put pathways* to inserts in the right order (1st insert on the 1st place, 
%the 2nd - on the 2nd place, etc.)
inserts = {'cmv.gb', 'gfp-corrected.gb'}; 
%4. Put spacers in the right order of junctions (1st: vector-insert1 junction, 
%2nd: insert1-insert2 junction, etc.); specify 'none' if there is no any
spacers = {' none','gccacc','none'}; 
%5. Press Run and Hope for the better!



%  *If a file in the working directory, put only the file's name


%% The code is here (+ there are some separate functions for convenience)
%       Processing inputs:
%Get sequences from files
[seq_vect] = get_data_from_file(vector);
[seq_inserts] = get_data_from_file(inserts);
%Delete occasional spaces in 'none'-specified spacers:
for ii = 1:length(spacers)
    spacers{ii}(isspace(spacers{ii})) = [];
end


%                   "Input Control"
%      (proceed assembly or break code with an error)
%Check if the vector has 2 recognition sites for TypeIIS enzyme
[f,r,s,n] = iis_enzymes(enzyme);
if length(regexpi(seq_vect{1}, f)) + length(regexpi(seq_vect{1}, r)) ~= 2
        error(['Vector must contain 2 binding sites for ', n])
end
%Check if the sites directed in the opposite sides (so that after 
%digestion no binding sites for the iis restrictase left)
if (length(regexpi(seq_vect{1}, f)) == 1 && length(regexpi(seq_vect{1}, r)) == 1) ~= 1
    error(['Wrong vector design: one binding site for ',n,' will still '...
            'be present in the vector after digestion'])
end
%Check if inserts have 0 or 2 recognition sites for TypeIIS enzyme
%and bifurcate processing of 0 and 2 sites containing inserts
just_assemble = logical.empty(length(inserts),0); 
for ii = 1: length(inserts)
    if (length(regexpi(seq_inserts{ii}, f)) + length(regexpi(seq_inserts{ii}, r))) == 2
        if (length(regexpi(seq_inserts{ii}, f)) == 1 && length(regexpi(seq_inserts{ii}, r))== 1) == 1
            just_assemble(ii) = true; %the further code will digest 
            %the insert and try to assemble if it is possible
        else
            error(['Wrong insert ',inserts{ii},' design: one binding site for ',n,' will still '...
            'be present in the insert after digestion']);
        end
    elseif length(regexpi(seq_inserts{ii}, f)) + length(regexpi(seq_inserts{ii}, r)) == 1
        error(['Insert ', inserts{ii}, ' has one ', n, ' recognition site']);
    elseif length(regexpi(seq_inserts{ii}, f)) + length(regexpi(seq_inserts{ii}, r)) == 0
        just_assemble(ii) = false; %the further code will create primers 
        %to extend inserts with the enzyme recognition sites
    else
        error(['Insert ', inserts{ii}, ' has 3 or more ', n, ' recognition sites']);
    end
end



%                              Digestion
%Note: all the sequences are written in 5'-3' direction 
% Vector digestion
cut_site_f = regexpi(seq_vect{1}, f) + length(f) + s; 
f_sticky = seq_vect{1}(cut_site_f:(cut_site_f+3)); 
cut_site_r = regexpi(seq_vect{1}, r) - s - 4; %4 is the length of sticky ends
r_sticky = seq_vect{1}(cut_site_r:(cut_site_r+3));

%Check the relative position of two binding sites in the overall seq-file
%(because matlab treats it like linear seq, not circular)
seq_after_r = extractAfter(seq_vect{1}, cut_site_r-1);    
if contains(seq_after_r, f) == 1 %the seq for cut out is fully inside the
    %overall sequence
    v_digested_R = eraseBetween(seq_vect{1},cut_site_r+5,cut_site_f-1); 
    %keep place at cut_site_r+4 for a unique character (start erasing
    %from cut_site_r+4(sticky_ends nucle)+1(reserved place))
    v_digested_R(cut_site_r+4) = 'R'; %"reserved" for insert (unique char
    %for a DNA seq - easy to manipulate)    
else
    %the seq for cut out starts in the end of overall sequence and continues in
    %the start
    v_digested_R = seq_vect{1};
    v_digested_R((cut_site_r+4):end) = [];
    v_digested_R(1:(cut_site_f-2)) = [];
    v_digested_R(1) = 'R';
end

%  Inserts digestion || primer design + "PCR" + digestion
i_digested = seq_inserts; 
[if_sticky, ir_sticky] = deal(cell(1, length(inserts)));
Primers = {'Primer', 'Seq', 'Length Whole',...
        'Length Anneal', 'Tm Anneal', 'GC% Anneal', 'Sticky ends'};
for ii = 1: length(inserts)
    %if inserts already have 2 binding sites for the restriction enzyme):
    if just_assemble(ii) == true 
        %Define existing sticky ends of the insert based on the cut site
        icut_site_f = regexpi(seq_inserts{ii}, f) + length(f) + s; 
        if_sticky{ii} = seq_inserts{ii}(icut_site_f:(icut_site_f+3)); 
        icut_site_r = regexpi(seq_inserts{ii}, r) - s - 4; 
        ir_sticky{ii} = seq_inserts{ii}(icut_site_r:(icut_site_r+3));      
        %Check complementarity of the sticky ends before digestion 
        %and break the code with an error if they are not complementary
        if ii == 1 %first insert after the vector
            if strcmpi(r_sticky, if_sticky{1}) == 0
                error(['The user-defined sticky ends in the junction of ',...
                    vector{1},' and ',inserts{1}, ' (the first junction)',...
                    ' are not complementary']);
            end
        
        elseif (ii > 1) && (ii < length(inserts)) %inserts from 2 to the penultimate
            %check comp with the previous insert 
            if strcmpi(ir_sticky{ii-1}, if_sticky{ii}) == 0
                error(['The user-defined sticky ends in the junction of ',...
                    inserts{ii-1},' and ',inserts{ii}, ' are not complementary']);
            end 
               
        elseif ii == length(inserts) %the last insert (before vector)
             if strcmpi(f_sticky, ir_sticky{end}) == 0
                 error(['The user-defined sticky ends in the junction of ',...
                    vector{1},' and ',inserts{end}, ' (the last junction)',...
                    'are not complementary']);
             end
             if strcmpi(ir_sticky{ii-1}, if_sticky{ii}) == 0
                error(['The user-defined sticky ends in the junction of ',...
                    inserts{ii-1},' and ',inserts{ii}, ' are not complementary']);
             end 
        end
        
        %Digest the insert:
        seq_after_r = extractAfter(seq_inserts{ii}, icut_site_r-1);
        if contains(seq_after_r, f) == 1 
            i_digested{ii} = eraseBetween(seq_inserts{ii},icut_site_r+5,icut_site_f-1); 
            i_digested{ii}(icut_site_r+4) = 'V'; %unique symbol for DNA > easy
            %to manipulate
            i_digested{ii} = circshift(i_digested{ii}, -(strfind(i_digested{ii}, 'V')-1));  
            %to bring to the convenient look: V-r(ii)-insert-f(ii)
        else
            i_digested{ii}(icut_site_r+4:end) = [];
            i_digested{ii}(1:icut_site_f-2) = [];
            i_digested{ii}(1) = 'V';
        end
        
    %Now let's process inserts which do not have binding sites yet
    PCR = cell(length(inserts),2); %memory preal
    elseif just_assemble(ii) == false 
        %Create forward primer: coding region
        Tm = 0; nn = 10; %nn is starting number of nucleotides in the primer
        while Tm < 60 %prolong primer until its Tm >= 60
            if nn == length(seq_inserts{ii})
                error(['The insert ', inserts{ii}, ' is too short, or '...
                    'GC-poor, or both. I could not create primers with '...
                    'Tm > 60C for it. Try putting this insert as a '...
                    'spacer, or extend ',inserts{ii},...
                    ' insert.']);
            end
            F_primer = seq_inserts{ii}(1:nn);
            F_properties = oligoprop(F_primer);
            Tm = F_properties.Tm(5); %Nearest-neighbor (SantaLucia Jr., 1998)
            nn = nn + 1; %add one more nucleotide
        end
        %give name to the primer
        Primers{size(Primers,1)+1,1} = [erase(inserts{ii}, {'.ape', ...
            '.gb', '.txt'}), '-Forward']; %'size' changes every loop
        %Record: coding region of the primer, anneal seq length, Tm,
        %GC-content of the annealing part
        Primers{size(Primers,1),2} = upper(F_primer);
        Primers{size(Primers,1),4} = nn;
        Primers{size(Primers,1),5} = Tm;
        Primers{size(Primers,1),6} = F_properties.GC;
        
        %Forward primer: overhang
          %   Add spacer (before this insert) if there is any:
        if strcmpi(spacers{ii}, 'none') == 0
            %Check if the entered spacer consists only of nucleotides:
            for hh = 1:length(spacers{ii}) %compare each symbol to 4 allowed nucleotides
                xx = strcmpi(spacers{ii}(hh), {'g','t','c','a'}); %it will give 1x4 logical array
                if sum(xx) == 0 %if the symbol is different from all 4 nucleotides
                    error(['The spacer ', spacers{ii}, ' must include only'...
                        ' DNA nucleotides (C, G, T, A)']);
                end
            end
            %If the input is okay, add it before the coding region:
            Primers{size(Primers,1),2} = [lower(spacers{ii}), Primers{size(Primers,1),2}];
        end
          %   Add sticky ends (which are determined by the end of the coding
          %region of the previous part = previous part's sticky ends)
          if ii == 1 %first insert after the vector
              Primers{size(Primers,1),2} = [lower(r_sticky), Primers{size(Primers,1),2}];
              Primers{size(Primers,1),7} = upper(r_sticky); %record sticky ends
              if_sticky{ii} = r_sticky; %record the variable for the proper
              %working of the previous loop
          else %inserts from 2 to the last
              Primers{size(Primers,1),2} = [lower(ir_sticky{ii-1}), Primers{size(Primers,1),2}];
              Primers{size(Primers,1),7} = upper(ir_sticky{ii-1});
              if_sticky{ii} = ir_sticky{ii-1};
          end
          %Record the initial seq for a quick "PCR+digest@ (next part):
          PCR{ii,1} = Primers{size(Primers,1),2}; 
          PCR{ii,1} = PCR{ii,1}(1:(regexpi(PCR{ii,1}, F_primer)-1));
          %PCR array: col 1 – forward part, col 2 – reverse part which is 
          %already rev comp and ready to be concatenate with the coding seq
          %Add 1-2 (depending on the enzyme's cutting shift, s) random
          %nucleotides:
          Primers{size(Primers,1),2} = [lower(randseq(s)), Primers{size(Primers,1),2}];
          %Add recognition site:
          Primers{size(Primers,1),2} = [lower(f), Primers{size(Primers,1),2}];
          %Add 6 extra nucleotides, so the enzyme could physically bind:
          Primers{size(Primers,1),2} = [lower(randseq(6)), Primers{size(Primers,1),2}];
          %Count the whole length:
          Primers{size(Primers,1),3} = length(Primers{size(Primers,1),2});
            
          %Reverse primer design:
          %Coding region
           Tm = 0; nn = 10; 
           while Tm < 60 %prolong primer until its Tm >= 60
               R_primer = seq_inserts{ii}(end+1-nn:end);
               R_primer = seqrcomplement(R_primer); %into 5'-3' seq-view
               R_properties = oligoprop(R_primer);
               Tm = R_properties.Tm(5); %Nearest-neighbor (SantaLucia Jr., 1998)
               nn = nn + 1; %add one more nucleotide-next loop
           end
           Primers{size(Primers,1)+1,1} = [erase(inserts{ii}, {'.ape', ... 
               '.gb', '.txt'}), '-Reverse']; %name
           Primers{size(Primers,1),2} = upper(R_primer); %seq
           Primers{size(Primers,1),4} = nn; %length
           Primers{size(Primers,1),5} = Tm;
           Primers{size(Primers,1),6} = R_properties.GC; %GC-content
           if ii == length(inserts) %(last insert) add sticky ends from the vector
               Primers{size(Primers,1),2} = [lower(seqrcomplement(f_sticky)), Primers{size(Primers,1),2}];
               Primers{size(Primers,1),7} = f_sticky;
               ir_sticky{ii} = f_sticky;
           else
           Primers{size(Primers,1),7} = upper(seqrcomplement(R_primer(1:4))); %sticky ends
           ir_sticky{ii} = Primers{size(Primers,1),7};
           end
          %Record the initial seq for a quick "PCR+digest" (next part):
          PCR{ii,2} = seqrcomplement(Primers{size(Primers,1),2}); 
          PCR{ii,2} = PCR{ii,2}(((regexpi(PCR{ii,2}, seqrcomplement(R_primer), 'end')+1):end));
          %PCR array: col 1 – forward part, col 2 – reverse part which is 
          %already rev comp and ready to be concatenate with the coding seq
          
          %Reverse primer: overhang
          %Add 1-2 (depending on the enzyme's cutting shift, s) random
          %nucleotides:
          Primers{size(Primers,1),2} = [lower(randseq(s)), Primers{size(Primers,1),2}];
          %Add recognition site:
          Primers{size(Primers,1),2} = [lower(f), Primers{size(Primers,1),2}];
          %Add 6 extra nucleotides, so the enzyme could physically bind:
          Primers{size(Primers,1),2} = [lower(randseq(6)), Primers{size(Primers,1),2}];
          %Count the whole length:
          Primers{size(Primers,1),3} = length(Primers{size(Primers,1),2});
          
    end
end

%        "Amplicons' digest" (my 'light' version with PCR variables)
for jj = 1:size(PCR,1)
    i_digested{1,jj} = ['V' PCR{jj,1} i_digested{jj} PCR{jj,2}];
    %V is a unique char for DNA seq, easy to manipulate in ligation stage
end  

%                           Ligation
%All the digested inserts have been previously brought to the format:
%V-r(i)-coding_seq-f(i) where is V is a special char, r and f - sticky ends
%(up and down in 5'-3' duplex, respectively); the vector has the format:
%seq-R-seq, where R is the site where we want to put all inserts.
%I am going to 1)ligate all inserts together and 2) put them instead of R.
for ii = 1:length(inserts) % (1)
i_digested{1,ii} = i_digested{1,ii}(6:end); %erase V-r(i)
if ii == length(inserts) %last
    i_digested{1,ii} = i_digested{1,ii}(1:(end-4)); %erase f(i) for the last
end
end
big_insert = strcat(i_digested{:}); %concat all inserts
assembly_seq = [v_digested_R(1:strfind(v_digested_R, 'R')) big_insert ...
    v_digested_R(strfind(v_digested_R, 'R'):end)]; % step (2)
assembly_seq = erase(assembly_seq, 'R');

%Record all non-complementary sticky-ends in GG reaction:
sticky_array = [if_sticky(:)' f_sticky];

%Create assembly name based on initial file names
name = strcat('Assembly_', erase(vector{1}, {'.ape', '.gb', '.txt'}));
for ii = 1:length(inserts)
inserts{ii} = erase(inserts{ii}, {'.ape', '.gb', '.txt'});
name = strcat(name, '–', inserts{ii});
end

clearvars -except name assembly_seq Primers sticky_array

%Check if sticky ends of different junctions can be used in one-pot GG
%reaction (for more information, please address check_sticky.m function)
[compatability, non_comp_sticky] = check_sticky(sticky_array, 0);
if compatability == 1
    fprintf('Sticky ends in this assembly are compatible:\n');
    disp(sticky_array);
    fprintf('\n');
else
    fprintf('Sticky ends in this assembly are NOT compatible:\n');
    disp(sticky_array);
    fprintf('Try to decrease the threshold for fidelity. \n\n');
end


%Export 
writecell({assembly_seq}, name);
writecell(Primers,'Primers.xls');

fprintf('Please find primers for PCR in Primers.xls and the assembled sequence in the file\n%s.txt.\n\n', name);

restriction_analysis(assembly_seq, 400); %400 is a thresh for light bands, 
%please address the function for more details




