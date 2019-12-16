%The function checks compatability of the sticky ends in the assembly based
%on the screening data in the table FileS03_T4_18h_25C.
%It outputs compatability = 1 if compatible and compatability = 0 if some
%of the input sticky ends are not compatible; in the last case it also
%outputs which sticky ends are not compatible – pairs go in raws.

%Inputs: 
%1) sticky ends to check (as a cell array): 5'-3' of one chain. E.g.:
% 5'-TACA-3'         5'-GCCA-3'
%    ||||     and       ||||
% 3'-ATGT-5'         3'-CGGT-5'
%should be written as: {'TACA','GCCA'} (which is standart designation)
%2) threshold which is the number of non-specific ligation events: 
%0 is for the highest fidelity, whereas e.g. specific ligation gives 
%3-5k ligation events

%Example: how the function works.
%sticky_array = {'TACA', 'CTAA', 'GGAA', 'GCCA'};
%highest_thres = 0;
%[compatability, non_comp_sticky] = check_sticky(sticky_array, highest_thres);
% The function will output:
%compatability =
%     0
%non_comp_sticky =
%  1?2 cell array
%    {'TACA'}    {'GCCA'}
%Why did we get such a result?
%Because 
% 5'-GCCA-3' 
%      ||
% 3'-ATGT-5'
%Will give us ~0.1% non-specific ligation events (Please address the Reference
%table, so you could find 1 ligation event in the intersection of GCCA and
%TGTA (5'>3') overhangs; 1 ~= highest_thres). If such accuracy is
%redundant, decrease threshold to 1 (or any value, but I would recommend 
%the range [0; 50] for aprx 0-0.5% non-specific ligation events):
%thres = 1;
%[compatability, non_comp_sticky] = check_sticky(sticky_array, thres)
%compatability =
%     1
%non_comp_sticky =
%  0?0 empty cell array

%Reference for the table FileS03_T4_18h_25C (attached in the 'Reference' folder):
%Potapov V, Ong JL, Kucera RB, et al. Comprehensive Profiling of Four Base Overhang
%Ligation Fidelity by T4 DNA Ligase and Application to DNA Assembly. ACS Synthetic
%Biology. 2018;7(11):2665-2674. (Supplementary materials > sb8b00333_si_002.zip
%> FileS03_T4_18h_25C.xlsx)
%In brief,the table contains information on the specificity of the sticky ends
%ligation, practically obtained with T4 ligase at 25C


function [compatability, non_comp_sticky] = check_sticky(sticky_array, threshold)
T = load('FileS03_T4_18h_25C.mat'); T = T.T; %load processed to a cell array
%(to have doubles and rows in dif dimensions) table
comp_mat = zeros(length(sticky_array),length(sticky_array))-1; %memory preallocation

%  Create pairwise comparison matrix:
for ii = 1:length(sticky_array) %ii - columns of the comp matrix
    seq1 = sticky_array{ii};
    for jj = 1:length(sticky_array) %jj -rows of the comp matrix
        seq2 = seqrcomplement(sticky_array{jj}); %rev compl table view
        %Find the number of ligation events in FileS03_T4_18h_25C.mat:
        a = strcmpi(T(1,:), seq1); %get logical array where the seq1 (col)
        b = strcmpi(T(:,1), seq2); %get logical array where the seq2 (row)
        c = cell2mat(T(a,b)); %get intersection (normalized number of cases when ligation occurs)
        if c <= threshold %put threshold (== 0 for the highest fidelity) 
            comp_mat(jj,ii) = 1; %sticky ends are compatible for GG
        else
            comp_mat(jj,ii) = 0; %sticky ends cann't be used in one-pot reaction
        end
    end
end
comp_mat(logical(eye(size(comp_mat,1),size(comp_mat,2)))) = 1; %change
%diagonal values (derived from one and the same sticky-ends pair) to 1
if sum(ismember(comp_mat, 0), 'all') ~= 0 %check if there are zeros
    compatability = 0; %non-compatible
    [row, col] = find(~comp_mat); %index zeros
    non_comp_sticky = cell(length(row),2); %table of non-comp sticky ends in pairs
    for mm = 1:length(row)
        non_comp_sticky{mm, 1} = sticky_array{col(mm)};
        non_comp_sticky{mm, 2} = sticky_array{row(mm)};%back to initial 5'-3'
    end
    %Remove repeated (inversed pairs like: TTTT + AATA and AATA + TTTT) 
    ll = 1; 
    while size(non_comp_sticky, 1) > nchoosek(length(sticky_array),2) %C(n,k)
        line = cell(1,2); line(1,1) = non_comp_sticky(ll,2);
        line(1,2) = non_comp_sticky(ll,1); %switch columns
        line = strcat(line(1,1), line(1,2)); %cat for comparison
        array = strcat(non_comp_sticky(:,1), non_comp_sticky(:,2));
        indices = strfind(array, line);
        indices(cellfun('isempty',indices)) = {0}; %remove empty cells
        indices = cell2mat(indices);
        if sum(indices, 'all') == 1 %break the loop if no more repeats
            break
        end
        indices(find(indices,1,'first')) = 0; %keep first encountered pair
        non_comp_sticky(logical(indices),:) = []; %remove the repeated one
        ll = ll + 1;
    end 
else
    compatability = 1;
    non_comp_sticky = {};
end











