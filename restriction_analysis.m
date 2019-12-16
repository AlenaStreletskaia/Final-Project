%The function outputs two optimal typeIIP enzymes for restriction analysis
%and bands (size, bp) that are supposed to be seen after the test digest
%Inputs:
%1) plasmid sequence as a char array;
%2) thres is a threshold to get rid of 'light' bands (e.g. bands < 400 bp 
%are poorly visible in 1.2% agarose gel after the mid-long run > we 
%want bands > 400 bp, therefore: thres = 400)

function restriction_analysis(seq, thres)
tde = readtable('Test_digest_enzymes.xlsx','basic',true); %'basic' 
%(non-formatted) reading saves time

%Detect if there are cutting sites for these enzymes on the seq and record
%if any
ind = false(size(tde.Enzyme,1),1); f = 0; 
for ii = 1:size(tde.Enzyme,1)
    site = erase(tde.Seq{ii}, '/'); %binding seq 5'-3'
    if isempty(regexpi(seq, site)) == 0 %since seq are palindromic – no need
        %to search for the rev compl sites (type IIP enzymes)
        ind(ii) = true; %replace zeros with ones if there is a binding site(s)
        f = f + 1;
        cut_sites{f} = {regexpi(seq, site)+(regexpi(tde.Seq{ii}, '/')-2)};
    end
end
if sum(ind, 'all') == 0
    error('No binding sites found. Extend the library for TypeIIP enzymes.')
end
%Proceed only with enzymes that can bind to the seq:
tde = tde(ind,:); 
%Add cut sites to the table for convenience:
tde = [tde cut_sites']; tde.Properties.VariableNames{3} = 'Cut_sites';
%'Digest' in pairs:
options = {};
for jj = 1:(size(tde,1)-1)
    all_cut_sites = sort([tde.Cut_sites{jj} tde.Cut_sites{jj+1}]);
    %Get bands between cut sites:
    bands = zeros(1, length(all_cut_sites));
    for nn = 1:(length(bands)-1)
        bands(nn) = all_cut_sites(nn+1) - all_cut_sites(nn);
    end
    %Get also a band between the last and first cut sites:
    bands(end) = length(seq) - max(all_cut_sites) + min(all_cut_sites);
    bands = sort(bands, 'descend');
    %Now check if bands < threshold are present:
    if isempty(bands(bands < thres)) == 1 %continue only with enzymes that
    %give heavier bands
        %Apply a "distance-optimum" eq:
        %(x is the distances between closest bands)
        func = @(x) -(20/(0.45*x + 0.75)^8 - 20/(0.8*x+0.5)^2);
        opt = zeros(1,length(bands)-1); 
        for i = 1:length(opt)
            x = (bands(i) - bands(i+1))/1000; %distance into kb
            if func(x) > 0.94 %chosen threshold for y(x) value
                opt(i) = 1;
            end
        end
        if sum(ismember(opt, 0), 'all') == 0 %check if there are no zeros 
        %left (optimal separation)
        options((size(options,1)+1), 1:3) = [tde.Enzyme(jj) tde.Enzyme(jj+1) {bands}]; 
        end
    end
end
    for w = 1:(size(options,1))
        if size(options,1) == 1
            fprintf('Only one test-digest option is found.\n');
        else
        fprintf('Option %d\n', w);
        end
        fprintf('Use %s and %s to get bands:\n', options{w,1}, options{w,2});
        for ww = 1:length(options{w,3})
            fprintf('%d\n', options{w,3}(ww));
        end
    end
    if isempty(options) == 1
        fprintf('No test digest options found. Extend enzyme database. \n');
    end
end