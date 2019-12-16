%Run this script if you just need a test digest for you plasmid (not GG
%assembly). 
%Steps: 
%1. Put your .ape or .gb file in the working directory
%2. Put its name in the brackets below (instead of the sample 
%pROC078 +insert_Y+insert_pA.gb name)
%3. Press 'Run'

seq = get_data_from_file({'pROC078 +insert_Y+insert_pA.gb'}); 
seq = seq{1};
restriction_analysis(seq, 400) %put threshold for light bands, e.g. 400 
%stands for no bands < 400 bp are allowed


%If no options found. Try to extend the enzyme database
%(add some enzymes to Test_digest_enzymes.xlsx). Use this web-site for info: 
%https://international.neb.com/tools-and-resources/usage-guidelines/nebuffer-performance-chart-with-restriction-enzymes

