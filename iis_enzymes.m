%Get the recognition site and cut shift of a type IIs restrictase 
%(this mini-database can be easily extended if needed)
%The information on the typeIIs restriction enzyme is on the NEB site.

function [rec_site_f, rec_site_r, shift, name] = iis_enzymes(enzyme)
if enzyme == 1 
    rec_site_f = 'GGTCTC'; rec_site_r = seqrcomplement(rec_site_f);
    shift = 1; % #nucleotides between cutting sites and sticky ends
    name = 'BsaI';
elseif enzyme == 2
    rec_site_f = 'GAAGAC'; rec_site_r = seqrcomplement(rec_site_f); 
    shift = 2;
    name = 'BbsI';
elseif enzyme == 3
    rec_site_f = 'CGTCTC'; rec_site_r = seqrcomplement(rec_site_f);
    shift = 1;
    name = 'Esp3I';
end
    