function nlj = str2nlj(str)
% str = '4f-';
expr_lterm = '[spd-z]';
expr = ['(\d+)(',expr_lterm,')(-?)'];
terms = regexp(str,expr,'tokens','once');

n = str2double(terms{1});
lterms = 'spdfghijklmnoqrtuvwxyz';
l = strfind(lterms,terms{2}) - 1;
jminusl = isempty(terms{3}) - 1/2;
j = jminusl + l;
nlj = [n,l,j];
end