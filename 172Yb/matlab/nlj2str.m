function str = nlj2str(nlj)
n = nlj(1); l = nlj(2); j = nlj(3);
lterms = 'spdfghijklmnoqrtuvwxyz';
lterm = lterms(l+1);
if j - l < 0 
    jterm = '-';
else
    jterm = '';
end
str = sprintf('%u%s%s',n,lterm,jterm);
