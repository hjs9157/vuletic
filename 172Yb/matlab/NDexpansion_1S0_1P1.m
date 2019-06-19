numberDensity_fermions_1S0_1P1;
drawnow;
ND = dND;

%%

% nuclear size
r_n_fm = 5.3; % fm; nuclear size 
a_0 = 52.917721067; % pm; Bohr radius = 5.2917721067(12)×10−11 m (2014 CODATA)
r_n = r_n_fm/a_0/1e3;

% Fermi distribution parameters
c_fm = 6.227;
a_fm = 2.18/4/log(3); % for 172Yb
c = c_fm/a_0/1e3; a = a_fm/a_0/1e3; 

order_fit = 3;
InxCap_fit = find(R > r_n*1e-2 & R < r_n*1e1);
InxBottom_fit = 10;

[S,NDexp_fig] = getNDexpansion(R,ND,order_fit,r_n,InxCap_fit,InxBottom_fit);