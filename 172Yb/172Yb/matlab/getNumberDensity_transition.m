function [r,dND,varargout] = getNumberDensity_transition(stateinfo1,stateinfo2,DHForCI)
%% calculate
fprintf('%s:\n',stateinfo1{1})
[r,ND1,ND_core1,ND_pop1,ND_corr1] = getNumberDensity_fermions(stateinfo1{:},DHForCI);
fprintf('%s:\n',stateinfo2{1})
[r,ND2,ND_core2,ND_pop2,ND_corr2] = getNumberDensity_fermions(stateinfo2{:},DHForCI);

%% output
dND = ND1 - ND2;
varargout{1} = ND_core1 - ND_core2; % dY_core
varargout{2} = ND_pop1 - ND_pop2; % dY_pop
varargout{3} = ND_corr1 - ND_corr2; % dY_corr
end