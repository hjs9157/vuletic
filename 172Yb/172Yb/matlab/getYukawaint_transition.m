function [dY,varargout] = getYukawaint_transition(stateinfo1,stateinfo2,DHForCI,m)
%% calculate
fprintf('%s:\n',stateinfo1{1})
[Y1,Y_core1,Y_pop1,Y_corr1] = getYukawaInt_fermions(stateinfo1{:},DHForCI,m);
fprintf('%s:\n',stateinfo2{1})
[Y2,Y_core2,Y_pop2,Y_corr2] = getYukawaInt_fermions(stateinfo2{:},DHForCI,m);

%% output
dY = Y1 - Y2;
varargout{1} = Y_core1 - Y_core2; % dY_core
varargout{2} = Y_pop1 - Y_pop2; % dY_pop
varargout{3} = Y_corr1 - Y_corr2; % dY_corr
end