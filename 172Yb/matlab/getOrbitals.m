function T_orbital = getOrbitals(stateID,DHForCI)
% stateID = '2S12';

% load summary data
S = readSum(stateID,DHForCI);
sumtable1 = S.sumtable1;
sumtable2 = S.sumtable2;
T_orbital = [sumtable1,sumtable2];

% Load radial wavefunctions
subshell = sumtable1.Properties.RowNames;
% getRwfnDat(stateID)
% rwfn = cell(length(subshell),1);
% for shinx = 1:length(subshell)
%     rwfn_dat = dlmread(sprintf('dat/rwfn_%s_%s.dat',subshell{shinx},stateID));
%     T_rwfn = array2table(rwfn_dat);
%     T_rwfn.Properties.VariableNames = {'r', 'P', 'Q'};
%     T_rwfn.Properties.VariableUnits = {'a_0','',''};
%     
%     rwfn{shinx} = T_rwfn;
% end
% T_orbital.rwfn = rwfn;
% T_orbital = T_orbital(:,[end,1:end-1]);

T_rwfn = getRwfn(stateID);
% T_orbital = join(T_rwfn,T_orbital);
T_orbital = [T_rwfn(:,1:5),T_orbital];


% T_orbital.Properties.Description = stateID; % put stateID in table description

fname_save = ['dat/T_orbital_',stateID,'_',DHForCI,'.mat'];
save(fname_save,'stateID','T_orbital')
end