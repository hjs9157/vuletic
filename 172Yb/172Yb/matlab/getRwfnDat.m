function getRwfnDat(stateID)

fname_load = ['../octave_',stateID,'.m'];
% fname_save = ['orbitalda',stateID,'.dat'];

%% load every line in octave_.m file
fstr = {};
i = 0;
fid = fopen(fname_load);
while ~feof(fid)
    i = i+1;
    fstr{i} = fgetl(fid);
%     disp(fstr{i})
end
fclose(fid);

n_line = length(fstr);

%% extract name of orbitals
str = fstr{end}; % 'legend(' 1s  ',....
% temp = regexp(str,''' [^,]+  ''','match');
% orbitals = cell(1,length(temp));
orbitals = regexp(str,'(?<=''\s*)[^,\s]+(?=\s*'')','match','all');
fprintf('orbitals: ')
for i = 1:length(orbitals)
    fprintf('%s, ',orbitals{i})
end
fprintf('\b\b\n')

%% save each wavefunction
li = 0; % index of line
oi = 0;% index of orbital
while li < n_line
    li = li+1;
    
    if strcmp(fstr{li},' P = [')
        oi = oi+1;
        fprintf(orbitals{oi})
        fname_save = ['dat/rwfn_',orbitals{oi},'_',stateID,'.dat'];
        fid = fopen(fname_save,'w');
        li = li+1;
        while ~strcmp(fstr{li},' ];')
            fprintf(fid,'%s\n',fstr{li}(5:end));
            li = li+1;
        end
        fclose(fid);
        fprintf(' saved.\n')
    end
end