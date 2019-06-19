function S = readSum(stateID,DHForCI)
DHForCI = 'DHF'; % no use of CI for now
% fname_load = '2S12.sum';
% fname_save = '2S12.sum.mat';
if strcmp(DHForCI,'DHF')
    fname_load = ['../',stateID,'.sum'];
    fname_save = ['dat/',stateID,'.sum.mat'];
elseif strcmp(DHForCI,'CI')
    fname_load = ['../',stateID,'.csum'];
    fname_save = ['dat/',stateID,'.csum.mat'];
else
    error('Wrong input for DHForCI.')
end

%% load every line in file
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

%% find first table
li_header = NaN;
for li = 1:n_line
    str_find = 'Subshell';
    if length(fstr{li}) >= length(str_find)
        if strcmp(fstr{li}(1:length(str_find)),str_find)
            li_header = li;
            break;
        end
    end
end
li_start = li_header + 2;
for li = li_start:n_line
    if isempty(fstr{li})
        li_end = li - 1;
        break;
    end
end

% get table headers
% header_dlm = split(fstr{li_header},' ');
% header = header_dlm(~cellfun('isempty',header_dlm));
if strcmp(DHForCI,'DHF')
    header = {'Subshell','e','p0','gamma','P2','Q2','SelfConsistency','MTP'};
else
    header = {'Subshell','e','p0','gamma','P2','Q2','MTP'};
end

% get orbital ID's and values
datacell = cell(li_end - li_start + 1,length(header));
i = 0;
for li = li_start:li_end
    i = i+1;
    str_dlm = split(fstr{li},' ');
    datacell(i,:) = str_dlm(~cellfun('isempty',str_dlm));
end
subshell = datacell(:,1);
valuestr = cellfun(@(str) replace(str,'D','E'),datacell(:,2:end),'UniformOutput',false);
value = str2double(valuestr);

% construct Table variable
sumtable1 = array2table(value);
sumtable1.Properties.VariableNames = header(2:end);
sumtable1.Properties.RowNames = subshell;

%% find second table 
li_header = NaN;
for li = li_end:n_line
    str_find = 'Subshell';
    if length(fstr{li}) >= length(str_find)
        if strcmp(fstr{li}(1:length(str_find)),str_find)
            li_header = li;
            break;
        end
    end
end
li_start = li_header + 2;
for li = li_start:n_line
    if isempty(fstr{li})
        li_end = li - 1;
        break;
    end
end

% get table headers
header = {'Subshell','mrm3','mrm1','mr1','mr2','mr4','occupation'};

% get orbital ID's and values
datacell = cell(li_end - li_start + 1,length(header));
i = 0;
for li = li_start:li_end
    i = i+1;
    str_dlm = split(fstr{li},' ');
    datacell(i,:) = str_dlm(~cellfun('isempty',str_dlm));
end
% subshell = datacell(:,1);
valuestr = cellfun(@(str) replace(str,'D','E'),datacell(:,2:end),'UniformOutput',false);
value = str2double(valuestr);

% construct Table variable
sumtable2 = array2table(value);
sumtable2.Properties.VariableNames = header(2:end);
sumtable2.Properties.RowNames = subshell;

%% return and save
save(fname_save,'fstr','sumtable1','sumtable2');
S = load(fname_save);
end