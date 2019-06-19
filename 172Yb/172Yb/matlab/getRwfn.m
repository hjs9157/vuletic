function T_orbital = getRwfn(stateID)
% stateID = '1S0';

direc_load = '../';
fname_load = [stateID,'.w.readrwf'];

%% load every line in file
fstr = {};
i = 0;
fid = fopen([direc_load,fname_load]);
while ~feof(fid)
    i = i+1;
    fstr{i} = fgetl(fid);
%     disp(fstr{i})
end
fclose(fid);
fstr = fstr.';

n_line = length(fstr);

% check if it is right file
if ~strcmp(fstr{1},'G92RWF')
    error('wrong .w.readrwf file.')
end

%% read data
% each header corresponds to a subshell
expr_header = '\s+(?<n>\d+)\s+(?<k>-?\d+)\s+(?<e>0\.\d+D(+|-)\d\d)\s+(?<MTP>\d+)';
% extract header only
temp = regexp(fstr,expr_header,'names');
lInx_header = find(cellfun(@(c) ~isempty(c),temp));
header = temp(lInx_header);
n_subshell = length(header);

% construct table
% varNames = {'n','l','j','k','rwfn','e','p0','MTP'};
% varTypes = {'double','double','double','double','table','double','double','double'};
% T_orbital = table('Size',[n_subshell,length(varNames)],...
%     'VariableNames',varNames,...
%     'VariableTypes',varTypes);
T_orbital_S = struct;
subshell = cell(n_subshell,1);
lterms={'s','p','d','f','g','h','i','j','k'};
for shi = 1:n_subshell % run over subshells
    % get info from header
    n = str2double(header{shi}.n);
    k = str2double(header{shi}.k);
    e = str2double(replace(header{shi}.e,'D','E'));
    MTP = str2double(header{shi}.MTP);
    
    % read after header (p0 & wavefunction)
    li_start = lInx_header(shi)+1;
    if shi == length(header)
        li_end = n_line;
    else
        li_end = lInx_header(shi+1) - 1;
    end
    
    temp = fstr(li_start:li_end);
    temp = replace(temp,'D','E');
    temp = replace(temp,' ','');
    
    p0 = str2double(temp{1});
    
    P = str2double(temp(2:(MTP+1)));
    Q = str2double(temp((MTP+2):(MTP*2+1)));
    r = str2double(temp((MTP*2+2):(MTP*3+1)));
    
    rwfn = table(r,P,Q); % r: Bohr radius
    
    % for saving
    Pstr = temp(2:(MTP+1));
    Qstr = temp((MTP+2):(MTP*2+1));
    rstr = temp((MTP*2+2):(MTP*3+1));
    
    % get state
    j = abs(k) - 1/2;
    l = j + sign(k)/2;
    subshell{shi} = sprintf('%u%s',n,lterms{l+1});
    if sign(k) > 0
        subshell{shi} = [subshell{shi},'-'];
    end
     
    % add line to table
    T_orbital_S(shi).n = n;
    T_orbital_S(shi).l = l;
    T_orbital_S(shi).j = j;
    T_orbital_S(shi).k = k;
    T_orbital_S(shi).rwfn = rwfn;
    T_orbital_S(shi).e = e;
    T_orbital_S(shi).p0 = p0;
    T_orbital_S(shi).MTP = MTP;
    
    % save to .dat file (as output of 'plotmcdf' in Grasp2k)
    direc_save = 'dat/';
    fname_save = ['rwfn_',stateID,'_',subshell{shi},'.dat'];
    fid = fopen([direc_save,fname_save],'w');
    for li = 1:MTP
        if ~strcmp(rstr{li}(1),'-'), rhead = ' '; else,  rhead = ''; end
        if ~strcmp(Pstr{li}(1),'-'), Phead = ' '; else,  Phead = ''; end
        if ~strcmp(Qstr{li}(1),'-'), Qhead = ' '; else,  Qhead = ''; end
        fprintf(fid,'%s%s   %s%s   %s%s\n',...
            rhead,rstr{li},Phead,Pstr{li},Qhead,Qstr{li});
    end
    fclose(fid);
end

T_orbital = struct2table(T_orbital_S);
T_orbital.Properties.RowNames = subshell;

end