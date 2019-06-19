function [mixc,csf] = readRmix(stateID,DHForCI,varargin)

% stateID = '1S0';

P_cutoff = 0; % omit negligible CSFs once papulation is accumulated enough 
if length(varargin) >= 1, P_cutoff = varargin{1}; end

%% load every line in file

direc_load = '../';
if strcmp(DHForCI,'DHF')
    fname_load = [stateID,'.rmix'];
elseif strcmp(DHForCI,'CI')
    fname_load = [stateID,'.crmix'];
else
    error('Wrong input for DHForCI.')
end

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
if ~strcmp(fstr{1},' RMIXEXTRACT')
    error('wrong .rmix file.')
end

%% read CSFs and mixing coeffs
for li = 1:n_line
    if length(fstr{li}) >= 7
        if strcmp(fstr{li}(2:7),'Energy')
            li_start = li + 2;
            break;
        end
    end
end

expr_lterm = '[spd-z]';
    % expression = ['(\d+)(',exp_orbital,')\((\d+)\)'];
%     expression = ['(?<n>\d+)(?<l>',expr_lterm,')\((?<occupn>\d+)\)'];
%     expr_csf = ['(?<n>\d+)(?<l>',expr_lterm,'-?)\s?\(\s?(?<occupn>\d+)\)'];
expr_csf = ['(?<subshell>\d+',expr_lterm,'-?)\s?\(\s?(?<occupn>\d+)\)'];

datastr = fstr(li_start:n_line-1);
n_csf = (length(datastr))/4;
mixc = nan(n_csf,1); % mixing coefficient
csf = cell(n_csf,1); % csfs
if round(n_csf) ~= n_csf, error('wrong number of lines in datastr.'); end
P_remain = 1; % 1 - accumulated papulation (assuming CSFs are in order of population)
% P_cutoff = 0; % omit negligible CSFs once papulation is accumulated enough 
for csi = 1:n_csf % run over csfs
    % get mixing coefficient
    mixc(csi) = str2double(datastr{(csi-1)*4+1}(15:23));
    P_remain = P_remain - mixc(csi)^2;
    % get list of subshells in csf and their occupations
    S_csf = regexp(datastr{(csi-1)*4+2},expr_csf,'names');
    cell_subshell = {S_csf.subshell};
    cell_occupn = num2cell(str2double({S_csf.occupn}));
    csf{csi} = [cell_subshell; cell_occupn];
    
    % cutoff negeligible CSFs
%     if P_remain < P_cutoff, break; end
end
n_omitted = n_csf - csi;
n_csf = n_csf - n_omitted;
mixc = mixc(1:n_csf);
csf = csf(1:n_csf);

end