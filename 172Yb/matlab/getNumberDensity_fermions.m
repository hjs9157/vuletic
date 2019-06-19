function [r,ND,varargout] = getNumberDensity_fermions(stateID,confstr_core,confstr_valence,DHForCI)

% % stateID = '1S0';
% DHForCI = 'DHF';

% % electron configuration of state
% confstr_core = ['1s(2)2s(2)2p-(2)2p(4)',...
%     '3s(2)3p-(2)3p(4)3d-(4)3d(6)',...
%     '4s(2)4p-(2)4p(4)4d-(4)4d(6)',...
%     '5s(2)5p-(2)5p(4)4f-(6)4f(8)'];
% confstr_valence = '6s(2)';

% grid
NNNP = 1990;
RNT = 2.857142857143e-8;
H = 0.05;
f_r = @(I) RNT*(exp((I-1)*H)-1);
r = f_r((1:NNNP).');

% load .w.readrwf and .sum file of subshells
T_orbital = getOrbitals(stateID,DHForCI);
fprintf('   - Wavefunctions loaded.\n')

% load ASF (mixing coeff & csfs) from .rmix (or .crmix) file
[mixc,csf] = readRmix(stateID,DHForCI); % m: m_j=-j,..,j
n_csf = length(mixc);
fprintf('   - ASF loaded; n_csf=%u\n',n_csf)

%% Core
fprintf('   - Core calculation started.\n')
% get list of subshells and occupations from str
expr_lterm = '[spd-z]';
expr_conf = ['(?<subshell>\d+',expr_lterm,'-?)\((?<occupn>\d+)\)'];
conf_core = regexp(confstr_core,expr_conf,'tokens');
conf_core = vertcat(conf_core{:});
conf_core(:,2) = cellfun(@(c) str2double(c),conf_core(:,2),'UniformOutput',false);
N_core = sum([conf_core{:,2}]); % number of core electrons
n_subshell_core = size(conf_core,1);
fprintf('\t- N_core=%u, n_subshell=%u\n',N_core,n_subshell_core)
% calculate Yukawa integration
ND_core_subshell = cell(n_subshell_core,1);
ND_core = zeros(NNNP,1);
fprintf('\t- Integration for core...')
for shi=1:n_subshell_core
    subshell = conf_core{shi,1};
    occupn = conf_core{shi,2};
    rwfn = T_orbital.rwfn{subshell};
    
    % get number density as a function of r
    MTP = size(rwfn,1);
    P = [rwfn.P; zeros(NNNP-MTP,1)];
    Q = [rwfn.Q; zeros(NNNP-MTP,1)];
    rho = P.^2 + Q.^2;
    
    ND = rho;
    %     ND=0;
    ND_core = ND_core + ND*occupn;
    ND_core_subshell{shi} = ND*occupn;
end
fprintf(' Done\n')

%% Valence -- diagonal, population, single-electron density terms
fprintf('   - Valence calculation started.\n')
conf_valence = regexp(confstr_valence,expr_conf,'tokens');
conf_valence = vertcat(conf_valence{:});
conf_valence(:,2) = cellfun(@(c) str2double(c),conf_valence(:,2),'UniformOutput',false);
N = sum([conf_valence{:,2}]); % number of valence electrons
fprintf('\t- N_valence=%u\n',N)

fprintf('\t- Population terms\n')
calList_pop = cell(0,3); % list of subshells, its occupation and mixc to avoid redundant integrations
cli = 0;
for csi = 1:n_csf
    n_subshell = size(csf{csi},2);
    NDtemp = 0;
    for shi=1:n_subshell
        subshell = csf{csi}{1,shi};
        occupn = csf{csi}{2,shi};
        rwfn = T_orbital.rwfn{subshell};
        % add subshell to Ylist
        cli = cli + 1;
        calList_pop(cli,:) = {subshell,occupn,mixc(csi)};
    end
end
calList_pop_subshells = unique(calList_pop(:,1),'stable');
n_valence_subshells = length(calList_pop_subshells);
calList_pop_unique = cell(n_valence_subshells,2); % {subshell,weight}
for cli = 1:n_valence_subshells
    subshell = calList_pop_subshells{cli};
    ylInx = find(strcmp(subshell,calList_pop(:,1)));
    occupn_temp = vertcat(calList_pop{ylInx,2});
    mixc_temp = vertcat(calList_pop{ylInx,3});
    weight = sum(occupn_temp.*mixc_temp.^2);
    calList_pop_unique(cli,:) = {subshell,weight};
end
fprintf('\t\t-- List of integration generated.\n')
fprintf('\t\t-- Performing integration...')
ND_pop = zeros(NNNP,1);
ND_pop_subshell = cell(n_valence_subshells,1);
% Y_pop_csf = nan(n_csf,1);
for cli = 1:n_valence_subshells
    subshell = calList_pop_unique{cli,1};
    weight = calList_pop_unique{cli,2};
    rwfn = T_orbital.rwfn{subshell};
    
    % get number density as a function of r
    MTP = size(rwfn,1);
    P = [rwfn.P; zeros(NNNP-MTP,1)];
    Q = [rwfn.Q; zeros(NNNP-MTP,1)];
    rho = P.^2 + Q.^2;
    
    ND = rho;
    %     ND = 0;
    ND_pop = ND_pop + ND*weight;
    ND_pop_subshell{cli} = ND*weight;
end
fprintf(' Done.\n')

%% Valence --  correlation terms
fprintf('\t- Correlation terms\n')
% list CSFs in nljm form; fill m from lowest value = -j
csf_nljm = cell(n_csf,1);
for csi = 1:n_csf
    csf_nljm{csi} = cell(3,N); % {[n,l,j,m];subshell;cmpstr}
    n_subshell = size(csf{csi},2);
    ei = 1;
    for shi = 1:n_subshell
        subshell = csf{csi}{1,shi};
        nlj = str2nlj(subshell);
        occupn = csf{csi}{2,shi};
        j = nlj(3);
        mj = (1:occupn) - j -1;
        temp = mat2cell([repmat(nlj,occupn,1),mj.'],ones(1,occupn)).';
        csf_nljm{csi}(1,ei:(ei+occupn-1)) = temp;
        csf_nljm{csi}(2,ei:(ei+occupn-1)) = repmat({subshell},1,occupn);
        csf_nljm{csi}(3,ei:(ei+occupn-1)) = ...
            compose([subshell,'%u'],(1:occupn));
        
        ei = ei + occupn;
    end
end

% calculate integration over all pair csi > csj
fprintf('\t\t-- List of integration being generated...\n')
cli = 0;
calList_corr = cell(0,4);
for csi = 1:n_csf
    for csj = 1:csi - 1
        % check if two CSF has minimum difference
        cmpstri = csf_nljm{csi}(3,:);
        cmpstrj = csf_nljm{csj}(3,:);
        ei_diffi = find(~ismember(cmpstri,cmpstrj));
        n_diff = length(ei_diffi);
        if n_diff == 0
            error('the same states: \n%u: %s\n%u: %s',...
                csi,strjoin(cmpstri,' '),csj,strjoin(cmpstrj,' '))
        elseif n_diff == 1
            ei_diffj = find(~ismember(cmpstrj,cmpstri));
            %             fprintf('%u in csfi and %u in csfj\n',shell_diff1,shell_diff2)
        else
            %             fprintf('more difference.\n')
            continue;
        end
        %         fprintf('\t\t\t -- test1 passed.\n')
        dispcsfflag = false;
        % check if (l,j,m) == (l',j',m')
        lji = csf_nljm{csi}{1,ei_diffi}(2:3);
        ljj = csf_nljm{csj}{1,ei_diffj}(2:3);
        if ~isequal(lji,ljj)
            continue;
        end
        %         fprintf('\t\t\t -- test2 passed.\n')
        dispcsfflag = false;
        if dispcsfflag
            fprintf('\t\t\t -- %u.\n',cli+1)
            csfstri = strjoin(cellfun(@(c) sprintf('%s(%u)',c{1},c{2}),mat2cell(csf{csi},2,ones(1,size(csf{csi},2))),'UniformOutput',false),' ');
            csfstrj = strjoin(cellfun(@(c) sprintf('%s(%u)',c{1},c{2}),mat2cell(csf{csj},2,ones(1,size(csf{csj},2))),'UniformOutput',false),' ');
            fprintf('\t\t\tcsfi: %s\n',csfstri)
            fprintf('\t\t\tcsfj: %s\n',csfstrj)
            fprintf('\t\t\tmi*mj= %g\n',mixc(csi)*mixc(csj))
        end
        
        % run over mixed states of all possible combinations of {mj} with
        % equal probabilities
        subshelli = csf{csi}(1,:);
        subshellj = csf{csj}(1,:);
        subshell_diff1 = csf_nljm{csi}{2,ei_diffi};
        subshell_diff2 = csf_nljm{csj}{2,ei_diffj};
        occupn1i = csf{csi}{2};
        temp = find(strcmp(subshell_diff1,subshellj));
        if isempty(temp), occupn1j = 0;
        else, occupn1j = csf{csj}{2,temp}; end
        temp = find(strcmp(subshell_diff2,subshelli));
        if isempty(temp), occupn2i = 0;
        else, occupn2i = csf{csi}{2,temp}; end
        occupn2j = csf{csj}{2};
        nlj1 = csf_nljm{csi}{1,ei_diffi}(1:3);
        nlj2 = csf_nljm{csj}{1,ei_diffj}(1:3);
        j1 = nlj1(3);
        j2 = nlj2(3);
        
        %         % for test
        %         j1 = 1.5; j2 = 2.5;
        %         occupn1i = 3; occupn2i = 1; occupn1j = 2; occupn2j = 2;
        %         j1 = 0.5; j2 = 0.5;
        %         occupn1i = 2; occupn2i = 0; occupn1j = 1; occupn2j = 1;
        
        N_m1i = nchoosek(2*j1+1,occupn1i);
        N_m1j = nchoosek(2*j1+1,occupn1j);
        N_m2i = nchoosek(2*j2+1,occupn2i);
        N_m2j = nchoosek(2*j2+1,occupn2j);
        comb_2m1i = -2*j1-2 + 2*nchoosek(1:2*j1+1,occupn1i); % all comb array in order
        comb_2m1j = -2*j1-2 + 2*nchoosek(1:2*j1+1,occupn1j);
        comb_2m2i = -2*j2-2 + 2*nchoosek(1:2*j2+1,occupn2i);
        comb_2m2j = -2*j2-2 + 2*nchoosek(1:2*j2+1,occupn2j);
        
        N_comb_m = N_m1i*N_m1j*N_m2i*N_m2j;
        % phase = ...
        % generate csf_mj for two subhsells of difference with all possible combinations of mj's
        N_comb_m_match = 0;
        for cm1i = 1:N_m1i
            for cm1j = 1:N_m1j
                for cm2i = 1:N_m2i
                    for cm2j = 1:N_m2j
                        csf_2mj1i = comb_2m1i(cm1i,:);
                        csf_2mj2i = comb_2m2i(cm2i,:);
                        csf_2mj1j = comb_2m1j(cm1j,:);
                        csf_2mj2j = comb_2m2j(cm2j,:);
                        %
                        %                         csf_diffi = [compose('%s|%d',subshell_diff1,csf_2mj1i),compose('%s|%d',subshell_diff2,csf_2mj2i)];
                        %                         csf_diffj = [compose('%s|%d',subshell_diff1,csf_2mj1j),compose('%s|%d',subshell_diff2,csf_2mj2j)];
                        %                         csf_diffi, csf_diffj
                        
                        ek_diff1i = find(~ismember(csf_2mj1i,csf_2mj1j));
                        ek_diff1j = find(~ismember(csf_2mj1j,csf_2mj1i));
                        if ~isempty(ek_diff1i) && ~isempty(ek_diff1j)
                            %                 unique(            fprintf('More than one difference in subshell 1.\n')
                            continue;
                        end
                        twom1 = [csf_2mj1i(ek_diff1i),csf_2mj1j(ek_diff1j)];
                        ek_diff2i = find(~ismember(csf_2mj2i,csf_2mj2j));
                        ek_diff2j = find(~ismember(csf_2mj2j,csf_2mj2i));
                        if ~isempty(ek_diff2i) && ~isempty(ek_diff2j)
                            %                             fprintf('More than one difference in subshell 2.\n')
                            continue;
                        end
                        twom2 = [csf_2mj2i(ek_diff2i),csf_2mj2j(ek_diff2j)];
                        %                         if isempty(ek_diffi) % csf_j has one more electron in subhsell1
                        %                             ek_diff1j = find(~ismember(csf_mj1j,csf_mj1i));
                        %                         else
                        %                         [twom1,twom2]
                        if twom1 ~= twom2
                            %                             fprintf('Different mj''s.\n')
                            continue;
                        end
                        %                         fprintf('All passed.\n')
                        N_comb_m_match = N_comb_m_match + 1;
                    end
                end
            end
        end
        if dispcsfflag
            fprintf('\t\t\tP = %u/%u\n',N_comb_m_match,N_comb_m)
        end
        P = N_comb_m_match/N_comb_m; % stat. prob of (ljm)=(l'j'm')
        
        % add subshell to Ylist
        cli = cli + 1;
        phase = (-1)^(ei_diffi-ei_diffj);
        weight = 2*phase*mixc(csi)*mixc(csj)*P;
        %         fprintf('\t\t\tphase= %+d\n',phase)
        subshellpair = {csf_nljm{csi}{2,ei_diffi},csf_nljm{csj}{2,ei_diffi}};
        subshellpair = sort(subshellpair);
        cmpstr = strjoin(subshellpair,' ');
        calList_corr(cli,:) = ...
            {subshellpair{:},cmpstr,weight}; % {subshell1,subshell2,cmpstr,weight}
    end
end
calList_corr_cmpstrs = unique(calList_corr(:,3),'stable');
n_corr = length(calList_corr_cmpstrs);
calList_corr_unique = cell(n_corr,3);
for cli=1:n_corr
    cmpstr = calList_corr_cmpstrs{cli};
    ylInx = find(strcmp(cmpstr,calList_corr(:,3)));
    subshellpair = calList_corr(ylInx(1),1:2);
    weight = sum([calList_corr{ylInx,4}]);
    calList_corr_unique(cli,:) = {subshellpair{:},weight}; % {subshell1,subshell2,weight}
end

[temp,I] = sort(abs([calList_corr_unique{:,3}]),'descend');
calList_corr_unique = calList_corr_unique(I,:);

fprintf('\t\t-- Performing integration...')
ND_corr = zeros(NNNP,1);
ND_corr_subshellpair = cell(n_corr,1);
% Y_corr_csfpair = nan(n_csf,n_csf);
for cli = 1:n_corr
    subshell1 = calList_corr_unique{cli,1};
    subshell2 = calList_corr_unique{cli,2};
    weight = calList_corr_unique{cli,3};
    rwfn1 = T_orbital.rwfn{subshell1};
    rwfn2 = T_orbital.rwfn{subshell2};
    
    % get number correlation density as a function of r
    MTP1 = size(rwfn1,1);
    P1 = [rwfn1.P; zeros(NNNP-MTP1,1)];
    Q1 = [rwfn1.Q; zeros(NNNP-MTP1,1)];
    MTP2 = size(rwfn2,1);
    P2 = [rwfn2.P; zeros(NNNP-MTP2,1)];
    Q2 = [rwfn2.Q; zeros(NNNP-MTP2,1)];
    rho = P1.*P2 + Q1.*Q2;
    
    ND = rho;
%     ND = 0;
    ND_corr = ND_corr + ND*weight;
    ND_corr_subshellpair{cli} = ND*weight;
end
fprintf(' Done\n')

%% Normalizations
fprintf('\t- Normalizations\n')
Int_core = trapz(r,ND_core);
norm_err_core = Int_core - N_core;
Int_valence = trapz(r,ND_pop + ND_corr);
norm_err_valence = Int_valence - N;
fprintf('\t\t- norm_err_core = %g (%.3g%%)\n',norm_err_core,norm_err_core/N_core*100)
fprintf('\t\t- norm_err_valence = %g (%.3g%%)\n',norm_err_valence,norm_err_valence/N*100)
ND_core = ND_core*N_core/Int_core;
ND_pop = ND_pop*N/Int_valence;
ND_corr = ND_corr*N/Int_valence;

%% Total Yukawa integral and output

ND = ND_core + ND_pop + ND_corr;

varargout{1} = ND_core;
varargout{2} = ND_pop;
varargout{3} = ND_corr;
varargout{4} = ND_pop_subshell;
varargout{5} = ND_corr_subshellpair;
end
