function [Iidf Mt Iexp l T t Nsc Nxc N_tsv_gates N_tsvs_1d Nc Nb N_ns_t N_ns_b t2 t3] = compare_calculations(use_corrected, Ng, N_strata,r,p,fo,alpha,k,max_area_ratio,w)

%% Debug stuff
DEBUG = 1;
fname = 'compare_calculations';

%%
Nt = Ng; % total number of gates
S = N_strata; % number of strata
Ns = Nt/S; % Avg number of gates per stratum
Nx = sqrt(Ns);

% get corrected mt_intra
N_tsv_gates = Ns* max_area_ratio/(1-max_area_ratio);

N_tsv_gates_1d = floor(sqrt(N_tsv_gates));
N_tsvs_1d = floor(N_tsv_gates_1d/w);

% corrected gate counts (length of chip in gates, as opposed to actual
% number of gates)
Nsc = ceil(Ns + N_tsv_gates);
Nxc = floor(sqrt(Nsc));
Ntc = S*Nsc;

T = floor(Nxc/N_tsvs_1d);
t = floor(1/2*(T-w));

lmax_2d = 2*sqrt(Ns); % max 2D wirelength (assuming optimal routing)
lmax_3d = lmax_2d + (N_strata-1)*r; % max 3D wirelength (assuming optimal routing)
lmax = lmax_3d;

if(use_corrected == 1)
    lmax_2d = 2*ceil(sqrt(Nsc));   % max 2D wirelength (assuming optimal routing)
	lmax_3d = lmax_2d + (N_strata-1)*r;      % max 3D wirelength (assuming optimal routing)
	lmax = floor(lmax_3d);
    
    if (exist('DEBUG','var') == 1)
        if(DEBUG == 1) % debug is a global variable
            str = sprintf('%s::\tlmax_2d %d\tlmax_3d %d\tlmax %d',fname,lmax_2d,lmax_3d,lmax);
            disp(str);
        end
    end
    
    [Mt_intra t2 t3] = calc_Mtcorr_newcalc(Nsc,Nxc,t,w,T,N_tsvs_1d);
    
    if (exist('DEBUG','var') == 1)
        if(DEBUG == 1) % debug is a global variable
            str = sprintf('%s::\tNs: %d\tNsc: %d\tNx: %d\tNxc: %d',fname,Ns,Nsc,Nx,Nxc);
            disp(str);
        end
    end
else
    Nxc = Nx;
    Nsc = Ns;
    Ntc = Nt;
    
    lmax_2d = floor(2*sqrt(Nsc));     % max 2D wirelength (assuming optimal routing)
	lmax_3d = lmax_2d + (N_strata-1)*r;         % max 3D wirelength (assuming optimal routing)
	lmax = floor(lmax_3d);

	Mt_intra = getMt_intra_vec(Ns,lmax_2d);
end

lengths = 0:lmax;
Mt = getMt_vec(Nsc,lmax,r,S,Mt_intra);
[Nc, N_ns_t, N_ns_b] = getNc_vec(Ntc, Nsc, lmax, Mt, S, r);
Nb = getNb_vec(Nc);

Iexp = getIexp_vec(alpha,k,p,Nb,Nc);
Iidf = Iexp.*Mt;

% 
% 
% %% calculate everything!
% for i=1:length(l)
%     Mt(i) = getMt(Nsc,l(i),r,S,Mt_intra_corr,use_corrected);
%     
%     Nb(i) = sum(Nc) - Nc(1); % want the sum of every Nc from l=1 to l=i-1
%     
%     N_nonstart_top(i) = getN_nonstart_top(Nsc, l(i));
%     N_nonstart_bot(i) = getN_nonstart_bot(Nsc, l(i), S, r);
%     Nstart(i) = Nt - N_nonstart_top(i) - N_nonstart_bot(i);
%     
%     Nc(i) = getNc(Ntc, Nsc, l(i), Mt(i), S, r);
% 
%     Iexp(i) = getIexp(alpha,k,p, Nb(i), Nc(i) );
%     Iidf(i) = Iexp(i)*Mt(i);
%     
% end

%% Cleanup for returning stuff
l = 0:lmax_3d;

