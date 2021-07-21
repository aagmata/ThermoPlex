function [F,G,H,I,J,K,L] = ThermoDhyb(Sequence1,Sequence2,spectemp,Mgconc)

%   Function for approximating the ensemble gibbs free energy for two 
%   hybridizing nucleic acid based on nearest-neighbor thermodynamic 
%   parameters (SantaLucia et al., 1998) with salt correction from Owczarzy 
%   et al. (2008). 
%
%   Input arguments:
%      Sequence1 - Sequence 1 in 5' -> 3' 
%      Sequence2 - Sequence 2 in 5' -> 3'
%      spectemp - specified temperature of hybridization
%      Mgconc - concentration of divalent Mg2+ ions in the solution.
%      
%   Output arguments 
%      F - ensemble Gibbs free energy of interaction
%      G - minimum Gibbs free energy
%      H - enables plotting of minimum free energy structure
%      I - enables user control to view all suboptimal structures
%      J - P structure containing all the info about each doublet pairing
%      K - P' structure containing all the info about each doublet pairing
%          propagated opposite the direction of P


% =========================================================================
%% Main Program
% =========================================================================

% Table anatomy:

%   3'- S e q u e n c e 2 -5'
% 5'   ___________________
% |   |                   |
% S   |                   |
% e   |                   |
% q   |                   |
% u   |                   |
% e   |                   |
% n   |                   |
% c   |                   |
% e   |                   |
% 1   |                   |
% |   |___________________|
% 3'
% 
% index list:
% i = row of table; doublet of Sequence 1
% j = column of table; doublet of Sequence 2
% -------------------------------------------------------------------------
%% Bidirectional filling of energy table
% -------------------------------------------------------------------------
if nargout>3
    xbar=waitbar(0,'Initiating parameters');
    titleHandle = get(findobj(xbar,'Type','axes'),'Title');
    set(titleHandle,'FontSize',9)
end
% retrieving the energy table for the given temperature and salt
% concentration conditions

fgc = interactgc(Sequence1,seqrcomplement(Sequence2));
np = (sum([length(Sequence1) length(Sequence2)])/2)-1;

[aaa,bbb,ccc,ddd,eee] = dnaenergies(spectemp,Mgconc,fgc,np);
[tabint,tabterm] = dnatablearrange(aaa,bbb,ccc,ddd,eee);

% data for bulge and hairpin loop penalty

spectempK = spectemp+273.15;
                    %1   2   3   4   5   6   7   8   9   10
bulgeparaminitdG = [4.0 2.9 3.1 3.2 3.3 3.5 3.7 3.9 4.1 4.3 0.2*(1:1000)+4.3];
bulgeparamdS = (bulgeparaminitdG.*1000)./310.15;
                   %2   4   6   8  10  12  14
iloopparaminitdG = [0 3.6 0.8 0.4 0.1 0.3 0.2 zeros(1,50)];
iloopparamdS = (iloopparaminitdG.*1000)./310.15;
bulgeparam = bulgeparamdS.*(spectempK)./1000;
iloopparam = iloopparamdS.*(spectempK)./1000;

% declaring global variables
global P Pp prime_state dyn_coords
if nargout>3
    waitbar(0.30,xbar,'Calculating P Matrix (30%)');
end
P = Hybridize(Sequence1,Sequence2);
if nargout>3
    waitbar(0.60,xbar,'Calculating P'' Matrix (60%)');
end
Pp = Hybridize(Sequence2,Sequence1);

% Rearrange P' tables
Pp.NNdG = rot90(Pp.NNdG',2);
Pp.dGmn = rot90(Pp.dGmn',2);
Pp.budG = rot90(Pp.budG',2);
Pp.bldG = rot90(Pp.bldG',2);
Pp.casescore = rot90(Pp.casescore',2);
Pp.terminal = rot90(Pp.terminal',2);
Pp.pathcase = rot90(Pp.pathcase',2);
Pp.bulgecumU = rot90(Pp.bulgecumU',2);
Pp.bulgecumL = rot90(Pp.bulgecumL',2);
Pp.loopcum = rot90(Pp.loopcum',2);
Pp.loopstate = rot90(Pp.loopstate',2);
Pp.doublet1 = rot90(Pp.doublet1',2);
Pp.doublet2 = rot90(Pp.doublet2',2);
Pp.doublet3 = rot90(Pp.doublet3',2);

[A,initiation,ATterm] = NNdG('GA','TC',tabint,tabterm,1);

% -------------------------------------------------------------------------
%% Consolidating and finalizing P & P' energies
% -------------------------------------------------------------------------
if nargout>3
    waitbar(0.80,xbar,'Consolidating P and P'' energies(60%)')
end
[ox,oy] = find(Pp.pathcase == 4);
logic_indx = logical(zeros(1,length(ox)));

mincoords_path = {};

% Tracing
for i = 1:length(ox)
    prime_state = 0;
    dyn_coords = [];
    pathtrace(ox(i),oy(i),[]);
    loop_path = dyn_coords;
    if isempty(loop_path) == 0
        if sum(loop_path(1,:) == loop_path(end,:)) == 2
            loop_path = loop_path(1:ceil(end/2),:);
            mincoords_path{i} = loop_path;
            logic_indx(i) = 1;
        end
    end
end

if isempty(mincoords_path)
    F = 0;
    G = 0;
    J = 0;
    return
else
    mincoords_path = mincoords_path(logic_indx);
end

% Resolving terminal energies

cscr_mat = zeros(length(mincoords_path));
Pp_NNdG_term = zeros(length(mincoords_path),1);
P_NNdG_term = zeros(length(mincoords_path),1);

for i = 1:length(mincoords_path)
    cscr_mat(i) = P.casescore{mincoords_path{i}(1,1)...
        ,mincoords_path{i}(1,2)}(1);
    Pp_NNdG_term(i) = ...
        Pp.NNdG(mincoords_path{i}(1,1),mincoords_path{i}(1,2));
    P_NNdG_term(i) = ...
        P.NNdG(mincoords_path{i}(1,1),mincoords_path{i}(1,2));
end

Energy_optima0 = cscr_mat(:,1) + ...
    Pp_NNdG_term - P_NNdG_term + initiation;

Energy_optima = Energy_optima0(Energy_optima0 < 0);
mincoords_path = mincoords_path(Energy_optima0 < 0);

% Resolving AT terminal penalties
if isempty(Energy_optima)
    Energy_optima = 0;
else
    for i = 1:length(Energy_optima)
        if sum(P.doublet1{mincoords_path{i}(1,1),mincoords_path{i}(1,2)}...
            (:,1) == ((['A';'T']))) == 2 || ...
            sum(P.doublet1{mincoords_path{i}(1,1),mincoords_path{i}(1,2)}...
            (:,1) == ((['T';'A']))) == 2 || ...
            sum(P.doublet1{mincoords_path{i}(1,1),mincoords_path{i}(1,2)}...
            (:,2) == ((['A';'T']))) == 2 || ...
            sum(P.doublet1{mincoords_path{i}(1,1),mincoords_path{i}(1,2)}...
            (:,2) == ((['T';'A']))) == 2

            Energy_optima(i) = Energy_optima(i) + ATterm;
        end
        
        if sum(P.doublet1{mincoords_path{i}(end,1),mincoords_path{i}(end,2)}...
            (:,1) == ((['A';'T']))) == 2 || ...
            sum(P.doublet1{mincoords_path{i}(end,1),mincoords_path{i}(end,2)}...
            (:,1) == ((['T';'A']))) == 2 || ...
            sum(P.doublet1{mincoords_path{i}(end,1),mincoords_path{i}(end,2)}...
            (:,2) == ((['A';'T']))) == 2 || ...
            sum(P.doublet1{mincoords_path{i}(end,1),mincoords_path{i}(end,2)}...
            (:,2) == ((['T';'A']))) == 2
        
            Energy_optima(i) = Energy_optima(i) + ATterm;
        end
    end
end

% -------------------------------------------------------------------------
%% Partition function and ensemble free energy
% -------------------------------------------------------------------------

Z = sum(exp(-(Energy_optima.*1000)./(1.986*(spectemp+273.15))));           % partition function of interaction
dGensemble = -(1.986*(spectemp+273.15))*log(Z)./1000;


[Enopt_sort,sort_add] = sort(Energy_optima);

F = dGensemble;
G = min(Energy_optima);
J = Enopt_sort;
L = Pp;

% -------------------------------------------------------------------------
%% Plotting the energy optima structures
% -------------------------------------------------------------------------

if nargout == 3
    if Energy_optima == 0
        fprintf('\n')
        fprintf('**** NO STABLE STRUCTURES PREDICTED ****')
        fprintf('\n')
        fprintf('\n')
        H = [];
        return
    end
    mfe_path = mincoords_path{Energy_optima == G};
    H = plot_strand(mfe_path,Sequence1,Sequence2,P,G);
elseif nargout > 3
    waitbar(0.90,xbar,'Plotting energy optima structures (90%)')
    if Energy_optima == 0
        fprintf('\n')
        fprintf('**** NO STABLE STRUCTURES PREDICTED ****')
        fprintf('\n')
        fprintf('\n')
        H = [];
        return
    end
    mfe_path = mincoords_path{Energy_optima == G};
    parfor k = 1:size(mincoords_path,2);
        l = sort_add(k);
        I(k) = plot_strand(mincoords_path{l},...
            Sequence1,Sequence2,P,Energy_optima(l));
    end
    H = plot_strand(mfe_path,Sequence1,Sequence2,P,G);
    close(xbar)
end
% =========================================================================
%% Hybridize Subfunction
% =========================================================================

function X = Hybridize(Sequence1,Sequence2)

Sequence1 = ['-' Sequence1 '-'];
Sequence2 = ['-' Sequence2 '-'];
Sequence2 = seqreverse(Sequence2);
    
% Establishing sizes
M = length(Sequence1);
N = length(Sequence2);

% Preallocation

NNdG_tab = zeros(M-1,N-1);                                                 % Nearest neighbor doublet values
dGmn_tab = NNdG_tab;                                                       % Minimum case cummulative values
budG_tab = NNdG_tab;                                                       % Bulge dG up values 
bldG_tab = NNdG_tab;                                                       % Bulge dG left values
cscr_tab = cell(M-1,N-1);                                                  % Casescore array
term_tab = NNdG_tab;                                                       % Terminal for computation
path_tab = NNdG_tab;                                                       % path case
bcmU_tab = NNdG_tab;                                                       % bulge penalty accumulation (up)
bcmL_tab = NNdG_tab;                                                       % bulge penalty accumulation (left)
lpcm_tab = NNdG_tab;                                                       % loop penalty accumulation
lpst_tab = NNdG_tab;                                                       % logical loop state
doublet1 = cell(M-1,N-1);                                                  % doublet comparison diagonal case
doublet2 = cell(M-1,N-1);                                                  % doublet comparison up case
doublet3 = cell(M-1,N-1);                                                  % doublet comparison left case


for m = 1:M-1
    for n = 1:N-1
        
        % Terminal case for table boundaries
        if m == 1 || n == 1
            term_tab(m,n) = 1;
        elseif dGmn_tab(m-1,n-1) == 0
            term_tab(m,n)  = 1;
        else
            term_tab(m,n)  = 0;
        end
        
        
        % Doublet pairs for all case
        doublet1{m,n} = [Sequence1(m:m+1);Sequence2(n:n+1)];
        if m ~= 1   % upcase
            gapstring2 = num2str([]);
            gapstring2(1:bcmU_tab(m-1,n)+1) = '-';
            doublet2{m,n} = [Sequence1(m-bcmU_tab(m-1,n)-1:m+1);...
                Sequence2(n),gapstring2,Sequence2(n+1)];
        else
            doublet2{m,n} = ['---';'---'];
        end
        if n ~= 1   % leftcase
            gapstring3 = num2str([]);
            gapstring3(1:bcmL_tab(m,n-1)+1) = '-';
            doublet3{m,n} = [Sequence1(m),gapstring3,Sequence1(m+1);...
                Sequence2(n-bcmL_tab(m,n-1)-1:n+1)];
            
        else
            doublet3{m,n} = ['---';'---'];
        end
        
        
        % Calculation of NNdG, budG and bldG
        NNdG_tab(m,n) = NNdG(doublet1{m,n}(1,:),doublet1{m,n}(2,:),...
             tabint,tabterm,term_tab(m,n));
        
        auxdoub2 = [doublet2{m,n}(:,1),doublet2{m,n}(:,end)];
        auxdoub3 = [doublet3{m,n}(:,1),doublet3{m,n}(:,end)];
        
        if  m ~= 1 && bcmU_tab(m-1,n) == 0
            budG_tab(m,n) = NNdG(auxdoub2(1,:),auxdoub2(2,:),tabint,tabterm,...
                0);
        else
            budG_tab(m,n) = 0;
        end
        
        if  n ~= 1 && bcmL_tab(m,n-1) == 0
            bldG_tab(m,n) = NNdG(auxdoub3(1,:),auxdoub3(2,:),tabint,tabterm,...
                0);
        else
            bldG_tab(m,n) = 0;
        end
        
        % internal loop initiation and accumulation
        
        ddGL = 0;
        if NNdG_tab(m,n) == 0.0001
            lpcm_tab(m,n) = lpcm_tab(m-1,n-1) + 1;
            if lpst_tab(m-1,n-1) ~= 1
                ddGL = NNdG(doublet1{m-1,n-1}(1,:),...
                    doublet1{m-1,n-1}(2,:),tabint,tabterm,1) - ...
                    NNdG_tab(m-1,n-1);
            end
        end
        
        if m ~= 1 && n ~= 1
            if lpcm_tab(m,n) == 1 
                lpst_tab(m,n) = 1;
            else lpst_tab(m,n) = lpst_tab(m-1,n-1);
            end
        end
        
       
        % Computing case energy
        if m == 1 || n == 1
            cscr_tab{m,n} = [NNdG_tab(m,n) 0 0 0];
            dGmn_tab(m,n) = min(cscr_tab{m,n});
            if dGmn_tab(m,n) < 0
                path_tab(m,n) = 4;
            else
                path_tab(m,n) = 0;
            end
            
        else
            if m == 2 && n == 2
                cscr_tab{m,n} = [NNdG_tab(m,n)+dGmn_tab(m-1,n-1)+ddGL+...
                         iloopparam(lpcm_tab(m,n)+1),...
                         0,...
                         0,...
                         0];
                     
             elseif m ~= 2 && n == 2
                
                cscr_tab{m,n} = [NNdG_tab(m,n)+dGmn_tab(m-1,n-1)+ddGL+...
                         iloopparam(lpcm_tab(m,n)+1),...
                         min([budG_tab(m,n),0])+...
                         bulgeparam(bcmU_tab(m-1,n)+1)+...
                         dGmn_tab(m-1-bcmU_tab(m-1,n)-1,n-1),...
                         0,...
                         0];
                     
            elseif m == 2 && n ~= 2
                
                cscr_tab{m,n} = [NNdG_tab(m,n)+dGmn_tab(m-1,n-1)+ddGL+...
                         iloopparam(lpcm_tab(m,n)+1),...
                         0,...
                         min([bldG_tab(m,n),0])+...
                         bulgeparam(bcmL_tab(m,n-1)+1)+...
                         dGmn_tab(m-1,n-1-bcmL_tab(m,n-1)-1),...
                         0];
            
            else
                cscr_tab{m,n} = [NNdG_tab(m,n)+dGmn_tab(m-1,n-1)+ddGL+...
                         iloopparam(lpcm_tab(m,n)+1),...
                         min([budG_tab(m,n),0])+...
                         bulgeparam(bcmU_tab(m-1,n)+1)+...
                         dGmn_tab(m-1-bcmU_tab(m-1,n)-1,n-1),...
                         min([bldG_tab(m,n),0])+...
                         bulgeparam(bcmL_tab(m,n-1)+1)+...
                         dGmn_tab(m-1,n-1-bcmL_tab(m,n-1)-1),...
                         0];
            
            end
            
            
            
            % internal loop termination
            if lpcm_tab(m-1,n-1) >= 1
                ddGR = NNdG(doublet1{m,n}(1,:),doublet1{m,n}(2,:),...
                    tabint,tabterm,1)-NNdG_tab(m,n);
                cscr_tab{m,n}(1) = cscr_tab{m,n}(1)+ddGR;
            end
            

            dGmn_tab(m,n) = min(cscr_tab{m,n});
            
            
            switch dGmn_tab(m,n)
                case cscr_tab{m,n}(1)
                    if cscr_tab{m,n}(1) ~= 0
                        path_tab(m,n) = 1;
                        if dGmn_tab(m-1,n-1) == 0 && dGmn_tab(m,n) < 0
                            path_tab(m,n) = 4;
                        end
                        bcmU_tab(m,n) = 0;
                        bcmL_tab(m,n) = 0;
                    else
                        path_tab(m,n) = 0;
                        bcmU_tab(m,n) = 0;
                        bcmL_tab(m,n) = 0;    
                    
                    end
                    
                case cscr_tab{m,n}(2)
                    if cscr_tab{m,n}(2) ~= 0
                        path_tab(m,n) = 2;
                        bcmU_tab(m,n) = bcmU_tab(m-1,n)+1;
                    else
                        path_tab(m,n) = 0;
                        bcmU_tab(m,n) = 0;
                        bcmL_tab(m,n) = 0;
                    end
                    lpcm_tab(m,n) = 0;
                    
                case cscr_tab{m,n}(3)
                    lpcm_tab(m,n) = 0;
                    if cscr_tab{m,n}(3) ~= 0
                        path_tab(m,n) = 3;
                        bcmL_tab(m,n) = bcmL_tab(m,n-1)+1;
                    else
                        path_tab(m,n) = 0;
                        bcmU_tab(m,n) = 0;
                        bcmL_tab(m,n) = 0;
                    end
                    
                case 0
                    path_tab(m,n) = 0;
                    bcmU_tab(m,n) = 0;
                    bcmL_tab(m,n) = 0;
                    lpcm_tab(m,n) = 0;
            end
            
               
        end
         
        
    end
end


X.NNdG = NNdG_tab;
X.dGmn = dGmn_tab;
X.budG = budG_tab;
X.bldG = bldG_tab;
X.casescore = cscr_tab;
X.terminal = term_tab;
X.pathcase = path_tab;
X.bulgecumU = bcmU_tab;
X.bulgecumL = bcmL_tab;
X.loopcum = lpcm_tab;
X.loopstate = lpst_tab;
X.doublet1 = doublet1;
X.doublet2 = doublet2;
X.doublet3 = doublet3;

end

% =========================================================================
%% Path-Tracing Subfunction
% =========================================================================

function pathtrace(seed_coordx,seed_coordy,coord_mat)
    
if prime_state == 0
    Path_mat = P.pathcase;
elseif prime_state == 1
    Path_mat = Pp.pathcase;
end

switch Path_mat(seed_coordx,seed_coordy)
    case 0
        coord_mat = [coord_mat;seed_coordx,seed_coordy];
        dyn_coords = coord_mat;
        return
    case 1
        move_vec = [-1,-1];
    case 2
        if prime_state == 0
            move_vec = [-1-P.bulgecumU(seed_coordx,seed_coordy),-1];
        else
            move_vec = [-1,-1-Pp.bulgecumU(seed_coordx,seed_coordy)];
        end
    case 3
        if prime_state == 0
            move_vec = [-1,-1-P.bulgecumL(seed_coordx,seed_coordy)];
        else
            move_vec = [-1-Pp.bulgecumL(seed_coordx,seed_coordy),-1];
        end
    case 4
        if prime_state == 0
            move_vec = [0,0];
            prime_state = 1;
        else
            coord_mat = [coord_mat;seed_coordx,seed_coordy];
            dyn_coords = coord_mat;
            return
        end
    case 5
        return
end

if prime_state == 1
    move_vec = -move_vec;
end

coord_mat = [coord_mat;seed_coordx,seed_coordy];

pathtrace(seed_coordx+move_vec(1),seed_coordy+move_vec(2),coord_mat);


end

end

% =========================================================================
%% Strand-Plotting Subfunction
% =========================================================================

function f = plot_strand(coord_path,Sequence1,Sequence2,P,Energy)

Csize = size(coord_path,1);
cp = flipud(coord_path);
clear plot_structL plot_structR
clf
h1 = figure('doublebuffer','off','Visible','Off');
set(0,'defaultfigurecolor',[1 1 1]);
set(gcf,'units','normalized','Position',[0.25,0,0.5,1]);
% i-1 = coord_path index

loopc = 0;
termstrand_0 = 0;
termstrand_end = 0;
auxdoubloop = [];
termstrand_0 = 0;
termstrand_end = 0;
term0countL = 0;
term0countR = 0;
termendcountL = 0;
termendcountR = 0;
plot_struct0L = struct('coord',[],'nbase',[],'bond',[],'color',[]);
plot_struct0R = struct('coord',[],'nbase',[],'bond',[],'color',[]);
plot_structendL = struct('coord',[],'nbase',[],'bond',[],'color',[]);
plot_structendR = struct('coord',[],'nbase',[],'bond',[],'color',[]);


%% Setting coordinates for internal pairing strands
for i = 2:Csize+1
    if i == 2
        if length(find(P.doublet1{cp(1,1),cp(1,2)}(1,:)...
                == seqcomplement(P.doublet1{cp(1,1),cp(1,2)}(2,:)))) == 2
            plot_structL(i-1).coord = [0 0];
            plot_structL(i-1).nbase = P.doublet1{cp(1,1),cp(1,2)}(1,1);
            plot_structL(i-1).bond = 1;
            
            plot_structL(i).coord = [0 1];
            plot_structL(i).nbase = P.doublet1{cp(1,1),cp(1,2)}(1,2);
            plot_structL(i).bond = 1;

            plot_structR(i-1).coord = [1 0];
            plot_structR(i-1).nbase = P.doublet1{cp(1,1),cp(1,2)}(2,1);
            plot_structR(i-1).bond = 1;
            
            plot_structR(i).coord = [1 1];
            plot_structR(i).nbase = P.doublet1{cp(1,1),cp(1,2)}(2,2);
            plot_structR(i).bond = 1;
            
        elseif P.doublet1{cp(1,1),cp(1,2)}(1,1)...
                == seqcomplement(P.doublet1{cp(1,1),cp(1,2)}(2,1))
            plot_structL(i-1).coord = [0 0];
            plot_structL(i-1).nbase = P.doublet1{cp(1,1),cp(1,2)}(1,1);
            plot_structL(i-1).bond = 1;
            
            plot_structR(i-1).coord = [1 0];
            plot_structR(i-1).nbase = P.doublet1{cp(1,1),cp(1,2)}(2,1);
            plot_structR(i-1).bond = 1;
            
            Coords = compu_circ(6,...
                    [plot_structR(end).coord;plot_structL(end).coord]);
                
            plot_structL(i).coord = Coords(3,:);
            plot_structL(end).nbase = P.doublet1{cp(1,1),cp(1,2)}(1,2);
            plot_structR(i).coord = Coords(end,:);
            plot_structR(end).nbase = P.doublet1{cp(1,1),cp(1,2)}(2,2);
            i = 3;
            plot_structL(i).coord = Coords(4,:);
            plot_structL(end).nbase = P.doublet1{cp(2,1),cp(2,2)}(1,2);
            plot_structR(i).coord = Coords(end-1,:);
            plot_structR(end).nbase = P.doublet1{cp(2,1),cp(2,2)}(2,2);
            
            
        elseif P.doublet1{cp(1,1),cp(1,2)}(1,2)...
                == seqcomplement(P.doublet1{cp(1,1),cp(1,2)}(2,2))
            
            plot_structL(i-1).coord = [0 0];
            plot_structL(i-1).nbase = P.doublet1{cp(1,1),cp(1,2)}(1,2);
            plot_structL(i-1).bond = 1;
            
            plot_structR(i-1).coord = [1 0];
            plot_structR(i-1).nbase = P.doublet1{cp(1,1),cp(1,2)}(2,2);
            plot_structR(i-1).bond = 1;
            
            termstrand_0 = 2;
            
            term0countL = 1;
            term0countR = 1;
            
        end
        
        
    else
        switch P.pathcase(cp(i-1,1),cp(i-1,2))
            case 1
                if length(find(P.doublet1{cp(i-1,1),cp(i-1,2)}(1,:)...
                 ==seqcomplement(P.doublet1{cp(i-1,1),cp(i-1,2)}(2,:))))==2
                    
                    Coords = compu_circ(4,...
                    [plot_structR(end).coord;plot_structL(end).coord]);
                    
                    plot_structL(end+1).coord = Coords(end-1,:);
                    plot_structL(end).nbase = P.doublet1{cp(i-1,1),cp(i-1,2)}(1,2);
                    plot_structL(end).bond = 1;
                    plot_structR(end+1).coord = Coords(end,:); 
                    plot_structR(end).nbase = P.doublet1{cp(i-1,1),cp(i-1,2)}(2,2);
                    plot_structR(end).bond = 1;
                    loopc = 0;
                else
                    loopc = 1 + loopc;
                    auxdoubloop = [auxdoubloop,...
                        P.doublet1{cp(i-1,1),cp(i-1,2)}(:,2)];
                    
                    if P.doublet1{cp(i-1,1),cp(i-1,2)}(1,2)...
                       ==seqcomplement(P.doublet1{cp(i-1,1),cp(i-1,2)}(2,2))
                       Coords = compu_circ(loopc*2+2,[plot_structR(end).coord;...
                           plot_structL(end).coord]);
                       Coords = Coords(3:end,:);
                       
                       for j = 1:loopc
                           plot_structL(end+1).coord = Coords(j,:);
                           plot_structL(end).nbase = P.doublet1...
                               {cp(i-loopc+j-1,1),cp(i-loopc+j-1,2)}(1,2);
                           plot_structR(end+1).coord = Coords(end-j+1,:);
                           plot_structR(end).nbase = P.doublet1...
                               {cp(i-loopc+j-1,1),cp(i-loopc+j-1,2)}(2,2);
                       end
                       
                       plot_structL(end).bond = 1;
                       plot_structR(end).bond = 1;
                       
                       
                    elseif i == Csize+1
                        termendcountL = 1;
                        termendcountR = 1;
                        termstrand_end = 2;
                    end
                    
                end
                
                
            case 2
                
                Coords = compu_circ(size(P.doublet2{cp(i-1,1),cp(i-1,2)},2)+2,...
                    [plot_structR(end).coord;plot_structL(end).coord]);
                
                for j = 1:size(Coords,1)-3
                    if j == size(Coords,1)-3
                        plot_structR(end+1).coord = Coords(end,:);
                        plot_structR(end).nbase = P.doublet2{cp(i-1,1),cp(i-1,2)}(2,end);
                    else
                        plot_structR(end+1).coord = Coords(1,:);
                        plot_structR(end).nbase = P.doublet2{cp(i-1,1),cp(i-1,2)}(2,1);
                    end
                    plot_structL(end+1).coord = Coords(2+j,:);
                    plot_structL(end).nbase = P.doublet2{cp(i-1,1),cp(i-1,2)}(1,j+1);
                end
                
                plot_structL(end).bond = 1;
                plot_structR(end).bond = 1;
                
            case 3
                Coords = compu_circ(size(P.doublet3{cp(i-1,1),cp(i-1,2)},2)+2,...
                    [plot_structR(end).coord;plot_structL(end).coord]);
                
                for j = 1:size(Coords,1)-3
                    if j == size(Coords,1)-3
                        plot_structL(end+1).coord = Coords(3,:);
                        plot_structL(end).nbase = P.doublet3{cp(i-1,1),cp(i-1,2)}(1,end);
                    else
                        plot_structL(end+1).coord = Coords(2,:);
                        plot_structL(end).nbase = P.doublet3{cp(i-1,1),cp(i-1,2)}(1,1);
                    end
                    plot_structR(end+1).coord = Coords(end+1-j,:);
                    plot_structR(end).nbase = P.doublet3{cp(i-1,1),cp(i-1,2)}(2,j+1);
                end
                
                plot_structL(end).bond = 1;
                plot_structR(end).bond = 1;
                
                
            case 4
                return
        end
        
        
        
    end
end

%% Setting coordinates for terminal non-pairing strands

termstrand_0 = 2+termstrand_0+cp(1,1)+cp(1,2);
termstrand_end = 6+termstrand_end+(length(Sequence1)-cp(end,1))...
    +(length(Sequence2)-cp(end,2));
coords0 = compu_circ(termstrand_0,...
    [plot_structL(1).coord;plot_structR(1).coord]);
coordsend = compu_circ(termstrand_end,...
    [plot_structR(end).coord;plot_structL(end).coord]);


term0countL = term0countL+cp(1,1)-2;
term0countR = term0countR+cp(1,2)-2;
termendcountL = termendcountL+(length(Sequence1)-cp(end,1));
termendcountR = termendcountR+(length(Sequence2)-cp(end,2));


% Setting termstrand_0

% Left
for L = 1:term0countL
    plot_struct0L(L).coord = [coords0(end-(term0countL-L),1),...
                              coords0(end-(term0countL-L),2)];
    plot_struct0L(L).nbase = P.doublet1{L,1}(1,2);
    plot_struct0L(L).bond = 0;
    
    switch plot_struct0L(L).nbase
        case 'A'
            plot_struct0L(L).color = [1 0.502 0];
        case 'T'
            plot_struct0L(L).color = [0.0549 0.3020 0.5725];
        case 'G'
            plot_struct0L(L).color = [0 0.6 0.298];
        case 'C'
            plot_struct0L(L).color = [1 1 0];
    end

end

% Right
for R = 1:term0countR
    plot_struct0R(R).coord = [coords0(term0countR+3-R,1),...
                              coords0(term0countR+3-R,2)];
    plot_struct0R(R).nbase = P.doublet1{1,R}(2,2);
    plot_struct0R(R).bond = 0;
    
    switch plot_struct0R(R).nbase
        case 'A'
            plot_struct0R(R).color = [1 0.502 0];
        case 'T'
            plot_struct0R(R).color = [0.0549 0.3020 0.5725];
        case 'G'
            plot_struct0R(R).color = [0 0.6 0.298];
        case 'C'
            plot_struct0R(R).color = [1 1 0];
    end
  
end



% Setting termstrand_end

% Left
for L = 1:termendcountL
    plot_structendL(L).coord = [coordsend(termendcountL+3-L,1),...
                              coordsend(termendcountL+3-L,2)];
    plot_structendL(L).nbase = P.doublet1{end+1-L,1}(1,1);
    plot_structendL(L).bond = 0;
    
    switch plot_structendL(L).nbase
        case 'A'
            plot_structendL(L).color = [1 0.502 0];
        case 'T'
            plot_structendL(L).color = [0.0549 0.3020 0.5725];
        case 'G'
            plot_structendL(L).color = [0 0.6 0.298];
        case 'C'
            plot_structendL(L).color = [1 1 0];
    end

end

% Right
for R = 1:termendcountR
    plot_structendR(R).coord = [coordsend(end+1-R,1),...
                              coordsend(end+1-R,2)];
    plot_structendR(R).nbase = P.doublet1{1,end-termendcountR+R}(2,1);
    plot_structendR(R).bond = 0;
    
    switch plot_structendR(R).nbase
        case 'A'
            plot_structendR(R).color = [1 0.502 0];
        case 'T'
            plot_structendR(R).color = [0.0549 0.3020 0.5725];
        case 'G'
            plot_structendR(R).color = [0 0.6 0.298];
        case 'C'
            plot_structendR(R).color = [1 1 0];
    end
    
end


%% Plotting of DNA strand

hold on

% Plotting of Phosphate backbone

plot_structendL = fliplr(plot_structendL);

phos_coordL = cell2mat({plot_struct0L.coord,...
                        plot_structL.coord,...
                        plot_structendL.coord}');
phos_coordR = cell2mat({plot_struct0R.coord,...
                        plot_structR.coord,...
                        plot_structendR.coord}');

plot(phos_coordL(:,1),phos_coordL(:,2),'LineWidth',3,...
     'Color',[0.3 0.3 0.3])
plot(phos_coordR(:,1),phos_coordR(:,2),'LineWidth',3,...
     'Color',[0.3 0.3 0.3])


axis equal
ax = gca;
axwid = ax.XLim(2)-ax.XLim(1);
sfac = 100/axwid;
set(ax,'XLim',[ax.XLim(1)-sfac/15,ax.XLim(2)+sfac/15]);
set(ax,'YLim',[ax.YLim(1)-sfac/15,ax.YLim(2)+sfac/15]);

plot(phos_coordL(:,1),phos_coordL(:,2),'LineWidth',0.4*sfac,...
     'Color',[0.3 0.3 0.3])
plot(phos_coordR(:,1),phos_coordR(:,2),'LineWidth',0.4*sfac,...
     'Color',[0.3 0.3 0.3])

% Plotting H-bonding

parfor b = 1:size(plot_structL,2)
    if plot_structL(b).bond == 1
        plot([plot_structL(b).coord(1);plot_structR(b).coord(1)],...
             [plot_structL(b).coord(2);plot_structR(b).coord(2)],...
             'LineWidth',0.4*sfac,...
             'Color',[0.3 0.3 0.3])
    end
end

for i = 1:size(plot_structL,2)
    
    switch plot_structL(i).nbase
        case 'A'
            plot_structL(i).color = [1 0.502 0];
        case 'T'
            plot_structL(i).color = [0.0549 0.3020 0.5725];
        case 'G'
            plot_structL(i).color = [0 0.6 0.298];
        case 'C'
            plot_structL(i).color = [1 1 0];
    end
    
    switch plot_structR(i).nbase
        case 'A'
            plot_structR(i).color = [1 0.502 0];
        case 'T'
            plot_structR(i).color = [0.0549 0.3020 0.5725];
        case 'G'
            plot_structR(i).color = [0 0.6 0.298];
        case 'C'
            plot_structR(i).color = [1 1 0];
    end
    
    
    scatter(plot_structL(i).coord(1),plot_structL(i).coord(2),40*sfac,...
        plot_structL(i).color,'filled','MarkerEdgeColor',[0.3 0.3 0.3],...
        'LineWidth',0.25*sfac)
	scatter(plot_structR(i).coord(1),plot_structR(i).coord(2),40*sfac,...
        plot_structR(i).color,'filled','MarkerEdgeColor',[0.3 0.3 0.3],...
        'LineWidth',0.25*sfac)
    
end

% plotting termstrand_0

% Left
for L = 1:term0countL
    scatter(plot_struct0L(L).coord(1),plot_struct0L(L).coord(2),40*sfac,...
        plot_struct0L(L).color,'filled','MarkerEdgeColor',[0.3 0.3 0.3],...
        'LineWidth',0.25*sfac)
end

% Right
for R = 1:term0countR
    scatter(plot_struct0R(R).coord(1),plot_struct0R(R).coord(2),40*sfac,...
        plot_struct0R(R).color,'filled','MarkerEdgeColor',[0.3 0.3 0.3],...
        'LineWidth',0.25*sfac)
    
end


% plotting termstrand_end

% Left
for L = 1:termendcountL
    scatter(plot_structendL(L).coord(1),plot_structendL(L).coord(2),40*sfac,...
        plot_structendL(L).color,'filled','MarkerEdgeColor',[0.3 0.3 0.3],...
        'LineWidth',0.25*sfac)
end

% Right
for R = 1:termendcountR
    scatter(plot_structendR(R).coord(1),plot_structendR(R).coord(2),40*sfac,...
        plot_structendR(R).color,'filled','MarkerEdgeColor',[0.3 0.3 0.3],...
        'LineWidth',0.25*sfac)
end




g = gscatter([NaN,NaN,NaN,NaN],...
    [NaN,NaN,NaN,NaN],['A';'T';'G';'C']...
    ,[1 0.502 0;0.0549 0.3020 0.5725;0 0.6 0.298;1 1 0],'.',50);
lgd = legend('show');

text(coords0(end-(term0countL),1),coords0(end-(term0countL),2),...
    '5''','HorizontalAlignment','Center',...
    'Fontsize',3*sfac);
text(coords0(term0countR+3,1),coords0(term0countR+3,2),...
    '3''','HorizontalAlignment','Center',...
    'Fontsize',3*sfac);
text(coordsend(termendcountL+3,1),coordsend(termendcountL+3,2),...
    '3''','HorizontalAlignment','Center',...
    'Fontsize',3*sfac);
text(coordsend(end-termendcountR,1),coordsend(end-termendcountR,2),...
    '5''','HorizontalAlignment','Center',...
    'Fontsize',3*sfac);

t = text(0,coords0(end-(term0countL),2)-2,...
    ['ΔG: ',num2str(Energy),' kcal/mol'],'HorizontalAlignment'...
    ,'center','Fontsize',20);
t.Units = 'normalized';
t.Position = [0.5 -0.1 0];

lgd.Units = 'normalized';
lgd.FontSize = 20;

axis off
f = getframe(gcf);


end