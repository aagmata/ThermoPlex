function [S, E, M] = rnafold(seq, varargin)
%RNAFOLD predicts minimum free energy secondary structure of RNA sequence
%
%   S = RNAFOLD(SEQ) returns the minimum free energy secondary structure S
%   in bracket notation.
%
%   [S, E] = RNAFOLD(SEQ) returns the minimum free energy value E (in
%   kcal/mol) for predicted secondary structure S.
%
%   [S, E, M] = RNAFOLD(SEQ) returns the minimum free energy secondary
%   structure as a connectivity matrix M. M is an upper triangular matrix
%   where M(i,j) = 1 if and only if the i-th residue in SEQ is paired with
%   the j-th residue of RNA sequence SEQ.
% 
%   RNAFOLD(SEQ, ..., 'MINLOOPSIZE', M) specifies the minimum size of the
%   substructure (in bases) to be considered when computing the free
%   energy. Default is 3.
%
%   RNAFOLD(SEQ, ..., 'NOGU', true) specifies whether GU or UG pairs are
%   forbidden to form (true) or not (false). Default is false. 
%
%   RNAFOLD(SEQ, ..., 'PROGRESS', true) specifies whether the folding
%   progress bar is displayed or not while minimum free energy secondary
%   structure is being computed. Default is false.
%
%   Examples:
%       % Find minimum free energy secondary structure 
%       seq = 'ACCCCCUCCUUCCUUGGAUCAAGGGGCUCAA';
%       [ss, energy] = rnafold(seq)
%
%   See also RNACONVERT, RNADEMO, RNAPLOT.

%   References: 
%   [1] Wuchty, S., Fontana, W., Hofacker I., Schuster P., Biopolymers (1999),
%   49:145-165.  
%   [2] Mathews, D., Sabina, J., Zuker, M., Turner, D., J. Mol. Biol. (1999)
%   288:911-940.

%   Copyright 2007-2012 The MathWorks, Inc.



noClosingGU = false; % if true, forbids GU closing pairs
m = 3;               % minimum size of loop [i,j]
progressBar = false; % if true, folding progress bar is displayed until folding is complete

%--------------------------------------------------------------------------
%  input checking
%--------------------------------------------------------------------------

% if a structure, try to extract the sequence
if isstruct(seq)
    seq = bioinfoprivate.seqfromstruct(seq);
end
seq = upper(seq);


%check that sequence is not empty
if isempty(seq)
    error(message('bioinfo:rnafold:EmptySequence'));
end

% check that seq is a vector 1xn
if isvector(seq)
    if size(seq, 1) ~= 1
        seq = seq';
    end
else
    error(message('bioinfo:rnafold:InputNotVector'));
end

% convert into numeric sequence and check validity of sequence
if ~isnumeric(seq)
    if any(~isletter(seq)) || ~bioinfoprivate.isnt(seq, 'ACGTUOnly', true)
        error(message('bioinfo:rnafold:InvalidInputNtSequence'));
    else
        s = nt2int(seq);
    end
else 
    s = seq;
    if (~all(s==1|s==2|s==3|s==4))
        error(message('bioinfo:rnafold:InvalidInputNumSequence'));
    end
end
   
% check arguments
if  nargin > 1
    if rem(nargin,2) == 0
        error(message('bioinfo:rnafold:IncorrectNumberOfArguments', mfilename));
    end
    okargs = {'noGU', 'minloopsize', 'progress'};
    
    for j=1:2:nargin-2
        pname = varargin{j};
        pval = varargin{j+1};
        k = find(strcmpi(pname, okargs));
        if isempty(k)
            error(message('bioinfo:rnafold:UnknownParameterName', pname));
        elseif length(k)>1
            error(message('bioinfo:rnafold:AmbiguousParameterName', pname));
        else
            switch(k)
                case 1 % noGU
                    if isempty(bioinfoprivate.opttf(pval))
                        error(message('bioinfo:rnafold:InvalidNoGUParameter', num2str( pval ), pname));
                    elseif (isnumeric(pval) && (~isscalar(pval) || (pval ~= 0 && pval ~= 1)))
                        error(message('bioinfo:rnafold:InvalidNoGUNumParameter', num2str( pval ), pname));
                    else 
                        noClosingGU = pval;
                    end
                case 2 % minloopsize
                    if (isnumeric(pval) && isscalar(pval) && pval > 0  && pval < length(s))
                            m = pval;
                    else
                        error(message('bioinfo:rnafold:InvalidMinLoopParameter'));
                    end
                case 3 % progress bar
                    if isempty(bioinfoprivate.opttf(pval))
                        error(message('bioinfo:rnafold:InvalidProgressParameter', num2str( pval ), pname));
                    elseif (isnumeric(pval) && (~isscalar(pval)) || (pval ~= 0 && pval ~= 1))
                        error(message('bioinfo:rnafold:InvalidProgressNumParameter', num2str( pval ), pname));
                    else
                        progressBar = bioinfoprivate.opttf(pval);
                    end
            end
        end
    end
end
    
%--------------------------------------------------------------------------
% load free energy data
%--------------------------------------------------------------------------
[internal37, bulge37, hairpin37, tetra, tetraloop37, stack37] = rnafoldenergies;
[tstacki37, tstackh37] = rnaterminalstack;
int1x1loop37 = rnaint1x1;
int2x1loop37 = rnaint2x1;
int2x2loop37 = rnaint2x2;
[dangling3, dangling5] = rnadangle;

%--------------------------------------------------------------------------
% setup constants
%--------------------------------------------------------------------------
% base pair types (AU=1; CG=2; GC=3; UA=4; GU=5; UG=6)
bpairs=[ 0 0 0 1; 0 0 2 0; 0 3 0 5; 4 0 6 0 ];
    
logc = 1.0785; %log constant for energy extrapolation, logc=RT*1.75, R=1.987 and T=temp(310.15K, 37C)   
terminalAU = 0.50; % terminal AU penalty for bulges > 1 base 
ninio = [0.5 0.5 0.5 0.5]; maxNinio = 3; % Ninio correction
bonusGGG = -2.2; % bonus for GGG hairpin
bonusC = [0.3; 1.6; 1.4];% penalty for poly-C tetraloops 
Mc = 3.4; Mi = 0.4; Mb = 0; % constant for multiloop decomposition

%--------------------------------------------------------------------------
%  initialization
%--------------------------------------------------------------------------
seqLength=length(s); 
C = ones(seqLength,seqLength)* inf;   % C(i,j) = min free energy for substring [i,j] with i<->j (closing)
FM = ones(seqLength,seqLength)* inf;  % FM(i,j) = min free energy of multi-loop for substring [i,j]
FM1 = ones(seqLength,seqLength)* inf; % FM1(i,j) = min free energy of rightmost stem of multiloop for substring [i,j]
F = ones(seqLength,1)* inf;           % F(j) = minimum free energy for substring [1,j] 

% auxiliary indices for backtracking 
CIndex = ones(seqLength, seqLength, 'single')* -1; 
CType = ones(seqLength, seqLength, 'single')* -1;
FMIndex = ones(seqLength, seqLength, 'single')* -1;
FM1Index = ones(seqLength, seqLength, 'single')* -1; 
FIndex = ones(seqLength, 1, 'single')* -1;
FMType = ones(seqLength, 1, 'single')* -1;
FType = ones(seqLength, 1, 'single')* 2;

%--------------------------------------------------------------------------
% precompute frequently used values
%--------------------------------------------------------------------------
% tables with propensity to pair and canonical pairs
[CP, TP] = findPairs(s, noClosingGU); % CP(i,j)=1 iff i<-->j, TP(i,j)=k iff ij = k

% energies for generic loops and bulges (extrapolated where necessary)
loopTable = [internal37(1:30)' internal37(30) + 1.079 * log((31:seqLength)/30)]; 
bulgeTable = [bulge37(1:30)' bulge37(30) + 1.079 * log((31:seqLength)/30)];


% dangling end contributions
D5 = zeros(1, seqLength); % D5(i,j)= 5' dangling end, i<-->j, i-1
D3 = zeros(1, seqLength); % D3(i,j)= 3' dangling end, i<-->j, j+1
for j = 2+m+1:seqLength
    D5(2:j-m-1,j) = dangling5(s(j), s(1:j-m-2) + 4 * (s(2:j-m-1)-1));
end
for j = 1+m+1:seqLength-1
    D3(1:j-m-1,j) = dangling3(s(j), s(j+1) + 4 * (s(1:j-m-1)-1));
end

%--------------------------------------------------------------------------
% main DP
%--------------------------------------------------------------------------
 
if progressBar 
    fprintf('Folding in progress:')
end

for i = seqLength-m-1:-1:1
    
    % progress bar, useful for long computations
    if progressBar
        if(mod(i,10)==0)
            fprintf('.')
        end
        if (mod(i,300)==0)
            fprintf('\n')
        end
    end
    
    for j = i+m+1:seqLength
        ij = TP(i,j);
        if ij % i <--> j (base-pair)
            if (ij > 4 && noClosingGU)
                C(i,j) = Inf; 
            else
                [minLoopEnergy, myp, myq] = getMinInteriorEnergy(i,j);
                minHairpinEnergy = hairpinEnergy(i,j);
                [minClosedMLoopEnergy, myk] = getMinClosedMLoopEnergy(i,j);
                [C(i,j), CType(i,j)] = min([minHairpinEnergy minLoopEnergy minClosedMLoopEnergy]); 
                
                if (CType(i,j) == 2)
                    CIndex(i,j)= myp; 
                    CIndex(j,i)= myq;
                elseif (CType(i,j) == 3)
                    CIndex(i,j) = myk;
                end
            end
        end
        [FM1(i,j), FM1Index(i,j)] = getMinStemLoopEnergy(i,j);
        [FM(i,j), FMIndex(i,j), FMType(i,j)] = getMinMLoopEnergy(i,j);
    end
end

%--------------------------------------------------------------------------
% Min free energy for segment [1,j]
%--------------------------------------------------------------------------

F(1:m,1) = 0; 
for j = m+1:seqLength
    [F(j) FIndex(j)] = getMinClosed5PrimeEnergy(j);
    if F(j) >= F(j-1)
        F(j) = F(j-1);
        FIndex(j) = j-1;
        FType(j) = 1;
    end
end

%--------------------------------------------------------------------------
% Traceback
%--------------------------------------------------------------------------

M = zeros(seqLength, seqLength); 
S = '';
S(1:seqLength)='.';
topstack = 1;
bstack(1,:) = [1 seqLength 1];

while (topstack)
    i = bstack(topstack,1); %starting index of substring to consider
    j = bstack(topstack,2); %ending index of substring to consider
    t = bstack(topstack,3); %type of array to consider in backtracking 
    topstack = topstack-1;
    
    if(j~=i && i>0 && j>0)
        topstack = topstack+1;
        switch (t)
            case 1 % backtrack in array F
                if (FType(j) == 1)
                    bstack(topstack,1:3) = [1 j-1 1]; % 1,j-1,F
                else % backtrack in C and F
                    bstack(topstack,1:3) = [FIndex(j) j 2]; % l,j,C
                    if (FIndex(j)>1)
                        topstack = topstack+1;
                        bstack(topstack,1:3) = [1 FIndex(j)-1 1]; % 1,j-1,F
                    end
                end
            case 2 % backtrack in array C
                S(i)='('; S(j)=')';
                M(i,j)=1; M(j,i)=1;
                if (CType(i,j)==2) % interior loop
                    bstack(topstack,1:3) = [CIndex(i,j) CIndex(j,i) 2]; %p q C
                elseif (CType(i,j)==3) % multiloop
                    bstack(topstack,1:3) = [i+1 CIndex(i,j)-1 4]; %i+1,k-1,FM
                    topstack = topstack+1;
                    bstack(topstack,1:3) = [CIndex(i,j) j-1 3];   %k,j-1,FM1
                else 
                    topstack = topstack - 1;
                end
            case 3 % backtrack in array FM1
                bstack(topstack,1:3) = [i FM1Index(i,j) 2]; % i,l,C
            case 4 % backtrack in array FM 
                if (FMType(i,j)==1)
                    bstack(topstack,1:3) = [i FMIndex(i,j)-1 4]; %i,k-1,FM
                    topstack = topstack+1;
                    bstack(topstack,1:3) = [FMIndex(i,j) j 3];   %k,j,FM1
                else
                    bstack(topstack,1:3) = [FMIndex(i,j) j 3];
                end
            otherwise
                error(message('bioinfo:rnafold:IncorrectBackTrackingType'));
        end
    end
end


if nargout ~= 0
   E = F(seqLength); % min free energy of [1,seqLength] in kcal/mol
   M = triu(M); % (upper triangular) connectivity matrix
end 

%-------------------------------------------------------------------------
% nested functions
%-------------------------------------------------------------------------

%--------------------------------------------------------------------------
% MINIMUM FREE ENERGY for 5' FRAGMENT 
%--------------------------------------------------------------------------   
function [minClosedF5 hstar] = getMinClosed5PrimeEnergy(j)
% Min free energy of fragment [s(1)...s(j)] with j paired 

dangleEnergy = zeros(j-m-2,1);
d3 = 0;
dangleEnergy = dangleEnergy + D5(2:j-m-1,j);
if (j < seqLength)
    dangleEnergy = dangleEnergy + D3(2:j-m-1,j);
    d3 = D3(1,j);
end
[minClosedF5 hstar] = min([C(1,j)+ d3; F(1:j-m-2) + C(2:j-m-1,j) + dangleEnergy]); 

if minClosedF5 == Inf
    hstar = j; 
end
end % end getMinClosed5PrimeEnergy


%--------------------------------------------------------------------------
% MINIMUM FREE ENERGY of MULTILOOP 
%--------------------------------------------------------------------------
function [minMLoopEnergy rstar type] = getMinMLoopEnergy(i,j)
% Min free energy for multiloop [i,j] consisting of a rightmost
% single stem and another multiloop component. i and j are not assumed to
% pair.

[minMLoopEnergy1, r1x] = min(FM(i,(i+m+1:j-m-1)-1)' + FM1((i+m+1:j-m-1),j));
[minMLoopEnergy2, r2x] = min(FM1((i:j-m-1),j) + Mb*((i:j-m-1)-i)');
if minMLoopEnergy1 < minMLoopEnergy2
    minMLoopEnergy = minMLoopEnergy1;
    rstar = r1x + i + m;
    type = 1;
else
    minMLoopEnergy = minMLoopEnergy2;
    rstar = r2x + i - 1;
    type = 2;
end

if isempty(minMLoopEnergy) || minMLoopEnergy == Inf
   minMLoopEnergy = Inf;
   rstar = i; 
   type = 0;
end

end % end function getMinMLoopEnergy

%--------------------------------------------------------------------------
% MINIMUM FREE ENERGY for SINGLE-STEM MULTILOOP 
%--------------------------------------------------------------------------
function [minStemLoopEnergy lstar] = getMinStemLoopEnergy(i,j)
% Min free energy for multiloop [i,j] consisting of a single stem plus
% unpaired bases at the right side. i and j are not assumed to pair.

dangleEnergy = zeros(1,j-i-m); % for dangling ends

energy = C(i,(i+m+1:j)) + Mb*(j-(i+m+1:j)) + Mi;
if isempty(energy(CP(i,(i+m+1:j))))
    minStemLoopEnergy = Inf;
    lstar = j;
else
    penaltyAU = (TP(i,(i+m+1:j))< 2 | TP(i,(i+m+1:j)) > 3) * terminalAU; 
    penaltyAU = penaltyAU .* CP(i,(i+m+1:j)); % apply only to possible canonical bps
    if i > 1
        dangleEnergy = dangleEnergy + D5(i,(i+m+1:j));
    end
    if j < seqLength
        dangleEnergy = dangleEnergy + D3(i,(i+m+1:j));
    end
    energy = energy .* CP(i,(i+m+1:j)) + penaltyAU + dangleEnergy;

    [minStemLoopEnergy, lstar] = min(energy);
    if minStemLoopEnergy == Inf
        lstar = j;
    else
        lstar = lstar + i + m;
    end
end
end % end minStemLoopEnergy function

%--------------------------------------------------------------------------
% MINIMUM FREE ENERGY for CLOSED MULTILOOP 
%--------------------------------------------------------------------------

function [minClosedMLoopEnergy kstar] = getMinClosedMLoopEnergy(i,j)
% Min free energy of multiloop structure [i,j] with i<->j (pairing) and
% kstar as branching base.

dangleEnergy = 0;
if i > 0
    dangleEnergy  = dangleEnergy + D5(i,j);
end
if (j <seqLength)
    dangleEnergy = dangleEnergy + D3(i,j);
end

[minClosedMLoopEnergy, k] = min(FM(i+1,(i+1:j-m-2)-1)' + FM1((i+1:j-m-2),j-1) + ...
   + Mc + dangleEnergy); %FIXME (see formula for 16)

if (isempty(minClosedMLoopEnergy) || minClosedMLoopEnergy == Inf)
    minClosedMLoopEnergy = Inf;
    kstar = j;
else
    kstar = k + i;
end
end % end getMinClosedMLoopEnergy

%--------------------------------------------------------------------------
% MINIMUM FREE ENERGY FOR INTERIOR LOOP 
%--------------------------------------------------------------------------

function [minLoopEnergy, pstar, qstar] = getMinInteriorEnergy(i,j)
% Min free energy of loop structure [i,j] with i<->j and pstar<-->qstar

%=== Special cases 
energy = ones(1,8) * Inf;
allBest = zeros(1,9);

if CP(i+1,j-1) % stacked basepairs
    energy(1) = stack37(s(i+1) + (s(i)-1)*4, s(j-1) + (s(j)-1)*4) + C(i+1,j-1);
end
if CP(i+2,j-2) % 1x1 loop
    pq =  bpairs(s(i+2),s(j-2));
    energy(2) = int1x1loop37(s(i+1) + (ij-1)*4, s(j-1) + (pq -1)*4) + C(i+2,j-2);
end
if CP(i+2,j-3) % 1x2 loop
    pq =  bpairs(s(i+2),s(j-3));
    energy(3) = int2x1loop37(s(i+1) + (s(j-3+1)-1)*4 + (ij-1)*16, s(j-1) + ...
        (pq-1)*4) + + C(i+2,j-3);
end
if CP(i+3,j-2) % 2x1 loop
    qp = bpairs(s(j-2),s(i+3));
    ji = bpairs(s(j),s(i));
    energy(4) = int2x1loop37(s(j-1) + (s(i+1)-1)*4 + (qp-1)*16, s(i+3-1) + ...
        (ji-1)*4) + C(i+3,j-2);
end
if CP(i+3,j-3) % 2x2 loop
    pq =  bpairs(s(i+3),s(j-3));
    energy(5) = int2x2loop37(s(j-1)+ (s(i+1)-1)*4 + (pq-1)*16 +(ij-1)* 96) + ...
        + C(i+3,j-3);
end
   
%=== Bulges 
q = find(CP(i+1,i+m+2:j-2)); % valid candidates q (relative indexing)
if ~isempty(q)
    sizes = j - i - q - m - 2; % size = j-q-1, q=q+i+m+1 (absolute indexing)
    if (sizes) 
        [energy(6), best] = min(bulgeTable(sizes) + C(i+1,q));
        allBest(6) = q(best)+ i + m + 1;
    end
end
p = find(CP(i+2:j-m-2,j-1)); % valid candidates p (relative indexing)
if ~isempty(p)
    sizes = p; % size = p-i-1, p=p+i+1 in absolute indexing
    if (sizes)
        [energy(7), best] = min(bulgeTable(sizes)'+ C(p,j-1));
        allBest(7) = p(best) + i + 1;
    end
end

%=== Generic loops (do not consider first row and last columns -> bulges)
[p,q] = find (CP(i+2:j-m-2, i+m+2:j-2)); % valid candidates p and q (relative indexing)
if (~isempty(p))
    sizes = j - q - i - m - 2 + p; % num of unpaired bases for all internal loops between i and j
    % size = (p-i-1) + (j-q-1); p=p+i+1, q=q+i+m+1 (in absolute indexing)
    
    ninioDiff = j - q - i - m - 2 - p; % (j-(q+i+m+1)-1)-((p+i+1)-i-1)
    ninioInd = find(ninioDiff < 0);
    if ~isempty(ninioInd)
        ninioDiff2 = p - j + q + i + m + 2;
        ninioDiff(ninioInd) = ninioDiff2(ninioInd);
    end
    %ninioCorr = min(maxNinio, ninioDiff * 0.3);
    ninioCorr = min(maxNinio, ninioDiff .* ninio(2));
    
    % two terminal mismatches: i<-->j and i+1:j-1, p<-->q and p-1:q+1
    termMismatchInd = sub2ind(size(tstacki37), s(q+i+m+2)+ (s(q+i+m+1)-1)*4,  s(p+i) + (s(p+i+1)-1)*4);
    termMismatch = tstacki37(termMismatchInd) + tstacki37(s(i+1) + (s(i)-1)*4, s(j-1) + (s(j)-1)*4);
    matchEnergy = C(sub2ind(size(C), p+i+1, q+i+m+1));

    % energy of generic loop
    [energy(8), best] = min(loopTable(sizes) + ninioCorr' + termMismatch + matchEnergy');
    allBest(8) = p(best) + i + 1;
    allBest(9) = q(best) + i + m + 1;
end

%=== minimum free energy of interior loop structure
[minLoopEnergy, star] = min(energy);
if minLoopEnergy == Inf
    pstar = i;
    qstar = j;
else 
    switch star
        case 1 % stacked bases
            pstar = i+1;
            qstar = j-1;
        case 2 % 1x1 loop
            pstar = i+2;
            qstar = j-2;
        case 3 % 1x2 loop
            pstar = i+2;
            qstar = j-3;
        case 4 % 2x1 loop
            pstar = i+3;
            qstar = j-2;
        case 5 % 2x2 loop
            pstar = i+3;
            qstar = j-3;
        case 6 % bulge on j
            pstar = i+1;
            qstar = allBest(6);
        case 7 % bulge on i
            pstar = allBest(7);
            qstar = j-1;
        case 8
            pstar = allBest(8);
            qstar = allBest(9);
    end
end
end % getMinInteriorEnergy


%---------------------------
% HAIRPIN ENERGY
%---------------------------

function energy = hairpinEnergy(i,j)
% Min free energy of a hairpin loop between i and j

size = j-i-1; 
energy = Inf;
if (size > 2) 
    if (noClosingGU && ij > 4) % GU closing pair forbidden
        energy = Inf;
    else
        if (size > 30) % j-i-1 = number of unpaired bases
            energy = hairpin37(30) + logc * log((size)/30);
        else
            energy = hairpin37(size);
            if ((size == 4) && (~isempty(strmatch(seq(i:j), tetra, 'exact')))) % special tetraloop
                energy = energy + tetraloop37(strmatch(seq(i:j), tetra, 'exact'));
            end
        end
        if (size ~= 3)% terminal mismatch (except for triloops)
            energy = energy + tstackh37(s(i+1) + (s(i)-1)*4, s(j-1) + (s(j)-1)*4); 
        end

        % special hairpins
        if (i>2 && all(s(i-2:i)==3 & s(j)==4)) % GGG hairpin
            energy = energy + bonusGGG;
        elseif (all(s(i+1:j-1)==2)) % poly-C hairpin
            if ( size ==3)
                energy = energy + bonusC(3);
            else
                energy = energy + bonusC(2) + bonusC(1) * size;
            end
        end
    end
end
end % end hairpinEnergy

%-------------------------------------------------------------------------
end % end rnafold function


%--------------------------------------------------------------------------
% subfunctions
%--------------------------------------------------------------------------
   
function [Pairs, PType] = findPairs(seq, flag)
% Determine pairing propensity for each pair of bases (P) and assign
% base-pair type.  Canonical basepairs are considered: AU, UA, CG, GC,
% including GU and UG if flag is false.

cseq = 5 - seq;
[X,Y]= meshgrid(seq, cseq);
Pairs = triu(X == Y);
PType = zeros(size(Pairs));

AU = (X == 4 & Pairs); PType(AU) = 1;
CG = (X == 3 & Pairs); PType(CG) = 2;
GC = (X == 2 & Pairs); PType(GC) = 3;
UA = (X == 1 & Pairs); PType(UA) = 4;

if (~flag)
    Pairs = triu(Pairs | (X + X' == 7));
    UG = X == 3 & (X + X' == 7); PType(UG) = 6;
    GU = X == 4 & (X + X' == 7); PType(GU) = 5;
end

PType = triu(PType);
end

