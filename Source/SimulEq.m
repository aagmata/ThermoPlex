function [f,g] = SimulEq(T,multiplexset,positionvector,revprim,template,Mgconc,c0vector)


% Function that simulates equilibrium scenario in the annealing step of a
% multiplex PCR assay by predicting multi-reaction equilibrium product
% distribution curves as a function of temperature. The scenario is
% simulated given one target template and the set of multiplex PCR primers
% amplifying the target. 
% 
% Input:
%   multiplexset - the set of multiplex PCR primers to be simulated
%   positionvector - most probable binding position of the multiplex PCR
%                    primers in the target template. Can be set to [] if
%                    inhibition scenario due to polymerase extension 
%                    blockage is disregarded
%   revprim - reverse primer sequence
%   template - sequence of the target template to be amplified
%   Mgconc - Magnesium Chloride (MgCl2) concentration of the solution under
%            which the simulation will be carried on
%   c0vector - row vector of initial primer concentration in uM with order: 
%              [<template> ----<forward primers>---- <reverse primer>]
%   fracspecthresh - fractional specificity threshold which is the ratio of
%                    target to non-target fractional converstion
% 
% Output:
%   f - Product concentration given multiplexed primers with the given 
%       template
%   g - Product concentration given single-plex primer with the given
%       template

% =========================================================================
%% Main Program
% =========================================================================
%

%% Initiation
tic;
xbar=waitbar(0,'Initiating parameters');
titleHandle = get(findobj(xbar,'Type','axes'),'Title');
set(titleHandle,'FontSize',9)
% Consider polymerase extension blockage or not
if isempty(positionvector)                                                                           
    c0vector = c0vector;
else
    c0temp(1:length(multiplexset)) = c0vector(1);
    c0vector = [c0temp,c0vector(2:end)];
end

% Initiate vector length parameters
li = length(c0vector);
lk = factorial(li+1)/(factorial(li-1)*2);
lw = (2*li)+lk;
np = length(multiplexset);

% Product coordinate upper triangular matrix 
for k = 1:lk
    l = li - floor(sqrt(-8*(k-1) + 4*(li+1)*(li)-7)/2.0 - 0.5);
    m = k + l - 1 - (li+1)*(li)/2 + (li-l+2)*((li-l+2)-1)/2;
    upmat(l,m) = k;

end
prodcoord = diag(upmat,length(c0temp));
prodcoord = prodcoord(1:np);

%% Calculation of pairwise interaction matrices per temperature increment 
fprintf('\n')
fprintf('\n')
fprintf('Preparing pairwise interaction matrices, patience...')
fprintf('\n')
fprintf('\n')

for i = 1:length(T)
    waitbar(i/length(T)*0.9,xbar,strcat('Calculating pairwise interaction matrices, patience (',...
        string(floor(i/length(T)*90)),'%)'));
    spectemp = T(i);
    dGvector(i,:) = dGvectcomp(template,...
                multiplexset,positionvector,revprim,spectemp,Mgconc);
end

%% Multi-reaction equilibrium thermodynamics

fprintf('Solving system of equations per temperature increment...')
fprintf('\n')
fprintf('\n')

for i = 1:length(T)
    waitbar(0.9+i/length(T)*0.1,xbar,strcat('Solving equilibrium system of equations (',...
        string(90+floor(i/length(T)*10)),'%)'));
    spectemp = T(i);
    eqnKeq = eqnsetKeq(c0vector,dGvector(i,:),spectemp);
    reactconcKeq = NewtonR(eqnKeq,c0vector,0);
    prodconcKeq = prodsolKeq(reactconcKeq,dGvector(i,:),40);
    initguess = [sqrt(abs(reactconcKeq').*1e6),sqrt(abs(prodconcKeq').*1e6),ones(1,li)];
    
    eqnsys = eqnsetLUM(c0vector,dGvector(i,:),spectemp);
    solution = NewtonR(eqnsys,initguess,0)';
    solsq = solution(li+1:li+lk).^2;
    targetprod = solsq(prodcoord);
    ProdconcT(i,:) = targetprod./c0vector(1);
    
    yyp = c0vector(length(c0temp)+1).*1e-6;
    ct = c0vector(1).*1e-6;
    for j = 1:np
        Xt = @(G,T,Cr,Ct) 0.5.*((1+((Cr+exp((G.*1000)./...
            (1.987.*(273.15+T))))/Ct))-sqrt((1+((Cr+exp((G.*1000)./...
            (1.987.*(273.15+T))))/Ct)).^2-(4.*(Cr/Ct))));
        idealX(1,j) = Xt(dGvector(i,prodcoord(j)),spectemp,yyp,ct);
        [dGvector(i,prodcoord(j)),spectemp];
    end
    idealconc(i,:) = idealX(:);
    

end

% Template product concentrations
Prodconcfinal = ProdconcT;

for i = 1:np-1
    Prodconcfinal(:,np-i) = ProdconcT(:,np-i).*...
        (prod(1-ProdconcT(:,np-i+1:end),2));
end

f = Prodconcfinal;
g = idealconc;
close(xbar)
end