function [a,b,c,d,e] = dnaenergies(spectemp,Mgconc,fgc,np)

%Parameters for Mg2+ ion correction (Owczarzy et al., 2008)
%Assumed monovalent concentration in PCR reaction buffer (KCl) = 50mM

aa = 3.92e-5*(0.843-(0.352*sqrt(50e-3)*log(50e-3)));
bb = -9.11e-6;
cc = 6.26e-5;
dd = 1.42e-5*(1.279-(4.03e-3*log(50e-3))-(8.03e-3*((log(50e-3))^2)));
ee = -4.82e-4;
ff = 5.25e-4;
gg = 8.31e-5*(0.486-(0.258*log(50e-3))+(5.25e-3*((log(50e-3))^3)));

%Thermodynamic parameters based on SantaLucia et al. (1998)

%Perfectly complementary doublet

NNWCHtab = [-7.6; -7.2; -7.2; -8.5; -8.4; -7.8; -8.2; -10.6; -9.8; -8.0; 0.2; 2.2; 4.4];
NNWCStab = [-21.3; -20.4; -21.3; -22.7; -22.4; -21.0; -22.2; -27.2; -24.4; -19.9; -5.7; 6.9; 13.8];
NNWCtabcomp = ((NNWCHtab.*1000) - (273.15+spectemp).*NNWCStab)/1000;
corrfac = NNWCHtab.*((273.15+spectemp).*(aa+(bb.*log(Mgconc))+(fgc.*(cc+(dd.*log(Mgconc))))+((ee+(ff.*log(Mgconc))+(gg.*(log(Mgconc).^2)))./(2.*(20-1)))));

% corrfac = 0; %for testing purposes  %%%%%%%%%%%%%%%

NNWCtabcomp = NNWCtabcomp-corrfac;
NNWCtab = {['AA';'TT'] NNWCtabcomp(1);    
           ['AT';'TA'] NNWCtabcomp(2);
           ['TA';'AT'] NNWCtabcomp(3);
           ['CA';'GT'] NNWCtabcomp(4);
           ['GT';'CA'] NNWCtabcomp(5);
           ['CT';'GA'] NNWCtabcomp(6);
           ['GA';'CT'] NNWCtabcomp(7);
           ['CG';'GC'] NNWCtabcomp(8);
           ['GC';'CG'] NNWCtabcomp(9);
           ['GG';'CC'] NNWCtabcomp(10);
           ['I'] NNWCtabcomp(11);
           ['ST'] NNWCtabcomp(12);
           ['DT'] +0.10;};

%single internal mismatch
deltaStab = [-9.8 14.2 -1.0 0; -3.8 8.9 0 5.4; 3.2 0 -15.8 10.4; 0 13.5 -12.3 -8.4;
             -4.2 3.7 -2.3 0; -0.6 -7.2 0 -4.5; -13.2 0 -15.3 -11.7; 0 -6.1 -8.0 -15.8;
             1.7 4.6 -2.3 0; 14.6 -4.4 0 0.2; -2.3 0 -9.5 0.9; 0 -6.2 -8.3 -10.8;
             12.9 8.0 0.7 0; 20.2 16.4 0 0.7; 7.4 0 3.6 -1.7; 0 0.7 -5.3 -1.5];
         
deltaHtab = [-2.9 5.2 -0.6 0; -0.7 3.6 0 2.3; 0.5 0 -6.0 3.3; 0 5.2 -4.4 -2.2;
             -0.9 1.9 -0.7 0; 0.6 -1.5 0 -0.8; -4.0 0 -4.9 -4.1; 0 -1.5 -2.8 -5.0;
             1.2 2.3 -0.6 0; 5.3 0 0 0.7; -0.7 0 -3.1 1.0; 0 -1.2 -2.5 -2.7;
             4.7 3.4 0.7 0; 7.6 6.1 0 1.2; 3.0 0 1.6 -0.1; 0 1.0 -1.4 0.2;];

deltaGtab2 = ((deltaHtab.*1000) - ((273.15+spectemp).*deltaStab))/1000;
corrfac = deltaHtab.*(273.15+spectemp).*((aa+(bb.*log(Mgconc))+(fgc.*(cc+(dd.*log(Mgconc))))+((ee+(ff.*log(Mgconc))+(gg.*(log(Mgconc)).^2))./(np-2))));

% corrfac = 0; %for testing purposes %%%%%%%%%%%%%

deltaGtab2 = deltaGtab2-corrfac;
deltaGtab2((deltaGtab2==-corrfac))=0;

deltaGtab = {['GX';'CY'] deltaGtab2(1:4,1:4);  
             ['CX';'GY'] deltaGtab2(5:8,1:4);
             ['AX';'TY'] deltaGtab2(9:12,1:4);
             ['TX';'AY'] deltaGtab2(13:16,1:4)};

%terminal mismatch
deltaStabterm = [-22.7 -13.9 -11.3 0; -7.1 -10.6 0 -13.5; -11.4 0 0.8 -16; 0 -7.8 -16.1 -21.1;  
             -10.6 -5.9 -9.7 0; -6 -5.1 0 -8.1; -15.4 0 -9.6 -9.2; 0 -10.6 -18.7 -16.9;
             -14.4 -10.6 -11 0; -10.4 -6.2 0 -8.4; -12.5 0 -8.7 -15.4; 0 -12.8 -16.1 -12.8;
             0.2 0.4 1.7 0; 0.4 -0.5 0 5.1; 0.5 0 -5.2 -4.3; 0 5.2 -3.3 -2.4];
         
deltaHtabterm = [-7.9 -5 -4.3 0; -3.2 -3.9 0 -4.9; -4.6 0 -0.7 -5.7; 0 -3 -5.9 -7.4;
             -4.3 -2.6 -3.9 0; -2.7 -2.1 0 -3.2; -6 0 -3.8 -3.8; 0 -3.9 -6.6 -6.1;
             -5.1 -3.6 -3.9 0; -3.8 -2.1 0 -2.9; -4.5 0 -3.1 -5.2; 0 -4.3 -5.5 -4.4;
             0.1 0.3 0.6 0; 0.3 -0.1 0 1.9; 0.2 0 1.5 -1.3; 0 1.9 -1 -0.6];

deltaGtab2term = ((deltaHtabterm.*1000) - ((273.15+spectemp).*deltaStabterm))/1000;
corrfac = deltaHtabterm.*(273.15+spectemp).*((aa+(bb.*log(Mgconc))+(fgc.*(cc+(dd.*log(Mgconc))))+((ee+(ff.*log(Mgconc))+(gg.*(log(Mgconc)).^2))./(np-2))));

% corrfac = 0; %for testing purposes %%%%%%%%%%%%%%%%

deltaGtab2term = deltaGtab2term-corrfac;
deltaGtab2term((deltaGtab2term==-corrfac))=0;
         
deltaGtabterm = {['GX';'CY'] deltaGtab2term(1:4,1:4);   
                 ['CX';'GY'] deltaGtab2term(5:8,1:4);
                 ['AX';'TY'] deltaGtab2term(9:12,1:4);
                 ['TX';'AY'] deltaGtab2term(13:16,1:4)};

         
         



DangendHtab5 = [0.2 0.6 -1.1 -6.9;
                -6.3 -4.4 -5.1 -4.0;
                -3.7 -4.0 -3.9 -4.9;
                -2.9 -4.1 -4.2 -0.2];
            
DangendHtab3 = [-0.5 4.7 -4.1 -3.8;
                -5.9 -2.6 -3.2 -5.2;
                -2.1 -0.2 -3.9 -4.4;
                -0.7 4.4 -1.6 2.9];
            
DangendStab5 = [2.3 3.3 -1.6 -20.0;
                -17.1 -12.6 -14.0 -10.9;
                -10.0 -11.9 -10.9 -13.8;
                -7.6 -13.0 -15.0 -0.5;];

DangendStab3 = [-1.1  14.2 -13.1 -12.6;
                -16.5 -7.4 -10.4 -15.0;
                -3.9 -0.1 -11.2 -13.1;
                -0.8 14.9 -3.6 10.4];

DangendGtab5p = ((DangendHtab5.*1000) - ((273.15+spectemp).*DangendStab5))/1000;
DangendGtab3p = ((DangendHtab3.*1000) - ((273.15+spectemp).*DangendStab3))/1000;
corrfac5 = DangendHtab5.*((273.15+spectemp).*(aa+(bb.*log(Mgconc))+(fgc.*(cc+(dd.*log(Mgconc))))+((ee+(ff.*log(Mgconc))+(gg.*(log(Mgconc).^2)))./(2.*(20-1)))));
corrfac3 = DangendHtab3.*((273.15+spectemp).*(aa+(bb.*log(Mgconc))+(fgc.*(cc+(dd.*log(Mgconc))))+((ee+(ff.*log(Mgconc))+(gg.*(log(Mgconc).^2)))./(2.*(20-1)))));

% corrfac5 = 0; %for testing purposes %%%%%%%%%%%%%%
% corrfac3 = 0; %for testing purposes %%%%%%%%%%%%%%

DangendGtab5 = DangendGtab5p-corrfac5;
DangendGtab3 = DangendGtab3p-corrfac3;

Dangend5 = {['XA';'-T'] DangendGtab5(1,:);     
            ['XC';'-G'] DangendGtab5(2,:);
            ['XG';'-C'] DangendGtab5(3,:);
            ['XT';'-A'] DangendGtab5(4,:)};

Dangend3 = {['AX';'T-'] DangendGtab3(1,:);     %%%%%%%%%%%%%%
            ['CX';'G-'] DangendGtab3(2,:);
            ['GX';'C-'] DangendGtab3(3,:);
            ['TX';'A-'] DangendGtab3(4,:)};

a = deltaGtab;
b = deltaGtabterm;
c = NNWCtab;
d = Dangend5;
e = Dangend3;
        
        
end
            
            
