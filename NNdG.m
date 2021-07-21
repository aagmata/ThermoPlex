  function [f,i,a] = NNdG(doub1,doub2,energytableint,energytableterm,terminal)

%   Function for calculating Gibbs free energy for two doublet pairs
%   based on the two-state nearest-neighbor thermodynamic parameters 
%   (SantaLucia et al.) with salt correction from Owczarzy et al. (2008).
%   
%
%   Input arguments are doub1 (Sequence 1 in 5' -> 3'); doub2 (Sequence 2 in 3' 
%   -> 5'); spectemp (specific temperature of hybridization); Mgconc
%   (concentration of divalent Mg2+ ions in the solution).
%
%   The output will be the Gibbs free energy for the interaction.
%
%

%Retrieving table values and assigning array elements to corresponding
%varibles

strnumval = ['ACGT-'];
numval = 1:5;

doubstr2num = [numval(doub1(1) == strnumval) numval(doub1(2) == strnumval);...
               numval(doub2(1) == strnumval) numval(doub2(2) == strnumval)];

if terminal == 0
    dG = energytableint{doubstr2num(1),doubstr2num(2)}...
        (doubstr2num(3),doubstr2num(4));
elseif terminal == 1
    dG = energytableterm{doubstr2num(1),doubstr2num(2)}...
        (doubstr2num(3),doubstr2num(4));
end
           
           
f = dG;
i = energytableint{end,1};
a = energytableint{end,2};

end

