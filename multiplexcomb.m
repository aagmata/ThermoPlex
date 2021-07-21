function [f,g,h] = multiplexcomb(primcandidatemeta,ampliconresolution)


Positionarray = {};
Vectorarray = {};
primcandidatemeta(1,:) = [];
nsp = size(primcandidatemeta,2);

for i = 1:nsp
    Primerarray = primcandidatemeta{1,i};
    Primerarray = Primerarray(2:end,1);
    Primerarray2(:,i) = {Primerarray};
    Positionarray(:,i) = {cell2mat(primcandidatemeta{1,i}(2:end,end))};
    Vectorarray(:,i) = {[1:size(Positionarray{i},1)]};
end

Combimat = combvec(Vectorarray{:})';
mataddress = 0;
Truecombimat = [];
Trueposmat = [];
% Producing the combination matrix

for i = 1:size(Combimat,1)                              %   Primer index
    for j = 1:size(Combimat,2)                          %   Sample index
        Positioneval(1,j) = Positionarray{j}(Combimat(i,j));
    end
    for k = 1:length(Positioneval)
        evalmat(k,:) = abs(Positioneval-Positioneval(k))<=ampliconresolution;%   Amplicon length resolution threshold
    end
    evalsum = sum(sum(evalmat')');

    if evalsum == nsp
        mataddress = mataddress+1;
        Truecombimat(mataddress,:) = Combimat(i,:);
        Trueposmat(mataddress,:) = Positioneval;
    end
    
    evalsum = 0;
end

% Primer sequence combination matrix
Primcombiarray = {};

if isempty(Truecombimat)
   f = [];
   g = [];
   h = [];
   return
end
for i = 1:size(Truecombimat,1)                              %   Primer index
    for j = 1:size(Truecombimat,2)                          %   Sample index
        Primcombiarray(i,j) = Primerarray2{j}(Truecombimat(i,j));
    end
end

f = Primcombiarray;
g = Trueposmat;
h = Truecombimat;

end