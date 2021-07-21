function f = Keqmat(tempseq,primercombi,revprim,spectemp,Mgconc)

pairingmat = zeros(size(primercombi,2)+2);

s = 1;

seqvector = {tempseq,primercombi{:},revprim};

for i = 1:size(primercombi,2)+2
    for j = s:size(primercombi,2)+2
        if i ~= 1 || j ~= 1
            pairingmat(i,j) = mfethermohyb(seqvector{1,i},seqvector{1,j},spectemp,Mgconc);
        end
    end
    s = s+1;
end

pairingmat;
pairingmatn = pairingmat + pairingmat';
pairingmatn(1:size(pairingmatn,1)+1:end) = diag(pairingmat);

Keqmatrix = exp(-pairingmatn./(1.987e-3*(273.15+spectemp)));
Keqmatrix(1,1) = 0;


f = Keqmatrix;


end