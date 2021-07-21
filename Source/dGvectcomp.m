function [f,g] = dGvectcomp(tempseq,primercombi,positionvector,revprim,spectemp,Mgconc)

s = 1;

if isempty(positionvector)
    temparray = tempseq;
else
    [positionvector,mm] = sort(positionvector);
    primercombi = primercombi(mm);
    
    for j = 1:size(primercombi,2)        
        bindpos = positionvector(j);
        tempseqnrc = seqrcomplement(tempseq);
        tempseqnrc = tempseqnrc(bindpos-2:bindpos+22);
        tempseqn = seqrcomplement(tempseqnrc);
        temparray{j} = tempseqn;
    end
end

seqvector = {temparray{:},primercombi{:},revprim};
seqvector{:};
pairingmat = zeros(size(seqvector,2));

for i = 1:size(seqvector,2)
    parfor j = s:size(seqvector,2)
        if isempty(positionvector)
            if i == 1 && j ~=1 && j ~= size(primercombi,2)+2
                bindpos = positionvector(j-1);
                tempseqnrc = seqrcomplement(tempseq);
                tempseqnrc = tempseqnrc(bindpos-2:bindpos+25);
                tempseqn = seqrcomplement(tempseqnrc);

                pairingmat(i,j) = ThermoDhyb(seqvector{1,j},tempseqn,...
                    spectemp,Mgconc);
            elseif i == 1 && j == size(primercombi,2)+2
                [a,b] = swalign(seqvector{1,j},tempseq,...
                    'gapopen',20);
                bindpos = strfind(tempseq,b(1,:));
                tempseqnrc = seqrcomplement(tempseq);
                tempseqnrc = tempseqnrc(bindpos-2:bindpos);
                tempseqn = seqrcomplement(tempseqnrc);

                pairingmat(i,j) = ThermoDhyb(seqvector{1,j},tempseqn,...
                    spectemp,Mgconc);
            elseif i ~= 1 || j ~= 1
                pairingmat(i,j) = ThermoDhyb(seqvector{1,j},seqvector{1,i},...
                    spectemp,Mgconc);
            end
        else
            pairingmat(i,j) = ThermoDhyb(seqvector{1,i},seqvector{1,j},...
                spectemp,Mgconc);
        end
        
        
    end

    s = s+1;
end

dGmatrix = pairingmat;


n = length(dGmatrix);
kindex = factorial(n+1)/(factorial(n-1)*2);

for k = 1:kindex
    l = n - floor(sqrt(-8*(k-1) + 4*(n+1)*(n)-7)/2.0 - 0.5);
    m = k + l - 1 - (n+1)*(n)/2 + (n-l+2)*((n-l+2)-1)/2;
    dGvec(k) = dGmatrix(l,m);
end


f = dGvec;
g = primercombi;

end