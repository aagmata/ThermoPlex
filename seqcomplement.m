function c = seqcomplement(sequence)

    complemat = ['ATCG-YRWSKMDVHBXN';...
                 'TAGC-RYWSMKHBDVXN'];
    lseq = length(sequence);
    
    for l = 1:lseq;
        sequencec(l) = complemat(2,sequence(l) == complemat(1,:));
    end
    c = sequencec;
end