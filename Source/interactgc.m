function f = interactgc(Seq1,Seq2)
gc = 0;
[value alignment] =nwalign(Seq1,Seq2);
for p = 1:size(alignment,2)
    if (alignment(1,p)=='G' || alignment(1,p)=='C')
        if alignment(2,p)=='|'
            gc = gc+1;
        end
    end
end
f=gc./length(alignment(2,alignment(2,:) == '|'));