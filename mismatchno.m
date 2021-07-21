function f = mismatchno(Seq1,Seq2)

mismatch = 0;
[a,b] = nwalign(Seq1,Seq2);

Seq1 = b(1,:);
Seq2 = b(3,:);

for i = 1:length(Seq1)
    if Seq1(i) ~= '-' && Seq2(i) ~= '-' && (Seq1(i) ~= Seq2(i))
        mismatch = mismatch + 1;
    end
end
f = mismatch;