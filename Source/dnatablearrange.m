function [f,g] = dnatablearrange(aaa,bbb,ccc,ddd,eee)

scafftab0 = zeros(5);
scafftab = scafftab0 + 0.0001;
scafftab(5,:) = 0;
scafftab(:,5) = 0;

for i = 1:5
    for j = 1:5
        if i == 5 || j == 5
            dnaentable{i,j} = scafftab0;
            dnaentableterm{i,j} = scafftab0;
        else
            dnaentable{i,j} = scafftab;
            dnaentableterm{i,j} = scafftab;
        end
    end
end

dnaentable{4,1}(1:4,1:4) = aaa{4,2};
dnaentable{3,2}(1:4,1:4) = aaa{1,2};
dnaentable{2,3}(1:4,1:4) = aaa{2,2};
dnaentable{1,4}(1:4,1:4) = aaa{3,2};

dnaentableterm{4,1}(1:4,1:4) = bbb{4,2};
dnaentableterm{3,2}(1:4,1:4) = bbb{1,2};
dnaentableterm{2,3}(1:4,1:4) = bbb{2,2};
dnaentableterm{1,4}(1:4,1:4) = bbb{3,2};

NNwccoordtab = [ccc{1,2} ccc{5,2}  ccc{6,2}  ccc{2,2};...
                ccc{4,2} ccc{10,2} ccc{8,2}  ccc{6,2};...
                ccc{7,2} ccc{9,2}  ccc{10,2} ccc{5,2};...
                ccc{3,2} ccc{7,2}  ccc{4,2}  ccc{1,2};];

for i = 1:4
    for j = 1:4
        dnaentable{i,5-i}(j,5-j) = NNwccoordtab(i,j);
        dnaentableterm{i,5-i}(j,5-j) = NNwccoordtab(i,j);
    end
end

for l = 1:4
    for m = 1:4
        for n = 1:4
            dnaentable{m,n}(l,5-l) = dnaentable{5-l,l}(n,m);
            dnaentableterm{m,n}(l,5-l) = dnaentableterm{5-l,l}(n,m);
        end
    end
end

for i = 1:4
    for j = 1:4
        dnaentableterm{i,5}(j,5-j) = ddd{j,2}(i);
        dnaentableterm{5,i}(j,5-j) = eee{5-j,2}(i);
        dnaentableterm{i,5-i}(5,j) = ddd{5-i,2}(j);
        dnaentableterm{i,5-i}(j,5) = eee{i,2}(j);
    end
end

dnaentable(end+1,1:2) = ccc(11:12,2);

f = dnaentable;
g = dnaentableterm;

end




