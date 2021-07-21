function f = seqpairplot(Seq1,Seq2,subseqalignment,m,n,mfe)



clf
set(0,'defaultfigurecolor',[1 1 1]);
Seq2r = seqreverse(Seq2);

Seq1 = Seq1(2:end-1);
Seq2 = Seq2(2:end-1);



gapposmat1 = find(subseqalignment(1,:) == '-');
gapposmat2 = find(subseqalignment(3,:) == '-');

gapposmat1(gapposmat1 == 1) = [];
gapposmat1(gapposmat1 == length(subseqalignment(1,:))) = [];
gapposmat2(gapposmat2 == 1) = [];
gapposmat2(gapposmat2 == length(subseqalignment(3,:))) = [];

ngap1 = length(gapposmat1);
ngap2 = length(gapposmat2);

startpt1 = m-length(subseqalignment)+1+ngap1;
startpt2 = n-length(subseqalignment)+1+ngap2;
endpt1 = m;
endpt2 = n;

% directional paths

mismatchcase1 = 0;
gapcase1 = 0;
bulgecase1 = 0;

for k = 1:size(subseqalignment,2)
    if subseqalignment(1,k) == seqcomplement(subseqalignment(3,k))
        pathSeq1(k) = 0;
        mismatchcase1 = 0;
        gapcase1 = 0;
        bulgecase1 = 0;
    elseif subseqalignment(1,k) == '-'
        gapcase1 = gapcase1 + 1;
        pathSeq1((k-gapcase1+1):k) = gapcase1+100;
    elseif subseqalignment(3,k) == '-'
        bulgecase1 = bulgecase1 + 1;
        pathSeq1((k-bulgecase1+1):k) = -bulgecase1;
    else 
        mismatchcase1 = mismatchcase1+1;
        pathSeq1((k-mismatchcase1+1):k) = mismatchcase1; 
    end
end
pathSeq1 = pathSeq1';
if pathSeq1(1) ~= 0
    pathSeq1(1) = [];
    startpt1 = startpt1+1;
end

mismatchcase2 = 0;
gapcase2 = 0;
bulgecase2 = 0;

for k = 1:size(subseqalignment,2)
    if subseqalignment(3,k) == seqcomplement(subseqalignment(1,k))
        pathSeq2(k) = 0;
        mismatchcase2 = 0;
        gapcase2 = 0;
        bulgecase2 = 0;
    elseif subseqalignment(3,k) == '-'
        gapcase2 = gapcase2 + 1;
        pathSeq2((k-gapcase2+1):k) = gapcase2+100;
    elseif subseqalignment(1,k) == '-'
        bulgecase2 = bulgecase2 + 1;
        pathSeq2((k-bulgecase2+1):k) = -bulgecase2;
    else 
        mismatchcase2 = mismatchcase2+1;
        pathSeq2((k-mismatchcase2+1):k) = mismatchcase2;
    end
end
pathSeq2 = pathSeq2';
if pathSeq2(1) ~= 0
    pathSeq2(1) = [];
    startpt2 = startpt2+1;
end

global x1 y1 x2 y2 r h1 k1 h2 k2
options = optimoptions('fsolve','Display','off');

% Substarting coordinates
r = 0.5/(sind((360/(4+(2*(mean([startpt1-1 startpt2-1])))))/2));

if startpt1 > 1
    addvec = [0 -1];
%     r = 0.5/(sind((360/(4+(2*(max([startpt1-1 startpt2-1])))))/2));
    for k = 2:startpt1
        if k == 2
            startseq1coord(1,:) = [-0.5 0];
            perpendvec = cross([0 0 1],[addvec 0]);
            perpendvec = perpendvec(1:2);
            x1 = startseq1coord(k-1,1);
            y1 = startseq1coord(k-1,2);
            x2 = startseq1coord(k-1,1)+perpendvec(1);
            y2 = startseq1coord(k-1,2)+perpendvec(2);
            x0 = [startseq1coord(k-1,:)+addvec];
            centercoord = fsolve(@centercirc,x0,options);
            h1 = centercoord(1);
            k1 = centercoord(2);
            h2 = startseq1coord(k-1,1);
            k2 = startseq1coord(k-1,2);
            startseq1coord(k,:) = fsolve(@coordcirc,x0,options);
            
        else
            x0 = [startseq1coord(k-1,:)+addvec];
            h1 = centercoord(1);
            k1 = centercoord(2);
            h2 = startseq1coord(k-1,1);
            k2 = startseq1coord(k-1,2);
            startseq1coord(k,:) = fsolve(@coordcirc,x0,options);
        end

    end
    startseq1coord(1,:) = [];
    startseq1coord = flipud(startseq1coord);
    
    
else startseq1coord = [];    
end

if startpt2 > 1
    addvec = [0 -1];
    for k = 2:startpt2
        if k == 2
            startseq2coord(1,:) = [0.5 0];
            perpendvec = cross([addvec 0],[0 0 1]);
            perpendvec = perpendvec(1:2);
            x1 = startseq2coord(k-1,1);
            y1 = startseq2coord(k-1,2);
            x2 = startseq2coord(k-1,1)+perpendvec(1);
            y2 = startseq2coord(k-1,2)+perpendvec(2);
            x0 = [startseq2coord(k-1,:)+addvec];
            centercoord = fsolve(@centercirc,x0,options);
            h1 = centercoord(1);
            k1 = centercoord(2);
            h2 = startseq2coord(k-1,1);
            k2 = startseq2coord(k-1,2);
            startseq2coord(k,:) = fsolve(@coordcirc,x0,options);
        else
            x0 = [startseq2coord(k-1,:)+addvec];
            h1 = centercoord(1);
            k1 = centercoord(2);
            h2 = startseq2coord(k-1,1);
            k2 = startseq2coord(k-1,2);
            startseq2coord(k,:) = fsolve(@coordcirc,x0,options);
        end

    end
    startseq2coord(1,:) = [];
    startseq2coord = flipud(startseq2coord);
    

else startseq2coord = [];    
end


% Plot vectors for Seq1

for k = 1:size(pathSeq1)
    if k == 1
        subseq1coord(k,:) = [-0.5 0];
        addvec1 = [0 1];
        seqbondvec1(k,:) = subseq1coord(k,:);
        loopin = 0;
        bulgein = 0;
    else
        if pathSeq1(k) == 0 && pathSeq1(k-1) == 0 && pathSeq1(k-1) < 100
            subseq1coord(k,:) = subseq1coord(k-1,:) + addvec1;
            seqbondvec1(k,:) = subseq1coord(k,:);
            loopin = 0;
            bulgein = 0;
        elseif pathSeq1(k) > 0 && loopin == 0 && pathSeq1(k) < 100;
            loopin = 1;
            perpendvec = cross([addvec1 0],[0 0 1]);
            perpendvec = perpendvec(1:2);
            x1 = subseq1coord(k-1,1);
            y1 = subseq1coord(k-1,2);
            x2 = subseq1coord(k-1,1)+perpendvec(1);
            y2 = subseq1coord(k-1,2)+perpendvec(2);
            r = 0.5/(sind((360/(4+(2*pathSeq1(k))))/2));
            x0 = [subseq1coord(k-1,:)+addvec1];
            centercoord = fsolve(@centercirc,x0,options);
            h1 = centercoord(1);
            k1 = centercoord(2);
            h2 = subseq1coord(k-1,1);
            k2 = subseq1coord(k-1,2);
            subseq1coord(k,:) = fsolve(@coordcirc,x0,options);
            
        elseif pathSeq1(k) > 0 && loopin ~= 0 && pathSeq1(k) < 100;
            
            x0 = [subseq1coord(k-1,:)+addvec1];
            h1 = centercoord(1);
            k1 = centercoord(2);
            h2 = subseq1coord(k-1,1);
            k2 = subseq1coord(k-1,2);
            subseq1coord(k,:) = fsolve(@coordcirc,x0,options);
            
        elseif pathSeq1(k) == 0 && pathSeq1(k-1) > 0 && pathSeq1(k-1) < 100
            loopin = 0;
            bulgein = 0;
            x0 = [subseq1coord(k-1,:)+addvec1];
            h1 = centercoord(1);
            k1 = centercoord(2);
            h2 = subseq1coord(k-1,1);
            k2 = subseq1coord(k-1,2);
            subseq1coord(k,:) = fsolve(@coordcirc,x0,options);
            seqbondvec1(k,:) = subseq1coord(k,:);
            
           
        elseif pathSeq1(k) < 0  && bulgein == 0;
            degdiv = -pathSeq1(k);
            bulgein = 1;
            perpendvec = cross([addvec1 0],[0 0 1]);
            perpendvec = perpendvec(1:2);
            x1 = subseq1coord(k-1,1);
            y1 = subseq1coord(k-1,2);
            x2 = subseq1coord(k-1,1)+perpendvec(1);
            y2 = subseq1coord(k-1,2)+perpendvec(2);
            r = 0.5/(sind((360/(4+(degdiv)))/2));
            x0 = [subseq1coord(k-1,:)+addvec1];
            centercoord = fsolve(@centercirc,x0,options);
            h1 = centercoord(1);
            k1 = centercoord(2);
            h2 = subseq1coord(k-1,1);
            k2 = subseq1coord(k-1,2);
            subseq1coord(k,:) = fsolve(@coordcirc,x0,options);
            
        elseif pathSeq1(k) < 0 && bulgein == 1;
            x0 = [subseq1coord(k-1,:)+addvec1];
            h1 = centercoord(1);
            k1 = centercoord(2);
            h2 = subseq1coord(k-1,1);
            k2 = subseq1coord(k-1,2);
            subseq1coord(k,:) = fsolve(@coordcirc,x0,options);
            
            
        elseif pathSeq1(k) == 0 && pathSeq1(k-1) < 0
            bulgein = 0;
            loopin = 0;
            x0cross = cross([addvec1 0],[0 0 1]);
            subseq1coord(k-1,:);
            x0 = [subseq1coord(k-1,:)+x0cross(1:2)];
            h1 = centercoord(1);
            k1 = centercoord(2);
            h2 = subseq1coord(k-1,1);
            k2 = subseq1coord(k-1,2);
            subseq1coord(k,:) = fsolve(@coordcirc,x0,options);
            seqbondvec1(k,:) = subseq1coord(k,:);
            
            x0 = [subseq1coord(k,:)+x0cross(1:2)];
            h2 = subseq1coord(k,1);
            k2 = subseq1coord(k,2);
            complebondcoord = fsolve(@coordcirc,x0,options);
            perpendvec1 = complebondcoord-subseq1coord(k,:);
            addvec1 = cross([perpendvec1 0],[0 0 -1]);
            addvec1 = addvec1(1:2);
            
        elseif pathSeq1(k) > 100
            subseq1coord(k,:) = [77 77];
            
            
        elseif pathSeq1(k) == 0 && pathSeq1(k-1) > 100
            bulgein = 0;
            loopin = 0;
            degdiv = pathSeq1(k-1)-100;
            perpendvec = cross([addvec1 0],[0 0 1]);
            perpendvec = perpendvec(1:2);
            x1 = subseq1coord(k-(pathSeq1(k-1)-99),1);
            y1 = subseq1coord(k-(pathSeq1(k-1)-99),2);
            x2 = subseq1coord(k-(pathSeq1(k-1)-99),1)+perpendvec(1);
            y2 = subseq1coord(k-(pathSeq1(k-1)-99),2)+perpendvec(2);
            r = 0.5/(sind((360/(4+(degdiv)))/2));
            x0 = [subseq1coord(k-2,:)+addvec1];
            centercoord = fsolve(@centercirc,x0,options);
            
            h1 = centercoord(1);
            k1 = centercoord(2);
            h2 = subseq1coord(k-(pathSeq1(k-1)-99),1);
            k2 = subseq1coord(k-(pathSeq1(k-1)-99),2);
            subseq1coord(k,:) = fsolve(@coordcirc,x0,options);
            seqbondvec1(k,:) = subseq1coord(k,:);
            
            h2 = subseq1coord(k,1);
            k2 = subseq1coord(k,2);
            complebondcoord = fsolve(@coordcirc,x0,options);
            perpendvec1 = complebondcoord-subseq1coord(k,:);
            addvec1 = cross([perpendvec1 0],[0 0 -1]);
            addvec1 = addvec1(1:2);

            
        end
        
    end
end
if pathSeq1(end,:) ~= 0
            subseq1coord(end,:) = [];
            endpt1 = endpt1-1;
end

subseq1coord((subseq1coord(:,1) == 77),:) = [];


% Plot vectors for Seq2
for k = 1:size(pathSeq2)
    if k == 1
        subseq2coord(k,:) = [0.5 0];
        seqbondvec2(k,:) = subseq2coord(k,:);
        addvec2 = [0 1];
        loopin = 0;
        bulgein = 0;
    else
        if pathSeq2(k) == 0 && pathSeq2(k-1) == 0 && pathSeq2(k-1) < 100
            subseq2coord(k,:) = subseq2coord(k-1,:) + addvec2;
            seqbondvec2(k,:) = subseq2coord(k,:);
            loopin = 0;
            bulgein = 0;
        elseif pathSeq2(k) > 0 && loopin == 0 && pathSeq2(k) < 100;
            loopin = 1;
            perpendvec = cross([addvec2 0],[0 0 -1]);
            perpendvec = perpendvec(1:2);
            x1 = subseq2coord(k-1,1);
            y1 = subseq2coord(k-1,2);
            x2 = subseq2coord(k-1,1)+perpendvec(1);
            y2 = subseq2coord(k-1,2)+perpendvec(2);
            r = 0.5/(sind((360/(4+(2*pathSeq2(k))))/2));
            x0 = [subseq2coord(k-1,:)+addvec2];
            centercoord = fsolve(@centercirc,x0,options);
            
            h1 = centercoord(1);
            k1 = centercoord(2);
            h2 = subseq2coord(k-1,1);
            k2 = subseq2coord(k-1,2);
            subseq2coord(k,:) = fsolve(@coordcirc,x0,options);
            
        elseif pathSeq2(k) > 0 && loopin ~= 0 && pathSeq2(k) < 100;
            
            x0 = [subseq2coord(k-1,:)+addvec2];
            h1 = centercoord(1);
            k1 = centercoord(2);
            h2 = subseq2coord(k-1,1);
            k2 = subseq2coord(k-1,2);
            subseq2coord(k,:) = fsolve(@coordcirc,x0,options);
            
        elseif pathSeq2(k) == 0 && pathSeq2(k-1) > 0 && pathSeq2(k-1) < 100
            bulgein = 0;
            loopin = 0;
            x0 = [subseq2coord(k-1,:)+addvec2];
            h1 = centercoord(1);
            k1 = centercoord(2);
            h2 = subseq2coord(k-1,1);
            k2 = subseq2coord(k-1,2);
            subseq2coord(k,:) = fsolve(@coordcirc,x0,options);
            seqbondvec2(k,:) = subseq2coord(k,:);
            
        elseif pathSeq2(k) < 0 && bulgein == 0
            degdiv = -pathSeq2(k);
            bulgein = 1;
            perpendvec = cross([addvec2 0],[0 0 -1]);
            perpendvec = perpendvec(1:2);
            x1 = subseq2coord(k-1,1);
            y1 = subseq2coord(k-1,2);
            x2 = subseq2coord(k-1,1)+perpendvec(1);
            y2 = subseq2coord(k-1,2)+perpendvec(2);
            r = 0.5/(sind((360/(4+(degdiv)))/2));
            x0 = [subseq2coord(k-1,:)+addvec2];
            centercoord = fsolve(@centercirc,x0,options);
            h1 = centercoord(1);
            k1 = centercoord(2);
            h2 = subseq2coord(k-1,1);
            k2 = subseq2coord(k-1,2);
            subseq2coord(k,:) = fsolve(@coordcirc,x0,options);
            
        elseif pathSeq2(k) < 0 && bulgein == 1;
            x0 = [subseq2coord(k-1,:)+addvec2];
            h1 = centercoord(1);
            k1 = centercoord(2);
            h2 = subseq2coord(k-1,1);
            k2 = subseq2coord(k-1,2);
            subseq2coord(k,:) = fsolve(@coordcirc,x0,options);
            
        elseif pathSeq2(k) == 0 && pathSeq2(k-1) < 0
            bulgein = 0;
            loopin = 0;
            x0cross = cross([addvec2 0],[0 0 -1]);
            x0 = [subseq2coord(k-1,:)+x0cross(1:2)];
            h1 = centercoord(1);
            k1 = centercoord(2);
            h2 = subseq2coord(k-1,1);
            k2 = subseq2coord(k-1,2);
            subseq2coord(k,:) = fsolve(@coordcirc,x0,options);
            seqbondvec2(k,:) = subseq2coord(k,:);
            
            x0 = [subseq2coord(k,:)+x0cross(1:2)];
            h2 = subseq2coord(k,1);
            k2 = subseq2coord(k,2);
            complebondcoord = fsolve(@coordcirc,x0,options);
            perpendvec2 = complebondcoord-subseq2coord(k,:);
            addvec2 = cross([perpendvec2 0],[0 0 1]);
            addvec2 = addvec2(1:2);
            
            
        elseif pathSeq2(k) > 100
            subseq2coord(k,:) = [77 77];
            
            
        elseif pathSeq2(k) == 0 && pathSeq2(k-1) > 100
            bulgein = 0;
            loopin = 0;
            degdiv = pathSeq2(k-1)-100;
            perpendvec = cross([addvec2 0],[0 0 -1]);
            perpendvec = perpendvec(1:2);
            x1 = subseq2coord(k-(pathSeq2(k-1)-99),1);
            y1 = subseq2coord(k-(pathSeq2(k-1)-99),2);
            x2 = subseq2coord(k-(pathSeq2(k-1)-99),1)+perpendvec(1);
            y2 = subseq2coord(k-(pathSeq2(k-1)-99),2)+perpendvec(2);
            r = 0.5/(sind((360/(4+(degdiv)))/2));
            x0 = [subseq2coord(k-(pathSeq2(k-1)-99),:)+addvec2];
            centercoord = fsolve(@centercirc,x0,options);
            
            h1 = centercoord(1);
            k1 = centercoord(2);
            h2 = subseq2coord(k-(pathSeq2(k-1)-99),1);
            k2 = subseq2coord(k-(pathSeq2(k-1)-99),2);
            subseq2coord(k,:) = fsolve(@coordcirc,x0,options);
            seqbondvec2(k,:) = subseq2coord(k,:);
            
            x0 = [subseq2coord(k,:)+addvec2];
            h2 = subseq2coord(k,1);
            k2 = subseq2coord(k,2);
            complebondcoord = fsolve(@coordcirc,x0,options);
            perpendvec2 = complebondcoord-subseq2coord(k,:);
            addvec2 = cross([perpendvec2 0],[0 0 1]);
            addvec2 = addvec2(1:2);
            
            
            
        end
    end
    
end

if pathSeq2(end,:) ~= 0
        subseq2coord(end,:) = [];
        endpt2 = endpt2-1;
end

subseq2coord((subseq2coord(:,1) == 77),:) = [];

% Superending coordinates

danglend1 = length(Seq1) - endpt1;
danglend2 = length(Seq2) - endpt2;
r = 0.5/(sind((360/(4+(2*(mean([danglend1 danglend2])))))/2));

if endpt1 < length(Seq1)
    for k = 2:1+danglend1
        if k == 2
            endseq1coord(1,:) = [subseq1coord(end,:)];
            perpendvec = cross([addvec1 0],[0 0 1]);
            perpendvec = perpendvec(1:2);
            x1 = endseq1coord(k-1,1);
            y1 = endseq1coord(k-1,2);
            x2 = endseq1coord(k-1,1)+perpendvec(1);
            y2 = endseq1coord(k-1,2)+perpendvec(2);
            x0 = [endseq1coord(k-1,:)+addvec1];
            centercoord = fsolve(@centercirc,x0,options);
            h1 = centercoord(1);
            k1 = centercoord(2);
            h2 = endseq1coord(k-1,1);
            k2 = endseq1coord(k-1,2);
            endseq1coord(k,:) = fsolve(@coordcirc,x0,options);
        else
            x0 = [endseq1coord(k-1,:)+addvec2];
            h1 = centercoord(1);
            k1 = centercoord(2);
            h2 = endseq1coord(k-1,1);
            k2 = endseq1coord(k-1,2);
            endseq1coord(k,:) = fsolve(@coordcirc,x0,options);
        end

    end
    endseq1coord(1,:) = [];

else endseq1coord = [];    
end

if endpt2 < length(Seq2)
    for k = 2:1+danglend2
        if k == 2
            endseq2coord(1,:) = [subseq2coord(end,:)];
            perpendvec = cross([addvec2 0],[0 0 -1]);
            perpendvec = perpendvec(1:2);
            x1 = endseq2coord(k-1,1);
            y1 = endseq2coord(k-1,2);
            x2 = endseq2coord(k-1,1)+perpendvec(1);
            y2 = endseq2coord(k-1,2)+perpendvec(2);
            x0 = [endseq2coord(k-1,:)+addvec2];
            centercoord = fsolve(@centercirc,x0,options);
            h1 = centercoord(1);
            k1 = centercoord(2);
            h2 = endseq2coord(k-1,1);
            k2 = endseq2coord(k-1,2);
            endseq2coord(k,:) = fsolve(@coordcirc,x0,options);
        else
            x0 = [endseq2coord(k-1,:)+addvec2];
            h1 = centercoord(1);
            k1 = centercoord(2);
            h2 = endseq2coord(k-1,1);
            k2 = endseq2coord(k-1,2);
            endseq2coord(k,:) = fsolve(@coordcirc,x0,options);
        end

    end
    endseq2coord(1,:) = [];

else endseq2coord = [];    
end

Seq1coord = [startseq1coord;subseq1coord;endseq1coord];
Seq2coord = [startseq2coord;subseq2coord;endseq2coord];

Seq1coordx = Seq1coord(2:end-1,1);
Seq1coordy = Seq1coord(2:end-1,2);
Seq2coordx = Seq2coord(2:end-1,1);
Seq2coordy = Seq2coord(2:end-1,2);



seqbondvecx = [seqbondvec1(:,1),seqbondvec2(:,1)];
seqbondvecy = [seqbondvec1(:,2),seqbondvec2(:,2)];

dotsize = max(abs([Seq1coordx(1)-Seq1coordx(end),Seq1coordy(1)-Seq1coordy(end),...
         Seq2coordx(1)-Seq2coordx(end),Seq2coordy(1)-Seq2coordy(end)]));

% plotting H-bonding
for i = 1:size(seqbondvecx,1)
    f = plot(seqbondvecx(i,:),seqbondvecy(i,:),'black','LineWidth',(50/dotsize));
    hold on
end


% plotting phospate backbone 


f = plot(Seq1coordx,Seq1coordy,'Color',[0.3 0.3 0.3],'LineWidth',(50/dotsize));
hold on
f = plot(Seq2coordx,Seq2coordy,'Color',[0.3 0.3 0.3],'LineWidth',(50/dotsize));
hold on

% plotting sequence

f = gscatter([Seq1coordx;Seq2coordx],...
    [Seq1coordy;Seq2coordy],[Seq1(2:end-1)';Seq2(2:end-1)']...
    ,[0.4 0.4 0.4;0.4 0.4 0.4;0.4 0.4 0.4;0.4 0.4 0.4],'.',...
    (950/dotsize)+(150/dotsize));

hold on

f = gscatter([0;0;0;0;0;Seq1coordx;Seq2coordx],...
    [-1;-1;-1;-1;-1;Seq1coordy;Seq2coordy],['A';'T';'G';'C';'-';Seq1(2:end-1)';Seq2(2:end-1)']...
    ,[1 0.502 0;0.0549 0.3020 0.5725;0 0.6 0.298;1 1 0;1 1 1],'.',950/dotsize);

hold on


text(Seq1coordx(1,:)-1,Seq1coordy(1,:)-1,'5''','HorizontalAlignment'...
    ,'center','Fontsize',300/dotsize);
text(Seq1coordx(end,:)+addvec1(1),Seq1coordy(end,:)+addvec1(2),'3'''...
    ,'HorizontalAlignment','center','Fontsize',300/dotsize);
text(Seq2coordx(1)+1,Seq2coordy(1)-1,'3''','HorizontalAlignment'...
    ,'center','Fontsize',300/dotsize);
text(Seq2coordx(end,:)+addvec2(1),Seq2coordy(end,:)+addvec2(2),'5'''...
    ,'HorizontalAlignment','center','Fontsize',300/dotsize);
t = text(0,min([Seq1coordy(1) Seq2coordy(1)]-2),...
    ['Minimum Free Energy: ',num2str(mfe),' kcal/mol'],'HorizontalAlignment'...
    ,'center','Fontsize',400/dotsize);
t.Units = 'normalized';
t.Position = [0.5 -0.1 0];


lgd = legend('show');
lgd.String(end) = [];
lgd.Units = 'normalized';
lgd.FontSize = 300/dotsize;

set(gcf,'Position',[300 60 700 700]);
axis tight
axis equal
axis off

end

function F = centercirc(h)

global x1 y1 x2 y2 r
F(1) = ((x1-h(1))^2) + ((y1-h(2))^2) - (r^2);
F(2) = ((x2-h(1))^2) + ((y2-h(2))^2) - (r^2);

end

function F = coordcirc(x)

global h1 k1 h2 k2 r
F(1) = ((x(1)-h1).^2) + ((x(2)-k1).^2) - (r.^2);
F(2) = (x(1)-h2)^2 + (x(2)-k2)^2 - 1;
end