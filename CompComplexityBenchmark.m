%% Randomize Samples
n = 15:100;

for i = 1:5
    parfor j = 1:length(n)
        seq1 = randseq(n(j));
        seq2 = randseq(n(j));
        
        tic;
        ThermoDhyb(seq1,seq2,37,2e-3);
        tEnd = toc;
        
        time(i,j) = tEnd;
    end
end

%% Curve Fitting
t = mean(time,1);
x1 = @(c) c(1)*n + c(2) - t;
x2 = @(c) c(1)*n.^2 + c(2) - t;
x3 = @(c) c(1)*n.^3 + c(2) - t;
x4 = @(c) c(1)*n.^4 + c(2) - t;

c1 = lsqnonlin(x1,[0 0]);
c2 = lsqnonlin(x2,[0 0]);
c3 = lsqnonlin(x3,[0 0]);
c4 = lsqnonlin(x4,[0 0]);

y1 = @(n) c1(1)*n + c1(2);
y2 = @(n) c2(1)*n.^2 + c2(2);
y3 = @(n) c3(1)*n.^3 + c3(2);
y4 = @(n) c4(1)*n.^4 + c4(2);

%% Plotting
%%
hold on
plot(n,y1(n),'LineWidth',1,'color',[1 0.502 0])
boxplot(time,'PlotStyle','Compact','Colors',[0.0549 0.3020 0.5725],'positions',n,...
    'Labels',n,'LabelOrientation','horizontal')
set(gca,'FontSize',15)
xlim([15 100])
xticks(10:10:100)
xticklabels({'10','20','30','40','50','60','70','80','90','100'})
title('ax+b')

%%
hold on
plot(n,y2(n),'LineWidth',1,'color',[1 0.502 0])
boxplot(time,'PlotStyle','Compact','Colors',[0.0549 0.3020 0.5725],'positions',n,...
    'Labels',n,'LabelOrientation','horizontal')
set(gca,'FontSize',15)
xlim([15 100])
xticks(10:10:100)
xticklabels({'10','20','30','40','50','60','70','80','90','100'})
title('ax^{2}+b')
%%
hold on
plot(n,y3(n),'LineWidth',1,'color',[1 0.502 0])
boxplot(time,'PlotStyle','Compact','Colors',[0.0549 0.3020 0.5725],'positions',n,...
    'Labels',n,'LabelOrientation','horizontal')
set(gca,'FontSize',15)
xlim([15 100])
xticks(10:10:100)
xticklabels({'10','20','30','40','50','60','70','80','90','100'})
title('ax^{3}+b')
%%
hold on
plot(n,y4(n),'LineWidth',1,'color',[1 0.502 0])
boxplot(time,'PlotStyle','Compact','Colors',[0.0549 0.3020 0.5725],'positions',n,...
    'Labels',n,'LabelOrientation','horizontal')
set(gca,'FontSize',15)
xlim([15 100])
xticks(10:10:100)
xticklabels({'10','20','30','40','50','60','70','80','90','100'})
title('ax^{4}+b')

%% Saving
h = gcf;
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(h,'filename','-dpdf','-r0')