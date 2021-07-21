%% Monte Carlo Simulation (3 variables)
n = 1000;
T = randi([45 75],n,1);
MgCl = 1e-3+1e-2*rand(n,1);
alpha = 10.^(10*rand(n,1));
X = zeros(n,1);

parfor i = 1:n
    X(i) = sens_fxn(T(i),MgCl(i),alpha(i));
end

%% Monte Carlo Simulation (2 variables)
T2 = 62;
MgCl2 = 1e-3+1e-2*rand(n,1);
alpha2 = 10.^(10*rand(n,1));
X2 = zeros(n,1);

parfor i = 1:n
    X2(i) = sens_fxn(T2,MgCl2(i),alpha2(i));
end


%% Plotting
c = ones(n,1)+double(alpha2<30);
colorfill = [.3 .3 .3];


subplot(2,3,1); scatter(alpha,X,10,'filled','MarkerEdgeColor',colorfill), set(gca,'xscale','log'),xlabel('\alpha'), ylabel('\chi')
subplot(2,3,2); scatter(MgCl,X,10,'filled','MarkerEdgeColor',colorfill),xlabel('MgCl2'), ylabel('\chi')
subplot(2,3,3); scatter(T,X,10,'filled','MarkerEdgeColor',colorfill),xlabel('T (°C)'), ylabel('\chi')
subplot(2,3,4); scatter(alpha2,X2,10,c,'filled','MarkerEdgeColor',colorfill), set(gca,'xscale','log'),xlabel('\alpha'), ylabel('\chi')
subplot(2,3,5); scatter(MgCl2,X2,10,c,'filled','MarkerEdgeColor',colorfill),xlabel('MgCl2'), ylabel('\chi')
