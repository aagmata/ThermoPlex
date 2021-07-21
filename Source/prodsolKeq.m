function N = prodsolKeq(reactantconc,dGvector,spectemp)

n = length(reactantconc);
k = factorial(n+1)/(factorial(n-1)*2);
R = 1.9872;
T = spectemp+273.15;

for r = 1:k
    i = n - floor(sqrt(-8*(r-1) + 4*(n+1)*(n)-7)/2.0 - 0.5);
    j = r + i - 1 - (n+1)*(n)/2 + (n-i+2)*((n-i+2)-1)/2;
    N(r,1) = exp((-1000*dGvector(r))/(R*T))*reactantconc(i)*reactantconc(j);
    
end






end