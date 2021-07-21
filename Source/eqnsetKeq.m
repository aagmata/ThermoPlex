function F = eqnsetKeq(c0vector,dGvector,spectemp)

n = length(c0vector);
R = 1.9862;
T = spectemp+273.15;
c0vector = c0vector*1e-6;

if isrow(c0vector)
    c0vector = c0vector';
end
Zvec = exp(-1000*dGvector./(R*T));

F = @eqnsys;

% combination with reptition formula: (n+k-1)!/((n-1)!*k!)
kindex = factorial(n+1)/(factorial(n-1)*2);

    function f = eqnsys(N)
        if isrow(N)
            N = N';
        end
        
        for i = 1:n
            
            upmat = zeros(n);

            for k = 1:kindex
                l = n - floor(sqrt(-8*(k-1) + 4*(n+1)*(n)-7)/2.0 - 0.5);
                m = k + l - 1 - (n+1)*(n)/2 + (n-l+2)*((n-l+2)-1)/2;
                upmat(l,m) = k;

            end

            r = unique([upmat(i,:) upmat(:,i)']);
            r = r(r~=0);
            l = n - floor(sqrt(-8*(r-1) + 4*(n+1)*(n)-7)/2.0 - 0.5);
            m = r + l - 1 - (n+1)*(n)/2 + (n-l+2).*((n-l+2)-1)./2;
  
            f(i,1) = c0vector(i)-N(i)-(Zvec(r(i))'.*N(l(i)).*(N(l(i))))...
                -sum(Zvec(r)'.*N(l).*N(m));
            
        end
        
        
    end


end