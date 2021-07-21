function f = eqnsetLUM(c0vector,dGvector,spectemp)


n = length(c0vector);
spectemp = spectemp+273.15;

if isrow(c0vector)
    c0vector = c0vector';
end

f = @eqnsys;

    function F = eqnsys(N)
    
    if isrow(N)
        N = N';
    end    
        
    % combination with reptition formula: (n+k-1)!/((n-1)!*k!)
    kindex = factorial(n+1)/(factorial(n-1)*2);
    
    % material balance eqns
    
    upmat = zeros(n);
        
    for k = 1:kindex
        l = n - floor(sqrt(-8*(k-1) + 4*(n+1)*(n)-7)/2.0 - 0.5);
        m = k + l - 1 - (n+1)*(n)/2 + (n-l+2)*((n-l+2)-1)/2;
        upmat(l,m) = k;

    end
    for i = 1:n
        r = unique([upmat(i,:) upmat(:,i)']);
        r = n+r(r~=0);
        
        F(1,i) = c0vector(i) - (N(i).^2)- N(r(i)).^2- ...
            sum(N(r).^2);

    end

    % single stranded equilibrium eqns
    for i = 1:n
        F(1,n+i) = log((N(i).^2)./c0vector(i))...
            + ((N(n+kindex+i)));

    end


    % double stranded equilibrium eqns
    for k = 1:kindex
        i = n - floor(sqrt(-8*(k-1) + 4*(n+1)*(n)-7)/2.0 - 0.5);
        j = k + i - 1 - (n+1)*(n)/2 + (n-i+2)*((n-i+2)-1)/2;
        
        F(1,(2*n)+k) = ((1000.*dGvector(k))/(1.9872.*spectemp)) ...
            +log((((N(n+k).^2)*1e6)./(c0vector(i).*c0vector(j)))) ...
            +((N(n+kindex+i))+(N(n+kindex+j)));
    end
    end



end
