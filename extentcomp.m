function F = extentcomp(KeqP,E,tolerance,n,v,ntot,vtot)



for k = 1
    for i = 1:6
        y(k,i) = (n(i)+(v(i,:)*E(k,:)'))/(ntot+(vtot*E(k,:)'));
    end
    Eqn{1}(k,:) = (y(k,4)*y(k,5))-(KeqP(k,1)*y(k,3)*y(k,6));
    Eqn{2}(k,:) = y(k,6)-(y(k,1)*y(k,5)*KeqP(k,2));
    Eqn{3}(k,:) = (y(k,4)^2)-((y(k,2)^2)*y(k,3)*KeqP(k,3));
    Eqn{4}(k,:) = y(k,6)-(y(k,2)*y(k,3)*KeqP(k,4));
end

for k = 1
    for l = 1:4
        Eh = E;
        Eh(:,l) = E(:,l)+tolerance;
        for i = 1:6
            yh(k,i) = (n(i)+(v(i,:)*Eh(k,:)'))/(ntot+(vtot*Eh(k,:)'));
        end
        Eqnh{1}(l,:) = (yh(k,4)*yh(k,5))-(KeqP(k,1)*yh(k,3)*yh(k,6));
        Eqnh{2}(l,:) = yh(k,6)-(yh(k,1)*yh(k,5)*KeqP(k,2));
        Eqnh{3}(l,:) = (yh(k,4)^2)-((yh(k,2)^2)*yh(k,3)*KeqP(k,3));
        Eqnh{4}(l,:) = yh(k,6)-(yh(k,2)*yh(k,3)*KeqP(k,4));
    end
        yh = 0;
        Eh = E;
 
    for j = 1:4
        for l = 1:4
            J{1,k}(j,l) = (Eqnh{j}(l)-Eqn{j}(k))/tolerance;
        end
    end
    EDelta(k,:) = J{1,k}\[Eqn{1}(k) Eqn{2}(k) Eqn{3}(k) Eqn{4}(k)]';
    Eprim(k,:) = E(k,:)-EDelta(k,:);
    F = Eprim;
end
end