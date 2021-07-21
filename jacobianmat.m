function f = jacobianmat(x0vector,h,eqnsys)

% A function that outputs the jacobian matrix of the given input equations,
% the given initial independent variable values and the forward distance h 
% (in percentage). The partial differentials are computed through 
% approximation from the definition f'(x0) = (f(x1)-f(x0))/h

n = length(x0vector);

for i = 1:n                                                % variable index
        
    % input matrix for partial differentials
    inpmat0 = x0vector;
    inpmat1 = inpmat0;
    inpmat2 = inpmat0;
    inpmat1(i) = inpmat1(i)-(inpmat0(i).*h);
    inpmat2(i) = inpmat2(i)+(inpmat0(i).*h);
     
    cdy = eqnsys(inpmat2)-eqnsys(inpmat1);
    
    % Jacobian matrix from dy/dx
    J(:,i) = cdy./(2*abs(inpmat0(i).*h));
       
end




f = J;
end