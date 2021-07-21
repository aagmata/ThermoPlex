function f = NewtonR(eqnsys,initguess,J)

% Newton-Raphson method for solving systems of non-linear equations
% 
% =========================================================================
% Input arguments:
% eqnarray = array of equations in anonymous function form
% initguess = vector of initial guess (n = number of equations)
%
% Output arguments:
% f = solution of the system of equations
% 
% =========================================================================

warning('off', 'MATLAB:nearlySingularMatrix');

if isrow(initguess)
    vector0 = initguess';
else vector0 = initguess;
end

stopper = 0;
tolerance0 = 1e10;
iter = 0;
while stopper == 0
    
    iter = iter + 1;
    
    if J == 1
        [fxnvec,anaJac] = eqnsys(vector0);
        Jacobian = anaJac;                                                 % analytical Jacobian
        
    else
        fxnvec = eqnsys(vector0);
        Jacobian = jacobianmat(vector0,0.0001,eqnsys);                     % finite difference Jacobian
    end
    
    if isrow(fxnvec)
        fxnvec = fxnvec';
    end
    
    delx = Jacobian\fxnvec;
    vectorn = vector0 - delx;
    fval = max(abs(eqnsys(vectorn)));
    
    
    if fval <= 1e-10
        stopper = stopper +1;
    elseif iter == 100
        stopper = stopper +1;
        fprintf('\n')
        fprintf('\n')
        fprintf('Solver stopped: reached maximum iteration = 200')
        fprintf('\n')
        fprintf('\n')
    end
        
    tolerance0 = fval;
    vector0 = abs(vectorn);
    
end

f = vectorn;

end