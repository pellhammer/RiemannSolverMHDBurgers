function [S1,S2] = RiemannSolverMHDBurgers(uL,vL,uR,vR,B,c,evalPts)
% FUNCTION
% S = RiemannSolverMHDBurgers(uL,uR,vL,vR,kappa,evalPts)
% Background viscosity dependent exact Riemann solver for the MHD-Burgers model.
%
% DESCRIPTION:
% This function computes the values of the unique solution to the Riemann problem with initial
% data ((uL,uR),(vL,vR)) at the points evalPts in space and the time level
% t=1 such that every discontonuity has a viscous profile w.r.t. the symmetric, positive definite
% viscosity matrix B.
%
% INPUT:
%         uL - (array) left values in the u-variable
%         uR - (array) right values in the u-variable
%         vL - (array) left values in the v-variable
%         vR - (array) right values in the v-variable
%         B - (2x2 array) symmetric positive definite viscosity matrix
%         evalPts - (array) evaluation points
% 
% OUTPUT: 
%         S1 - (array) of the u components of the solution(s) evaluated at evalPts and t=1
%         S2 - (array) of the v components of the solution(s) evaluated at evalPts and t=1
% 
% ASSUMPTIONS AND LIMITATIONS: 
%         Only the following combinations of input data are possible:
%         1. uL,uR,vL,vR are numbers and evalPts is a matrix:
%           The output contains the values of the solution of the single Riemann
%           problem with data ((uL,uR),(vL,vR)).
%         2. uL,uR,vL,vR,evalPts are matrices of the same size:
%           The (i,j)-th entry of the output contains the value of the
%           solution to the Riemann problem with data ((uL(i,j),uR(i,j)),(vL(i,j),vR(i,j)))
%           evaluated at evalPts(i,j)
%         3. uL,uR,vL,vR are matrices of the same size and evalPts is a number:
%           The (i,j)-th entry of the output contains the value of the
%           solution to the Riemann problem with data ((uL(i,j),uR(i,j)),(vL(i,j),vR(i,j)))
%           evaluated at evalPts
%
% AUTHOR:
% Valentin Pellhammer
% Department of Mathematics and Statistics,
% University of Konstanz, 78457 Konstanz
% email adress: valentin.pellhammer@uni-konstanz.de
% homepage: http://www.math.uni-konstanz.de/~pellhammer/




    if ~isequal(size(B),[2,2]) || B(1,1)<=0 || det(B)<=0
        error('Fifth argument must be a 2x2 symmteric, positive definite matrix');       
    end
    if ~isscalar(c) || c<=0
        error('Sixth argument must be positive scalar');       
    end


    % first transform the Matrix B to standart form

    zeta = B(1,1)/B(2,2);
    kappa = B(1,2)/B(2,2);


    % handle 3. in ASSUMPTIONS AND LIMITATIONS
    if length(evalPts) == 1
       evalPts = ones(size(uL)).*evalPts; 
    end
   
    % speeds of possible shockwaves
    uS = (1/2).*(uL+uR);
    vS = (1/2).*c*(vL+vR);
    

    % scalar rarefaction wave functions
    rhoU = @(um,up,sigma) ((sigma<um) .* um  + (sigma>=um).*(sigma<=up).*sigma  + (sigma>up) .*up);
    rhoV = @(vm,vp,sigma) ((sigma<c*vm) .* vm  + (sigma>=c*vm).*(sigma<=c*vp).*sigma./c  + (sigma>c*vp) .*vp);
  
    

    
    if kappa<=0  % Riemann solver in the classical case

        
        % set values for u-component            
        S1 = (uL<=uR) .* rhoU(uL,uR,evalPts)...
            +(uL>uR) .*  ((evalPts<=uS) .* uL + (evalPts > uS).*uR);
        %set values for v-component
        S2 = (vL<=vR) .* rhoV(vL,vR,evalPts)...
             +(vL>vR) .*  ((evalPts<=vS) .* vL + (evalPts > vS).*vR);

    else % Riemann solver in the non-classical case
    
    % define visc angles
    [omega1,omega2] = ViscAngles(kappa,zeta,c);


    Z_u = @(um,up,v) ((1/2).*(um-up)) <= (-(1/omega2)* abs(v-(1/(2*c)).*(um+up)));
    Z_v = @(vm,vp,u) ((1/2).*(vm-vp)) <= (-omega1* abs(u-(1/2)*c.*(vm+vp)));


    % use logical matrices to avoid if-statements and handle multible Riemann problems at once
    C1 = (uL<=uR) .* (vL<=vR);
    C2 = (uL<=uR) .* (vL>vR) .* Z_v(vL,vR,rhoU(uL,uR,vS));
    C3 = (uL>uR) .* (vL<=vR) .* Z_u(uL,uR,rhoV(vL,vR,uS));
    C41 = (uL>uR) .* (vL>vR) .* ((uL+uR) < c*(vL+vR)) .*  Z_u(uL,uR,vL) .* Z_v(vL,vR,uR);
    C42 = (uL>uR) .* (vL>vR) .* ((uL+uR) >= c*(vL+vR)) .* Z_u(uL,uR,vR) .* Z_v(vL,vR,uL);

    NCv = (uL>uR) .* ((-(1/omega2).*(vL-(1/(2*c)).*(uL+uR))) < ((1/2)*(uL-uR))) .* (((1/omega2)*(vR-(1/(2*c))*(uL+uR))) < ((1/2)*(uL-uR)));
    NCu = (vL>vR) .* ((-omega1*(uL-(1/2)*c*(vL+vR))) < ((1/2)*(vL-vR))) .* ((omega1*(uR-(1/2)*c*(vL+vR))) < ((1/2)*(vL-vR)));
    
    
    %%% values for the u-variable (S1) %%%
    
    % setting speed of the undercompresive shock as well as the
    % intermediate states
    ucS = vS;
    um1 =  (1/(2*omega1))*(vL-vR) + (1/2)*c*(vL+vR);
    um2 = -(1/(2*omega1))*(vL-vR) + (1/2)*c*(vL+vR);
    
    % speed of possible intermediate shocks
    s1 = (1/2).*(uL +um1);
    s2 = (1/2).*(um2+uR);
        
    RSR = (uL <= um1) .* (uR >= um2);
    SSR = (uL > um1) .* (uR >= um2);
    SSS = (uL > um1) .* (uR < um2);
    RSS = (uL <= um1) .* (uR < um2);
    
    % setting values for the solution in the u-variable for the classical and non-classical solution (S1)
    S1class = C1 .* ((evalPts <= uL) .* uL   +   (evalPts > uL).*(evalPts <= uR).*evalPts  +   (evalPts > uR).*uR)...
            + C2 .* ((evalPts <= uL) .* uL   +   (evalPts > uL).*(evalPts <= uR).*evalPts  +   (evalPts > uR).*uR)...
            + C3 .* ((evalPts <= uS).*uL  +  (evalPts > uS).*uR)...
            + C41.* ((evalPts <= uS).*uL  +  (evalPts > uS).*uR)...
            + C42.* ((evalPts <= uS).*uL  +  (evalPts > uS).*uR);
    S1nonclass = NCv .* ((evalPts <= uS).*uL  +  (evalPts > uS).*uR)...
                +NCu .* (RSR .* ((evalPts <= uL) .* uL   +   (evalPts > uL).*(evalPts <= um1).*evalPts...  
                              + (evalPts > um1).*(evalPts < ucS).*um1  +  (evalPts >= ucS).*(evalPts <= um2).*um2...
                              + (evalPts > um2).*(evalPts < uR).*evalPts + (evalPts >= uR) .* uR)...
                       + SSR .* ((evalPts <= s1) .* uL   +   (evalPts > s1).*(evalPts <= ucS).*um1...  
                              + (evalPts > ucS).*(evalPts < um2).*um2  +  (evalPts >= um2).*(evalPts <= uR).*evalPts...
                              + (evalPts > uR).*uR)...
                       + SSS .* ((evalPts <= s1) .* uL   +   (evalPts > s1).*(evalPts <= ucS).*um1...  
                              + (evalPts > ucS).*(evalPts < s2).*um2  +  (evalPts >= s2).*uR)...
                       + RSS .* ((evalPts <= uL) .* uL   +   (evalPts > uL).*(evalPts <= um1).*evalPts...  
                              + (evalPts > um1).*(evalPts < ucS).*um1  +  (evalPts >= ucS).*(evalPts <= s2).*um2...
                              + (evalPts > s2).*uR) );
        
    S1 = S1class+S1nonclass;
    
    
    %%% values for the v-variable (S2) %%%
    
    % setting speed of the undercompresive shock as well as the
    % intermediate states
    ucS = uS;
    vm1 =  omega2*(1/2)*(uL-uR) + (1/(2*c))*(uL+uR);
    vm2 = -omega2*(1/2)*(uL-uR) + (1/(2*c))*(uL+uR);
    
    s1 = (1/2).*c*(vL + vm1);
    s2 = (1/2).*c*(vm2 + vR);
        
    RSR = (vL <= vm1) .* (vR >= vm2);
    SSR = (vL > vm1) .* (vR >= vm2);
    SSS = (vL > vm1) .* (vR < vm2);
    RSS = (vL <= vm1) .* (vR < vm2);
 
    % setting values for the solution in the v-variable for the classical and non-classical solution (S1)
    S2class =  C1 .* ((evalPts <= c*vL) .* vL   +   (evalPts > c*vL).*(evalPts <= c*vR).*evalPts./c  +   (evalPts > c*vR).*vR)...
        + C2 .* ((evalPts <= vS).*vL  +  (evalPts > vS).*vR)...
        + C3 .* ((evalPts <= c*vL) .* vL   +   (evalPts > c*vL).*(evalPts <= c*vR).*evalPts./c  +   (evalPts > c*vR).*vR)...
        + C41.* ((evalPts <= vS).*vL  +  (evalPts > vS).*vR)...
        + C42.* ((evalPts <= vS).*vL  +  (evalPts > vS).*vR);
    S2nonclass = NCu .* ((evalPts <= vS).*vL  +  (evalPts > vS).*vR)...
                +NCv .* (RSR .* ((evalPts <= c*vL) .* vL   +   (evalPts > c*vL).*(evalPts <= c*vm1).*evalPts./c...  
                              + (evalPts > c*vm1).*(evalPts < ucS).*vm1  +  (evalPts >= ucS).*(evalPts <= c*vm2).*vm2...
                              + (evalPts > c*vm2).*(evalPts < c*vR).*evalPts./c + (evalPts >= c*vR) .* vR)...
                       + SSR .* ((evalPts <= s1) .* vL   +   (evalPts > s1).*(evalPts <= ucS).*vm1...  
                              + (evalPts > ucS).*(evalPts < (c*vm2)).*vm2  +  (evalPts >= (c*vm2)).*(evalPts <= (c*vR)).*evalPts./c...
                              + (evalPts > c*vR).*vR)...
                       + SSS .* ((evalPts <= s1) .* vL   +   (evalPts > s1).*(evalPts <= ucS).*vm1...  
                              + (evalPts > ucS).*(evalPts < s2).*vm2  +  (evalPts >= s2).*vR)...
                       + RSS .* ((evalPts <= c*vL) .* vL   +   (evalPts > c*vL).*(evalPts <= c*vm1).*evalPts./c...  
                              + (evalPts > c*vm1).*(evalPts < ucS).*vm1  +  (evalPts >= ucS).*(evalPts <= s2).*vm2...
                              + (evalPts > s2).*vR)); 
        S2 = S2class + S2nonclass;

    end
end