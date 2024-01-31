function [S1,S2] = riemannSolverBurgersSquared(uL,vL,uR,vR,kappa,evalPts)
% FUNCTION
% S = riemannSolverBurgersSquared(uL,uR,vL,vR,kappa,evalPts)
% Background viscosity dependent exact Riemann solver for two inviscid Burgers equations.
%
% DESCRIPTION:
% This function computes the values of the unique solution to the Riemann problem with initial
% data ((uL,uR),(vL,vR)) at the points evalPts in space and the time level
% t=1 such that every discontonuity has a viscous profile w.r.t. the
% viscosity matrix [1,kappa;kappa,1] for -1 < kappa < 1.
% This code refers to the paper "Dependence on the background viscosity of 
% solutions to a prototypical non-strictly hyperbolic system of conservation laws"
% written by Heinrich Freistühler and the author of this code.
%
% INPUT:
%         uL - (array) left values in the u-variable
%         uR - (array) right values in the u-variable
%         vL - (array) left values in the v-variable
%         vR - (array) right values in the v-variable
%         kappa - (double) viscosity parameter. Number between -1 and 1
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
%
% Date: March 2020
% Last revision: March 2020


    if kappa<=-1 || kappa>=1
        error('The third argument must be a value in (-1,1)');       
    end


    % handle 3. in ASSUMPTIONS AND LIMITATIONS
    if length(evalPts) == 1
       evalPts = ones(size(uL)).*evalPts; 
    end
   
    % speeds of possible shockwaves
    uS = (1/2).*(uL+uR);
    vS = (1/2).*(vL+vR);
    
    % handle dependency of kappa
    if kappa > 0
        f = (1+kappa+sqrt((1+kappa).^2-4*kappa.^2))./(2*kappa);  
        Z = @(um,up,v) ((1/2).*(um-up)) <= (f.*abs(v-(1/2).*(um+up)));
    else
        f = Inf;  
        Z = @(um,up,v) ones(size(um)); 
    end
    % scalar rarefaction wave function
    rho = @(um,up,sigma) ((sigma<um) .* um  + (sigma>=um).*(sigma<=up).*sigma  + (sigma>up) .*up);
   
    % use logical matrices to avoid if-statements and handle multible Riemann problems at once
    C1 = (uL<=uR) .* (vL<= vR);
    C2 = (uL<=uR) .* (vL> vR) .* Z(vL,vR,rho(uL,uR,vS));
    C3 = (uL>uR) .* (vL<= vR) .* Z(uL,uR,rho(vL,vR,uS));
    C41 = (uL>uR) .* (vL> vR) .* ((uL+uR) < (vL+vR)) .* Z(uL,uR,vL) .* Z(vL,vR,uR);
    C42 = (uL>uR) .* (vL> vR) .* ((uL+uR) >= (vL+vR)) .* Z(uL,uR,vR) .* Z(vL,vR,uL);

    NCv = (uL>uR) .* (vL < ((1/2)*(uL+uR+(1/f)*(uL-uR)))) .* (vR > ((1/2)*(uL+uR-(1/f)*(uL-uR))));
    NCu = (vL>vR) .* (uL < ((1/2)*(vL+vR+(1/f)*(vL-vR)))) .* (uR > ((1/2)*(vL+vR-(1/f)*(vL-vR))));
    
    
    %%% values for the u-variable (S1) %%%
    
    % setting speed of the undercompresive shock as well as the
    % intermediate states
    ucS = vS;
    um1 = -abs((1/f)*((1/2)*(vL-vR)))+(1/2)*(vL+vR);
    um2 = abs((1/f)*((1/2)*(vL-vR)))+(1/2)*(vL+vR);
    
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
    vm1 = -abs((1/f)*((1/2)*(uL-uR)))+(1/2)*(uL+uR);
    vm2 = abs((1/f)*((1/2)*(uL-uR)))+(1/2)*(uL+uR);
    
    s1 = (1/2).*(vL +vm1);
    s2 = (1/2).*(vm2+vR);
        
    RSR = (vL <= vm1) .* (vR >= vm2);
    SSR = (vL > vm1) .* (vR >= vm2);
    SSS = (vL > vm1) .* (vR < vm2);
    RSS = (vL <= vm1) .* (vR < vm2);
 
    % setting values for the solution in the v-variable for the classical and non-classical solution (S1)
    S2class =  C1 .* ((evalPts <= vL) .* vL   +   (evalPts > vL).*(evalPts <= vR).*evalPts  +   (evalPts > vR).*vR)...
        + C2 .* ((evalPts <= vS).*vL  +  (evalPts > vS).*vR)...
        + C3 .* ((evalPts <= vL) .* vL   +   (evalPts > vL).*(evalPts <= vR).*evalPts  +   (evalPts > vR).*vR)...
        + C41.* ((evalPts <= vS).*vL  +  (evalPts > vS).*vR)...
        + C42.* ((evalPts <= vS).*vL  +  (evalPts > vS).*vR);
    S2nonclass = NCu .* ((evalPts <= vS).*vL  +  (evalPts > vS).*vR)...
                +NCv .* (RSR .* ((evalPts <= vL) .* vL   +   (evalPts > vL).*(evalPts <= vm1).*evalPts...  
                              + (evalPts > vm1).*(evalPts < ucS).*vm1  +  (evalPts >= ucS).*(evalPts <= vm2).*vm2...
                              + (evalPts > vm2).*(evalPts < vR).*evalPts + (evalPts >= vR) .* vR)...
                       + SSR .* ((evalPts <= s1) .* vL   +   (evalPts > s1).*(evalPts <= ucS).*vm1...  
                              + (evalPts > ucS).*(evalPts < vm2).*vm2  +  (evalPts >= vm2).*(evalPts <= vR).*evalPts...
                              + (evalPts > vR).*vR)...
                       + SSS .* ((evalPts <= s1) .* vL   +   (evalPts > s1).*(evalPts <= ucS).*vm1...  
                              + (evalPts > ucS).*(evalPts < s2).*vm2  +  (evalPts >= s2).*vR)...
                       + RSS .* ((evalPts <= vL) .* vL   +   (evalPts > vL).*(evalPts <= vm1).*evalPts...  
                              + (evalPts > vm1).*(evalPts < ucS).*vm1  +  (evalPts >= ucS).*(evalPts <= s2).*vm2...
                              + (evalPts > s2).*vR)); 
    S2 = S2class+S2nonclass;
    S = [S1;S2];
end