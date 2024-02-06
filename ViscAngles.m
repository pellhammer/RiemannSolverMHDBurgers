function [Omega1,Omega2] = ViscAngles(kappa,zeta,c)

% Compute the critical slopes of lines of the two saddle omega_{alpha,beta}
% as the zeros of V
% 
% INPUT:
%         kappa,kappa - (scalars) ,viscosity matrix B=[zeta,kappa;kappa,1]
%
%
% AUTHOR:
% Valentin Pellhammer
% Department of Mathematics and Statistics,
% University of Konstanz, 78457 Konstanz
% email adress: valentin.pellhammer@uni-konstanz.de
% homepage: http://www.math.uni-konstanz.de/~pellhammer/



    Vcoeffs = [kappa,zeta,-1/c,-kappa/c];

    Rts = sort(roots(Vcoeffs));

    Omega1 = Rts(1);
    Omega2 = Rts(2);



end