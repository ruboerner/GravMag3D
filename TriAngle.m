function ang = TriAngle(A,B,C)
% Triangle computes solid angle spanned by a triangle in R3 as seen from from origin
%
% ang = TriAngle(A,B,C)
% 
% Input:
% =====
% A, B, C: (3 x 1) arrays of triangle vertices ordered counterclockwise
% 
% Output:
% ======
%
% ang: Solid angle in rad
%
% Implemented after Osteroom and Stracke (1983) (https://doi.org/10.1109%2FTBME.1983.325207)
%
arguments
    A (3, 1)
    B (3, 1)
    C (3, 1)
end
An = norm(A);
Bn = norm(B);
Cn = norm(C);
crsBC = cross(B,C);
deter = A(1) * crsBC(1) + A(2) * crsBC(2) + A(3) * crsBC(3);
dotab = A(1) * B(1) + A(2) * B(2) + A(3) * B(3);
dotac = A(1) * C(1) + A(2) * C(2) + A(3) * C(3);
dotbc = B(1) * C(1) + B(2) * C(2) + B(3) * C(3);
tan2ang = deter / (An * Bn * Cn + dotab * Cn + dotac * Bn + dotbc * An);
ang = 2 * atan(tan2ang);

end