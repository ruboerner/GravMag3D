function ang = TriAngle(A,B,C)
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