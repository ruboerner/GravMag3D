function [Hx, Hy, Hz, gx, gy, gz] = get_H(Face, cor, Un, M, density)
% get_H computes the magnetic and gravity anomaly 
% 
% get_H computes the magnetic and gravity anomaly at the origin of a 
% right-handed coordinate system. The anomalies are caused by a solid body
% of homogeneous magnetization and/or density.
% The shape of the body is approximated by a sufficiently large number of
% vertices. The convex (or concave) hull is provided by a Delaunay
% triangulation.
%
% [Bx, By, Bz, gx, gy, gz] = get_H(face, vtx, un, m, density)
%
% Input:
% =====
%
% faces: (nf x 3) array of indices into the array of triangles vtx
% vtx: (nc x 3) array of triangular vertex coordinates
% un: (nf x 3) array of unit vectors normal to the triangular faces
% m: (3 x 1) array of induced magnetization of the body enclosed by triangulation, given in A/m
% density: mass density of the body, given in kg/m^3
%
% Output:
% ======
%
% Bx, By, Bz: Cartesian components of the magnetic anomaly B at the origin of
% a right-handed coordinate system, given in nanoTesla (nT)
% gx, gy, gz: Cartesian components of the gravitational attraction g at the
% origin of a right-handed coordinate system, given in SI units of m/s^2
%
%
arguments
    Face (:, 3)
    cor (:, 3)
    Un (:, 3)
    M (3,1) {mustBeNumeric}
    density double {mustBeNumeric}
end

[Hx, Hy, Hz] = deal(0.0, 0.0, 0.0);
[gx, gy, gz] = deal(0.0, 0.0, 0.0);
rhof = density * 6.6732e-11;
Nf = size(Face, 1);
cor = cor';
M = 1e-7 * M;
for f = 1:Nf
    idx = Face(f, :);
    crs = cor(:, idx);
    A = crs(:, 1);
    B = crs(:, 2);
    C = crs(:, 3);
    N = Un(f, :);
    Omega = TriAngle(C, B, A);

    di = dot(A, N);
    if di < 0
        Omega = -sign(di) * Omega;
    end

    PQR = get_PQR(crs);

    l = N(1);
    m = N(2);
    n = N(3);
    p = PQR(1);
    q = PQR(2);
    r = PQR(3);

    hx = l * Omega + n * q - m * r;
    hy = m * Omega + l * r - n * p;
    hz = n * Omega + m * p - l * q;

    Pd = N * M;
    Hx = Hx + Pd * hx;
    Hy = Hy + Pd * hy;
    Hz = Hz + Pd * hz;

    di = rhof * di;
    gx = gx - di * hx;
    gy = gy - di * hy;
    gz = gz - di * hz;
end

end
