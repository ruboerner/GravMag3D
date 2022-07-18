function PQR = get_PQR(crs)
% get_PQR computes the contributons of each edge of a triangle formed by vertices
% P, Q, and R
%
% PQR = get_PQR(VTX)
% 
% Input:
%
% VTX: (3 x 3) array of the vertex coordinates of a triangle oriented ccw with
% respect to the outer normal direction.
%
% Output:
% 
% The contribution of all three edges to the surface integral formed by the
% area of the triangle PQR.
%
PQR = [0; 0; 0];
edge_index = [1; 2; 3; 1];
for t = 1:3
    p1 = crs(:, edge_index(t));
    p2 = crs(:, edge_index(t + 1));
    v = p2 - p1;
    L = norm(v);
    b = 2 * dot(p1, v);
    r1 = norm(p1);
    denom = r1 + b / 2 / L;
    if abs(denom) < eps()
        I = 1 / L * log(abs(L - r1) / r1);
    else
        I = 1 / L * log( ...
            (sqrt(L^2 + b + r1^2) + L + b / 2 / L) ...
            / denom);
    end
    PQR = PQR + I * v;
end
end