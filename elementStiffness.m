function [K_e] = ElementStiffness(A, E, x1, y1, z1, x2, y2, z2)

    % Calculation of Direction Cosines
    L = sqrt( (x1-x2)^2 + (y1-y2)^2 + (z1-z2)^2 );
    l = (x1 - x2) / L;
    m = (y1 - y2) / L;
    n = (z1 - z2) / L;
    
    T = [l, m, n, 0, 0, 0; 0, 0, 0, l, m, n];

    K_m = (A*E/L) * [1, -1; -1, 1];
    
    K_e = T'* K_m * T;
end