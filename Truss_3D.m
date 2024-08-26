clc
clear

load coord.txt;
load connect.txt;
load force.txt;
load displacement.txt;

A = 200 * 1e-6;
E = 70 * 1e9;

num_members = length(connect);
num_nodes = length(coord);
dof = 3 * num_nodes;

% Global stiffness matrix
K_g = zeros(dof, dof);
for e = 1 : num_members
    x_1 = coord(connect(e,1),1);
    y_1 = coord(connect(e,1),2);
    z_1 = coord(connect(e,1),3);
    x_2 = coord(connect(e,2),1);
    y_2 = coord(connect(e,2),2);
    z_2 = coord(connect(e,2),3);
    K_e = elementStiffness(A, E, x_1, y_1, z_1, x_2, y_2, z_2);
    
    for i = 1 : 2
        i_g = connect(e, i);
        for j = 1 : 2
            j_g = connect(e, j);
            K_g(3*i_g-2 : 3*i_g, 3*j_g-2 : 3*j_g) = K_g(3*i_g-2 : 3*i_g, 3*j_g-2 : 3*j_g) + K_e(3*i-2: 3*i, 3*j-2: 3*j);
        end
    end
end

% displacements for all dofs
U = zeros(3*num_nodes, 1);
 
% reaction forces for all dofs
P = zeros(3*num_nodes, 1);

% assigning given loads to P
for i = 1 : length(force(:,1))
    n = force(i, 1);
    f_x = force(i, 2);
    f_y = force(i, 3);
    f_z = force(i, 4);
    P(3*n-2) = f_x;
    P(3*n-1) = f_y;
    P(3*n) = f_z;
end

% known dispacements matrix
U_known = zeros(length(displacement), 1); 

% insert known displacements in U_known
for j = 1 : length(displacement) 
    n = displacement(j, 1);
    a = displacement(j, 2);
    if (a == 1)
        U_known(j) = 3*n-2;
    elseif (a == 2)
        U_known(j) = 3*n-1;
    else
        U_known(j) = 3*n;
    end
end

%Initializing partitioned stiffness matrix
K11 = zeros(dof - length(U_known), dof - length(U_known));
K21 = zeros(length(U_known), dof-length(U_known));

% known reaction forces
P_known = zeros(dof - length(U_known), 1);

% unknown displacement
U_unknown = zeros(dof - length(U_known), 1); 
i = 1;

% Adding unknown displacement dofs in U_unknown
for j = 1 : dof
    flag = 0;
    for k = 1 : length(U_known)
        if(U_known(k) == j)
            flag = 1;
        end
    end

    if(flag == 0)
        U_unknown(i) = j;
        i = i+1;
    end
end

% partitioning the global stiffness matrix
for j = 1 : length(U_unknown)
    for i = 1 : length(U_unknown)
        K11(i,j) = K_g(U_unknown(i),U_unknown(j));
    end
    for i = 1 : length(U_known)
        K21(i,j) = K_g(U_known(i),U_unknown(j));
    end

    P_known(j) = P(U_unknown(j));
end

% Solving for unknown displacements
U_u = linsolve(K11, P_known);

% Solving for unknown recations
P_u = K21 * U_u;

% Assembling displacements & reactions
for i = 1 : length(U_unknown)
    U(U_unknown(i)) = U_u(i); 
end

for i = 1 : length(U_known)
    P(U_known(i)) = P_u(i);
end

% Calculating member forces
F = zeros(num_members, 1);

for e = 1 : num_members
    i = connect(e, 1);
    j = connect(e, 2);
    x_1 = coord(j,1);
    y_1 = coord(j,2);
    z_1 = coord(j,3);
    x_2 = coord(i,1);
    y_2 = coord(i,2);
    z_2 = coord(i,3);
    L = sqrt((x_1-x_2)^2 + (y_1-y_2)^2 + (z_1-z_2)^2);
    lx = (x_2 - x_1) / L;
    ly = (y_2 - y_1) / L;
    lz = (z_2 - z_1) / L;
    K_m = [3*i-2, 3*i-1, 3*i, 3*j-2, 3*j-1, 3*j];
    U_m = [U(K_m(1)), U(K_m(2)), U(K_m(3)), U(K_m(4)), U(K_m(5)), U(K_m(6))];
    F(e) = -(A*E/L) * [lx, ly, lz, -lx, -ly, -lz] * U_m';
end

disp("Reaction Forces:");
disp(P);
disp("Displacements:");
disp(U);
disp("Member Forces:");
disp(F);

% Plotting the deformed shape
P_new = zeros(num_nodes, 3);
for i = 1 : num_nodes
    P_new(i,:) = coord(i,:) + 1e+3*[U(3*i-2) U(3*i-1) U(3*i)];
end

for e = 1 : length(connect)
    node1 = connect(e,1);
    node2 = connect(e,2);
    x1 = P_new(node1, 1);
    y1 = P_new(node1, 2);
    z1 = P_new(node1, 3);
    x2 = P_new(node2, 1);
    y2 = P_new(node2, 2);
    z2 = P_new(node2, 3);

    x = [coord(node1, 1) coord(node2, 1)];
    y = [coord(node1, 2) coord(node2, 2)];
    z = [coord(node1, 3) coord(node2, 3)];

    
    L = sqrt((x1-x2)^2 + (y1-y2)^2 + (z1-z2)^2);
    plot3(x, y, z, 'b');         % plotting original shape
    X = (x2 - x1) / 3;
    Y = (y2 - y1) / 3;
    Z = (z2 - z1) / 3;
    
    if F(e) < 0     %compression
        quiver3(x1, y1, z1, X, Y, Z, 'm','linewidth', 1.5); % plotting arrows indicating memberforce
        hold on;
        quiver3(x2, y2, z2, -X, -Y, -Z, 'm', 'linewidth', 1.5);
    else
        quiver3((x1+x2)/2, (y1+y2)/2, (z1+z2)/2, -X, -Y, -Z, 'm', 'linewidth', 1.5);
        hold on;
        quiver3((x1+x2)/2, (y1+y2)/2, (z1+z2)/2, X, Y, Z, 'm', 'linewidth', 1.5);
    end

    ux = [x1 x2];
    uy = [y1 y2];
    uz = [z1 z2];

    plot3(ux, uy, uz, 'r'); % plotting deformed shape
    hold on;

end

title("Truss");