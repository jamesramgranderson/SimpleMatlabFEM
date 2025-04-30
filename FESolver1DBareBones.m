% Barebones FESolver1D 

% Written by James Anderson, April 2025
%% Prep Workspace
clear
clf
close all
clc

%% Define The Problem
% Model Definition
tot_len = 10; % Total length of problem
num_el = 10; %# of elements desired
A = [1,5]; % cross-sectional area - add more if discontinuous
A_bounds = [0,5;5,10]; % input boundaries of each cross-sectional area as row in matrix

% Boundary Conditions
pos_fixed = [0]; % Position of fixed nodes - Make sure on real nodes!
pos_loads = []; % Position of loads
load_vals = []; % load values

% Material Conditions
E = [4]; % Young's Modulus - more if discontinuous, order left to right
E_bounds = [0,10]; % input boundaries of each modulus as row in matrix

%% Create Mesh
num_nod = num_el +1; % Number of nodes

% Define the nodes struct
nodes.num = [1:num_nod]; % node numbers (id's)
nodes.pos = [0:tot_len/num_el:tot_len]; % node positions
nodes.fixed = zeros(1,num_nod);
nodes.load = zeros(1,num_nod);

% Assigning Fixed BCs to appropriate nodes
for i = 1:length(pos_fixed)
    fixed = find(nodes.pos == pos_fixed(i));
    nodes.fixed(fixed) = 1;
end

% Assigning loads to appropriate nodes
for i = 1:length(load_vals)
    loaded = find(nodes.pos == pos_loads(i));
    nodes.load(loaded) = load_vals(i);
end

% Defining Elements Struct
elem.num = [1:num_el];
elem.E = zeros(1,num_el);
elem.A = zeros(1,num_el);
elem.k = zeros(1,num_el);
elem.len = zeros(1,num_el);
elem.mid = zeros(1,num_el); % midpoint of each element
elem.n1 = nodes.num(1:num_nod-1);
elem.n2 = nodes.num(2:num_nod);

% Defining element midpoints
for i = 1:num_el
    elem.mid(i) = (nodes.pos(i) + nodes.pos(i+1))/2;
end

% Determining element lengths
for i = 1:num_el
    elem.len(i) = nodes.pos(elem.n2(i))-nodes.pos(elem.n1(i));
end

% Apply E values to according elements
for i = 1:num_el
    for j = 1:length(E)
        if elem.mid(i) >= E_bounds(j,1) && elem.mid(i) <= E_bounds(j,2)
            elem.E(i) = E(j);
        end
    end
end

% Apply A values to appropriate elements
for i = 1:num_el
    for j = 1:length(A)
        if elem.mid(i) >= A_bounds(j,1) && elem.mid(i) <= A_bounds(j,2)
            elem.A(i) = A(j);
        end
    end
end

% Determine k value for each element
for i = 1:num_el
    elem.k(i) = elem.A(i) * elem.E(i) / elem.len(i);
end

%% Create Stiffness Matrix AE/L

% Create k, the global stiffness matrix
k = zeros(num_nod);
for i = 1:num_nod-1
    for j = 1:num_nod-1
        if i == j
            k(i,j) = k(i,j) + elem.k(i);
            k(i+1,j) = -elem.k(i);
            k(i,j+1) = -elem.k(i);
            k(i+1,j+1) = elem.k(i);
        end
    end
end

%% Finalize F and U Matrices, and apply BCs in U, F, & k - Remember F=k*U

% F matrix - this is the one that holds the load values
F = nodes.load';

% U matrix - This is the one with the displacements
U = zeros(1,num_nod);

% Apply fixed boundary conditions by removing rows/columns
free_nod = find(nodes.fixed == 0); % Unconstrained DOFs
fixed_nod = find(nodes.fixed == 1);
k_reduced = k(free_nod, free_nod);
khold_reduced = khold(free_nod, free_nod);
F_reduced = F(free_nod);

% Solve only for free DOFs
U(free_nod) = k_reduced \ F_reduced;


%% Postprocessing

% New node positions
for i = 1:num_nod
    nodes.ppos(i) = nodes.pos(i) + U(i);
end

% New element lengths
for i = 1:num_el
    elem.plen(i) = nodes.ppos(elem.n2(i))-nodes.ppos(elem.n1(i));
    if elem.plen(i) <= 0
        error('Negative jacobian at element %f0.0',i)
    end
end

% Calculating stress -> remember stress = E*epsilon
elem.stress = zeros(1, num_el);  % Preallocate
for i = 1:num_el
    u1 = U(elem.n1(i));
    u2 = U(elem.n2(i));
    elem.stress(i) = elem.E(i) * (u2 - u1) / elem.len(i);
end