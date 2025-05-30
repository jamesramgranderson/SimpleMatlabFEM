%% 1D Mechanical FE Solver
% This code uses linear elements and the direct stiffness method to solve a
% simple or complex 1D problem, using linear interpolation to provide
% continuity in solutions

% Written by James Anderson with only a little help from ChatGPT
% May 2025

%% Prep Workspace
clear
clf
close all
clc

%% Define The Problem

% Total length of problem
tot_len = 3; 
% # of elements desired
num_el = 8; 

% % cross-sectional area - add more if discontinuous
% A = [pi*.0225,pi*.01,pi*.0025]; 
% % input boundaries of each cross-sectional area as row in matrix
% A_bounds = [0,1;1,2;2,3]; 


% Area for Tapered Bar Problem
rightr = .05; % Right radius
leftr = .15; % Left radius
dr = (rightr-leftr)/num_el;
el_len = tot_len/num_el;
for i = 1:num_el
    A(i) = pi*(leftr+i*dr)^2;
    A_bounds(i,1) = (i-1)*el_len;
    A_bounds(i,2) = i*el_len;
end


% Boundary Conditions
pos_fixed = [0]; % Position of fixed nodes - Make sure on real nodes!
pos_loads = [3]; % Position of loads
load_vals = [100]; % load values
pos_const_mo = [0]; % Constrained motion positions
const_mo = [0];

% Material Conditions
E = [269770]; % Young's Modulus - more if discontinuous, order left to right
E_bounds = [0,3]; % input boundaries of each modulus as row in matrix

%% Create Mesh
num_nod = num_el +1; % Number of nodes

% Define the nodes struct
nodes.num = [1:num_nod]; % node numbers (id's)
nodes.pos = [0:tot_len/num_el:tot_len]; % node positions
nodes.fixed = zeros(1,num_nod);
nodes.load = zeros(1,num_nod);
nodes.const = zeros(1,num_nod);

% Error codes to verify no human error
for i = 1:length(pos_fixed) 
    if isempty(find(nodes.pos == pos_fixed(i)))
        error('fixed BC doesn''t line up with actual nodes')
    end
end
for i = 1:length(pos_loads)
    if isempty(find(nodes.pos == pos_loads(i)))
        error('load BC doesn''t line up with actual nodes')
    end
end
for i = 1:length(E_bounds) 
    if isempty(find(nodes.pos == E_bounds(i)))
        error('E_bounds doesn''t line up with actual nodes')
    end
end
% for i = 1:length(A_bounds)
%     if isempty(find(nodes.pos == A_bounds(i), 1))
%         error('A_bounds doesn''t line up with actual nodes')
%     end
% end
for i = 1:length(const_mo)
    if isempty(find(nodes.pos == pos_const_mo(i)));
        error('pos_const_mo doesn''t line up with actual nodes')
    end
end
if length(pos_loads) ~= length(load_vals)
    error('pos_loads and load_vals different lengths')
end
if isempty(pos_const_mo) && isempty(pos_fixed)
    error('Must constrain motion at at least one point')
end

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
    
% Assigning Constrained motions to appropriate nodes
for i = 1:length(const_mo)
    constrained = find(nodes.pos == pos_const_mo(i));
    nodes.const(constrained) = const_mo(i);
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
%% Visualize Mesh
figure
hold on

kmin = min(elem.k); % Used for drawing elements (lines) the right color
kmax = max(elem.k);
Emin = min(elem.E); % Used for drawing elements (rectangles) the right color
Emax = max(elem.E);

% Draws elements with color associated with the stiffness of the element
for i = 1:num_nod
    if i <= num_el % for nodes with elements to their right
        if nodes.fixed(i) == 1
            plot(nodes.pos(i),0,'or') % fixed nodes red
        elseif nodes.const(i) ~= 0
            plot(nodes.pos(i),0,'ok') % constrained nodes black
        else
            plot(nodes.pos(i),0,'og') % free nodes green
        end 
%         if kmax - kmin < 1e-3
        if Emax-Emin < 1e-3
            colorval = 0.5;
        else 
%             colorval = (elem.k(i) - min(elem.k))/(max(elem.k)-min(elem.k));
              colorval = (elem.E(i) - min(elem.E))/(max(elem.E)-min(elem.E));
        end
%         line([nodes.pos(i),nodes.pos(i+1)],[0,0],'Color',[.5,colorval,.5])
        rectangle('Position',[nodes.pos(i),-A(i)/2,el_len,A(i)],'FaceColor',[.5,colorval,.5])
    else % for last node only
        if nodes.fixed(i) == 1
            plot(nodes.pos(i),0,'or') % red if fixed
        elseif nodes.const(i) ~= 0
            plot(nodes.pos(i),0,'ok') % black if constrained
        else
            plot(nodes.pos(i),0,'og') % green if free
        end 
    end
end
xlim([-1,tot_len + 1])
ylim([-max(A),max(A)])

% Draws arrow for applied forces
drawArrow = @(x,y) quiver( x(1),y(1),x(2)-x(1),y(2)-y(1),0,'Color','k'); % Defines arrow-drawing function
for i = 1:num_nod
    if nodes.load(i) > 0
        drawArrow([nodes.pos(i)-0.5,nodes.pos(i)+0.5] , [0.1,0.1]);
    elseif nodes.load(i)<0
        drawArrow([nodes.pos(i)+0.5,nodes.pos(i)-0.5] , [0.1,0.1]);
    end
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

% F matrix - this is the one that with the loads
F = nodes.load';

% U matrix - This is the one with the displacements
U = zeros(1,num_nod);

% hold onto original k for later use
khold = k;

% Apply motion constraints by modifying k and F
constrained_dofs = find(nodes.const ~= 0);
for i = constrained_dofs
    k(i,:) = 0;         % Zero out row
%     k(:,i) = 0;         % Zero out column
    k(i,i) = 1;         % Place 1 on the diagonal
    F(i) = nodes.const(i);  % Set F = prescribed displacement
end

% Apply fixed boundary conditions by removing rows/columns
free_nod = find(nodes.fixed == 0); % Unconstrained DOFs
k_reduced = k(free_nod, free_nod);
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

%% Visualize Postprocessing

% Create colorscale for stress viz
maxstress = max(elem.stress);
minstress = min(elem.stress);

% Plots new node positions with lines between them of a color scaled to the
% stress they experience - red = positive, green = negative
figure
hold on
for i = 1:num_nod
    if i <= num_el
        if elem.stress(i) > 0
            plot(nodes.ppos(i),0,'o','Color',[.5,.5,.5])
%             line([nodes.ppos(i),nodes.ppos(i+1)],[0,0],'Color',[abs(elem.stress(i)/maxstress),0,0])
            rectangle('Position',[nodes.ppos(i),-A(i)/2,elem.plen(i),A(i)],'FaceColor',[abs(elem.stress(i)/maxstress),0,0])
        else
            plot(nodes.ppos(i),0,'o','Color',[.5,.5,.5])
%             line([nodes.ppos(i),nodes.ppos(i+1)],[0,0],'Color',[0,abs(elem.stress(i)/minstress),0])
            rectangle('Position',[nodes.ppos(i),-A(i)/2,elem.plen(i),A(i)],'FaceColor',[0,abs(elem.stress(i)/minstress,0)])
        end
    else
        plot(nodes.ppos(i),0,'o')  
    end
end
xlim([-1,tot_len+1])
ylim([-max(A),max(A)])


%% Final Value Query
% This will make a popup menu and then use it and the Command Window to
% evaluate a specified value at a specified location
Query = menu('would you like to query a value?','yes','no');

if Query == 1
    qtype = menu('what would you like to query?','displacement','final position','strain','stress') % Makes popup

    x = input('Input Location:'); % Asks for location
    if x > nodes.pos(num_nod)
        disp('please select a value within the original length')
    else 
        if qtype == 1 % If displacement selected
            for i = 1:num_nod
                if x > nodes.pos(i) && x < nodes.pos(i+1) % Finding the right element in question
                    u1 = u(i); % Displacement of first node of element
                    u2 = u(i+1); % Displacement of 2nd node of element
                    N1 = (nodes.pos(i+1)-x)/el_len; % First Shape function
                    N2 = (x-nodes.pos(i))/el_len; % 2nd Shape function
                    displacementx = N1*u1+N2*u2; % Multiplying to interpolate values. Could also use interp1 or interp2
                end
            end
            fprintf('Displacement at x= %2.1f is %2.3f',x,output)
        elseif qtype == 2 % if final position selected
            for i = 1:num_nod
                if x > nodes.pos(i) && x < nodes.pos(i+1)
                    x1 = nodes.ppos(i);
                    x2 = nodes.ppos(i+1);
                    N1 = (nodes.pos(i+1)-x)/el_len;
                    N2 = (x-nodes.pos(i))/el_len;
                    output = N1*x1+N2*x2;
                end
            end
            fprintf('final position of what was originally at x = %2.1f is %2.3f',x,output)
        elseif qtype ==3 % strain - notice no shape functions because it's an element value, not nodal. Could do smoothing with 
            for i = 1:num_nod
                if x > nodes.pos(i) && x < nodes.pos(i+1)
                    x1 = nodes.ppos(i);
                    x2 = nodes.ppos(i+1);
                    u1 = U(i); % Displacement of first node of element
                    u2 = U(i+1);
                    strainx = (u2-u1)/(x2-x1);
                end
            end
            fprintf('strain at x = %2.1f is %2.3f \n',x,strainx)
            disp('This strain does not use nodal averaging - it''s the same strain in the whole element')
        elseif qtype ==4 % stress - uses nodal averaging - averages element stress values at nodes, then uses shape functions to interpolate between those values - could have done this with strain, too, but wanted to show two options
            if x > nodes.pos(num_nod-1)
                fprintf('stress at x = %2.1f is %2.3f (not averaged because it''s in the last element',x,elem.stress(num_el))
            else
                for i = 1:num_nod
                    if x > nodes.pos(i) && x < nodes.pos(i+1)
                        el1 = i-1;
                        el2 = i;
                        el3 = i+1; % these el_ values needed for nodal averaging
                        ns1 = (elem.stress(el1)+elem.stress(el2))/2;
                        ns2 = (elem.stress(el2)+elem.stress(el3))/2; 
                        N1 = (nodes.pos(i+1)-x)/el_len;
                        N2 = (x-nodes.pos(i))/el_len;
                        output = N1*ns1+N2*ns2;
                    end
                end
                fprintf('stress at x = %2.1f is %2.3f (using nodal averaging)',x,output)
            end
        end
    end
else
    disp('ok, thank you') % if didn't want to query a value
end










