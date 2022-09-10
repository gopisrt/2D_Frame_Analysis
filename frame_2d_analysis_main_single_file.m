%% giving input
clear
clc
close all

input_file_name = 'input_2d_frame_sample';

load(input_file_name, 'nodal_coordinates', 'member_connectivity', 'sectional_properties', ...
    'material_properties', 'member_geometry', 'number_load_cases', ...
    'restraint_and_support_settlements', 'joint_loads', 'member_loads_local_point_load', ...
    'member_loads_local_distributed_load', 'load_combinations', 'x')

%% description of variables

% descreption of variables of input
% nodal_coordinates : descibes the position of each coordinate (x and z)
% x is the horizontal coordinate and z is the vertical coordinate
% member_connectivity : describes the index of two ends of the member
% sectional_properties : types of cross-sections and their breadth and depth
% material_properties : types of materials and their youngs modulii
% member_geometry : assigning the geometrical and material properties to each member
% number_load_cases : number of load cases
% restraint_and_support_settlements : assigning the support conditions and the possible suport settlement associated with a particular load case
% joint_loads : joint loads at the specified joints associated with a particular load case
% member_loads_local_point_load : point member loads at the specified members associated with a particular load case
% member_loads_local_distributed_load : distributed member loads at the specified members associated with a particular load case
% load_combinations : load combinations and the associated scaling

ndata = nodal_coordinates{:,:};
mdata = member_connectivity{:,:};

n_nodes = size(ndata, 1);
n_members = size(mdata, 1);

%% connectivity & DOF

node_1 = mdata(:, 2);
node_2 = mdata(:, 3);
dof = [3*node_1-2 3*node_1-1 3*node_1 3*node_2-2 3*node_2-1 3*node_2]';

x1 = ndata(node_1, 2);
y1 = ndata(node_1, 3);
x2 = ndata(node_2, 2);
y2 = ndata(node_2, 3);

elem_length = sqrt((x2-x1).^2 + (y2-y1).^2);

%% section properties

n_sections = size(sectional_properties, 1);
Area = zeros(n_members, 1);
I_moment = zeros(n_members, 1);

for i1  = 1 : n_sections

    Area(member_geometry.section_type == i1) = ...
        sectional_properties.breadth(i1)*sectional_properties.depth(i1);
    I_moment(member_geometry.section_type == i1) = ...
        sectional_properties.breadth(i1)*(sectional_properties.depth(i1)^3)/12;

end

%% geometric properties

n_materials = size(material_properties, 1);
Youngs_mod = zeros(n_members, 1);

for i2  = 1 : n_materials

    Youngs_mod(member_geometry.material_type == i2) = ...
        material_properties.youngs_modulus(i2);

end

%% computing transformation & stiffness matrices

Trans_matrix = zeros(6, 6 ,n_members);
K_e_local= zeros(6, 6, n_members);
K_e_global= zeros(6, 6, n_members);
K_assembled = zeros(3*n_nodes, 3*n_nodes);

for i3 = 1 : n_members

    Trans_matrix(:, :, i3) = transformation_matrix(elem_length(i3), x1(i3), y1(i3), x2(i3), y2(i3));

    K_e_local(:, :, i3) = stiffness_matrix(elem_length(i3), Area(i3), I_moment(i3), Youngs_mod(i3));
    K_e_global(:, :, i3) = Trans_matrix(:,:,i3)'*K_e_local(:, :, i3)*Trans_matrix(:, :, i3);

    K_assembled(dof(:, i3), dof(:, i3)) = K_assembled(dof(:, i3), dof(:, i3)) + K_e_global(:, :, i3);

end

%% joint loads

F_nodal = zeros(3*n_nodes, number_load_cases);

n_joint_loads = size(joint_loads, 1);

for i4 = 1 : n_joint_loads

    lc = joint_loads.load_case(i4);
    jl = joint_loads.node_number(i4);
    F_nodal(3*jl-2 : 3*jl, lc) = ...
        [joint_loads.force_x(i4); joint_loads.force_z(i4); ...
        joint_loads.moment_y(i4)];

end

%% member loads

% point load

Fj_e_local= zeros(6, n_members, number_load_cases);

n_point_loads = size(member_loads_local_point_load, 1);

for i5 = 1 : n_point_loads

    pl = member_loads_local_point_load.member_number(i5);
    lc = member_loads_local_point_load.load_case(i5);

    p_axial = member_loads_local_point_load.Point_load_axial(i5);
    a_axial = member_loads_local_point_load.location_axial(i5);

    p_transverse = member_loads_local_point_load.Point_load_transverse(i5);
    a_transverse = member_loads_local_point_load.location_transverse(i5);

    p_moment = member_loads_local_point_load.point_moment(i5);
    a_moment = member_loads_local_point_load.location_moment(i5);

    Fj_e_local(:, pl, lc) = equivalent_point_member_load(Fj_e_local(:, pl, lc), ...
        elem_length(pl), p_axial, ...
        a_axial, p_transverse, a_transverse, p_moment, a_moment);

end

% distributed load

n_distributed_loads = size(member_loads_local_distributed_load, 1);

for i6 = 1 : n_distributed_loads

    dl = member_loads_local_distributed_load.member_number(i6);
    lc = member_loads_local_distributed_load.load_case(i6);

    trapezoidal_axial_x_1 = ...
        member_loads_local_distributed_load.trapezoidal_load_x_1(i6);
    trapezoidal_axial_x_2 = ...
        member_loads_local_distributed_load.trapezoidal_load_x_2(i6);
    trapezoidal_axial_start = ...
        member_loads_local_distributed_load.trapezoidal_load_x_location_from_start_node(i6);
    trapezoidal_axial_length = ...
        member_loads_local_distributed_load.trapezoidal_load_x_length(i6);

    trapezoidal_transverse_z_1 = ...
        member_loads_local_distributed_load.trapezoidal_load_z_1(i6);
    trapezoidal_transverse_z_2 = ...
        member_loads_local_distributed_load.trapezoidal_load_z_2(i6);
    trapezoidal_transverse_start = ...
        member_loads_local_distributed_load.trapezoidal_load_z_location_from_start_node(i6);
    trapezoidal_transverse_length = ...
        member_loads_local_distributed_load.trapezoidal_load_z_length(i6);

    trapezoidal_moment_y_1 = ...
        member_loads_local_distributed_load.trapezoidal_moment_y_1(i6);
    trapezoidal_moment_y_2 = ...
        member_loads_local_distributed_load.trapezoidal_moment_y_2(i6);
    trapezoidal_moment_start = ...
        member_loads_local_distributed_load.trapezoidal_moment_y_location_from_start_node(i6);
    trapezoidal_moment_length = ...
        member_loads_local_distributed_load.trapezoidal_moment_y_length(i6);

    Fj_e_local(:, dl, lc) = equivalent_distributed_member_load(Fj_e_local(:, dl, lc), ...
        x, elem_length(dl), trapezoidal_axial_x_1, trapezoidal_axial_x_2, ...
        trapezoidal_axial_start, trapezoidal_axial_length, ...
        trapezoidal_transverse_z_1, trapezoidal_transverse_z_2, ...
        trapezoidal_transverse_start, trapezoidal_transverse_length, ...
        trapezoidal_moment_y_1, trapezoidal_moment_y_2, ...
        trapezoidal_moment_start, trapezoidal_moment_length);

end

Fj_e_global= zeros(6, n_members, number_load_cases);
Fj_assembled = zeros(3*n_nodes, number_load_cases);

for i11 = 1 : number_load_cases

    for i7 = 1 : n_members

        Fj_e_global(:, i7, i11) = Trans_matrix(:, :, i7)'*Fj_e_local(:, i7, i11);
        Fj_assembled(dof(:, i7), i11) = Fj_assembled(dof(:, i7), i11) + ...
            Fj_e_global(:, i7, i11);

    end

end

%% considering support settlement & nodal load

restraints = zeros(3*n_nodes, 1);
support_settlement = zeros(3*n_nodes, number_load_cases);
n_supports = size(restraint_and_support_settlements, 1);

for i8 = 1 : n_supports

    ss = restraint_and_support_settlements.node_number(i8);
    lc = restraint_and_support_settlements.load_case(i8);

    restraints(3*ss - 2 : 3*ss) = [restraint_and_support_settlements.restraint_translation_x(i8); ...
        restraint_and_support_settlements.restraint_translation_z(i8); ...
        restraint_and_support_settlements.restraint_rotation_y(i8)];

    support_settlement(3*ss - 2 : 3*ss, lc) = [restraint_and_support_settlements.support_settlement_x(i8); ...
        restraint_and_support_settlements.support_settlement_z(i8); ...
        restraint_and_support_settlements.support_settlement_y(i8)];

end

%% identifying restraining degrees of freedom

total_dof = (1 : 3*n_nodes)';
restrained_dof = find(restraints == 1);
unrestrained_dof = setdiff(total_dof,restrained_dof);

%% getting active & restrained part of stiffness matrix & force vector

Kaa = K_assembled(unrestrained_dof,unrestrained_dof);
Kra = K_assembled(restrained_dof,unrestrained_dof);
Kar = K_assembled(unrestrained_dof,restrained_dof);
Krr = K_assembled(restrained_dof,restrained_dof);

F_active = F_nodal(unrestrained_dof, :);
F_fixed_active = Fj_assembled(unrestrained_dof, :);

D_restrained = support_settlement(restrained_dof, :);

%% solving for displacement

n_unrestrained_dof = length(unrestrained_dof);
D_member_local = zeros(6, n_members, number_load_cases);
F_member = zeros(6, n_members, number_load_cases);
D_active = zeros(n_unrestrained_dof, number_load_cases);
D = zeros(3*n_nodes, number_load_cases);

for i9 = 1 : number_load_cases

    D_active(:, i9) = Kaa\(F_active(:, i9) + F_fixed_active(:, i9) - ...
        Kar*D_restrained(:, i9));

    D(:, i9) = support_settlement(:, i9);
    D(unrestrained_dof, i9) = D_active(:, i9);

    %% obtaining member forces

    for i10 = 1 : n_members

        D_member_local(:, i10, i9) = Trans_matrix(:, :, i10)*D(dof(:, i10), i9);
        F_member(:, i10, i9) = -Fj_e_local(:, i10, i9) + ...
            K_e_local(:, :, i10)*D_member_local(:, i10, i9);

    end

end

n_load_combinations = max(load_combinations.combination);
D_member_local_combinations = zeros(6, n_members, n_load_combinations);
F_member_combinations = zeros(6, n_members, n_load_combinations);

for i11 = 1 : n_load_combinations

    lcns = find(load_combinations.combination == i11);
    lcs = load_combinations.load_case(lcns);
    scl = load_combinations.scale(lcns);

    for i12 = 1 : n_members

        D_member_local_combinations(:, i12, i11) = sum(bsxfun(@times, ...
            reshape(D_member_local(:, i12, lcs), 6, length(lcs)), scl'), 2);
        F_member_combinations(:, i12, i11) = sum(bsxfun(@times, ...
            reshape(F_member(:, i12, lcs), 6, length(lcs)), scl'), 2);
    end

end

save('output_2d_frame')

%% supplementary functions

%% local stiffness matrix

function [K_local] = stiffness_matrix(l, A, I, E)

a = E*A/l;
b = 12*E*I/(l^3);
c = 6*E*I/(l^2);
d = 4*E*I/l;
e = 2*E*I/l;

K_local = [a 0 0 -a 0 0;0 b c 0 -b c;0 c d 0 -c e;-a 0 0 a 0 0;0 -b -c 0 b -c;0 c e 0 -c d];

end

%% transformation matrix

function [T] = transformation_matrix(l, x1, y1, x2, y2)

cost = (x2-x1)/l;
sint = (y2-y1)/l;

k = [cost sint 0;-sint cost 0;0 0 1];
y = zeros(3);

T = [k y;y k];

end

%% converting point span loading to equivalent point member load

function [F_ep] = equivalent_point_member_load(F_ep, ll, p_axial, ...
    a_axial, p_transverse, a_transverse, p_moment, a_moment)

F_axial = p_axial*[1-a_axial/ll; 0; 0; a_axial/ll; 0; 0];

F_transverse = p_transverse*[0;...
    ((1 - a_transverse/ll)^2)*(1+2*a_transverse/ll); ...
    a_transverse*(1 - a_transverse/ll)^2; ...
    0; ...
    1-((1 - a_transverse/ll)^2)*(1+2*a_transverse/ll); ...
    -(ll - a_transverse)*(a_transverse/ll)^2];

F_moment = p_moment*[0; ...
    -6*a_moment*(ll - a_moment)/(ll^3);...
    -(ll - a_moment)*(3*a_moment - ll)/(ll^2);...
    0;...
    6*a_moment*(ll - a_moment)/(ll^3);...
    -a_moment*(2*ll - 3*a_moment)/(ll^2)];

F_ep = F_ep + F_axial + F_transverse + F_moment;

end

%% converting distributed span (moment) loading to equivalent distributed member load

function [fd] = force_distributed_moment(x, qm, le, ad, bd)

[~, phid] = shapefunctions_beam(x, le);

temp = qm*phid';

fd = int(sym(temp),ad,le-bd);

fd = double(fd);

end

%% converting distributed span loading to equivalent distributed member load

function [F_ep] = equivalent_distributed_member_load(F_ep, x, ll, ...
    trapezoidal_axial_x_1, trapezoidal_axial_x_2, trapezoidal_axial_start, ...
    trapezoidal_axial_length, trapezoidal_transverse_z_1, ...
    trapezoidal_transverse_z_2, trapezoidal_transverse_start, ...
    trapezoidal_transverse_length, trapezoidal_moment_y_1, ...
    trapezoidal_moment_y_2, trapezoidal_moment_start, ...
    trapezoidal_moment_length)

if trapezoidal_axial_length > 0
    
    t_a_x = trapezoidal_axial_x_1.*...
        (1-(x-trapezoidal_axial_start)./trapezoidal_axial_length) + ...
        trapezoidal_axial_x_2.*(x-trapezoidal_axial_start)./...
        trapezoidal_axial_length;
    
    Fdx = force_distributed_axial(x, t_a_x, ll, trapezoidal_axial_start, ...
        trapezoidal_axial_start + trapezoidal_axial_length);
    F_axial = [Fdx(1); 0; 0; Fdx(2); 0; 0];
    
else
    
    F_axial = zeros(6, 1);
    
end

if trapezoidal_transverse_length > 0
    
    t_t_z = trapezoidal_transverse_z_1.*...
        (1-(x-trapezoidal_transverse_start)./trapezoidal_transverse_length) + ...
        trapezoidal_transverse_z_2.*(x-trapezoidal_transverse_start)./...
        trapezoidal_transverse_length;
    
    Fdz = force_distributed_transverse(x, t_t_z, ll, trapezoidal_transverse_start, ...
        ll - trapezoidal_transverse_start - trapezoidal_transverse_length);
    F_transverse = [0; Fdz(1); Fdz(2); 0; Fdz(3); Fdz(4)];
    
else
    
    F_transverse = zeros(6, 1);
    
end

if trapezoidal_moment_length > 0
    
    t_m_y = trapezoidal_moment_y_1.*...
        (1-(x-trapezoidal_moment_start)./trapezoidal_moment_length) + ...
        trapezoidal_moment_y_2.*(x-trapezoidal_moment_start)./...
        trapezoidal_moment_length;
    
    Fdm = force_distributed_moment(x, t_m_y, ll, trapezoidal_moment_start, ...
        ll- trapezoidal_moment_start - trapezoidal_moment_length);
    F_moment = [0; Fdm(1); Fdm(2); 0; Fdm(3); Fdm(4)];
    
else
    
    F_moment = zeros(6, 1);
    
end

F_ep = F_ep + F_axial + F_transverse + F_moment;

end

%% converting distributed span (axial) loading to equivalent distributed member load

function fd = force_distributed_axial(x, q, l, a1, a2)

px = int(sym(q),a1,l-a2);
px_m = int(sym(q).*x,a1,l-a2);

if px ~= 0
    
    ax = px_m/px;
    fd = double(px*[1-ax/l ax/l]);
    
else
    
    fd = [0; 0];
    
end
end

%% converting distributed span (transverse) loading to equivalent distributed member load

function [fd] = force_distributed_transverse(x, q, le, ad, bd)

[phi, ~] = shapefunctions_beam(x, le);

temp = q*phi';

fd = int(sym(temp),ad,le-bd);

fd = double(fd);

end

%% shape functions

function [phi, phid] = shapefunctions_beam(x, le)

phi1 = 1 - 3*(x/le)^2 + 2*(x/le)^3;
phi2 = x*(1-(x/le))^2;
phi3 = 3*(x/le)^2 - 2*(x/le)^3;
phi4 = x^2/le*((x/le)-1);

phi = [phi1 phi2 phi3 phi4];

phi1d = diff(phi1);
phi2d = diff(phi2);
phi3d = diff(phi3);
phi4d = diff(phi4);

phid = [phi1d phi2d phi3d phi4d];

end
