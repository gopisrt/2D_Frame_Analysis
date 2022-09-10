clear
clc
close all

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

%% declaring symbolic variable for shape functions

syms x

%% steps for giving input to the main script

%% step 1: enter information of nodal coordinates

% nodal_coordinates : descibes the position of each coordinate (x and z)
% x is the horizontal coordinate and z is the vertical coordinate

nodal_coordinates  = cell2table(cell(1, 3), 'VariableNames', ...
    {'node_number'	'x_coordinate'	'z_coordinate'});

%% step 2: enter information of member connectivity

% member_connectivity : describes the index of two ends of the member
member_connectivity  = cell2table(cell(1, 3), 'VariableNames', ...
    {'member_number'	'node_1'	'node_2'});

%% step 3: enter information of sectional properties

% sectional_properties : types of cross-sections and their breadth and depth
sectional_properties  = cell2table(cell(1, 3), 'VariableNames', ...
    {'section_type'	'breadth'	'depth'});

%% step 4: enter information of material properties

% material_properties : types of materials and their youngs modulii
material_properties  = cell2table(cell(1, 2), 'VariableNames', ...
    {'material_type'	'youngs_modulus'});

%% step 5: assign geometrical and material properties to members

% member_geometry : assigning the geometrical and material properties to each member
member_geometry  = cell2table(cell(1, 3), 'VariableNames', ...
    {'member_number'	'section_type'	'material_type'});

%% step 6: calculating the number of load cases

% number_load_cases : number of load cases
number_load_cases  = 1;

%% step 7: assign the support conditions and the possible suport settlements

% restraint_and_support_settlements : assigning the support conditions and the possible suport settlement associated with a particular load case
restraint_and_support_settlements  = cell2table(cell(1, 8), 'VariableNames', ...
    {'node_number'	'restraint_translation_x'	'restraint_translation_z'	...
    'restraint_rotation_y'	'support_settlement_x'	'support_settlement_z'	...
    'support_settlement_y'	'load_case'});

%% step 8: assign joint loads

% joint_loads : joint loads at the specified joints associated with a particular load case
joint_loads  = cell2table(cell(1, 5), 'VariableNames', ...
    {'node_number'	'force_x'	'force_z'	'moment_y'	'load_case'});

%% step 9: assign span point load

% member_loads_local_point_load : point member loads at the specified members associated with a particular load case
member_loads_local_point_load  = cell2table(cell(1, 8), 'VariableNames', ...
    {'member_number'	'Point_load_axial'	'location_axial'	...
    'Point_load_transverse'	'location_transverse'	'point_moment'	...
    'location_moment'	'load_case'});

%% step 10: assign span distributed load

% member_loads_local_distributed_load : distributed member loads at the specified members associated with a particular load case
member_loads_local_distributed_load  = cell2table(cell(1, 14), 'VariableNames', ...
    {'member_number'	'trapezoidal_load_x_1'	'trapezoidal_load_x_2'	...
    'trapezoidal_load_x_location_from_start_node'	'trapezoidal_load_x_length'	...
    'trapezoidal_load_z_1'	'trapezoidal_load_z_2'	...
    'trapezoidal_load_z_location_from_start_node'	...
    'trapezoidal_load_z_length'	'trapezoidal_moment_y_1'	...
    'trapezoidal_moment_y_2'	'trapezoidal_moment_y_location_from_start_node'	...
    'trapezoidal_moment_y_length'	'load_case'});

%% step 11: assign load combinations

% load_combinations : load combinations and the associated scaling
load_combinations  = cell2table(cell(1, 3), 'VariableNames', ...
    {'combination'	'load_case'	'scale'});

%%

save 'input_2d_frame' 





