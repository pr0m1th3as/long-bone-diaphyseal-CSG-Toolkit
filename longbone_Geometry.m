% Copyright (C) 2018-2020 Andreas Bertsatos <abertsatos@biol.uoa.gr>
%
% This program is free software; you can redistribute it and/or modify it under
% the terms of the GNU General Public License as published by the Free Software
% Foundation; either version 3 of the License, or (at your option) any later
% version.
%
% This program is distributed in the hope that it will be useful, but WITHOUT
% ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
% FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
% details.
%
% You should have received a copy of the GNU General Public License along with
% this program; if not, see <http://www.gnu.org/licenses/>.
%
%
function [CS_Geometry, SMoA, polyline] = longbone_Geometry(varargin)
  % function [CS_Geometry, SMoA, polyline] = longbone_Geometry(filename,type,points)
  %
  % This function slices the provided mesh of a humerus, femur or tibia bone
  % at 20%, 35%, 50%, 65% and 80% along the bone's maximum length and returns
  % the geometric properties for each planar cross section.
  %
  % This function loads a triangular mesh stored in Wavefront OBJ file along
  % with the respective Meshlab Point file that contains the two points that
  % define the mediolateral axis of the bone's anatomical orientation. Both
  % files must be present in the working directory and share a common name.
  % Alternatively, the user may provide the 2 points for the mediolateral axis
  % as a third optional argument, if a Meshlab Point file is not present.
  % 
  % example:    points = [x1, y2, z3; x2, y2, z2];
  %
  % The function requires at least 2 input arguments, the first be the filename
  % of the $$.obj including the file extension, i.e. 'ID_humerus.obj' and the
  % second argument a char string defining which bone is analyzed, and may be
  % any of three available options: "Humerus", "Femur" or "Tibia".
  %
  % The mediolateral axis alignment points for the humerus need to be positioned
  % anteriorly on the trochlea and capitulum at the distal end of the humerus.
  % For the femur, the points need to be located posteriorly on the medial and
  % lateral condyles at the distal end of the femur, whereas for the tibia, the
  % points should be positioned posteriorly of the medial and lateral condyles
  % at the proximal end of the tibia. The function will automatically optimize
  % the local extremal points of the mediolateral axis based on the initial user
  % defined points and will use the optimized coronal plane orientation for
  % the calculations of the second moments of area.
  %
  % When providing orientation points through the respective Meshlab point file,
  % caution should be taken to that 'longbone_Geometry.m' function will only
  % use the first two points available in the file. Furthermore, the optimized
  % points are appended in the original .pp file after the user defined points.
  % So, if more points are present in the original .pp file, these are disregarded
  % and eventually deleted from the final .pp file after processing.
  %
  % Independently of the number of input arguments, 'longbone_Geometry.m'
  % returns 3 structure arrays as output arguments with the following fields.
  %
  %
  %
  %   'CS_Geometry(x)'      is an array structure, where x=[1:5] defines the
  %                         order of section from proximal 20% to distal 80%
  %                         of the bone's max length and contains the fields:
  %         'CS_Geometry(x).Area'     : calculated in mm2
  %         'CS_Geometry(x).Perimeter': calculated in mm
  %         'CS_Geometry(x).Centroid' : containing the (x,y,z) coordinates in R3
  %         'CS_Geometry(x).Slice_n'  : containing the (x,y,z) normal vector of
  %                                     the slicing plane
  %         'CS_Geometry(x).Coronal_n': containing the (x,y,z) normal vector of
  %                                     the coronal plane facing anteriorly
  %
  %   'SMoA(x)'             is an array structure containing the fields:
  %         'SMoA(x).Ix'      : 2nd moment of area with respect to x axis,
  %                             which is collinear with the intersection of
  %                             the coronal plane and the sectioning plane
  %         'SMoA(x).Iy'      : 2nd moment of area with respect to y axis,
  %                             which is collinear with the intersection of
  %                             the sagital plane and the sectioning plane
  %         'SMoA(x).Ixy'     : product 2nd moment of area calculated in mm4
  %         'SMoA(x).Imin'    : minimum 2nd moment of area calculated in mm4
  %         'SMoA(x).Imax'    : maximum 2nd moment of area calculated in mm4
  %         'SMoA(x).theta'   : angle of rotation of the principal axis of
  %                             2nd moment of area with respect to x axis
  %                             expressed in degrees
  %
  %   'polyline(x)'         is a structure array containing the fields:
  %         'polyline.poly2D' : Nx2 matrix containing the (x,y) coordinates of
  %                             the cross section on the 2D local axes of the
  %                             sectioning plane ordered counter-clockwise
  %         'polyline.poly3D' : Nx3 matrix containing the original 3D coordinates
  %                             of the cross section ordered counter clockwise
  %
  % The function requires the 'io', 'matgeom' & 'statistics' packages to be loaded.
  % It also relies on the functions 'longbone_maxDistance', 'simple_polygon3D', 
  % 'slice_Mesh_Plane', 'read_MeshlabPoints', 'write_MeshlabPoints,' and 'readObj',
  % which must be present in the working directory.
  
  % check for number of input variables
  if nargin < 2 || nargin > 3
    printf("invalid number of input arguments\n");
    return;
  endif
  % check first two arguments are strings
  if !ischar(varargin{1}(:)') || !ischar(varargin{2}(:)')
    printf("first two input arguments should be char strings\n");
    return;
  else
    filename = varargin{1}(:)';
    type = varargin{2}(:)';
  endif
  % check second argument properly defines a bone
  if !(strcmp(type, "Humerus") || strcmp(type, "Femur") || strcmp(type, "Tibia"))
    printf("bone should be defined either as humerus, femur of tibia\n");
    return;
  endif
	
  % load vertices and faces of triangular mesh from obj file
  [v,f] = readObj(filename);
  % find the maximum distance of the bone
  [maxDistance, maxd_V1, maxd_V2] = longbone_maxDistance(v);
  % calculate the normal vector and the 5 points in R3 that define the 5 slicing
  % planes at 20, 35, 50, 65 and 80% of the maximum length of the bone
  normal = (maxd_V2 - maxd_V1) ./ sqrt(sum((maxd_V2 - maxd_V1).^2));
  point_1 = maxd_V1 + (maxd_V2 - maxd_V1) .* 0.20;
  point_2 = maxd_V1 + (maxd_V2 - maxd_V1) .* 0.35;
  point_3 = maxd_V1 + (maxd_V2 - maxd_V1) .* 0.50;
  point_4 = maxd_V1 + (maxd_V2 - maxd_V1) .* 0.65;
  point_5 = maxd_V1 + (maxd_V2 - maxd_V1) .* 0.80;

  % calculate the sectioning points for each plane
  plane_1 = slice_Mesh_Plane(v,f,point_1,normal);
  plane_2 = slice_Mesh_Plane(v,f,point_2,normal);
  plane_3 = slice_Mesh_Plane(v,f,point_3,normal);
  plane_4 = slice_Mesh_Plane(v,f,point_4,normal);
  plane_5 = slice_Mesh_Plane(v,f,point_5,normal);
	
  % for each cross section calculate the 3D coordinates of its centroid and area
  CS_Geometry(1) = simple_polygon3D(plane_1, normal);
  CS_Geometry(2) = simple_polygon3D(plane_2, normal);
  CS_Geometry(3) = simple_polygon3D(plane_3, normal);
  CS_Geometry(4) = simple_polygon3D(plane_4, normal);
  CS_Geometry(5) = simple_polygon3D(plane_5, normal);
	
  % store the centroids of the initial cross-sectional areas
  Centroid_1 = CS_Geometry(1).Centroid;
  Centroid_2 = CS_Geometry(2).Centroid;
  Centroid_3 = CS_Geometry(3).Centroid;
  Centroid_4 = CS_Geometry(4).Centroid;
  Centroid_5 = CS_Geometry(5).Centroid;

  % check if orientation points are provided as a third input argument
  % otherwise reuse the mesh filename by appending the required extension
  % for Meshlab point files
  if nargin != 3
    filenamePP = filename([1:length(filename)-4]);
    extension = ".pp";
    filenamePP = strcat(filenamePP, extension);
    % load Meshlab points for mediolateral axis from pp file
    MLA_points = read_MeshlabPoints(filenamePP);
    MLA_points(:,1) = [];
  elseif nargin == 3
    MLA_points = varargin{3};
  endif
  % calculate the proximal distal axis of the bone and use it together with
  % the mediolateral axis vector to define the coronal plane and calculate
  % its normal. Additionally, check if the centroids' locations progress from
  % promimal (20%) to distal (80%) and if required reverse their order to
  % comply with this standard.
  proximal = distancePoints(Centroid_1, MLA_points(1,:));
  distal = distancePoints(Centroid_5, MLA_points(1,:));
  if (strcmp(type, "Humerus") || strcmp(type, "Femur")) && (proximal < distal)
    Centroid_1 = CS_Geometry(5).Centroid;
    Centroid_2 = CS_Geometry(4).Centroid;
    Centroid_3 = CS_Geometry(3).Centroid;
    Centroid_4 = CS_Geometry(2).Centroid;
    Centroid_5 = CS_Geometry(1).Centroid;
  elseif (strcmp(type, "Tibia")) && (proximal > distal)
    Centroid_1 = CS_Geometry(5).Centroid;
    Centroid_2 = CS_Geometry(4).Centroid;
    Centroid_3 = CS_Geometry(3).Centroid;
    Centroid_4 = CS_Geometry(2).Centroid;
    Centroid_5 = CS_Geometry(1).Centroid;
  endif

  % calculate the mediolateral axis vector
  MLA_vector = MLA_points(1,:) - MLA_points(2,:);
  % calculate proximal to distal vector
  PDA_vector = Centroid_1 - Centroid_5;
  % calculate the unit normal of the transverse plane facing upwards
  TransPlane_normal = PDA_vector ./ sqrt(sum(PDA_vector.^2));
  % calculate the normal of the coronal plane
  CorPlane_normal = cross(PDA_vector, MLA_vector);
  % normalize the normal vector of the coronal plane
  CorPlane_normal = CorPlane_normal ./ sqrt(sum(CorPlane_normal.^2));
  
  
  % find if coronal plane normal points in the right direction (towards the
  % front) by checking the dot product of the normal with the vector between 
  % the first point of MLA vector and the nearest centroid. If not, reverse
  % the normal so that it points to the front
  %
  % for humerus bone
  if strcmp(type, "Humerus")
    MLP1_C5 = Centroid_5 - MLA_points(1,:);
    if sum(MLP1_C5 .* CorPlane_normal) > 0      % MLA_points should be in front
      CorPlane_normal = CorPlane_normal * -1;
    endif
    % find the points that lie on or below the transverse plane of the distal centroid
    TP_dotP = sum(TransPlane_normal .* (v - Centroid_5), 2);
    v_below = v(find(TP_dotP <= 0),:);
    % find the midpoint between the user defined points and project it to the
    % vertical vector defined by centroids 1 and 5 to find the closest point along
    % the vertical axis
    MLA_midpoint = (MLA_points(1,:) + MLA_points(2,:)) ./ 2;
    C1_midpoint = MLA_midpoint - Centroid_1;
    C1_C5 = Centroid_5 - Centroid_1;
    midpoint_proj = Centroid_1 + (dot(C1_midpoint,C1_C5)/dot(C1_C5,C1_C5)) * C1_C5;
    % recalculate the coronal plane unit vector by normalizing the vector from
    % midpoint_proj to MLA_midpoint
    CorPlane_normal = MLA_midpoint - midpoint_proj;
    CorPlane_normal = CorPlane_normal ./ sqrt(sum(CorPlane_normal.^2));
    % calculate the unit normal of the sagital plane facing to the right side
    SagPlane_normal = cross(CorPlane_normal, TransPlane_normal);
    % find the points that lie in front of the coronal plane of the midpoint_proj point
    CP_dotP = sum(CorPlane_normal .* (v_below - midpoint_proj), 2);
    v_below_front = v_below(find(CP_dotP > 0),:);
    % find the points that lie on one or the other side
    % of the sagital plane of the MLA_midpoint
    ML_dotP = sum(SagPlane_normal .* (v_below_front - MLA_midpoint), 2);
    v_side1 = v_below_front(find(ML_dotP <= 0),:);
    v_side2 = v_below_front(find(ML_dotP >= 0),:);
    % initialize 5 iterations for finding the optimal points
    for i=1:5
      % calculate the distance for each point of each side from the coronal plane
      % passing through the midpoint_proj point
      d_side1 = sum(CorPlane_normal .* (v_side1 - midpoint_proj),2);
      d_side2 = sum(CorPlane_normal .* (v_side2 - midpoint_proj),2);
      % find new points for optimized mediolateral axis vector that have the
      % maximum distance from the initial coronal plane (optimize local extrema)
      [d1,i1]=max(d_side1);
      MLA_optimal_point1 = v_side1(i1,:);
      [d2,i2]=max(d_side2);
      MLA_optimal_point2 = v_side2(i2,:);
      MLA_vector = MLA_optimal_point1 - MLA_optimal_point2;
      % calculate the optimized normal of the coronal plane
      CorPlane_normal = cross(PDA_vector, MLA_vector);
      % normalize the optimized normal vector of the coronal plane
      CorPlane_normal = CorPlane_normal ./ sqrt(sum(CorPlane_normal.^2));
      % check the optimized CorPlane_normal points towards the right direction
      MLP1_C5 = Centroid_5 - MLA_optimal_point1;
      if sum(MLP1_C5 .* CorPlane_normal) > 0      % MLA_points should be in front
        CorPlane_normal = CorPlane_normal * -1;
      endif
    endfor
  % for femur bone
  elseif strcmp(type, "Femur")
    MLP1_C5 = Centroid_5 - MLA_points(1,:);
    if sum(MLP1_C5 .* CorPlane_normal) < 0      % MLA_points should be behind
      CorPlane_normal = CorPlane_normal * -1;
    endif
    % find the points that lie on or below the transverse plane of the distal centroid
    TP_dotP = sum(TransPlane_normal .* (v - Centroid_5), 2);
    v_below = v(find(TP_dotP <= 0),:);
    % find the midpoint between the user defined points and project it to the
    % vertical vector defined by centroids 1 and 5 to find the closest point along
    % the vertical axis
    MLA_midpoint = (MLA_points(1,:) + MLA_points(2,:)) ./ 2;
    C1_midpoint = MLA_midpoint - Centroid_1;
    C1_C5 = Centroid_5 - Centroid_1;
    midpoint_proj = Centroid_1 + (dot(C1_midpoint,C1_C5)/dot(C1_C5,C1_C5)) * C1_C5;
    % recalculate the coronal plane unit vector by normalizing the vector from
    % MLA_midpoint to midpoint_proj
    CorPlane_normal = midpoint_proj - MLA_midpoint;
    CorPlane_normal = CorPlane_normal ./ sqrt(sum(CorPlane_normal.^2));
    % calculate the unit normal of the sagital plane facing to the right side
    SagPlane_normal = cross(CorPlane_normal, TransPlane_normal);
    % find the points that lie behind the coronal plane of the midpoint_proj point
    CP_dotP = sum(CorPlane_normal .* (v_below - midpoint_proj), 2);
    v_below_behind = v_below(find(CP_dotP < 0),:);
    % find the points that lie on one or the other side
    % of the sagital plane of the MLA_midpoint
    ML_dotP = sum(SagPlane_normal .* (v_below_behind - MLA_midpoint), 2);
    v_side1 = v_below_behind(find(ML_dotP <= 0),:);
    v_side2 = v_below_behind(find(ML_dotP >= 0),:);
    % initialize 5 iterations for finding the optimal points
    for i=1:5
      % calculate the distance for each point of each side from the coronal plane
      % passing through the midpoint_proj point
      d_side1 = sum(CorPlane_normal .* (v_side1 - midpoint_proj),2);
      d_side2 = sum(CorPlane_normal .* (v_side2 - midpoint_proj),2);
      % find new points for optimized mediolateral axis vector that have the
      % maximum distance from the initial coronal plane (optimize local extrema)
      [d1,i1]=min(d_side1);
      MLA_optimal_point1 = v_side1(i1,:);
      [d2,i2]=min(d_side2);
      MLA_optimal_point2 = v_side2(i2,:);
      MLA_vector = MLA_optimal_point1 - MLA_optimal_point2;
      % calculate the optimized normal of the coronal plane
      CorPlane_normal = cross(PDA_vector, MLA_vector);
      % normalize the optimized normal vector of the coronal plane
      CorPlane_normal = CorPlane_normal ./ sqrt(sum(CorPlane_normal.^2));
      % check the optimized CorPlane_normal points towards the right direction
      MLP1_C5 = Centroid_5 - MLA_optimal_point1;
      if sum(MLP1_C5 .* CorPlane_normal) < 0      % MLA_points should be behind
        CorPlane_normal = CorPlane_normal * -1;
      endif
    endfor
  % for tibia bone
  elseif strcmp(type, "Tibia")
    MLP1_C1 = Centroid_1 - MLA_points(1,:);
    if sum(MLP1_C1 .* CorPlane_normal) < 0      % MLA_points should be behind
      CorPlane_normal = CorPlane_normal * -1;
    endif
    % find the points that lie on or above the transverse plane of the proximal centroid
    TP_dotP = sum(TransPlane_normal .* (v - Centroid_1), 2);
    v_above = v(find(TP_dotP >= 0),:);
    % find the midpoint between the user defined points and project it to the
    % vertical vector defined by centroids 5 and 1 to find the closest point along
    % the vertical axis
    MLA_midpoint = (MLA_points(1,:) + MLA_points(2,:)) ./ 2;
    C5_midpoint = MLA_midpoint - Centroid_5;
    C5_C1 = Centroid_1 - Centroid_5;
    midpoint_proj = Centroid_5 + (dot(C5_midpoint,C5_C1)/dot(C5_C1,C5_C1)) * C5_C1;
    % recalculate the coronal plane unit vector by normalizing the vector from
    % MLA_midpoint to midpoint_proj
    CorPlane_normal = midpoint_proj - MLA_midpoint;
    CorPlane_normal = CorPlane_normal ./ sqrt(sum(CorPlane_normal.^2));
    % calculate the unit normal of the sagital plane facing to the right side
    SagPlane_normal = cross(CorPlane_normal, TransPlane_normal);
    % find the points that lie behind the coronal plane of the midpoint_proj point
    CP_dotP = sum(CorPlane_normal .* (v_above - midpoint_proj), 2);
    v_above_behind = v_above(find(CP_dotP < 0),:);
    % find the points that lie on one or the other side
    % of the sagital plane of the MLA_midpoint
    ML_dotP = sum(SagPlane_normal .* (v_above_behind - MLA_midpoint), 2);
    v_side1 = v_above_behind(find(ML_dotP <= 0),:);
    v_side2 = v_above_behind(find(ML_dotP >= 0),:);
    % initialize 5 iterations for finding the optimal points
    for i=1:5
      % calculate the distance for each point of each side from the coronal plane
      % passing through the midpoint_proj point
      d_side1 = sum(CorPlane_normal .* (v_side1 - midpoint_proj),2);
      d_side2 = sum(CorPlane_normal .* (v_side2 - midpoint_proj),2);
      % find new points for optimized mediolateral axis vector that have the
      % maximum distance from the initial coronal plane (optimize local extrema)
      [d1,i1]=min(d_side1);
      MLA_optimal_point1 = v_side1(i1,:);
      [d2,i2]=min(d_side2);
      MLA_optimal_point2 = v_side2(i2,:);
      MLA_vector = MLA_optimal_point1 - MLA_optimal_point2;
      % calculate the optimized normal of the coronal plane
      CorPlane_normal = cross(PDA_vector, MLA_vector);
      % normalize the optimized normal vector of the coronal plane
      CorPlane_normal = CorPlane_normal ./ sqrt(sum(CorPlane_normal.^2));
      % check the optimized CorPlane_normal points towards the right direction
      MLP1_C1 = Centroid_1 - MLA_optimal_point1;
      if sum(MLP1_C1 .* CorPlane_normal) < 0      % MLA_points should be behind
        CorPlane_normal = CorPlane_normal * -1;
      endif
    endfor
  endif
  
	% flush the screen output to display the results during iterations through multiple mesh files
  page_screen_output(0);
  page_output_immediately(1);
  % print initial user defined MLA_points
  printf("\nUser defined MLA points A and B are:\n");
  printf("Ax: %f Ay: %f Az: %f and Bx: %f By: %f Bz: %f\n", MLA_points(1,:), MLA_points(2,:));
  
	% find which optimized MLA points correspond to the user defined pair
  d1 = distancePoints(MLA_points(1,:), MLA_optimal_point1);
  d2 = distancePoints(MLA_points(1,:), MLA_optimal_point2);
  if d1 > d2
    MLA_opt_point_A = MLA_optimal_point2;
    MLA_opt_point_B = MLA_optimal_point1;
  else
    MLA_opt_point_A = MLA_optimal_point1;
    MLA_opt_point_B = MLA_optimal_point2;
  endif
  % print optimized MLA_points
  printf("\nOptimized MLA points A and B are:\n");
  printf("Ax: %f Ay: %f Az: %f and Bx: %f By: %f Bz: %f\n", MLA_opt_point_A, MLA_opt_point_B);
  % if orientation points were retrived from a Meshlab point file
  % save user defined and optimized MLA point back to the Meshlab Point file
  if nargin != 3
    MLP = [MLA_points([1:2],:); MLA_opt_point_A; MLA_opt_point_B];
    write_MeshlabPoints(filenamePP, filename, MLP);
  endif
	
  % calculate new normals for each final cross section at 20, 35, 50, 65 and 80%
  % all normals should be pointing upwards, that is from distal towards proximal
  n1 = (Centroid_1 - Centroid_2)./sqrt(sum((Centroid_1 - Centroid_2).^2));
  n2 = (Centroid_1 - Centroid_3)./sqrt(sum((Centroid_1 - Centroid_3).^2));
  n3 = (Centroid_2 - Centroid_4)./sqrt(sum((Centroid_2 - Centroid_4).^2));
  n4 = (Centroid_3 - Centroid_5)./sqrt(sum((Centroid_3 - Centroid_5).^2));
  n5 = (Centroid_4 - Centroid_5)./sqrt(sum((Centroid_4 - Centroid_5).^2));
	
  % calculate the sectioning points for each plane
  section_1 = slice_Mesh_Plane(v,f,Centroid_1,n1);
  section_2 = slice_Mesh_Plane(v,f,Centroid_2,n2);
  section_3 = slice_Mesh_Plane(v,f,Centroid_3,n3);
  section_4 = slice_Mesh_Plane(v,f,Centroid_4,n4);
  section_5 = slice_Mesh_Plane(v,f,Centroid_5,n5);
	
  % for each cross section calculate the 3D coordinates of its centroid and area
  [CS_Geometry(1), SMoA(1), polyline(1)] = simple_polygon3D(section_1, n1, CorPlane_normal);
  [CS_Geometry(2), SMoA(2), polyline(2)] = simple_polygon3D(section_2, n2, CorPlane_normal);
  [CS_Geometry(3), SMoA(3), polyline(3)] = simple_polygon3D(section_3, n3, CorPlane_normal);
  [CS_Geometry(4), SMoA(4), polyline(4)] = simple_polygon3D(section_4, n4, CorPlane_normal);
  [CS_Geometry(5), SMoA(5), polyline(5)] = simple_polygon3D(section_5, n5, CorPlane_normal);

  % store the normals for each sectioning plane and the normal of the coronal plane
  CS_Geometry(1).Slice_n = n1; CS_Geometry(1).Coronal_n = CorPlane_normal;
  CS_Geometry(2).Slice_n = n2; CS_Geometry(2).Coronal_n = CorPlane_normal;
  CS_Geometry(3).Slice_n = n3; CS_Geometry(3).Coronal_n = CorPlane_normal;
  CS_Geometry(4).Slice_n = n4; CS_Geometry(4).Coronal_n = CorPlane_normal;
  CS_Geometry(5).Slice_n = n5; CS_Geometry(5).Coronal_n = CorPlane_normal;
    
  % print results for centroids and cross sectional areas
  printf("\n%s in %s has a maximum distance of %f mm\n\n", type, filename, maxDistance);
  printf("Cross section at 20%% has an area of %f mm2, perimeter of %f mm \n\
and centroid coordinates are: x:%f y:%f z:%f\n\n",...
          CS_Geometry(1).Area, CS_Geometry(1).Perimeter, CS_Geometry(1).Centroid);
  printf("Cross section at 35%% has an area of %f mm2, perimeter of %f mm \n\
and centroid coordinates are: x:%f y:%f z:%f\n\n",...
          CS_Geometry(2).Area, CS_Geometry(2).Perimeter, CS_Geometry(2).Centroid);
  printf("Cross section at 50%% has an area of %f mm2, perimeter of %f mm \n\
and centroid coordinates are: x:%f y:%f z:%f\n\n",...
          CS_Geometry(3).Area, CS_Geometry(3).Perimeter, CS_Geometry(3).Centroid);
  printf("Cross section at 65%% has an area of %f mm2, perimeter of %f mm \n\
and centroid coordinates are: x:%f y:%f z:%f\n\n",...
          CS_Geometry(4).Area, CS_Geometry(4).Perimeter, CS_Geometry(4).Centroid);
  printf("Cross section at 80%% has an area of %f mm2, perimeter of %f mm \n\
and centroid coordinates are: x:%f y:%f z:%f\n\n",...
          CS_Geometry(5).Area, CS_Geometry(5).Perimeter, CS_Geometry(5).Centroid);
endfunction