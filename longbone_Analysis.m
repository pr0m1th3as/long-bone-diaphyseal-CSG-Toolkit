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

% This script reads the available triangular meshes stored in .obj format that
% are present in the working folder and utilizes the 'longbone_Geometry.m' function
% to analyze their geometric properties, which are subsequently stored in .csv
% files named following the name convention of their initial mesh filename.
%
% For example, for a triangular mesh named as 'ID_humerus.obj', the following
% files are produced:
%
%   geometry-ID_humerus.csv     # containing the properties of area, perimeter,
%                               # centroid, slicing and orientation normals for
%                               # each cross section stored as a row vector in
%                               # the order (:,[1],[2],[3:5],[6:8],[9:11]) from
%                               # proximal to distal cross sections ([1:5],:)
%                               
%   inertia-ID_humerus.csv      # containing the properties of Ix, Iy, Ixy, Imin,
%                               # Imax and theta angle for each cross section as
%                               # a row vector (:,[1:6]) from proximal to distal
%                               # cross sections ([1:5],:)
%                               
%   polyline2D-ID_humerus.csv   # containing the 2D coordinates for each cross
%                               # section as an Nx2 matrix, where N is the number
%                               # of points for each polygon. The polygons are
%                               # ordered as (:,[1:2],[3:4],[5:6],[7:8],[9:10])
%                               # from proximal to distal cross sections
%
%   polyline3D-ID_humerus.csv   # containing the 3D coordinates for each cross
%                               # section as an Nx3 matrix, where N is the number
%                               # of points for each polygon. The polygons are
%                               # ordered as (:,[1:3],[4:6],[7:9],[10:12],[13:15])
%                               # from proximal to distal cross sections
%
% The script requires that 'statistcs', 'geometry' and 'io' packages are
% installed and for each $$.obj file there is a corresponding $$.pp Meshlab Point
% file, which contains two points that provide guidance for the anatomical
% orientation of the bone represented in the triangular mesh.
%
% The user is promted with the list of longbones and must choose the appropriate
% bone, that is 'Humerus', 'Femur' or 'Tibia', contained in the obj files so
% that the underlying functions can properly calculate its orientation. Only one
% type of bone can be analyzed each time, so the user must ensure that only
% meshes of a particular bone are present in the working folder during an analysis.
% Of course, both sides of the same bone can be analyzed together, but each
% .obj file must explicitly contain a single bone. Coordinates of the triangular
% meshes are considered to be in mm.
%
% This script requires the 'io', 'matgeom' & 'statistics' packages to be loaded.
% It also relies on the functions 'longbone_Geometry', 'longbone_maxDistance',
% 'slice_Mesh_Plane', 'simple_polygon3D', 'read_MeshlabPoints', 'write_MeshlabPoints',
% and 'readObj', which must be present in the working directory.

% define the type of bone contained in the mesh files
bone_type = questdlg("Select bone", "Bone Selection", "Humerus", "Femur", "Tibia", "Humerus");

% load required packages
pkg load statistics
pkg load matgeom
pkg load io;

% list the filenames with .obj extension in the working folder
filenames = dir("*.obj");

% calculate the geometric properties for each 
for i = 1:length(filenames)
  filename = strcat(filenames(i).name);

  [CS_Geometry, SMoA, polyline] = longbone_Geometry(filename, bone_type);

  geometry(:,1) = [CS_Geometry(1).Area; CS_Geometry(2).Area;...
    CS_Geometry(3).Area; CS_Geometry(4).Area; CS_Geometry(5).Area];
  geometry(:,2) = [CS_Geometry(1).Perimeter; CS_Geometry(2).Perimeter;...
    CS_Geometry(3).Perimeter; CS_Geometry(4).Perimeter; CS_Geometry(5).Perimeter];
  geometry(:,[3:5]) = [CS_Geometry(1).Centroid; CS_Geometry(2).Centroid;...
    CS_Geometry(3).Centroid; CS_Geometry(4).Centroid; CS_Geometry(5).Centroid];
  geometry(:,[6:8]) = [CS_Geometry(1).Slice_n; CS_Geometry(2).Slice_n;...
    CS_Geometry(3).Slice_n; CS_Geometry(4).Slice_n; CS_Geometry(5).Slice_n];
  geometry(:,[9:11]) = [CS_Geometry(1).Coronal_n; CS_Geometry(2).Coronal_n;...
    CS_Geometry(3).Coronal_n; CS_Geometry(4).Coronal_n; CS_Geometry(5).Coronal_n];

  inertia(1,:) = [SMoA(1).Ix, SMoA(1).Iy, SMoA(1).Ixy, SMoA(1).Imin, SMoA(1).Imax, SMoA(1).theta];
  inertia(2,:) = [SMoA(2).Ix, SMoA(2).Iy, SMoA(2).Ixy, SMoA(2).Imin, SMoA(2).Imax, SMoA(2).theta];
  inertia(3,:) = [SMoA(3).Ix, SMoA(3).Iy, SMoA(3).Ixy, SMoA(3).Imin, SMoA(3).Imax, SMoA(3).theta];
  inertia(4,:) = [SMoA(4).Ix, SMoA(4).Iy, SMoA(4).Ixy, SMoA(4).Imin, SMoA(4).Imax, SMoA(4).theta];
  inertia(5,:) = [SMoA(5).Ix, SMoA(5).Iy, SMoA(5).Ixy, SMoA(5).Imin, SMoA(5).Imax, SMoA(5).theta];

  polygon2D([1:length(polyline(1).poly2D)],[1:2]) = polyline(1).poly2D;
  polygon2D([1:length(polyline(2).poly2D)],[3:4]) = polyline(2).poly2D;
  polygon2D([1:length(polyline(3).poly2D)],[5:6]) = polyline(3).poly2D;
  polygon2D([1:length(polyline(4).poly2D)],[7:8]) = polyline(4).poly2D;
  polygon2D([1:length(polyline(5).poly2D)],[9:10]) = polyline(5).poly2D;

  polygon3D([1:length(polyline(1).poly3D)],[1:3]) = polyline(1).poly3D;
  polygon3D([1:length(polyline(2).poly3D)],[4:6]) = polyline(2).poly3D;
  polygon3D([1:length(polyline(3).poly3D)],[7:9]) = polyline(3).poly3D;
  polygon3D([1:length(polyline(4).poly3D)],[10:12]) = polyline(4).poly3D;
  polygon3D([1:length(polyline(5).poly3D)],[13:15]) = polyline(5).poly3D;

  starting = "geometry-";
  name = filename([1:length(filename)-4]);
  endfile = ".csv";
  filename = strcat(starting,name,endfile);
  csvwrite(filename,geometry);
  starting = "inertia-";
  filename = strcat(starting,name,endfile);
  csvwrite(filename,inertia);
  starting = "polyline2D-";
  filename = strcat(starting,name,endfile);
  csvwrite(filename,polygon2D);
  starting = "polyline3D-";
  filename = strcat(starting,name,endfile);
  csvwrite(filename,polygon3D);

  clear geometry; clear inertia; clear polygon2D; clear polygon3D; 
  clear CS_Geometry; clear SMoA; clear polyline;
endfor