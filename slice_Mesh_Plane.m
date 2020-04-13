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
function section = slice_Mesh_Plane(v,f,point,normal)
  % function [section] = slice_Mesh_Plane(v, f, point, normal)
  %
  % This function slices a triangular mesh defined by 'v' (vertices) and 'f' (faces)
  % with a plane defined by its normal and a point that lies on the slicing plane
  % and returns the intersection points of the face edges between the vertices
  % that lie on opposite sides of the sectioning plane. Duplicate points due to
  % adjacent faces are removed. Note that the sectioning points are returned in
  % random order.
  % 'v' and 'f' need be (Nx3) matrices containing 3D coordinates and vertex indices
  % respectively. 'point' and 'normal' should be defined as row vectors containing
  % x, y, z coordinates in R3. 'section' is the return variable Nx3 in size,
  % where N is the number of unique intersection points of the 3D mesh represented
  % by 'v' and 'f' input arguments.

  % calculate the dot products for each vertex of the mesh so we can search which
  % faces are intersecting the plane and which are not
  section = sum((v - point) .* normal, 2);
  % search for faces intersecting each sectioning plane
  face = 0;
  for i=1:length(f)
    face_V1 = f(i,1);
    face_V2 = f(i,2);
    face_V3 = f(i,3);
    % plane at 20%
    if !((section(face_V1) < 0 && section(face_V2) < 0 && section(face_V3) < 0)...
        || (section(face_V1) > 0 && section(face_V2) > 0 && section(face_V3) > 0))
      face += 1;
      section_faces(face,:) = [face_V1, face_V2, face_V3];
    endif
  endfor
  % calculate the intersection points of the face edges between the face
  % vertices that lie on opposite sides of the sectioning plane
  sp = 0;
  for i=1:length(section_faces)
    % check if vertices A and B lie on opposite sides
    if (section(section_faces(i,1)) * section(section_faces(i,2))) < 0
      sp += 1;
      P_0 = v(section_faces(i,1),:);
      P_1 = v(section_faces(i,2),:);
      u = P_1 - P_0;
      w = P_0 - point;
      D = sum(normal .* u);
      N = -sum(normal .* w);
      sI = N / D;
      section_CS(sp,:) = P_0 + (sI .* u);
    endif
    % check if vertices A and C lie on opposite sides
    if (section(section_faces(i,1)) * section(section_faces(i,3))) < 0
      sp += 1;
      P_0 = v(section_faces(i,1),:);
      P_1 = v(section_faces(i,3),:);
      u = P_1 - P_0;
      w = P_0 - point;
      D = sum(normal .* u);
      N = -sum(normal .* w);
      sI = N / D;
      section_CS(sp,:) = P_0 + (sI .* u);
    endif
    % check if vertices B and C lie on opposite sides
    if (section(section_faces(i,2)) * section(section_faces(i,3))) < 0
      sp += 1;
      P_0 = v(section_faces(i,2),:);
      P_1 = v(section_faces(i,3),:);
      u = P_1 - P_0;
      w = P_0 - point;
      D = sum(normal .* u);
      N = -sum(normal .* w);
      sI = N / D;
      section_CS(sp,:) = P_0 + (sI .* u);
    endif
  endfor
  % clear the duplicate points from each cross section
  section = unique(section_CS,"rows");
endfunction

