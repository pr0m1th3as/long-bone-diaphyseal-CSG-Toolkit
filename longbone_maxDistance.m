% Copyright (C) 2018 Andreas Bertsatos <andreas.bertsatos@gmail.com>
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
function [varargout] = longbone_maxDistance(v)
  % function [maxDistance, maxD_v1, maxD_v2] = longbone_maxDistance(v)
  %
  % This function calculates the maximum distance of a long bone as represented
  % by its mesh vertices. It requires an Nx3 matrix as input argument, which
  % contains the (x,y,z) coordinates in R3. The function also returns the point
  % coordinates that correspond to the maximum distance.
  %
  % The function 'longbone_maxDistance.m' may be called with varying number of
  % output arguments.
  % If 1 output argument is defined, then maximum distance is returned.
  %
  %     e.g. maxDistance = longbone_maxDistance(v)
  %
  % If 2 output arguments are defined, then the corresponding points of
  % maximum distance are returned.
  %
  %     e.g. [maxD_v1, maxD_v2] = longbone_maxDistance(v)
  %
  % If 3 output arguments are defined, then the maximum distance along with
  % its corresponding points are returned.
  %
  %     e.g. [maxDistance, maxD_v1, maxD_v2] = longbone_maxDistance(v)
  %
  % The function requires the 'geometry' and 'statistics' packages to be loaded.
  
  % find the extreme vertices of the mesh for each axis
  [w,iw] = min(v);
  for i=1:3
    extreme_V(i,:) = v(iw(i),:);
    extreme_Vi(i) = iw(i);
  endfor
  [w,iw] = max(v);
  for i=1:3
    extreme_V(i+3,:) = v(iw(i),:);
    extreme_Vi(i+3) = iw(i);
  endfor

  % find max distance between pair of vertices
  d = squareform(pdist(extreme_V));
  max_d = max(max(d));
  % find corresponding vertices
  [q1, q2] = find(d==max_d);
  % check if indices match
  if q1(1)==q2(2) && q1(2)==q2(1)
    maxd_V1 = extreme_V(q1(1),:);
    maxd_V2 = extreme_V(q1(2),:);
  else
    printf("There is something wrong! Check your input data.");
  endif
  % find vertices that lie outside the boundaries of initial most distant vertices
  % iterate three times to ensure that maximum distance is found
  len = distancePoints(maxd_V1,v);
  [w,iw1] = max(len);
  maxd_V2 = v(iw1,:);
  max_d = distancePoints(maxd_V1, maxd_V2);
  len = distancePoints(maxd_V2,v);
  [w,iw2] = max(len);
  maxd_V1 = v(iw2,:);
  max_d = distancePoints(maxd_V1, maxd_V2);
  len = distancePoints(maxd_V1,v);
  [w,iw1] = max(len);
  maxd_V2 = v(iw1,:);
  max_d = distancePoints(maxd_V1, maxd_V2);
  len = distancePoints(maxd_V2,v);
  [w,iw2] = max(len);
  maxd_V1 = v(iw2,:);
  max_d = distancePoints(maxd_V1, maxd_V2);
  len = distancePoints(maxd_V1,v);
  [w,iw1] = max(len);
  maxd_V2 = v(iw1,:);
  max_d = distancePoints(maxd_V1, maxd_V2);
  len = distancePoints(maxd_V2,v);
  [w,iw2] = max(len);
  maxd_V1 = v(iw2,:);
  max_d = distancePoints(maxd_V1, maxd_V2);
  if (nargout==0)
    printf("Bone maximum distance is %f and is found between vertices %d and %d.\n",...
            max_d, iw1, iw2);
    printf("3D coordinates of most distant vertices are:\n\nV1 = %f %f %f\nV2 = %f %f %f\n\n",...
            maxd_V1, maxd_V2);
  elseif (nargout==1)
    varargout{1} = max_d;
  elseif (nargout==2)
    varargout{1} = maxd_V1;
    varargout{2} = maxd_V2;
  else
    varargout{1} = max_d;
    varargout{2} = maxd_V1;
    varargout{3} = maxd_V2;
  endif
endfunction



