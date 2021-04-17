% Copyright (C) 2018-2021 Andreas Bertsatos <abertsatos@biol.uoa.gr>
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
function [varargout] = visualize_CrossSections(varargin)
  % function visualize_CrossSections("bone_id", CS_Geometry, SMoA, polyline)
  % function [CS_Geometry, SMoA, polyline] = visualize_CrossSections("bone_id")
  %
  % This function reads the 'geometry-$$.csv', 'inertia-$$.csv' and 'polyline2D-$$.csv'
  % files produced with the 'longbone_Analysis.m' script, where $$ is the 'bone_id'
  % char string for the required bone, and plots the 2D cross sectional polygons
  % from proximal (top) to distal (bottom) along with certain information regarding
  % cross sectional area, perimeter, and second moments of area.
  %
  % The centroid of each cross section is centered at origin and the x-axis of
  % the plots represents the frontal axis of the bone from left to right side,
  % whereas the y-axis is aligned to the sagital axis with the top pointing
  % towards the front side (anterior). Both axes retain equal size so that the
  % actual shape of the cross sectional area is preserved. However, each figure
  % has a different scaling factor. Consequently, size may not be visually
  % proportional among different figures, but the axes in each figure preserve
  % the values of the actual size for the displayed cross section.
  %
  % All three .csv files for the respective bone defined in the 'bone_id' input
  % argument must be present in the working directory. Alternatively, the data
  % structures 'CS_Geometry', 'SMoA' and 'polyline' may be parsed to the function
  % instead. In such case, the fields of each structure need to comply with the
  % the standards used by the 'longbone_geometry.m' function.
  %
  % If output arguments are defined then the function returns the aforementioned
  % data structures with the data retrieved from the .csv files. For 1 output
  % argument 'CS_Geometry' is returned, for 2 output arguments 'CS_Geometry' and
  % 'SMoA' are returned. When 3 output arguments are defined, then 'polyline'
  % structure is also returned, but contains only the 'poly2D' field with the
  % 2D coordinates of the planar cross sections.
  %
  % Function 'visualize_CrossSection.m' relies on 'io' and 'matgeom' packages
  % that need to be installed.
  
  % load io and geometry packages
  pkg load io;
  pkg load matgeom;
  % check if only 'bone_id' is given to the function and if this is the case
  % scan the files of the working directory for the required .csv files and
  % check if all three .csv files are present.
  if (nargin == 1)
    bone_id = varargin{1};
    filenames = readdir(pwd);
    % merge filenames
    g_filename = strcat("geometry-", bone_id, ".csv");
    i_filename = strcat("inertia-", bone_id, ".csv");
    p_filename = strcat("polyline2D-", bone_id, ".csv");
    % read files if available
    if (sum(strcmp(filenames, g_filename)) == 1)
      geometry = csvread(filenames(strcmp(filenames, g_filename)){:});
    endif
    if (sum(strcmp(filenames, i_filename)) == 1)
      inertia = csvread(filenames(strcmp(filenames, i_filename)){:});
    endif
    if (sum(strcmp(filenames, p_filename)) == 1)
      poly = csvread(filenames(strcmp(filenames, p_filename)){:});
    endif
    % check if all files present
    if !(exist("geometry") == 1 && exist("geometry") == 1 && exist("poly") == 1)
      error 'not all .csv files are present';
    else
      % store the matrices from the .csv files to the corresponding
      % structures in the same format as returned from 'longbone_Geometry.m'
      for i=1:5
        CS_Geometry(i).Area = geometry(i,1);
        CS_Geometry(i).Perimeter = geometry(i,2);
        SMoA(i).Ix = inertia(i,1);
        SMoA(i).Iy = inertia(i,2);
        SMoA(i).Ixy = inertia(i,3);
        SMoA(i).Imin = inertia(i,4);
        SMoA(i).Imax = inertia(i,5);
        SMoA(i).theta = inertia(i,6);
        % remove the trailing zeros from each columns
        poly_x = poly(:,(i*2)-1);
        poly_y = poly(:,(i*2));
        lnz_x = find(poly_x, 1, 'last');
        lnz_y = find(poly_y, 1, 'last');
        last_index = max(lnz_x, lnz_y);
        polyX = poly_x([1:last_index],:);
        polyY = poly_y([1:last_index],:);
        polyline(i).poly2D = [polyX, polyY];
        clear poly_x; clear poly_y; clear lnz_x; clear lnz_y;
      endfor
    endif
  % alternatively, check if three input arguments are parsed into the function
  % and check that they are structures
  elseif (nargin == 4)
    bone_id = varargin{1};
    CS_Geometry = varargin{2};
    SMoA = varargin{3};
    polyline = varargin{4};
    if !(isstruct(CS_Geometry) &&  isstruct(SMoA) && isstruct(polyline))
      error 'all three input arguments need to be structures';
    endif
  else
    error 'invalid number of input arguments';
  endif
  a = {" 20%", " 35%", " 50%", " 65%", " 80%"};
  for i=1:5
    figure(i, "name", bone_id, "numbertitle", "off");
    drawPolygon(polyline(i).poly2D(:,1),polyline(i).poly2D(:,2));
    axis("image");
    title(strcat("Cross section at ",a(i)), 'FontSize', 16);
    xlabel('Frontal axis (mm)','FontSize', 10);
    ylabel('Sagital axis (mm)','FontSize', 10);
    grid("on");
    legend(["Perimeter = ", num2str(CS_Geometry(i).Perimeter,'%5.1f'), " mm",...
            "\nArea = ", num2str(CS_Geometry(i).Area,'%5.1f'), " mm^2",...
            "\nIx = ", num2str(SMoA(i).Ix,'%5.1f'), " mm^4",...
            "\nIy = ", num2str(SMoA(i).Iy,'%5.1f'), " mm^4",...
            "\nIxy = ", num2str(SMoA(i).Ixy,'%5.1f'), " mm^4",...
            "\nImin = ", num2str(SMoA(i).Imin,'%5.1f'), " mm^4",...
            "\nImax = ", num2str(SMoA(i).Imax,'%5.1f'), " mm^4",...
            "\ntheta = ", num2str(SMoA(i).theta,'%5.1f'), "^o"],...
            'FontSize', 14, 'box', "off", 'location', 'eastoutside');
    pos(i,:) = get(gcf,'Position');
    pos(i,1) = pos(i,1) + (i-1) * 150;
    set(gcf, 'Position', [pos(i,1), pos(i,2), pos(i,3), pos(i,4)]);
  endfor
  % check if output arguments are defined and return the data structures
  if nargout == 1
    varargout{1} = CS_Geometry;
  elseif nargout == 2
    varargout{1} = CS_Geometry;
    varargout{2} = SMoA;
  elseif nargout == 3
    varargout{1} = CS_Geometry;
    varargout{2} = SMoA;
    varargout{3} = polyline;
  endif
endfunction