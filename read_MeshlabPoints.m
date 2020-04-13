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
function [varargout] = read_MeshlabPoints(filename)
  % function [MLP, name_list] = read_MeshlabPoints(filename)
  %
  % This function reads a .pp MeshLab Point file and returns and Nx4 matrix
  % containing the name for each point along with the corresponding x y z 
  % coordinates in each row. If name is given as a numerical value then it is
  % stored as numerical value. If the name of a point is alpharithmetic or 0
  % then NaN is stored in the relevant row of the MPL array and the name is
  % stored in a separate cell array as a character string.
  % 
  % If one output argument is provided then the function returns the Nx4 matrix,
  % if two output arguments are defined then an Nx3 matrix containing only the
  % x, y, z coordinates and a cell array containing the names of each point are
  % returned.
  %
  % Input:
  %  a char string of the filename of Meshlab's Point file
  %
  % Output:
  %   MLP  number of Points x 4 array in the format
  %   MLP(1,:) = [name, x, y, z]
  %
  % or
  %   MLP  number of Points x 3 array in the format
  %   MLP(1,:) = [x, y, z]
  %   name_list(1) = "name of first point"
  %
  %
  MLP = zeros(1,4);
  name_list = {''};
  
  
  point_index = 1;
  fid = fopen(filename,'rt');
  line = fgets(fid);
  while ischar(line)
    if strfind(line, "<point")
      x_start = strfind(line,"x=\"") + 3;
      y_start = strfind(line,"y=\"") + 3;
      z_start = strfind(line,"z=\"") + 3;
      name_start = strfind(line,"name=\"") + 6;
      quotes = strfind(line,"\"");
      x_end = min(quotes(find(quotes>x_start))) - 1;
      y_end = min(quotes(find(quotes>y_start))) - 1;
      z_end = min(quotes(find(quotes>z_start))) - 1;
      name_end = min(quotes(find(quotes>name_start))) - 1;
      x = str2num(line(x_start:x_end));
      y = str2num(line(y_start:y_end));
      z = str2num(line(z_start:z_end));
      if (isnumeric(str2num(line(name_start:name_end))))
        name = str2num(line(name_start:name_end));
        name_list(point_index) = str2num(line(name_start:name_end));
      else
        name = NaN;
        name_list(point_index) = line(name_start:name_end);
      endif
      MLP(point_index,:) = [name x y z];
      point_index += 1;
    endif
    line = fgets(fid);
  end
  fclose(fid);
  if (nargout==1)
    varargout{1} = MLP;
  elseif (nargout==2)
    varargout{1} = MLP(:,[2:4]);
    varargout{2} = name_list;
  endif
endfunction