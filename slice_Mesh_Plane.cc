/*
Copyright (C) 2018-2020 Andreas Bertsatos <abertsatos@biol.uoa.gr>

This program is free software; you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free Software
Foundation; either version 3 of the License, or (at your option) any later
version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
details.

You should have received a copy of the GNU General Public License along with
this program; if not, see <http://www.gnu.org/licenses/>.
*/

#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <octave/oct.h>
#include <octave/parse.h>

using namespace std;

struct VCoord						//for point and normal coordinates
{
	double x, y, z;
};
struct Mesh             //for triangular mesh object data
{
	double V1x, V1y, V1z;
	double V2x, V2y, V2z;
	double V3x, V3y, V3z;
};
// function for calculating DotProduct between two vectors
double dotProduct(VCoord A, VCoord B)
{
	double dotprod = A.x * B.x + A.y * B.y + A.z * B.z;
	return dotprod;
}
// function for removing duplicate points
vector<VCoord> unique_rows(vector<VCoord> sectionPoints)
{
	vector<VCoord> uniquePoints;
	for(int i = 0; i < sectionPoints.size() - 1; ++i)
	{
		for(int j = i + 1; j < sectionPoints.size(); ++j)
		{
			double x_i = sectionPoints[i].x; double x_j = sectionPoints[j].x;
			double y_i = sectionPoints[i].y; double y_j = sectionPoints[j].y;
			double z_i = sectionPoints[i].z; double z_j = sectionPoints[j].z;
			if(x_i == x_j && y_i == y_j && z_i == z_j)
			{
				sectionPoints.erase(sectionPoints.begin() + j);
				sectionPoints.push_back({x_i, y_i, z_i});
			}
		}
	}
	bool duplicate_found = false;
	int i = 0;
	while(!duplicate_found)
	{
		for(int j = i + 1; j < sectionPoints.size(); ++j)
		{
			double x_i = sectionPoints[i].x; double x_j = sectionPoints[j].x;
			double y_i = sectionPoints[i].y; double y_j = sectionPoints[j].y;
			double z_i = sectionPoints[i].z; double z_j = sectionPoints[j].z;
			if(x_i == x_j && y_i == y_j && z_i == z_j)
			{
				sectionPoints.erase(sectionPoints.begin() + j, sectionPoints.end());
				duplicate_found = true;
			}
		}
		i++;
	}
	uniquePoints = sectionPoints;
	return uniquePoints;
}
// function for extracting intersection points ginen a plane specified by its normal and a point
vector<VCoord> sliceMesh(vector<Mesh> Mesh3D, VCoord point, VCoord normal)
{
	// initialize vectors and variables
	int iFaces = 0;
	vector<VCoord> iPoints;
	VCoord V1, V2, V3;
	double dotV1, dotV2, dotV3;
	// for each triplet of vertices find which face is intersected by specified plane
	for(int i = 0; i < Mesh3D.size(); ++i)
	{
		V1 = {Mesh3D[i].V1x - point.x, Mesh3D[i].V1y - point.y, Mesh3D[i].V1z - point.z};
		V2 = {Mesh3D[i].V2x - point.x, Mesh3D[i].V2y - point.y, Mesh3D[i].V2z - point.z};
		V3 = {Mesh3D[i].V3x - point.x, Mesh3D[i].V3y - point.y, Mesh3D[i].V3z - point.z};
		dotV1 = dotProduct(V1, normal);
		dotV2 = dotProduct(V2, normal);
		dotV3 = dotProduct(V3, normal);
		if(!((dotV1 < 0 && dotV2 < 0 && dotV3 < 0) || (dotV1 > 0 && dotV2 > 0 && dotV3 > 0)))   // true if intersecting
		{
			if(dotV1 * dotV2 < 0)       // vertices 1 and 2 of ith face lie on opposite sides
			{
				VCoord PL = {point.x - Mesh3D[i].V1x, point.y - Mesh3D[i].V1y, point.z - Mesh3D[i].V1z};
				VCoord L = {Mesh3D[i].V2x - Mesh3D[i].V1x, Mesh3D[i].V2y - Mesh3D[i].V1y, Mesh3D[i].V2z - Mesh3D[i].V1z};
				double d = dotProduct(PL, normal) / dotProduct(L, normal);
				VCoord CSv1v2 = {Mesh3D[i].V1x + d * L.x, Mesh3D[i].V1y + d * L.y, Mesh3D[i].V1z + d * L.z};
				iPoints.push_back(CSv1v2);
			}
			if(dotV1 * dotV3 < 0)       // vertices 1 and 3 of ith face lie on opposite sides
			{
				VCoord PL = {point.x - Mesh3D[i].V1x, point.y - Mesh3D[i].V1y, point.z - Mesh3D[i].V1z};
				VCoord L = {Mesh3D[i].V3x - Mesh3D[i].V1x, Mesh3D[i].V3y - Mesh3D[i].V1y, Mesh3D[i].V3z - Mesh3D[i].V1z};
				double d = dotProduct(PL, normal) / dotProduct(L, normal);
				VCoord CSv1v3 = {Mesh3D[i].V1x + d * L.x, Mesh3D[i].V1y + d * L.y, Mesh3D[i].V1z + d * L.z};
				iPoints.push_back(CSv1v3);
			}
			if(dotV2 * dotV3 < 0)       // vertices 2 and 3 of ith face lie on opposite sides
			{
				VCoord PL = {point.x - Mesh3D[i].V2x, point.y - Mesh3D[i].V2y, point.z - Mesh3D[i].V2z};
				VCoord L = {Mesh3D[i].V3x - Mesh3D[i].V2x, Mesh3D[i].V3y - Mesh3D[i].V2y, Mesh3D[i].V3z - Mesh3D[i].V2z};
				double d = dotProduct(PL, normal) / dotProduct(L, normal);
				VCoord CSv2v3 = {Mesh3D[i].V2x + d * L.x, Mesh3D[i].V2y + d * L.y, Mesh3D[i].V2z + d * L.z};
				iPoints.push_back(CSv2v3);
			}
			if(dotV1 == 0)      // vertex 1 of ith face lies on the cross-section plane
			{
				VCoord CSv1;
				CSv1 = {Mesh3D[i].V1x, Mesh3D[i].V1y, Mesh3D[i].V1z};
				iPoints.push_back(CSv1);
			}
			if(dotV2 == 0)      // vertex 2 of ith face lies on the cross-section plane
			{
				VCoord CSv2;
				CSv2 = {Mesh3D[i].V2x, Mesh3D[i].V2y, Mesh3D[i].V2z};
				iPoints.push_back(CSv2);
			}
			if(dotV3 == 0)      // vertex 2 of ith face lies on the cross-section plane
			{
				VCoord CSv3;
				CSv3 = {Mesh3D[i].V3x, Mesh3D[i].V3y, Mesh3D[i].V3z};
				iPoints.push_back(CSv3);
			}
		}
	}
	return iPoints;
}
DEFUN_DLD(slice_Mesh_Plane, args, nargout, 
          "-*- texinfo -*-\n\
@deftypefn{function} @var{cross_section} = slice_Mesh_Plane(@var{v},@var{f},@var{point},@var{normal})\n\
\n\
\n\
This function loads the vertices @var{v} and faces @var{f} of a triangular 3D Mesh along with a\
slicing plane defined by its @var{normal} and a @var{point} that lies on the slicing plane and\
returns the intersection points of the face edges between the vertices that lie on opposite sides\
of the sectioning plane. Duplicate points due to adjacent faces are removed.\n\n\
@var{v} and @var{f} need be (Nx3) matrices containing 3D coordinates and vertex indices respectively. \
@var{normal} and @var{point} should be defined as row vectors containing x, y, z coordinates in R3. \
@var{cross_section} is the return variable Nx3 in size, where N is the number of unique intersection\
of the 3D mesh represented by 'v' and 'f' input arguments.\n\
@end deftypefn")
{
  // count the number of input arguments and store their values into the appropriate variables
  // check for invalid number of input arguments
  if (args.length() != 4)
  {
    std::cout << "Invalid number of input arguments.\n";
    return octave_value_list();
  }
  // check all arguments being real matrices
  if (!args(0).is_matrix_type() || !args(1).is_matrix_type() || !args(2).is_matrix_type() || !args(3).is_matrix_type())
  {
    std::cout << "All input arguments should be real matrices.\n";
    return octave_value_list();
  }
  // check vertices containing 3D coordinates
	if (args(0).columns() != 3)
	{
		std::cout << "Vertex matrix should be Nx3 containing x,y,z coordinates.\n";
    return octave_value_list();
  }
	// check mesh for being triangular
	if (args(1).columns() != 3)
	{
		std::cout << "Mesh should be triangular. Face matrix should be Nx3 containing three vertices.\n";
    return octave_value_list();
  }
	// check point and normal being 1x3 arrays
	if (args(2).rows() != 1 || args(2).columns() != 3)
	{
		std::cout << "point should be 1x3 array containing x y z coordinates.\n";
    return octave_value_list();
  }
	if (args(3).rows() != 1 || args(3).columns() != 3)
	{
		std::cout << "normal should be 1x3 array containing x y z coordinates.\n";
    return octave_value_list();
  }
	// store vertices and faces in Mesh data structure
  Matrix V = args(0).array_value();
  Matrix F = args(1).array_value();
	vector<Mesh> Mesh3D;
  // find number of vertices and faces
  octave_idx_type V_rows = args(0).rows();
  octave_idx_type F_rows = args(1).rows();
	// for each face store its corresponding vertex coordinates
	double tmpV1x, tmpV1y, tmpV1z;
	double tmpV2x, tmpV2y, tmpV2z;
	double tmpV3x, tmpV3y, tmpV3z;
	for (octave_idx_type i = 0; i < F_rows; i++)
	{
		tmpV1x = V(F(i,0) - 1, 0); tmpV1y = V(F(i,0) - 1, 1); tmpV1z = V(F(i,0) - 1, 2);
		tmpV2x = V(F(i,1) - 1, 0); tmpV2y = V(F(i,1) - 1, 1); tmpV2z = V(F(i,1) - 1, 2);
		tmpV3x = V(F(i,2) - 1, 0); tmpV3y = V(F(i,2) - 1, 1); tmpV3z = V(F(i,2) - 1, 2);
		Mesh temp_mesh = {tmpV1x, tmpV1y, tmpV1z, tmpV2x, tmpV2y, tmpV2z, tmpV3x, tmpV3y, tmpV3z};
		Mesh3D.push_back(temp_mesh);
	}	
	// store point and normal from input arguments
	Matrix P = args(2).array_value();
	Matrix N = args(3).array_value();
	VCoord point = {P(0), P(1), P(2)};
	VCoord normal = {N(0), N(1), N(2)};
	// calculate cross-sectioning points
	vector<VCoord> section_points;
	section_points = sliceMesh(Mesh3D, point, normal);
	vector<VCoord> unique_points;
	unique_points = unique_rows(section_points);
	// cast the cross-sectioning points into an octave type Matrix
	int n = unique_points.size();	
	Matrix cross_section (n, 3);
	for (octave_idx_type i = 0; i < n; i++)
	{
		cross_section(i,0) = unique_points[i].x;
		cross_section(i,1) = unique_points[i].y;
		cross_section(i,2) = unique_points[i].z;
	}
	// define return value list and return the cross-sectioning points
  octave_value_list retval;
	retval(0) = cross_section;
  return retval;
}