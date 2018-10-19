# long-bone-diaphyseal-CSG-Toolkit
A GNU Octave toolkit for analyzing diaphyseal long bone cross sectional geometry

The present set of GNU Octave functions provides a novel and robust algorithm for
analyzing the diaphyseal cross-sectional geometric properties of long bones, which
can be applied to any 3D digital model of a humerus, femur or tibia bone represented
as a triangular mesh in a Wavefront OBJ file format.

The CSG Toolkit requires the 'io', 'statistical' and 'geometry' packages to be
installed. Reading the .obj mesh files is performed with the readObj.oct file,
which may be used directly as a binary. In case this is not possible (due to
different version of GNU Octave installed or else), the user may directly recompile
it from source with the following command

>> mkoctfile readObj.cc

The CSG Toolkit can be called in batch processing mode using the 'longbone_Analysis.m'
script which can handle all available 3D models found in OBJ format in the working
directory as long as a specific bone, i.e. humerus, femur or tibia, is selected. Type

>> help longbone_Analysis

for information about its usage. The CSG Toolkit also provides functionality for
graphical representation of the cross-sectional contours and their respective CSG
properties, which can be accessed with the 'visualize_CrossSections.m' function by
reading the output results of the longbone_Analysis stored in relevant CSV files.

The present toolkit has been extensively tested and these results are presented
in the relevant validation study available at doi:10.5281/zenodo.1466135.
The testing dataset, which comprises three 3D mesh models of real humans bones 
(a humerus, a femur and a tibia, which are part of the Athens modern reference
skeletal collection), is also freely available at doi:10.5281/zenodo.1466962 and 
may be used with the CSG Toolkit for demonstrating its operation.
