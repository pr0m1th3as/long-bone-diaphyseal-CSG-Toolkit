# long-bone-diaphyseal-CSG-Toolkit
A GNU Octave toolkit for analyzing diaphyseal long bone cross sectional geometry

The present set of GNU Octave functions provides a novel and robust algorithm for
analyzing the diaphyseal cross-sectional geometric properties of long bones, which
can be applied to any 3D digital model of a humerus, femur or tibia bone represented
as a triangular mesh in a Wavefront OBJ file format.

The CSG Toolkit requires the 'io', 'statistical' and 'geometry' packages to be
installed. Reading the .obj mesh files is performed with the readObj.oct file,
which may be used directly as a binary. In case this is not achieved (due to
different version of GNU Octave installed or else), the user may directly recompile
it from source with the following command

>> mkoctfile readObj.cc

