# long-bone-diaphyseal-CSG-Toolkit
A GNU Octave toolkit for analyzing diaphyseal long bone cross sectional geometry

The present set of GNU Octave functions provides a novel and robust algorithm for
analyzing the diaphyseal cross-sectional geometric properties of long bones, which
can be applied to any 3D digital model of a humerus, femur or tibia bone represented
as a triangular mesh in a Wavefront OBJ file format.

The CSG Toolkit requires the 'io >= 2.4.12', 'statistical >= 1.4.1' and 'matgeom >= 1.2.1'
packages to be installed. Reading the .obj mesh files and calculating the cross sections
relies on the 'readObj', 'slice_Mesh_Plane', and 'longbone_Registration' functions,
which should be compiled before their initial use.
These can be compiled from within GNU Octave environment by issuing the following commands

	mkoctfile readObj.cc
	mkoctfile slice_Mesh_Plane.cc
	mkoctfile longbone_Registration.cc

Alternatively, you may download and install the CSG Toolkit as a GNU Octave package from
the archive 'csg-toolkit-1.2.1.tar.gz' which corresponds to the latest version. When installing
the packaged format, compiling and loaded paths are handled automatically. Download and run

	pkg install csg-toolkit-1.2.1.tar.gz

to install from a locally downloaded archive, or issue the following command

	pkg install "https://raw.githubusercontent.com/pr0m1th3as/long-bone-diaphyseal-CSG-Toolkit/master/csg-toolkit-1.2.1.tar.gz"

to automatically download and install it from online repository. Prior to using the CSG Toolkit,
remember to run

	pkg load csg-toolkit

Dependencies on aforementioned packages are still valid prior to installation , but they are
automatically loaded along with the csg-toolkit. Happy long bone analysis!

The CSG Toolkit can be called in batch processing mode using the 'longbone_Analysis'
script which can handle all available 3D models found in OBJ format in the working
directory or any other folder path. Although this functionality is still available,
there is no longer need to work with a specific bone being present in a folder for batch
proccessing. The user may select a specific bone, i.e. Humerus, or may choose to process
different bones in a single batch process. Not selected bones are not processed.
The need for MeshLab .pp side car files has been also relaxed. If they are present along
with their OBJ counterparts, they are utilized, if not, the initial alignment points are
automatically registered with the 'longbone_Registration' function. Type

	help longbone_Analysis

for information about its usage. The CSG Toolkit also provides functionality for
graphical representation of the cross-sectional contours and their respective CSG
properties, which can be accessed with the 'visualize_CrossSections.m' function by
reading the output results of the longbone_Analysis stored in relevant CSV files. When
numerous samples have been analyzed in batch mode, the function 'inspect_CSG.m' can also
be used to facilitate their visual inspection according to the files present in the working
directory. Furthermore, 'inspect_CSG.m' will assemble all calculated values in tabular
form and save them in a CSV file. Type 

	help inspect_CSG

for more information about how to use 'inspect_CSG.m' and the layout of the resulting CSV file.

The present toolkit has been extensively tested and these results are presented
in the relevant validation study available at doi:10.5281/zenodo.1466135.
The testing dataset, which comprises three 3D mesh models of real humans bones 
(a humerus, a femur and a tibia, which are part of the Athens modern reference
skeletal collection), is also freely available at doi:10.5281/zenodo.1466962 and 
may be used with the CSG Toolkit for demonstrating its operation. Please, when using this toolkit,
cite the following reference:

Bertsatos A, Chovalopoulou M-E. 2019. A novel method for analyzing long bone diaphyseal 
cross-sectional geometry. A GNU Octave CSG Toolkit. Forensic Science International 297: 65–71. 
DOI: 10.1016/j.forsciint.2019.01.041

