# The Intel MPI compilers have problems when a C-header is included before the
# MPI header. This script inserts an include of the MPI header into every .cpp
# file, and guards it by #ifndef PRECICE_NO_MPI. This is sufficient to cirmcumvent
# the problem. The script is meant to be executed on the fly and the changes
# it makes are meant to be included only temporarily.
#
# How to use this script: Execute it in the preCICE root folder, i.e., above the
# src folder.

import os
import glob
import shutil

header0 = "#ifndef PRECICE_NO_MPI\n";
header1 = "#include \"mpi.h\"\n";
header2 = "#endif\n";

filenames = os.listdir(".")
#filenames = [ glob("*.hpp"), glob("*.cpp")]
for filename in filenames:
   #print "Looking at", filename
   if os.path.isdir(filename):
      nestedfiles = os.listdir(filename + "/")
      for i in range(len(nestedfiles)):
         nestedfiles[i] = filename + "/" + nestedfiles[i]
      filenames.extend(nestedfiles)
   else:
      insert_cpp = filename[-4:] == ".cpp"
      #insert_cpp = insert_cpp or filename[-4:] ==".hpp" 
      #insert_cpp = insert_cpp or filename[-5:] == ".cpph"
      
      #insert_c = filename[-2:] == ".h" 
      #insert_c = insert_c or filename[-2:] == ".c"
      
      if insert_cpp:
         file = open(filename, 'r+')
         if file.readline() != header0:
            print "Insert mpi include in file \"" + filename + "\""
            file.seek(0)
            text = file.read()
            file.seek(0)
            file.write(header0)
            file.write(header1)
            file.write(header2)
            file.write(text)
         file.close()         
