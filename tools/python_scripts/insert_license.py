import os
import glob
import shutil

cpp_license0 = "// Copyright (C) 2011 Technische Universitaet Muenchen\n"
cpp_license1 = "// This file is part of the preCICE project. For conditions of distribution and\n"
cpp_license2 = "// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License\n"

c_license0 = "/* Copyright (C) 2011 Technische Universitaet Muenchen\n"
c_license1 = " * This file is part of the preCICE project. For conditions of distribution and\n"
c_license2 = " * use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License */\n"

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
      insert_cpp = insert_cpp or filename[-4:] ==".hpp" 
      insert_cpp = insert_cpp or filename[-5:] == ".cpph"
      
      insert_c = filename[-2:] == ".h" 
      insert_c = insert_c or filename[-2:] == ".c"
      
      if insert_cpp:
         file = open(filename, 'r+')
         if file.readline() != cpp_license0:
            print "Insert cpp license in file \"" + filename + "\""
            file.seek(0)
            text = file.read()
            file.seek(0)
            file.write(cpp_license0)
            file.write(cpp_license1)
            file.write(cpp_license2)
            file.write(text)
         file.close()
         
      elif insert_c:
         file = open(filename, 'r+')
         if file.readline() != c_license0:
            print "Insert c license in file \"" + filename + "\""
            file.seek(0)
            text = file.read()
            file.seek(0)
            file.write(c_license0)
            file.write(c_license1)
            file.write(c_license2)
            file.write(text)
         file.close()
