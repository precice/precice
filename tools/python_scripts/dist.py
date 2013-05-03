import os
import tarfile
import shutil

# Howto use this script:
#
# 1) Adapt the name of the archive to be built (typically the version number)
# 2) Check the files to be added
# 3) Go to the root folder of the preCICE project and execute it from there

##### Determine path to peano root from environment variable
#
tarchSrc = os.getenv ('TARCH_SRC')
if ( tarchSrc == None ):
   print 'ERROR: Environment variable TARCH_SRC not defined!'
   sys.exit(1)
else:
   print 'Using environment variable TARCH_SRC =', tarchSrc

shutil.copytree(tarchSrc + "/tarch", "src/tarch")

name = "preCICE_v0_8_1.tar.gz"
print "Creating archive \"" + name + "\" ..."
tar = tarfile.open(name, "w:gz")
addCount = 0
excludeCount = 0
filenames = ["src", "tools/livegraph", "tools/solverproxies", "SConstruct", "changelog.txt",
             "test-config.xml", "integration-test-config.xml"]
for filename in filenames:
   #print "File", filename
   lastPart = filename[filename.rfind("/")+1:]
   #print "lastPart:", lastPart
   if lastPart[:1] == ".":
      #print "Exclude \"" + filename + "\""
      excludeCount += 1
   elif os.path.isdir(filename):
      nestedfiles = os.listdir(filename + "/")
      for i in range(len(nestedfiles)):
         nestedfiles[i] = filename + "/" + nestedfiles[i]
      filenames.extend(nestedfiles)
   else:
      #print "Adding file \"" + filename + "\""
      addCount += 1
      tar.add(filename)
tar.close()
shutil.rmtree("src/tarch")
print "Done. Added " + str(addCount) + ", excluded " + str(excludeCount)
