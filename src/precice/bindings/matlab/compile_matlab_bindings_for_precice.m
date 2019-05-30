% Script for compiling the MATLAB bindings for preCICE

% Get the path to the binding files. Use mfilename to get the full path of
% the current file. Thus, the script will work even if the user calls it
% from a different folder
path = string(fileparts(mfilename('fullpath')));
path_Interface = strjoin([path,"+precice","@SolverInterface","private","preciceGateway"],filesep);
path_Constants = strjoin([path,"+precice","@Constants","private","preciceGateway"],filesep);

% Get the flags for linking to preCICE
[status,flags] = system('pkg-config --cflags --libs libprecice');
if status==1
    error("pkg-config was unable to determine the compiler flags for preCICE.")
end
flags = strsplit(flags);

% Run mex commands to compile
mex(strcat(path_Interface,".cpp"),"-output",path_Interface,flags{:});
mex(strcat(path_Constants,".cpp"),"-output",path_Constants,flags{:});

%mex ~/precice/src/precice/bindings/matlab/+precice/@SolverInterface/private/preciceGateway.cpp -output ...
%    ~/precice/src/precice/bindings/matlab/+precice/@SolverInterface/private/preciceGateway -lprecice;
%mex ~/precice/src/precice/bindings/matlab/+precice/@Constants/private/preciceConstants.cpp -output ...
%    ~/precice/src/precice/bindings/matlab/+precice/@Constants/private/preciceConstants -lprecice;