% MATLAB solverdummy.
% To use it, don't forget to install the matlab bindings and add them to
% the path
function solverdummy(configFileName,solverName,meshName)
    if nargin~=3
        disp('Usage: solverdummy configFile solverName meshName');
        disp('');
        disp('Parameter description');
        disp('  configurationFile: Path and filename of preCICE configuration');
        disp('  solverName:        SolverDummy participant name in preCICE configuration');
        disp('  meshName:          Mesh in preCICE configuration that carries read and write data');
        return;
    end
    
    interface = precice.SolverInterfaceOOP(solverName);
    interface.configure(configFileName);
    
    meshID = interface.getMeshID(meshName);
    dims = interface.getDimensions();
    
    dataIndices = interface.setMeshVertex(meshID,zeros(dims,1));
    
    dt = interface.initialize();
    
    while(interface.isCouplingOngoing())
        if(interface.isActionRequired('write-iteration-checkpoint'))
            disp('DUMMY: Writing iteration checkpoint.')
            interface.fulfilledAction('write-iteration-checkpoint');
        end
        
        dt = interface.advance(dt);
        
        if(interface.isActionRequired('read-iteration-checkpoint'))
            disp('DUMMY: Reading iteration checkpoint.')
            interface.fulfilledAction('read-iteration-checkpoint');
        else
            disp('DUMMY: Advaning in time.')
        end
    end
    
    interface.finalize();
    disp('DUMMY: Closing MATLAB solverdummy.')
end