clear; close all; clc;

% Initialize and configure preCICE
interface = precice.SolverInterface("SolverU", "precice-config.xml", 0, 1);
cowid = precice.Constants.actionWriteInitialData(); % Required for data initialization
coric = precice.Constants.actionReadIterationCheckpoint(); % For implicit coupling
cowic = precice.Constants.actionWriteIterationCheckpoint(); % For implicit coupling

% Geometry IDs. As it is a 0-D simulation, only one vertex is necessary.
meshID = interface.getMeshID("MeshU");
vertex_ID = interface.setMeshVertex(meshID, [0 0]);

% Data IDs
I_ID = interface.getDataID("I", meshID);
U_ID = interface.getDataID("U", meshID);

% Simulation parameters and initial condition
C = 2;                      % Capacitance
L = 1;                      % Inductance
t0 = 0;                     % Initial simulation time
t_max = 10;                 % End simulation time
Io = 1;                     % Initial current
phi = 0;                    % Phase of the signal

w0 = 1/sqrt(L*C);           % Resonant frequency
I0 = Io*cos(phi);           % Initial condition for I
U0 = -w0*L*Io*sin(phi);     % Initial condition for U

f_U = @(t, U, I) -I/C;      % Time derivative of U

% Initialize simulation
dt = interface.initialize();
if (interface.isActionRequired(cowid))
    interface.writeScalarData(U_ID, vertex_ID, U0);
    interface.fulfilledAction(cowid)
end
interface.initializeData();

% Start simulation
t = t0 + dt;
while interface.isCouplingOngoing()

    % Record checkpoint if necessary
    if interface.isActionRequired(cowic)
        I0_checkpoint = I0;
        U0_checkpoint = U0;
        interface.fulfilledAction(cowic)
    end

    % Make Simulation Step
    [t_ode, U_ode] = ode45(@(t, y) f_U(t, y, I0), [t0 t], U0);
    U0 = U_ode(end);

    % Exchange data
    interface.writeScalarData(U_ID, vertex_ID, U0);
    dt = interface.advance(dt);

    % Recover checkpoint if not converged, else finish time step
    if interface.isActionRequired(coric)
        I0 = I0_checkpoint;
        U0 = U0_checkpoint;
        interface.fulfilledAction(coric)
    else
        I0 = interface.readScalarData(I_ID, vertex_ID);
        t0 = t;
        t = t0 + dt;
    end

end

% Stop coupling
interface.finalize();