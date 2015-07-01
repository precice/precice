-- Configuration file for Ateles --


-- This is a configuration file for the Finite Volume / Discontinuous Galerkin Solver ATELES.
-- It provides a testcase for the simulation of Euler equations in a homogenous media. The simulation domain
-- is a periodic cube with edge length 2.0. Therefore this is a very good way to verify your algorithmic implementations,
-- since this testcase does not involve any boundary conditions.
-- The testcase simulates the temporal development of Gaussian pulse in density. Since we
-- are considering a very simple domain, an analytic solution is well known and given as Lua functions in this script.
-- Therefore we suggest this testcase to verify one of the following topics
-- ... algorihtmic correctness
-- ... spatial and temporal order of your scheme
-- ... diffusion and dispersion relations of your scheme
-- ... and many more.
-- This testcase can be run in serial (only one execution unit) or in parallel (with multiple mpi ranks).
-- To specify the number of ranks please modify nprocs variable. To calculate a grid convergence behavior please modify the
-- level variable. An increment of one will half the radius of your elements.

timestep_info = 1

-- Check for Nans and unphysical values
check =  {
           interval = 1,
         }

cubeLength = 2.0
-- global simulation options
simulation_name='ateles_left'
fin_time = 0.008
sim_control = {
             time_control = {
                  min = 0,
                  max = {sim=fin_time}, -- final simulation time
                  interval = {iter=1}
                }
}

-- table for preCICE
precice = {
           accessor = 'Ateles_left',
           configFile ='../precice-config.xml',
          }

-- Mesh definitions --
mesh = 'mesh/16384/'

---- Restart settings
estart = {
--            -- file to restart from
--            read = './restart/twoway/left/simulation_lastHeader.lua',
--            -- folder to write restart data to
            write = './restart/',
            -- temporal definition of restart write
            time_control = { min = 0, max = fin_time, interval = {iter=1} }
          }

-- timing settings (i.e. output for performance measurements, this table is otional)
timing = {
          folder = './',                  -- the folder for the timing results
          filename = 'timing_left.res'         -- the filename of the timing results
         }

-- Equation definitions --
equation = {
    name   = 'euler',
    therm_cond = 2.555e-02,
    isen_coef = 1.4,
    r      = 296.0,
}
-- (cv) heat capacity and (r) ideal gas constant
equation["cv"] = equation["r"] / (equation["isen_coef"] - 1.0)

-- Scheme definitions --
scheme = {
    -- the spatial discretization scheme
    spatial =  {
               name = 'modg',           -- we use the modal discontinuous Galerkin scheme
               modg_space = 'Q',        -- the polynomial space Q or P
               m = 32,                   -- the maximal polynomial degree for each spatial direction
               },
    -- the temporal discretization scheme
    temporal = {
               name = 'explicitRungeKutta',  --'explicitEuler',
               steps = 4,
               -- how to control the timestep
               control = {
                          name = 'cfl',   -- the name of the timestep control mechanism
                          cfl  = 0.8,     -- Courant–Friedrichs–Lewy number
                         },
               },
}

function dens(x,y,z)
  return x
end

projection = {
              kind = 'fpt',  -- 'fpt' or 'l2p', default 'l2p'
              -- for fpt the  nodes are automatically 'chebyshev'
              -- for lep the  nodes are automatically 'gauss-legendre'
           -- lobattoPoints = false  -- if lobatto points should be used, default = false
              factor = 1.0,          -- dealising factpr for fpt, oversampling factor for l2p, float, default 1.0
           -- blocksize = 32,        -- for fpt, default -1
           -- fftMultiThread = false -- for fpt, logical, default false
             }

-- This is a very simple example to define constant boundary condtions.
-- Transport velocity of the pulse in x direction.
velocityX = 250
initial_condition = {-- density = 3.0,
                             --  {
                             --    predefined='gausspulse',
                             --    center={1.8, 1.8, 1.8},
                             --    halfwidth=0.20,
                             --    amplitude=2.0,
                             --    background=1.225
                             --   },
                      density = {
                                    predefined='gausspulse',
                                    center={-1.0, 0.0, 0.0},
                                    halfwidth=0.2,
                                    amplitude=2.0,
                                    background=1.225
                                 },
                      pressure = 100000,
                      velocityX = velocityX,
                      velocityY = 0.0,
                      velocityZ = 0.0,
                    }


-- Tracking
racking = {
             label = 'track_momentum_left',
             folder = './',
             variable = {'energy'},
             shape = {kind = 'canoND', object= { origin ={-1.5, 0., 0.} } },
             time = {
                      useIterations = true,
                      min = 0,
                      max = sim_control.time_control.max,
                      interval = 1,
                    },
             format = 'ascii'
           }

 -- Boundary definitions
 boundary_condition = {
                         {
                           label = 'wall_1',
                           kind = 'outflow',
                           pressure = 100000,
                           density = 'extrapolate'
                         }
                         ,
                         {
                           label = 'wall_2',
                           kind = 'outflow',
                           pressure = 100000,
                           density = 'extrapolate'
                         }
                         ,
                         {
                           label = 'wall_3',
                             kind = 'inflow',
                             pressure = 'extrapolate',
                             density = 1.225,
                             v_x = velocityX,
                             v_y = 0.0,
                             v_z = 0.0
                         }
                         ,
                         {
                           label = 'wall_5',
                           kind = 'outflow',
                           pressure = 100000,
                           density = 'extrapolate'
                         }
                         ,
                         {
                           label = 'wall_6',
                           kind = 'outflow',
                           pressure = 100000,
                           density = 'extrapolate'
                         }
                         ,
                         {
                           label = 'precice_leftmesh',
                           kind = 'precice',
                           precice_mesh = 'AcousticSurface_left',
                           provide_mesh = true,
                           exchange_data_write =  {'Density_left',
                                                   'Momentum_X_left',
                                                   'Momentum_Y_left',
                                                   'Momentum_Z_left',
                                                   'Energy_left'
                                                  },
                           exchange_data_read =  {'Density_right',
                                                  'Momentum_X_right',
                                                  'Momentum_Y_right',
                                                  'Momentum_Z_right',
                                                  'Energy_right'
                                                  }
                         }
                       }

precice_write = {  label = {'precice_leftmesh'}
               }
-- DEBUG OPTIONS
debug = {
         debugMode = true,        -- default= false
         verbose = 99,             -- default= 0
         debugFiles = true,       -- default= false
         -- What to dump into debugFiles? --
           dumpTreeIDs = true,      -- default= false
           dumpPropBits = true,     -- default= false
           dumpAuxLists = true,     -- default= false
           dumpDependencies = true, -- default= false
           dumpState = true,        -- default= false
           dumpHaloState = true,    -- default= false
         --  end debugFiles  --
         debugDependencies = true, -- default= false
         checkDependencies = true, -- default= false
         checkNans = true,         -- default= false
         checkSteps = true,        -- default= false
         debugMesh = 'dbg/mesh_',  -- default= ''
         debugSource = true,       -- default= false
         debugRestart = true,      -- default= false
         traceMemory = true,       -- default= false
        }
