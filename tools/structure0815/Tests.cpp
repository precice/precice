#include "Tests.hpp"
#include "Globals.hpp"
#include "Structure0815.hpp"

#define validate(booleanExpr, msg) if (!(booleanExpr)) { \
    std::cerr << "  boolean test failed " << std::endl \
              << "  file: " << __FILE__ << " \t line: " << __LINE__ << std::endl \
              << "  statement: " << #booleanExpr \
              << "  msg: " << msg \
              << std::endl; \
    return;\
  }

void Tests:: run()
{
  STRUCTURE_DEBUG("Running tests...");
  test2D();
  test3D();
  STRUCTURE_DEBUG("...done running tests.");
}

void Tests:: test2D()
{
  STRUCTURE_DEBUG("  test2D()...");
  using precice::utils::DynVector;
  using precice::utils::Vector2D;
  int dim = 2;
  double density = 1;
  DynVector gravity(2, 0.0);
  DynVector vertices(4*dim,0.0);
  tarch::la::DynamicVector<int> faces(4*dim,0);

  vertices = 0.0, 0.0,
             1.0, 0.0,
             1.0, 1.0,
             0.0, 1.0;
  faces = 0, 1, 1, 2, 2, 3, 3, 0;

  Structure0815 structure(dim, density, gravity, vertices, faces);

  DynVector expected(4*dim, 0.0);
  DynVector& forces = structure.forces();
  const DynVector& displacements = structure.displacements();
  const DynVector& displDeltas = structure.displacementDeltas();
  const DynVector& velocities = structure.velocities();
  const DynVector& velDeltas = structure.velocityDeltas();

  validate(forces.size() == 8, forces.size());
  validate(displacements.size() == 8, displacements.size());
  validate(displDeltas.size() == 8, displDeltas.size());
  validate(velocities.size() == 8, velocities.size());
  validate(velDeltas.size() == 8, velDeltas.size());

  assign(expected) = 0.0;
  validate(tarch::la::equals(forces, expected), forces);
  validate(tarch::la::equals(displacements, expected), displacements);
  validate(tarch::la::equals(velocities, expected), velocities);
  validate(tarch::la::equals(displDeltas, expected), displDeltas);
  validate(tarch::la::equals(velDeltas, expected), velDeltas);

  validate(structure.getCenterOfGravity().size() == 2, structure.getCenterOfGravity().size());
  validate(tarch::la::equals(structure.getCenterOfGravity(), Vector2D(0.5, 0.5)),
           structure.getCenterOfGravity());

  double dt = 0.0;
  structure.iterate(dt);

  validate(tarch::la::equals(forces, expected), forces);
  validate(tarch::la::equals(displacements, expected), displacements);
  validate(tarch::la::equals(velocities, expected), velocities);
  validate(tarch::la::equals(displDeltas, expected), displDeltas);
  validate(tarch::la::equals(velDeltas, expected), velDeltas);

  forces = 1.0, 0.0,
           1.0, 0.0,
           1.0, 0.0,
           1.0, 0.0;

  dt = 1.0;
  structure.iterate(dt);

  expected = forces;
  validate(tarch::la::equals(forces, expected), forces);
  expected = 4.0, 0.0,
             4.0, 0.0,
             4.0, 0.0,
             4.0, 0.0;
  validate(tarch::la::equals(displacements, expected), displacements);
  validate(tarch::la::equals(velocities, expected), velocities);
  validate(tarch::la::equals(displDeltas, expected), displDeltas);
  validate(tarch::la::equals(velDeltas, expected), velDeltas);

  structure.timestep(dt);

  validate(tarch::la::equals(displacements, expected), displacements);
  validate(tarch::la::equals(velocities, expected), velocities);
  validate(tarch::la::equals(displDeltas, expected), displDeltas);
  validate(tarch::la::equals(velDeltas, expected), velDeltas);

  forces = 0.0, 1.0,
           0.0, 1.0,
           0.0, 1.0,
           0.0, 1.0;

  structure.iterate(dt);
  expected = 8.0, 4.0,
             8.0, 4.0,
             8.0, 4.0,
             8.0, 4.0;
  validate(tarch::la::equals(displacements, expected), displacements);
  expected = 4.0, 4.0,
             4.0, 4.0,
             4.0, 4.0,
             4.0, 4.0;
  validate(tarch::la::equals(velocities, expected), velocities);
  expected = 4.0, 4.0,
             4.0, 4.0,
             4.0, 4.0,
             4.0, 4.0;
  validate(tarch::la::equals(displDeltas, expected), displDeltas);
  expected = 0.0, 4.0,
             0.0, 4.0,
             0.0, 4.0,
             0.0, 4.0;
  validate(tarch::la::equals(velDeltas, expected), velDeltas);

  forces = -1.0, 0.0,
           -1.0, 0.0,
           -1.0, 0.0,
           -1.0, 0.0;

  structure.iterate(dt);
  expected = 4.0, 0.0,
             4.0, 0.0,
             4.0, 0.0,
             4.0, 0.0;
  validate(tarch::la::equals(displacements, expected), displacements);
  expected = 0.0, 0.0,
             0.0, 0.0,
             0.0, 0.0,
             0.0, 0.0;
  validate(tarch::la::equals(velocities, expected), velocities);
  expected = 0.0, 0.0,
             0.0, 0.0,
             0.0, 0.0,
             0.0, 0.0;
  validate(tarch::la::equals(displDeltas, expected), displDeltas);
  expected = -4.0, 0.0,
             -4.0, 0.0,
             -4.0, 0.0,
             -4.0, 0.0;
  validate(tarch::la::equals(velDeltas, expected), velDeltas);

  structure.timestep(dt);
  structure.iterate(dt);

  expected = 0.0, 0.0,
             0.0, 0.0,
             0.0, 0.0,
             0.0, 0.0;
  validate(tarch::la::equals(displacements, expected), displacements);
  expected = -4.0, 0.0,
             -4.0, 0.0,
             -4.0, 0.0,
             -4.0, 0.0;
  validate(tarch::la::equals(velocities, expected), velocities);
  expected = -4.0, 0.0,
             -4.0, 0.0,
             -4.0, 0.0,
             -4.0, 0.0;
  validate(tarch::la::equals(displDeltas, expected), displDeltas);
  expected = -4.0, 0.0,
             -4.0, 0.0,
             -4.0, 0.0,
             -4.0, 0.0;
  validate(tarch::la::equals(velDeltas, expected), velDeltas);

  structure.timestep(dt);

  forces = 1.0, 0.0,
           1.0, 0.0,
           1.0, 0.0,
           1.0, 0.0;

  structure.iterate(dt);

  expected = 0.0, 0.0,
             0.0, 0.0,
             0.0, 0.0,
             0.0, 0.0;
  validate(tarch::la::equals(displacements, expected), displacements);
  expected = 0.0, 0.0,
             0.0, 0.0,
             0.0, 0.0,
             0.0, 0.0;
  validate(tarch::la::equals(velocities, expected), velocities);
  expected = 0.0, 0.0,
             0.0, 0.0,
             0.0, 0.0,
             0.0, 0.0;
  validate(tarch::la::equals(displDeltas, expected), displDeltas);
  expected = 4.0, 0.0,
             4.0, 0.0,
             4.0, 0.0,
             4.0, 0.0;
  validate(tarch::la::equals(velDeltas, expected), velDeltas);

  structure.timestep(dt);

  forces = 1.0, 0.0,
           0.0, 0.0,
           -1.0, 0.0,
           0.0, 0.0;

  structure.iterate(dt);

  double rot = 0.5;
  expected = rot, -rot,
             rot, rot,
             -rot, rot,
             -rot, -rot;
  validate(tarch::la::equals(displacements, expected), displacements);

  STRUCTURE_DEBUG("  ...done test2D()");
}

void Tests:: test3D()
{
  STRUCTURE_DEBUG("  test3D()...");
  using precice::utils::DynVector;
  using precice::utils::Vector3D;
  int dim = 3;
  double density = 1;
  DynVector gravity(3, 0.0);
  DynVector vertices(8*dim,0.0); // Cube vertex coordinates
  tarch::la::DynamicVector<int> faces(6*2*dim,0); // Triangle vertex indices
  double eps = 1e-13;

  vertices = 0.0, 0.0, 0.0, // 0
             1.0, 0.0, 0.0, // 1
             1.0, 1.0, 0.0, // 2
             0.0, 1.0, 0.0, // 3
             0.0, 0.0, 1.0, // 4
             1.0, 0.0, 1.0, // 5
             1.0, 1.0, 1.0, // 6
             0.0, 1.0, 1.0; // 7
  faces = 0, 1, 2, // z=0
          0, 2, 3,
          4, 5, 6, // z=1
          4, 6, 7,
          0, 3, 7, // x=0
          0, 7, 4,
          1, 2, 6, // x=1
          1, 6, 5,
          0, 5, 4, // y=0
          0, 1, 5,
          3, 6, 7, // y=1
          3, 2, 6;

  Structure0815 structure(dim, density, gravity, vertices, faces);

  DynVector expected(8*dim, 0.0);
  DynVector& forces = structure.forces();
  const DynVector& displacements = structure.displacements();
  const DynVector& displDeltas = structure.displacementDeltas();
  const DynVector& velocities = structure.velocities();
  const DynVector& velDeltas = structure.velocityDeltas();

  validate(forces.size() == 24, forces.size());
  validate(displacements.size() == 24, displacements.size());
  validate(displDeltas.size() == 24, displDeltas.size());
  validate(velocities.size() == 24, velocities.size());
  validate(velDeltas.size() == 24, velDeltas.size());

  assign(expected) = 0.0;
  validate(tarch::la::equals(forces, expected), forces);
  validate(tarch::la::equals(displacements, expected), displacements);
  validate(tarch::la::equals(velocities, expected), velocities);
  validate(tarch::la::equals(displDeltas, expected), displDeltas);
  validate(tarch::la::equals(velDeltas, expected), velDeltas);

  validate(structure.getCenterOfGravity().size() == 3, structure.getCenterOfGravity().size());
  validate(tarch::la::equals(structure.getCenterOfGravity(), Vector3D(0.5, 0.5, 0.5)),
           structure.getCenterOfGravity());

  double dt = 0.0;
  structure.iterate(dt);

  validate(tarch::la::equals(forces, expected), forces);
  validate(tarch::la::equals(displacements, expected), displacements);
  validate(tarch::la::equals(velocities, expected), velocities);
  validate(tarch::la::equals(displDeltas, expected), displDeltas);
  validate(tarch::la::equals(velDeltas, expected), velDeltas);

  forces = 0.0, 0.0, 1.0,
           0.0, 0.0, 1.0,
           0.0, 0.0, 1.0,
           0.0, 0.0, 1.0,
           0.0, 0.0, 1.0,
           0.0, 0.0, 1.0,
           0.0, 0.0, 1.0,
           0.0, 0.0, 1.0;
  dt = 2.0;
  structure.iterate(dt);

  expected = 0.0, 0.0, 32.0,
             0.0, 0.0, 32.0,
             0.0, 0.0, 32.0,
             0.0, 0.0, 32.0,
             0.0, 0.0, 32.0,
             0.0, 0.0, 32.0,
             0.0, 0.0, 32.0,
             0.0, 0.0, 32.0;
  validate(tarch::la::equals(displacements, expected, eps), displacements);
  validate(tarch::la::equals(displDeltas, expected, eps), displDeltas);
  expected = 0.0, 0.0, 16.0,
             0.0, 0.0, 16.0,
             0.0, 0.0, 16.0,
             0.0, 0.0, 16.0,
             0.0, 0.0, 16.0,
             0.0, 0.0, 16.0,
             0.0, 0.0, 16.0,
             0.0, 0.0, 16.0;
  validate(tarch::la::equals(velocities, expected, eps), velocities);
  validate(tarch::la::equals(velDeltas, expected, eps), velDeltas);

  forces = 0.0, 1.0, 0.0,
           0.0, 1.0, 0.0,
           0.0, 1.0, 0.0,
           0.0, 1.0, 0.0,
           0.0, 1.0, 0.0,
           0.0, 1.0, 0.0,
           0.0, 1.0, 0.0,
           0.0, 1.0, 0.0;
  dt = 1.0;
  structure.iterate(dt);

  expected = 0.0, 8.0, 0.0,
             0.0, 8.0, 0.0,
             0.0, 8.0, 0.0,
             0.0, 8.0, 0.0,
             0.0, 8.0, 0.0,
             0.0, 8.0, 0.0,
             0.0, 8.0, 0.0,
             0.0, 8.0, 0.0;
  validate(tarch::la::equals(displacements, expected, eps), displacements);
  validate(tarch::la::equals(displDeltas, expected, eps), displDeltas);
  validate(tarch::la::equals(velocities, expected, eps), velocities);
  validate(tarch::la::equals(velDeltas, expected, eps), velDeltas);

  forces = 1.0, 0.0, 0.0,
           1.0, 0.0, 0.0,
           1.0, 0.0, 0.0,
           1.0, 0.0, 0.0,
           1.0, 0.0, 0.0,
           1.0, 0.0, 0.0,
           1.0, 0.0, 0.0,
           1.0, 0.0, 0.0;
  dt = 1.0;
  structure.iterate(dt);

  expected = 8.0, 0.0, 0.0,
             8.0, 0.0, 0.0,
             8.0, 0.0, 0.0,
             8.0, 0.0, 0.0,
             8.0, 0.0, 0.0,
             8.0, 0.0, 0.0,
             8.0, 0.0, 0.0,
             8.0, 0.0, 0.0;
  validate(tarch::la::equals(displacements, expected, eps), displacements);
  validate(tarch::la::equals(displDeltas, expected, eps), displDeltas);
  validate(tarch::la::equals(velocities, expected, eps), velocities);
  validate(tarch::la::equals(velDeltas, expected, eps), velDeltas);

  forces = 1.0, 0.0, 0.0,
           1.0, 0.0, 0.0,
           1.0, 0.0, 0.0,
           1.0, 0.0, 0.0,
           -1.0, 0.0, 0.0,
           -1.0, 0.0, 0.0,
           -1.0, 0.0, 0.0,
           -1.0, 0.0, 0.0;
  dt = 1.0;
  structure.iterate(dt);

  expected = 2.0, 0.0, -2.0,
             2.0, 0.0, 2.0,
             2.0, 0.0, 2.0,
             2.0, 0.0, -2.0,
             -2.0, 0.0, -2.0,
             -2.0, 0.0, 2.0,
             -2.0, 0.0, 2.0,
             -2.0, 0.0, -2.0;
  validate(tarch::la::equals(displacements, expected, eps), displacements);

  forces = 0.0, 1.0, 0.0,
           0.0, 1.0, 0.0,
           0.0, 1.0, 0.0,
           0.0, 1.0, 0.0,
           0.0, -1.0, 0.0,
           0.0, -1.0, 0.0,
           0.0, -1.0, 0.0,
           0.0, -1.0, 0.0;
  dt = 1.0;
  structure.iterate(dt);

  expected = 0.0, 2.0, -2.0,
             0.0, 2.0, -2.0,
             0.0, 2.0, 2.0,
             0.0, 2.0, 2.0,
             0.0, -2.0, -2.0,
             0.0, -2.0, -2.0,
             0.0, -2.0, 2.0,
             0.0, -2.0, 2.0;
  validate(tarch::la::equals(displacements, expected, eps), displacements);

  forces = 0.0, 0.0, 1.0,
           0.0, 0.0, -1.0,
           0.0, 0.0, -1.0,
           0.0, 0.0, 1.0,
           0.0, 0.0, 1.0,
           0.0, 0.0, -1.0,
           0.0, 0.0, -1.0,
           0.0, 0.0, 1.0;
  dt = 1.0;
  structure.iterate(dt);

  expected = -2.0, 0.0, 2.0,
             -2.0, 0.0, -2.0,
             -2.0, 0.0, -2.0,
             -2.0, 0.0, 2.0,
             2.0, 0.0, 2.0,
             2.0, 0.0, -2.0,
             2.0, 0.0, -2.0,
             2.0, 0.0, 2.0;
  validate(tarch::la::equals(displacements, expected, eps), displacements);

  STRUCTURE_DEBUG("  ...done test3D()");
}

