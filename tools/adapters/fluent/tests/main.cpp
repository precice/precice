#include <iostream>
#include <cmath>
extern "C" {
#include "wave_profile_udf.h"
}

#define print(message) \
  std::cout << message << std::endl

#define validate(statement) \
  test_count++; \
  if (statement == false){ \
    print("Validate fails for statement: " << #statement << " (line " << __LINE__ << ")"  ); \
    errors++; \
  }

bool equals(double a, double b, double eps = 1e-14)
{
  return std::abs(a - b) < eps;
}

int test_count = 0;
int errors = 0;

int main(int argc, char** argv)
{
  print("Running fluent tests");

  double coords[6];
  double* a = &coords[0];
  double* b = &coords[2];
  double* c = &coords[4];

  a[0] = 0.0; a[1] = 0.0;
  b[0] = 1.0; b[1] = 0.0;
  c[0] = 1.0; c[1] = 1.0;
  double area = compute_triangle_area_2d(coords);
  validate(equals(area, 0.5));
  c[0] = 0.5; c[1] = 1.0;
  area = compute_triangle_area_2d(coords);
  validate(equals(area, 0.5));

  double height = 2.0;
  double ratio = get_ratio_triangle_is_covered(coords, height);
  validate(equals(ratio, 1.0));
  height = 1.0;
  ratio = get_ratio_triangle_is_covered(coords, height);
  validate(equals(ratio, 1.0));
  height = -1.0;
  ratio = get_ratio_triangle_is_covered(coords, height);
  validate(equals(ratio, 0.0));
  height = 0.0;
  ratio = get_ratio_triangle_is_covered(coords, height);
  validate(equals(ratio, 0.0));
  height = 0.5;
  ratio = get_ratio_triangle_is_covered(coords, height);
  validate(equals(ratio, 3.0/4.0));
  height = 0.25;
  ratio = get_ratio_triangle_is_covered(coords, height);
  validate(equals(ratio, 7.0/16.0));
  a[0] = 0.0; a[1] = -0.5;
  b[0] = 1.0; b[1] = 0.0;
  c[0] = 0.0; c[1] = 0.5;
  height = 0.0;
  ratio = get_ratio_triangle_is_covered(coords, height);
  validate(equals(ratio, 0.5));
  height = -0.25;
  ratio = get_ratio_triangle_is_covered(coords, height);
  validate(equals(ratio, 1.0/8.0));
  height = 0.25;
  ratio = get_ratio_triangle_is_covered(coords, height);
  validate(equals(ratio, 7.0/8.0));

  print("Finished running " << test_count << " tests with " << errors << " errors");
  return 0;
}

