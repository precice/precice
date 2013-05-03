#ifdef TEST
#define real double
#include "wave_profile_udf.h"
#include "stdio.h"
#else
#include "udf.h"

real get_ratio_triangle_is_covered (
  real coords[6],
  real wave_height );

real compute_triangle_area_2d(
  real coords[6]);

const real PI = 3.1415926535897932;
const real GRAVITY = 9.80665;
const real SHIFT = 3.1415926535897932*0.5;
int counter = 0;
real wave_amplitude = 0.0;
real omega = 0.0;
real k = 0.0;
real surface_level = 0.0;

/**
 * Returns the phase ratio (by water) of a boundary cell.
 */
real get_ratio_edge_is_covered (
  face_t  face,
  Thread* face_thread,
  real    wave_height );

DEFINE_INIT(init,domain)
{
  Message ( "\nEntering INIT\n" );
  Domain* subdomain;
  Thread* cell_thread; /* linked list of cells */
  cell_t cell;    /* one cell */
  int phase_domain_index;
  real cell_centroid[ND_ND]; /* n dimensional cell center coordinates */
  int flipzones = 0;
  int phase_sign = 0;
  real upper_phase = 0.0;
  real lower_phase = 0.0;

  int i;
  real coords[6];
  int node_index;
  Node* node;
  real phase;

# if !RP_NODE /* Serial or Host */
  surface_level = RP_Get_Float("udf/surfacelevel");
  flipzones = RP_Get_Integer("udf/flipzones");
  wave_amplitude = 0.5 * RP_Get_Float("udf/waveheight");
  omega = (2.0 * PI) / RP_Get_Float("udf/period");
  k = (2.0 * PI) / RP_Get_Float("udf/wavelength");
# endif
  host_to_node_real_4(surface_level, wave_amplitude, omega, k);
  host_to_node_int_1(flipzones);

  Message("surface level = %f\n", surface_level);
  Message("flipzones = %d\n", flipzones);
  Message("wave_amplitude = %f\n", wave_amplitude);
  Message("omega = %f\n", omega);
  Message("k = %f\n", k);
  
# if !RP_HOST /* serial or node */
  sub_domain_loop (subdomain, domain, phase_domain_index){ 
    if (phase_domain_index == 0){ /* primary phase */
      if (! flipzones){
        phase_sign = 1;
        upper_phase = 1.0;
        lower_phase = 0.0;
      }
      else {
        phase_sign = -1;
        upper_phase = 0.0;
        lower_phase = 1.0;
      }
    }
    else if (phase_domain_index == 1){
      if (! flipzones){
        phase_sign = -1;
        upper_phase = 0.0;
        lower_phase = 1.0;
      }
      else {
        phase_sign = 1;
        upper_phase = 1.0;
        lower_phase = 0.0;
      }
    }
    Message("Setting phase for upper=%f, lower=%f\n", upper_phase, lower_phase);
    thread_loop_c (cell_thread, subdomain){
      begin_c_loop_all (cell, cell_thread){
#       if ND_ND == 3
        C_CENTROID(cell_centroid, cell, cell_thread);
        if (cell_centroid[1] < surface_level){
          C_VOF(cell, cell_thread) = lower_phase;
        }
        else {
          C_VOF(cell, cell_thread) = upper_phase;
        }
#       endif /* ND_ND == 3 */
#       if ND_ND == 2
        i = 0;
        c_node_loop(cell, cell_thread, node_index){
          node = C_NODE(cell, cell_thread, node_index);
          coords[i++] = NODE_X(node);
          coords[i++] = NODE_Y(node);
        }
        /*if (i != 6){
          Message("Error: found cell with size nodes != 3\n");
          exit(1);
        }*/
        phase = phase_sign == 1
            ? 1.0 - get_ratio_triangle_is_covered(coords, surface_level)
            : get_ratio_triangle_is_covered(coords, surface_level);
        /*if (phase < 0.0 || phase > 1.0){
          Message("Error: phase = %f\n", phase);
        }*/

        C_CENTROID(cell_centroid, cell, cell_thread);
        /*if (cell_centroid[1] > 0.5 && cell_centroid[1] < 0.55){
          Message("At %f,%f: %f (coords = %f, %f, %f, %f, %f, %f)\n", cell_centroid[0], cell_centroid[1], phase,
                  coords[0], coords[1], coords[2], coords[3], coords[4], coords[5]);
        }*/
        C_VOF(cell, cell_thread) = phase;
#       endif /* ND_ND == 2 */
      } end_c_loop_all(cell_type, cell_thread)
    }
  }
# endif /* serial or node */
  Message ( "Leaving INIT\n" );
}

DEFINE_PROFILE(MP_wave_x_velocity,face_thread,var_index)
{
  /*Message ( "\nEntering DEFINE_PROFILE(wave_x_velocity)\n");*/
  real time = CURRENT_TIME;
  real centroid_coords[ND_ND];
  face_t face;
  real height = 0.0;
  real vel = 0.0;
  real wave_height = wave_amplitude * cos(-1.0 * omega * time - SHIFT);
  real vel_amplitude = (wave_amplitude * GRAVITY * k * cos(-1.0 * omega * time - SHIFT)) / (cosh(k * surface_level) * omega);

# if !RP_HOST /* Serial or node */
  begin_f_loop(face, face_thread){
    if (PRINCIPAL_FACE_P(face, face_thread)){
      F_CENTROID(centroid_coords, face, face_thread);
      height = centroid_coords[1] - surface_level;
      if (height <= wave_height){
        vel = vel_amplitude * cosh(k * (height + surface_level));
      }
      else {
        vel = vel_amplitude * cosh(k * (wave_height + surface_level));
      }
      /*Message("Setting x vel = %f at face centroid_y = %f\n", vel, height);*/
      F_PROFILE(face, face_thread, var_index) = vel;
    }
  } end_f_loop(face, face_thread)
# endif
  /*Message ( "Leaving DEFINE_PROFILE(wave_x_velocity)\n" );*/
}

DEFINE_PROFILE(MP_wave_y_velocity, face_thread, var_index)
{
  /*Message ( "\nEntering DEFINE_PROFILE(wave_y_velocity)\n" );*/
  real time = CURRENT_TIME;
  real centroid_coords[ND_ND];
  face_t face;
  real height = 0.0;
  real vel = 0.0;
  real wave_height = wave_amplitude * cos(-1.0 * omega * time - SHIFT);
  real vel_amplitude = (wave_amplitude * GRAVITY * k * sin(-1.0 * omega * time - SHIFT)) / (cosh(k * surface_level) * omega);

# if !RP_HOST /* Serial or node */
  begin_f_loop(face, face_thread){
    if (PRINCIPAL_FACE_P(face, face_thread)){
      F_CENTROID(centroid_coords, face, face_thread);
      height = centroid_coords[1] - surface_level;
      if (height <= wave_height){
        vel = vel_amplitude * sinh(k * (height + surface_level));
        /*vel = 0.2 * sin(-1.0 * omega * time  - SHIFT);*/
      }
      else {
        vel = vel_amplitude * sinh(k * (wave_height + surface_level));
      }
      /*Message("Setting y vel = %f at face centroid_y = %f\n", vel, height);*/
      F_PROFILE(face, face_thread, var_index) = vel;
    }
  } end_f_loop(face, face_thread)
# endif
  /*Message ( "Leaving DEFINE_PROFILE(wave_y_velocity)\n" );*/
}

DEFINE_PROFILE(MP_wave_volume_fraction, face_thread, var_index)
{
  /*Message ( "\nEntering DEFINE_PROFILE(wave_volume_fraction)\n" );*/
  real time = CURRENT_TIME;
  face_t face;
  real water_fraction = 0.0;
  real wave_height = 0.0;

  real coords[6];
  int i;
  int node_index;
  Node *node;

  wave_height = wave_amplitude * cos(-1.0 * omega * time - SHIFT);

# if !RP_HOST /* Serial or node */
  if (ND_ND == 2){
    begin_f_loop(face, face_thread){
      if (PRINCIPAL_FACE_P(face, face_thread)){
        water_fraction = get_ratio_edge_is_covered(face, face_thread, wave_height);
        /*Message("Setting water fraction = %f at face centroid_y = %f\n", water_fraction, centroid_y);*/
        F_PROFILE(face, face_thread, var_index) = water_fraction;
      }
    } end_f_loop(face, face_thread)
  }
  else { /* ND_ND == 3 */
    begin_f_loop(face, face_thread){
      if (PRINCIPAL_FACE_P(face, face_thread)){
        i = 0;
        f_node_loop(face, face_thread, node_index){
          node = F_NODE(face, face_thread, node_index);
          coords[i++] = NODE_Z(node);
          coords[i++] = NODE_Y(node);
        }
        water_fraction = get_ratio_triangle_is_covered(coords, wave_height);
        /*Message("Setting water fraction = %f at face centroid_y = %f\n", water_fraction, centroid_y);*/
        F_PROFILE(face, face_thread, var_index) = water_fraction;
      }
    } end_f_loop(face, face_thread)
  }
# endif
  /*Message ( "Leaving DEFINE_PROFILE(wave_volume_fraction)\n" );*/
}

real get_ratio_edge_is_covered (
  face_t  face,
  Thread* face_thread,
  real    wave_height )
{
  real centroid_coords[ND_ND];
  real area_vector[ND_ND];
  real area = 0.0;
  real centroid_y = 0.0;
  real upper_y = 0.0;
  real lower_y = 0.0;

  F_CENTROID(centroid_coords, face, face_thread);
  F_AREA(area_vector, face, face_thread);
  area = NV_MAG(area_vector);
  centroid_y = centroid_coords[1];
  upper_y = centroid_y + area/2.0;
  lower_y = centroid_y - area/2.0;
  if (lower_y > surface_level + wave_height){
     return 0.0;
  }
  else if (upper_y < surface_level + wave_height){
     return 1.0;
  }
  return (surface_level + wave_height - lower_y) / area;
}

#endif /* not TEST */

real get_ratio_triangle_is_covered (
  real coords[6],
  real wave_height )
{
  real ratio = 0.0;
  real min = coords[1];
  real max = coords[1];
  int above = 0;
  int below = 0;
  int above_count = 0;
  int i, j;
  real coords_leg_vertices[4];
  real lambda;
  real x_intersected_1;
  real x_intersected_2;
  real coords_small[6];
  real area, area_small;

  for (i=3; i < 6; i+=2){
    if (min > coords[i]) min = coords[i];
    if (max < coords[i]) max = coords[i];
  }

  if (min >= wave_height) return 0.0;
  if (max <= wave_height) return 1.0;

  /* Determine, whether one or two vertices are above interface */
  for (i=0; i < 3; i++){
    if (coords[i*2 +1] > wave_height){
      above_count++;
      above = i;
    }
    else {
      below = i;
    }
  }
  /*printf("above_count = %d, above = %d, below = %d\n", above_count, above, below);*/

  j=0;
  if (above_count == 1){ /* Select leg vertices below interface */
    coords_small[0] = coords[above*2];
    coords_small[1] = coords[above*2+1];
    for (i=0; i < 3; i++){
      if (i != above){
        coords_leg_vertices[j++] = coords[i*2];
        coords_leg_vertices[j++] = coords[i*2+1];
      }
    }
  }
  else { /* Select leg vertices above interface */
    coords_small[0] = coords[below*2];
    coords_small[1] = coords[below*2+1];
    for (i=0, j=0; i < 3; i++){
      if (i != below){
        coords_leg_vertices[j++] = coords[i*2];
        coords_leg_vertices[j++] = coords[i*2+1];
      }
    }
  }

  /*printf("coords head: %f, %f\n", coords_small[0], coords_small[1]);
  printf("leg vertex 1: %f, %f\n", coords_leg_vertices[0], coords_leg_vertices[1]);
  printf("leg vertex 2: %f, %f\n", coords_leg_vertices[2], coords_leg_vertices[3]);*/

  /* compute intersection points of interface with legs */
  lambda = (wave_height - coords_small[1]) / (coords_leg_vertices[1] - coords_small[1]);
  coords_small[2] = coords_small[0] + lambda*(coords_leg_vertices[0] - coords_small[0]);
  /*printf("lambda leg 1: %f, x leg 1: %f\n", lambda, coords_small[2]);*/
  lambda = (wave_height - coords_small[1]) / (coords_leg_vertices[3] - coords_small[1]);
  coords_small[4] = coords_small[0] + lambda*(coords_leg_vertices[2] - coords_small[0]);
  /*printf("lambda leg 2: %f, x leg 2: %f\n", lambda, coords_small[4]);*/
  coords_small[3] = wave_height;
  coords_small[5] = wave_height;

  area = compute_triangle_area_2d(coords);
  area_small = compute_triangle_area_2d(coords_small);
  /*printf("area = %f, area_small = %f\n", area, area_small);*/
  if (above_count == 1){
    ratio = (area - area_small) / area;
  }
  else {
    ratio = area_small / area;
  }
  if (ratio > 1.0 || ratio < 0.0){
    printf("Wrong triangle covered ratio of %f!\n", ratio);
    printf("area = %f, area_small = %f, above_count = %d, above = %d\n", area, area_small, above_count, above);
    printf("coords = (%f, %f) (%f, %f) (%f, %f); coords small = (%f, %f) (%f, %f) (%f, %f)\n",
           coords[0], coords[1], coords[2], coords[3], coords[4], coords[5],
           coords_small[0], coords_small[1], coords_small[2], coords_small[3], coords_small[4], coords_small[5]);
    exit(1);
  }
  return ratio;
}

real compute_triangle_area_2d(real coords[6])
{
  real ab[] = {coords[2] - coords[0], coords[3] - coords[1]};
  real ac[] = {coords[4] - coords[0], coords[5] - coords[1]};
  return fabs(0.5 * (ab[0]*ac[1] - ab[1]*ac[0]));
}
