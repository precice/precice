#include "fsi.h"
#include "wave_profile.h"

DEFINE_INIT(init,domain)
{
  wave_profile_init(domain);
  fsi_init(domain);
}

DEFINE_PROFILE(MP_wave_x_velocity,face_thread,var_index)
{
  wave_profile_x_velocity(face_thread, var_index);
}

DEFINE_PROFILE(MP_wave_y_velocity, face_thread, var_index)
{
  wave_profile_y_velocity(face_thread, var_index);
}

DEFINE_PROFILE(MP_wave_volume_fraction, face_thread, var_index)
{
  wave_profile_volume_fraction(face_thread, var_index);
}

DEFINE_ON_DEMAND(write_and_advance)
{
  fsi_write_and_advance();
}

DEFINE_GRID_MOTION(gridmotions,domain,dt,time,dtime)
{
  fsi_grid_motion(domain, dt, time, dtime);
}
