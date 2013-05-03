#include "wave_profile.h"
#include "udf.h"

DEFINE_INIT(init,domain)
{
  Message ( "\nEntering INIT\n" );
  wave_profile_init(domain);
  Message ( "Leaving INIT\n" );
}

DEFINE_PROFILE(MP_wave_x_velocity,face_thread,var_index)
{
  wave_profile_x_velocity(face_thread, var_index);
}

DEFINE_PROFILE(MP_wave_y_velocity, face_thread, var_index)
{  wave_profile_y_velocity(face_thread, var_index);
}

DEFINE_PROFILE(MP_wave_volume_fraction, face_thread, var_index)
{
  wave_profile_volume_fraction(face_thread, var_index);
}
