#ifndef WAVE_PROFILE_H_
#define WAVE_PROFILE_H_

#include "udf.h"

void wave_profile_init(Domain* domain);

void wave_profile_x_velocity(Thread* face_thread, int var_index);
void wave_profile_y_velocity(Thread* face_thread, int var_index);
void wave_profile_volume_fraction(Thread* face_thread, int var_index);


#endif /* WAVE_PROFILE_H_ */
