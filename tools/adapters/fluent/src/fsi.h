#ifndef FSI_H_
#define FSI_H_

#include "udf.h"
#include "dynamesh_tools.h"

void fsi_init(Domain* domain);
void fsi_write_and_advance();
void fsi_grid_motion(Domain* domain, Dynamic_Thread* dt, real time, real dtime);
void fsi_plot_coords();

#endif /* FSI_H_ */
