#include "fsi.h"

DEFINE_INIT(init,domain)
{
  fsi_init(domain);
}

DEFINE_ON_DEMAND(write_and_advance)
{
  fsi_write_and_advance();
}

DEFINE_GRID_MOTION(gridmotions,domain,dt,time,dtime)
{
  fsi_grid_motion(domain, dt, time, dtime);
}

DEFINE_ON_DEMAND(plot_coords)
{
  fsi_plot_coords();
}
