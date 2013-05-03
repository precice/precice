#include "udf.h"

DEFINE_INIT(init,domain)
{
  Message ( "\nEntering INIT\n" );
  Domain* subdomain;
  Thread* cell_thread; /* linked list of cells */
  cell_t cell;    /* one cell */
  int phase_domain_index;
  real cell_centroid[ND_ND]; /* n dimensional cell center coordinates */
  real surface_level = RP_Get_Float("udf/surfacelevel");
  int flipzones = RP_Get_Integer("udf/flipzones");
  real upper_phase = 0.0;
  real lower_phase = 0.0;
  
  sub_domain_loop (subdomain, domain, phase_domain_index){ 
    if (phase_domain_index == 0){ /* primary phase */
      if (! flipzones){
        upper_phase = 1.0;
        lower_phase = 0.0;
      }
      else {
        upper_phase = 0.0;
        lower_phase = 1.0;
      }
    }
    else if (phase_domain_index == 1){
      if (! flipzones){
        upper_phase = 0.0;
        lower_phase = 1.0;
      }
      else {
        upper_phase = 1.0;
        lower_phase = 0.0;
      }
    }
    Message("Setting phasfor upper=%f, lower=%f\n", upper_phase, lower_phase);
    thread_loop_c (cell_thread, subdomain){
      /*Message ( "-- 2\n" );*/
      begin_c_loop_all (cell, cell_thread){
        /*Message ( "-- 3\n" );*/
        C_CENTROID(cell_centroid, cell, cell_thread);
        /*Message ( "-- 4\n" );*/
        if (cell_centroid[1] < surface_level){
          C_VOF(cell, cell_thread) = lower_phase; /* returns volume fraction for c cell and t phase thread */
          /*Message("Setting water\n");*/
        }
        else {
          C_VOF(cell, cell_thread) = upper_phase;
          /*Message("Setting air\n");*/
        }
        /*Message ( "-- 5\n" ); */
      } end_c_loop_all(cell_type, cell_thread)
    }
  }
  Message ( "Leaving INIT\n" );
}

DEFINE_GRID_MOTION(gridrotation,domain,dt,time,dtime)
{
  Thread *tf = DT_THREAD(dt);
  face_t f;
  Node *v;
  int n;
  real amplitude = 0.0;
  real period = 0.0;
  real x, dx;
  real pi = 3.14;
  real oldtime = time - dtime;

# if !RP_NODE /* Serial or Host */
  amplitude = RP_Get_Float("udf/amplitude");
  period = RP_Get_Float("udf/period");
# endif
  host_to_node_real_2(amplitude, period);

  /* set deforming flag on adjacent cell zone */
  SET_DEFORMING_THREAD_FLAG(THREAD_T0(tf));

  x = amplitude * (1.0 - exp((-5.0*oldtime)/(2.0*period))) * sin(2.0 * pi * oldtime / period);
  dx = (amplitude * (1.0 - exp((-5.0*time)/(2.0*period))) * sin(2.0 * pi * time / period)) - x;
  /*NV_S(vecdx, =, dx);*/

  Message ("time = %f, dtime = %f, x = %f, dx = %f\n", time, dtime, x, dx);

  begin_f_loop(f,tf){
    if (PRINCIPAL_FACE_P(f,tf)){
      f_node_loop(f,tf,n) {
        v = F_NODE(f,tf,n);
        if (NODE_POS_NEED_UPDATE (v)){
          NODE_POS_UPDATED(v);
          NODE_COORD(v)[0] += dx;
        }
      }
    }
  } end_f_loop(f,tf);
}
