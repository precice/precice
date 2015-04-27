#include "udf.h"
#include "dynamesh_tools.h"
#include <float.h>
#include <math.h>
#include <string.h>
#include <stdio.h>

#ifdef BOOL_TRUE
#error BOOL_TRUE already defined!
#endif

#ifdef BOOL_FALSE
#error BOOL_TRUE already defined!
#endif

#define BOOL_TRUE  1
#define BOOL_FALSE 0

DEFINE_GRID_MOTION(gridmotions,domain,dt,time,dtime)
{
  Message("Setting displacements...\n");
  Thread* face_thread  = DT_THREAD(dt);
  face_t face;
  int n;
  Node* node;

  begin_f_loop (face, face_thread){
    if (PRINCIPAL_FACE_P(face,face_thread)){
      f_node_loop (face, face_thread, n){
        node = F_NODE(face, face_thread, n);
        if (NODE_POS_NEED_UPDATE(node)){
          NODE_POS_UPDATED(node);
          NODE_COORD(node)[0] += 0.1;
        }
      }
    }
  } end_f_loop (face, face_thread);
  Message("Leaving GRID_MOTION\n");
}
