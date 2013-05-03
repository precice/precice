(rpsetvar 'dynamesh/update-in-timestep/interval 1000000000)
(ti-menu-load-string "solve dual-time-iterate 1 0")

(do () ((= (%rpgetvar 'udf/ongoing) 0))
   (if (= (%rpgetvar 'udf/iterate) 0) 
   ( begin 
      (ti-menu-load-string "solve iterate 201")
   )
   ( begin
      (do () ((= (%rpgetvar 'udf/convergence) 1))
         (ti-menu-load-string "solve iterate 202")
         (ti-menu-load-string "define/user-defined/execute-on-demand \"write_and_advance::libudf.so\"")
         (if (= (%rpgetvar 'udf/subcycling) 0)
         ( begin 
            (if (= (%rpgetvar 'udf/convergence) 0)
            ( begin 
               (rpsetvar 'dynamesh/update-in-timestep/interval 1)
               (ti-menu-load-string "solve iterate 1")
               (rpsetvar 'dynamesh/update-in-timestep/interval 1000000000)
               (rpsetvar 'udf/subcycling 1)
            ))
         ))
      )
   ))
   (ti-menu-load-string "solve dual-time-iterate 1 0")
)
