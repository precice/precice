(rpsetvar 'dynamesh/update-in-timestep/interval 1000000000)
(ti-menu-load-string "solve dual-time-iterate 1 0")

(do () ((= (%rpgetvar 'udf/ongoing) 0))
   (if (= (%rpgetvar 'udf/iterate) 0) 
   ( begin 
      (ti-menu-load-string "solve iterate 101")
   )
   ( begin
      (do () ((= (%rpgetvar 'udf/convergence) 1))
         (ti-menu-load-string "solve iterate 22")
         (ti-menu-load-string "define/user-defined/execute-on-demand \"write_and_advance::libudf\"")
         (if (= (%rpgetvar 'udf/convergence) 0)
         ( begin 
            (rpsetvar 'dynamesh/update-in-timestep/interval 1)
            (ti-menu-load-string "solve iterate 1")
            (rpsetvar 'dynamesh/update-in-timestep/interval 1000000000)
         ))
      )
   ))
   (if (= (%rpgetvar 'udf/checkpoint) 1)
   ( begin
      (rpsetvar 'udf/appendcheckpoint 1)
      (ti-menu-load-string "file/write-case-data fluent_checkpoint.cas yes")
      (rpsetvar 'udf/appendcheckpoint 0)
   ))
   (ti-menu-load-string "solve dual-time-iterate 1 0")
)
