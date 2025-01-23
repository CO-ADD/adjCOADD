update dplate.testplate tp
set n_reads = (
  select count(*)
  from dplate.testwell tw
  where tw.plate_id = tp.plate_id 
   and cardinality(tw.readouts) > 0
)


update dplate.testplate tp
set n_sample = (
  select count(*)
  from dplate.testwell tw
  where tw.plate_id = tp.plate_id 
   and cardinality(tw.cmpbatch_lst) > 0
)
