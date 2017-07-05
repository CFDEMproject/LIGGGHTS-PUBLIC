Purpose
===================

This case is to test correct inter-process communication when using LIGGGHTS and ParScale. The case simulates a large number of individual runs.

To start, use

  mpirun -np 2 liggghts < in.liggghts_init
> mpirun -np 2 liggghts < in.liggghts_run

Clean with

> ./cleanCase.sh


