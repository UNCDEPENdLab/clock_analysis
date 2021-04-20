#!/bin/bash

#encode="TRUE FALSE"
#rt_predict="TRUE FALSE"
#alignment="clock RT"

#rt aligned
sbatch --export=alignment=RT,rt_predict=TRUE,encode=FALSE sbatch_meg_mixed_by.bash
sbatch --export=alignment=RT,rt_predict=FALSE,encode=TRUE sbatch_meg_mixed_by.bash

#clock aligned
sbatch --export=alignment=clock,rt_predict=FALSE,encode=TRUE sbatch_meg_mixed_by.bash
sbatch --export=alignment=clock,rt_predict=TRUE,encode=FALSE sbatch_meg_mixed_by.bash

