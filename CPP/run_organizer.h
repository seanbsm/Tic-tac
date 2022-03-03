#ifndef RUN_ORGANIZER_H
#define RUN_ORGANIZER_H

#include <iostream>
#include <vector>

#include "constants.h"
#include "type_defs.h"
#include "disk_io_routines.h"
#include "General_functions/kinetic_conversion.h"

void find_on_shell_bins(solution_configuration& solve_config,
                        swp_statespace swp_states,
                        run_params run_parameters);

void find_deuteron_channels(solution_configuration& solve_config,
                            pw_3N_statespace pw_states);

#endif // RUN_ORGANIZER_H