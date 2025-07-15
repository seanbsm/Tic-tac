#ifndef RUN_ORGANIZER_H
#define RUN_ORGANIZER_H

#include <iostream>
#include <vector>

#include "utils/constants.h"
#include "utils/type_defs.h"
#include "utils/disk_io_routines.h"
#include "utils/kinetic_conversion.h"
#include "core/make_pw_symm_states.h"

void find_on_shell_bins(solution_configuration& solve_config,
						channel_os_indexing&	solve_config_subchn,
						pw_3N_statespace pw_states,
                        fwp_statespace fwp_states,
                        swp_statespace swp_states,
                        run_params run_parameters);

void find_deuteron_channels(solution_configuration& solve_config,
                            pw_3N_statespace pw_states);

#endif // RUN_ORGANIZER_H