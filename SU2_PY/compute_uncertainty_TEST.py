#!/usr/bin/env python

## \file compute_uncertainty.py
#  \brief Python script for performing model-form UQ for SST turbulence model
#  \author J. Mukhopadhaya
#  \version 8.0.0 "Harrier"
#
# SU2 Project Website: https://su2code.github.io
#
# The SU2 Project is maintained by the SU2 Foundation
# (http://su2foundation.org)
#
# Copyright 2012-2023, SU2 Contributors (cf. AUTHORS.md)
#
# SU2 is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.
#
# SU2 is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with SU2. If not, see <http://www.gnu.org/licenses/>.
# imports
import numpy as np
from optparse import OptionParser
import os
import sys
import shutil
import copy
import os.path

sys.path.append(os.environ["SU2_RUN"])
import SU2


def main():
    # Command Line Options
    parser = OptionParser()
    parser.add_option(
        "-f", "--file", dest="filename", help="read config from FILE", metavar="FILE"
    )
    parser.add_option(
        "-n",
        "--partitions",
        dest="partitions",
        default=1,
        help="number of PARTITIONS",
        metavar="PARTITIONS",
    )
    parser.add_option(
        "-u",
        "--underRelaxation",
        dest="uq_urlx",
        default=0.1,
        help="under relaxation factor",
        metavar="UQ_URLX",
    )
    parser.add_option(
        "-b",
        "--deltaB",
        dest="uq_delta_b",
        default=1.0,
        help="magnitude of perturbation",
        metavar="UQ_DELTA_B",
    )

    (options, args) = parser.parse_args()
    options.partitions = int(options.partitions)
    # check the typecasting
    options.beta_delta = float(options.uq_delta_b)
    options.urlx = float(options.uq_urlx)

    # load config, start state
    config = SU2.io.Config(options.filename)
    state = SU2.io.State()

    # find solution files if they exist
    state.find_files(config)

    # prepare config
    config.NUMBER_PART = options.partitions
    config.NZONES = 1
    config.SST_OPTIONS = "UQ"
    config.UQ_DELTA_B = options.beta_delta
    config.UQ_URLX = options.urlx
    config.UQ_PERMUTE = "NO"
    
    perf_key = config.CREDIBILITY
    perf_vec = []

    # perform eigenvalue perturbations
    for comp in range(1, 4):
        print(
            "\n\n =================== Performing "
            + str(comp)
            + "  Component Perturbation =================== \n\n"
        )

        # make copies
        konfig = copy.deepcopy(config)
        ztate = copy.deepcopy(state)

        # set componentality
        konfig.UQ_COMPONENT = comp

        # send output to a folder
        folderName = str(comp) + "c/"
        if os.path.isdir(folderName):
            shutil.rmtree(folderName)
        os.mkdir(folderName)
        sendOutputFiles(konfig, folderName)

        # run su2
        SU2.run.CFD(konfig)

        # History reading
        info = historyReading(config, konfig)
        ztate.update(info)

        perf_vec.append(ztate.FUNCTIONS[perf_key])


    print(
        "\n\n =================== Performing p1c1 Component Perturbation =================== \n\n"
    )

    # make copies
    konfig = copy.deepcopy(config)
    ztate = copy.deepcopy(state)

    # set componentality
    konfig.UQ_COMPONENT = 1
    konfig.UQ_PERMUTE = "YES"

    # send output to a folder
    folderName = "p1c1/"
    if os.path.isdir(folderName):
        shutil.rmtree(folderName)
    os.mkdir(folderName)
    sendOutputFiles(konfig, folderName)

    # run su2
    SU2.run.CFD(konfig)

    # History reading
    info = historyReading(config, konfig)
    ztate.update(info)

    perf_vec.append(ztate.FUNCTIONS[perf_key])

    print(
        "\n\n =================== Performing p1c2 Component Perturbation =================== \n\n"
    )

    # make copies
    konfig = copy.deepcopy(config)
    ztate = copy.deepcopy(state)

    # set componentality
    konfig.UQ_COMPONENT = 2
    konfig.UQ_PERMUTE = "YES"

    # send output to a folder
    folderName = "p1c2/"
    if os.path.isdir(folderName):
        shutil.rmtree(folderName)
    os.mkdir(folderName)
    sendOutputFiles(konfig, folderName)

    # run su2
    SU2.run.CFD(konfig)

    # History reading
    info = historyReading(config, konfig)
    ztate.update(info)

    perf_vec.append(ztate.FUNCTIONS[perf_key])

    print('Performace vector =\t', perf_vec)
    credibility = SU2.util.ordered_bunch()
    credibility[config.CREDIBILITY + '_' + 'CREDIBILITY'] = credibility_indicator(perf_vec)
    state.FUNCTIONS.update(credibility)
    print(state)
    




def sendOutputFiles(config, folderName=""):
    config.CONV_FILENAME = folderName + config.CONV_FILENAME
    # config.BREAKDOWN_FILENAME = folderName + config.BREAKDOWN_FILENAME
    config.RESTART_FILENAME = folderName + config.RESTART_FILENAME
    config.VOLUME_FILENAME = folderName + config.VOLUME_FILENAME
    config.SURFACE_FILENAME = folderName + config.SURFACE_FILENAME


def historyReading(config, konfig=None):

    if konfig is None:
        konfig = copy.deepcopy(config)

    # multizone cases
    multizone_cases = SU2.io.get_multizone(konfig)

    # merge
    konfig["SOLUTION_FILENAME"] = konfig["RESTART_FILENAME"]
    if "FLUID_STRUCTURE_INTERACTION" in multizone_cases:
        konfig["SOLUTION_FILENAME"] = konfig["RESTART_FILENAME"]

    # filenames
    plot_format = konfig.get("TABULAR_FORMAT", "CSV")
    plot_extension = SU2.io.get_extension(plot_format)

    # adapt the history_filename, if a restart solution is chosen
    # check for 'RESTART_ITER' is to avoid forced restart situation in "compute_polar.py"...
    if konfig.get("RESTART_SOL", "NO") == "YES" and konfig.get("RESTART_ITER", 1) != 1:
        if konfig.get("CONFIG_LIST", []) != []:
            konfig[
                "CONV_FILENAME"
            ] = "config_CFD"  # master cfg is always config_CFD. Hardcoded names are prob nt ideal.
        restart_iter = "_" + str(konfig["RESTART_ITER"]).zfill(5)
        history_filename = konfig["CONV_FILENAME"] + restart_iter + plot_extension
    else:
        if konfig.get("CONFIG_LIST", []) != []:
            konfig["CONV_FILENAME"] = "config_CFD"
        history_filename = konfig["CONV_FILENAME"] + plot_extension

    special_cases = SU2.io.get_specialCases(konfig)

    # averaging final iterations
    final_avg = config.get("ITER_AVERAGE_OBJ", 0)
    # get chosen windowing function, default is square
    wnd_fct = config.get("WINDOW_FUNCTION", "SQUARE")

    # get history and objectives
    history = SU2.io.read_history(history_filename, config.NZONES)
    aerodynamics = SU2.io.read_aerodynamics(
        history_filename, config.NZONES, special_cases, final_avg, wnd_fct
    )

    # update super config
    config.update({"MATH_PROBLEM": konfig["MATH_PROBLEM"]})

    # info out
    info = SU2.io.State()
    info.FUNCTIONS.update(aerodynamics)
    info.FILES.DIRECT = konfig["RESTART_FILENAME"]
    if "INV_DESIGN_CP" in special_cases:
        info.FILES.TARGET_CP = "TargetCp.dat"
    if "INV_DESIGN_HEATFLUX" in special_cases:
        info.FILES.TARGET_HEATFLUX = "TargetHeatFlux.dat"
    info.HISTORY.DIRECT = history

    """If WINDOW_CAUCHY_CRIT is activated and the time marching converged before the final time has been reached,
       store the information for the adjoint run"""
    if config.get("WINDOW_CAUCHY_CRIT", "NO") == "YES" and config.TIME_MARCHING != "NO":
        konfig["TIME_ITER"] = int(
            info.HISTORY.DIRECT.Time_Iter[-1] + 1
        )  # update the last iteration
        if konfig["UNST_ADJOINT_ITER"] > konfig["TIME_ITER"]:
            konfig["ITER_AVERAGE_OBJ"] = max(
                0,
                konfig["ITER_AVERAGE_OBJ"]
                - (konfig["UNST_ADJOINT_ITER"] - konfig["TIME_ITER"]),
            )
            konfig["UNST_ADJOINT_ITER"] = konfig["TIME_ITER"]

        info["WND_CAUCHY_DATA"] = {
            "TIME_ITER": konfig["TIME_ITER"],
            "UNST_ADJOINT_ITER": konfig["UNST_ADJOINT_ITER"],
            "ITER_AVERAGE_OBJ": konfig["ITER_AVERAGE_OBJ"],
        }

    SU2.run.merge(konfig)

    return info

def credibility_indicator(performance_vec):
    delta_p = np.max(performance_vec) - np.min(performance_vec)
    return delta_p


if __name__ == "__main__":
    main()
