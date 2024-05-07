#!/usr/bin/env python

## \file shape_optimization.py
#  \brief Python script for performing the surrogate initialization.
#  \author A. Perlini
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

import os, sys, shutil, glob, time
from optparse import OptionParser
from warnings import warn

sys.path.append(os.environ["SU2_RUN"])
import SU2

# -------------------------------------------------------------------
#  Main
# -------------------------------------------------------------------


def main():

    parser = OptionParser()
    parser.add_option(
        "-f", "--file", dest="filename", help="read config from FILE", metavar="FILE"
    )
    parser.add_option(
        "-r",
        "--name",
        dest="projectname",
        default="",
        help="try to restart from project file NAME",
        metavar="NAME",
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
        "-i",
        "--initialization",
        dest="initialization",
        default="LHS",
        help="INITIALIZATION techique (LHS)",
        metavar="INITIALIZATION",
    )
    parser.add_option(
        "-q",
        "--quiet",
        dest="quiet",
        default="True",
        help="True/False Quiet all SU2 output (optimizer output only)",
        metavar="QUIET",
    )
    parser.add_option(
        "-z",
        "--zones",
        dest="nzones",
        default="1",
        help="Number of Zones",
        metavar="ZONES",
    )

    (options, args) = parser.parse_args()

    # process inputs
    options.partitions = int(options.partitions)
    options.quiet = options.quiet.upper() == "TRUE"
    options.nzones = int(options.nzones)

    sys.stdout.write(
        "\n-------------------------------------------------------------------------\n"
    )
    sys.stdout.write(
        "|    ___ _   _ ___                                                      |\n"
    )
    sys.stdout.write(
        '|   / __| | | |_  )   Release 8.0.0 "Harrier"                           |\n'
    )
    sys.stdout.write(
        "|   \\__ \\ |_| |/ /                                                      |\n"
    )
    sys.stdout.write(
        "|   |___/\\___//___|   Surrogate Initialization Script                   |\n"
    )
    sys.stdout.write(
        "|                                                                       |\n"
    )
    sys.stdout.write(
        "-------------------------------------------------------------------------\n"
    )
    sys.stdout.write(
        "| SU2 Project Website: https://su2code.github.io                        |\n"
    )
    sys.stdout.write(
        "|                                                                       |\n"
    )
    sys.stdout.write(
        "| The SU2 Project is maintained by the SU2 Foundation                   |\n"
    )
    sys.stdout.write(
        "| (http://su2foundation.org)                                            |\n"
    )
    sys.stdout.write(
        "-------------------------------------------------------------------------\n"
    )
    sys.stdout.write(
        "| Copyright 2012-2023, SU2 Contributors (cf. AUTHORS.md)                |\n"
    )
    sys.stdout.write(
        "|                                                                       |\n"
    )
    sys.stdout.write(
        "| SU2 is free software; you can redistribute it and/or                  |\n"
    )
    sys.stdout.write(
        "| modify it under the terms of the GNU Lesser General Public            |\n"
    )
    sys.stdout.write(
        "| License as published by the Free Software Foundation; either          |\n"
    )
    sys.stdout.write(
        "| version 2.1 of the License, or (at your option) any later version.    |\n"
    )
    sys.stdout.write(
        "|                                                                       |\n"
    )
    sys.stdout.write(
        "| SU2 is distributed in the hope that it will be useful,                |\n"
    )
    sys.stdout.write(
        "| but WITHOUT ANY WARRANTY; without even the implied warranty of        |\n"
    )
    sys.stdout.write(
        "| MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU      |\n"
    )
    sys.stdout.write(
        "| Lesser General Public License for more details.                       |\n"
    )
    sys.stdout.write(
        "|                                                                       |\n"
    )
    sys.stdout.write(
        "| You should have received a copy of the GNU Lesser General Public      |\n"
    )
    sys.stdout.write(
        "| License along with SU2. If not, see <http://www.gnu.org/licenses/>.   |\n"
    )
    sys.stdout.write(
        "-------------------------------------------------------------------------\n"
    )

    surrogate_initialization(
        options.filename,
        options.projectname,
        options.partitions,
        options.initialization,
        options.quiet,
        options.nzones,
    )


#: main()


def surrogate_initialization(
    filename,
    projectname="",
    partitions=0,
    initialization="LHS",
    quiet=False,
    nzones=1,
):
    # Config
    config = SU2.io.Config(filename)
    config.NUMBER_PART = partitions
    config.NZONES = int(nzones)
    if quiet:
        config.CONSOLE = "CONCISE"

    bound_upper = float(
        config.OPT_BOUND_UPPER
    )  # variable bound to be scaled by the line search
    bound_lower = float(
        config.OPT_BOUND_LOWER
    )  # variable bound to be scaled by the line search
    relax_factor = float(config.OPT_RELAX_FACTOR)  # line search scale
    def_dv = config.DEFINITION_DV  # complete definition of the desing variable
    n_dv = sum(def_dv["SIZE"])  # number of design variables
    n_samples = int(config.SAMPLES_NUMBER)  # optimizer accuracy
    x0 = [0.0] * n_dv  # initial design
    xb_low = [
        float(bound_lower) / float(relax_factor)
    ] * n_dv  # lower dv bound it includes the line search acceleration factor
    xb_up = [
        float(bound_upper) / float(relax_factor)
    ] * n_dv  # upper dv bound it includes the line search acceleration fa
    xb = list(zip(xb_low, xb_up))  # design bounds

    # State
    state = SU2.io.State()
    state.find_files(config)

    # add restart files to state.FILES
    if (
        config.get("TIME_DOMAIN", "NO") == "YES"
        and config.get("RESTART_SOL", "NO") == "YES"
    ):
        restart_name = config["RESTART_FILENAME"].split(".")[0]
        restart_filename = (
            restart_name + "_" + str(int(config["RESTART_ITER"]) - 1).zfill(5) + ".dat"
        )
        if not os.path.isfile(
            restart_filename
        ):  # throw, if restart files does not exist
            sys.exit("Error: Restart file <" + restart_filename + "> not found.")
        state["FILES"]["RESTART_FILE_1"] = restart_filename

        # use only, if time integration is second order
        if config.get("TIME_MARCHING", "NO") == "DUAL_TIME_STEPPING-2ND_ORDER":
            restart_filename = (
                restart_name
                + "_"
                + str(int(config["RESTART_ITER"]) - 2).zfill(5)
                + ".dat"
            )
            if not os.path.isfile(
                restart_filename
            ):  # throw, if restart files does not exist
                sys.exit("Error: Restart file <" + restart_filename + "> not found.")
            state["FILES"]["RESTART_FILE_2"] = restart_filename

    # Project

    if os.path.exists(projectname):
        project = SU2.io.load_data(projectname)
        project._design_folder = 'INITIALIZATION/INIT_*'
        project.config = config
    else:
        project = SU2.opt.Project(config, state)
        project._design_folder = 'INITIALIZATION/INIT_*'
        # removing old designs
        pull, link = [], []
        folder = "./"
        with SU2.io.redirect_folder(folder, pull, link, force=True):
            folders = glob.glob(project._design_folder)
            if len(folders) > 0:
                sys.stdout.write("Removing old designs in 10s.")
                sys.stdout.flush()
                if warn:
                    time.sleep(10)
                sys.stdout.write(" Done!\n\n")
                for f in folders:
                    shutil.rmtree(f)
            #: if existing designs

            # save project
            SU2.io.save_data(project.filename, project)


    # Initialize
    if initialization == 'LHS':
        SU2.opt.LHS(project,x0,xb,n_samples)


    # rename project file
    if projectname:
        shutil.move("project.pkl", projectname)

    # rename history file
    output_format = config.TABULAR_FORMAT
    if output_format == "CSV":
        shutil.move("history_project.csv", "history_init.csv")
    else:
        shutil.move("history_project.dat", "history_init.dat")

    return project


#: shape_optimization()


# -------------------------------------------------------------------
#  Run Main Program
# -------------------------------------------------------------------

# this is only accessed if running from command prompt
if __name__ == "__main__":
    main()
