#!/usr/bin/env python

## \file surrogate.py
#  \brief python package for performing surrogate-based gradient estimation for credibility constraints
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

# ----------------------------------------------------------------------
#  Imports
# ----------------------------------------------------------------------

import copy
import os
import numpy as np

from .. import io as su2io

from smt.surrogate_models import KRG
from smt.utils.design_space import DesignSpace, FloatVariable

# ----------------------------------------------------------------------
#  Surrogate Modeling
# ----------------------------------------------------------------------


def surrogate(config):
    """info = SU2.run.surrogate(config)

    Runs the interface with Surrogate Modeling Toolbox
    Copyright (c) 2017, SMT developers
    All rights reserved.

    Assumptions:
        Does not rename restart filename to solution filename
        Adds 'epm' suffix to convergence filename
        config has UQ entries

    Outputs:
        info - SU2 State with keys:
            FUNCTIONS


    Updates:
        config.MATH_PROBLEM

    Executes in:
        ./
    """

    # local copy
    konfig = copy.deepcopy(config)

    # checking the function
    func_out_name = konfig["OBJECTIVE_FUNCTION"]

    # design space definition
    def_dv = konfig.DEFINITION_DV
    n_dv = sum(def_dv["SIZE"])
    bound_upper = float(konfig.OPT_BOUND_UPPER)   
    bound_lower = float(konfig.OPT_BOUND_LOWER)  
    relax_factor = float(config.OPT_RELAX_FACTOR)
    xb_low = [float(bound_lower) / float(relax_factor)] * n_dv
    xb_up = [float(bound_upper) / float(relax_factor) ] * n_dv  

    # reading the data from past EPM
    pull = []
    link = []
    i = 1

    xt = np.empty((0, n_dv))
    yt = np.array([])
    
    isPath = True
    while isPath:
        path = "../../DSN_" + str(i).zfill(3) + "/EPM"
        isPath = os.path.exists(path)

        with su2io.redirect_folder(path, pull, link) as push:
            print(os.getcwd())

            func, dv_vec = su2io.tools.read_epm("epm.dat", func_out_name)
        
        xt = np.vstack((xt, np.array(dv_vec)))
        yt = np.append(yt, func)

        i += 1
    
    """increasing the characteristic lengthscales, if the training points are less than n_dv
    the GP is kept isotropic, if they are greater or equal than n_dv the process becomes
    anisotropic
    """
    n_data = yt.shape[0]

    theta0 = [1e-2]
    if n_data >= n_dv:
        theta0 = [1e-2] * n_dv 
    
    design_space = DesignSpace([FloatVariable(xb_low, xb_up)])

    sm = KRG(design_space=design_space, theta0=theta0)
    sm.set_training_values(xt, yt)

    if n_data == 1:
        grad = [0.0] * n_dv

    else:
        sm.train()

        # querying the derivative at the last xt
        xtest = xt[-1,:]
        grad = []
        for i in n_dv:
            grad.append(float(sm.predict_derivatives(xtest, i).item()))
    
    # files out
    surr_title = "SURR_GRAD_" + func_out_name
    # info out
    info = su2io.State()
    info.HISTORY[surr_title] = grad

    return info
        




