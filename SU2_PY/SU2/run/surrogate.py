#!/usr/bin/env python

## \file adjoint.py
#  \brief python package for performing surrogate-based gradient estimation for credibility constraints
#  \author T. Lukaczyk, F. Palacios
#  \version 7.5.1 "Blackbird"
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

import scipy.optimize as sci_opt
import sklearn.gaussian_process as sklgp
import warnings

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
    objective = konfig["OBJECTIVE_FUNCTION"]

    # verbose
    if konfig.get("CONSOLE", "VERBOSE") in ["QUIET", "CONCISE"]:
        warnings.filterwarnings('ignore') 

 

    # design space definition and scaling factors
    def_dv = konfig.DEFINITION_DV
    n_dv = sum(def_dv["SIZE"])
    dv_scales = def_dv["SCALE"]
    dv_size   = def_dv["SIZE"]
    global_factor = float(konfig["OPT_GRADIENT_FACTOR"])
    bound_upper = float(konfig.OPT_BOUND_UPPER)   
    bound_lower = float(konfig.OPT_BOUND_LOWER)  
    relax_factor = float(config.OPT_RELAX_FACTOR)
    xb_low = [float(bound_lower) / float(relax_factor)] * n_dv
    xb_up = [float(bound_upper) / float(relax_factor) ] * n_dv  

    # reading the data from past EPM
    pull = []
    link = []


    xt = np.empty((0, n_dv))
    yt = np.array([])

    # looping on INIT folders to read EPM data
    i_init = 1
    isPath = True
    while isPath:
        path = "../../../INITIALIZATION/INIT_" + str(i_init).zfill(3) + "/EPM"
        isPath = os.path.exists(path)

        if isPath:

            with su2io.redirect_folder(path, pull, link) as push:

                func, dv_vec = su2io.tools.read_epm("epm.dat", objective)
                
            
            xt = np.vstack((xt, np.array(dv_vec)))
            yt = np.append(yt, func)

        i_init += 1

    
    # looping on DSN folders to read EPM data
    i_dsn = 1    
    isPath = True
    while isPath:
        path = "../../DSN_" + str(i_dsn).zfill(3) + "/EPM"
        isPath = os.path.exists(path)

        if isPath:

            with su2io.redirect_folder(path, pull, link) as push:

                func, dv_vec = su2io.tools.read_epm("epm.dat", objective)
                
            
            xt = np.vstack((xt, np.array(dv_vec)))
            yt = np.append(yt, func)

        i_dsn += 1
    

    n_data = yt.shape[0]

    theta0 = np.array([1e-1]).reshape(1,1)
    # if n_data >= n_dv:
    #     theta0 = np.array([1e-1] * n_dv).reshape(n_dv,1) 
    

    # GP instantiation
    kernel = sklgp.kernels.ConstantKernel(1.0, (1e-1**3, 1e1**2)) * sklgp.kernels.RBF(length_scale=theta0, length_scale_bounds=(1e-2, 1e2))
    gp = sklgp.GaussianProcessRegressor(kernel=kernel, optimizer='fmin_l_bfgs_b', n_restarts_optimizer=50, alpha=1e-7)

    if n_data == 1:

        raw_gradient = [0.0] * n_dv

        # rand = np.random.rand(n_dv)
        # raw_gradient = rand / np.linalg.norm(rand, 2)

        # # making the random first guess have a 1E-6 magnitude
        # raw_gradient = 1E-6 * raw_gradient
        # k = 0
        # for i_dv,dv_scl in enumerate(dv_scales):
        #     for i_grd in range(dv_size[i_dv]):
        #         raw_gradient[k] = raw_gradient[k] * dv_scl / global_factor 
        #         k = k + 1

    elif n_data > 1 and n_data < n_dv:
        gp.fit(xt, yt)

        # derivatives computation
        def fun_prediction(x, gp):
            fun_pred = gp.predict(x.reshape(1, -1))
            return fun_pred

        # querying the derivative at the last xt
        raw_gradient = sci_opt.approx_fprime(xt[-1,:], fun_prediction, [np.sqrt(np.finfo(float).eps)]*n_dv, gp)

        # finding the reduction factor necessary for having a raw_gradient with a norm smaller that 5e-6 after the scalings 
        scaled_grad = raw_gradient
        k = 0
        for i_dv,dv_scl in enumerate(dv_scales):
            for i_grd in range(dv_size[i_dv]):
                scaled_grad[k] = scaled_grad[k] / dv_scl 
                k = k + 1
        
        i = 0
        reduction_factor = 0.5**i
        while np.linalg.norm(scaled_grad, 2) > 5e-6 / global_factor:
            i += 1
            reduction_factor = 0.5**i
            scaled_grad = scaled_grad * 0.5
        
        raw_gradient = raw_gradient * reduction_factor


    else:
        gp.fit(xt, yt)

        # derivatives computation
        def fun_prediction(x, gp):
            fun_pred = gp.predict(x.reshape(1, -1))
            return fun_pred

        # querying the derivative at the last xt
        raw_gradient = sci_opt.approx_fprime(xt[-1,:], fun_prediction, [np.sqrt(np.finfo(float).eps)]*n_dv, gp)
    

    # info out
    info = su2io.State()


    gradients = {objective: raw_gradient}
    info.GRADIENTS.update(gradients)

    # writing the surrogate
    su2io.write_surrogate("surrogate.dat", xt, yt, raw_gradient)

    return info
        
