#!/usr/bin/env python

## \file initialize.py
#  \brief tools for sampling for surrogate initialization
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

# -------------------------------------------------------------------
#  Imports
# -------------------------------------------------------------------

import sys

from .. import eval as su2eval
from numpy import array, zeros
import numpy as np


# -------------------------------------------------------------------
#  Initialization
# -------------------------------------------------------------------


def lhs_initialize(project,x0=None,xb=None,n_samples=10):
    """ result = lhs_initialize(project,x0=[],xb=[],n_samples=10)

        Runs the surrogate initialization with
        an SU2 project

        Inputs:
            project     - an SU2 project
            x0          - optional, initial guess
            xb          - optional, design variable bounds
            n_samples  - max number of samples, default = 10

    Outputs:
       result - the outputs from lhs sampling
    """


    # handle input cases
    if x0 is None:
        x0 = []
    if xb is None:
        xb = []


    # number of design variables
    dv_size = project.config["DEFINITION_DV"]["SIZE"]
    n_dv = sum(dv_size)
    project.n_dv = n_dv

    # Initial guess
    if not x0:
        x0 = [0.0] * n_dv

    # prescale x0
    dv_scales = project.config["DEFINITION_DV"]["SCALE"]
    k = 0
    for i, dv_scl in enumerate(dv_scales):
        for j in range(dv_size[i]):
            x0[k] = x0[k] / dv_scl
            k = k + 1

    # scale accuracy
    obj = project.config["OPT_OBJECTIVE"]
    obj_scale = []
    for this_obj in obj.keys():
        obj_scale = obj_scale + [obj[this_obj]["SCALE"]]


    # scale accuracy
    eps = 1.0e-04

    # optimizer summary
    sys.stdout.write("Initialization parameters:\n")
    sys.stdout.write(
        "Number of design variables: " + str(len(dv_size)) + " ( " + str(n_dv) + " ) \n"
    )
    sys.stdout.write("Objective function scaling factor: " + str(obj_scale) + "\n")
    sys.stdout.write("Number of samples: " + str(n_samples) + "\n")
    sys.stdout.write("Initial guess for the independent variable(s): " + str(x0) + "\n")
    sys.stdout.write(
        "Lower and upper bound for each independent variable: " + str(xb) + "\n\n"
    )

    # Latin hypercube sampling and 
    xx = lhs_sampling(n_samples, n_dv)


    for i_smp in range(n_samples): 

        x_norm = np.asarray(xx[i_smp, :]).reshape(1,n_dv)
        x = (norm_to_real(xb, x_norm)).reshape(n_dv,)

        obj_f(x, project)
    

    # inserire una funzione che entra in tutte le cartelle INIT e le va a leggere e poi salva i risultati


    # Done
    return 





def obj_f(x, project):
    """obj = obj_f(x,project)

    Objective Function
    SU2 Project interface to scipy.fmin_slsqp

    su2:         minimize f(x), list[nobj]
    scipy_slsqp: minimize f(x), float
    """

    obj_list = project.obj_f(x)
    obj = 0
    for this_obj in obj_list:
        obj = obj + this_obj

    return obj


def obj_df(x, project):
    """dobj = obj_df(x,project)

    Objective Function Gradients
    SU2 Project interface to scipy.fmin_slsqp

    su2:         df(x), list[nobj x dim]
    scipy_slsqp: df(x), ndarray[dim]
    """

    dobj_list = project.obj_df(x)
    dobj = [0.0] * len(dobj_list[0])

    for this_dobj in dobj_list:
        idv = 0
        for this_dv_dobj in this_dobj:
            dobj[idv] = dobj[idv] + this_dv_dobj
            idv += 1
    dobj = array(dobj)

    return dobj


def con_ceq(x, project):
    """cons = con_ceq(x,project)

    Equality Constraint Functions
    SU2 Project interface to scipy.fmin_slsqp

    su2:         ceq(x) = 0.0, list[nceq]
    scipy_slsqp: ceq(x) = 0.0, ndarray[nceq]
    """

    cons = project.con_ceq(x)

    if cons:
        cons = array(cons)
    else:
        cons = zeros([0])

    return cons


def con_dceq(x, project):
    """dcons = con_dceq(x,project)

    Equality Constraint Gradients
    SU2 Project interface to scipy.fmin_slsqp

    su2:         dceq(x), list[nceq x dim]
    scipy_slsqp: dceq(x), ndarray[nceq x dim]
    """

    dcons = project.con_dceq(x)

    dim = project.n_dv
    if dcons:
        dcons = array(dcons)
    else:
        dcons = zeros([0, dim])

    return dcons


def con_cieq(x, project):
    """cons = con_cieq(x,project)

    Inequality Constraints
    SU2 Project interface to scipy.fmin_slsqp

    su2:         cieq(x) < 0.0, list[ncieq]
    scipy_slsqp: cieq(x) > 0.0, ndarray[ncieq]
    """

    cons = project.con_cieq(x)

    if cons:
        cons = array(cons)
    else:
        cons = zeros([0])

    return -cons


def con_dcieq(x, project):
    """dcons = con_dcieq(x,project)

    Inequality Constraint Gradients
    SU2 Project interface to scipy.fmin_slsqp

    su2:         dcieq(x), list[ncieq x dim]
    scipy_slsqp: dcieq(x), ndarray[ncieq x dim]
    """

    dcons = project.con_dcieq(x)

    dim = project.n_dv
    if dcons:
        dcons = array(dcons)
    else:
        dcons = zeros([0, dim])

    return -dcons


def lhs_sampling(ndata, dim):
    """
    Function to perform Latin Hypercube Sampling.

    Inputs:
    - ndata is the number of required data
    - dim is the dimension of the domain

    Output:
    - V_set is the set of sampled points
    """

    # Generating dim random indpendent partitions of the set {1,...,ndata} uniformly distributed on the interval (0,ndata!)
    perm = np.zeros((dim,ndata))
    set =  np.linspace(1,ndata,ndata)
    for k in range(0,dim):
        perm[k,:] = np.random.permutation(set)
    # Vector LHS
    V_set = np.zeros((ndata,dim))
    for j in range(0,ndata):
        for k in range(0,dim):
            V_set[j,k] = (np.random.rand() + perm[k,j] - 1) / ndata

    return V_set


def norm_to_real(xb, x):
    """
    Function to transform from unit hypercube to real domain.

    Inputs:
    - xb is the zip of the domain bounds
    - x is the row vector containing the coordinates of the point in the [0,1] hypercube

    Output:
    - x_real is the row vector in the real domain
    """
    xb = list(zip(*xb))
    bndlw = np.array(xb[0])
    bndup = np.array(xb[1])
    # Transformation matrix and constant vector
    dim = len(bndup)
    T_vec = np.zeros(dim)
    b_vec = np.zeros(dim)
    for i in range(0, dim):
        T_vec[i] = bndup[i] - bndlw[i]
        b_vec[i] = bndlw[i]
    T_mat = np.diag(T_vec)
    # Transformation 
    x_real = T_mat @ x.T + b_vec.reshape(dim,1)
    # Output has to be a row vector
    x_real = x_real.T
    return x_real