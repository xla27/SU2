#!/usr/bin/env python

## \file scipy_tools.py
#  \brief tools for interfacing with scipy
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

# -------------------------------------------------------------------
#  Imports
# -------------------------------------------------------------------

import sys

from .. import eval as su2eval
from numpy import array, zeros


# -------------------------------------------------------------------
#  Initialization
# -------------------------------------------------------------------

def lhs_initialize(project,x0=None,xb=None,its=100):
    """ result = lhs_initialize(project,x0=[],xb=[],its=100)

        Runs the surrogate initialization with
        an SU2 project

        Inputs:
            project - an SU2 project
            x0      - optional, initial guess
            xb      - optional, design variable bounds
            its     - max outer iterations, default 100

        Outputs:
           result - the outputs from scipy.fmin_slsqp
    """

    # import sampling method
    from SU2.util import lhc_unif

    # handle input cases
    if x0 is None: x0 = []
    if xb is None: xb = []

    # function handles
    func           = obj_f
    f_eqcons       = con_ceq
    f_ieqcons      = con_cieq

    # number of design variables
    dv_size = project.config['DEFINITION_DV']['SIZE']
    n_dv = sum( dv_size)
    project.n_dv = n_dv

    # Initial guess
    if not x0: x0 = [0.0]*n_dv

    # prescale x0
    dv_scales = project.config['DEFINITION_DV']['SCALE']
    k = 0
    for i, dv_scl in enumerate(dv_scales):
        for j in range(dv_size[i]):
            x0[k] =x0[k]/dv_scl;
            k = k + 1

    # scale accuracy
    obj = project.config['OPT_OBJECTIVE']
    obj_scale = []
    for this_obj in obj.keys():
        obj_scale = obj_scale + [obj[this_obj]['SCALE']]


    # scale accuracy
    eps = 1.0e-04

    # optimizer summary
    sys.stdout.write('Initialization parameters:\n')
    sys.stdout.write('Number of design variables: ' + str(len(dv_size)) + ' ( ' + str(n_dv) + ' ) \n' )
    sys.stdout.write('Objective function scaling factor: ' + str(obj_scale) + '\n')
    sys.stdout.write('Number of samples: ' + str(its) + '\n')
    sys.stdout.write('Initial guess for the independent variable(s): ' + str(x0) + '\n')
    sys.stdout.write('Lower and upper bound for each independent variable: ' + str(xb) + '\n\n')

    xx = lhc_unif(xb, its)

    for i_smp in range(0, xx.shape[0]):
        
        x = list(xx[i,:])

        obj_f(x, project)
    

    # inserire una funzione che entra in tutte le cartelle INIT e le va a leggere e poi salva i risultati


    # Done
    return 



def obj_f(x,project):
    """ obj = obj_f(x,project)

        Objective Function
        SU2 Project interface to scipy.fmin_slsqp

        su2:         minimize f(x), list[nobj]
        scipy_slsqp: minimize f(x), float
    """

    obj_list = project.obj_f(x)
    obj = 0
    for this_obj in obj_list:
        obj = obj+this_obj

    print('obj_f =\t', obj)

    return obj

def obj_df(x,project):
    """ dobj = obj_df(x,project)

        Objective Function Gradients
        SU2 Project interface to scipy.fmin_slsqp

        su2:         df(x), list[nobj x dim]
        scipy_slsqp: df(x), ndarray[dim]
    """

    dobj_list = project.obj_df(x)
    dobj=[0.0]*len(dobj_list[0])

    for this_dobj in dobj_list:
        idv=0
        for this_dv_dobj in this_dobj:
            dobj[idv] = dobj[idv]+this_dv_dobj;
            idv+=1
    dobj = array( dobj )

    print('obj_df =\t', dobj)

    return dobj

def con_ceq(x,project):
    """ cons = con_ceq(x,project)

        Equality Constraint Functions
        SU2 Project interface to scipy.fmin_slsqp

        su2:         ceq(x) = 0.0, list[nceq]
        scipy_slsqp: ceq(x) = 0.0, ndarray[nceq]
    """

    cons = project.con_ceq(x)

    if cons: cons = array(cons)
    else:    cons = zeros([0])

    print('con_ceq =\t', cons)

    return cons

def con_dceq(x,project):
    """ dcons = con_dceq(x,project)

        Equality Constraint Gradients
        SU2 Project interface to scipy.fmin_slsqp

        su2:         dceq(x), list[nceq x dim]
        scipy_slsqp: dceq(x), ndarray[nceq x dim]
    """

    dcons = project.con_dceq(x)

    dim = project.n_dv
    if dcons: dcons = array(dcons)
    else:     dcons = zeros([0,dim])

    print('con_dceq =\t', dcons)

    return dcons

def con_cieq(x,project):
    """ cons = con_cieq(x,project)

        Inequality Constraints
        SU2 Project interface to scipy.fmin_slsqp

        su2:         cieq(x) < 0.0, list[ncieq]
        scipy_slsqp: cieq(x) > 0.0, ndarray[ncieq]
    """

    cons = project.con_cieq(x)

    if cons: cons = array(cons)
    else:    cons = zeros([0])

    print('-con_cieq =\t', -cons)

    return -cons

def con_dcieq(x,project):
    """ dcons = con_dcieq(x,project)

        Inequality Constraint Gradients
        SU2 Project interface to scipy.fmin_slsqp

        su2:         dcieq(x), list[ncieq x dim]
        scipy_slsqp: dcieq(x), ndarray[ncieq x dim]
    """

    dcons = project.con_dcieq(x)

    dim = project.n_dv
    if dcons: dcons = array(dcons)
    else:     dcons = zeros([0,dim])

    print('-con_dcieq =\t', -dcons)

    return -dcons
