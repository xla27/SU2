#!/usr/bin/env python

## \file mmg.py
#  \brief python script for running mesh adaptation using the MMG Inria library
#  \author Victorien Menier, Brian Mungu\'ia
#  \version 7.3.0 "Blackbird"
#
# SU2 Project Website: https://su2code.github.io
#
# The SU2 Project is maintained by the SU2 Foundation
# (http://su2foundation.org)
#
# Copyright 2012-2022, SU2 Contributors (cf. AUTHORS.md)
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

import os, shutil, copy, time

from .. import io as su2io
from .. import adap as su2adap
from ..run.interface import CFD as SU2_CFD

def mmg(config):
    """
    Runs the a mesh adaptation loop with the MMG library.

    Inputs:
        config - an SU2 config object
    """

    print('SU2-MMG Anisotropic Mesh Adaptation')

    #--- Check config options related to mesh adaptation

    pyadap_options = [ 'ADAP_SIZES', 'ADAP_SUBITER', 'ADAP_HGRAD', 'ADAP_RESIDUAL_REDUCTION', 
                      'ADAP_FLOW_ITER', 'ADAP_ADJ_ITER', 'ADAP_CFL', 'ADAP_HAUSD' ]
    required_options = [ 'ADAP_SIZES', 'ADAP_SUBITER', 'ADAP_HMAX', 'ADAP_HMIN', 'MESH_FILENAME', 
                        'RESTART_SOL', 'MESH_OUT_FILENAME' ]

    if not all (opt in config for opt in required_options):
        err = '\n\n## ERROR : Missing options: \n'
        for opt in required_options:
            if not opt in config:
                err += opt + '\n'
        raise AttributeError(err)
    
    #--- NEMO solver check
    if 'NEMO' in config.SOLVER:
        nemo = True
    else:
        nemo = False

    #--- Print adap options

    print(su2adap.print_adap_options(config))

    #--- Target mesh sizes and subiterations at each size

    mesh_sizes = su2adap.get_mesh_sizes(config)
    sub_iter   = su2adap.get_sub_iterations(config)

    if len(mesh_sizes) != len(sub_iter):
        raise ValueError(f'Inconsistent number of mesh sizes and sub-iterations. {len(mesh_sizes)} mesh sizes and {len(sub_iter)} sub-iterations provided.')

    #--- Solver iterations/ residual reduction param for each size level

    flow_iter = su2adap.get_flow_iter(config)
    adj_iter  = su2adap.get_adj_iter(config)
    flow_cfl  = su2adap.get_flow_cfl(config)

    adap_sensors = su2adap.get_adap_sensors(config)
    sensor_avail = ['GOAL', 'MACH', 'PRESSURE', 'TEMPERATURE', 'ENERGY', 'DENSITY']

    for sensor in adap_sensors:
        if sensor not in sensor_avail:
            raise ValueError(f'Unknown adaptation sensor {sensor}. Available options are {sensor_avail}.')

    gol = 'GOAL' in adap_sensors

    #--- Change current directory

    warn = True
    base_dir = os.getcwd()
    adap_dir = './adap'

    if os.path.exists(adap_dir):
        print('./adap exists. Removing old mesh adaptation in 10s.')
        if warn : time.sleep(10)
        shutil.rmtree(adap_dir)
        print(f'The {adap_dir} folder was deleted.')

    dir = f'{adap_dir}/ite0'
    os.makedirs(dir)
    os.chdir(dir)
    os.symlink(os.path.join(base_dir, config.MESH_FILENAME), config.MESH_FILENAME)

    meshfil = config['MESH_FILENAME']

    #--- Format of history file

    history_format = config.TABULAR_FORMAT
    if (history_format == 'TECPLOT'):
        history_filename = os.path.join(base_dir, 'history_adap.dat')
    else:
        history_filename = os.path.join(base_dir, 'history_adap.csv')

    #--- Get mesh dimension

    dim = su2adap.get_su2_dim(meshfil)
    if ( dim != 2 and dim != 3 ):
        raise ValueError('Wrong dimension number.')

    #--- MMG parameters

    config_mmg = su2adap.get_mmg_config(config, dim)

    #--- Compute initial solution if needed, else link current files

    config_cfd = copy.deepcopy(config)
    config_cfd_ad = copy.deepcopy(config)
    for opt in pyadap_options:
        config_cfd.pop(opt, None)
        config_cfd_ad.pop(opt, None)

    #--- Check config for filenames if restarting
    restart = config['RESTART_SOL'] == 'YES'
    if restart:
        required_options = ['SOLUTION_FILENAME']
        if gol: required_options.append('SOLUTION_ADJ_FILENAME')
        if not all (opt in config for opt in required_options):
            err = 'RESTART_SOL is set to YES, but the solution is missing:\n'
            for opt in required_options:
                if not opt in config:
                    err += opt + '\n'
            raise ValueError(err)

        os.symlink(os.path.join(base_dir, config.SOLUTION_FILENAME), config.SOLUTION_FILENAME)

        print('\nInitial CFD solution is provided.')

    else:
        print('\nRunning initial CFD solution.')

    #--- Only allow ASCII restarts for file conversion
    if not gol:
        sol_ext_cfd = '.csv'
        config_cfd.OUTPUT_FILES = ['RESTART_ASCII','PARAVIEW','SURFACE_PARAVIEW']
    if gol:
        sol_ext_cfd = '.dat'
        config_cfd.OUTPUT_FILES = ['RESTART','PARAVIEW','SURFACE_PARAVIEW']
        sol_ext_cfd_ad = '.csv'
        config_cfd_ad.OUTPUT_FILES = ['RESTART_ASCII','PARAVIEW','SURFACE_PARAVIEW']


    meshfil = config['MESH_FILENAME']
    solfil  = f'restart_flow{sol_ext_cfd}'
    su2adap.set_flow_config_ini(config_cfd, solfil, adap_sensors, mesh_sizes[0])

    try: # run with redirected outputs
        #--- Run a single iteration of the flow if restarting to get history info
        if restart:
            config_cfd.ITER = 1
            config_cfd.RESTART_CFL = 'YES'

        with su2io.redirect.output('su2.out'): SU2_CFD(config_cfd)

        if restart:
            os.remove(solfil)
            os.symlink(os.path.join(base_dir, config.SOLUTION_FILENAME), solfil)

        #--- Set RESTART_SOL=YES for runs after adaptation
        if not nemo:
            config_cfd.RESTART_SOL = 'YES' 
            config_cfd.RESTART_CFL = 'YES'

        if gol:
            adjsolfil = f'restart_adj{sol_ext_cfd_ad}'
            su2adap.set_adj_config_ini(config_cfd_ad, solfil, adjsolfil, mesh_sizes[0])

            #--- If restarting, check for the existence of an adjoint restart
            if restart:
                adjsolfil_ini = config_cfd_ad.SOLUTION_ADJ_FILENAME
                func_name          = config.OBJECTIVE_FUNCTION
                suffix             = su2io.get_adjointSuffix(func_name)
                adjsolfil_ini = su2io.add_suffix(adjsolfil_ini, suffix)

                #--- Run an adjoint if the solution file doesn't exist
                if not (os.path.exists(os.path.join(base_dir, adjsolfil_ini))):
                    config_cfd_ad.ITER        = config.ITER
                    config_cfd_ad.RESTART_SOL = 'NO'

                    print('Running initial adjoint CFD solution.')

                #--- Otherwise just compute the metric
                else:
                    os.symlink(os.path.join(base_dir, adjsolfil_ini), adjsolfil_ini)
                    config_cfd_ad.ITER = 0

                    print('Initial adjoint CFD solution is provided.')

            else:
                print('Running initial adjoint CFD solution.')

            with su2io.redirect.output('su2.out'): SU2_CFD(config_cfd_ad)

            func_name      = config.OBJECTIVE_FUNCTION
            suffix         = su2io.get_adjointSuffix(func_name)
            adjsolfil = su2io.add_suffix(adjsolfil, suffix)

            #--- Set RESTART_SOL=YES for runs after adaptation  
            config_cfd_ad.RESTART_SOL = 'YES'    

    except:
        raise

    #--- Check existence of initial mesh, solution

    required_files = [meshfil, solfil]

    if not all (os.path.exists(fil) for fil in required_files):
        err = "Can't find the following files:\n"
        for fil in required_files:
            if not os.path.exists(fil):
                err += fil + '\n'
        raise Exception(err)

    #--- Start adaptive loop

    global_iter = 0

    #--- Print convergence history

    npoin = su2adap.get_su2_npoin(meshfil)
    su2adap.plot_results(history_format, history_filename, global_iter, npoin)

    print('\nStarting mesh adaptation process.\n')

    nSiz = len(mesh_sizes)
    for iSiz in range(nSiz):
        nSub = int(sub_iter[iSiz])
        for iSub in range(nSub):

            global_iter += 1

            #--- Instantiating the mesh converter
            fileconverter = su2adap.MeshSolConverter()

            #--- Load and read .su2 mesh and dump .mesh file

            fileconverter.SU2ToMeditMesh(meshfil, meshfil.replace('.su2', '.mesh'))

            #--- Load and read .csv solution and dump .sol file

            fileconverter.SU2ToMeditSol(solfil.replace(sol_ext_cfd,'.csv'), 
                                        solfil.replace('.csv', '.sol'))

            mesh_size = int(mesh_sizes[iSiz])
            if iSub == nSub-1 and iSiz != nSiz-1: mesh_size = int(mesh_sizes[iSiz+1])
            config_mmg['size'] = mesh_size

            # if gol:

            #     #--- Read and merge adjoint solution to be interpolated

            #     fileconverter.SU2ToMeditSol(adjsolfil, adjsolfil.replace('.csv', '.sol'))


            #--- Adapt mesh with MMG
            meshin = config_cfd['MESH_FILENAME'].replace('.su2','.mesh')
            meshout = config_cfd['MESH_FILENAME'].replace('.su2','_adap.mesh')
            if gol:
                adjsolfile = config_cfd_ad['RESTART_ADJ_FILENAME'].replace('.csv','.sol')
                solfile = config_cfd['RESTART_FILENAME'].replace(sol_ext_cfd,'.sol')
                su2adap.call_mmg(meshin, meshout, solfile, config_mmg)
            else:
                solfile = config_cfd['RESTART_FILENAME'].replace('.csv','.sol')
                su2adap.call_mmg(meshin, meshout, solfile, config_mmg)

            mesh_new = fileconverter.ReadMeshMedit(meshout)

            #--- Dumping a copy of the adapted mesh 
            fileconverter.WriteMeshSU2(meshout.replace('.mesh','.su2'))

            #--- Reading the new output mesh

            # extra_files=['']
            # for file in extra_files:
            #     try:
            #         os.remove(file)
            #     except OSError:
            #         pass


            #--- Print mesh sizes
            su2adap.print_adap_table(iSiz, mesh_sizes, iSub, nSub, mesh_new)

            dir = f'./ite{global_iter}'
            os.makedirs(os.path.join('..',dir))
            os.chdir(os.path.join('..',dir))

            meshfil = config_cfd['MESH_FILENAME']
            #solfil  = f'flo{sol_ext}'

            fileconverter.WriteMeshSU2(meshfil)

            # if gol:
            #     adjsolfil = f'adj{sol_ext}'
            #     sol_adj = su2adap.split_adj_sol(mesh_new)
            #     su2adap.write_sol(adjsolfil, sol_adj)

            # meshfil_gmf    = 'flo_itp.meshb'
            # solfil_gmf     = 'flo_itp.solb'
            # su2adap.write_mesh_and_sol(meshfil_gmf, solfil_gmf, mesh_new)

            del mesh_new

            # if gol:
            #     solfil_gmf_adj = 'adj_itp.solb'
            #     su2adap.write_sol(solfil_gmf_adj, sol_adj)
            #     del sol_adj

            #--- Run su2

            try: # run with redirected outputs

                solfil_ini = f'flo_ini{sol_ext_cfd}'
                # solfil_ini  = f'restart_flow{sol_ext}'
                # os.rename(solfil, solfil_ini)

                su2adap.update_flow_config(config_cfd, meshfil, solfil, solfil_ini,
                                          flow_iter[iSiz], flow_cfl[iSiz], adap_sensors, mesh_size)

                with su2io.redirect.output('su2.out'): SU2_CFD(config_cfd)

                if not os.path.exists(solfil) :
                    raise RuntimeError('SU2_CFD failed.\n')

                #--- Print convergence history

                npoin = su2adap.get_su2_npoin(meshfil)
                su2adap.plot_results(history_format, history_filename, global_iter, npoin)

                if gol:

                    adjsolfil_ini = f'adj_ini{sol_ext_cfd_ad}'
                    adjsolfil_ini = su2io.add_suffix(adjsolfil_ini, suffix)
                    # os.rename(adjsolfil, adjsolfil_ini)
                    adjsolfil_ini = f'adj_ini{sol_ext_cfd_ad}'

                    su2adap.update_adj_config(config_cfd_ad, meshfil, solfil, adjsolfil,
                                             adjsolfil_ini, adj_iter[iSiz], mesh_size)

                    with su2io.redirect.output('su2.out'): SU2_CFD(config_cfd_ad)

                    adjsolfil = su2io.add_suffix(adjsolfil, suffix)

                    if not os.path.exists(adjsolfil) :
                        raise RuntimeError('SU2_CFD_AD failed.\n')

            except:
                raise

            del fileconverter

    #--- Write final files

    fileconverter = su2adap.MeshSolConverter()
    fileconverter.SU2ToMeditMesh(meshfil, meshfil.replace('.su2', '.mesh'))

    os.rename(solfil, os.path.join(base_dir, config.RESTART_FILENAME))
    os.rename(meshfil, os.path.join(base_dir, config.MESH_OUT_FILENAME))

    pad_nul = ' '*15
    print('\nMesh adaptation successfully ended.')
    print(f'Results files: {config.MESH_OUT_FILENAME}\n{pad_nul}{config.RESTART_FILENAME}')
