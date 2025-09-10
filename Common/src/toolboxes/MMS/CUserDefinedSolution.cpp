/*!
 * \file CUserDefinedSolution.cpp
 * \brief Implementations of the member functions of CUserDefinedSolution.
 * \author T. Economon, E. van der Weide
 * \version 8.2.0 "Harrier"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2025, SU2 Contributors (cf. AUTHORS.md)
 *
 * SU2 is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * SU2 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with SU2. If not, see <http://www.gnu.org/licenses/>.
 */

#include "../../../include/toolboxes/MMS/CUserDefinedSolution.hpp"

CUserDefinedSolution::CUserDefinedSolution() : CVerificationSolution() {}

CUserDefinedSolution::CUserDefinedSolution(unsigned short val_nDim, unsigned short val_nVar, unsigned short val_iMesh,
                                           CConfig* config)
    : CVerificationSolution(val_nDim, val_nVar, val_iMesh, config) {

  /*--- Write a message that the solution is initialized for the manufactured
    solution for the incompressible Euler equations. ---*/

  if ((rank == MASTER_NODE) && (val_iMesh == MESH_0)) {
    cout << endl;
    cout << "Warning: Fluid properties and solution are being " << endl;
    cout << "         initialized for the manufactured solution " << endl;
    cout << "         of the incompressible Euler equations!!!" << endl;
    cout << endl << flush;
  }

  /*--- Coefficients, needed to determine the solution. ---*/
  Density = config->GetDensity_FreeStreamND();
  Temperature = config->GetTemperature_FreeStreamND();

  /*--- Constants, which describe this manufactured solution. This is a
    solution where the primitive variables vary as a combination
    of sine and cosine functions, while pressure has boundary layers ---*/

  P_0 = 1.0;
  u_0 = 1.0;
  v_0 = 1.0;
  epsilon = 0.05;

  /*--- Perform some sanity and error checks for this solution here. ---*/
  if (config->GetTime_Marching() != TIME_MARCHING::STEADY)
    SU2_MPI::Error("Steady mode must be selected for the MMS incompressible Euler case", CURRENT_FUNCTION);

  if (Kind_Solver != MAIN_SOLVER::INC_EULER && Kind_Solver != MAIN_SOLVER::INC_NAVIER_STOKES &&
      Kind_Solver != MAIN_SOLVER::INC_RANS)
    SU2_MPI::Error("Incompressible flow equations must be selected for the MMS incompressible Euler case",
                    CURRENT_FUNCTION);

  if (Kind_Solver != MAIN_SOLVER::INC_EULER)
    SU2_MPI::Error("Euler equations must be selected for the MMS incompressible Euler case", CURRENT_FUNCTION);

  if (config->GetKind_FluidModel() != CONSTANT_DENSITY)
    SU2_MPI::Error("Constant density fluid model must be selected for the MMS incompressible Euler case",
                    CURRENT_FUNCTION);

  if (config->GetEnergy_Equation())
    SU2_MPI::Error("Energy equation must be disabled (isothermal) for the MMS incompressible Euler case",
                    CURRENT_FUNCTION);
}

CUserDefinedSolution::~CUserDefinedSolution() = default;

void CUserDefinedSolution::GetBCState(const su2double* val_coords, const su2double val_t,
                                      su2double* val_solution) const {
  /*--- The exact solution is prescribed on the boundaries. ---*/
  GetSolution(val_coords, val_t, val_solution);
}

void CUserDefinedSolution::GetSolution(const su2double* val_coords, const su2double val_t,
                                       su2double* val_solution) const {
  /* Easier storage of the x- and y-coordinates. */
  const su2double x = val_coords[0];
  const su2double y = val_coords[1];

  const su2double pi = 3.14159265358979323846;
  /* Compute the primitives from the defined solution. */
  const su2double u = u_0 * (sin(x * x + y * y));
  const su2double v = v_0 * (cos(x * x + y * y));
  const su2double p = P_0 * (1 - exp(-x/epsilon)) * (1 - exp(-y/epsilon));

  /* For the incompressible solver, we return the primitive variables
    directly, as they are used for the working variables in the solver.
    Note that the implementation below is valid for both 2D and 3D. */
  val_solution[0] = p;
  val_solution[1] = u;
  val_solution[2] = v;
  val_solution[3] = 0.0;
  val_solution[nVar - 1] = Temperature;
}

void CUserDefinedSolution::GetMMSSourceTerm(const su2double* val_coords, const su2double val_t,
                                            su2double* val_source) const {
  /*--- Easier storage of the x- and y-coordinates. ---*/
  const su2double x = val_coords[0];
  const su2double y = val_coords[1];

  const su2double pi = 3.14159265358979323846;

  /*--- The expressions for the source terms are generated
    automatically by the sympy package in python.---*/
  val_source[0] =  2 * Density * (u_0 * x * cos(x*x + y*y) - v_0 * y * sin(x*x + y*y));
  val_source[1] = (2 * Density * epsilon * u_0 * (u_0 * x * sin(2 * (x*x + y*y)) - v_0 * y * pow(sin(x*x + y*y), 2) + v_0 * y * pow(cos(x*x + y*y), 2)) * exp((x + y) / epsilon) - P_0 * (1 - exp(y / epsilon))) * exp(-(x + y) / epsilon) / epsilon;
  val_source[2] = (2 * Density * epsilon * v_0 * (-u_0 * x * pow(sin(x*x + y*y), 2) + u_0 * x * pow(cos(x*x + y*y), 2) - v_0 * y * sin(2 * (x*x + y*y))) * exp((x + y) / epsilon) - P_0 * (1 - exp(x / epsilon))) * exp(-(x + y) / epsilon) / epsilon;
  val_source[3] = 0.0;
  val_source[nVar - 1] = 0.0;
}

bool CUserDefinedSolution::IsManufacturedSolution() const { return true; }
