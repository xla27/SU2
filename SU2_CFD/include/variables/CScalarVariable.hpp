/*!
 * \file CScalarVariable.hpp
 * \brief Base class for defining the shared variables of scalar solvers.
 * \author F. Palacios, T. Economon
 * \version 8.0.1 "Harrier"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2024, SU2 Contributors (cf. AUTHORS.md)
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

#pragma once

#include "CVariable.hpp"

/*!
 * \class CScalarVariable
 * \brief Base class for defining the shared variables of scalar solvers.
 */
class CScalarVariable : public CVariable {
 protected:
  MatrixType HB_Source; /*!< \brief Harmonic Balance source term. */

  CVectorOfMatrix& Gradient_Reconstruction; /*!< \brief Reference to the gradient of the primitive variables for MUSCL
                                               reconstruction for the convective term */
  CVectorOfMatrix
      Gradient_Aux; /*!< \brief Auxiliary structure to store a second gradient for reconstruction, if required. */

  MatrixType &Primitive;     /*!< \brief Primitive variables, equivalent to Solution for non-SST. */
  MatrixType Primitive_Aux;  /*!< \brief Primitive variables, equivalent to Solution for non-SST. */


 public:
  /*!
   * \brief Constructor of the class.
   * \param[in] npoint - Number of points/nodes/vertices in the domain.
   * \param[in] ndim - Number of dimensions of the problem.
   * \param[in] nvar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CScalarVariable(unsigned long npoint, unsigned long ndim, unsigned long nvar, const CConfig* config);

  /*!
   * \brief Get the array of the reconstruction variables gradient at a node.
   * \param[in] iPoint - Index of the current node.
   * \return Array of the reconstruction variables gradient at a node.
   */
  inline CMatrixView<su2double> GetGradient_Reconstruction(unsigned long iPoint) final {
    return Gradient_Reconstruction[iPoint];
  }

  /*!
   * \brief Get the reconstruction gradient for primitive variable at all points.
   * \return Reference to variable reconstruction gradient.
   */
  inline CVectorOfMatrix& GetGradient_Reconstruction() final { return Gradient_Reconstruction; }
  inline const CVectorOfMatrix& GetGradient_Reconstruction() const final { return Gradient_Reconstruction; }

  /*!
   * \brief Set the harmonic balance source term.
   * \param[in] iPoint - Point index.
   * \param[in] iVar - Index of the variable.
   * \param[in] source - Value of the harmonic balance source term. for the index <i>iVar</i>.
   */
  inline void SetHarmonicBalance_Source(unsigned long iPoint, unsigned long iVar, su2double source) final {
    HB_Source(iPoint, iVar) = source;
  }

  /*!
   * \brief Get the harmonic balance source term.
   * \param[in] iPoint - Point index.
   * \param[in] iVar - Index of the variable.
   * \return Value of the harmonic balance source term for the index <i>val_var</i>.
   */
  inline su2double GetHarmonicBalance_Source(unsigned long iPoint, unsigned long iVar) const final {
    return HB_Source(iPoint, iVar);
  }

  /*!
   * \brief Get the value of the mass diffusivity
   * \param[in] iPoint - Point index.
   * \param[in] val_ivar - eqn. index to the mass diffusivity.
   * \return Value of the mass diffusivity
   */
  inline virtual su2double GetDiffusivity(unsigned long iPoint, unsigned short val_ivar) const { return 0.0; }

  /*!
   * \brief Get the primitive variables for all points.
   * \return Reference to primitives.
   */
  inline const MatrixType& GetPrimitive(void) const final { return Primitive; }

  /*!
   * \brief Get the primitive variables.
   * \param[in] iVar - Index of the variable.
   * \return Value of the primitive variable for the index <i>iVar</i>.
   */
  inline su2double GetPrimitive(unsigned long iPoint, unsigned long iVar) const final { return Primitive(iPoint,iVar); }

  /*!
   * \brief Set the value of the primitive variables.
   * \param[in] iVar - Index of the variable.
   * \param[in] iVar - Index of the variable.
   * \return Set the value of the primitive variable for the index <i>iVar</i>.
   */
  inline void SetPrimitive(unsigned long iPoint, unsigned long iVar, su2double val_prim) final { Primitive(iPoint,iVar) = val_prim; }

  /*!
   * \brief Set the value of the primitive variables.
   * \param[in] val_prim - Primitive variables.
   * \return Set the value of the primitive variable for the index <i>iVar</i>.
   */
  inline void SetPrimitive(unsigned long iPoint, const su2double *val_prim) final {
    for (unsigned long iVar = 0; iVar < nPrimVar; iVar++)
      Primitive(iPoint,iVar) = val_prim[iVar];
  }

  /*!
   * \brief Get the primitive variables of the problem.
   * \return Pointer to the primitive variable vector.
   */
  inline su2double *GetPrimitive(unsigned long iPoint) final {return Primitive[iPoint]; }
};
