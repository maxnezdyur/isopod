//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

// Moose Includes
#include "DiracKernel.h"
#include "ReporterInterface.h"

/**

 */
class MaterialReporterPointSource : public DiracKernel, public ReporterInterface
{
public:
  static InputParameters validParams();
  MaterialReporterPointSource(const InputParameters & parameters);
  virtual void addPoints() override;

protected:
  virtual Real computeQpResidual() override;

  void errorCheck(const std::string & input_name, std::size_t reporterSize);

  /// bool if data format read in is points
  const bool _read_in_points;

  /// convenience vectors (these are not const because reporters can change their size)
  std::vector<Real> _ones_vec;
  std::vector<Real> _zeros_vec;
  std::vector<Point> _zeros_pts;
  /// x coordinate
  const std::vector<Real> & _coordx;
  /// y coordinate
  const std::vector<Real> & _coordy;
  ///z coordinate
  const std::vector<Real> & _coordz;
  ///xyz point
  const std::vector<Point> & _point;
  /// source material
  const MaterialProperty<Real> & _value;
  const Real _coeff;
};
