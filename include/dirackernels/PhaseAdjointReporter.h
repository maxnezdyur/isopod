
//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html
#pragma once

#include "ReporterPointSource.h"

/**
 *
 */
class PhaseAdjointReporter : public ReporterPointSource
{
public:
  static InputParameters validParams();

  PhaseAdjointReporter(const InputParameters & parameters);


  virtual Real computeQpResidual() override;

  protected:

  Real computeRealResidual();

  Real computeImaginaryResidual();

  Real computePhase();
  Real computePhaseOffset();
  Real computeRealDerivative();
  Real computeImaginaryDerivative();
    Real computeRealOffsetDerivative();
  Real computeImaginaryOffsetDerivative();



  const bool _real_component;

  const VariableValue & _real;

  const VariableValue & _imaginary;
  const PostprocessorValue * _offset_real;
  const PostprocessorValue * _offset_imaginary;
  const bool _offset;

};
