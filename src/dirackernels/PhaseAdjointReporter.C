//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "MathUtils.h"
#include "PhaseAdjointReporter.h"
#include <cmath>

registerMooseObject("isopodApp", PhaseAdjointReporter);

InputParameters
PhaseAdjointReporter::validParams()
{
  InputParameters params = ReporterPointSource::validParams();
  params.addClassDescription("");

  params.addCoupledVar("real_variable", "The name of the material property for the source term");
  params.addCoupledVar("imaginary_variable",
                       "The name of the material property for the source term");

  params.addRequiredParam<bool>("real_adjoint_component", "");
  params.addParam<PostprocessorName>("real_postprocessor", "");
  params.addParam<PostprocessorName>("imaginary_postprocessor", "");
  return params;
}

PhaseAdjointReporter::PhaseAdjointReporter(const InputParameters & parameters)
  : ReporterPointSource(parameters),
    _real_component(getParam<bool>("real_adjoint_component")),
    _real(coupledValue("real_variable")),
    _imaginary(coupledValue("imaginary_variable")),
    _offset_real(isParamValid("real_postprocessor") ? &getPostprocessorValue("real_postprocessor") : NULL),
    _offset_imaginary(isParamValid("imaginary_postprocessor") ? &getPostprocessorValue("imaginary_postprocessor")
                                                              : NULL),
    _offset((isParamValid("real_postprocessor") && isParamValid("imaginary_postprocessor") ? true
                                                                                           : false))
{

  if (isParamValid("real_postprocessor") != isParamValid("imaginary_postprocessor"))
    mooseError("Must supply both or none of the postprocessor offsets.");
}

Real
PhaseAdjointReporter::computeQpResidual()
{

  Real derivative = _real_component ? computeRealResidual() : computeImaginaryResidual();
  return -_test[_i][_qp] * derivative;
}

Real
PhaseAdjointReporter::computeRealResidual()
{
  auto measurement = libmesh_map_find(_point_to_weightedValue, _current_point);

  auto phase = computePhase();
  if (_offset)
    phase -= computePhaseOffset();

  auto phase_misfit = phase - measurement;

  auto real_phase_derivative = computeRealDerivative();

  if (_offset)
    real_phase_derivative -= computeRealOffsetDerivative();

  return phase_misfit * real_phase_derivative;
}

Real
PhaseAdjointReporter::computeImaginaryResidual()
{
   auto measurement = libmesh_map_find(_point_to_weightedValue, _current_point);

  auto phase = computePhase();
  if (_offset)
    phase -= computePhaseOffset();

  auto phase_misfit = phase - measurement;

  auto imaginary_phase_derivative = computeImaginaryDerivative();

  if (_offset)
    imaginary_phase_derivative -= computeImaginaryOffsetDerivative();

  return phase_misfit * imaginary_phase_derivative;
}

Real
PhaseAdjointReporter::computePhase()
{
  return atan2(_imaginary[_qp], _real[_qp]);
}

Real
PhaseAdjointReporter::computePhaseOffset()
{
  return atan2(*_offset_imaginary, *_offset_real);
}

Real
PhaseAdjointReporter::computeRealDerivative()
{
  Real numerator = -_imaginary[_qp];
  Real denominator = Utility::pow<2>(_imaginary[_qp]) + Utility::pow<2>(_real[_qp]);
  return numerator / denominator;
}

Real
PhaseAdjointReporter::computeImaginaryDerivative()
{
  Real numerator = _real[_qp];
  Real denominator = Utility::pow<2>(_imaginary[_qp]) + Utility::pow<2>(_real[_qp]);
  return numerator / denominator;
}

Real
PhaseAdjointReporter::computeRealOffsetDerivative()
{
  Real numerator = -(*_offset_imaginary);
  Real denominator = Utility::pow<2>(*_offset_imaginary) + Utility::pow<2>(*_offset_real);
  return numerator / denominator;
}

Real
PhaseAdjointReporter::computeImaginaryOffsetDerivative()
{
  Real numerator = *_offset_real;
  Real denominator = Utility::pow<2>(*_offset_imaginary) + Utility::pow<2>(*_offset_real);
  return numerator / denominator;
}
