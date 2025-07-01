//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "VariableReporterPointSource.h"
#include "MooseUtils.h"

registerMooseObject("isopodApp", VariableReporterPointSource);

InputParameters
VariableReporterPointSource::validParams()
{
  InputParameters params = DiracKernel::validParams();

  params.addClassDescription("Apply a point load defined by Reporter.");

  params.addCoupledVar(
      "source_variable", "The name of the material property for the source term");
  params.addParam<ReporterName>(
      "x_coord_name",
      "reporter x-coordinate name.  This uses the reporter syntax <reporter>/<name>.");
  params.addParam<ReporterName>(
      "y_coord_name",
      "reporter y-coordinate name.  This uses the reporter syntax <reporter>/<name>.");
  params.addParam<ReporterName>(
      "z_coord_name",
      "reporter z-coordinate name.  This uses the reporter syntax <reporter>/<name>.");
  params.addParam<ReporterName>("point_name",
                                "reporter point name.  This uses the reporter syntax "
                                "<reporter>/<name>.");
  params.addParam<Real>("coeff", 1.0, "Coeff");
  return params;
}

VariableReporterPointSource::VariableReporterPointSource(const InputParameters & parameters)
  : DiracKernel(parameters),
    ReporterInterface(this),
    _read_in_points(isParamValid("point_name")),
    _coordx(isParamValid("x_coord_name")
                ? getReporterValue<std::vector<Real>>("x_coord_name", REPORTER_MODE_REPLICATED)
                : _zeros_vec),
    _coordy(isParamValid("y_coord_name")
                ? getReporterValue<std::vector<Real>>("y_coord_name", REPORTER_MODE_REPLICATED)
                : _zeros_vec),
    _coordz(isParamValid("z_coord_name")
                ? getReporterValue<std::vector<Real>>("z_coord_name", REPORTER_MODE_REPLICATED)
                : _zeros_vec),
    _point(_read_in_points
               ? getReporterValue<std::vector<Point>>("point_name", REPORTER_MODE_REPLICATED)
               : _zeros_pts),
    _value(coupledValue("source_variable")),
    _coeff(getParam<Real>("coeff"))
{
  if (isParamValid("point_name") == (isParamValid("x_coord_name") || isParamValid("y_coord_name") ||
                                     isParamValid("z_coord_name")))
    paramError("Either supply x,y, and z reporters or a point reporter.");
}

void
VariableReporterPointSource::addPoints()
{

  if (_read_in_points)
    for (const auto & i : index_range(_point))
      addPoint(_point[i]);
  else
  {
    if ((_coordx.size() != _coordy.size()) || (_coordx.size() != _coordz.size()))
      mooseError("Coords don't match.");
    std::vector<Point> points(_coordx.size());
    for (const auto i : index_range(_coordx))
    {
      const Point point = Point(_coordx[i], _coordy[i], _coordz[i]);
      addPoint(point);
    }
  }
}

Real
VariableReporterPointSource::computeQpResidual()
{
  // This is negative because it's a forcing function that has been brought over to the left side
  return -_coeff * _test[_i][_qp] * _value[_qp];
}
