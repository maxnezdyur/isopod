//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "VectorSmoothnessRegularization.h"

#include "libmesh/int_range.h"

registerMooseObject("isopodApp", VectorSmoothnessRegularization);

InputParameters
VectorSmoothnessRegularization::validParams()
{
  InputParameters params = GeneralReporter::validParams();
  params.addClassDescription(
      "First-difference (smoothness) regularization of an optimization parameter vector. Declares "
      "a scalar 'value' equal to J_reg = 0.5 * beta * sum over the non-skipped adjacent pairs i of "
      "(v[i+1] - v[i])^2, and its exact analytic 'gradient' with respect to each parameter. "
      "Adjacent pairs listed in 'skip_pairs' are excluded from the penalty so a known physical "
      "interface jump between two parameters is not penalized.");
  params.addRequiredParam<ReporterName>(
      "vector_reporter", "Reporter name of the vector of parameter values to regularize.");
  params.addRequiredParam<Real>(
      "beta",
      "Penalty coefficient scaling the regularization; a value of 0 makes the reporter inert.");
  params.addParam<std::vector<unsigned int>>(
      "skip_pairs",
      "0-based indices i of adjacent pairs (i, i+1) excluded from the penalty, e.g. the pair "
      "straddling a known material interface.");
  return params;
}

VectorSmoothnessRegularization::VectorSmoothnessRegularization(const InputParameters & parameters)
  : GeneralReporter(parameters),
    _beta(getParam<Real>("beta")),
    _skip_pairs(getParam<std::vector<unsigned int>>("skip_pairs")),
    _vector_reporter(
        getReporterValueByName<std::vector<Real>>(getParam<ReporterName>("vector_reporter"))),
    _value(declareValueByName<Real>("value", REPORTER_MODE_ROOT)),
    _gradient(declareValueByName<std::vector<Real>>("gradient", REPORTER_MODE_ROOT))
{
}

void
VectorSmoothnessRegularization::execute()
{
  const std::size_t n = _vector_reporter.size();

  // With fewer than two entries there are no adjacent pairs to penalize.
  if (n < 2)
  {
    _value = 0.0;
    _gradient.assign(n, 0.0);
    return;
  }

  const std::size_t num_pairs = n - 1;

  // Flag which adjacent pairs contribute to the penalty.
  std::vector<bool> active(num_pairs, true);
  for (const auto pair : _skip_pairs)
  {
    if (pair >= num_pairs)
      paramError("skip_pairs",
                 "Pair index ",
                 pair,
                 " is out of range; the input vector has ",
                 n,
                 " entries, so valid pair indices are 0 to ",
                 num_pairs - 1,
                 ".");
    active[pair] = false;
  }

  // First differences d_i = v[i+1] - v[i].
  std::vector<Real> d(num_pairs);
  for (const auto i : make_range(num_pairs))
    d[i] = _vector_reporter[i + 1] - _vector_reporter[i];

  // Objective: J_reg = 0.5 * beta * sum over active pairs of d_i^2.
  Real sum = 0.0;
  for (const auto i : make_range(num_pairs))
    if (active[i])
      sum += d[i] * d[i];
  _value = 0.5 * _beta * sum;

  // Analytic gradient: g[j] = beta * (active(j-1) * d_{j-1} - active(j) * d_j),
  // dropping any term whose referenced pair is out of range or skipped.
  _gradient.assign(n, 0.0);
  for (const auto j : make_range(n))
  {
    const Real left = (j >= 1 && active[j - 1]) ? d[j - 1] : 0.0;
    const Real right = (j + 1 < n && active[j]) ? d[j] : 0.0;
    _gradient[j] = _beta * (left - right);
  }
}
