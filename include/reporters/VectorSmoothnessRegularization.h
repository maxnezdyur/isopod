//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "GeneralReporter.h"

/**
 * Reporter computing first-difference (smoothness) regularization of an optimization
 * parameter vector held in another Reporter.
 *
 * Declares two Reporter values consumed by the calling input file:
 *  - "value" (Real): the regularization objective
 *        J_reg = 0.5 * beta * sum over active pairs i of (v[i+1] - v[i])^2
 *  - "gradient" (std::vector<Real>): the exact analytic gradient of J_reg,
 *        g[j] = beta * (active(j-1) * d_{j-1} - active(j) * d_j),
 *    where d_i = v[i+1] - v[i] and active(i) is zero when pair i is skipped or out of
 *    range and one otherwise.
 *
 * Adjacent pairs listed in "skip_pairs" are excluded from the penalty so a known physical
 * interface jump between two parameters is not penalized.
 */
class VectorSmoothnessRegularization : public GeneralReporter
{
public:
  static InputParameters validParams();

  VectorSmoothnessRegularization(const InputParameters & parameters);

  virtual void initialize() override {}
  virtual void execute() override;
  virtual void finalize() override {}

private:
  /// Penalty coefficient scaling the regularization; a value of 0 makes the reporter inert
  const Real _beta;

  /// 0-based indices of adjacent pairs (i, i+1) excluded from the penalty
  const std::vector<unsigned int> & _skip_pairs;

  /// Input vector of parameter values being regularized
  const std::vector<Real> & _vector_reporter;

  /// Output regularization objective value
  Real & _value;

  /// Output gradient of the regularization objective with respect to each parameter
  std::vector<Real> & _gradient;
};
