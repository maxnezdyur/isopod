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
#include "MooseEnum.h"

#include <utility>

/**
 * Reporter computing connectivity-based first-difference (smoothness) regularization of a
 * parameter-mesh field.
 *
 * Unlike VectorSmoothnessRegularization, which assumes the parameter vector is an ordered 1-D
 * chain and needs a hand-supplied list of pairs to exempt a known interface, this reporter derives
 * the neighbor pairs from the parameter mesh's element connectivity. The parameter mesh is built
 * exactly as the optimization module builds it (see ParameterMesh / ParameterMeshFunction) so the
 * DOF ordering, and therefore the parameter indexing, is identical. With a CONSTANT MONOMIAL field
 * each element carries a single DOF, so an interior face between two elements is one adjacent pair;
 * higher-order fields are not supported.
 *
 * Rather than exempting a known interface, an edge-preserving penalty lets a genuine material
 * discontinuity survive wherever the data places it:
 *  - l2:            rho(d) = 0.5 d^2                                    (quadratic everywhere)
 *  - huber:         rho(d) = 0.5 d^2 for |d| <= delta, else
 *                            delta (|d| - 0.5 delta)                    (quadratic then linear)
 *  - geman_mcclure: rho(d) = 0.5 d^2 / (1 + (d/delta)^2)               (quadratic then saturating)
 * where d = v[j] - v[i] over a neighbor pair (i, j).
 *
 * Per-pair geometric weighting ('weighting' parameter):
 *  - none:               every pair weighs 1 (the historical behavior; exact on any mesh whose
 *                        pairs are geometrically interchangeable, e.g. the uniform 1-D chains).
 *  - face_over_distance: w_ij = |interior face| / |centroid distance|, normalized so the mean
 *                        weight over all pairs is 1. This is the finite-volume-consistent
 *                        anisotropy correction for 2-D/3-D grids with unequal element pitches:
 *                        differences across the tighter direction represent steeper physical
 *                        gradients and are weighted up by (pitch ratio)^2, while delta keeps its
 *                        physical meaning as a raw jump scale in the parameter's units. The mean
 *                        normalization keeps beta on the same per-pair scale as 'none'.
 *
 * Declares two Reporter values consumed by the calling input file:
 *  - "value" (Real):             J_reg = beta * sum over pairs of w rho(d)
 *  - "gradient" (std::vector<Real>): its exact analytic derivative, with
 *        gradient[j] += beta w rho'(d) and gradient[i] -= beta w rho'(d) for each pair.
 *
 * The l2 penalty on a 1-D chain mesh reproduces VectorSmoothnessRegularization with an empty
 * skip_pairs list; that equivalence is the validation of this reporter.
 */
class ParameterMeshSmoothnessRegularization : public GeneralReporter
{
public:
  static InputParameters validParams();

  ParameterMeshSmoothnessRegularization(const InputParameters & parameters);

  virtual void initialize() override {}
  virtual void execute() override;
  virtual void finalize() override {}

private:
  /// Per-edge penalty shape
  enum class PenaltyType
  {
    L2,
    HUBER,
    GEMAN_MCCLURE
  };

  /// Per-pair geometric weighting scheme
  enum class WeightingType
  {
    NONE,
    FACE_OVER_DISTANCE
  };

  /// Map the 'penalty_type' MooseEnum onto the internal enum
  static PenaltyType penaltyType(const MooseEnum & penalty_type);

  /// Map the 'weighting' MooseEnum onto the internal enum
  static WeightingType weightingType(const MooseEnum & weighting);

  /// Per-edge penalty rho(d) for the selected penalty type (excludes beta)
  Real penalty(Real d) const;
  /// Derivative rho'(d) of the per-edge penalty (excludes beta)
  Real penaltyDerivative(Real d) const;

  /// Penalty coefficient scaling the regularization; a value of 0 makes the reporter inert
  const Real _beta;

  /// Selected per-edge penalty shape
  const PenaltyType _penalty;

  /// Selected per-pair geometric weighting scheme
  const WeightingType _weighting;

  /// Transition scale separating smooth variation from an interface jump (huber, geman_mcclure)
  const Real _delta;

  /// Input vector of current parameter values being regularized
  const std::vector<Real> & _vector_reporter;

  /// Unique neighbor DOF-index pairs derived from the parameter mesh connectivity
  std::vector<std::pair<dof_id_type, dof_id_type>> _pairs;

  /// Per-pair weight (all 1 for 'none'; mean-normalized face/distance otherwise), aligned with
  /// _pairs
  std::vector<Real> _weights;

  /// Number of parameter DOFs, i.e. the expected parameter vector length
  std::size_t _n_dofs = 0;

  /// Output regularization objective value
  Real & _value;

  /// Output gradient of the regularization objective with respect to each parameter
  std::vector<Real> & _gradient;
};
