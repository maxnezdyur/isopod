//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ParameterMeshSmoothnessRegularization.h"

#include "AddVariableAction.h"
#include "MooseError.h"
#include "ParameterMesh.h"

#include "libmesh/system.h"
#include "libmesh/dof_map.h"
#include "libmesh/elem.h"
#include "libmesh/int_range.h"

#include <algorithm>
#include <cmath>

namespace
{
/**
 * Thin adapter that reuses ParameterMesh's exact mesh and equation-system construction (so DOF
 * index equals parameter index) and exposes the neighbor DOF-index pairs implied by the element
 * connectivity. Kept private to this translation unit because it is not user facing.
 */
class ParameterMeshConnectivity : public ParameterMesh
{
public:
  using ParameterMesh::ParameterMesh;

  /// Unique neighbor DOF-index pairs, one per interior face (recorded from the lower element id).
  /// Also fills raw_weights with the matching |interior face| / |centroid distance| for each pair.
  std::vector<std::pair<dof_id_type, dof_id_type>>
  neighborDofPairs(std::vector<Real> & raw_weights) const
  {
    const libMesh::DofMap & dof_map = _sys->get_dof_map();
    const unsigned int var = _sys->variable_number("_parameter_mesh_var");

    std::vector<std::pair<dof_id_type, dof_id_type>> pairs;
    raw_weights.clear();
    std::vector<dof_id_type> elem_dofs;
    std::vector<dof_id_type> neighbor_dofs;
    for (const auto & elem : _mesh.active_element_ptr_range())
    {
      dof_map.dof_indices(elem, elem_dofs, var);
      if (elem_dofs.size() != 1)
        mooseError("ParameterMeshSmoothnessRegularization requires a CONSTANT MONOMIAL parameter "
                   "field with exactly one DOF per element; element ",
                   elem->id(),
                   " has ",
                   elem_dofs.size(),
                   " DOFs. Higher-order parameter fields are not supported.");

      for (const auto side : make_range(elem->n_sides()))
      {
        const libMesh::Elem * const neighbor = elem->neighbor_ptr(side);
        if (!neighbor)
          continue;
        // Record each interior face once, from the element with the smaller id.
        if (elem->id() < neighbor->id())
        {
          dof_map.dof_indices(neighbor, neighbor_dofs, var);
          if (neighbor_dofs.size() != 1)
            mooseError("ParameterMeshSmoothnessRegularization requires a CONSTANT MONOMIAL "
                       "parameter field with exactly one DOF per element; element ",
                       neighbor->id(),
                       " has ",
                       neighbor_dofs.size(),
                       " DOFs. Higher-order parameter fields are not supported.");
          pairs.emplace_back(elem_dofs[0], neighbor_dofs[0]);

          const Real face = elem->build_side_ptr(side)->volume();
          const Real dist = (elem->vertex_average() - neighbor->vertex_average()).norm();
          if (dist <= 0.0)
            mooseError("ParameterMeshSmoothnessRegularization found coincident centroids between "
                       "neighboring elements ",
                       elem->id(),
                       " and ",
                       neighbor->id(),
                       "; the parameter mesh is degenerate.");
          raw_weights.push_back(face / dist);
        }
      }
    }
    return pairs;
  }
};
}

registerMooseObject("isopodApp", ParameterMeshSmoothnessRegularization);

InputParameters
ParameterMeshSmoothnessRegularization::validParams()
{
  InputParameters params = GeneralReporter::validParams();
  params.addClassDescription(
      "Connectivity-based first-difference (smoothness) regularization of a parameter-mesh field. "
      "Neighbor pairs are derived from the parameter mesh's element connectivity (CONSTANT "
      "MONOMIAL, one DOF per element), so it applies to any layout, dimension, or ordering. The "
      "edge-preserving 'huber' and 'geman_mcclure' penalties let a genuine material discontinuity "
      "survive without telling the reporter where it is; the 'l2' penalty on a 1-D chain mesh "
      "reproduces VectorSmoothnessRegularization with an empty skip_pairs list.");

  params.addRequiredParam<FileName>(
      "parameter_mesh",
      "Exodus file of the parameter mesh whose element connectivity defines the neighbor pairs.");

  MooseEnum families = AddVariableAction::getNonlinearVariableFamilies();
  families = "MONOMIAL";
  params.addParam<MooseEnum>(
      "family",
      families,
      "Family of the FE shape functions for the parameter field; must match how the optimization "
      "declares the parameters because the DOF ordering defines the parameter indexing.");

  MooseEnum orders = AddVariableAction::getNonlinearVariableOrders();
  orders = "CONSTANT";
  params.addParam<MooseEnum>(
      "order",
      orders,
      "Order of the FE shape functions for the parameter field; must match how the optimization "
      "declares the parameters because the DOF ordering defines the parameter indexing.");

  params.addRequiredParam<ReporterName>(
      "vector_reporter", "Reporter name of the vector of current parameter values to regularize.");
  params.addRequiredParam<Real>(
      "beta",
      "Penalty coefficient scaling the regularization; a value of 0 makes the reporter inert.");

  params.addParam<MooseEnum>(
      "penalty_type",
      MooseEnum("l2 huber geman_mcclure", "huber"),
      "Per-edge penalty shape: 'l2' (quadratic everywhere), 'huber' (quadratic below 'delta', "
      "linear beyond), or 'geman_mcclure' (quadratic below 'delta', saturating beyond so large "
      "jumps are nearly unpenalized).");
  params.addParam<Real>(
      "delta",
      0.0,
      "Transition scale separating smooth variation from an interface jump for the 'huber' and "
      "'geman_mcclure' penalties; must be positive for those types and is ignored for 'l2'.");

  params.addParam<MooseEnum>(
      "weighting",
      MooseEnum("none face_over_distance", "none"),
      "Per-pair geometric weight: 'none' weighs every pair 1 (historical behavior); "
      "'face_over_distance' weighs each pair by its interior-face measure divided by the "
      "neighboring elements' centroid distance, normalized to a mean weight of 1 — the "
      "finite-volume-consistent anisotropy correction for 2-D/3-D grids with unequal pitches.");

  return params;
}

ParameterMeshSmoothnessRegularization::ParameterMeshSmoothnessRegularization(
    const InputParameters & parameters)
  : GeneralReporter(parameters),
    _beta(getParam<Real>("beta")),
    _penalty(penaltyType(getParam<MooseEnum>("penalty_type"))),
    _weighting(weightingType(getParam<MooseEnum>("weighting"))),
    _delta(getParam<Real>("delta")),
    _vector_reporter(
        getReporterValueByName<std::vector<Real>>(getParam<ReporterName>("vector_reporter"))),
    _value(declareValueByName<Real>("value", REPORTER_MODE_ROOT)),
    _gradient(declareValueByName<std::vector<Real>>("gradient", REPORTER_MODE_ROOT))
{
  if ((_penalty == PenaltyType::HUBER || _penalty == PenaltyType::GEMAN_MCCLURE) && _delta <= 0.0)
    paramError("delta", "'delta' must be positive for the 'huber' and 'geman_mcclure' penalties.");

  // Build the parameter mesh exactly as the optimization module does so that DOF index equals
  // parameter index, then cache the neighbor pairs and the parameter count from its connectivity.
  const ParameterMeshConnectivity mesh(AddVariableAction::feType(parameters),
                                       getParam<FileName>("parameter_mesh"));
  std::vector<Real> raw_weights;
  _pairs = mesh.neighborDofPairs(raw_weights);
  _n_dofs = mesh.size();

  if (_weighting == WeightingType::FACE_OVER_DISTANCE && !_pairs.empty())
  {
    Real mean = 0.0;
    for (const auto w : raw_weights)
      mean += w;
    mean /= raw_weights.size();
    _weights.resize(raw_weights.size());
    for (const auto p : index_range(raw_weights))
      _weights[p] = raw_weights[p] / mean;
  }
  else
    _weights.assign(_pairs.size(), 1.0);
}

ParameterMeshSmoothnessRegularization::PenaltyType
ParameterMeshSmoothnessRegularization::penaltyType(const MooseEnum & penalty_type)
{
  if (penalty_type == "l2")
    return PenaltyType::L2;
  if (penalty_type == "huber")
    return PenaltyType::HUBER;
  return PenaltyType::GEMAN_MCCLURE;
}

ParameterMeshSmoothnessRegularization::WeightingType
ParameterMeshSmoothnessRegularization::weightingType(const MooseEnum & weighting)
{
  if (weighting == "none")
    return WeightingType::NONE;
  return WeightingType::FACE_OVER_DISTANCE;
}

Real
ParameterMeshSmoothnessRegularization::penalty(Real d) const
{
  switch (_penalty)
  {
    case PenaltyType::L2:
      return 0.5 * d * d;
    case PenaltyType::HUBER:
    {
      const Real a = std::abs(d);
      return a <= _delta ? 0.5 * d * d : _delta * (a - 0.5 * _delta);
    }
    case PenaltyType::GEMAN_MCCLURE:
    {
      const Real s = d / _delta;
      return 0.5 * d * d / (1.0 + s * s);
    }
  }
  mooseError("Unhandled penalty type in ParameterMeshSmoothnessRegularization::penalty().");
}

Real
ParameterMeshSmoothnessRegularization::penaltyDerivative(Real d) const
{
  switch (_penalty)
  {
    case PenaltyType::L2:
      return d;
    case PenaltyType::HUBER:
      return std::clamp(d, -_delta, _delta);
    case PenaltyType::GEMAN_MCCLURE:
    {
      const Real s = d / _delta;
      const Real denom = 1.0 + s * s;
      return d / (denom * denom);
    }
  }
  mooseError(
      "Unhandled penalty type in ParameterMeshSmoothnessRegularization::penaltyDerivative().");
}

void
ParameterMeshSmoothnessRegularization::execute()
{
  if (_vector_reporter.size() != _n_dofs)
    mooseError("The parameter vector has ",
               _vector_reporter.size(),
               " entries but the parameter mesh has ",
               _n_dofs,
               " DOFs.");

  _gradient.assign(_n_dofs, 0.0);
  Real objective = 0.0;
  for (const auto p : index_range(_pairs))
  {
    const auto & [i, j] = _pairs[p];
    const Real w = _weights[p];
    const Real d = _vector_reporter[j] - _vector_reporter[i];
    objective += w * penalty(d);

    const Real dp = _beta * w * penaltyDerivative(d);
    _gradient[j] += dp;
    _gradient[i] -= dp;
  }
  _value = _beta * objective;
}
