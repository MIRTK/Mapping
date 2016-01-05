// =============================================================================
// Project: Image Registration Toolkit (IRTK)
// Package: Registration
//
// Copyright (c) 2015 Imperial College London
// Copyright (c) 2015 Andreas Schuh
// =============================================================================

#include <mirtkCommon.h>

// =============================================================================
// Help
// =============================================================================

// -----------------------------------------------------------------------------
void PrintHelp(const char* name)
{
  using namespace mirtk;
  cout << endl;
  cout << "usage: " << name << " <target> <output> [options]" << endl;
  cout << "       " << name << " [<target>] -eval <file> [options]" << endl;
  cout << endl;
  cout << "This tool parameterizes an shape embedded in 3D space. The generated output" << endl;
  cout << "map assigns each point in the input domain a point in the output domain." << endl;
  cout << endl;
  cout << "The input domain is defined either by the convex hull of an input point" << endl;
  cout << "cloud or a surface mesh delineating the boundary of the domain, or a volumetric" << endl;
  cout << "data set such as a tetrahedral mesh. Furthermore, this tool can parameterize" << endl;
  cout << "common primitive shapes such as a sphere by either triangulating its surface" << endl;
  cout << "or creating a tetrahedral mesh of its interior. A triangulated boundary surface" << endl;
  cout << "generated with this tool can subsequently be mapped to the boundary of the desired" << endl;
  cout << "output domain. Such boundary map can then be extended into the interior of the" << endl;
  cout << "input domain to generate a volumetric map." << endl;
  cout << endl;
  cout << "The output is by default either the tesselation of the target point set" << endl;
  cout << "or, when a volumetric map to a source domain is computed, the target point" << endl;
  cout << "set (or tesselation of the target volume) with the corresponding volumetric" << endl;
  cout << "map coordinates as point data. When the -map option is given, the coefficients" << endl;
  cout << "of the volumetric map itself are written to a binary output file. This map can" << endl;
  cout << "later applied to the target or another point set using the -apply option." << endl;
  cout << endl;
  cout << "Required arguments:" << endl;
  cout << "  target                      Point set delineating the boundary of the input domain" << endl;
  cout << "                              or name of primitive shape such as: \"sphere\"." << endl;
  cout << "  output                      File path of output tesselation, volumetric map," << endl;
  cout << "                              or reparameterized target point set (see -apply)." << endl;
  cout << endl;
  cout << "Optional arguments:" << endl;
  cout << "  -radius <float>             Radius of sphere to parameterize. (default: 1)" << endl;
  cout << "  -surface                    Tesselate surface of the shape. (default: off)" << endl;
  cout << "  -volume                     Tesselate or parameterize the interior. (default: on)" << endl;
  cout << "  -source <file>              Point set delineating the boundary of the output domain." << endl;
  cout << "  -array <name>               Name of point data array used as boundary map. (default: tcoords)" << endl;
  cout << "  -iterations <n>             No. of discrete solver iterations. (default: 500)" << endl;
  cout << "  -map                        Write volumetric map to binary output file. (default: off)" << endl;
  cout << "  -apply <file>               Apply a previously computed volumetric map to the" << endl;
  cout << "                              points of the input target point set. (default: off)" << endl;
  cout << "  -eval <file>                Evaluate a previously computed volumetric map." << endl;
  cout << "                              If a target point set is given, the volumetric map" << endl;
  cout << "                              is evaluated at each target point and stored as tcoords" << endl;
  cout << "                              of the output point set. (default: off)" << endl;
  cout << "  -harmonic-energy [<name>]   Evaluate harmonic energy of the volumetric map." << endl;
  cout << "                              Optionally write harmonic energy at each point" << endl;
  cout << "                              of a volumetric grid of the input domain to the" << endl;
  cout << "                              named image file. (default: off)" << endl;
  cout << "  -lattice <file>             Lattice attributes used to discretize the input" << endl;
  cout << "                              domain on a regular grid are read from the given" << endl;
  cout << "                              image file. (default: derived from input domain)" << endl;
  cout << endl;
  cout << "Output map type:" << endl;
  cout << "  -harmonic                   Harmonic volumetric map." << endl;
  cout << "  -biharmonic                 Biharmonic volumetric map." << endl;
  cout << "  -acap                       As-conformal-as-possible volumetric map." << endl;
  PrintCommonOptions(cout);
  cout << endl;
}

// =============================================================================
// Includes
// =============================================================================

#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <vtkPointData.h>
#include <vtkDataArray.h>
#include <vtkFloatArray.h>
#include <vtkSphereSource.h>
#include <vtkImplicitPolyDataDistance.h>

#include <mirtkGenericImage.h>
#include <mirtkGradientImageFilter.h>
#include <mirtkPointSetUtils.h>

#include <mirtkDiscreteMap.h>
#include <mirtkHarmonicTetrahedralVolumeParameterizer.h>
#include <mirtkHarmonicFundamentalVolumeParameterizer.h>
//#include <mirtkBiharmonicFundamentalVolumeParameterizer.h>
#include <mirtkAsConformalAsPossibleVolumeParameterizer.h>

using namespace mirtk;

// =============================================================================
// Auxiliaries
// =============================================================================

// -----------------------------------------------------------------------------
/// Enumeration of available boundary and/or volumetric maps
enum MapType
{
  MAP_Spherical,     ///< Map target surface to a sphere
  MAP_Spectral,      ///< Map boundary surfaces using spectral matching
  MAP_Harmonic,      ///< Compute harmonic volumetric map
  MAP_HarmonicFEM,   ///< Compute tetrahedral harmonic volumetric map
  MAP_HarmonicMFS,   ///< Compute harmonic volumetric map using MFS
  MAP_BiharmonicMFS, ///< Compute biharmonic volumetric map using MFS
  MAP_ACAP           ///< Compute as-conformal-as-possible volumetric map
};

// -----------------------------------------------------------------------------
/// Get boundary map array from target point set
vtkDataArray *GetBoundaryMap(vtkPointSet *target, const char *name)
{
  vtkPointData *targetPD     = target->GetPointData();
  vtkDataArray *boundary_map = targetPD->GetArray(name);
  if (!boundary_map) {
    if      (strcmp(name, "tcoords") == 0) boundary_map = targetPD->GetTCoords();
    else if (strcmp(name, "vectors") == 0) boundary_map = targetPD->GetVectors();
    else if (strcmp(name, "scalars") == 0) boundary_map = targetPD->GetScalars();
  }
  if (!boundary_map) {
    cerr << "Error: Target point set has no point data array named " << name << endl;
    exit(1);
  }
  return boundary_map;
}

// -----------------------------------------------------------------------------
/// Create tesselation of the surface of the sphere
vtkSmartPointer<vtkPolyData> Sphere(double radius)
{
  vtkSmartPointer<vtkSphereSource> source = vtkSmartPointer<vtkSphereSource>::New();
  source->SetRadius(radius);
  source->Update();
  return source->GetOutput();
}

// -----------------------------------------------------------------------------
/// Create or read named point set
vtkSmartPointer<vtkPointSet> GetPointSet(const char *name, double radius = 1.0)
{
  vtkSmartPointer<vtkPointSet> pointset;
  if (strcmp(name, "sphere") == 0) {
    pointset = Sphere(radius);
  } else {
    if (verbose) cout << "Reading point set from " << name << "...", cout.flush();
    pointset = ReadPointSet(name);
    if (pointset->GetNumberOfPoints() == 0) {
      if (verbose) cout << " failed" << endl;
      cerr << "Error: Failed to read point set from " << name << " or it has no points" << endl;
      exit(1);
    }
    if (verbose) cout << " done" << endl;
    if (pointset->GetNumberOfCells() == 0) {
      if (verbose) cout << "Computing convex hull...", cout.flush();
      pointset = ConvexHull(pointset);
      if (verbose) cout << " done" << endl;
    }
  }
  return pointset;
}

// -----------------------------------------------------------------------------
/// Map target surface onto sphere
vtkSmartPointer<vtkDataArray> SphericalMap(vtkSmartPointer<vtkPointSet> target,
                                           double                       radius = 1.0)
{
  cerr << "SphericalMap: not implemented" << endl;
  exit(1);

  // Initialize output boundary map
  vtkSmartPointer<vtkDataArray> boundary_map;
  boundary_map = vtkSmartPointer<vtkFloatArray>::New();
  boundary_map->SetName("BoundaryMap");
  boundary_map->SetNumberOfComponents(3);
  boundary_map->SetNumberOfTuples(target->GetNumberOfPoints());

  // TODO: Implement spherical map such as spherical Isomap (cf. MATLAB implementation)

  return boundary_map;
}

// -----------------------------------------------------------------------------
/// Map target surface onto sphere
vtkSmartPointer<vtkDataArray> SpectralMatch(vtkSmartPointer<vtkPointSet> target,
                                            vtkSmartPointer<vtkPointSet> source)
{
  cerr << "SpectralMatch: not implemented" << endl;
  exit(1);

  // Initialize output boundary map
  vtkSmartPointer<vtkDataArray> boundary_map;
  boundary_map = vtkSmartPointer<vtkFloatArray>::New();
  boundary_map->SetName("BoundaryMap");
  boundary_map->SetNumberOfComponents(3);
  boundary_map->SetNumberOfTuples(target->GetNumberOfPoints());

  // TODO: Implement spectral matching (cf. mirtkSpectralDecomposition.h)

  return boundary_map;
}

// -----------------------------------------------------------------------------
/// Map target surface onto source surface
vtkSmartPointer<vtkDataArray> BoundaryMap(vtkSmartPointer<vtkPointSet> target,
                                          vtkSmartPointer<vtkPointSet> source,
                                          double                       radius,
                                          MapType                      type)
{
  // Extract target surface
  vtkSmartPointer<vtkPolyData>  target_surface = DataSetSurface(target, true);
  vtkSmartPointer<vtkDataArray> surface_map;

  // Compute spherical map
  if (type == MAP_Spherical) {
    surface_map = SphericalMap(target, radius);
  // Perform surface matching
  } else {
    // Extract source surface
    vtkSmartPointer<vtkPolyData> source_surface = DataSetSurface(source);

    // Compute spectral match
    if (type == MAP_Spectral) {
      surface_map = SpectralMatch(target, source);
    }
  }

  // Convert surface map to output array
  vtkSmartPointer<vtkDataArray> boundary_map;
  if (surface_map) {
    if (surface_map->GetNumberOfTuples() == target->GetNumberOfPoints()) {
      boundary_map = surface_map;
    } else {
      boundary_map = vtkSmartPointer<vtkFloatArray>::New();
      boundary_map->SetNumberOfComponents(surface_map->GetNumberOfComponents());
      boundary_map->SetNumberOfTuples(target->GetNumberOfPoints());
      vtkDataArray *origPtIds = target_surface->GetPointData()->GetArray("vtkOriginalPointIds");
      for (vtkIdType ptId = 0, origPtId; ptId < surface_map->GetNumberOfTuples(); ++ptId) {
        origPtId = static_cast<vtkIdType>(origPtIds->GetComponent(ptId, 0));
        boundary_map->SetTuple(origPtId, surface_map->GetTuple(ptId));
      }
    }
    boundary_map->SetName("BoundaryMap");
  }

  return boundary_map;
}

// -----------------------------------------------------------------------------
/// Parameterize interior of target point set
VolumetricMap *ComputeVolumetricMap(vtkSmartPointer<vtkPointSet>  target,
                                    vtkSmartPointer<vtkDataArray> boundary_map,
                                    vtkSmartPointer<vtkDataArray> boundary_mask,
                                    MapType                       type,
                                    int                           no_of_iterations)
{
  VolumetricMap *map = NULL;
  if (type == MAP_Harmonic) {
    if (IsTetrahedralMesh(target)) type = MAP_HarmonicFEM;
    else                           type = MAP_HarmonicMFS;
  }
  switch (type) {
    case MAP_HarmonicFEM: {
      if (verbose) cout << "Computing tetrahedral harmonic map...", cout.flush();
      HarmonicTetrahedralVolumeParameterizer filter;
      filter.InputSet(target);
      filter.InputMap(boundary_map);
      filter.InputMask(boundary_mask);
      filter.NumberOfIterations(no_of_iterations);
      filter.Run();
      map = filter.GetOutputMap();
    } break;
    case MAP_HarmonicMFS: {
      if (verbose) cout << "Computing harmonic map using MFS...", cout.flush();
      HarmonicFundamentalVolumeParameterizer filter;
      filter.InputSet(target);
      filter.InputMap(boundary_map);
      filter.Run();
      map = filter.GetOutputMap();
    } break;
    case MAP_BiharmonicMFS: {
      cerr << "Error: Biharmonic mapping not implemented" << endl;
      exit(1);
//      if (verbose) cout << "Computing biharmonic map using MFS...", cout.flush();
//      BiharmonicFundamentalVolumeParameterizer filter;
//      filter.InputSet(target);
//      filter.InputMap(boundary_map);
//      filter.Run();
//      map = filter.GetOutputMap();
    } break;
    case MAP_ACAP: {
      if (verbose) cout << "Computing as-conformal-as-possible map...", cout.flush();
      AsConformalAsPossibleVolumeParameterizer filter;
      filter.InputSet(target);
      filter.InputMap(boundary_map);
      filter.Run();
      map = filter.GetOutputMap();
    } break;
    default:
      cerr << "Error: Invalid volumetric map type: " << type << endl;
      exit(1);
  }
  if (verbose) cout << " done" << endl;
  return map;
}

// -----------------------------------------------------------------------------
/// Normalize target and source models to fit into the unit box
template <class Real>
void NormalizeCoordinates(GenericImage<Real>          &map,
                          vtkSmartPointer<vtkPointSet> target,
                          vtkSmartPointer<vtkPointSet> source)
{
  double bounds[6];

  // Get homogeneous transformation matrix to map target points to unit box
  Matrix T(4, 4);
  target->GetBounds(bounds);
  T(0, 0) = 1.0 / (bounds[1] - bounds[0]);
  T(1, 1) = 1.0 / (bounds[3] - bounds[2]);
  T(2, 2) = 1.0 / (bounds[5] - bounds[4]);
  T(0, 3) = - bounds[0] / (bounds[1] - bounds[0]);
  T(1, 3) = - bounds[2] / (bounds[3] - bounds[2]);
  T(2, 3) = - bounds[4] / (bounds[5] - bounds[4]);
  T(3, 3) = 1.0;

  // Get homogeneous transformation matrix to map source points to unit box
  Matrix S(4, 4);
  source->GetBounds(bounds);
  S(0, 0) = 1.0 / (bounds[1] - bounds[0]);
  S(1, 1) = 1.0 / (bounds[3] - bounds[2]);
  S(2, 2) = 1.0 / (bounds[5] - bounds[4]);
  S(0, 3) = - bounds[0] / (bounds[1] - bounds[0]);
  S(1, 3) = - bounds[2] / (bounds[3] - bounds[2]);
  S(2, 3) = - bounds[4] / (bounds[5] - bounds[4]);
  S(3, 3) = 1.0;

  // Adjust attributes of discrete target domain lattice
  map.PutAffineMatrix(T, true);

  // Map source domain to unit box
  const int nvox = map.NumberOfSpatialVoxels();
  Real *x = map.Data(), *y = x + nvox, *z = y + nvox;
  for (int vox = 0; vox < nvox; ++vox, ++x, ++y, ++z) {
    Transform(S, *x, *y, *z);
  }
}

// -----------------------------------------------------------------------------
/// Evaluate harmonic energy of discretized volumetric map
///
/// Li et al. (2009). Meshless harmonic volumetric mapping using fundamental solution methods.
/// IEEE Transactions on Automation Science and Engineering, 6(3), 409–422.
template <class Real>
double EvaluateHarmonicEnergy(const GenericImage<Real> &map,
                              GenericImage<Real>       *energy = NULL)
{
  typedef GradientImageFilter<Real> GradientFilter;

  const ImageAttributes &lattice = map.Attributes();
  if (energy) energy->Initialize(lattice, 1);

  // Compute squared norm of gradient vectors for each component
  GenericImage<Real> squared_gradient_norm;
  GradientFilter gradient(GradientFilter::GRADIENT_DOT_PRODUCT);
  gradient.Input(&map);
  gradient.Output(&squared_gradient_norm);
  gradient.UseVoxelSize(true);
  gradient.UseOrientation(false); // R = I
  gradient.PaddingValue(numeric_limits<Real>::quiet_NaN());
  gradient.Run();

  // Evaluate harmonic energy
  const double vol  = lattice._dx * lattice._dy * lattice._dz;
  const int    nvox = lattice.NumberOfSpatialPoints();
  double value, harmonic_energy = .0;
  for (int vox = 0; vox < nvox; ++vox) {
    value = .0;
    for (int j = 0; j < lattice._t; ++j) {
      value += squared_gradient_norm(vox + j * nvox);
    }
    value *= vol;
    if (energy) energy->Put(vox, value);
    harmonic_energy += value;
  }

  return harmonic_energy;
}

// -----------------------------------------------------------------------------
/// Evaluate deformation energy of discretized volumetric map
///
/// Li et al. (2009). Meshless harmonic volumetric mapping using fundamental solution methods.
/// IEEE Transactions on Automation Science and Engineering, 6(3), 409–422.
template <class Real>
double EvaluateDeformationEnergy(const GenericImage<Real> &map,
                                 GenericImage<Real>       *energy = NULL,
                                 double lambda = .0335, double mu = .0224)
{
  typedef GradientImageFilter<Real> GradientFilter;

  if (map.T() != 3) {
    cerr << "Error: Can compute deformation energy only from 3D -> 3D volumetric map" << endl;
    exit(1);
  }

  const ImageAttributes &lattice = map.Attributes();
  if (energy) energy->Initialize(lattice, 1);

  // Compute squared norm of gradient vectors for each component
  GenericImage<Real> jac;
  GradientFilter gradient(GradientFilter::GRADIENT_VECTOR);
  gradient.Input(&map);
  gradient.Output(&jac);
  gradient.UseVoxelSize(true);
  gradient.UseOrientation(false); // R = I
  gradient.PaddingValue(numeric_limits<Real>::quiet_NaN());
  gradient.Run();

  // Evaluate deformation energy
  const double vol = lattice._dx * lattice._dy * lattice._dz;

  double stress, value, deformation_energy = .0;
  Vector3D<Real> ga, gb;
  Matrix strain(3, 3);

  for (int k = 0; k < lattice._z; ++k)
  for (int j = 0; j < lattice._y; ++j)
  for (int i = 0; i < lattice._x; ++i) {
    // Compute strain tensor, \epsilon
    for (int b = 0; b < 3; ++b) {
      gb._x = jac(i, j, k, b    ); // dq_x / dp_b
      gb._y = jac(i, j, k, b + 3); // dq_y / dp_b
      gb._z = jac(i, j, k, b + 6); // dq_z / dp_b
      for (int a = 0; a < 3; ++a) {
        ga._x = jac(i, j, k, a    ); // dq_x / dp_a
        ga._y = jac(i, j, k, a + 3); // dq_y / dp_a
        ga._z = jac(i, j, k, a + 6); // dq_z / dp_a
        strain(a, b) = ga.DotProduct(gb);
        if (a == b) strain(a, b) -= 1.0;
      }
    }
    // Compute elastic potential, \eta
    value = .0;
    for (int b = 0; b < 3; ++b)
    for (int a = 0; a < 3; ++a) {
      stress = 2.0 * mu * strain(a, b);
      if (a == b) {
        stress += lambda * strain(0, 0);
        stress += lambda * strain(1, 1);
        stress += lambda * strain(2, 2);
      }
      value += stress * strain(a, b);
    }
    value *= .5;
    // Multiply by cube volume
    value *= vol;
    // Add to total deformation energy
    if (energy) energy->Put(i, j, k, value);
    deformation_energy += value;
  }

  return deformation_energy;
}

// -----------------------------------------------------------------------------
/// Count number of lattice points mapped outside the output domain
template <class Real>
int NumberOfPointsOutside(const GenericImage<Real>     &map,
                          vtkSmartPointer<vtkPointSet>  source,
                          RealImage                    *dfield = NULL)
{
  const int nvox = map.NumberOfSpatialVoxels();

  int    n = 0;
  double p[3], d;

  vtkSmartPointer<vtkPolyData> surface = DataSetSurface(source);

  vtkSmartPointer<vtkImplicitPolyDataDistance> dist;
  dist = vtkSmartPointer<vtkImplicitPolyDataDistance>::New();
  dist->SetInput(surface);

  if (dfield) dfield->Initialize(map.Attributes(), 1);

  const Real *x = map.Data(), *y = x + nvox, *z = y + nvox;
  for (int vox = 0; vox < nvox; ++vox, ++x, ++y, ++z) {
    p[0] = *x, p[1] = *y, p[2] = *z;
    d = dist->EvaluateFunction(p);
    if (dfield) dfield->Put(vox, d);
    if (d > 0) ++n;
  }

  return n;
}

// -----------------------------------------------------------------------------
/// Compute distance of each mapped point to the output domain boundary
vtkSmartPointer<vtkPointSet> DistanceField(const DiscreteMap           *map,
                                           vtkSmartPointer<vtkPointSet> source)
{
  vtkSmartPointer<vtkPolyData> surface = DataSetSurface(source);

  vtkSmartPointer<vtkImplicitPolyDataDistance> dist;
  dist = vtkSmartPointer<vtkImplicitPolyDataDistance>::New();
  dist->SetInput(surface);

  vtkSmartPointer<vtkDataArray> darray = vtkSmartPointer<vtkFloatArray>::New();
  darray->SetName("Distance");
  darray->SetNumberOfComponents(1);
  darray->SetNumberOfTuples(map->Input()->GetNumberOfPoints());

  vtkSmartPointer<vtkPointSet> dfield;
  dfield = vtkSmartPointer<vtkPointSet>::NewInstance(map->Input());
  dfield->ShallowCopy(map->Input());
  dfield->GetPointData()->Initialize();
  dfield->GetPointData()->SetScalars(darray);

  double p[3], d;
  for (vtkIdType ptId = 0; ptId < dfield->GetNumberOfPoints(); ++ptId) {
    map->Values()->GetTuple(ptId, p);
    d = dist->EvaluateFunction(p);
    darray->SetTuple(ptId, &d);
  }

  return dfield;
}

// =============================================================================
// Main
// =============================================================================

// -----------------------------------------------------------------------------
int main(int argc, char *argv[])
{
  // ---------------------------------------------------------------------------
  // Parse arguments
  REQUIRES_POSARGS(0);

  if (argc == 1 || NUM_POSARGS > 2) {
    PrintHelp(EXECNAME);
    exit(1);
  }

  const char *target_name = NULL;
  const char *output_name = NULL;

  if (NUM_POSARGS >= 1) target_name = POSARG(1);
  if (NUM_POSARGS >= 2) output_name = POSARG(2);

  int         volumetric_map          = -1;
  const char *source_name             = NULL;
  const char *boundary_map_name       = NULL;
  const char *boundary_mask_name      = NULL;
  const char *volumetric_map_name     = NULL;
  const char *lattice_name            = NULL;
  double      radius                  = 1.0;
  MapType     output_map_type         = MAP_Harmonic;
  int         no_of_iterations        = 500;
  bool        write_volumetric_map    = false;
  bool        apply_volumetric_map    = false;
  bool        eval_harmonic_energy    = false;
  const char *harmonic_energy_name    = NULL;
  bool        eval_deformation_energy = false;
  const char *deformation_energy_name = NULL;
  bool        eval_outside            = false;
  const char *outside_name            = NULL;
  const char *distance_name           = NULL;

  for (ALL_OPTIONS) {
    if      (OPTION("-surface"))      volumetric_map = false;
    else if (OPTION("-volume" ))      volumetric_map = true;
    else if (OPTION("-array"))        boundary_map_name    = ARGUMENT;
    else if (OPTION("-mask"))         boundary_mask_name   = ARGUMENT;
    else if (OPTION("-map"))          write_volumetric_map = true;
    else if (OPTION("-apply")) volumetric_map_name = ARGUMENT, apply_volumetric_map = true;
    else if (OPTION("-eval"))  volumetric_map_name = ARGUMENT, apply_volumetric_map = false;
    else if (OPTION("-source"))       source_name = ARGUMENT;
    else if (OPTION("-lattice"))      lattice_name  = ARGUMENT;
    else if (OPTION("-radius"))       radius = atof(ARGUMENT);
    else if (OPTION("-harmonic"))     output_map_type = MAP_Harmonic;
    else if (OPTION("-harmonic-fem")) output_map_type = MAP_HarmonicFEM;
    else if (OPTION("-harmonic-mfs")) output_map_type = MAP_HarmonicMFS;
    else if (OPTION("-biharmonic"))   output_map_type = MAP_BiharmonicMFS;
    else if (OPTION("-acap"))         output_map_type = MAP_ACAP;
    else if (OPTION("-iter") || OPTION("-iterations")) {
      const char *arg = ARGUMENT;
      if (!FromString(arg, no_of_iterations)) {
        cerr << "Error: Invalid -iterations argument: " << arg << endl;
        exit(1);
      }
    }
    else if (OPTION("-harmonic-energy")) {
      eval_harmonic_energy = true;
      if (HAS_ARGUMENT) harmonic_energy_name = ARGUMENT;
    }
    else if (OPTION("-deformation-energy")) {
      eval_deformation_energy = true;
      if (HAS_ARGUMENT) deformation_energy_name = ARGUMENT;
    }
    else if (OPTION("-outside")) {
      eval_outside = true;
      if (HAS_ARGUMENT) outside_name = ARGUMENT;
    }
    else if (OPTION("-distance")) {
      distance_name = ARGUMENT;
    }
    else HANDLE_COMMON_OR_UNKNOWN_OPTION();
  }

  if (!target_name && (!volumetric_map_name || apply_volumetric_map)) {
    cerr << "Error: Target point set required unless the -eval option is given" << endl;
    exit(1);
  }

  if (volumetric_map_name && eval_outside && !source_name) {
    cerr << "Error: Source point set required by -outside option" << endl;
    exit(1);
  }

  // Read target point set
  vtkSmartPointer<vtkPointSet> target;
  if (target_name) {
    target = GetPointSet(target_name, radius);
    if (!volumetric_map_name && volumetric_map == -1) {
      volumetric_map = IsTetrahedralMesh(target);
    }
  }

  // Read source point set
  vtkSmartPointer<vtkPointSet> source;
  if (source_name) source = GetPointSet(source_name, radius);

  // ---------------------------------------------------------------------------
  // Apply/evaluate volumetric map
  if (volumetric_map_name) {

    if (volumetric_map == 0) {
      cerr << "Warning: Ignoring -surface option when running in -apply mode" << endl;
    }

    // Read input volumetric map
    if (verbose) cout << "Read volumetric map from " << volumetric_map_name << "...", cout.flush();
    unique_ptr<VolumetricMap> map(VolumetricMap::New(volumetric_map_name));
    DiscreteMap *dmap = dynamic_cast<DiscreteMap *>(map.get());
    if (verbose) cout << " done" << endl;

    if (target && output_name) {
      vtkPoints * const points = target->GetPoints();
      double p[3];

      // Apply volumetric map to points of target point set
      if (apply_volumetric_map) {

        if (verbose) cout << "Apply volumetric map...", cout.flush();
        if (map->NumberOfComponents() != 3) {
          if (verbose) cout << " failed" << endl;
          cerr << "Error: Volumetric map must have output dimension 3" << endl;
          exit(1);
        }
        for (vtkIdType ptId = 0; ptId < points->GetNumberOfPoints(); ++ptId) {
          points->GetPoint(ptId, p);
          map->Evaluate(p, p);
          points->SetPoint(ptId, p);
        }
        if (verbose) cout << " done" << endl;

      // Evaluate volumetric map at points of target point set
      } else {

        if (verbose) cout << "Evaluate volumetric map...", cout.flush();
        vtkSmartPointer<vtkDataArray> discrete_map = vtkSmartPointer<vtkFloatArray>::New();
        discrete_map->SetName("VolumetricMap");
        discrete_map->SetNumberOfComponents(map->NumberOfComponents());
        discrete_map->SetNumberOfTuples(target->GetNumberOfPoints());
        double *v = new double[map->NumberOfComponents()];
        for (vtkIdType ptId = 0; ptId < target->GetNumberOfPoints(); ++ptId) {
          target->GetPoint(ptId, p);
          if (!map->Evaluate(v, p)) {
            cerr << "Warning: Volumetric map undefined at point ("
                 << p[0] << ", " << p[1] << ", " << p[2]
                 << ") with ID " << ptId << endl;
          }
          discrete_map->SetTuple(ptId, v);
        }
        delete[] v;
        if (discrete_map->GetNumberOfComponents() == 1) {
          target->GetPointData()->SetScalars(discrete_map);
        } else if (discrete_map->GetNumberOfComponents() == 3) {
          target->GetPointData()->SetTCoords(discrete_map);
        } else {
          target->GetPointData()->AddArray(discrete_map);
        }
        if (verbose) cout << " done" << endl;

      }

      // Write output point set
      if (!WritePointSet(output_name, target)) {
        cerr << "Error: Failed to write point set to " << output_name << endl;
        exit(1);
      }
    }

    // Evaluate volumetric map at lattice points
    RealImage discrete_map;
    if (eval_harmonic_energy || eval_deformation_energy || eval_outside || (distance_name && !dmap)) {
      if (verbose) cout << "Discretize volumetric map", cout.flush();
      if (lattice_name) {
        discrete_map.Read(lattice_name);
        discrete_map.Initialize(discrete_map.Attributes(), map->NumberOfComponents());
      } else {
        discrete_map.Initialize(map->Attributes(128, 128, 128), map->NumberOfComponents());
      }
      if (verbose > 1) {
        cout << " (N = " << discrete_map.X()
             <<    " x " << discrete_map.Y()
             <<    " x " << discrete_map.Z() << ")";
      }
      cout << "...";
      cout.flush();
      map->Evaluate(discrete_map, 0, target);
      if (debug) discrete_map.Write("discrete_map.nii.gz");
      if (verbose) cout << " done" << endl;
    }

    // Evaluate distance of mapped points to output domain boundary
    if (distance_name) {
      if (dmap) {
        vtkSmartPointer<vtkPointSet> dfield = DistanceField(dmap, source);
        if (!WritePointSet(distance_name, dfield)) {
          cerr << "Error: Failed to write mapped distance field to " << distance_name << endl;
          exit(1);
        }
      } else {
        // TODO: Modify NumberOfPointsOutside function
        cerr << "Error: Can compute -distance field currently only for tetrahedral volumetric maps" << endl;
        exit(1);
      }
    }

    double rms_boundary_error, min_boundary_error, max_boundary_error;
    double harmonic_energy;
    double deformation_energy;
    int    noutside;

    // Evaluate error of volumetric map at boundary
    if (boundary_map_name) {
      if (verbose) {
        cout << "Evaluate boundary map error...";
        if (verbose > 1) cout << "\n";
        cout.flush();
      }
      vtkDataArray *boundary_map = GetBoundaryMap(target, boundary_map_name);
      if (map->NumberOfComponents() != boundary_map->GetNumberOfComponents()) {
        cerr << "Error: Boundary map and volumetric map have differing output domoin dimension" << endl;
        exit(1);
      }
      int    num = 0;
      double p[3], error, sum = .0;
      min_boundary_error = +numeric_limits<double>::infinity();
      max_boundary_error = -numeric_limits<double>::infinity();
      double *b = new double[map->NumberOfComponents()];
      double *v = new double[map->NumberOfComponents()];
      for (vtkIdType ptId = 0; ptId < target->GetNumberOfPoints(); ++ptId) {
        target->GetPoint(ptId, p);
        if (!map->Evaluate(v, p)) {
          cerr << "Warning: Volumetric map undefined at point ("
                 << p[0] << ", " << p[1] << ", " << p[2]
                 << ") with ID " << ptId << endl;
          continue;
        }
        error = .0;
        for (int j = 0; j < map->NumberOfComponents(); ++j) {
          error += pow(v[j] - boundary_map->GetComponent(ptId, j), 2);
        }
        if (verbose > 1) {
          cout << "Point " << setw(7) << (ptId + 1)
               << ": Boundary map error = " << sqrt(error) << endl;
        }
        if (error < min_boundary_error) min_boundary_error = error;
        if (error > max_boundary_error) max_boundary_error = error;
        sum += error;
      }
      min_boundary_error = sqrt(min_boundary_error);
      max_boundary_error = sqrt(max_boundary_error);
      rms_boundary_error = sqrt(sum / num);
      delete[] v;
      delete[] b;
      if (verbose) {
        if (verbose > 1) cout << "Evaluate boundary map error...";
        cout << " done";
        if (verbose > 1) cout << "\n";
        cout << endl;
      }
    }

    // Determine how many points are mapped outside the output domain
    if (eval_outside || (distance_name && !dmap)) {
      if (verbose) cout << "Evaluate distance of mapped points to output domain boundary...", cout.flush();
      if (outside_name || distance_name) {
        RealImage dfield;
        noutside = NumberOfPointsOutside(discrete_map, source, &dfield);
        if (outside_name ) dfield.Write(outside_name);
        if (distance_name) dfield.Write(distance_name);
      } else {
        noutside = NumberOfPointsOutside(discrete_map, source);
      }
      if (verbose) cout << " done" << endl;
    }

    // Normalize models to unit box before evaluation energy measures
    if ((eval_harmonic_energy || eval_deformation_energy) && target && source) {
      NormalizeCoordinates(discrete_map, target, source);
    }

    // Evaluate harmonic energy of volumetric map
    if (eval_harmonic_energy) {
      if (verbose) cout << "Evaluate harmonic energy...", cout.flush();
      if (harmonic_energy_name) {
        RealImage energy;
        harmonic_energy = EvaluateHarmonicEnergy(discrete_map, &energy);
        energy.Write(harmonic_energy_name);
      } else {
        harmonic_energy = EvaluateHarmonicEnergy(discrete_map);
      }
      if (verbose) cout << " done" << endl;
    }

    // Evaluate deformation energy of volumetric map
    if (eval_deformation_energy) {
      if (verbose) cout << "Evaluate deformation energy...", cout.flush();
      if (deformation_energy_name) {
        RealImage energy;
        deformation_energy = EvaluateDeformationEnergy(discrete_map, &energy);
        energy.Write(deformation_energy_name);
      } else {
        deformation_energy = EvaluateDeformationEnergy(discrete_map);
      }
      if (verbose) cout << " done" << endl;
    }

    if (boundary_map_name) {
      cout << "Boundary RMS error     = " << rms_boundary_error << endl;
      cout << "Minimum boundary error = " << min_boundary_error << endl;
      cout << "Maximum boundary error = " << max_boundary_error << endl;
    }
    if (eval_harmonic_energy   ) cout << "Harmonic energy        = " << harmonic_energy << endl;
    if (eval_deformation_energy) cout << "Deformation energy     = " << deformation_energy << endl;
    if (eval_outside           ) cout << "No. of points outside  = " << noutside << endl;

  // ---------------------------------------------------------------------------
  // Compute volumetric map
  } else if (volumetric_map) {

    vtkSmartPointer<vtkDataArray> boundary_map, boundary_mask;

    // Read or compute boundary map
    if (boundary_map_name) {
      boundary_map = GetBoundaryMap(target, boundary_map_name);
    } else {
      boundary_map = target->GetPointData()->GetTCoords();
    }
    if (source) {
      if (boundary_map) {
        cerr << "Warning: Ignoring -source point set, using input -map instead" << endl;
      } else {
        if (verbose) cout << "Computing boundary map...", cout.flush();
        boundary_map = BoundaryMap(target, source, radius, output_map_type);
        if (!boundary_map) {
          if (verbose) cout << " failed" << endl;
          cerr << "Error: Failed to compute boundary map" << endl;
          exit(1);
        }
        if (verbose) cout << " done" << endl;
      }
    }
    // Get boundary mask
    if (boundary_mask_name) {
      boundary_mask = target->GetPointData()->GetArray(boundary_mask_name);
      if (!boundary_mask) {
        cerr << "Error: Target point set has no point data array named " << boundary_mask_name << endl;
        exit(1);
      }
    }
    // Compute volumetric map given boundary map
    if (boundary_map) {
      unique_ptr<VolumetricMap> map;
      map.reset(ComputeVolumetricMap(target, boundary_map, boundary_mask, output_map_type, no_of_iterations));
      if (write_volumetric_map) {
        if (!map->Write(output_name)) {
          cerr << "Error: Failed to write volumetric map to " << output_name << endl;
          exit(1);
        }
      } else {
        if (verbose) cout << "Evaluating volumetric map...", cout.flush();
        map->Initialize();
        vtkSmartPointer<vtkDataArray> discrete_map = vtkSmartPointer<vtkFloatArray>::New();
        discrete_map->SetName("VolumetricMap");
        discrete_map->SetNumberOfComponents(map->NumberOfComponents());
        discrete_map->SetNumberOfTuples(target->GetNumberOfPoints());
        double p[3], *v = new double[map->NumberOfComponents()];
        for (vtkIdType ptId = 0; ptId < target->GetNumberOfPoints(); ++ptId) {
          target->GetPoint(ptId, p);
          if (!map->Evaluate(v, p)) {
            cerr << "Warning: Volumetric map undefined at point ("
                 << p[0] << ", " << p[1] << ", " << p[2]
                 << ") with ID " << ptId << endl;
          }
          discrete_map->SetTuple(ptId, v);
        }
        delete[] v;
        target->GetPointData()->SetTCoords(discrete_map);
        if (verbose) cout << " done" << endl;
        if (!WritePointSet(output_name, target)) {
          cerr << "Error: Failed to write point set with volumetric map to " << output_name << endl;
          exit(1);
        }
      }
    // Tesselate interior of target point set
    } else {
      if (verbose) cout << "Computing tetrahedral mesh...", cout.flush();
      vtkSmartPointer<vtkPointSet> mesh = Tetrahedralize(target);
      if (!WritePointSet(output_name, mesh)) {
        if (verbose) cout << " failed" << endl;
        cerr << "Error: Failed to write tetrahedral mesh to " << output_name << endl;
        exit(1);
      }
      if (verbose) cout << " done" << endl;
    }

  // ---------------------------------------------------------------------------
  // Compute surface map
  } else {

    // Map target surface to source surface
    if (source) {

      if (boundary_map_name) {
        cerr << "Warning: Ignoring input -map, computing a new boundary map instead" << endl;
      }
      if (verbose) cout << "Computing boundary map...", cout.flush();
      vtkSmartPointer<vtkDataArray> boundary_map;
      boundary_map = BoundaryMap(target, source, radius, output_map_type);
      if (!boundary_map) {
        if (verbose) cout << " failed" << endl;
        cerr << "Error: Failed to compute boundary map" << endl;
        exit(1);
      }
      target->GetPointData()->AddArray(boundary_map);
      if (!WritePointSet(output_name, target)) {
        if (verbose) cout << " failed" << endl;
        cerr << "Error: Failed to write target point set with computed boundary map to " << output_name << endl;
        exit(1);
      }
      if (verbose) cout << " done" << endl;

    // Tesselate surface
    } else {
      if (boundary_map_name) {
        cerr << "Warning: Ignoring input -map, tesselating target surface instead" << endl;
      }
      if (verbose) cout << "Computing surface mesh...", cout.flush();
      vtkSmartPointer<vtkPolyData> mesh = DataSetSurface(target);
      if (!WritePolyData(output_name, mesh)) {
        if (verbose) cout << " failed" << endl;
        cerr << "Error: Failed to write surface mesh to " << output_name << endl;
        exit(1);
      }
      if (verbose) cout << " done" << endl;
    }

  }

  return 0;
}
