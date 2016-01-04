/*
 * Medical Image Registration ToolKit (MMIRTK)
 *
 * Copyright 2013-2015 Imperial College London
 * Copyright 2013-2015 Andreas Schuh
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#include <mirtkVolumetricMap.h>

#include <mirtkMath.h>
#include <mirtkCfstream.h>
#include <mirtkBaseImage.h>
#include <mirtkVoxelFunction.h>
#include <mirtkDiscreteMap.h>
#include <mirtkHarmonicMap.h>
#include <mirtkBiharmonicMap.h>
#include <mirtkPointSetUtils.h>

#include <vtkSmartPointer.h>
#include <vtkPointSet.h>
#include <vtkPolyData.h>
#include <vtkImageData.h>


namespace mirtk {


// =============================================================================
// Auxiliary functors
// =============================================================================

namespace VolumetricMapUtils {


// -----------------------------------------------------------------------------
/// Evaluate volumetric map at lattice points
class EvaluateMap : public VoxelFunction
{
  const VolumetricMap *_Map;
  vtkImageData        *_Domain;
  const BaseImage     *_Output;
  const int            _NumberOfVoxels;
  const int            _l1, _l2;

public:

  EvaluateMap(const VolumetricMap *map,
              vtkImageData        *domain,
              const BaseImage     *output,
              int l1, int l2)
  :
    _Map(map),
    _Domain(domain),
    _Output(output),
    _NumberOfVoxels(output->NumberOfSpatialVoxels()),
    _l1(l1), _l2(l2)
  {}

  template <class T>
  void operator ()(int i, int j, int k, int, T *v) const
  {
    if (!_Domain || _Domain->GetScalarComponentAsFloat(i, j, k, 0) != .0) {
      double x = i, y = j, z = k;
      _Output->ImageToWorld(x, y, z);
      double *f = new double[_Map->NumberOfComponents()];
      _Map->Evaluate(f, x, y, z);
      for (int l = _l1; l < _l2; ++l, v += _NumberOfVoxels) {
        *v = f[l];
      }
      delete[] f;
    } else {
      for (int l = _l1; l < _l2; ++l, v += _NumberOfVoxels) {
        *v = numeric_limits<T>::quiet_NaN();
      }
    }
  }
};


} // namespace VolumetricMapUtils
using namespace VolumetricMapUtils;

// =============================================================================
// Construction/Destruction
// =============================================================================

// -----------------------------------------------------------------------------
VolumetricMap *VolumetricMap::New(const char *fname)
{
  VolumetricMap *map = NULL;

  const size_t max_name_len = 32;
  char         map_type_name[max_name_len];

  Cifstream is(fname);
  is.ReadAsChar(map_type_name, max_name_len);

  if (strncmp(map_type_name, HarmonicMap::NameOfType(), max_name_len) == 0) {
    map = new HarmonicMap();
  } else if (strncmp(map_type_name, BiharmonicMap::NameOfType(), max_name_len) == 0) {
    map = new BiharmonicMap();
  } else {
    map = new DiscreteMap();
  }

  map->Read(fname);
  return map;
}

// -----------------------------------------------------------------------------
void VolumetricMap::CopyAttributes(const VolumetricMap &other)
{
  _OutsideValue = other._OutsideValue;
}

// -----------------------------------------------------------------------------
VolumetricMap::VolumetricMap()
:
  _OutsideValue(numeric_limits<double>::quiet_NaN())
{
}

// -----------------------------------------------------------------------------
VolumetricMap::VolumetricMap(const VolumetricMap &other)
:
  Object(other)
{
  CopyAttributes(other);
}

// -----------------------------------------------------------------------------
VolumetricMap &VolumetricMap::operator =(const VolumetricMap &other)
{
  if (this != &other) {
    Object::operator =(other);
    CopyAttributes(other);
  }
  return *this;
}

// -----------------------------------------------------------------------------
void VolumetricMap::Initialize()
{
}

// -----------------------------------------------------------------------------
VolumetricMap::~VolumetricMap()
{
}

// =============================================================================
// Input domain
// =============================================================================

// -----------------------------------------------------------------------------
ImageAttributes VolumetricMap::Attributes(int nx, int ny, int nz) const
{
  if (ny <= 0) ny = nx;
  if (nz <= 0) nz = nx;

  double x1, y1, z1, x2, y2, z2;
  this->BoundingBox(x1, y1, z1, x2, y2, z2);
  const double lx = x2 - x1;
  const double ly = y2 - y1;
  const double lz = z2 - z1;

  ImageAttributes lattice;
  lattice._xorigin = x1 + .5 * lx;
  lattice._yorigin = y1 + .5 * ly;
  lattice._zorigin = z1 + .5 * lz;
  lattice._x       = nx;
  lattice._y       = ny;
  lattice._z       = nz;
  lattice._dx      = (lattice._x == 0 ? .0 : lx / nx);
  lattice._dy      = (lattice._y == 0 ? .0 : ly / ny);
  lattice._dz      = (lattice._z == 0 ? .0 : lz / nz);
  if (lattice._x == 0) lattice._x = 1;
  if (lattice._y == 0) lattice._y = 1;
  if (lattice._z == 0) lattice._z = 1;

  return lattice;
}

// -----------------------------------------------------------------------------
ImageAttributes VolumetricMap::Attributes(double dx, double dy, double dz) const
{
  double x1, y1, z1, x2, y2, z2;
  this->BoundingBox(x1, y1, z1, x2, y2, z2);
  const double lx = x2 - x1;
  const double ly = y2 - y1;
  const double lz = z2 - z1;

  if (dx <= .0) dx = sqrt(lx*lx + ly*ly + lz*lz) / 256;
  if (dy <= .0) dy = dx;
  if (dz <= .0) dz = dx;

  ImageAttributes lattice;
  lattice._xorigin = x1 + .5 * lx;
  lattice._yorigin = y1 + .5 * ly;
  lattice._zorigin = z1 + .5 * lz;
  lattice._x       = iceil(lx / dx);
  lattice._y       = iceil(ly / dy);
  lattice._z       = iceil(lz / dz);
  lattice._dx      = (lattice._x == 0 ? .0 : dx);
  lattice._dy      = (lattice._y == 0 ? .0 : dy);
  lattice._dz      = (lattice._z == 0 ? .0 : dz);
  if (lattice._x == 0) lattice._x = 1;
  if (lattice._y == 0) lattice._y = 1;
  if (lattice._z == 0) lattice._z = 1;

  return lattice;
}

// =============================================================================
// Evaluation
// =============================================================================

// -----------------------------------------------------------------------------
void VolumetricMap::Evaluate(GenericImage<float> &f, int l, vtkSmartPointer<vtkPointSet> domain) const
{
  ImageAttributes lattice = f.Attributes();
  lattice._dt = .0;

  if (l >= NumberOfComponents() || l + lattice._t > NumberOfComponents()) {
    cerr << "VolumetricMap::Evaluate: Component index out of range" << endl;
    exit(1);
  }

  vtkSmartPointer<vtkImageData> mask;
  if (domain) {
    mask = NewVtkMask(lattice._x, lattice._y, lattice._z);
    ImageStencilToMask(ImageStencil(mask, WorldToImage(domain, &f)), mask);
  }

  ParallelForEachVoxel(EvaluateMap(this, mask, &f, l, l + lattice._t), lattice, f);
}

// -----------------------------------------------------------------------------
void VolumetricMap::Evaluate(GenericImage<double> &f, int l, vtkSmartPointer<vtkPointSet> domain) const
{
  ImageAttributes lattice = f.Attributes();
  lattice._dt = .0;

  if (l >= NumberOfComponents() || l + lattice._t > NumberOfComponents()) {
    cerr << "VolumetricMap::Evaluate: Component index out of range" << endl;
    exit(1);
  }

  vtkSmartPointer<vtkImageData> mask;
  if (domain) {
    mask = NewVtkMask(lattice._x, lattice._y, lattice._z);
    ImageStencilToMask(ImageStencil(mask, WorldToImage(domain, &f)), mask);
  }

  ParallelForEachVoxel(EvaluateMap(this, mask, &f, l, l + lattice._t), lattice, f);
}

// =============================================================================
// I/O
// =============================================================================

// -----------------------------------------------------------------------------
bool VolumetricMap::Read(const char *fname)
{
  Cifstream is(fname);
  const size_t max_name_len = 32;
  char         map_type_name[max_name_len];
  is.ReadAsChar(map_type_name, max_name_len);
  if (strncmp(map_type_name, this->NameOfClass(), max_name_len) != 0) {
    return false;
  }
  this->ReadMap(is);
  this->Initialize();
  return true;
}

// -----------------------------------------------------------------------------
bool VolumetricMap::Write(const char *fname) const
{
  Cofstream os(fname);
  const size_t max_name_len = 32;
  char         map_type_name[max_name_len] = {0};
  strncpy(map_type_name, this->NameOfClass(), max_name_len);
  os.WriteAsChar(map_type_name, max_name_len);
  this->WriteMap(os);
  return true;
}

// -----------------------------------------------------------------------------
void VolumetricMap::ReadMap(Cifstream &)
{
  cerr << this->NameOfClass() << "::ReadMap not implemented" << endl;
  exit(1);
}

// -----------------------------------------------------------------------------
void VolumetricMap::WriteMap(Cofstream &) const
{
  cerr << this->NameOfClass() << "::WriteMap not implemented" << endl;
  exit(1);
}


} // namespace mirtk
