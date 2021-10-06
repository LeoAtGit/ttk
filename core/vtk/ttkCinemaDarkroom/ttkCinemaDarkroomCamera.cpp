#include <ttkCinemaDarkroomCamera.h>

#include <vtkCamera.h>
#include <vtkDoubleArray.h>
#include <vtkInformation.h>
#include <vtkIntArray.h>
#include <vtkObjectFactory.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>

#include <vtkCompositeDataSet.h>

#include <vtkPythonInterpreter.h>

#include <ttkUtils.h>

const int nLinesPerFrustum{16};

template <typename T, int DataType>
T *prepArray(vtkFieldData *fieldData,
             const std::string name,
             const int nTuples,
             const int nComponents) {
  auto array = vtkSmartPointer<vtkDataArray>::Take(
    vtkDataArray::CreateDataArray(DataType));
  array->SetName(name.data());
  array->SetNumberOfComponents(nComponents);
  array->SetNumberOfTuples(nTuples);

  fieldData->AddArray(array);

  return ttkUtils::GetPointer<T>(array);
}

vtkStandardNewMacro(ttkCinemaDarkroomCamera);

ttkCinemaDarkroomCamera::ttkCinemaDarkroomCamera() {
  this->setDebugMsgPrefix("CinemaDarkroomCamera");

  this->SetNumberOfInputPorts(3);
  this->SetNumberOfOutputPorts(2);
}

ttkCinemaDarkroomCamera::~ttkCinemaDarkroomCamera() {
}

int ttkCinemaDarkroomCamera::FillInputPortInformation(int port,
                                                      vtkInformation *info) {
  if(port >= 0 && port < this->GetNumberOfInputPorts()) {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkDataObject");
    info->Set(vtkAlgorithm::INPUT_IS_OPTIONAL(), 1);
    return 1;
  }
  return 0;
}

int ttkCinemaDarkroomCamera::FillOutputPortInformation(int port,
                                                       vtkInformation *info) {
  if(port == 0 || port == 1) {
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkPolyData");
    return 1;
  }
  return 0;
}

int ttkCinemaDarkroomCamera::SyncParaViewCamera() {
  ttk::Timer timer;
  const std::string msg{"Updating ParaView Camera Parameters"};
  this->printMsg(msg, 0, 0, 1, ttk::debug::LineMode::REPLACE);

  std::string code(R"(
from paraview.simple import GetActiveSource
from paraview.simple import GetActiveCamera

drCamera = GetActiveSource()
pvCamera = GetActiveCamera()

if not drCamera or drCamera.__class__.__name__!='TTKDarkroomCamera':
  print("")
  print("ERROR: Active Source must be a 'TTKDarkroomCamera'")
  raise IOError

if not pvCamera:
  print("")
  print("ERROR: No active PV Camera")
  raise IOError

pvCamera.SetPosition(drCamera.Position)
pvCamera.SetViewUp(drCamera.Up)
pvCamera.SetFocalPoint([
  drCamera.Position[0] + drCamera.Direction[0],
  drCamera.Position[1] + drCamera.Direction[1],
  drCamera.Position[2] + drCamera.Direction[2]
])
pvCamera.SetParallelProjection(drCamera.Orthographic)
pvCamera.SetViewAngle(drCamera.Angle)
pvCamera.SetParallelScale(drCamera.Scale)

)");

  vtkPythonInterpreter::RunSimpleString(code.data());
  this->printMsg(msg, 1, timer.getElapsedTime(), 1);
  return 1;
}

int ttkCinemaDarkroomCamera::SyncDarkroomCamera() {
  ttk::Timer timer;
  const std::string msg{"Updating Darkroom Camera Parameters"};
  this->printMsg(msg, 0, 0, 1, ttk::debug::LineMode::REPLACE);

  std::string code(R"(
from paraview.simple import GetActiveSource
from paraview.simple import GetActiveCamera
from paraview.simple import GetActiveView

drCamera = GetActiveSource()
pvCamera = GetActiveCamera()
pvView = GetActiveView()

if not drCamera or drCamera.__class__.__name__!='TTKDarkroomCamera':
  print("")
  print("ERROR: Active Source must be a 'TTKDarkroomCamera'")
  raise IOError

if not pvCamera or not pvView:
  print("")
  print("ERROR: No active PV Camera and/or View")
  raise IOError

drCamera.Resolution = pvView.ViewSize

drCamera.Position = pvCamera.GetPosition()

fp = pvCamera.GetFocalPoint()

drCamera.Direction = [
  fp[0]-drCamera.Position[0],
  fp[1]-drCamera.Position[1],
  fp[2]-drCamera.Position[2]
]
drCamera.Up = pvCamera.GetViewUp()
drCamera.Orthographic = pvCamera.GetParallelProjection()
drCamera.Angle = pvCamera.GetViewAngle()
drCamera.Scale = pvCamera.GetParallelScale()

)");

  vtkPythonInterpreter::RunSimpleString(code.data());
  this->printMsg(msg, 1, timer.getElapsedTime(), 1);
  return 1;
}

ttkCinemaDarkroomCamera::Calibration::Calibration(vtkPointSet *cameras) {
  auto camerasPD = cameras->GetPointData();

  this->Resolution.Array = camerasPD->GetArray("CamResolution");
  this->Position.Array = cameras->GetPoints()->GetData();
  this->Direction.Array = camerasPD->GetArray("CamDirection");
  this->Up.Array = camerasPD->GetArray("CamUp");
  this->NearFar.Array = camerasPD->GetArray("CamNearFar");
  this->Scale.Array = camerasPD->GetArray("CamScale");
  this->Angle.Array = camerasPD->GetArray("CamAngle");

  this->Orthographic = this->Scale.Array;
  bool valid = this->Resolution.Array
               && this->Resolution.Array->GetNumberOfComponents() == 2
               && this->Direction.Array
               && this->Direction.Array->GetNumberOfComponents() == 3
               && this->Up.Array && this->Up.Array->GetNumberOfComponents() == 3
               && this->NearFar.Array
               && this->NearFar.Array->GetNumberOfComponents() == 2
               && ((this->Orthographic
                    && this->Scale.Array->GetNumberOfComponents() == 1)
                   || (!this->Orthographic && this->Angle.Array
                       && this->Angle.Array->GetNumberOfComponents() == 1));

  if(valid) {
    this->nCameras = cameras->GetNumberOfPoints();
  }
}

std::vector<double>
  ttkCinemaDarkroomCamera::Calibration::Parameter::Get(const int idx) {
  std::vector<double> a(this->Array->GetNumberOfComponents());
  this->Array->GetTuple(idx, a.data());
  return a;
}

template <typename T0, typename T1, int n>
void copy(T0 *des, const T1 *src) {
  for(int i = 0; i < n; i++)
    des[i] = (T0)src[i];
}

double length(const double vec[3]) {
  return std::sqrt(vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2]);
};

void add(double a[3], const double b[3]) {
  a[0] += b[0];
  a[1] += b[1];
  a[2] += b[2];
};

void multiply(double vec[3], const double s) {
  vec[0] *= s;
  vec[1] *= s;
  vec[2] *= s;
};

void norm(double vec[3]) {
  multiply(vec, 1.0 / length(vec));
};

void cross(double c[3], const double a[3], const double b[3]) {
  c[0] = a[1] * b[2] - a[2] * b[1];
  c[1] = a[2] * b[0] - a[0] * b[2];
  c[2] = a[0] * b[1] - a[1] * b[0];
};

int computeBoundingBox(double bounds[6],
                       double center[3],
                       double &diameter,
                       vtkDataObject *obj) {
  if(obj->IsA("vtkCompositeDataSet"))
    static_cast<vtkCompositeDataSet *>(obj)->GetBounds(bounds);
  else if(obj->IsA("vtkDataSet"))
    static_cast<vtkDataSet *>(obj)->GetBounds(bounds);
  else
    return 0;

  const double d[3]{
    bounds[1] - bounds[0], bounds[3] - bounds[2], bounds[5] - bounds[4]};

  diameter = length(d);
  center[0] = bounds[0] + 0.5 * d[0];
  center[1] = bounds[2] + 0.5 * d[1];
  center[2] = bounds[4] + 0.5 * d[2];

  return 1;
}

int addPlaneCorners(vtkPoints *points,
                    const double distance,
                    const bool orthographic,
                    const double angleOrScale,
                    const double *position,
                    const double *dir,
                    const double *up,
                    const double *resolution) {
  double dv = orthographic
                ? angleOrScale
                : std::tan(angleOrScale / 2.0 / 180.0 * 3.1415) * distance;

  double u[3]{up[0], up[1], up[2]};
  norm(u);
  multiply(u, dv);

  double r[3]{0, 0, 0};
  cross(r, u, dir);
  norm(r);
  multiply(r, dv * resolution[0] / resolution[1]);

  // move to center of plane
  double t[3]{0, 0, 0};
  add(t, dir);
  norm(t);
  multiply(t, distance);
  add(t, position);

  // top right
  add(t, u);
  add(t, r);
  points->InsertNextPoint(t);

  // bottom right
  multiply(u, -2.0);
  add(t, u);
  points->InsertNextPoint(t);

  // bottom left
  multiply(r, -2.0);
  add(t, r);
  points->InsertNextPoint(t);

  // top left
  multiply(u, -1.0);
  add(t, u);
  points->InsertNextPoint(t);

  return 1;
}

int ttkCinemaDarkroomCamera::ComputeCalibrations(
  vtkPolyData *cameras,
  vtkPointSet *positions,
  vtkDataObject *nearFarBoundingBoxObject,
  vtkDataObject *focusObject) {
  ttk::Timer timer;
  const std::string msg{"Computing Calibrations"};
  this->printMsg(msg, 0, 0, 1, ttk::debug::LineMode::REPLACE);

  // points
  {
    auto cameraPoints = vtkSmartPointer<vtkPoints>::New();
    cameraPoints->SetDataTypeToDouble();

    if(positions) {
      auto inputPoints = positions->GetPoints();
      if(inputPoints->GetDataType() == VTK_DOUBLE) {
        cameraPoints->ShallowCopy(inputPoints);
      } else {
        auto inputPositionsData
          = ttkUtils::GetPointer<float>(inputPoints->GetData());
        const int n = inputPoints->GetNumberOfPoints();
        cameraPoints->SetNumberOfPoints(n);
        auto outputPositionsData
          = ttkUtils::GetPointer<double>(cameraPoints->GetData());

        for(int i = 0; i < n * 3; i++)
          outputPositionsData[i] = (double)inputPositionsData[i];
      }
    } else {
      cameraPoints->InsertNextPoint(
        this->Position[0], this->Position[1], this->Position[2]);
    }
    cameras->SetPoints(cameraPoints);
  }
  const int nPoints = cameras->GetNumberOfPoints();

  // cells
  {
    auto cells = vtkSmartPointer<vtkCellArray>::New();
    auto offsetArray = vtkSmartPointer<vtkIntArray>::New();
    {
      offsetArray->SetNumberOfTuples(nPoints + 1);
      auto data = ttkUtils::GetPointer<int>(offsetArray);
      for(int i = 0; i <= nPoints; i++)
        data[i] = i;
    }

    auto connectivityArray = vtkSmartPointer<vtkIntArray>::New();
    {
      connectivityArray->SetNumberOfTuples(nPoints);
      auto data = ttkUtils::GetPointer<int>(connectivityArray);
      for(int i = 0; i < nPoints; i++)
        data[i] = i;
    }
    cells->SetData(offsetArray, connectivityArray);
    cameras->SetVerts(cells);
  }

  auto cameraPD = cameras->GetPointData();

  // compute point data arrays
  auto resolutionData
    = prepArray<int, VTK_INT>(cameraPD, "CamResolution", nPoints, 2);

  auto positionArray = cameras->GetPoints()->GetData();
  auto directionData
    = prepArray<double, VTK_DOUBLE>(cameraPD, "CamDirection", nPoints, 3);
  auto upData = prepArray<double, VTK_DOUBLE>(cameraPD, "CamUp", nPoints, 3);
  auto nearFarData
    = prepArray<double, VTK_DOUBLE>(cameraPD, "CamNearFar", nPoints, 2);
  auto scaleData = prepArray<double, VTK_DOUBLE>(
    cameraPD, this->Orthographic ? "CamScale" : "CamAngle", nPoints, 1);

  double nearFarBBCenter[3];
  double nearFarBBDiameter;
  if(nearFarBoundingBoxObject) {
    double bounds[6];
    if(!computeBoundingBox(
         bounds, nearFarBBCenter, nearFarBBDiameter, nearFarBoundingBoxObject))
      return !this->printErr(
        "Automatic Parameter Tuning requires 'vtkCompositeDataSet' or "
        "'vtkDataSet' input data obejcts.");
  }

  double focusObjectCenter[3];
  if(focusObject) {
    double bounds[6];
    double diameter{0.0};
    if(!computeBoundingBox(bounds, focusObjectCenter, diameter, focusObject))
      return !this->printErr(
        "Automatic Parameter Tuning requires 'vtkCompositeDataSet' or "
        "'vtkDataSet' input data obejcts.");
  }

  for(int p = 0; p < nPoints; p++) {
    copy<int, int, 2>(&resolutionData[p * 2], this->Resolution);

    scaleData[p] = this->Orthographic ? this->Scale : this->Angle;

    double position[3];
    positionArray->GetTuple(p, position);

    if(nearFarBoundingBoxObject) {
      const double dirToNFBBC[3]{nearFarBBCenter[0] - position[0],
                                 nearFarBBCenter[1] - position[1],
                                 nearFarBBCenter[2] - position[2]};
      const double distanceToNFBBC = length(dirToNFBBC);

      double nearFar[2];
      nearFar[0] = std::max(0.01, distanceToNFBBC - nearFarBBDiameter / 2.0);
      nearFar[1] = nearFar[0] + nearFarBBDiameter;
      copy<double, double, 2>(&nearFarData[p * 2], nearFar);
    } else {
      copy<double, double, 2>(&nearFarData[p * 2], this->NearFar);
    }

    if(focusObject) {
      double direction[3]{0, 0, 0};

      add(direction, position);
      multiply(direction, -1.0);
      add(direction, focusObjectCenter);
      norm(direction);
      copy<double, double, 3>(&directionData[p * 3], direction);
    } else {
      copy<double, double, 3>(&directionData[p * 3], this->Direction);
    }
    norm(&directionData[p * 3]);

    // make sure up vector is at right angle to dir
    double right[3]{0, 0, 0};
    cross(right, &directionData[p * 3], this->Up);
    cross(&upData[p * 3], right, &directionData[p * 3]);
    norm(&upData[p * 3]);
  }

  this->printMsg(msg, 1, timer.getElapsedTime(), 1);

  return 1;
}

int ttkCinemaDarkroomCamera::ComputeFrustums(vtkPointSet *cameras,
                                             vtkPolyData *frustums) {
  ttk::Timer timer;
  const std::string msg{"Computing Frustums"};
  this->printMsg(msg, 0, 0, 1, ttk::debug::LineMode::REPLACE);

  const int nEdgesPerFrustum = 16;
  const int nPointsPerFrustum = 9;

  Calibration calibrations(cameras);
  if(calibrations.nCameras < 0)
    return !this->printErr("Invalid Camera Calibrations.");

  // points
  auto points = vtkSmartPointer<vtkPoints>::New();
  frustums->SetPoints(points);

  // points
  for(int c = 0; c < calibrations.nCameras; c++) {
    points->InsertNextPoint(calibrations.Position.Get(c).data());
    for(int nf = 0; nf < 2; nf++) {
      addPlaneCorners(
        points, calibrations.NearFar.Get(c)[nf], calibrations.Orthographic,
        (calibrations.Orthographic ? calibrations.Scale.Get(c)
                                   : calibrations.Angle.Get(c))[0],
        calibrations.Position.Get(c).data(),
        calibrations.Direction.Get(c).data(), calibrations.Up.Get(c).data(),
        calibrations.Resolution.Get(c).data());
    }
  }

  auto frustumId
    = prepArray<int, VTK_INT>(frustums->GetPointData(), "FrustumId",
                              calibrations.nCameras * nPointsPerFrustum, 1);
  for(int c = 0, q = 0; c < calibrations.nCameras; c++) {
    for(int i = 0; i < nPointsPerFrustum; i++)
      frustumId[q++] = c;
  }

  // cells
  auto cells = vtkSmartPointer<vtkCellArray>::New();
  auto offsetArray = vtkSmartPointer<vtkIntArray>::New();
  auto connectivityArray = vtkSmartPointer<vtkIntArray>::New();

  const int nEdges = calibrations.nCameras * nEdgesPerFrustum;

  // offsets
  {
    offsetArray->SetNumberOfTuples(nEdges + 1);
    auto data = ttkUtils::GetPointer<int>(offsetArray);
    for(int i = 0; i <= nEdges; i++)
      data[i] = i * 2;
  }

  {
    connectivityArray->SetNumberOfTuples(2 * nEdges);
    auto data = ttkUtils::GetPointer<int>(connectivityArray);

    for(int c = 0, q = 0; c < calibrations.nCameras; c++) {
      const int offset = c * nPointsPerFrustum;

      // lines to near plane corners
      for(int i = 1; i < 5; i++) {
        data[q++] = offset + 0;
        data[q++] = offset + i;
      }

      // lines between near and far plane corners
      for(int i = 1; i < 5; i++) {
        data[q++] = offset + i;
        data[q++] = offset + i + 4;
      }

      // near plane
      for(int i = 1; i < 5; i++) {
        data[q++] = offset + i;
        data[q++] = offset + (i > 3 ? 1 : i + 1);
      }

      // far plane
      for(int i = 5; i < 9; i++) {
        data[q++] = offset + i;
        data[q++] = offset + (i > 7 ? 5 : i + 1);
      }
    }
  }

  cells->SetData(offsetArray, connectivityArray);
  frustums->SetLines(cells);

  this->printMsg(msg, 1, timer.getElapsedTime(), 1);

  return 1;
}

int ttkCinemaDarkroomCamera::RequestData(vtkInformation *request,
                                         vtkInformationVector **inputVector,
                                         vtkInformationVector *outputVector) {

  auto cameras = vtkPolyData::GetData(outputVector, 0);

  int status = 0;

  status
    = this->ComputeCalibrations(cameras, vtkPointSet::GetData(inputVector[0]),
                                vtkDataObject::GetData(inputVector[1]),
                                vtkDataObject::GetData(inputVector[2]));
  if(!status)
    return 0;

  status
    = this->ComputeFrustums(cameras, vtkPolyData::GetData(outputVector, 1));
  if(!status)
    return 0;

  return 1;
}
