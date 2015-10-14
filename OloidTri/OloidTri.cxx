//
// Oloid - a ruled surface (a developable) based on two unit circles.
//   One circle is in the xy plane and centred at (0,0,0) the other
//   is centred at (1,0,0) in the xz plane. The ruled surface is
//   developed systematically in 4 phases, starting at the x axis
//   of the xy circle with straight lines traced to the other circle.
//   100 STEPS are taken over 2/3rds of a half circle in each phase.
//
//                                                 -=:LogicMonkey:=-

#include <vtkVersion.h>
#include <vtkSmartPointer.h>

#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkCellArray.h>
#include <vtkTriangleStrip.h>
#include <vtkTransform.h>
#include <vtkTransformPolyDataFilter.h>
#include <vtkAppendPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkProperty.h>
#include <vtkCamera.h>
#include <vtkTriangleFilter.h>
#include <vtkSTLWriter.h>

#define PI         3.141592653
#define VERTICES 20000
#define RES       25

int main(int, char *[]) {

  double hx, hy, vx, vz;

  // trace 4 arcs each 2/3 of a half circle (pi radians)
  // sub arc length 2/3 * Pi / number of steps
  double sub_arc_length = 2.0/3.0 * PI / (VERTICES/2-2);

  // first quarter

  double ht = PI;              // Pi decreasing to Pi/3
  double vt = 2.0f/3.0f * PI;  // 2/3 Pi decreasing to 0

  hx = cos(ht);
  hy = sin(ht);

  vx = cos(vt);
  vz = sin(vt);

  vtkSmartPointer<vtkPoints>
    points = vtkSmartPointer<vtkPoints>::New();
  points->InsertPoint(0, hx, hy, 0);
  points->InsertPoint(1, 1.0f + vx, 0, vz);

  vtkSmartPointer<vtkTriangleStrip>
    strip1= vtkSmartPointer<vtkTriangleStrip>::New();
  strip1->GetPointIds()->SetNumberOfIds(VERTICES);
  strip1->GetPointIds()->SetId(0,0);
  strip1->GetPointIds()->SetId(1,1);

  for( int i=1; i<VERTICES/2; i++ ) {

    ht -= sub_arc_length;
    vt -= sub_arc_length;

    hx = cos(ht);
    hy = sin(ht);

    vx = cos(vt);
    vz = sin(vt);

    points->InsertPoint(2*i, hx, hy, 0);
    points->InsertPoint(2*i+1, 1.0f + vx, 0, vz);

    strip1->GetPointIds()->SetId(2*i,2*i);
    strip1->GetPointIds()->SetId(2*i+1,2*i+1);
  }

  // second quarter

  ht = -PI;             // -Pi increasing to -Pi/3
  vt = 2.0f/3.0f * PI;  // 2/3 Pi decreasing to 0

  hx = cos(ht);
  hy = sin(ht);

  vx = cos(vt);
  vz = sin(vt);

  points->InsertPoint(VERTICES+0, hx, hy, 0);
  points->InsertPoint(VERTICES+1, 1.0f + vx, 0, vz);

  vtkSmartPointer<vtkTriangleStrip>
    strip2= vtkSmartPointer<vtkTriangleStrip>::New();
  strip2->GetPointIds()->SetNumberOfIds(VERTICES);
  strip2->GetPointIds()->SetId(0,VERTICES+0);
  strip2->GetPointIds()->SetId(1,VERTICES+1);

  for( int i=1; i<VERTICES/2; i++ ) {

    ht += sub_arc_length;
    vt -= sub_arc_length;

    hx = cos(ht);
    hy = sin(ht);

    vx = cos(vt);
    vz = sin(vt);

    points->InsertPoint(VERTICES+2*i, hx, hy, 0);
    points->InsertPoint(VERTICES+2*i+1, 1.0f + vx, 0, vz);

    strip2->GetPointIds()->SetId(2*i,VERTICES+2*i);
    strip2->GetPointIds()->SetId(2*i+1,VERTICES+2*i+1);
  }

  // third quarter

  ht = PI;              // Pi decreasing to Pi/3
  vt = -2.0f/3.0f * PI; // -2/3 Pi increasing to 0

  hx = cos(ht);
  hy = sin(ht);

  vx = cos(vt);
  vz = sin(vt);

  points->InsertPoint(2*VERTICES+0, hx, hy, 0);
  points->InsertPoint(2*VERTICES+1, 1.0f + vx, 0, vz);

  vtkSmartPointer<vtkTriangleStrip>
    strip3= vtkSmartPointer<vtkTriangleStrip>::New();
  strip3->GetPointIds()->SetNumberOfIds(VERTICES);
  strip3->GetPointIds()->SetId(0,2*VERTICES+0);
  strip3->GetPointIds()->SetId(1,2*VERTICES+1);

  for( int i=1; i<VERTICES/2; i++ ) {

    ht -= sub_arc_length;
    vt += sub_arc_length;

    hx = cos(ht);
    hy = sin(ht);

    vx = cos(vt);
    vz = sin(vt);

    points->InsertPoint(2*VERTICES+2*i, hx, hy, 0);
    points->InsertPoint(2*VERTICES+2*i+1, 1.0f + vx, 0, vz);

    strip3->GetPointIds()->SetId(2*i,2*VERTICES+2*i);
    strip3->GetPointIds()->SetId(2*i+1,2*VERTICES+2*i+1);
  }

  // fourth quarter

  ht = -PI;             // -Pi increasing to -Pi/3
  vt = -2.0f/3.0f * PI; // -2/3 Pi increasing to 0

  hx = cos(ht);
  hy = sin(ht);

  vx = cos(vt);
  vz = sin(vt);

  points->InsertPoint(3*VERTICES+0, hx, hy, 0);
  points->InsertPoint(3*VERTICES+1, 1.0f + vx, 0, vz);

  vtkSmartPointer<vtkTriangleStrip>
    strip4= vtkSmartPointer<vtkTriangleStrip>::New();
  strip4->GetPointIds()->SetNumberOfIds(VERTICES);
  strip4->GetPointIds()->SetId(0,3*VERTICES+0);
  strip4->GetPointIds()->SetId(1,3*VERTICES+1);

  for( int i=1; i<VERTICES/2; i++ ) {

    ht += sub_arc_length;
    vt += sub_arc_length;

    hx = cos(ht);
    hy = sin(ht);

    vx = cos(vt);
    vz = sin(vt);

    points->InsertPoint(3*VERTICES+2*i, hx, hy, 0);
    points->InsertPoint(3*VERTICES+2*i+1, 1.0f + vx, 0, vz);

    strip4->GetPointIds()->SetId(2*i,3*VERTICES+2*i);
    strip4->GetPointIds()->SetId(2*i+1,3*VERTICES+2*i+1);
  }

  vtkSmartPointer<vtkPolyData>
    polydata = vtkSmartPointer<vtkPolyData>::New();
  polydata->SetPoints(points);

  vtkSmartPointer<vtkCellArray>
    cells = vtkSmartPointer<vtkCellArray>::New();
  cells->InsertNextCell(strip1);
  cells->InsertNextCell(strip2);
  cells->InsertNextCell(strip3);
  cells->InsertNextCell(strip4);

  polydata->SetStrips(cells);

  // write out a STL file
  //   run a triangle filter on the mesh to create polys
  vtkSmartPointer<vtkTriangleFilter>
    tris = vtkTriangleFilter::New();
  tris->SetInputData(polydata);

  vtkSmartPointer<vtkSTLWriter>
    stlWriter = vtkSmartPointer<vtkSTLWriter>::New();
  stlWriter->SetFileName("OloidTri.stl");
  stlWriter->SetInputConnection(tris->GetOutputPort());
  stlWriter->Write();

  // now the usual VTK render and displey pipeline

  vtkSmartPointer<vtkPolyDataMapper>
    mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
  mapper->SetInputData(polydata);                      // data to mapper

  vtkSmartPointer<vtkActor>
    actor = vtkSmartPointer<vtkActor>::New();
  actor->SetMapper(mapper);                            // mapper to actor

  actor->GetProperty()->SetColor(0.89, 0.81, 0.34);    // goldish
//actor->GetProperty()->SetColor(0.75, 0.75, 0.75);    // silver, really?
  actor->GetProperty()->SetSpecular(1.0);

  // Add the actors to the renderer, set the background and size
  vtkSmartPointer<vtkRenderer>
    renderer = vtkSmartPointer<vtkRenderer>::New();

  vtkSmartPointer<vtkRenderWindow>
    renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
  renderWindow->AddRenderer(renderer);

  vtkSmartPointer<vtkRenderWindowInteractor>
    interactor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
  interactor->SetRenderWindow(renderWindow);

  renderer->AddActor(actor);
  renderer->SetBackground(.3, .4, .5);

  renderer->ResetCamera();
  renderer->GetActiveCamera()->Azimuth(60);
  renderer->GetActiveCamera()->Elevation(60);
  renderer->GetActiveCamera()->Dolly(1.2);
  renderer->ResetCameraClippingRange();

  renderWindow->Render();
  interactor->Start();

  return EXIT_SUCCESS;
}
