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

#include <vtkRuledSurfaceFilter.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkCellArray.h>
#include <vtkTransform.h>
#include <vtkTransformPolyDataFilter.h>
#include <vtkAppendPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkProperty.h>
#include <vtkCamera.h>
#include <vtkLine.h>
#include <vtkTriangleFilter.h>
#include <vtkSTLWriter.h>

#define PI    3.141592653
#define STEPS 10
#define RES   25

#define RADIUS 25.0

int main(int, char *[]) {

  double hx, hy, vx, vz;

  // trace 4 arcs each 2/3 of a half circle (pi radians)
  // sub arc length 2/3 * Pi / number of steps
  double sub_arc_length = 2.0/3.0 * PI / STEPS;

  vtkSmartPointer<vtkPoints>
    points = vtkSmartPointer<vtkPoints>::New();
  vtkSmartPointer<vtkLine>
    hsecant = vtkSmartPointer<vtkLine>::New();
  vtkSmartPointer<vtkLine>
    vsecant = vtkSmartPointer<vtkLine>::New();
  vtkSmartPointer<vtkCellArray>
    lines = vtkSmartPointer<vtkCellArray>::New();
  vtkSmartPointer<vtkPolyData>
    polydata = vtkSmartPointer<vtkPolyData>::New();
  vtkSmartPointer<vtkRuledSurfaceFilter>
    ruled = vtkSmartPointer<vtkRuledSurfaceFilter>::New();

  // first quarter

  double ht = PI;              // Pi decreasing to Pi/3
  double vt = 2.0f/3.0f * PI;  // 2/3 Pi decreasing to 0

  hx = RADIUS * (cos(ht) - 0.5);
  hy = RADIUS * sin(ht);

  vx = RADIUS * (cos(vt) + 0.5);
  vz = RADIUS * sin(vt);

  int n = 0;

  points->InsertPoint(n++, hx, hy, 0);
  points->InsertPoint(n++, vx, 0, vz);

  for( int i=0; i<STEPS; i++ ) {

    ht -= sub_arc_length;
    vt -= sub_arc_length;

    hx = RADIUS * (cos(ht) - 0.5);
    hy = RADIUS * sin(ht);

    vx = RADIUS * (cos(vt) + 0.5);
    vz = RADIUS * sin(vt);

    points->InsertPoint(n, hx, hy, 0);
    points->InsertPoint(n+1, vx, 0, vz);

    hsecant->GetPointIds()->SetId(0,n-2);
    hsecant->GetPointIds()->SetId(1,n);

    vsecant->GetPointIds()->SetId(0,n-1);
    vsecant->GetPointIds()->SetId(1,n+1);

    n += 2;

    lines->InsertNextCell(hsecant);
    lines->InsertNextCell(vsecant);
  }

  // second quarter

  ht = -PI;             // -Pi increasing to -Pi/3
  vt = 2.0f/3.0f * PI;  // 2/3 Pi decreasing to 0

  hx = RADIUS * (cos(ht) - 0.5);
  hy = RADIUS * sin(ht);

  vx = RADIUS * (cos(vt) + 0.5);
  vz = RADIUS * sin(vt);

  points->InsertPoint(n++, hx, hy, 0);
  points->InsertPoint(n++, vx, 0, vz);

  for( int i=0; i<STEPS; i++ ) {

    ht += sub_arc_length;
    vt -= sub_arc_length;

    hx = RADIUS * (cos(ht) - 0.5);
    hy = RADIUS * sin(ht);

    vx = RADIUS * (cos(vt) + 0.5);
    vz = RADIUS * sin(vt);

    points->InsertPoint(n, hx, hy, 0);
    points->InsertPoint(n+1, vx, 0, vz);

    hsecant->GetPointIds()->SetId(0,n-2);
    hsecant->GetPointIds()->SetId(1,n);

    vsecant->GetPointIds()->SetId(0,n-1);
    vsecant->GetPointIds()->SetId(1,n+1);

    n += 2;

    lines->InsertNextCell(hsecant);
    lines->InsertNextCell(vsecant);
  }

  // third quarter

  ht = PI;              // Pi decreasing to Pi/3
  vt = -2.0f/3.0f * PI; // -2/3 Pi increasing to 0

  hx = RADIUS * (cos(ht) - 0.5);
  hy = RADIUS * sin(ht);

  vx = RADIUS * (cos(vt) + 0.5);
  vz = RADIUS * sin(vt);

  points->InsertPoint(n++, hx, hy, 0);
  points->InsertPoint(n++, vx, 0, vz);

  for( int i=0; i<STEPS; i++ ) {

    ht -= sub_arc_length;
    vt += sub_arc_length;

    hx = RADIUS * (cos(ht) - 0.5);
    hy = RADIUS * sin(ht);

    vx = RADIUS * (cos(vt) + 0.5);
    vz = RADIUS * sin(vt);

    points->InsertPoint(n, hx, hy, 0);
    points->InsertPoint(n+1, vx, 0, vz);

    hsecant->GetPointIds()->SetId(0,n-2);
    hsecant->GetPointIds()->SetId(1,n);

    vsecant->GetPointIds()->SetId(0,n-1);
    vsecant->GetPointIds()->SetId(1,n+1);

    n += 2;

    lines->InsertNextCell(hsecant);
    lines->InsertNextCell(vsecant);
  }

  // fourth quarter

  ht = -PI;             // -Pi increasing to -Pi/3
  vt = -2.0f/3.0f * PI; // -2/3 Pi increasing to 0

  hx = RADIUS * (cos(ht) - 0.5);
  hy = RADIUS * sin(ht);

  vx = RADIUS * (cos(vt) + 0.5);
  vz = RADIUS * sin(vt);

  points->InsertPoint(n++, hx, hy, 0);
  points->InsertPoint(n++, vx, 0, vz);

  for( int i=0; i<STEPS; i++ ) {

    ht += sub_arc_length;
    vt += sub_arc_length;

    hx = RADIUS * (cos(ht) - 0.5);
    hy = RADIUS * sin(ht);

    vx = RADIUS * (cos(vt) + 0.5);
    vz = RADIUS * sin(vt);

    points->InsertPoint(n, hx, hy, 0);
    points->InsertPoint(n+1, vx, 0, vz);

    hsecant->GetPointIds()->SetId(0,n-2);
    hsecant->GetPointIds()->SetId(1,n);

    vsecant->GetPointIds()->SetId(0,n-1);
    vsecant->GetPointIds()->SetId(1,n+1);

    n += 2;

    lines->InsertNextCell(hsecant);
    lines->InsertNextCell(vsecant);

  }

  polydata->SetPoints(points);
  polydata->SetLines(lines);

  ruled->SetInputData(polydata);   // pick up the poly data to run filter
  ruled->SetResolution(RES,RES);
  ruled->SetRuledModeToResample(); // apply VTK ruled surface wizardry

  // write out a STL file
  //   run a triangle filter on the ruled surface to create polys
  vtkSmartPointer<vtkTriangleFilter>
    tris = vtkTriangleFilter::New();
  tris->SetInputConnection(ruled->GetOutputPort());

  vtkSmartPointer<vtkSTLWriter>
    stlWriter = vtkSmartPointer<vtkSTLWriter>::New();
  stlWriter->SetFileName("OloidRuled.stl");
  stlWriter->SetInputConnection(tris->GetOutputPort());
  stlWriter->Write();

  // now the usual VTK render and displey pipeline

  vtkSmartPointer<vtkPolyDataMapper>
    mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
  mapper->SetInputConnection(ruled->GetOutputPort());  // data to mapper

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
