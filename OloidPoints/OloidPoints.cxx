#include <vtkVersion.h>
#include <vtkSmartPointer.h>
#include <vtkPoints.h>
#include <vtkCellArray.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkProperty.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkRenderWindowInteractor.h>

#include <vtkLabeledDataMapper.h>
#include <vtkActor2D.h>

#define RADIUS 1.0
#define VMAX 8
#define PI atan(1)*4.0

int main(int, char *[])
{
  // Create the geometry of a point (the coordinate)
  vtkSmartPointer<vtkPoints>
    points = vtkSmartPointer<vtkPoints>::New();

  vtkIdType pid[2*(VMAX+1)];

  double dt = 4.0/3.0 * PI / VMAX;

  double t = -2.0/3.0 * PI;

  double p[2];

  for( vtkIdType i=0; i<=VMAX; i++ ) {

    p[0] = RADIUS * cos(t);
    p[1] = RADIUS * sin(t);

    // reflect in x
    pid[i] = points->InsertNextPoint( -p[0],  p[1], 0 );

    t += dt;
  }

  t = 2.0/3.0 * PI - VMAX/2*dt;

  for( vtkIdType i=VMAX+1; i<=(3*VMAX+2)/2; i++ ) {


    p[0] = RADIUS * (1.0 - cos(t)/(1.0+cos(t)));
    p[1] = RADIUS * (sqrt(1.0+2.0*cos(t))/(1.0+cos(t)));

    pid[i] = points->InsertNextPoint( p[0], 0,  p[1] );

    t += dt;
  }

  t = 2.0/3.0 * PI - dt;

  for( vtkIdType i=(3*VMAX+2)/2+1; i<=2*VMAX+1; i++ ) {

    p[0] = RADIUS * (1.0 - cos(t)/(1.0+cos(t)));
    p[1] = RADIUS * (sqrt(1.0+2.0*cos(t))/(1.0+cos(t)));

    pid[i] = points->InsertNextPoint( p[0], 0, -p[1] );

    t -= dt;
  }

  vtkSmartPointer<vtkCellArray>
    vertices = vtkSmartPointer<vtkCellArray>::New();
  vertices->InsertNextCell(2*(VMAX+1),pid);

  vtkSmartPointer<vtkPolyData>
    polydata = vtkSmartPointer<vtkPolyData>::New();
  polydata->SetPoints(points);
  polydata->SetVerts(vertices);

  vtkSmartPointer<vtkLabeledDataMapper>
    labelMapper = vtkSmartPointer<vtkLabeledDataMapper>::New();
  labelMapper->SetInputData(polydata);

  vtkSmartPointer<vtkActor2D>
    labelActor = vtkSmartPointer<vtkActor2D>::New();
  labelActor->SetMapper(labelMapper);

  // now the usual VTK render and display pipeline

  vtkSmartPointer<vtkPolyDataMapper>
    mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
  mapper->SetInputData(polydata);

  vtkSmartPointer<vtkActor>
    actor = vtkSmartPointer<vtkActor>::New();
  actor->SetMapper(mapper);                            // mapper to actor

  actor->GetProperty()->SetColor(0.89, 0.81, 0.34);    // goldish
  //actor->GetProperty()->SetSpecular(1.0);
  actor->GetProperty()->SetPointSize(5);

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
  renderer->AddActor(labelActor);
  renderer->SetBackground(.3, .4, .5);

  renderWindow->Render();
  interactor->Start();
  return EXIT_SUCCESS;
}
