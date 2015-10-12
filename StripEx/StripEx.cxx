#include <vtkVersion.h>
#include <vtkSmartPointer.h>
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
#include <vtkTriangleStrip.h>
#include <vtkSTLWriter.h>

#define VERTICES 100

int main(int, char *[]) {

  vtkSmartPointer<vtkPoints>
    points = vtkSmartPointer<vtkPoints>::New();
  points->InsertPoint(0, 0, 0, 0);
  points->InsertPoint(1, 0, 1, 0);

  vtkSmartPointer<vtkTriangleStrip>
    strip = vtkSmartPointer<vtkTriangleStrip>::New();
  strip->GetPointIds()->SetNumberOfIds(VERTICES);
  strip->GetPointIds()->SetId(0,0);
  strip->GetPointIds()->SetId(1,1);

  for( int i=1; i<VERTICES/2; i++ ) {

    points->InsertPoint(2*i,   (double) i, 0, 0);
    points->InsertPoint(2*i+1, (double) i, 1, 0);

    strip->GetPointIds()->SetId(2*i,2*i);
    strip->GetPointIds()->SetId(2*i+1,2*i+1);
  }

  vtkSmartPointer<vtkPolyData>
    polydata = vtkSmartPointer<vtkPolyData>::New();

  polydata->SetPoints(points);

  vtkSmartPointer<vtkCellArray>
    cells = vtkSmartPointer<vtkCellArray>::New();
  cells->InsertNextCell(strip);

  polydata->SetStrips(cells);

  // STL file output

  vtkSmartPointer<vtkSTLWriter>
    stlWriter = vtkSmartPointer<vtkSTLWriter>::New();
  stlWriter->SetFileName("StripEx.stl");
  stlWriter->SetInputData(polydata);
  stlWriter->Write();

  // now the usual VTK render and displey pipeline

  vtkSmartPointer<vtkPolyDataMapper>
    mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
  mapper->SetInputData(polydata);                    // data to mapper

  vtkSmartPointer<vtkActor>
    actor = vtkSmartPointer<vtkActor>::New();
  actor->SetMapper(mapper);                            // mapper to actor

  actor->GetProperty()->SetColor(0.89, 0.81, 0.34);
  actor->GetProperty()->SetRepresentationToWireframe();

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
