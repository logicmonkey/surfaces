//
// Oloid - a ruled surface (a developable) based on two unit circles.
//   One circle is in the xy plane and centred at (0,0,0) the other
//   is centred at (1,0,0) in the xz plane. At least, that would be
//   the case, but I shift everything left by 0.5 to position the
//   centre of mass at O.
//
//   This implementation generates all the required vertices on each
//   circle up front (geometry phase) and then stitches them together
//   into a single strip that closes on itself (toplology phase).
//
//                                                 -=:LogicMonkey:=-

#include <vtkVersion.h>
#include <vtkTransform.h>
#include <vtkProperty.h>

#include <vtkSmartPointer.h>              // for points
#include <vtkCellArray.h>                 // tri strip in cells
#include <vtkTransformPolyDataFilter.h>   // for polys
#include <vtkAppendPolyData.h>
#include <vtkPolyDataMapper.h>

#include <vtkCamera.h>                    // the visualisation
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>

#include <vtkTriangleFilter.h>            // needed for STL output
#include <vtkSTLWriter.h>                 // physical realisation

#define PI   4*atan(1)
#define VMAX 401                          // maximum vertex index

#define RADIUS 25.4                       // imperious :)

int main(int, char *[]) {

  double p[2];

  double sub_arc_length = 4.0/3.0 * PI / VMAX;

  vtkSmartPointer<vtkPoints>
    points = vtkSmartPointer<vtkPoints>::New();

  // GEOMETRY PHASE (create point set)

  // trace out all the vertices on the horizontal circle
  // just one loop of those nasty transcendental calcs

  double theta = -2.0/3.0 * PI;

  for( int i=0; i<=VMAX; i++ ) {

    p[0] = RADIUS * (cos( theta ) + 0.5);
    p[1] = RADIUS * sin( theta );

    points->InsertNextPoint( p[0], p[1], 0 );

    theta += sub_arc_length;
  }

  // copy each point on the horizontal circle to a point on
  // the vertical one, reflecting in x at the same time

  for( int i=0; i<=VMAX; i++ ) {
    points->GetPoint( i, p );
    points->InsertNextPoint( -p[0], 0, p[1] );
  }

  // TOPOLOGY PHASE (create point interconnect: tri strip)

  // separate the two stitching mechanisms for readability
  // the odd vertices case needs no special handling at the end of
  // each quarter. the even case has an extra triangle to be added
  // in order to fill a hole when adding vertices to the strip in pairs

  vtkIdType strip[4*VMAX+2];
  int j = 0;

  if( VMAX%2 == 0 ) { // odd vertices (VMAX+1), even steps

    vtkIdType t = VMAX+1 + VMAX/2;

    // s increasing, t increasing
    for( vtkIdType s=0; s<=VMAX/2; s++ ) {
      strip[j++] = s;
      strip[j++] = t++;
    }
    t--;  // hit N, post inc to N+1 so rewind to N again
    // s increasing, t decreasing
    for( vtkIdType s=VMAX/2+1; s<=VMAX; s++ ) {
      strip[j++] = s;
      strip[j++] = --t;
    }
    // s decreasing, t decreasing
    for( vtkIdType s=VMAX-1; s>=VMAX/2; s-- ) {
      strip[j++] = s;
      strip[j++] = --t;
    }
    // s decreasing, t increasing
    for( vtkIdType s=VMAX/2-1; s>=0; s-- ) {
      strip[j++] = s;
      strip[j++] = ++t;
    }

  } else { // even vertices, odd steps

    vtkIdType t = VMAX+1 + (VMAX+1)/2;

    // s increasing, t increasing
    for( vtkIdType s=0; s<=(VMAX+1)/2-1; s++ ) {
      strip[j++] = s;
      strip[j++] = t++;
    }
    strip[j++] = (VMAX+1)/2;
    t--;  // hit N, post inc to N+1 so rewind to N again
    // s increasing, t decreasing
    for( vtkIdType s=(VMAX+1)/2+1; s<=VMAX; s++ ) {
      strip[j++] = s;
      strip[j++] = --t;
    }
    strip[j++] = --t;
    // s decreasing, t decreasing
    for( vtkIdType s=VMAX-1; s>=(VMAX+1)/2; s-- ) {
      strip[j++] = s;
      strip[j++] = --t;
    }
    strip[j++] = (VMAX+1)/2-1;
    // s decreasing, t increasing
    for( vtkIdType s=(VMAX+1)/2-2; s>=0; s-- ) {
      strip[j++] = s;
      strip[j++] = ++t;
    }
    strip[j++] = ++t;

  }

  vtkSmartPointer<vtkCellArray>
    cells = vtkSmartPointer<vtkCellArray>::New();
  cells->InsertNextCell(j, strip);

  //for( vtkIdType s=0; s<j; s++ ) {
  //  std::cout << s << " : " << strip[s] << std::endl;
  //}

  vtkSmartPointer<vtkPolyData>
    polydata = vtkSmartPointer<vtkPolyData>::New();
  polydata->SetPoints(points);
  polydata->SetStrips(cells);

  // write out STL

  // need for this filter has been designed out in leading edge VTK
  // for vtk < 6.3.0 at least, you need it for non-empty STL files
  vtkSmartPointer<vtkTriangleFilter>
    tris = vtkTriangleFilter::New();
  tris->SetInputData(polydata);

  vtkSmartPointer<vtkSTLWriter>
    stlWriter = vtkSmartPointer<vtkSTLWriter>::New();
  stlWriter->SetFileName("OloidMesh.stl");
  stlWriter->SetInputConnection(tris->GetOutputPort());
  stlWriter->Write();

  // now the usual VTK render and display pipeline

  vtkSmartPointer<vtkPolyDataMapper>
    mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
  mapper->SetInputData(polydata);

  vtkSmartPointer<vtkActor>
    actor = vtkSmartPointer<vtkActor>::New();
  actor->SetMapper(mapper);                            // mapper to actor

  actor->GetProperty()->SetColor(0.89, 0.81, 0.34);    // goldish
  //actor->GetProperty()->SetSpecular(1.0);

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
