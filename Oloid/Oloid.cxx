//
//   Copyright (C) 2015 by Piers Barber
//
//   This program is free software; you can redistribute it and/or modify
//   it under the terms of the GNU General Public License as published by
//   the Free Software Foundation; either version 2 of the License, or
//   (at your option) any later version.
//
//   This program is distributed in the hope that it will be useful,
//   but WITHOUT ANY WARRANTY; without even the implied warranty of
//   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//   GNU General Public License for more details.
//
//   You should have received a copy of the GNU General Public License
//   along with this program; if not, write to the
//   Free Software Foundation, Inc.,
//   51 Franklin Street, Fifth Floor,Boston, MA 02110-1301 USA
//
// -----------------------------------------------------------------------
//
// Oloid - a ruled surface based on two unit circles.
//   One circle is in the xy plane and centred at (0,0,0) the other
//   is centred at (1,0,0) in the xz plane. The ruled surface is swept
//   out by a straight line from one circle to the other, over the arc
//   with range -2*pi/3 to 2*pi/3. The mid point of each arc connects
//   to the end points of the arc of the corresponding circle.
//
//   This implementation calculates all the required triangle mesh vertices
//   on each arc up front (geometry stage) and then stitches them together
//   into a single strip that closes on itself (toplology stage).
//   The triangle strip is generated in four phases as contiguous
//   sub-strips each with a turning point at the end of one arc. Vertices
//   are re-used between phases in order to ensure there are no holes.
//   The triangle mesh is created 2 triangles at at time by adding 2
//   vertices to an initial pair. A slightly different mesh structure
//   comes from using an odd or even vertex count. The odd or even
//   approaches are separated for readability.
//   The resulting structure is well behaved and suitable for 3D printing.
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

#define PI 4.0*atan(1)                    // pi for this computer

#define UMAX 128

#define RADIUS 25.4                       // imperious :)

int main(int, char *[]) {

  const int UMIN = 0;
  const int UMID = UMAX/2;
  const int VMIN = UMAX+1;
  const int VMID = (3*UMAX+2)/2;
  const int VMAX = 2*UMAX+1;

  double p[2];

  vtkSmartPointer<vtkPoints>
    points = vtkSmartPointer<vtkPoints>::New();

  // GEOMETRY PHASE (create point set)

  double dt = 4.0/3.0 * PI / UMAX;
  double t = -2.0/3.0 * PI;

  std::cout << "GEOMETRY" << std::endl;

  for( vtkIdType i=UMIN; i<=UMAX; i++ ) {

    p[0] = RADIUS * cos(t);
    p[1] = RADIUS * sin(t);

    // reflect in x
    points->InsertNextPoint( -p[0],  p[1], 0 );
    std::cout << i << " U " << p[0] << " : " << p[1] << std::endl;

    t += dt;
  }

  t = 2.0/3.0 * PI - UMID*dt;

  for( vtkIdType i=VMIN; i<=VMID-1; i++ ) {

    p[0] = RADIUS * (1.0 - cos(t)/(1.0+cos(t)));
    p[1] = RADIUS * (sqrt(1.0+2.0*cos(t))/(1.0+cos(t)));

    points->InsertNextPoint( p[0], 0,  p[1] );
    std::cout << i << " V " << p[0] << " : " << p[1] << std::endl;

    t += dt;
  }

  // avoid divide by zero at the midpoint of V
  points->InsertNextPoint( 2*RADIUS, 0, 0 );
  std::cout << VMID << " V " << p[0] << " : " << p[1] << std::endl;

  t = 2.0/3.0 * PI - dt;

  for( vtkIdType i=VMID+1; i<=VMAX; i++ ) {

    p[0] = RADIUS * (1.0 - cos(t)/(1.0+cos(t)));
    p[1] = RADIUS * (sqrt(1.0+2.0*cos(t))/(1.0+cos(t)));

    points->InsertNextPoint( p[0], 0, -p[1] );
    std::cout << i << " V " << p[0] << " : " << p[1] << std::endl;

    t -= dt;
  }

  // TOPOLOGY PHASE (create point interconnect: tri strip)

  // separate the two stitching mechanisms for readability
  // the odd vertices case needs no special handling at the end of
  // each quarter. the even case has an extra triangle to be added
  // in order to fill a hole when adding vertices to the strip in pairs

  vtkIdType strip1[UMAX+2];
  vtkIdType strip2[UMAX+2];
  vtkIdType strip3[UMAX+2];
  vtkIdType strip4[UMAX+2];

  int strip_idx = 0;
  vtkIdType v = VMID; // VMID to VMAX

  for( vtkIdType u=UMIN; u<=UMID; u++ ) {
    strip1[strip_idx++] = v++;
    strip1[strip_idx++] = u;
  }

  strip_idx = 0;
  v = VMIN;           // VMIN to VMID

  for( vtkIdType u=UMID; u>=UMIN; u-- ) {
    strip2[strip_idx++] = u;
    strip2[strip_idx++] = v++;
  }

  strip_idx = 0;
  v = VMID;           // VMID to VMAX

  for( vtkIdType u=UMAX; u>=UMID; u-- ) {
    strip3[strip_idx++] = v++;
    strip3[strip_idx++] = u;
  }

  strip_idx = 0;
  v = VMIN;           // VMIN to VMID

  for( vtkIdType u=UMID; u<=UMAX; u++ ) {
    strip4[strip_idx++] = u;
    strip4[strip_idx++] = v++;
  }

  vtkSmartPointer<vtkCellArray>
    cells = vtkSmartPointer<vtkCellArray>::New();
  cells->InsertNextCell( UMAX+2, strip1 );
  cells->InsertNextCell( UMAX+2, strip2 );
  cells->InsertNextCell( UMAX+2, strip3 );
  cells->InsertNextCell( UMAX+2, strip4 );

  std::cout << "TOPOLOGY" << std::endl;
  for( vtkIdType s=0; s<UMAX+2; s+=2 ) {
    std::cout << s << " : " << strip1[s] << " : " << strip1[s+1] << std::endl;
  }
  for( vtkIdType s=0; s<UMAX+2; s+=2 ) {
    std::cout << s << " : " << strip2[s] << " : " << strip2[s+1] << std::endl;
  }
  for( vtkIdType s=0; s<UMAX+2; s+=2 ) {
    std::cout << s << " : " << strip3[s] << " : " << strip3[s+1] << std::endl;
  }
  for( vtkIdType s=0; s<UMAX+2; s+=2 ) {
    std::cout << s << " : " << strip4[s] << " : " << strip4[s+1] << std::endl;
  }

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
  stlWriter->SetFileName("Oloid.stl");
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

  renderWindow->Render();
  interactor->Start();

  return EXIT_SUCCESS;
}
