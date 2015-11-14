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
// rOloid - an Oloid using rational parametrization.
//
// I consider the definitive text on the Oloid to be Dirnbock and Stachel's
// paper of 1997:
//
//     http://www.heldermann-verlag.de/jgg/jgg01_05/jgg0113.pdf
//
// In this software the following original alternative paramatrizaton
// is presented:
//
//              ,                      .
//             |   2t    -1   1-t^2     |
//         U = | ----- ,  - - ----- , 0 |
//             | 1+t^2    2   1+t^2     |
//              '                      '
//
//              ,                                  .
//             |         sqrt((3-t^2)(1+t^2))   t^2 |
//         V = | 0 , +/- -------------------- , --- |
//             |                   2             2  |
//              '                                  '
//
//         t:(-sqrt(3)..sqrt(3))
//
// This predominantly rational parametrization avoids the transcendental
// functions and singularities exhibited by Dirnbock and Stachel's solution.
//
// Derivation of this implementation and that of D&S is here:
//
// link to my PDF
//
//                                                         -=:LogicMonkey:=-

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

#define UMAX 28

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

  double dt = 2.0 * sqrt( 3.0 ) / UMAX;
  double t  = -sqrt( 3.0 );

  std::cout << "GEOMETRY" << std::endl;

  for( vtkIdType i=UMIN; i<=UMAX; i++ ) {

    p[0] = RADIUS * 2.0 * t/(1 + t*t);
    p[1] = RADIUS * ( -0.5 - (1-t*t)/(1 + t*t));

    // reflect in x
    points->InsertNextPoint( p[0],  p[1], 0 );
    std::cout << i << " U " << p[0] << " : " << p[1] << std::endl;

    t += dt;
  }

  t = sqrt( 3.0 ) - UMID*dt;

  for( vtkIdType i=VMIN; i<=VMID-1; i++ ) {

    p[0] = RADIUS * sqrt((3.0 - t*t)*(1.0 + t*t)) / 2.0;
    p[1] = RADIUS * t*t / 2.0;

    points->InsertNextPoint( 0, p[1], p[0] );
    std::cout << i << " V " << p[1] << " : " << p[0] << std::endl;

    t += dt;
  }

  // midpoint of V
  points->InsertNextPoint( 0, RADIUS * 1.5, 0 );
  std::cout << VMID << " V " << RADIUS * 1.5 << " : " << 0 << std::endl;

  t = sqrt( 3.0 ) - dt;

  for( vtkIdType i=VMID+1; i<=VMAX; i++ ) {

    p[0] = RADIUS * sqrt((3.0 - t*t)*(1.0 + t*t)) / 2.0;
    p[1] = RADIUS * t*t / 2.0;

    points->InsertNextPoint( 0, p[1], -p[0] );
    std::cout << i << " V " << p[1] << " : " << -p[0] << std::endl;

    t -= dt;
  }

  // TOPOLOGY PHASE (create point interconnect: tri strip)

  // Run the triangulation process over 4 sub-strips. This is the simplest
  // way to do it without incurring artefacts at the turning points as
  // a result of vertex sharing (local vertex fans)

  vtkIdType strip[4*(UMAX+2)];
  int strip_idx = 0;

  vtkIdType v = VMID; // VMID to VMAX

  for( vtkIdType u=UMIN; u<=UMID; u++ ) {
    strip[strip_idx++] = v++;
    strip[strip_idx++] = u;
  }

  v = VMIN;           // VMIN to VMID

  for( vtkIdType u=UMID; u>=UMIN; u-- ) {
    strip[strip_idx++] = u;
    strip[strip_idx++] = v++;
  }

  v = VMID;           // VMID to VMAX

  for( vtkIdType u=UMAX; u>=UMID; u-- ) {
    strip[strip_idx++] = v++;
    strip[strip_idx++] = u;
  }

  v = VMIN;           // VMIN to VMID

  for( vtkIdType u=UMID; u<=UMAX; u++ ) {
    strip[strip_idx++] = u;
    strip[strip_idx++] = v++;
  }

  vtkSmartPointer<vtkCellArray>
    cells = vtkSmartPointer<vtkCellArray>::New();
  cells->InsertNextCell( 4*(UMAX+2), strip );

  std::cout << "TOPOLOGY" << std::endl;
  for( vtkIdType s=0; s<4*(UMAX+2); s+=2 ) {
    std::cout << s << " : " << strip[s] << " : " << strip[s+1] << std::endl;
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
