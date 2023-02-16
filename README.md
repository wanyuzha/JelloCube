# Jello Cube

This is the assignment#1 of CSCI520: Computer Animation and Simulation.

## Description
This is the simulation of a jello cube inside a big bounding box. Users could interact with the jello cube by selecting 
the mass point with mouse.
### Overview
<table>
    <tr>
        <td><img src="https://github.com/wanyuzha/JelloCube/blob/main/img/jello.gif" width="320" height="240" alt=""></td>
        <td><img src="https://github.com/wanyuzha/JelloCube/blob/main/img/wireframes.gif" width="320" height="240" alt=""></td>
    </tr>
    <tr>
        <td align="center">jello surface model</td>
        <td align="center">jello wireframe model</td>
    </tr>
    <tr>
        <td><img src="https://github.com/wanyuzha/JelloCube/blob/main/img/select.gif" width="320" height="240" alt=""></td>
        <td><img src="https://github.com/wanyuzha/JelloCube/blob/main/img/surface.gif" width="320" height="240" alt=""></td>
    </tr>
    <tr>
        <td align="center">jello wireframe model + selection</td>
        <td align="center">jello surface model + selection (not very easy to select because the point is hidden)</td>
    </tr>
</table>

### Features
- Applied Hook Law and Damp Law to the mass spring system, including structural, shear, and bend springs
- Collision Detection with bounding box and inclined plane
- Integration: elastic, damping, external forces, and collision bounce back
- Interactive selection: User can choose mass point to pull 

## Install
```shell
mkdir build
cd build
cmake ..
make 
```
Or simply use jello Application in the directory
```shell
./jello world/jello.w
```
## Documentation
This documentation is to clarify the basic steps of simulator
### Integrator
* For each mass point inside the mass system
  * Call computeSprings to get back all the structural,shear,bend springs
  * Call computeHookLaw, computeDampingLaw function to compute Force F
  * Call computeCollisionSprings to get back all the collision springs with bounding box
  * Call computeHookLaw, computeDampingLaw function to compute Force F
  * Call computeInclinedCollision to get back all the collision springs with inclined plane
  * Call computeHookLaw, computeDampingLaw function to compute Force F
  * Apply the force field with tri-linear interpolation
$$ F = F_{hook} + F_{damp} + F{Force field} $$

$$ F_{hook} = F_{normal spring hook} + F_{collision spring hook} + F{inclined plane collision spring hook} $$
$$ F_{damp} = F_{normal spring damp} + F_{collision spring damp} + F{inclined plane collision spring damp} $$
### OpenGL Selection
Reference: [Selection and Feedback](https://www.glprogramming.com/red/chapter13.html) (Chapter 13)
#### Selection by gluPickMatrix
Generally, selection is applied by SELECT_MODE;
- Enter Select Mode
- Load Matrix gluPickMatrix -> gluPerspective
- glLoadName() give each entity an ID
- Enter Render Mode, and trace back the hit data from given area
- Post Process
```c++
    glSelectBuffer (BUFSIZE, selectBuf);
    (void) glRenderMode (GL_SELECT);

    glInitNames();
    glPushName(0);

    glMatrixMode (GL_PROJECTION);
    glPushMatrix ();
    glLoadIdentity ();

    gluPickMatrix ((GLdouble) x, (GLdouble) (viewport[3] - y),
                   5.0, 5.0, viewport);
    
    gluPerspective(60.0f, aspectRatio, 0.01f, 1000.0f);
    
    glLoadName();
    glDrawSomething;

    glMatrixMode (GL_PROJECTION);
    glPopMatrix ();
    glFlush ();

    hits = glRenderMode (GL_RENDER);
    processHits (hits, selectBuf); // process hit data
```
The select area is restricted to 5 X 5, which is a pretty small area around cursor.

#### Move Object Towards Mouse Direction

To give users a feeling of pulling a point, we have to move the point along the direction 
of mouse movement. However, the mouse coordinate is in screen space, we have to transform it back to
world coordinate.

```c++
void glUntransform(double x, double y, double z, double *objX, double *objY, double *objZ)
{
    GLdouble model[16];
    GLdouble proj[16];
    GLint view[4];
    glGetIntegerv(GL_VIEWPORT, view);
    glGetDoublev(GL_PROJECTION_MATRIX, proj);
    glGetDoublev(GL_MODELVIEW_MATRIX, model);
    gluUnProject(x, y, z, model, proj, view, objX, objY, objZ);
}

```
Getting all the matrices, we can use this function to transform input coordinates to world coordinates.
After that, make mass point in the world space move a little bit along the direction.