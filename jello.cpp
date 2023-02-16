/*

  USC/Viterbi/Computer Science
  "Jello Cube" Assignment 1 starter code

  Your name:
  Wanyu Zhang

*/

#include "jello.h"
#include "showCube.h"
#include "input.h"
#include "physics.h"
#include "iostream"
#include "ctime"

// camera parameters
double Theta = pi / 6;
double Phi = pi / 6;
double R = 6;

// mouse control
int g_iMenuId;
int g_vMousePos[2];
int g_iLeftMouseButton,g_iMiddleMouseButton,g_iRightMouseButton;

// number of images saved to disk so far
int sprite=0;

// these variables control what is displayed on screen
int shear=0, bend=0, structural=1, pause=0, viewingMode=0, saveScreenToFile=0;

// aspect ratio
double aspectRatio = 0;

point startMouse, endMouse;

struct world jello;

int windowWidth, windowHeight;

void myinit()
{
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  gluPerspective(90.0,1.0,0.01,1000.0);

  // set background color to grey
  glClearColor(0.5, 0.5, 0.5, 0.0);

  glCullFace(GL_BACK);
  glEnable(GL_CULL_FACE);

  glShadeModel(GL_SMOOTH);
  glEnable(GL_POLYGON_SMOOTH);
  glEnable(GL_LINE_SMOOTH);

  return; 
}

void reshape(int w, int h) 
{
  // Prevent a divide by zero, when h is zero.
  // You can't make a window of zero height.
  if(h == 0)
    h = 1;

  glViewport(0, 0, w, h);

  // Reset the coordinate system before modifying
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  // Set the perspective
  aspectRatio = 1.0 * w / h;
  gluPerspective(60.0f, aspectRatio, 0.01f, 1000.0f);

  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity(); 

  windowWidth = w;
  windowHeight = h;

  glutPostRedisplay();
}

void display()
{
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();

  // camera parameters are Phi, Theta, R
  gluLookAt(R * cos(Phi) * cos (Theta), R * sin(Phi) * cos (Theta), R * sin (Theta),
	        0.0,0.0,0.0, 0.0,0.0,1.0);


  /* Lighting */
  /* You are encouraged to change lighting parameters or make improvements/modifications
     to the lighting model . 
     This way, you will personalize your assignment and your assignment will stick out. 
  */

  // global ambient light
  GLfloat aGa[] = { 0.0, 0.0, 0.0, 0.0 };
  
  // light 's ambient, diffuse, specular
  GLfloat lKa0[] = { 0.0, 0.0, 0.0, 1.0 };
  GLfloat lKd0[] = { 1.0, 1.0, 1.0, 1.0 };
  GLfloat lKs0[] = { 1.0, 1.0, 1.0, 1.0 };

  GLfloat lKa1[] = { 0.0, 0.0, 0.0, 1.0 };
  GLfloat lKd1[] = { 1.0, 0.0, 0.0, 1.0 };
  GLfloat lKs1[] = { 1.0, 0.0, 0.0, 1.0 };

  GLfloat lKa2[] = { 0.0, 0.0, 0.0, 1.0 };
  GLfloat lKd2[] = { 1.0, 1.0, 0.0, 1.0 };
  GLfloat lKs2[] = { 1.0, 1.0, 0.0, 1.0 };

  GLfloat lKa3[] = { 0.0, 0.0, 0.0, 1.0 };
  GLfloat lKd3[] = { 0.0, 1.0, 1.0, 1.0 };
  GLfloat lKs3[] = { 0.0, 1.0, 1.0, 1.0 };

  GLfloat lKa4[] = { 0.0, 0.0, 0.0, 1.0 };
  GLfloat lKd4[] = { 0.0, 0.0, 1.0, 1.0 };
  GLfloat lKs4[] = { 0.0, 0.0, 1.0, 1.0 };

  GLfloat lKa5[] = { 0.0, 0.0, 0.0, 1.0 };
  GLfloat lKd5[] = { 1.0, 0.0, 1.0, 1.0 };
  GLfloat lKs5[] = { 1.0, 0.0, 1.0, 1.0 };

  GLfloat lKa6[] = { 0.0, 0.0, 0.0, 1.0 };
  GLfloat lKd6[] = { 1.0, 1.0, 1.0, 1.0 };
  GLfloat lKs6[] = { 1.0, 1.0, 1.0, 1.0 };

  GLfloat lKa7[] = { 0.0, 0.0, 0.0, 1.0 };
  GLfloat lKd7[] = { 0.0, 1.0, 1.0, 1.0 };
  GLfloat lKs7[] = { 0.0, 1.0, 1.0, 1.0 };

  // light positions and directions
  GLfloat lP0[] = { -1.999, -1.999, -1.999, 1.0 };
  GLfloat lP1[] = { 1.999, -1.999, -1.999, 1.0 };
  GLfloat lP2[] = { 1.999, 1.999, -1.999, 1.0 };
  GLfloat lP3[] = { -1.999, 1.999, -1.999, 1.0 };
  GLfloat lP4[] = { -1.999, -1.999, 1.999, 1.0 };
  GLfloat lP5[] = { 1.999, -1.999, 1.999, 1.0 };
  GLfloat lP6[] = { 1.999, 1.999, 1.999, 1.0 };
  GLfloat lP7[] = { -1.999, 1.999, 1.999, 1.0 };
  
  // jelly material color

  GLfloat mKa[] = { 0.0, 0.0, 0.0, 1.0 };
  GLfloat mKd[] = { 0.3, 0.3, 0.3, 1.0 };
  GLfloat mKs[] = { 1.0, 1.0, 1.0, 1.0 };
  GLfloat mKe[] = { 0.0, 0.0, 0.0, 1.0 };

  /* set up lighting */
  glLightModelfv(GL_LIGHT_MODEL_AMBIENT, aGa);
  glLightModelf(GL_LIGHT_MODEL_LOCAL_VIEWER, GL_TRUE);
  glLightModelf(GL_LIGHT_MODEL_TWO_SIDE, GL_FALSE);

  // set up cube color
  glMaterialfv(GL_FRONT, GL_AMBIENT, mKa);
  glMaterialfv(GL_FRONT, GL_DIFFUSE, mKd);
  glMaterialfv(GL_FRONT, GL_SPECULAR, mKs);
  glMaterialfv(GL_FRONT, GL_EMISSION, mKe);
  glMaterialf(GL_FRONT, GL_SHININESS, 120);
    
  // macro to set up light i
  #define LIGHTSETUP(i)\
  glLightfv(GL_LIGHT##i, GL_POSITION, lP##i);\
  glLightfv(GL_LIGHT##i, GL_AMBIENT, lKa##i);\
  glLightfv(GL_LIGHT##i, GL_DIFFUSE, lKd##i);\
  glLightfv(GL_LIGHT##i, GL_SPECULAR, lKs##i);\
  glEnable(GL_LIGHT##i)
  
  LIGHTSETUP (0);
  LIGHTSETUP (1);
  LIGHTSETUP (2);
  LIGHTSETUP (3);
  LIGHTSETUP (4);
  LIGHTSETUP (5);
  LIGHTSETUP (6);
  LIGHTSETUP (7);

  // enable lighting
  glEnable(GL_LIGHTING);    
  glEnable(GL_DEPTH_TEST);

  /*
   *    code starts here, test gl first
   *    hard-core modification test pass
   */
  if(pause == 0)
  RK4(&jello);

  // show the cube
  showCube(&jello);
  glFlush();

  glDisable(GL_LIGHTING);

  // show the bounding box
  showBoundingBox();

  showInclinedPlane();
 
  glutSwapBuffers();
}

double glReadDepth(double x,double y,double *per=NULL)                  // x,y [pixels], per[16]
{
    GLfloat _z=0.0; double m[16],z,zFar,zNear;
    if (per==NULL){ per=m; glGetDoublev(GL_PROJECTION_MATRIX,per); }    // use actual perspective matrix if not passed
    zFar =0.5*per[14]*(1.0-((per[10]-1.0)/(per[10]+1.0)));              // compute zFar from perspective matrix
    zNear=zFar*(per[10]+1.0)/(per[10]-1.0);                             // compute zNear from perspective matrix
    glReadPixels(x,y,1,1,GL_DEPTH_COMPONENT,GL_FLOAT,&_z);              // read depth value
    z=_z;                                                               // logarithmic
    z=(2.0*z)-1.0;                                                      // logarithmic NDC
    z=(2.0*zNear*zFar)/(zFar+zNear-(z*(zFar-zNear)));                   // linear <zNear,zFar>
    return -z;
}

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

void mouseButton(int button, int state, int x, int y)
{
    switch (button)
    {
        case GLUT_LEFT_BUTTON:
            g_iLeftMouseButton = (state==GLUT_DOWN);
            std::cout<<"x is "<<x<<std::endl;
            std::cout<<"y is "<<y<<std::endl;
            float z;
            glEnable (GL_DEPTH_TEST);
            glReadPixels(x, y, 1, 1, GL_DEPTH_COMPONENT, GL_FLOAT, &z);
            //z = glReadDepth(x, y);
            std::cout<<"z is "<<z<<std::endl;
            GLdouble model[16];
            GLdouble proj[16];
            GLint view[4];
            glGetIntegerv(GL_VIEWPORT, view);
            glGetDoublev(GL_PROJECTION_MATRIX, proj);
            glGetDoublev(GL_MODELVIEW_MATRIX, model);
            double objX, objY, objZ;
            gluUnProject(x, y, 0.5, model, proj, view, &objX, &objY, &objZ);
            startMouse = {objX, objY, objZ};
            std::cout<<"obj x is "<<objX<<std::endl;
            std::cout<<"obj y is "<<objY<<std::endl;
            std::cout<<"obj z is "<<objZ<<std::endl;
            pickPoint(x, y, &jello, aspectRatio);
            break;
        case GLUT_MIDDLE_BUTTON:
            g_iMiddleMouseButton = (state==GLUT_DOWN);
            break;
        case GLUT_RIGHT_BUTTON:
            g_iRightMouseButton = (state==GLUT_DOWN);
            break;
    }

    g_vMousePos[0] = x;
    g_vMousePos[1] = y;
}

/* converts mouse drags into information about rotation/translation/scaling */
void mouseMotionDrag(int x, int y)
{
    int vMouseDelta[2] = {x-g_vMousePos[0], y-g_vMousePos[1]};

    if (g_iRightMouseButton) // handle camera rotations
    {
        Phi += vMouseDelta[0] * 0.01;
        Theta += vMouseDelta[1] * 0.01;

        if (Phi>2*pi)
            Phi -= 2*pi;

        if (Phi<0)
            Phi += 2*pi;

        if (Theta>pi / 2 - 0.01) // dont let the point enter the north pole
            Theta = pi / 2 - 0.01;

        if (Theta<- pi / 2 + 0.01)
            Theta = -pi / 2 + 0.01;

        g_vMousePos[0] = x;
        g_vMousePos[1] = y;
    }
    else if(g_iLeftMouseButton)
    {
        double objX, objY, objZ;
        glUntransform(x, windowHeight - y, 1.0, &objX, &objY, &objZ);
        std::cout<<"new obj x is "<<objX<<std::endl;
        std::cout<<"new obj y is "<<objY<<std::endl;
        std::cout<<"new obj z is "<<objZ<<std::endl;
        endMouse = {objX, objY, objZ};
        point unit_vector = endMouse - startMouse;
        double unit_length = sqrt(unit_vector * unit_vector);
        unit_vector = 1/unit_length * unit_vector;
        double length = sqrt(vMouseDelta[0] * vMouseDelta[0] + vMouseDelta[1] * vMouseDelta[1]);
        // one deform function to call (vector of direction and magnitude (x,y)-())
        pullPoint(unit_vector, length, &jello);
    }
}

void doIdle()
{
  char s[20]="picxxxx.ppm";
  int i;
  
  // save screen to file
  s[3] = 48 + (sprite / 1000);
  s[4] = 48 + (sprite % 1000) / 100;
  s[5] = 48 + (sprite % 100 ) / 10;
  s[6] = 48 + sprite % 10;

  if (saveScreenToFile==1)
  {
    saveScreenshot(windowWidth, windowHeight, s);
    saveScreenToFile=0; // save only once, change this if you want continuos image generation (i.e. animation)
    sprite++;
  }

  if (sprite >= 300) // allow only 300 snapshots
  {
    exit(0);	
  }

  if (pause == 0)
  {
    // insert code which appropriately performs one step of the cube simulation:
  }

  glutPostRedisplay();
}

void computeInclinedPlane(world& jello)
{
    jello.incPlanePresent = 1;
    double x1 = 1.0, y1 = -2.0, z1 = -2.0;
    double x2 = -2.0, y2 = 1.0, z2 = -2.0;
    double x3 = -2.0, y3 = -2.0, z3 = 1.0;
    double a = (y2-y1)*(z3-z1)-(y3-y1)*(z2-z1);
    double b = (z2-z1)*(x3-x1)-(z3-z1)*(x2-x1);
    double c = (x2-x1)*(y3-y1)-(x3-x1)*(y2-y1);
    double d = - a * x1 - b * y1 - c * z1;
    jello.a = a;
    jello.b = b;
    jello.c = c;
    jello.d = d;
}

int main (int argc, char ** argv)
{
  if (argc<2)
  {  
    printf ("Oops! You didn't say the jello world file!\n");
    printf ("Usage: %s [worldfile]\n", argv[0]);
    exit(0);
  }

  readWorld(argv[1],&jello);

  //jello.dt = 0.0005;
  /*
   * Added by Wanyu Zhang
   * */
  std::cout<<jello.dt<<std::endl;
  std::cout<<jello.mass<<std::endl;
  std::cout<<jello.v[0][0][0].x<<std::endl;
  std::cout<<jello.v[0][0][0].y<<std::endl;
  std::cout<<jello.v[0][0][0].z<<std::endl;
  std::cout<<jello.resolution<<std::endl;
  /*
   * test operator overloading
   */
  point a = {1,2,3};
  point b = {2,3,4};
  std::cout<<(a - b).z<<std::endl;
  /*
   * call Inclined Plane function to compute a, b, c, d
   */
  computeInclinedPlane(jello);
  std::cout<<"a is"<<jello.a<<std::endl;
  std::cout<<"b is"<<jello.b<<std::endl;
  std::cout<<"c is"<<jello.c<<std::endl;
  std::cout<<"d is"<<jello.d<<std::endl;

  glutInit(&argc,argv);
  
  /* double buffered window, use depth testing, 640x480 */
  glutInitDisplayMode (GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
  
  windowWidth = 640;
  windowHeight = 480;
  glutInitWindowSize (windowWidth, windowHeight);
  glutInitWindowPosition (0,0);
  glutCreateWindow ("Jello cube");

  /* tells glut to use a particular display function to redraw */
  glutDisplayFunc(display);

  /* replace with any animate code */
  glutIdleFunc(doIdle);

  /* callback for mouse drags */
  glutMotionFunc(mouseMotionDrag);

  /* callback for window size changes */
  glutReshapeFunc(reshape);

  /* callback for mouse movement */
  glutPassiveMotionFunc(mouseMotion);

  /* callback for mouse button changes */
  glutMouseFunc(mouseButton);

  /* register for keyboard events */
  glutKeyboardFunc(keyboardFunc);

  /* do initialization */
  myinit();

  /* forever sink in the black hole */
  glutMainLoop();

  return(0);
}

