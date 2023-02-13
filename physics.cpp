/*

  USC/Viterbi/Computer Science
  "Jello Cube" Assignment 1 starter code

*/

#include "vector"
#include "jello.h"
#include "physics.h"

#define CUBE_MIN_INDEX 0
#define CUBE_MAX_INDEX 7
#define BOUNDING_BOX_LENGTH 4.0

double computeOriginalDistance(position& a, position& b)
{
    return sqrt((a.x - b.x) * (a.x - b.x)
                + (a.y - b.y) * (a.y - b.y)
                + (a.z - b.z) * (a.z - b.z))/7;
}

double computeCurrentDistance(point& a, point& b)
{
    return sqrt((a.x - b.x) * (a.x - b.x)
        + (a.y - b.y) * (a.y - b.y)
        + (a.z - b.z) * (a.z - b.z));
}

point computeHookLaw(double kHook, partitcle& a, partitcle &b)
{
    double L = computeCurrentDistance(a.p, b.p);
    double rest = computeOriginalDistance(a.pos, b.pos);
    point Lvector = a.p - b.p;
    double coefficient = - kHook * (L - rest) / L;
    point F = {.x = coefficient * Lvector.x, .y = coefficient * Lvector.y, .z = coefficient * Lvector.z};
    return F;
}

point computeDampingLaw(double kDamp, partitcle &a, partitcle& b)
{
    double L = computeCurrentDistance(a.p, b.p);
    point Lvector = a.p - b.p;
    double coefficient = - kDamp * ((a.v - b.v) * Lvector) / (L * L);
    point F = {.x = coefficient * Lvector.x, .y = coefficient * Lvector.y, .z = coefficient * Lvector.z};
    return F;
}

bool checkBoundary(int x, int y, int z, int i, int j, int k)
{
    // return true if inside boundary
    return (x + i >= CUBE_MIN_INDEX && x + i <= CUBE_MAX_INDEX)
        && (y + j >= CUBE_MIN_INDEX && y + j <= CUBE_MAX_INDEX)
        && (z + k >= CUBE_MIN_INDEX && z + k <= CUBE_MAX_INDEX);
}

bool checkOutsideBoundingBox(point* p)
{
    return (p->x<-2 || p->x>2) || (p->y<-2 || p->y>2) || (p->z<-2 || p->z>2);
}

std::vector<partitcle> computeSprings(world * jello, int x, int y, int z)
{
    std::vector<partitcle> springs;
    /*
     * compute structural, shear, bend springs
     */

    /* we can do structural and shear springs separately,
     * but a better method is to simply add all the available neighbours of the node into vector
     * since it either be a structural spring or a shear spring
     * For instance in a 3 X 3 grids, the center node has 6 structural nodes and 20 shear nodes,
     * adding them up plus one center node makes total 27 (3 X 3 X 3) nodes.
     */
    for(int i = -1; i < 2; i++)
    {
        for(int j = -1; j < 2; j++)
        {
            for(int k = -1; k < 2; k++)
            {
                if(i == 0 && j == 0 && k == 0)
                {
                    //itself skip
                    continue;
                }
                if(checkBoundary(x, y, z, i, j, k))
                {
                    position pos = {.x = x+i, .y = y+j, .z = z+k};
                    partitcle par = {.p = jello->p[x+i][y+j][z+k], .pos = pos, .v = jello->v[x+i][y+j][z+k]};
                    springs.push_back(par);
                }
            }
        }
    }

    //Now add bend springs, only 6 directions
    if(checkBoundary(x, y, z, -2, 0, 0)){
        position pos = {.x = x-2, .y = y, .z = z};
        partitcle par = {.p = jello->p[x-2][y][z], .pos = pos, .v = jello->v[x-2][y][z]};
        springs.push_back(par);
    }
    if(checkBoundary(x, y, z, 2, 0, 0)){
        position pos = {.x = x+2, .y = y, .z = z};
        partitcle par = {.p = jello->p[x+2][y][z], .pos = pos, .v = jello->v[x+2][y][z]};
        springs.push_back(par);
    }
    if(checkBoundary(x, y, z, 0, -2, 0)){
        position pos = {.x = x, .y = y-2, .z = z};
        partitcle par = {.p = jello->p[x][y-2][z], .pos = pos, .v = jello->v[x][y-2][z]};
        springs.push_back(par);
    }
    if(checkBoundary(x, y, z, 0, 2, 0)){
        position pos = {.x = x, .y = y+2, .z = z};
        partitcle par = {.p = jello->p[x][y+2][z], .pos = pos, .v = jello->v[x][y+2][z]};
        springs.push_back(par);
    }
    if(checkBoundary(x, y, z, 0, 0, -2)){
        position pos = {.x = x, .y = y, .z = z-2};
        partitcle par = {.p = jello->p[x][y][z-2], .pos = pos, .v = jello->v[x][y][z-2]};
        springs.push_back(par);
    }
    if(checkBoundary(x, y, z, 0, 0, 2)){
        position pos = {.x = x, .y = y, .z = z+2};
        partitcle par = {.p = jello->p[x][y][z+2], .pos = pos, .v = jello->v[x][y][z+2]};
        springs.push_back(par);
    }

    return springs;
}

std::vector<partitcle> computeCollisionSprings(world * jello, int x, int y, int z)
{
    std::vector<partitcle>  collisions;
    if(checkOutsideBoundingBox(&jello->p[x][y][z]))
    {
        if(jello->p[x][y][z].x < -2 || jello->p[x][y][z].x > 2)
        {
            partitcle xc = {.p = {.x = jello->p[x][y][z].x<-2.0?-2.0:2.0, .y = jello->p[x][y][z].y, .z = jello->p[x][y][z].z}, .pos = {.x = x, .y = y, .z = z}, .v = {0,0,0}};
            collisions.push_back(xc);
        }
        if(jello->p[x][y][z].y < -2 || jello->p[x][y][z].y > 2)
        {
            partitcle yc = {.p = {.x = jello->p[x][y][z].x, .y = jello->p[x][y][z].y<-2.0?-2.0:2.0, .z = jello->p[x][y][z].z}, .pos = {.x = x, .y = y, .z = z}, .v = {0,0,0}};
            collisions.push_back(yc);
        }
        if(jello->p[x][y][z].z < -2 || jello->p[x][y][z].z > 2)
        {
            partitcle zc = {.p = {.x = jello->p[x][y][z].x, .y = jello->p[x][y][z].y, .z = jello->p[x][y][z].z<-2.0?-2.0:2.0}, .pos = {.x = x, .y = y, .z = z}, .v = {0,0,0}};
            collisions.push_back(zc);
        }
    }
    return collisions;
}

int computeInterval(double p, double interval, double& proportion)
{
    if(p < -2.0 || p > 2.0 || isnan(p))
        return -1;
    int index = (int)((p+2.0)/interval);
    proportion = ((p+2.0)-interval*index)/interval;
    return index;
}


std::vector<partitcle> computeInclinedCollision(world * jello, int x, int y, int z)
{
    std::vector<partitcle> inclinedCollisions;
    double x1 = jello->p[x][y][z].x;
    double y1 = jello->p[x][y][z].y;
    double z1 = jello->p[x][y][z].z;
    if(jello->a * x1 + jello->b * y1 + jello->c * z1 + jello->d < 0)
    {
        double t = (jello->a * x1 + jello->b * y1 + jello->c * z1 + jello->d)
                /(jello->a * jello->a + jello->b * jello->b + jello->c * jello->c);
        partitcle ic = {.p = {.x = x1 - t * jello->a, .y = y1 - t * jello->b, .z = z1 - t * jello->c}, .pos = {x, y, z}, .v = {0,0,0}};
        inclinedCollisions.push_back(ic);
    }

    return inclinedCollisions;
}

/* Computes acceleration to every control point of the jello cube, 
   which is in state given by 'jello'.
   Returns result in array 'a'. */
void computeAcceleration(struct world * jello, struct point a[8][8][8])
{
      /* for you to implement ... */
      /*
       * Figure out how many points are connected to current point
       * Apply Hook's Law and damping law in 3D to every pair of spring
       * */
      std::vector<partitcle> springs;
      std::vector<partitcle> collisions;
      for(int x = 0; x < 8; x++)
      {
          for(int y = 0; y < 8; y++)
          {
              for(int z = 0; z < 8; z++)
              {
                  point F = {.x=0, .y=0, .z=0};
                  springs = computeSprings(jello, x, y, z);
                  partitcle cur = {.p = jello->p[x][y][z], .pos = {.x = x, .y = y, .z = z}, .v = jello->v[x][y][z]};
                  for(auto & p : springs)
                  {
                        point Fhook = computeHookLaw(200, cur, p);
                        point Fdamp = computeDampingLaw(0.2, cur, p);
                        F = F + Fhook + Fdamp;
                  }
                  // compute collisions
                  collisions = computeCollisionSprings(jello, x, y, z);
                  for(auto& p : collisions)
                  {
                      point FCollisionHook = computeHookLaw(2000, cur, p);
                      point FCollisionDamp = computeDampingLaw(0.2, cur, p);
                      F = F + FCollisionHook + FCollisionDamp;
                  }
                  // compute collisions with inclined plane
                  if(jello->incPlanePresent == 1)
                  {
                      std::vector<partitcle> inclinedCollisions = computeInclinedCollision(jello, x, y, z);
                      for(auto& p : inclinedCollisions)
                      {
                          point FCollisionHook = computeHookLaw(2000, cur, p);
                          point FCollisionDamp = computeDampingLaw(0.2, cur, p);
                          F = F + FCollisionHook + FCollisionDamp;
                      }
                  }
                  // interpolation
                  if(jello->resolution > 0)
                  {
                      double interval = BOUNDING_BOX_LENGTH / (jello->resolution-1);
                      int res = jello->resolution;
                      double xp, yp, zp;
                      int xInterval = computeInterval(jello->p[x][y][z].x, interval, xp);
                      int yInterval = computeInterval(jello->p[x][y][z].y, interval, yp);
                      int zInterval = computeInterval(jello->p[x][y][z].z, interval, zp);
                      if(xInterval >= 0 && xInterval < res
                      && yInterval >= 0 && yInterval < res
                      && zInterval >= 0 && zInterval < res)
                      {
                          /*F.x = F.x + (1.0-xp) * jello->forceField[xInterval*res*res+yInterval*res+zInterval].x
                                  + xp * jello->forceField[(xInterval+1)*res*res+yInterval*res+zInterval].x;
                          F.y = F.y + (1.0-yp) * jello->forceField[xInterval*res*res+yInterval*res+zInterval].y
                                  + yp * jello->forceField[xInterval*res*res+(yInterval+1)*res+zInterval].y;
                          F.z = F.z + (1.0-zp) * jello->forceField[xInterval*res*res+yInterval*res+zInterval].z
                                  + zp * jello->forceField[xInterval*res*res+yInterval*res+(zInterval+1)].z;*/
                          F = F + (1.0-xp)*(1.0-yp)*(1.0-zp)*jello->forceField[xInterval*res*res+yInterval*res+zInterval]
                                + xp*(1.0-yp)*(1.0-zp)*jello->forceField[(xInterval+1)*res*res+yInterval*res+zInterval]
                                + (1.0-xp)*yp*(1.0-zp)*jello->forceField[xInterval*res*res+(yInterval+1)*res+zInterval]
                                + xp*yp*(1.0-zp)*jello->forceField[(xInterval+1)*res*res+(yInterval+1)*res+zInterval]
                                + (1.0-xp)*(1.0-yp)*zp*jello->forceField[xInterval*res*res+yInterval*res+(zInterval+1)]
                                + xp*(1.0-yp)*zp*jello->forceField[(xInterval+1)*res*res+yInterval*res+(zInterval+1)]
                                + (1.0-xp)*yp*zp*jello->forceField[xInterval*res*res+(yInterval+1)*res+(zInterval+1)]
                                + xp*yp*zp*jello->forceField[(xInterval+1)*res*res+(yInterval+1)*res+(zInterval+1)];
                      }
                  }

                  a[x][y][z].x = F.x/jello->mass;
                  a[x][y][z].y = F.y/jello->mass;
                  a[x][y][z].z = F.z/jello->mass;
              }
          }
      }
}

/* performs one step of Euler Integration */
/* as a result, updates the jello structure */
void Euler(struct world * jello)
{
  int i,j,k;
  point a[8][8][8];

  computeAcceleration(jello, a);
  
  for (i=0; i<=7; i++)
    for (j=0; j<=7; j++)
      for (k=0; k<=7; k++)
      {
        jello->p[i][j][k].x += jello->dt * jello->v[i][j][k].x;
        jello->p[i][j][k].y += jello->dt * jello->v[i][j][k].y;
        jello->p[i][j][k].z += jello->dt * jello->v[i][j][k].z;
        jello->v[i][j][k].x += jello->dt * a[i][j][k].x;
        jello->v[i][j][k].y += jello->dt * a[i][j][k].y;
        jello->v[i][j][k].z += jello->dt * a[i][j][k].z;

      }
}

/* performs one step of RK4 Integration */
/* as a result, updates the jello structure */
void RK4(struct world * jello)
{
  point F1p[8][8][8], F1v[8][8][8], 
        F2p[8][8][8], F2v[8][8][8],
        F3p[8][8][8], F3v[8][8][8],
        F4p[8][8][8], F4v[8][8][8];

  point a[8][8][8];


  struct world buffer;

  int i,j,k;

  buffer = *jello; // make a copy of jello

  computeAcceleration(jello, a);

  for (i=0; i<=7; i++)
    for (j=0; j<=7; j++)
      for (k=0; k<=7; k++)
      {
         pMULTIPLY(jello->v[i][j][k],jello->dt,F1p[i][j][k]);
         pMULTIPLY(a[i][j][k],jello->dt,F1v[i][j][k]);
         pMULTIPLY(F1p[i][j][k],0.5,buffer.p[i][j][k]);
         pMULTIPLY(F1v[i][j][k],0.5,buffer.v[i][j][k]);
         pSUM(jello->p[i][j][k],buffer.p[i][j][k],buffer.p[i][j][k]);
         pSUM(jello->v[i][j][k],buffer.v[i][j][k],buffer.v[i][j][k]);
      }

  computeAcceleration(&buffer, a);

  for (i=0; i<=7; i++)
    for (j=0; j<=7; j++)
      for (k=0; k<=7; k++)
      {
         // F2p = dt * buffer.v;
         pMULTIPLY(buffer.v[i][j][k],jello->dt,F2p[i][j][k]);
         // F2v = dt * a(buffer.p,buffer.v);     
         pMULTIPLY(a[i][j][k],jello->dt,F2v[i][j][k]);
         pMULTIPLY(F2p[i][j][k],0.5,buffer.p[i][j][k]);
         pMULTIPLY(F2v[i][j][k],0.5,buffer.v[i][j][k]);
         pSUM(jello->p[i][j][k],buffer.p[i][j][k],buffer.p[i][j][k]);
         pSUM(jello->v[i][j][k],buffer.v[i][j][k],buffer.v[i][j][k]);
      }

  computeAcceleration(&buffer, a);

  for (i=0; i<=7; i++)
    for (j=0; j<=7; j++)
      for (k=0; k<=7; k++)
      {
         // F3p = dt * buffer.v;
         pMULTIPLY(buffer.v[i][j][k],jello->dt,F3p[i][j][k]);
         // F3v = dt * a(buffer.p,buffer.v);     
         pMULTIPLY(a[i][j][k],jello->dt,F3v[i][j][k]);
         pMULTIPLY(F3p[i][j][k],1.0,buffer.p[i][j][k]);
         pMULTIPLY(F3v[i][j][k],1.0,buffer.v[i][j][k]);
         pSUM(jello->p[i][j][k],buffer.p[i][j][k],buffer.p[i][j][k]);
         pSUM(jello->v[i][j][k],buffer.v[i][j][k],buffer.v[i][j][k]);
      }
         
  computeAcceleration(&buffer, a);


  for (i=0; i<=7; i++)
    for (j=0; j<=7; j++)
      for (k=0; k<=7; k++)
      {
         // F3p = dt * buffer.v;
         pMULTIPLY(buffer.v[i][j][k],jello->dt,F4p[i][j][k]);
         // F3v = dt * a(buffer.p,buffer.v);     
         pMULTIPLY(a[i][j][k],jello->dt,F4v[i][j][k]);

         pMULTIPLY(F2p[i][j][k],2,buffer.p[i][j][k]);
         pMULTIPLY(F3p[i][j][k],2,buffer.v[i][j][k]);
         pSUM(buffer.p[i][j][k],buffer.v[i][j][k],buffer.p[i][j][k]);
         pSUM(buffer.p[i][j][k],F1p[i][j][k],buffer.p[i][j][k]);
         pSUM(buffer.p[i][j][k],F4p[i][j][k],buffer.p[i][j][k]);
         pMULTIPLY(buffer.p[i][j][k],1.0 / 6,buffer.p[i][j][k]);
         pSUM(buffer.p[i][j][k],jello->p[i][j][k],jello->p[i][j][k]);

         pMULTIPLY(F2v[i][j][k],2,buffer.p[i][j][k]);
         pMULTIPLY(F3v[i][j][k],2,buffer.v[i][j][k]);
         pSUM(buffer.p[i][j][k],buffer.v[i][j][k],buffer.p[i][j][k]);
         pSUM(buffer.p[i][j][k],F1v[i][j][k],buffer.p[i][j][k]);
         pSUM(buffer.p[i][j][k],F4v[i][j][k],buffer.p[i][j][k]);
         pMULTIPLY(buffer.p[i][j][k],1.0 / 6,buffer.p[i][j][k]);
         pSUM(buffer.p[i][j][k],jello->v[i][j][k],jello->v[i][j][k]);
      }

  return;  
}
