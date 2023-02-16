/*

  USC/Viterbi/Computer Science
  "Jello Cube" Assignment 1 starter code

*/


#ifndef _SHOWCUBE_H_
#define _SHOWCUBE_H_

void showCube(struct world * jello);

void showBoundingBox();

void showInclinedPlane();

void pickPoint(int x, int y, world* jello, double aspectRatio);

void pullPoint(point unit_vector, double length, world* jello);

#endif
