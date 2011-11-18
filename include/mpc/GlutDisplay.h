/*
 * GlutDisplay.h
 *
 *  Created on: Nov 18, 2011
 *      Author: walz
 */

#ifndef GLUTDISPLAY_H_
#define GLUTDISPLAY_H_

#include <GL/glut.h>
#include <GL/freeglut.h>
#include <iostream>
#include <stdlib.h>
#include <vector>
#include "mpc/Candidate.h"

namespace mpc {

std::vector<Particle> states;

void display() {
	// clear color and depth buffers
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	// reset transformations
	glLoadIdentity();
	// set the camera
	gluLookAt(0.f, 0.f, 5.f, 0.f, 0.f, 0.f, 0.f, 1.f,  0.f);
	// draw trajectory points
	glColor3f(1.f, 1.f, 1.f); // white
	for (std::vector<Particle>::size_type i = 0; i < states.size(); i++) {
		Hep3Vector pos = states[i].getPosition() / Mpc;
		glPushMatrix();
		glTranslatef(pos.x(), pos.y(), pos.z());
		glutSolidSphere(0.015f, 8, 4);
		glPopMatrix();
	}
	// Swap hidden and visible buffer
	glutSwapBuffers();
}

void changeSize(int w1,int h1) {
	float ratio;
	ratio = 1.0f * w1 / h1;
	// reset the coordinate system before modifying
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	// set the viewport to be the entire window
    glViewport(0, 0, w1, h1);
	// set the clipping volume
	gluPerspective(45,ratio,0.1,1000);
	// always want to be in MODELVIEW mode
	glMatrixMode(GL_MODELVIEW);
}

void wait_for_key() {
    std::cout << "Press ENTER to continue...";
    std::cin.clear();
    std::cin.ignore(std::cin.rdbuf()->in_avail());
    std::cin.get();
    return;
}

class GlutDisplay {
public:
	int counter;
	int refresh;

	GlutDisplay(double r) {
		counter = 0;
		refresh = r;
		// set up GLUT
		int argc = 0;
		char *argv[1] = { 0 };
		// initialize GLUT and create window
		glutInit(&argc, argv);
		glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGBA);
		glutInitWindowPosition(50, 50);
		glutInitWindowSize(600, 600);
		glutCreateWindow("Trajectory Display");
		//
		glEnable(GL_DEPTH_TEST);
		// register callback functions
		glutDisplayFunc(display);
		glutReshapeFunc(changeSize);
	}

	void apply(Candidate &candidate) {
		if (counter % refresh == 0) {
			states.push_back(candidate.current);
			std::cout << candidate.current.getPosition() << std::endl;
			// mark the current window to be redisplayed
			glutPostRedisplay();

			glutMainLoopEvent();
			wait_for_key();
		}
		counter += 1;
		return;
	}
};

} // namspace mpc

#endif /* GLUTDISPLAY_H_ */







