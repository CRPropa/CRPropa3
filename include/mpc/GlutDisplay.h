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

/* Viewer state */
float sphi = 90.0, stheta = 45.0, sdepth = 10;
bool leftButton = false, middleButton = false;
int downX, downY;
float windowWidth = 500, windowHeight = 500;
bool continueDisplay;

/* keyboard */
void pressKey(unsigned char key, int x, int y) {
	switch (key) {
	case 27: // escape
		continueDisplay = 0;
		break;
	}
}

/* mouse */
void mouseButton(int button, int state, int x, int y) {
	downX = x;
	downY = y;
	leftButton = ((button == GLUT_LEFT_BUTTON) && (state == GLUT_DOWN));
	middleButton = ((button == GLUT_MIDDLE_BUTTON) && (state == GLUT_DOWN));
	glutPostRedisplay();
}

void mouseMove(int x, int y) {
	if (leftButton) {
		sphi += (float) (x - downX) / 4.0;
		stheta += (float) (downY - y) / 4.0;
	} // rotate
	if (middleButton) {
		sdepth += (float) (downY - y) / 10.0;
	} // scale
	downX = x;
	downY = y;
	glutPostRedisplay();
}

void drawTrajectoryPoints() {
	glColor3f(1., 1., 1.);
	for (std::vector<Particle>::size_type i = 0; i < states.size(); i++) {
		Hep3Vector pos = states[i].getPosition() / Mpc;
		glPushMatrix();
		glTranslatef(pos.x(), pos.y(), pos.z());
		glutSolidSphere(0.015, 10, 5);
		glPopMatrix();
	}
}

void displayCallback(void) {
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	// gluPerspective(64.0, aspect, zNear, zFar);
	gluPerspective(64.0, 1.f, 1.f, 100.f);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	glTranslatef(0.0, 0.0, -sdepth);
	glRotatef(-stheta, 1.0, 0.0, 0.0);
	glRotatef(sphi, 0.0, 0.0, 1.0);
	glutSwapBuffers();
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	drawTrajectoryPoints();
	glutSwapBuffers();
}

void reshapeCallback(int width, int height) {
	windowWidth = width;
	windowHeight = height;
	glViewport(0, 0, width, height);
	gluPerspective(45, 1.f, 0.1, 1000);
}

void init() {
	int argc = 0;
	char *argv[1] = { 0 };
	// initialize GLUT and create window
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGBA);
	glutInitWindowPosition(50, 50);
	glutInitWindowSize(windowWidth, windowHeight);
	glutCreateWindow("Trajectory Display");
	glutSetOption(GLUT_ACTION_ON_WINDOW_CLOSE, GLUT_ACTION_CONTINUE_EXECUTION);
	glEnable(GL_DEPTH_TEST);
	// register callback functions
	glutDisplayFunc(displayCallback);
	glutReshapeFunc(reshapeCallback);
	glutIgnoreKeyRepeat(1);
	glutKeyboardFunc(pressKey);
	glutMouseFunc(mouseButton);
	glutMotionFunc(mouseMove);
}

class GlutDisplay {
public:
	int counter;
	int refresh;

	GlutDisplay(double r) {
		counter = 0;
		refresh = r;
		init();
	}

	void apply(Candidate &candidate) {
		if (counter % refresh == 0) {
			states.push_back(candidate.current);
			std::cout << candidate.current.getPosition() << std::endl;
			// mark the current window to be redisplayed
			glutPostRedisplay();
			continueDisplay = 1;
			while (continueDisplay == 1) {
				glutMainLoopEvent();
			}
		}
		counter += 1;
		return;
	}
};

} // namspace mpc

#endif /* GLUTDISPLAY_H_ */

