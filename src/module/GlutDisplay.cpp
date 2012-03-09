#include <GL/glut.h>
#include <GL/freeglut.h>
#include <iostream>
#include <stdlib.h>
#include <vector>

#include "mpc/module/GlutDisplay.h"
#include "mpc/ParticleState.h"

namespace mpc {

std::vector<ParticleState> points;

// Viewer state
int refreshAfter = 5;
float cameraPhi = 45.0;
float cameraTheta = 90.0;
float cameraDepth = 10;
bool leftButton = false;
bool middleButton = false;
int downX;
int downY;
float windowWidth = 500;
float windowHeight = 500;
bool continueDisplay;

void keyboard(unsigned char key, int x, int y) {
	switch (key) {
	case 13: // enter
		continueDisplay = 0;
		break;
	case 27: // escape
		break;
	case 32: // space
		continueDisplay = 0;
		break;
	}
}

void mouseButton(int button, int state, int x, int y) {
	// mouse position at moment of button press
	downX = x;
	downY = y;
	// mouse buttons
	leftButton = ((button == 0) && (state == 0));
	middleButton = ((button == 1) && (state == 0));
	// scroll wheel
	if (button == 3)
		cameraDepth -= 0.5;
	if (button == 4)
		cameraDepth += 0.5;
	glutPostRedisplay();
}

void mouseMove(int x, int y) {
	// rotate
	if (leftButton) {
		cameraPhi += (float) (x - downX) / 4.0;
		cameraTheta += (float) (downY - y) / 4.0;
	}
	// zoom
	if (middleButton) {
		cameraDepth += (float) (downY - y) / 10.0;
	}
	downX = x;
	downY = y;
	glutPostRedisplay();
}

void drawTrajectoryPoints() {
	glColor3f(1., 1., 1.);
	for (std::vector<ParticleState>::size_type i = 0; i < points.size(); i++) {
		Vector3 pos = points[i].getPosition() / Mpc;
		glPushMatrix();
		glTranslatef(pos.x(), pos.y(), pos.z());
		glutSolidSphere(0.015, 10, 5);
		glPopMatrix();
	}
}

void display(void) {
	glLoadIdentity();
	glTranslatef(0, 0, -cameraDepth);
	glRotatef(-cameraTheta, 1.0, 0.0, 0.0);
	glRotatef(cameraPhi, 0.0, 0.0, 1.0);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	drawTrajectoryPoints();
	glutSwapBuffers();
}

void reshape(int width, int height) {
	windowWidth = width;
	windowHeight = height;
	glViewport(0, 0, width, height);
}

//int displayMenu, mainMenu;
//enum {
//	TIME, ENERGY, CHARGE
//};
//int displayMode = TIME;
//
//void setMainMenu(int value) {
//	switch (value) {
//	case 99:
//		exit(0);
//		break;
//	}
//}
//
//void setDisplayMenu(int value) {
//	displayMode = value;
//	switch (value) {
//	case TIME:
//		break;
//	case ENERGY:
//		break;
//	case CHARGE:
//		break;
//	}
//	glutPostRedisplay();
//}

void initGlutDisplay() {
	int argc = 0;
	char *argv[1] = { 0 };
	// initialize GLUT and create window
	glutInit(&argc, argv);
	glutInitWindowSize(windowWidth, windowHeight);
	glutCreateWindow("Trajectory Display");
	// register callback functions
	glutDisplayFunc(display);
	glutReshapeFunc(reshape);
	glutKeyboardFunc(keyboard);
	glutMouseFunc(mouseButton);
	glutMotionFunc(mouseMove);
	// camera view
	glMatrixMode(GL_PROJECTION);
	gluPerspective(64.0, 1., 0., 100.);
	glMatrixMode(GL_MODELVIEW);
//	// menu
//	displayMenu = glutCreateMenu(setDisplayMenu);
//	glutAddMenuEntry("Wireframe", WIREFRAME);
//	glutAddMenuEntry("Hidden Line", HIDDENLINE);
//	glutAddMenuEntry("Flat Shaded", FLATSHADED);
//	glutAddMenuEntry("Smooth Shaded", SMOOTHSHADED);
//	mainMenu = glutCreateMenu(setMainMenu);
//	glutAddSubMenu("Display", displayMenu);
//	glutAddMenuEntry("Exit", 99);
//	glutAttachMenu(GLUT_RIGHT_BUTTON);
}

GlutDisplay::GlutDisplay() {
	counter = 0;
	initGlutDisplay();
}

GlutDisplay::~GlutDisplay() {
	glutDestroyWindow(glutGetWindow());
}

void GlutDisplay::process(Candidate *candidate) const {

	if (counter % refreshAfter == 0) {
		// append trajectory point
		points.push_back(candidate->current);

		glutPostRedisplay();

		// wait for user to continue
		continueDisplay = 1;
		while (continueDisplay == 1) {
			glutMainLoopEvent();
		}
	}
	return;
}

std::string GlutDisplay::getDescription() const {
	return "GlutDisplay";
}

} // namspace mpc
