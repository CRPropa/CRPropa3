#ifndef CRPROPA_PROGRESSBAR_H
#define CRPROPA_PROGRESSBAR_H

#include <stdio.h>
#include <string>
#include <iostream>
#include <ctime>

namespace crpropa {

class ProgressBar {
private:
	unsigned long _steps;
	unsigned long _currentCount;
	unsigned long _maxbarLength;
	unsigned long _nextStep;
	unsigned long _updateSteps;
	time_t _startTime;
	std::string stringTmpl;
	std::string arrow;

public:
	/// Initialize a ProgressBar with [steps] number of steps, updated at [updateSteps] intervalls
	ProgressBar(unsigned long steps = 0, unsigned long updateSteps = 100);
	void start(const std::string &title);

	/// update the progressbar
	/// should be called steps times in a loop 
	void update(); 

	// sets the position of the bar to a given value
	void setPosition(unsigned long position);

	/// Mark the progressbar with an error
	void setError();
};

} // namespace crpropa

#endif // CRPROPA_PROGRESSBAR_H
