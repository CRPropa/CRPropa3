#include "crpropa/ProgressBar.h"

namespace crpropa {

/// Initialize a ProgressBar with [steps] number of steps, updated at [updateSteps] intervalls
ProgressBar::ProgressBar(unsigned long steps, unsigned long updateSteps) :
		_steps(steps), _currentCount(0), _maxbarLength(10), _updateSteps(
				updateSteps), _nextStep(1), _startTime(0) {
	if (_updateSteps > _steps)
		_updateSteps = _steps;
	arrow.append(">");
}

void ProgressBar::start(const std::string &title) {
	_startTime = time(NULL);
	std::string s = ctime(&_startTime);
	s.erase(s.end() - 1, s.end());
	stringTmpl = "  Started ";
	stringTmpl.append(s);
	stringTmpl.append(" : [%-10s] %3i%%    %s: %02i:%02i:%02i %s\r");
	std::cout << title << std::endl;

}
/// update the progressbar
/// should be called steps times in a loop
void ProgressBar::update() {
	_currentCount++;
	if (_currentCount == _nextStep || _currentCount == _steps
			|| _currentCount == 1000) {
				_nextStep += long(_steps / float(_updateSteps));
			setPosition(_currentCount);
			}
}

void ProgressBar::setPosition(unsigned long position) {
	int percentage = int(100 * (position / float(_steps)));
	time_t currentTime = time(NULL);
	if (position < _steps) {
		int j = 0;
		if (arrow.size() <= (_maxbarLength) * (position) / (_steps))
			arrow.insert(0, "=");
		float tElapsed = currentTime - _startTime;
		float tToGo = (_steps - position) * tElapsed / position;
		printf(stringTmpl.c_str(), arrow.c_str(), percentage, "Finish in",
				int(tToGo / 3600), (int(tToGo) % 3600) / 60,
				int(tToGo) % 60, "");
		fflush(stdout);
	} else {
		float tElapsed = currentTime - _startTime;
		std::string s = " - Finished at ";
		s.append(ctime(&currentTime));
		char fs[255];
		sprintf(fs, "%c[%d;%dm Finished %c[%dm", 27, 1, 32, 27, 0);
		printf(stringTmpl.c_str(), fs, 100, "Needed",
				int(tElapsed / 3600), (int(tElapsed) % 3600) / 60,
				int(tElapsed) % 60, s.c_str());
	}
}


/// Mark the progressbar with an error
void ProgressBar::setError() {
	time_t currentTime = time(NULL);
	_currentCount++;
	float tElapsed = currentTime - _startTime;
	std::string s = " - Finished at ";
	s.append(ctime(&currentTime));
	char fs[255];
	sprintf(fs, "%c[%d;%dm  ERROR   %c[%dm", 27, 1, 31, 27, 0);
	printf(stringTmpl.c_str(), fs, _currentCount, "Needed",
			int(tElapsed / 3600), (int(tElapsed) % 3600) / 60,
			int(tElapsed) % 60, s.c_str());
}

} // namespace crpropa
