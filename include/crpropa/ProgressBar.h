#ifndef CRPROPA_PROGRESSBAR_H
#define CRPROPA_PROGRESSBAR_H

#include <string>
#include <ctime>

namespace crpropa {

/**
 * \addtogroup Core
 * @{
 */

/** 
 @class ProgressBar
 @brief Track the evolution of the simulations with a progress bar
 */
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
	/** Constructor to initialize a progress bar
	 @param steps		number of steps
	 @param updateSteps	progress bar will be updated at the steps given by this parameter
	 */
	ProgressBar(unsigned long steps = 0, unsigned long updateSteps = 100);
	void start(const std::string &title);

	/** Update the progressbar
	 This should be called steps times in a loop.
	*/
	void update(); 

	/** Sets the position of the progress bar to a given value
	 @param position	current position of the progress bar
	 */
	void setPosition(unsigned long position);

	/** Mark the progress bar with an error
	 */
	void setError();
};
/** @}*/

} // namespace crpropa

#endif // CRPROPA_PROGRESSBAR_H
