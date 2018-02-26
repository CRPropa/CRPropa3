#ifndef CRPROPA_CLOCK_H
#define CRPROPA_CLOCK_H

namespace crpropa {
/**
 * \addtogroup Tools
 * @{
 */

/**
 @class Clock
 @brief Provides a clock / timer.

 */
class Clock {
private:
	class Impl;
	Impl *impl;
public:
	Clock();
	virtual ~Clock();

	void reset();
	double getSecond();
	double getMillisecond();
	static Clock &getInstance();
};
/** @} */

} // namespace crpropa

#endif // CRPROPA_CLOCK_H
