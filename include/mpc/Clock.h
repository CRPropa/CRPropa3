#ifndef CLOCK_H_
#define CLOCK_H_

namespace mpc {

//class ClockImpl;

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

} /* namespace scs */
#endif /* CLOCK_H_ */
