#ifndef KISS_UUID_H
#define KISS_UUID_H

#include <iostream>

namespace kiss {

class uuid {
private:
	unsigned int ints[4]; /// Storage for the actual 16-digit ID.

public:
	/// Constructor, creates a new uuid.
	uuid();
	uuid(unsigned int *p);
	uuid(unsigned int a, unsigned int b, unsigned int c, unsigned int d);

	/// Constructor from a string \p id.

	/// Sets the ID to 0 (i.e. all components to 0).
	void reset();

	/// check if all zero
	bool valid() const;

	/// Create a new uuid. Standard way of creating a new uuid.
	static uuid create();

	bool operator ==(const uuid& id) const;
	bool operator !=(const uuid& id) const;

	/// Less-than-operator, provides ordering of IDs.
	bool operator <(const uuid& op) const;

	/// Returns a string representation of the ID.
	//std::string str() const;
	void print(char *p) const;

	static uuid parse(const char *id);
	static uuid parse(const std::string& id);

	unsigned int &operator[](size_t i);
	const unsigned int &operator[](size_t i) const;
	const unsigned char *data() const;
};

} // namespace kiss

std::ostream& operator <<(std::ostream& os, const kiss::uuid &id);
std::istream& operator >>(std::istream& os, kiss::uuid &id);

#endif // KISS_UUID_H
