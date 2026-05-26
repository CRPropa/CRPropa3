#ifndef CRPROPA_REFERENCED_H
#define CRPROPA_REFERENCED_H

#include <cstddef>
#include <memory>

namespace crpropa {
/**
 * \addtogroup Core
 * @{
 */

/**
 @class ref_ptr
@brief Referenced pointer
*/
template<class T>
class ref_ptr {
	public:
	typedef T element_type;

	ref_ptr() : _raw_ptr(NULL), _shared_ptr(NULL) {}
	ref_ptr(T& obj) {
		_raw_ptr = &obj;
		_shared_ptr = NULL;
	}
	ref_ptr(T* ptr) {
		_raw_ptr = NULL;
		_shared_ptr = std::shared_ptr<T>(ptr);
	}
	ref_ptr(const ref_ptr& rp) {
		_raw_ptr = rp._raw_ptr;
		_shared_ptr = rp._shared_ptr;
	}
	template<class Other> ref_ptr(const ref_ptr<Other>& rp) {
		_raw_ptr = rp._raw_ptr;
		_shared_ptr = rp._shared_ptr;
	}
	template<class Other> ref_ptr(const std::shared_ptr<Other>& shared_ptr) {
		_raw_ptr = NULL;
		_shared_ptr = shared_ptr;
	}

	~ref_ptr() {
		_raw_ptr = NULL;  // do not delete, it is expected to be managed by user
		_shared_ptr = NULL;
	}

	ref_ptr& operator =(const ref_ptr& rp) {
		assign(rp);
		return *this;
	}

	template<class Other> ref_ptr& operator =(const ref_ptr<Other>& rp) {
		assign(rp);
		return *this;
	}

	inline ref_ptr& operator =(long int ptr) {
		_raw_ptr = NULL;
		_shared_ptr = NULL;
		return *this;
	}

	operator T*() const {
		return get();
	}

	T& operator*() const {
		if (_raw_ptr) return *_raw_ptr;
		return *_shared_ptr;
	}
	T* operator->() const {
		if (_raw_ptr) return _raw_ptr;
		return _shared_ptr.get();
	}

	// bool operator==(const ref_ptr& rp) const {
	// 	return get()==rp.get();
	// }

	T* get() const {
		if (_raw_ptr) return _raw_ptr;
		return _shared_ptr.get();
	}

	std::shared_ptr<T> get_shared() const{
		return _shared_ptr;
	}

	bool valid() const {
		if (_raw_ptr) return true;
		return !(_shared_ptr==NULL);
	}

	void release() {
		_raw_ptr = NULL;  // do not delete, it is expected to be managed by user
		_shared_ptr = NULL;
	}

	void swap(ref_ptr& rp) {
		if (_raw_ptr){
			T* tmp = _raw_ptr;
			_raw_ptr = rp._raw_ptr;
			rp._raw_ptr = tmp;
		} else {
			_shared_ptr.swap(rp._shared_ptr);
		}
	}

	private:

	template<class Other> void assign(const ref_ptr<Other>& rp) {
		_raw_ptr = rp._raw_ptr;
		_shared_ptr = rp._shared_ptr;
	}

	template<class Other> friend class ref_ptr;

	T* _raw_ptr = NULL;
	std::shared_ptr<T> _shared_ptr = NULL;
};

template<class T> inline
void swap(ref_ptr<T>& rp1, ref_ptr<T>& rp2) {
	rp1.swap(rp2);
}

template<class T, class Y> 
inline ref_ptr<T> static_pointer_cast(const ref_ptr<Y>& rp) {
	return std::static_pointer_cast<T>(rp.get_shared());
}

template<class T, class Y> 
inline ref_ptr<T> dynamic_pointer_cast(const ref_ptr<Y>& rp) {
	return std::dynamic_pointer_cast<T>(rp.get_shared());
}

template<class T, class Y> 
inline ref_ptr<T> const_pointer_cast(const ref_ptr<Y>& rp) {
	return std::const_pointer_cast<T>(rp.get_shared());
}

/** @}*/
} // namespace crpropa

#endif // CRPROPA_REFERENCED_H
