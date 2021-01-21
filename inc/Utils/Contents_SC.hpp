#ifndef CONTENTS_SC_HPP
#define CONTENTS_SC_HPP

template <typename B, typename R>
class Contents_SC
{
public:
	std::vector<R> lambda;			// alpha
	std::vector<B> s;				// beta
	bool           is_frozen_bit;

	explicit Contents_SC(int size) : lambda(size), s(size), is_frozen_bit(false) {}
	virtual ~Contents_SC() {}
};


#endif