#include "hongyulibcpp.h"

template<class Iter, class new_record>
Iter binary_find(Iter begin, Iter end, Iter val)
{
    // Finds the lower bound in at most log(last - first) + 1 comparisons
    Iter i = std::lower_bound(begin, end, val);

    if (i != end && !(val < *i))
        return i; // found
    else
        return end; // not found
}
