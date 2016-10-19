#ifndef TRISYCL_SYCL_ID_HPP
#define TRISYCL_SYCL_ID_HPP

/** \file The OpenCL SYCL id<>

    Ronan at Keryell point FR

    This file is distributed under the University of Illinois Open Source
    License. See LICENSE.TXT for details.
*/

#include <algorithm>
#include <cstddef>

#include "CL/sycl/detail/small_array.hpp"
#include "CL/sycl/range.hpp"

namespace cl {
namespace sycl {

template <std::size_t dims> class item;

/** \addtogroup parallelism Expressing parallelism through kernels
    @{
*/

/** Define a multi-dimensional index, used for example to locate a work
    item
*/
template <std::size_t dims = 1>
class id : public detail::small_array_123<std::size_t, id<dims>, dims> {

public:

  // Inherit from all the constructors
  using detail::small_array_123<std::size_t,
                                id<dims>,
                                dims>::small_array_123;


  /// Construct an id from the dimensions of a range
  id(const range<dims> &range_size)
    /** Use the fact we have a constructor of a small_array from a another
        kind of small_array
     */
    : detail::small_array_123<std::size_t, id<dims>, dims> { range_size } {}


  /// Construct an id from an item global_id
  id(const item<dims> &rhs)
    : detail::small_array_123<std::size_t, id<dims>, dims>
      { rhs.get() }
  {}

  /// Keep other constructors
  id() = default;

};


/** Implement a make_id to construct an id<> of the right dimension with
    implicit conversion from an initializer list for example.

    Cannot use a template on the number of dimensions because the implicit
    conversion would not be tried. */
inline auto make_id(id<1> i) { return i; }
inline auto make_id(id<2> i) { return i; }
inline auto make_id(id<3> i) { return i; }


/** Construct an id<> from a function call with arguments, like
    make_id(1, 2, 3) */
template<typename... BasicType>
auto make_id(BasicType... Args) {
  // Call constructor directly to allow narrowing
  return id<sizeof...(Args)>(Args...);
}

/// @} End the parallelism Doxygen group

}
}

/*
    # Some Emacs stuff:
    ### Local Variables:
    ### ispell-local-dictionary: "american"
    ### eval: (flyspell-prog-mode)
    ### End:
*/

#endif // TRISYCL_SYCL_ID_HPP
