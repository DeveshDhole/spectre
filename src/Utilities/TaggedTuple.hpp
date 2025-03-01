// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstddef>
#include <functional>
#include <initializer_list>
#include <ostream>
#include <stack>
#include <string>
#include <utility>

#include "Utilities/Overloader.hpp"
#include "Utilities/PrettyType.hpp"
#include "Utilities/PrintHelpers.hpp"
#include "Utilities/TMPL.hpp"

/// \cond
namespace PUP {
class er;
}  // namespace PUP
/// \endcond

/*!
 * \brief Contains utilities for working with tuples
 */
namespace tuples {

#if __cplusplus >= 201402L
#define TUPLES_LIB_CONSTEXPR_CXX_14 constexpr
#else
#define TUPLES_LIB_CONSTEXPR_CXX_14
#endif

namespace tuples_detail {

template <class T>
inline constexpr T&& forward(typename std::remove_reference<T>::type& t) {
  return static_cast<T&&>(t);
}

template <class T>
inline constexpr T&& forward(typename std::remove_reference<T>::type&& t) {
  static_assert(!std::is_lvalue_reference<T>::value,
                "cannot forward an rvalue as an lvalue");
  return static_cast<T&&>(t);
}

template <class T, T...>
struct value_list {};

template <class...>
struct typelist {};

template <bool... Bs>
using all = typename std::is_same<
    value_list<bool, Bs...>,
    value_list<bool, (static_cast<void>(Bs), true)...>>::type;

struct no_such_type {
  no_such_type() = delete;
  no_such_type(no_such_type const& /*unused*/) = delete;
  no_such_type(no_such_type&& /*unused*/) = delete;
  ~no_such_type() = delete;
  no_such_type& operator=(no_such_type const& /*unused*/) = delete;
  no_such_type operator=(no_such_type&& /*unused*/) = delete;
};

namespace detail {
using std::swap;

template <class T, class S,
          bool = not std::is_void<T>::value and not std::is_void<S>::value>
struct is_swappable_with {
  template <class L, class R>
  static auto test_swap(int)
      -> decltype(swap(std::declval<L&>(), std::declval<R&>()));
  template <class L, class R>
  static tuples::tuples_detail::no_such_type test_swap(...);

  static const bool value =
      not std::is_same<decltype(test_swap<T, S>(0)),
                       tuples::tuples_detail::no_such_type>::value and
      not std::is_same<decltype(test_swap<S, T>(0)),
                       tuples::tuples_detail::no_such_type>::value;
};

template <class T, class S>
struct is_swappable_with<T, S, false> : std::false_type {};
}  // namespace detail

template <class T, class S>
using is_swappable_with = detail::is_swappable_with<T, S>;

template <typename... Ts>
constexpr char expand_pack(Ts&&... /*unused*/) {
  return '0';
}
}  // namespace tuples_detail

namespace tuples_detail {
template <class Tag, bool Ebo = std::is_empty<typename Tag::type>::value &&
                                !__is_final(typename Tag::type)>
class TaggedTupleLeaf;

template <class T, bool B>
void swap(TaggedTupleLeaf<T, B>& lhs, TaggedTupleLeaf<T, B>& rhs) {
  using std::swap;
  swap(lhs.get_data(), rhs.get_data());
}

template <class Tag>
class TaggedTupleLeaf<Tag, false> {
  using value_type = typename Tag::type;
  value_type value_;

  template <class T>
  static constexpr bool can_bind_reference() {
    using rem_ref_value_type = typename std::remove_reference<value_type>::type;
    using rem_ref_T = typename std::remove_reference<T>::type;
    using is_lvalue_type = std::integral_constant<
        bool,
        std::is_lvalue_reference<T>::value or
            std::is_same<std::reference_wrapper<rem_ref_value_type>,
                         rem_ref_T>::value or
            std::is_same<std::reference_wrapper<typename std::remove_const<
                             rem_ref_value_type>::type>,
                         rem_ref_T>::value>;
    return not std::is_reference<value_type>::value or
           (std::is_lvalue_reference<value_type>::value and
            is_lvalue_type::value) or
           (std::is_rvalue_reference<value_type>::value and
            not std::is_lvalue_reference<T>::value);
  }

 public:
  // Tested in constexpr context in Unit.TaggedTuple.Ebo
  constexpr TaggedTupleLeaf() : value_() {
    static_assert(
        !std::is_reference<value_type>::value,
        "Cannot default construct a reference element in a TaggedTuple");
  }

  // clang-tidy: forwarding references can be bad
  template <
      class T,
      typename std::enable_if<
          !std::is_same<typename std::decay<T>::type, TaggedTupleLeaf>::value &&
          std::is_constructible<value_type, T>::value>::type* = nullptr>
  constexpr explicit TaggedTupleLeaf(T&& t)
      : value_(tuples_detail::forward<T>(t)) {
    static_assert(can_bind_reference<T>(),
                  "Cannot construct an lvalue reference with an rvalue");
  }

  constexpr TaggedTupleLeaf(TaggedTupleLeaf const& /*rhs*/) = default;
  constexpr TaggedTupleLeaf(TaggedTupleLeaf&& /*rhs*/) = default;
  constexpr TaggedTupleLeaf& operator=(TaggedTupleLeaf const& /*rhs*/) =
      default;
  constexpr TaggedTupleLeaf& operator=(TaggedTupleLeaf&& /*rhs*/) = default;

  ~TaggedTupleLeaf() = default;

  // Note: name get_data instead of get to enable structured binding support.
#if __cplusplus < 201402L
  value_type& get_data() { return value_; }
#else
  constexpr value_type& get_data() { return value_; }
#endif
  constexpr const value_type& get_data() const { return value_; }

  bool swap(TaggedTupleLeaf& t) {
    using std::swap;
    swap(*this, t);
    return false;
  }

  // clang-tidy: runtime-references
  void pup(PUP::er& p) { p | value_; }  // NOLINT
};

template <class Tag>
class TaggedTupleLeaf<Tag, true> : private Tag::type {
  using value_type = typename Tag::type;

 public:
  constexpr TaggedTupleLeaf() : value_type{} {}

  template <
      class T,
      typename std::enable_if<
          !std::is_same<typename std::decay<T>::type, TaggedTupleLeaf>::value &&
          std::is_constructible<value_type, T&&>::value>::type* = nullptr>
  constexpr explicit TaggedTupleLeaf(T&& t)
      : value_type(tuples_detail::forward<T>(t)) {}

  constexpr TaggedTupleLeaf(TaggedTupleLeaf const& /*rhs*/) = default;
  constexpr TaggedTupleLeaf(TaggedTupleLeaf&& /*rhs*/) = default;
  constexpr TaggedTupleLeaf& operator=(TaggedTupleLeaf const& /*rhs*/) =
      default;
  constexpr TaggedTupleLeaf& operator=(TaggedTupleLeaf&& /*rhs*/) = default;

  ~TaggedTupleLeaf() = default;

  // Note: name get_data instead of get to enable structured binding support.
#if __cplusplus < 201402L
  value_type& get_data() { return static_cast<value_type&>(*this); }
#else
  constexpr value_type& get_data() { return static_cast<value_type&>(*this); }
#endif

  constexpr const value_type& get_data() const {
    return static_cast<const value_type&>(*this);
  }

  bool swap(TaggedTupleLeaf& t) {
    using std::swap;
    swap(*this, t);
    return false;
  }

  // NOLINTNEXTLINE(google-runtime-references)
  void pup(PUP::er& p) { p | static_cast<typename Tag::type&>(*this); }
};

struct disable_constructors {
  static constexpr bool enable_default() { return false; }
  static constexpr bool enable_explicit() { return false; }
  static constexpr bool enable_implicit() { return false; }
};
}  // namespace tuples_detail

/*!
 * \ingroup UtilitiesGroup
 * \brief An associative container that is indexed by structs
 *
 * A data structure that is indexed by Tags. A Tag is a struct that contains
 * a type alias named `type`, which is the type of the object stored with
 * index Tag.
 *
 * \tparam Tags the tags of the objects to be placed in the tuple
 */
template <class... Tags>
class TaggedTuple;

template <class Tag, class... Tags>
constexpr const typename Tag::type& get(const TaggedTuple<Tags...>& t);
template <class Tag, class... Tags>
constexpr typename Tag::type& get(TaggedTuple<Tags...>& t);
template <class Tag, class... Tags>
constexpr const typename Tag::type&& get(const TaggedTuple<Tags...>&& t);
template <class Tag, class... Tags>
constexpr typename Tag::type&& get(TaggedTuple<Tags...>&& t);

/*!
 * \brief Returns the type of the Tag
 */
template <class Tag>
using tag_type = typename Tag::type;

// clang-tidy: class does not define copy or move assignment (it does)
template <class... Tags>
class TaggedTuple : private tuples_detail::TaggedTupleLeaf<Tags>... {  // NOLINT
  template <class... Args>
  struct pack_is_TaggedTuple : std::false_type {};
  template <class... Args>
  struct pack_is_TaggedTuple<TaggedTuple<Args...>> : std::true_type {};

  template <bool EnableConstructor, class Dummy = void>
  struct args_constructor : tuples_detail::disable_constructors {};

  template <class Dummy>
  struct args_constructor<true, Dummy> {
    static constexpr bool enable_default() {
      return tuples_detail::all<
          std::is_default_constructible<tag_type<Tags>>::value...>::value;
    }

    template <class... Ts>
    static constexpr bool enable_explicit() {
      return tuples_detail::all<std::is_constructible<
                 tuples_detail::TaggedTupleLeaf<Tags>, Ts>::value...>::value and
             not tuples_detail::all<
                 std::is_convertible<Ts, tag_type<Tags>>::value...>::value;
    }
    template <class... Ts>
    static constexpr bool enable_implicit() {
      return tuples_detail::all<std::is_constructible<
                 tuples_detail::TaggedTupleLeaf<Tags>, Ts>::value...>::value and
             tuples_detail::all<
                 std::is_convertible<Ts, tag_type<Tags>>::value...>::value;
    }
  };

  // C++17 Draft 23.5.3.2 Assignment - helper aliases
  using is_copy_assignable =
      tuples_detail::all<std::is_copy_assignable<tag_type<Tags>>::value...>;
  using is_nothrow_copy_assignable = tuples_detail::all<
      std::is_nothrow_copy_assignable<tag_type<Tags>>::value...>;
  using is_move_assignable =
      tuples_detail::all<std::is_move_assignable<tag_type<Tags>>::value...>;
  using is_nothrow_move_assignable = tuples_detail::all<
      std::is_nothrow_move_assignable<tag_type<Tags>>::value...>;

  // clang-tidy: redundant declaration
  template <class Tag, class... LTags>
  friend constexpr const typename Tag::type& get(  // NOLINT
      const TaggedTuple<LTags...>& t);
  template <class Tag, class... LTags>
  friend constexpr typename Tag::type& get(  // NOLINT
      TaggedTuple<LTags...>& t);
  template <class Tag, class... LTags>
  friend constexpr const typename Tag::type&& get(  // NOLINT
      const TaggedTuple<LTags...>&& t);
  template <class Tag, class... LTags>
  friend constexpr typename Tag::type&& get(  // NOLINT
      TaggedTuple<LTags...>&& t);

 public:
  using tags_list = tmpl::list<Tags...>;

  static constexpr size_t size() { return sizeof...(Tags); }

  // clang-tidy: runtime-references
  // NOLINTNEXTLINE(google-runtime-references)
  void pup(PUP::er& p) {
    static_cast<void>(std::initializer_list<char>{
        (tuples_detail::TaggedTupleLeaf<Tags>::pup(p), '0')...});
  }

  // C++17 Draft 23.5.3.1 Construction
  template <bool Dummy = true, typename std::enable_if<args_constructor<
                                   Dummy>::enable_default()>::type* = nullptr>
  constexpr TaggedTuple() {}

  TaggedTuple(TaggedTuple const& /*rhs*/) = default;
  TaggedTuple(TaggedTuple&& /*rhs*/) = default;

  /*!
   * \brief Construct a TaggedTuple with Args
   * \requires `std::is_convertible_v<Us, typename Tags::type>...` is `true`
   *
   * \example
   * \snippet Test_TaggedTuple.cpp construction_example
   */
  template <class... Us,
            typename std::enable_if<
                args_constructor<not pack_is_TaggedTuple<Us...>::value and
                                 sizeof...(Us) == sizeof...(Tags)>::
                    template enable_explicit<Us...>()>::type* = nullptr>
  constexpr explicit TaggedTuple(Us&&... us)
      : tuples_detail::TaggedTupleLeaf<Tags>(
            tuples_detail::forward<Us>(us))... {}

  /*!
   * \brief Construct a TaggedTuple with Args
   * \requires `std::is_convertible_v<Us, typename Tags::type>...` is `true`
   *
   * \example
   * \snippet Test_TaggedTuple.cpp construction_example
   */
  template <class... Us,
            typename std::enable_if<
                args_constructor<not pack_is_TaggedTuple<Us...>::value and
                                 sizeof...(Us) == sizeof...(Tags)>::
                    template enable_implicit<Us...>()>::type* = nullptr>
  // clang-tidy: mark explicit
  constexpr TaggedTuple(Us&&... us)
      : tuples_detail::TaggedTupleLeaf<Tags>(
            tuples_detail::forward<Us>(us))... {}

  template <
      class... UTags,
      typename std::enable_if<
          sizeof...(Tags) == sizeof...(UTags) and
          tuples_detail::all<std::is_constructible<
              tag_type<Tags>, const tag_type<UTags>&>::value...>::value and
          not tuples_detail::all<std::is_same<Tags, UTags>::value...>::value>::
          type* = nullptr>
  constexpr explicit TaggedTuple(TaggedTuple<UTags...> const& t)
      : tuples_detail::TaggedTupleLeaf<Tags>(get<UTags>(t))... {}

  template <class... UTags,
            typename std::enable_if<
                sizeof...(Tags) == sizeof...(UTags) and
                tuples_detail::all<std::is_constructible<
                    tag_type<Tags>, tag_type<UTags>&&>::value...>::value and
                not tuples_detail::all<std::is_same<Tags, UTags>::value...>::
                    value>::type* = nullptr>
  constexpr explicit TaggedTuple(TaggedTuple<UTags...>&& t)
      : tuples_detail::TaggedTupleLeaf<Tags>(std::move(get<UTags>(t)))... {}

  ~TaggedTuple() = default;

  // C++17 Draft 23.5.3.2 Assignment
  TaggedTuple& operator=(
      tmpl::conditional_t<is_copy_assignable::value, TaggedTuple,
                          tuples_detail::no_such_type> const& t) {
    static_cast<void>(
        tuples_detail::expand_pack((get<Tags>(*this) = get<Tags>(t))...));
    return *this;
  }

  TaggedTuple& operator=(
      tmpl::conditional_t<is_move_assignable::value, TaggedTuple,
                          tuples_detail::no_such_type>&& t) {
    static_cast<void>(tuples_detail::expand_pack(
        (get<Tags>(*this) =
             tuples_detail::forward<tag_type<Tags>>(get<Tags>(t)))...));
    return *this;
  }

  template <class... UTags,
            typename std::enable_if<
                sizeof...(Tags) == sizeof...(UTags) and
                tuples_detail::all<std::is_assignable<
                    tag_type<Tags>&,
                    tag_type<UTags> const&>::value...>::value>::type* = nullptr>
  TaggedTuple& operator=(TaggedTuple<UTags...> const& t) {
    static_cast<void>(
        tuples_detail::expand_pack((get<Tags>(*this) = get<UTags>(t))...));
    return *this;
  }

  template <
      class... UTags,
      typename std::enable_if<
          sizeof...(Tags) == sizeof...(UTags) and
          tuples_detail::all<std::is_assignable<
              tag_type<Tags>&, tag_type<UTags>&&>::value...>::value>::type* =
          nullptr>
  TaggedTuple& operator=(TaggedTuple<UTags...>&& t) {
    static_cast<void>(tuples_detail::expand_pack(
        (get<Tags>(*this) =
             tuples_detail::forward<tag_type<UTags>>(get<UTags>(t)))...));
    return *this;
  }

  // C++17 Draft 23.5.3.3 swap
  void swap(TaggedTuple& t) {
    tuples_detail::expand_pack(tuples_detail::TaggedTupleLeaf<Tags>::swap(
        static_cast<tuples_detail::TaggedTupleLeaf<Tags>&>(t))...);
  }
};

template <>
class TaggedTuple<> {
 public:
  using tags_list = tmpl::list<>;
  static constexpr size_t size() { return 0; }
  TaggedTuple() = default;
  void swap(TaggedTuple& /*unused*/) {}
  // clang-tidy: runtime-references
  void pup(PUP::er& /*p*/) {}  // NOLINT
};

// C++17 Draft 23.5.3.6 Tuple helper classes
template <class T>
struct tuple_size;

template <class... Tags>
struct tuple_size<TaggedTuple<Tags...>>
    : std::integral_constant<size_t, sizeof...(Tags)> {};
template <class... Tags>
struct tuple_size<const TaggedTuple<Tags...>>
    : tuple_size<TaggedTuple<Tags...>> {};
template <class... Tags>
struct tuple_size<volatile TaggedTuple<Tags...>>
    : tuple_size<TaggedTuple<Tags...>> {};
template <class... Tags>
struct tuple_size<const volatile TaggedTuple<Tags...>>
    : tuple_size<TaggedTuple<Tags...>> {};

// C++17 Draft 23.5.3.7 Element access
/// @{
/*!
 * \ingroup UtilitiesGroup
 * \brief Retrieve the element of `Tag` in the TaggedTuple
 */
template <class Tag, class... Tags>
inline constexpr const typename Tag::type& get(const TaggedTuple<Tags...>& t) {
  static_assert(std::is_base_of<tuples_detail::TaggedTupleLeaf<Tag>,
                                TaggedTuple<Tags...>>::value,
                "Could not retrieve Tag from TaggedTuple. See the first "
                "template parameter of the instantiation for what Tag is being "
                "retrieved and the remaining template parameters for what Tags "
                "are available.");
  return static_cast<const tuples_detail::TaggedTupleLeaf<Tag>&>(t).get_data();
}
template <class Tag, class... Tags>
inline constexpr typename Tag::type& get(TaggedTuple<Tags...>& t) {
  static_assert(std::is_base_of<tuples_detail::TaggedTupleLeaf<Tag>,
                                TaggedTuple<Tags...>>::value,
                "Could not retrieve Tag from TaggedTuple. See the first "
                "template parameter of the instantiation for what Tag is being "
                "retrieved and the remaining template parameters for what Tags "
                "are available.");
  return static_cast<tuples_detail::TaggedTupleLeaf<Tag>&>(t).get_data();
}
template <class Tag, class... Tags>
inline constexpr const typename Tag::type&& get(
    const TaggedTuple<Tags...>&& t) {
  static_assert(std::is_base_of<tuples_detail::TaggedTupleLeaf<Tag>,
                                TaggedTuple<Tags...>>::value,
                "Could not retrieve Tag from TaggedTuple. See the first "
                "template parameter of the instantiation for what Tag is being "
                "retrieved and the remaining template parameters for what Tags "
                "are available.");
  return static_cast<const typename Tag::type&&>(
      static_cast<const tuples_detail::TaggedTupleLeaf<Tag>&&>(t).get_data());
}
template <class Tag, class... Tags>
inline constexpr typename Tag::type&& get(TaggedTuple<Tags...>&& t) {
  static_assert(std::is_base_of<tuples_detail::TaggedTupleLeaf<Tag>,
                                TaggedTuple<Tags...>>::value,
                "Could not retrieve Tag from TaggedTuple. See the first "
                "template parameter of the instantiation for what Tag is being "
                "retrieved and the remaining template parameters for what Tags "
                "are available.");
  return static_cast<typename Tag::type&&>(
      static_cast<tuples_detail::TaggedTupleLeaf<Tag>&&>(t).get_data());
}
/// @}

template <size_t I, class... Tags>
inline constexpr typename tmpl::at_c<tmpl::list<Tags...>, I>::type&& get(
    TaggedTuple<Tags...>&& t) {
  return get<tmpl::at_c<tmpl::list<Tags...>, I>>(t);
}

template <size_t I, class... Tags>
inline constexpr const typename tmpl::at_c<tmpl::list<Tags...>, I>::type& get(
    const TaggedTuple<Tags...>& t) {
  return get<tmpl::at_c<tmpl::list<Tags...>, I>>(t);
}

template <size_t I, class... Tags>
inline constexpr typename tmpl::at_c<tmpl::list<Tags...>, I>::type& get(
    TaggedTuple<Tags...>& t) {
  return get<tmpl::at_c<tmpl::list<Tags...>, I>>(t);
}

// C++17 Draft 23.5.3.8 Relational operators
namespace tuples_detail {
struct equal {
  template <class T, class U>
  static TUPLES_LIB_CONSTEXPR_CXX_14 void apply(T const& lhs, U const& rhs,
                                                bool* result) {
    *result = *result and lhs == rhs;
  }
};

template <class... LTags, class... RTags>
TUPLES_LIB_CONSTEXPR_CXX_14 bool tuple_equal_impl(
    TaggedTuple<LTags...> const& lhs, TaggedTuple<RTags...> const& rhs) {
  bool equal = true;
  // This short circuits in the sense that the operator== is only evaluated if
  // the result thus far is true
  static_cast<void>(std::initializer_list<char>{
      (equal::apply(get<LTags>(lhs), get<RTags>(rhs), &equal), '0')...});
  return equal;
}
}  // namespace tuples_detail

template <class... LTags, class... RTags,
          typename std::enable_if<sizeof...(LTags) == sizeof...(RTags)>::type* =
              nullptr>
TUPLES_LIB_CONSTEXPR_CXX_14 bool operator==(TaggedTuple<LTags...> const& lhs,
                                            TaggedTuple<RTags...> const& rhs) {
  return tuples_detail::tuple_equal_impl(lhs, rhs);
}

template <class... LTags, class... RTags,
          typename std::enable_if<sizeof...(LTags) == sizeof...(RTags)>::type* =
              nullptr>
TUPLES_LIB_CONSTEXPR_CXX_14 bool operator!=(TaggedTuple<LTags...> const& lhs,
                                            TaggedTuple<RTags...> const& rhs) {
  return not(lhs == rhs);
}

namespace tuples_detail {
struct less {
  template <class T, class U>
  static TUPLES_LIB_CONSTEXPR_CXX_14 void apply(T const& lhs, U const& rhs,
                                                bool* last_rhs_less_lhs,
                                                bool* result) {
    if (*result or *last_rhs_less_lhs) {
      return;
    }
    *result = lhs < rhs;
    if (*result) {
      return;
    }
    *last_rhs_less_lhs = rhs < lhs;
  }
};

template <class... LTags, class... RTags>
TUPLES_LIB_CONSTEXPR_CXX_14 bool tuple_less_impl(
    TaggedTuple<LTags...> const& lhs, TaggedTuple<RTags...> const& rhs) {
  bool result = false;
  bool last_rhs_less_lhs = false;
  static_cast<void>(
      std::initializer_list<char>{(less::apply(get<LTags>(lhs), get<RTags>(rhs),
                                               &last_rhs_less_lhs, &result),
                                   '0')...});
  return result;
}
}  // namespace tuples_detail

template <class... LTags, class... RTags,
          typename std::enable_if<sizeof...(LTags) == sizeof...(RTags)>::type* =
              nullptr>
TUPLES_LIB_CONSTEXPR_CXX_14 bool operator<(TaggedTuple<LTags...> const& lhs,
                                           TaggedTuple<RTags...> const& rhs) {
  return tuples_detail::tuple_less_impl(lhs, rhs);
}

template <class... LTags, class... RTags,
          typename std::enable_if<sizeof...(LTags) == sizeof...(RTags)>::type* =
              nullptr>
TUPLES_LIB_CONSTEXPR_CXX_14 bool operator>(TaggedTuple<LTags...> const& lhs,
                                           TaggedTuple<RTags...> const& rhs) {
  return rhs < lhs;
}

template <class... LTags, class... RTags,
          typename std::enable_if<sizeof...(LTags) == sizeof...(RTags)>::type* =
              nullptr>
TUPLES_LIB_CONSTEXPR_CXX_14 bool operator<=(TaggedTuple<LTags...> const& lhs,
                                            TaggedTuple<RTags...> const& rhs) {
  return not(rhs < lhs);
}

template <class... LTags, class... RTags,
          typename std::enable_if<sizeof...(LTags) == sizeof...(RTags)>::type* =
              nullptr>
TUPLES_LIB_CONSTEXPR_CXX_14 bool operator>=(TaggedTuple<LTags...> const& lhs,
                                            TaggedTuple<RTags...> const& rhs) {
  return not(lhs < rhs);
}

// C++17 Draft 23.5.3.3 swap
template <
    class... Tags,
    typename std::enable_if<tuples_detail::all<tuples_detail::is_swappable_with<
        tuples_detail::TaggedTupleLeaf<Tags>,
        tuples_detail::TaggedTupleLeaf<Tags>>::value...>::value>::type* =
        nullptr>
void swap(TaggedTuple<Tags...>& lhs, TaggedTuple<Tags...>& rhs) {
  lhs.swap(rhs);
}

namespace TaggedTuple_detail {
template <typename T>
struct tagged_tuple_typelist_impl;

template <template <typename...> class List, typename... Tags>
struct tagged_tuple_typelist_impl<List<Tags...>> {
  using type = TaggedTuple<Tags...>;
};
}  // namespace TaggedTuple_detail

/// \ingroup UtilitiesGroup
template <typename T>
using tagged_tuple_from_typelist =
    typename TaggedTuple_detail::tagged_tuple_typelist_impl<T>::type;

namespace TaggedTuple_detail {
template <typename... InputTags, typename... OutputTags>
TaggedTuple<OutputTags...> reorder_impl(TaggedTuple<InputTags...>&& input,
                                        tmpl::list<OutputTags...> /*meta*/) {
  static_assert(
      std::is_same_v<tmpl::list_difference<tmpl::list<OutputTags...>,
                                           tmpl::list<InputTags...>>,
                     tmpl::list<>> and
          std::is_same_v<tmpl::list_difference<tmpl::list<InputTags...>,
                                               tmpl::list<OutputTags...>>,
                         tmpl::list<>>,
      "The input and output TaggedTuples must be the same except "
      "for ordering.");
  return TaggedTuple<OutputTags...>(std::move(get<OutputTags>(input))...);
}
}  // namespace TaggedTuple_detail

/// Given an input TaggedTuple, produce an output TaggedTuple
/// with the tags in a different order.  All tags must be the same
/// except for ordering.
/// \example
/// \snippet Test_TaggedTuple.cpp reorder_example
template <typename ReturnedTaggedTuple, typename... Tags>
ReturnedTaggedTuple reorder(TaggedTuple<Tags...> input) {
  return TaggedTuple_detail::reorder_impl(
      std::move(input), typename ReturnedTaggedTuple::tags_list{});
}

/// Stream operator for TaggedTuple
using ::operator<<;
template <class... Tags>
std::ostream& operator<<(std::ostream& os, const TaggedTuple<Tags...>& t) {
  os << "TaggedTuple:\n";
  const auto print_item = [&os, &t](auto tag_v) {
    using tag = tmpl::type_from<decltype(tag_v)>;
    using type = typename tag::type;
    os << "----------\n";
    os << "Name:  " << pretty_type::get_name<tag>() << "\n";
    os << "Type:  " << pretty_type::get_name<type>() << "\n";
    os << "Value: ";
    print_value(os, get<tag>(t));
    os << "\n";
  };
  tmpl::for_each<tmpl::list<Tags...>>(print_item);
  return os;
}

namespace TaggedTuple_detail {

template <typename F, typename... Tags, typename... ApplyTags>
constexpr decltype(auto) apply_impl(F&& f, const TaggedTuple<Tags...>& t,
                                    tmpl::list<ApplyTags...> /* meta */) {
  return std::forward<F>(f)(get<ApplyTags>(t)...);
}

}  // namespace TaggedTuple_detail

/// @{
/*!
 * \ingroup UtilitiesGroup
 * \brief Invoke `f` with the `ApplyTags` taken from `t` expanded in a parameter
 * pack
 *
 * `ApplyTags` defaults to the full list of tags in `t`.
 *
 * Here is an example how to use the function:
 *
 * \snippet Test_TaggedTuple.cpp expand_tuple_example
 *
 * This is the function being called in the above example:
 *
 * \snippet Test_TaggedTuple.cpp expand_tuple_example_function
 *
 * \see std::apply
 */
template <typename ApplyTags, typename F, typename... Tags>
constexpr decltype(auto) apply(F&& f, const TaggedTuple<Tags...>& t) {
  return TaggedTuple_detail::apply_impl(std::forward<F>(f), t, ApplyTags{});
}

template <typename F, typename... Tags>
constexpr decltype(auto) apply(F&& f, const TaggedTuple<Tags...>& t) {
  return TaggedTuple_detail::apply_impl(
      std::forward<F>(f), t, typename TaggedTuple<Tags...>::tags_list{});
}
/// @}

}  // namespace tuples

namespace std {
template <typename... Tags>
struct tuple_size<tuples::TaggedTuple<Tags...>>
    : std::integral_constant<int, sizeof...(Tags)> {};
template <size_t I, typename... Tags>
struct tuple_element<I, tuples::TaggedTuple<Tags...>> {
  using type = typename tmpl::at_c<tmpl::list<Tags...>, I>::type;
};
}  // namespace std
