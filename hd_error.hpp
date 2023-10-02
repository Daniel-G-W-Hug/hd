#ifndef HD_ERROR_HPP
#define HD_ERROR_HPP

#include <cassert> // assert() for run-time assertion checking
                   // static_assert() for compile-time assertion checking
#include <exception>
#include <iostream> // std::cerr, std::cout
#include <source_location>
#include <string>

namespace hd {

using namespace std::string_literals;

// simple location message including source file location
void file_loc_msg(std::string const& message = ""s,
                  const std::source_location location = std::source_location::current())
{
    std::cout << "file: " << location.file_name() << '\n'
              << "(line " << location.line() << ", column " << location.column() << ") "
              << "in function '" << location.function_name() << "':\n"
              << message << "\n";
}

// error handling as proposed by
// Stroustrup, Bjarne. Tour of C++, 3rd edition, p.48. Pearson Education. Kindle-Version.

enum class Error_action {
    ignore,
    throwing,
    terminating,
    logging
}; // error-handling alternatives

// default (=> user defined & evaluated at compile time)
constexpr Error_action default_Error_action = Error_action::logging;

enum class Error_code // individual errors
{
    range_error,
    length_error
}; // (extend as required)

// names of individual errors
std::string error_code_name[]{
    "range error"s,
    "length error"s}; // (extend as required)

//
// usage: pass condition via lambda to expect
//
// example:
//
// double& Vector::operator[](int i)
// {
//   expect([i,this]{return 0<=i && i<size();}, Error_code::range_error);
//   return elem[i];
// }
//

template <Error_action action = default_Error_action, class C>
constexpr void expect(C cond,
                      Error_code x,
                      const std::source_location location = std::source_location::current())
// take "action" if the expected condition "cond" doesn't hold
// if constexpr tests done at compile time depending on Error_action
// at most one runtime test of cond() is performed for each call of expect
{
    if constexpr (action == Error_action::logging)
        if (!cond)
            std::cerr
                << "file: "
                << location.file_name() << '\n'
                << "(line " << location.line() << ", column " << location.column()
                << ") in function '" << location.function_name() << "':\n"
                << "expect() failure (#" << int(x) << ", " << error_code_name[int(x)] << ")\n";
    if constexpr (action == Error_action::throwing)
        if (!cond)
            throw x;
    if constexpr (action == Error_action::terminating)
        if (!cond)
            std::terminate();
    // or no action (in case of Error_action::ignore)
}

// By setting default_Error_action a user can select an action
// suitable for a particular deployment of the program, such as terminating or logging.
// To support logging, a table of error_code_names needs to be defined.
// The logging information could be improved by using source_location (§16.5).

} // namespace hd

#endif // HD_ERROR_HPP