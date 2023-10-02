# hd

- provides header only support functions in namespace hd
- to use, include corresponding header (.hpp) in your code

- testing of functionality implemented with doctest in corresponding "..._test.cpp"
  Not needed otherwise, i.e. can be ignored for usage.
  cmake and corresponding CMakeLists.txt is just used to complile test cases, otherwise "header only" for usage.

Dependencies:
  - doctest (e.g. brew install doctest)
  - date (e.g. brew install howard-hinnant-date)