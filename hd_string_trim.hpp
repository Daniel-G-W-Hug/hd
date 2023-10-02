#ifndef HD_STRING_TRIM_HPP
#define HD_STRING_TRIM_HPP

#include <algorithm> // std::transform
#include <string>
namespace hd {

// trim from left
inline std::string& ltrim(std::string& s, const char* t = " \t\n\r\f\v")
{
    s.erase(0, s.find_first_not_of(t));
    return s;
}

// trim from right
inline std::string& rtrim(std::string& s, const char* t = " \t\n\r\f\v")
{
    s.erase(s.find_last_not_of(t) + 1);
    return s;
}

// trim from left & right
inline std::string& trim(std::string& s, const char* t = " \t\n\r\f\v")
{
    return ltrim(rtrim(s, t), t);
}

// copying versions

inline std::string ltrim_cp(std::string s, const char* t = " \t\n\r\f\v")
{
    return ltrim(s, t);
}

inline std::string rtrim_cp(std::string s, const char* t = " \t\n\r\f\v")
{
    return rtrim(s, t);
}

inline std::string trim_cp(std::string s, const char* t = " \t\n\r\f\v")
{
    return trim(s, t);
}

std::string wstring_to_string(const std::wstring& ws)
{
    // does only work for subset of wchar_t
    std::string outstr;
    std::transform(ws.begin(), ws.end(),
                   std::back_inserter(outstr),
                   [](wchar_t c) { return (char)c; });
    return outstr;
}

std::wstring string_to_wstring(const std::string& s)
{
    std::wstring outwstr;
    std::transform(s.begin(), s.end(),
                   std::back_inserter(outwstr),
                   [](char c) { return (wchar_t)c; });
    return outwstr;
}

} // namespace hd

#endif // HD_STRING_TRIM_HPP