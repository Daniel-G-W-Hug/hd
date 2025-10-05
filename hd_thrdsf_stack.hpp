#ifndef HD_THRDSF_STACK_HPP
#define HD_THRDSF_STACK_HPP
//
// threadsafe stack (based on A. Williams "Concurrency in Action")
//
// protect all access routines with a mutex to get a threadsafe stack
// (only pop - not separate pop & top - to avoid race condition)
//

#include <exception>
#include <memory>
#include <mutex>
#include <stack>

namespace hd {

struct empty_stack : std::exception {
    const char* what() const throw() { return "HD::thrdsf_stack: Empty stack!"; }
};

template <typename T> class thrdsf_stack {

  private:

    std::stack<T> data;
    mutable std::mutex m;

  public:

    thrdsf_stack() {}

    thrdsf_stack(const thrdsf_stack& other)
    {
        // copy in ctor body to be able to lock other.m
        std::lock_guard<std::mutex> lock(other.m);
        data = other.data;
    }

    thrdsf_stack& operator=(const thrdsf_stack& other) = delete;

    void push(T new_value)
    {
        std::lock_guard<std::mutex> lock(m);
        data.push(new_value);
    }

    std::shared_ptr<T> pop()
    {
        std::lock_guard<std::mutex> lock(m);
        if (data.empty())
            throw empty_stack(); // check for data.empty() before trying to top
        std::shared_ptr<T> const value(
            std::make_shared<T>(data.top())); // allocate return value before modifying
                                              // stack (safe in face of exceptions)
        data.pop();
        return value;
    }

    void pop(T& value)
    {
        std::lock_guard<std::mutex> lock(m);
        if (data.empty())
            throw empty_stack(); // check for data.empty() before trying to top
        value = data.top();      // storage for value provided by user
        data.pop();
    }

    bool empty() const
    {
        std::lock_guard<std::mutex> lock(m);
        return data.empty();
    }

}; // class thrdsf_stack

} // namespace hd

#endif // HD_THRDSF_STACK_HPP
