#ifndef HD_THRDSF_QUEUE_HPP
#define HD_THRDSF_QUEUE_HPP
//
// threadsafe queue (based on A. Williams "Concurrency in Action")
//

#include <condition_variable>
#include <exception>
#include <memory>
#include <mutex>
#include <queue>

namespace hd {

template <typename T> class thrdsf_queue {

  private:

    mutable std::mutex mtx;
    std::queue<T> data_queue;
    std::condition_variable data_cond;

  public:

    thrdsf_queue() {}

    void push(T new_value)
    {

        std::shared_ptr<T> value(make_shared<T>(std::move(new_value)));
        std::lock_guard<std::mutex> lk(mtx);
        data_queue.push(value);
        data_cond.notify_one();
    }

    void wait_and_pop(T& value)
    {
        std::unique_lock<std::mutex> lk(mtx);
        data_cond.wait(lk, [this] { return !data_queue.empty(); });
        value = std::move(*data_queue.front());
        data_queue.pop();
    }

    bool try_pop(T& value)
    {
        std::lock_guard<std::mutex> lk(mtx);
        if (data_queue.empty()) return false;
        value = std::move(*data_queue.front());
        data_queue.pop();
        return true;
    }

    std::shared_ptr<T> wait_and_pop()
    {
        std::unique_lock<std::mutex> lk(mtx);
        data_cond.wait(lk, [this] { return !data_queue.empty(); });
        std::shared_ptr<T> value = data_queue.front();
        data_queue.pop();
        return value;
    }

    std::shared_ptr<T> try_pop()
    {
        std::lock_guard<std::mutex> lk(mtx);
        if (data_queue.empty()) return std::shared_ptr<T>(); //  return nullptr
        std::shared_ptr<T> value = data_queue.front();
        data_queue.pop();
        return value;
    }

    bool empty() const
    {
        std::lock_guard<std::mutex> lk(mtx);
        return data_queue.empty();
    }

}; // class thrdsf_queue

} // namespace hd

#endif // HD_THRDSF_QUEUE_HPP