#ifndef BLOCKINGQUEUE_HPP_
#define BLOCKINGQUEUE_HPP_

#include <algorithm>
#include <mutex>
#include <condition_variable>
#include <memory>

// A simple parallel task: holds a bound callable as std::function<void()>
struct ParallelTask {
    std::function<void()> fn;
    bool poison = false;

    ParallelTask() = default;

    template <class F, class... Args>
    explicit ParallelTask(F&& f, Args&&... args) {
	// Bind f and args without invoking; execute later in consumer thread
	auto bound = std::tuple<std::decay_t<F>, std::decay_t<Args>...>(
	    std::forward<F>(f), std::forward<Args>(args)...
	    );
	fn = [t = std::move(bound)]() mutable {
	    std::apply([](auto&& f0, auto&&... a0) {
		std::invoke(std::forward<decltype(f0)>(f0),
			    std::forward<decltype(a0)>(a0)...);
	    }, t);
	};
    }

    // Factory for poison pill
    static ParallelTask poisonpill() {
	ParallelTask t;
	t.poison = true;
	return t;
    }

    bool is_poison() const noexcept { return this->poison; }

    void operator()() { if (fn) fn(); }
    bool valid() const noexcept { return static_cast<bool>(fn); }
};


/**
 * A fixed capacity, thread stable queue used to store
 * tasks in a producer/consumer setup
 */
template <typename T>
class BlockingQueue
{

private:
    std::size_t capacity;

    std::mutex mtx;

    std::condition_variable c_head;
    std::condition_variable c_tail;

    std::size_t tail;
    std::size_t head;

    std::unique_ptr<T[]> queue;

public:
    // BlockingQueue() {
    //     capacity = 100;
    //     head = 0;
    //     tail = 0;
    //     queue = std::make_unique<T[]>(capacity);
    // }
  
    BlockingQueue(std::size_t _capacity) {
        capacity = _capacity;
        head = 0;
        tail = 0;
        queue = std::unique_ptr<T[]>(new T[capacity]);
    }

    BlockingQueue(const BlockingQueue& other) : BlockingQueue(other.capacity) {
	head = other.head;
	tail = other.tail;
	std::copy_n(other.queue.get(), other.capacity, queue.get());
    }

    BlockingQueue& operator=(const BlockingQueue& other) {
	capacity = other.capacity;
	head = other.head;
	tail = other.tail;
	queue = std::unique_ptr<T[]>(new T[capacity]);
        std::copy_n(other.queue.get(), other.capacity, queue.get());
	return *this;
    }

    BlockingQueue(BlockingQueue&& other) : BlockingQueue(0) {
	std::swap(capacity, other.capacity);
	std::swap(head, other.head);
	std::swap(tail, other.tail);
	std::swap(queue, other.queue);
	// other.queue = NULL;
    }

    BlockingQueue& operator=(BlockingQueue&& other) {
        std::swap(capacity, other.capacity);
	std::swap(head, other.head);
	std::swap(tail, other.tail);
	std::swap(queue, other.queue);
	// other.queue = NULL;
	return *this;
    }

    ~BlockingQueue() = default;
    // {
    // 	if (queue != nullptr)
    // 	    delete[] queue;
    // }

    void push(T const &value)
    {
        std::unique_lock<std::mutex> lock(mtx);
        c_tail.wait(lock, [&] { return (tail - head) != capacity; });

        queue.get()[tail % capacity] = value;
        tail++;
        c_head.notify_one();
    }

    T pop()
    {
        std::unique_lock<std::mutex> lock(mtx);
        c_head.wait(lock, [&] { return head != tail; });
	
        T ret(std::move(queue.get()[head % capacity]));
        head++;
        c_tail.notify_one();
	
        return ret;
    }

    std::size_t getCapacity() { return capacity; }

    // TODO - lock?
    std::size_t getSize() {
	std::unique_lock<std::mutex> lock(mtx);
	return tail - head;
    }
    bool isEmpty() {
	std::unique_lock<std::mutex> lock(mtx);
	return head == tail;
    }
};

#endif /* BLOCKINGQUEUE_HPP_ */
