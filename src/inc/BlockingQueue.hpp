#ifndef BLOCKINGQUEUE_HPP_
#define BLOCKINGQUEUE_HPP_

#include <mutex>
#include <condition_variable>
#include <memory>
#include <fstream>

/**
 * A fixed capacity, thread stable queue used to store
 * tasks in a producer/consumer setup
 */
template <typename T>
class BlockingQueue
{

private:
    int capacity;

    std::mutex mtx;

    std::condition_variable c_head;
    std::condition_variable c_tail;

    unsigned long tail;
    unsigned long head;

    std::unique_ptr<T[]> queue;

public:
    // BlockingQueue() {
    //     capacity = 100;
    //     head = 0;
    //     tail = 0;
    //     queue = std::make_unique<T[]>(capacity);
    // }
  
    BlockingQueue(int _capacity) {
        capacity = _capacity;
        head = 0;
        tail = 0;
        queue = std::make_unique<T[]>(capacity);
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
	queue = std::make_unique<T[]>(capacity);
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
    // 	if (queue != NULL) delete[] queue;
    // }

    void push(T const &value)
    {
        std::unique_lock<std::mutex> lock(mtx);
        c_tail.wait(lock, [=] { return (tail - head) != capacity; });

        queue.get()[tail % capacity] = value;
        tail++;
        c_head.notify_one();
    }

    T pop()
    {
        std::unique_lock<std::mutex> lock(mtx);
        c_head.wait(lock, [=] { return head != tail; });
	
        T ret(std::move(queue.get()[head % capacity]));
        head++;
        c_tail.notify_one();
	
        return ret;
    }

    int getCapacity() { return capacity; }

    // TODO - lock?
    int getSize() {
	std::unique_lock<std::mutex> lock(mtx);
	return tail - head;
    }
    bool isEmpty() {
	std::unique_lock<std::mutex> lock(mtx);
	return head == tail;
    }
};

#endif /* BLOCKINGQUEUE_HPP_ */
