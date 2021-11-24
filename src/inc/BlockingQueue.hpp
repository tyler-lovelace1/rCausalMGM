#ifndef BLOCKINGQUEUE_HPP_
#define BLOCKINGQUEUE_HPP_

#include <mutex>
#include <condition_variable>
#include <list>
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

    T *queue;

public:
    BlockingQueue(int _capacity)
    {
        capacity = _capacity;
        head = 0;
        tail = 0;
        queue = new T[capacity];
    }

    ~BlockingQueue()
    {
        delete[] queue;
    }

    void push(T const &value)
    {
        std::unique_lock<std::mutex> lock(mtx);
        c_tail.wait(lock, [=] { return (tail - head) != capacity; });

        queue[tail % capacity] = value;
        tail++;
        c_head.notify_one();
    }

    T pop()
    {
        std::unique_lock<std::mutex> lock(mtx);
        c_head.wait(lock, [=] { return head != tail; });
	
        T ret(std::move(queue[head % capacity]));
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
