#ifndef BLOCKINGQUEUE_HPP_
#define BLOCKINGQUEUE_HPP_

#include <mutex>
#include <condition_variable>
#include <list>
// #include <boost/circular_buffer.hpp>

/**
 * A fixed capacity, thread stable queue used to store
 * tasks in a producer/consumer setup
 */ 
template <typename T>
class BlockingQueue {

private:

    int capacity;

    std::mutex m_head;
    std::mutex m_tail;

    // std::mutex mtx;
    // std::mutex f_mtx;

    std::condition_variable c_head;
    std::condition_variable c_tail;

    // std::condition_variable cv;
    // std::condition_variable full;

    unsigned long tail;
    unsigned long head;

    T* queue;

public:
    BlockingQueue(int _capacity) {
        capacity = _capacity;
	head = 0;
	tail = 0;
        queue = new T[capacity];
    }

    ~BlockingQueue() {
        delete[] queue;
    }

    void push(T const& value) {
	std::unique_lock<std::mutex> lock(m_tail);
	c_tail.wait(lock, [=]{ return getSize() != getCapacity(); });
	// {
	//     std::lock_guard<std::mutex> lock(mtx);
	//     Rcpp::Rcout << "push v = " << value << std::endl;
	//     Rcpp::Rcout << "push Q: " << head << " " << tail << std::endl;
	// }
	queue[tail % capacity] = value;
	tail++;
        c_head.notify_one();
    }
    T pop() {
        std::unique_lock<std::mutex> lock(m_head);
        c_head.wait(lock, [=]{ return !isEmpty(); });
	T ret(std::move(queue[head % capacity]));
	// {
	//     std::lock_guard<std::mutex> lock(mtx);
	//     Rcpp::Rcout << "pop v = " << ret << std::endl;
	//     Rcpp::Rcout << "pop Q: " << head << " " << tail << std::endl;
	// }
	head++;
	c_tail.notify_one();
        return ret;
    }

    int getCapacity() {return capacity;}

    // TODO - lock?
    int getSize() {return tail - head;} //queue.size();}
    bool isEmpty() {return head == tail; } //queue.empty();}

};

#endif /* BLOCKINGQUEUE_HPP_ */
