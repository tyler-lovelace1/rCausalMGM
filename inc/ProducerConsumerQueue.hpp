#ifndef PRODUCERCONSUMERQUEUE_HPP_
#define PRODUCERCONSUMERQUEUE_HPP_

#include <mutex>
#include <condition_variable>
#include <boost/circular_buffer.hpp>

/**
 * A fixed capacity, thread stable queue used to store
 * tasks in a producer/consumer setup
 */ 
template <typename T>
class ProducerConsumerQueue {

private:

    int capacity;

    std::mutex mu;

    std::condition_variable producerGo;
    std::condition_variable consumerGo;

    boost::circular_buffer<T> queue;

public:
    ProducerConsumerQueue(int _capacity) {
        capacity = _capacity;
        queue = boost::circular_buffer<T>(capacity);
    }

    void push(const T& t);
    T pop();

    int getCapacity() {return capacity;}

    // TODO - lock?
    int getSize() {return queue.size();}
    bool isEmpty() {return queue.empty();}

};

#endif /* PRODUCERCONSUMERQUEUE_HPP_ */