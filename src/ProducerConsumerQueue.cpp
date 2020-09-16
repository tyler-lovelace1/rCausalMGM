#include "ProducerConsumerQueue.hpp"

template <typename T>
void ProducerConsumerQueue<T>::push(const T& t) {
    // std::unique_lock<std::mutex> lk(mu);
    // producerGo.wait(lk, [&]() { return queue.size() < capacity; }); // Wait for space in buffer
    // queue.push_back(t);
    // lk.unlock();
    // consumerGo.notify_one();
}

template <typename T>
T ProducerConsumerQueue<T>::pop() {
    std::unique_lock<std::mutex> lk(mu);
    consumerGo.wait(lk, [&]() { return !queue.empty(); });
    T ret = queue.front();
    queue.pop_front();
    lk.unlock();
    producerGo.notify_one();
    return ret;
}