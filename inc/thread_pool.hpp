#ifndef THREAD_POOL_H
#define THREAD_POOL_H

#include <thread>
#include <mutex>
#include <condition_variable>
#include <functional>
#include <queue>

class thread_pool {
public:

    thread_pool(int threads);
    ~thread_pool();

    void do_job(std::function<void(void)> func);

    void wait();

protected:

    void enter_thread();

    bool in_destructor_;
    std::vector<std::thread> threads_;
    std::queue<std::function<void(void)>> jobs_;
    std::mutex job_mutex_;
    std::condition_variable job_sync_;
    int waiting_on_;
    std::mutex wait_mutex_;
    std::condition_variable wait_sync_;

};

#endif