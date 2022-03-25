#include "thread_pool.hpp"
#include <functional>

thread_pool::thread_pool(int threads) : in_destructor_(false), waiting_on_(0) {
    threads_.reserve(threads);
    for (int i = 0; i < threads; i++) threads_.emplace_back(std::bind(&thread_pool::enter_thread, this));
}

thread_pool::~thread_pool() {
    {
        auto lock = std::unique_lock(job_mutex_);

        in_destructor_ = true;
        job_sync_.notify_all();
    }

    for (auto& thread : threads_) thread.join();
}

void thread_pool::do_job(std::function<void(void)> func) {
    {
        auto lock = std::unique_lock(wait_mutex_);

        waiting_on_ += 1;
    }
    {
        auto lock = std::unique_lock(job_mutex_);

        jobs_.emplace(std::move(func));
        job_sync_.notify_one();
    }
}

void thread_pool::wait() {
    auto lock = std::unique_lock(wait_mutex_);

    wait_sync_.wait(
        lock,
        [&] {
            return waiting_on_ == 0;
        }
    );
}

void thread_pool::enter_thread() {
    std::function<void(void)> job;

    while (true) {
        {
            auto lock = std::unique_lock(job_mutex_);

            job_sync_.wait(
                lock,
                [&] {
                    return !(!in_destructor_ && jobs_.empty());
                }
            );

            if (jobs_.empty()) return;

            // consider what happens if we get to the destructor and have many jobs to do?

            job = std::move(jobs_.front());
            jobs_.pop();
        }

        job();

        {
            auto active_lock = std::unique_lock(wait_mutex_);

            waiting_on_ -= 1;
            wait_sync_.notify_all();
        }
    }
}