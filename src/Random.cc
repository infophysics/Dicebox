/** @file Random.cc
 *  @brief 
 *  @author Nicholas Carrara [ncarrara.physics@gmail.com]
 *  @date 
 *  @details
 */
#include "../include/Random.hh"

namespace dicebox
{
    Random Random::sInstance;
    std::mutex Random::sMutex;

    Random::Random()
    {}
    Random::~Random()
    {}

    Random& Random::Get()
    {
        // if (sInstance == nullptr)
        // {
        //     std::unique_lock<std::mutex> lock(sMutex);
        //     if (sInstance == nullptr)
        //     {
        //         sInstance = new (std::nothrow) Random;
        //     }
        // }
        return sInstance;
    }

    void Random::Delete()
    {
        // std::unique_lock<std::mutex> lock(sMutex);
        // if (sInstance)
        // {
        //     delete sInstance;
        //     sInstance = nullptr;
        // }
    }

    int Random::_uniform_int()
    {
    }
    double Random::_uniform_real()
    {
    }
    double Random::_normal()
    {
        return 0.5;
    }
}