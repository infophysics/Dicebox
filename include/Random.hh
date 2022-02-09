/** @file Random.hh
 *  @brief 
 *  @author Nicholas Carrara [ncarrara.physics@gmail.com]
 *  @date 
 *  @details
 */
#pragma once

#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <random>
#include <mutex>

namespace dicebox
{
    class Random
    {
    public:
        Random(const Random&)         = delete;     // delete copy constructor
        void operator=(const Random&) = delete;     // delete operator
        
        // get an instance of the generator
        static Random& Get();
        // delete an instance
        static void Delete();

        // generators
        static int UniformInt()     { return Get()._uniform_int(); }
        static double UniformReal() { return Get()._uniform_real(); }
        static double Normal()      { return Get()._normal(); }

    private:
        int _uniform_int();
        double _uniform_real();
        double _normal();

    private:
        Random();
        ~Random();
        // instance for this singleton
        static Random sInstance;
        static std::mutex sMutex;
        // Mersenne twister
        std::random_device sDevice;
        std::mt19937 sGenerator;
        // distributions
        // std::uniform_int_distribution<> sUniformInt;
        // std::uniform_real_distribution<> sUniformReal;

    };
}