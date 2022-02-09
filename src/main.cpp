/**
 * 
 */
#include "../include/Random.hh"
#include <stdio.h>

int main(int argc, char** argv)
{
    double number = dicebox::Random::Normal();

    std::cout << number << std::endl;
    std::cin.get();
}