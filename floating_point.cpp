/**
 * @author Chirag Kaushik
 * @date August 25, 2021
 * This file contains various tests to explore the limitations and 
 * boundaries of floating point numbers.
 */

#include <iostream>

/**
 * The main method does several floating point calculations 
 * and prints their results.
 */
int main() 
{
    std::cout << "nonzero / zero:   " << (3.0f/0.0f) << std::endl;

    std::cout << "zero / zero:   " << (0.0f/0.0f) << std::endl;

    float inf = 3.0f/0.0f;
    float neginf = -3.0f/0.0f;
    std::cout << "inf / 0:   " << (inf/0.0f) << std::endl;
    std::cout << "-inf / 0:   " <<( neginf/0.0f) << std::endl;

    std::cout << "0 / inf:   " << (0.0f/inf) << std::endl;
    std::cout << "0 / -inf:   " << (0.0f/neginf) << std::endl;

    std::cout << "0 * inf:   " << (0.0f * inf) << std::endl;
    std::cout << "0 * -inf:   " << (0.0f * neginf) << std::endl;

    std::cout << "inf * inf:   " << (inf * inf) << std::endl;
    std::cout << "inf * -inf:   " << (inf * neginf) << std::endl;
    std::cout << "-inf * inf:   " << (neginf * inf) << std::endl;
    std::cout << "-inf * -inf:   " << (neginf * neginf) << std::endl;

    std::cout << "inf / inf:   " << (inf / inf) << std::endl;
    std::cout << "inf / -inf:   " << (inf / neginf) << std::endl;
    std::cout << "-inf / inf:   " << (neginf / inf) << std::endl;
    std::cout << "-inf / -inf:   " << (neginf / neginf) << std::endl;
    
    float min = 1;
    while (2 * (min / 2) != 0) 
    {
        min /= 2;
    }
    std::cout << "min:    +/-" << min << std::endl;


    float max = 1;
    while ((max * 2) / 2 != inf) 
    {
        max *= 2;
    }
    std::cout << "max:    +/-" << max << std::endl;

    float eps = 1;
    while (1 + (eps / 2) != 1) 
    {
        eps /= 2;
    }
    std::cout << "eps:    +/-" << eps << std::endl;

    return(0);
} // int main()