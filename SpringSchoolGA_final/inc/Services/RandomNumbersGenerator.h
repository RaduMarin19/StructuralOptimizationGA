#pragma once

#include <vector>
#include <random>

class RandomNumbersGenerator
{
public:
	static int GenerateIntegerNumberInRange(int lowerBound, int upperBound);
	static double GenerateRealNumberInRange(double lowerBound, double upperBound);
	static std::vector<double> GenerateRealNumbers(double lowerBound, double upperBound, size_t size);

private:
	static std::mt19937& GetGenerator();
};
