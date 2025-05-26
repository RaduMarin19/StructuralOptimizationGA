#include <Services/RandomNumbersGenerator.h>

std::mt19937& RandomNumbersGenerator::GetGenerator()
{
	static std::random_device rd;
	static std::mt19937 gen(rd());
	return gen;
}

int RandomNumbersGenerator::GenerateIntegerNumberInRange(int lowerBound, int upperBound)
{
	std::uniform_int_distribution<int> dist(lowerBound, upperBound);
	return dist(GetGenerator());
}

double RandomNumbersGenerator::GenerateRealNumberInRange(double lowerBound, double upperBound)
{
	std::uniform_real_distribution<double> dist(lowerBound, upperBound);
	return dist(GetGenerator());
}

std::vector<double> RandomNumbersGenerator::GenerateRealNumbers(double lowerBound, double upperBound, size_t size)
{
	std::vector<double> numbers;
	numbers.reserve(size);

	for (size_t i = 0; i < size; ++i)
	{
		numbers.push_back(GenerateRealNumberInRange(lowerBound, upperBound));
	}

	return numbers;
}
