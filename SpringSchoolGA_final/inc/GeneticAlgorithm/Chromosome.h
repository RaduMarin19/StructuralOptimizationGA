#pragma once

#include <vector>
#include <Services/RandomNumbersGenerator.h>
#include <functional>

struct Chromosome {
	std::vector<bool> genes;

	bool operator==(const Chromosome& other) const {
		return genes == other.genes;
	}

	void Crossover(Chromosome& other)
	{
		size_t numberOfGenes = genes.size();
		int point = RandomNumbersGenerator::GenerateIntegerNumberInRange(1, numberOfGenes - 1);

		for (size_t i = point; i < numberOfGenes; ++i)
			std::swap(genes[i], other.genes[i]);
	}

	void TwoPointCrossover(Chromosome& other)
	{
		size_t numberOfGenes = genes.size();

		int indexOne = RandomNumbersGenerator::GenerateIntegerNumberInRange(0, numberOfGenes - 1);
		int indexTwo = indexOne;
		while (indexTwo == indexOne)
			indexTwo = RandomNumbersGenerator::GenerateIntegerNumberInRange(0, numberOfGenes - 1);

		if (indexOne > indexTwo)
			std::swap(indexOne, indexTwo);

		for (int i = indexOne; i <= indexTwo; ++i)
			std::swap(genes[i], other.genes[i]);
	}

	void UniformCrossover(Chromosome& other)
	{
		size_t numberOfGenes = genes.size();

		for (size_t i = 0; i < numberOfGenes; ++i)
		{
			if (RandomNumbersGenerator::GenerateRealNumberInRange(0.0, 1.0) < 0.5)
				std::swap(genes[i], other.genes[i]);
		}
	}

	void Mutation(double mutationProbability)
	{
		size_t numberOfGenes = genes.size();
		std::vector<double> randomNumbers = RandomNumbersGenerator::GenerateRealNumbers(0.0, 1.0, numberOfGenes);

		for (size_t i = 0; i < numberOfGenes; ++i)
		{
			if (randomNumbers[i] < mutationProbability)
				genes[i] = !genes[i];
		}
	}

	void InverseMutation(double mutationProbability)
	{
		size_t numberOfGenes = genes.size();

		if (numberOfGenes < 2)
			return;

		int indexOne = RandomNumbersGenerator::GenerateIntegerNumberInRange(0, numberOfGenes - 1);
		int indexTwo = indexOne;

		while (indexTwo == indexOne)
			indexTwo = RandomNumbersGenerator::GenerateIntegerNumberInRange(0, numberOfGenes - 1);

		if (indexOne > indexTwo)
			std::swap(indexOne, indexTwo);

		std::vector<double> randomNumbers = RandomNumbersGenerator::GenerateRealNumbers(0.0, 1.0, numberOfGenes);

		for (int i = indexOne; i <= indexTwo; ++i)
		{
			if (randomNumbers[i] < mutationProbability)
				genes[i] = !genes[i];
		}
	}

	void SwapMutation(double mutationProbability)
	{
		if (RandomNumbersGenerator::GenerateRealNumberInRange(0.0, 1.0) >= mutationProbability)
			return;

		int numberOfGenes = static_cast<int>(genes.size());

		int indexOne = RandomNumbersGenerator::GenerateIntegerNumberInRange(0, numberOfGenes - 1);
		int indexTwo = indexOne;

		while (indexTwo == indexOne)
			indexTwo = RandomNumbersGenerator::GenerateIntegerNumberInRange(0, numberOfGenes - 1);

		std::swap(genes[indexOne], genes[indexTwo]);
	}
};

namespace std {
	template <>
	struct hash<Chromosome> {
		size_t operator()(const Chromosome& c) const {
			size_t hash = 0;
			for (bool gene : c.genes) {
				hash <<= 1;
				hash ^= gene ? 1 : 0;
			}
			return hash;
		}
	};
}