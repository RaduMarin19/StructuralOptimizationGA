#pragma once

#include <functional>
#include <iostream>
#include <map>
#include <unordered_map>

#include <GeneticAlgorithm/IIndividual.h>
#include <GeneticAlgorithm/Chromosome.h>

#include <Services/RandomNumbersGenerator.h>
#include <Services/IOIndividualManager.h>
#include <Services/constants.h>

class GeneticAlgorithm
{
public:
	GeneticAlgorithm(
		std::function<IIndividual* ()> createIndividual,
		size_t populationSize,
		size_t numberOfEpochs,
		double crossoverProbabillity,
		double mutationProbability);

	GeneticAlgorithm(const GeneticAlgorithm& other) = delete;
	GeneticAlgorithm(GeneticAlgorithm&& other) = delete;

	GeneticAlgorithm& operator=(const GeneticAlgorithm& other) = delete;
	GeneticAlgorithm& operator=(GeneticAlgorithm&& other) = delete;

	~GeneticAlgorithm() = default;

	void Run();

	Chromosome GetWinnerIndividual();

private:
	void InitializePopulation();

	std::unordered_map<Chromosome, double> CalculateFitnessValues();

	double CalculateSumOfFitnessValues();
	std::vector<double> CalculateProbabilityOfSelection();
	//std::vector<double> CalculateProbabilityOfRankSelection();
	std::vector<double> CalcutateCumulativeProbabilityOfSelection();

	void Selection();
	//void TournamentSelection();

	void Crossover();
	void Mutation();

	bool IsGraterThan(double value, double lowerBound) const;
	bool IsLessThanOrEqualTo(double value, double upperBound) const;

	void WriteWinners(int epoch);

private:
	std::vector<std::shared_ptr<IIndividual>> m_population;
	std::vector<std::shared_ptr<IIndividual>> m_workingPopulation;
	std::vector<Chromosome> m_chromPopulation;
	std::vector<Chromosome> m_chromWorkingPopulation;
	std::unordered_map<Chromosome, double> m_calculatedFitnessValues;
	std::shared_ptr<IIndividual> m_individual;

	std::function<IIndividual* ()> m_createIndividual;

	//std::map<IIndividual*, double> m_fitnessValues;
	std::unordered_map<Chromosome, double> m_chromFitnessValues;

	size_t m_populationSize;
	size_t m_numberOfEpochs;

	double m_crossoverProbability;
	double m_mutationProbability;
};