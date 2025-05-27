#include <GeneticAlgorithm/GeneticAlgorithm.h>

GeneticAlgorithm::GeneticAlgorithm(
	std::function<IIndividual* ()> createIndividual,
	size_t populationSize, size_t numberOfEpochs,
	double crossoverProbabillity, double mutationProbability) :
	m_createIndividual{ createIndividual },
	m_populationSize{ populationSize },
	m_numberOfEpochs{ numberOfEpochs },
	m_crossoverProbability{ crossoverProbabillity },
	m_mutationProbability{ mutationProbability }
{
}

void GeneticAlgorithm::Run()
{
	InitializePopulation();

	for (int index = 0; index < 3; ++index)
	{
		std::cout << std::endl << "Epoch: " << index + 1 << std::endl;

		m_fitnessValues = CalculateFitnessValues();

		std::cout << "Calculated Fitness Values\n";

		if (SELECTION_TYPE == "roulette_wheel") Selection();
		else TournamentSelection();
		std::cout << "Selection done! - " << SELECTION_TYPE << "\n";
		Crossover();
		std::cout << "Crossover done! - " << CROSSOVER_TYPE << "\n";
		Mutation();
		std::cout << "Mutation  done! - " << MUTATION_TYPE << "\n";


		WriteWinners(index);
	}
}

IIndividual* GeneticAlgorithm::GetWinnerIndividual()
{
	double maxValue = 0.0;
	IIndividual* winner = nullptr;

	for (const auto value : m_fitnessValues)
		if (value.second > maxValue)
		{
			maxValue = value.second;
			winner = value.first;
		}

	return winner;
}

void GeneticAlgorithm::InitializePopulation()
{
	for (int index = 0; index < m_populationSize; ++index)
	{
		std::cout << "Created individual " << index + 1 << "\n";

		m_population.push_back(std::move(std::shared_ptr<IIndividual>(m_createIndividual())));
		m_workingPopulation.push_back(m_population[index]);
	}
}

std::map<IIndividual*, double> GeneticAlgorithm::CalculateFitnessValues()
{
	static int times = 0;
	std::map<IIndividual*, double> fitnessValues;
	for (const auto& individual : m_workingPopulation)
	{
		auto& genes = std::dynamic_pointer_cast<Individual>(individual)->GetChromosome();
		if (m_calculatedFitnessValues.find(genes) != m_calculatedFitnessValues.end()) {
			fitnessValues[individual.get()] = m_calculatedFitnessValues[genes];
			times++;
			std::cout << "Retrieved old individual fitness("<<times<<")"<< std::endl;
			continue;
		}
		std::cout << "Calculated fitness\n";
		double value = individual->Evaluate();
		fitnessValues[individual.get()] = value;
		m_calculatedFitnessValues[genes] = value;
	}

	return fitnessValues;
}

double GeneticAlgorithm::CalculateSumOfFitnessValues()
{
	double sum{};
	for (const auto& individual : m_workingPopulation)
	{
		sum += m_fitnessValues[individual.get()];
	}

	return sum;
}

std::vector<double> GeneticAlgorithm::CalculateProbabilityOfSelection()
{
	std::vector<double> probabilityOfSelectionVector;
	double sum = CalculateSumOfFitnessValues();

	for (const auto& individual : m_workingPopulation)
	{
		probabilityOfSelectionVector.emplace_back(m_fitnessValues[individual.get()] / sum);
	}

	return probabilityOfSelectionVector;
}

std::vector<double> GeneticAlgorithm::CalculateProbabilityOfRankSelection()
{
	std::vector<std::shared_ptr<IIndividual>> newPopulation;
	std::vector<std::pair < IIndividual*, double >> rankProbabilities{ m_fitnessValues.begin(), m_fitnessValues.end() };
	std::vector<double> probabilities;

	size_t ranksLength = rankProbabilities.size();
	probabilities.reserve(ranksLength);

	std::sort(rankProbabilities.begin(), rankProbabilities.end(), [](const auto& lhs, const auto& rhs) {
		return lhs.second > rhs.second;
		});

	int sumOfRanks = ranksLength * (ranksLength + 1) / 2;

	for (int index = 0; index < ranksLength; ++index)
	{
		rankProbabilities[index].second = index + 1;
	}

	for (auto& pair : rankProbabilities)
	{
		pair.second /= sumOfRanks;
	}

	std::transform(rankProbabilities.begin(), rankProbabilities.end(), std::back_inserter(probabilities),
		[](const std::pair<IIndividual*, double> p) { return p.second; });

	return probabilities;
}

std::vector<double> GeneticAlgorithm::CalcutateCumulativeProbabilityOfSelection()
{
	std::vector<double> cumulativeProbabilityOfSelectionVector;
	std::vector<double> probabilityOfSelectionVector;

	if (SELECTION_TYPE == "roulette_wheel") probabilityOfSelectionVector = CalculateProbabilityOfSelection();
	else probabilityOfSelectionVector = CalculateProbabilityOfRankSelection();

	for (int currentIndividualIndex = 0; currentIndividualIndex < m_populationSize; ++currentIndividualIndex)
	{
		double probability{};

		for (int index = 0; index <= currentIndividualIndex; ++index)
		{
			probability += probabilityOfSelectionVector[index];
		}

		cumulativeProbabilityOfSelectionVector.emplace_back(probability);
	}

	return cumulativeProbabilityOfSelectionVector;
}

void GeneticAlgorithm::Selection()
{
	std::vector<std::shared_ptr<IIndividual>> newPopulation;

	std::vector<double> cumulativeProbabilityOfSelectionVector = CalcutateCumulativeProbabilityOfSelection();
	std::vector<double> randomNumbers = RandomNumbersGenerator::GenerateRealNumbers(LOWER_BOUND, UPPER_BOUND, m_populationSize);

	for (const auto& randomNumber : randomNumbers)
	{
		if (IsGraterThan(randomNumber, LOWER_BOUND) &&
			IsLessThanOrEqualTo(randomNumber, cumulativeProbabilityOfSelectionVector[0]))
		{
			newPopulation.push_back(m_workingPopulation[0]);
			continue;
		}

		for (size_t probabilityIndex = 0; probabilityIndex < cumulativeProbabilityOfSelectionVector.size() - 1; ++probabilityIndex)
		{
			if (IsGraterThan(randomNumber, cumulativeProbabilityOfSelectionVector[probabilityIndex]) &&
				IsLessThanOrEqualTo(randomNumber, cumulativeProbabilityOfSelectionVector[probabilityIndex + 1]))
			{
				newPopulation.push_back(m_workingPopulation[probabilityIndex + 1]);
				break;
			}
		}
	}

	m_workingPopulation = newPopulation;
}

void GeneticAlgorithm::TournamentSelection()
{
	std::vector<std::shared_ptr<IIndividual>> newPopulation{};
	size_t populationSize = m_workingPopulation.size();
	size_t tournamentSize = populationSize / 10;

	for (int index = 0; index < populationSize; ++index)
	{
		std::vector<double> randomNumbers = RandomNumbersGenerator::GenerateRealNumbers(0, populationSize, tournamentSize);
		std::vector<int> tournamentIndices(tournamentSize);

		std::transform(randomNumbers.begin(), randomNumbers.end(), tournamentIndices.begin(),
			[](double d) { return static_cast<int>(d); });

		std::vector<std::pair<std::shared_ptr<IIndividual>, double>> tournamentPopulation{};

		for (auto index : tournamentIndices)
		{
			std::shared_ptr<IIndividual> individual = m_workingPopulation[index];
			double fitnessValue = m_fitnessValues[individual.get()];
			tournamentPopulation.push_back({ individual, fitnessValue });
		}

		double bestFitnessValue = tournamentPopulation[0].second;
		int bestIndex = 0;

		for (int index = 0; index < tournamentSize; ++index)
		{
			if (tournamentPopulation[index].second > bestFitnessValue)
			{
				bestFitnessValue = tournamentPopulation[index].second;
				bestIndex = index;
			}
		}

		newPopulation.push_back(tournamentPopulation[bestIndex].first);
	}

	m_workingPopulation = newPopulation;
}

void GeneticAlgorithm::Crossover()
{
	std::vector<std::shared_ptr<IIndividual>> selectedPopulationForCrossover;
	std::vector<std::shared_ptr<IIndividual>> newPopulation;

	std::vector<double> randomNumbers = RandomNumbersGenerator::GenerateRealNumbers(LOWER_BOUND, UPPER_BOUND, m_populationSize);

	for (size_t index = 0; index < m_populationSize; ++index)
	{
		if (randomNumbers[index] < m_crossoverProbability)
		{
			selectedPopulationForCrossover.push_back(m_workingPopulation[index]);
		}
	}

	if (selectedPopulationForCrossover.size() % 2 != 0)
	{
		selectedPopulationForCrossover.pop_back();
	}

	for (size_t index = 0; index < selectedPopulationForCrossover.size(); index += 2)
	{
		if (CROSSOVER_TYPE == "single_point") selectedPopulationForCrossover[index]->Crossover(*selectedPopulationForCrossover[index + 1]);
		else if (CROSSOVER_TYPE == "two_point") selectedPopulationForCrossover[index]->TwoPointCrossover(*selectedPopulationForCrossover[index + 1]);
		else selectedPopulationForCrossover[index]->UniformCrossover(*selectedPopulationForCrossover[index + 1]);
	}
}

void GeneticAlgorithm::Mutation()
{
	for (auto& individual : m_workingPopulation)
	{
		if (MUTATION_TYPE == "bit_flip") individual->Mutation(m_mutationProbability);
		else if (MUTATION_TYPE == "inverse") individual->InverseMutation(m_mutationProbability);
		else individual->SwapMutation(m_mutationProbability);
	}
}

bool GeneticAlgorithm::IsGraterThan(double value, double lowerBound) const
{
	return value > lowerBound;
}

bool GeneticAlgorithm::IsLessThanOrEqualTo(double value, double upperBound) const
{
	return value <= upperBound;
}

void GeneticAlgorithm::WriteWinners(int epoch)
{
	IIndividual* winner = GetWinnerIndividual();

	if (epoch == 0)
	{
		IOIndividualManager::WriteIndividualValueInFile(epoch + 1, m_fitnessValues[winner], false);
	}
	else
	{
		if (epoch == m_numberOfEpochs - 1)
		{
			IOIndividualManager::WriteIndividualDetailsInFile(winner);
		}
		IOIndividualManager::WriteIndividualValueInFile(epoch + 1, m_fitnessValues[winner], true);
	}
}
