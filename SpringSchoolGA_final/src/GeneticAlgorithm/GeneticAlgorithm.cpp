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

	for (int index = 0; index < m_numberOfEpochs; ++index)
	{
		std::cout << std::endl << "Epoch: " << index + 1 << std::endl;

		m_chromFitnessValues = CalculateFitnessValues();

		std::cout << "Calculated Fitness Values\n";

		/*if (SELECTION_TYPE == "roulette_wheel")*/ Selection();
		//else TournamentSelection();
		std::cout << "Selection done! - " << SELECTION_TYPE << "\n";
		Crossover();
		std::cout << "Crossover done! - " << CROSSOVER_TYPE << "\n";
		Mutation();
		std::cout << "Mutation  done! - " << MUTATION_TYPE << "\n";


		WriteWinners(index);
	}
}

Chromosome GeneticAlgorithm::GetWinnerIndividual()
{
	double maxValue = 0.0;
	Chromosome winner;

	for (const auto value : m_chromFitnessValues)
		if (value.second > maxValue)
		{
			maxValue = value.second;
			winner = value.first;
		}

	return winner;
}

void GeneticAlgorithm::InitializePopulation()
{
	m_individual = std::shared_ptr<IIndividual>(m_createIndividual());

	if (!m_individual) {
		std::cerr << "[FATAL] m_individual is nullptr after creation.\n";
		exit(EXIT_FAILURE);
	}

	auto individualPtr = std::dynamic_pointer_cast<Individual>(m_individual);
	size_t geneLength = individualPtr->GetBuilding()->GetCubesExistence().size();

	for (int index = 0; index < m_populationSize; ++index)
	{
		std::cout << "Created individual " << index + 1 << "\n";
		Chromosome chrom = individualPtr->GetChromosome();
		m_chromPopulation.push_back(chrom);
		m_chromWorkingPopulation.push_back(chrom);
	}
}

std::unordered_map<Chromosome, double> GeneticAlgorithm::CalculateFitnessValues()
{
	static int times = 0;
	std::unordered_map<Chromosome, double> fitnessValues;

	auto individualPtr = std::dynamic_pointer_cast<Individual>(m_individual);
	if (!individualPtr) {
		std::cerr << "[FATAL] m_individual cast to Individual failed.\n";
		exit(EXIT_FAILURE);
	}

	size_t expectedGeneSize = individualPtr->GetBuilding()->GetCubesExistence().size();
	int generationCount = 0;

	for (const auto& individual : m_chromWorkingPopulation)
	{
		if (m_calculatedFitnessValues.find(individual) != m_calculatedFitnessValues.end()) {
			fitnessValues[individual] = m_calculatedFitnessValues[individual];
			times++;
			std::cout << "Retrieved old individual fitness(" << times << ")\n";
			continue;
		}

		std::cerr << "Evaluating chromosome #" << generationCount++ << "\n";

		if (individual.genes.size() != expectedGeneSize) {
			std::cerr << "[FATAL] Gene size mismatch at generation " << times << ": "
				<< individual.genes.size() << " vs expected " << expectedGeneSize << "\n";
			exit(EXIT_FAILURE);
		}

		// Apply gene mask only after size is validated
		individualPtr->GetBuilding()->EliminateCubesBasedOnCubesExistence(individual.genes);
		individualPtr->GetBuilding()->AddCubesBasedOnCubesExistence(individual.genes);

		double value = individualPtr->Evaluate();
		fitnessValues[individual] = value;
		m_calculatedFitnessValues[individual] = value;
		std::cout << m_calculatedFitnessValues.size()<<"\n";
	}

	return fitnessValues;
}


double GeneticAlgorithm::CalculateSumOfFitnessValues()
{
	double sum{};
	for (const auto& individual : m_chromWorkingPopulation)
	{
		sum += m_chromFitnessValues[individual];
	}

	return sum;
}

std::vector<double> GeneticAlgorithm::CalculateProbabilityOfSelection()
{
	std::vector<double> probabilityOfSelectionVector;
	double sum = CalculateSumOfFitnessValues();

	for (const auto& individual : m_chromWorkingPopulation)
	{
		probabilityOfSelectionVector.emplace_back(m_chromFitnessValues[individual] / sum);
	}

	return probabilityOfSelectionVector;
}

//std::vector<double> GeneticAlgorithm::CalculateProbabilityOfRankSelection()
//{
//	std::vector<std::shared_ptr<IIndividual>> newPopulation;
//	std::vector<std::pair < IIndividual*, double >> rankProbabilities{ m_chromFitnessValues.begin(), m_chromFitnessValues.end() };
//	std::vector<double> probabilities;
//
//	size_t ranksLength = rankProbabilities.size();
//	probabilities.reserve(ranksLength);
//
//	std::sort(rankProbabilities.begin(), rankProbabilities.end(), [](const auto& lhs, const auto& rhs) {
//		return lhs.second > rhs.second;
//		});
//
//	int sumOfRanks = ranksLength * (ranksLength + 1) / 2;
//
//	for (int index = 0; index < ranksLength; ++index)
//	{
//		rankProbabilities[index].second = index + 1;
//	}
//
//	for (auto& pair : rankProbabilities)
//	{
//		pair.second /= sumOfRanks;
//	}
//
//	std::transform(rankProbabilities.begin(), rankProbabilities.end(), std::back_inserter(probabilities),
//		[](const std::pair<IIndividual*, double> p) { return p.second; });
//
//	return probabilities;
//}

std::vector<double> GeneticAlgorithm::CalcutateCumulativeProbabilityOfSelection()
{
	std::vector<double> cumulativeProbabilityOfSelectionVector;
	std::vector<double> probabilityOfSelectionVector;

	/*if (SELECTION_TYPE == "roulette_wheel")*/ probabilityOfSelectionVector = CalculateProbabilityOfSelection();
	//else probabilityOfSelectionVector = CalculateProbabilityOfRankSelection();

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
	std::vector<Chromosome> newPopulation;

	std::vector<double> cumulativeProbabilityOfSelectionVector = CalcutateCumulativeProbabilityOfSelection();
	std::vector<double> randomNumbers = RandomNumbersGenerator::GenerateRealNumbers(LOWER_BOUND, UPPER_BOUND, m_populationSize);

	for (const auto& randomNumber : randomNumbers)
	{
		if (IsGraterThan(randomNumber, LOWER_BOUND) &&
			IsLessThanOrEqualTo(randomNumber, cumulativeProbabilityOfSelectionVector[0]))
		{
			newPopulation.push_back(m_chromWorkingPopulation[0]);
			continue;
		}

		for (size_t probabilityIndex = 0; probabilityIndex < cumulativeProbabilityOfSelectionVector.size() - 1; ++probabilityIndex)
		{
			if (IsGraterThan(randomNumber, cumulativeProbabilityOfSelectionVector[probabilityIndex]) &&
				IsLessThanOrEqualTo(randomNumber, cumulativeProbabilityOfSelectionVector[probabilityIndex + 1]))
			{
				newPopulation.push_back(m_chromWorkingPopulation[probabilityIndex + 1]);
				break;
			}
		}
	}

	m_chromWorkingPopulation = newPopulation;
}


//void GeneticAlgorithm::TournamentSelection()
//{
//	std::vector<std::shared_ptr<IIndividual>> newPopulation{};
//	size_t populationSize = m_workingPopulation.size();
//	size_t tournamentSize = populationSize / 10;
//
//	for (int index = 0; index < populationSize; ++index)
//	{
//		std::vector<double> randomNumbers = RandomNumbersGenerator::GenerateRealNumbers(0, populationSize, tournamentSize);
//		std::vector<int> tournamentIndices(tournamentSize);
//
//		std::transform(randomNumbers.begin(), randomNumbers.end(), tournamentIndices.begin(),
//			[](double d) { return static_cast<int>(d); });
//
//		std::vector<std::pair<std::shared_ptr<IIndividual>, double>> tournamentPopulation{};
//
//		for (auto index : tournamentIndices)
//		{
//			Chromosome individual = m_chromWorkingPopulation[index];
//			double fitnessValue = m_chromFitnessValues[individual];
//			tournamentPopulation.push_back({ individual, fitnessValue });
//		}
//
//		double bestFitnessValue = tournamentPopulation[0].second;
//		int bestIndex = 0;
//
//		for (int index = 0; index < tournamentSize; ++index)
//		{
//			if (tournamentPopulation[index].second > bestFitnessValue)
//			{
//				bestFitnessValue = tournamentPopulation[index].second;
//				bestIndex = index;
//			}
//		}
//
//		newPopulation.push_back(tournamentPopulation[bestIndex].first);
//	}
//
//	m_workingPopulation = newPopulation;
//}

void GeneticAlgorithm::Crossover()
{
	std::vector<Chromosome> selectedPopulationForCrossover;
	std::vector<Chromosome> newPopulation;

	std::vector<double> randomNumbers = RandomNumbersGenerator::GenerateRealNumbers(LOWER_BOUND, UPPER_BOUND, m_populationSize);

	for (size_t index = 0; index < m_populationSize; ++index)
	{
		if (randomNumbers[index] < m_crossoverProbability)
		{
			selectedPopulationForCrossover.push_back(m_chromWorkingPopulation[index]);
		}
	}

	if (selectedPopulationForCrossover.size() % 2 != 0)
	{
		selectedPopulationForCrossover.pop_back();
	}

	for (size_t index = 0; index < selectedPopulationForCrossover.size(); index += 2)
	{
		if (CROSSOVER_TYPE == "single_point") selectedPopulationForCrossover[index].Crossover(selectedPopulationForCrossover[index + 1]);
		else if (CROSSOVER_TYPE == "two_point") selectedPopulationForCrossover[index].TwoPointCrossover(selectedPopulationForCrossover[index + 1]);
		else selectedPopulationForCrossover[index].UniformCrossover(selectedPopulationForCrossover[index + 1]);
	}
}

void GeneticAlgorithm::Mutation()
{
	for (auto& individual : m_chromWorkingPopulation)
	{
		if (MUTATION_TYPE == "bit_flip") individual.Mutation(m_mutationProbability);
		else if (MUTATION_TYPE == "inverse") individual.InverseMutation(m_mutationProbability);
		else individual.SwapMutation(m_mutationProbability);
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
	Chromosome winner = GetWinnerIndividual();

	if (epoch == 0)
	{
		IOIndividualManager::WriteIndividualValueInFile(epoch + 1, m_chromFitnessValues[winner], false);
	}
	else
	{
		if (epoch == m_numberOfEpochs - 1)
		{
			//IIndivi
			//IOIndividualManager::WriteIndividualDetailsInFile(winner);
		}
		IOIndividualManager::WriteIndividualValueInFile(epoch + 1, m_chromFitnessValues[winner], true);
	}
}
