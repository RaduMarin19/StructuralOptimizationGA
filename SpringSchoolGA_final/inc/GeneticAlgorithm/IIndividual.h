#pragma once

class IIndividual
{
public:
	virtual double Evaluate() = 0;

	virtual void Crossover(IIndividual& other) = 0;
	virtual void Mutation(double mutationProbability) = 0;

	virtual void TwoPointCrossover(IIndividual& other) = 0;
	virtual void InverseMutation(double mutationProbability) = 0;

	virtual void UniformCrossover(IIndividual& other) = 0;
	virtual void SwapMutation(double mutationProbability) = 0;
};