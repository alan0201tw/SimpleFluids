#pragma once

/**
 *
 *  Reference : https://github.com/finallyjustice/stablefluids/blob/master/GridStableFluid2D/GridStableSolver.cpp
 *
 */

#include <vector>

class StableFluidSimulator
{
public:
	StableFluidSimulator();
	~StableFluidSimulator();

private:
    const size_t rowSize;
	const size_t colSize;
	const size_t totalSize;

	const float h; // = 1.0f;

	// float minX;
	// float maxX;
	// float minY;
	// float maxY;
	/*
	minX = 1.0f;
    maxX = rowSize-1.0f;
    minY = 1.0f;
    maxY = colSize-1.0f;
	*/

	//params
	float visc;
	float diff;
	float vorticity;
	float timeStep;

	/*
    visc = 0.0f;
    diff = 0.0f;
    vorticity = 0.0f;
    timeStep = 1.0f;
	*/

	std::vector<float> vx;
	std::vector<float> vy;
	std::vector<float> vx0;
	std::vector<float> vy0;
	std::vector<float> d;
	std::vector<float> d0;
	std::vector<float> px;
	std::vector<float> py;
	std::vector<float> div;
	std::vector<float> p;
	//vorticity confinement
	std::vector<float> vort;
	std::vector<float> absVort;
	std::vector<float> gradVortX;
	std::vector<float> gradVortY;
	std::vector<float> lenGrad;
	std::vector<float> vcfx;
	std::vector<float> vcfy;
};