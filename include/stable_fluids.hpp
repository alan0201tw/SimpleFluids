#pragma once

/**
 *
 *  Reference : https://github.com/finallyjustice/stablefluids/blob/master/GridStableFluid2D/GridStableSolver.cpp
 *
 */

#include <vector>

enum class BoundaryType : short
{
	Density,
	VelocityAlongX,
	VelocityAlongY,
};

class StableFluidSimulator
{
public:
	StableFluidSimulator();

	void Reset();
    void VortConfinement();
    void AddSource();
    void AnimateVel();
    void AnimateDen();

	inline float GetBilinearFilteredDensity(size_t i, size_t j)
	{
		return (d[Get1DIndex(i-1, j-1)]
			 + d[Get1DIndex(i, j-1)] 
			 + d[Get1DIndex(i-1, j)] 
			 + d[Get1DIndex(i, j)]) * 0.25f;
	}

	void SetD0(size_t i, size_t j, float value)
	{
		d0[Get1DIndex(i, j)] = value;
	}

	void SetVX0(size_t i, size_t j, float value)
	{
		vx0[Get1DIndex(i, j)] = value;
	}

	void SetVY0(size_t i, size_t j, float value)
	{
		vy0[Get1DIndex(i, j)] = value;
	}

	void CleanBuffer()
	{
		for(size_t i = 0; i < totalSize; ++i)
		{
			vx0[i] = 0.0f;
			vy0[i] = 0.0f;
			d0[i] = 0.0f;
		}
	}

private:
    const size_t rowSize;
	const size_t colSize;
	const size_t totalSize;

	const float minX;
	const float maxX;
	const float minY;
	const float maxY;
	/*
	minX = 1.0f;
    maxX = rowSize-1.0f;
    minY = 1.0f;
    maxY = colSize-1.0f;
	*/

	//params
	const float visc;
	const float diff;
	const float vorticity;
	const float timeStep;
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

	inline size_t Get1DIndex(size_t i, size_t j) { return j * rowSize + i; }

	void SetBoundary(std::vector<float>& value, BoundaryType boundaryType);
	void Projection();
	void Advection(std::vector<float>& value, std::vector<float>& value0, 
		std::vector<float>& u, std::vector<float>& v, BoundaryType boundaryType);
    void Diffusion(std::vector<float>& value, std::vector<float>& value0, 
		float rate, BoundaryType boundaryType);
};