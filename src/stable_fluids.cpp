#include "stable_fluids.hpp"

#include <cmath>
#include <iostream>

/**
 *
 *  Reference : https://github.com/finallyjustice/stablefluids/blob/master/GridStableFluid2D/GridStableSolver.cpp
 *
 */

StableFluidSimulator::StableFluidSimulator() : 
	rowSize(512), 
	colSize(512), 
	totalSize(rowSize * colSize), 
	minX(1.0f),
    maxX(rowSize - 1.0f),
    minY(1.0f),
    maxY(colSize - 1.0f),
	visc(0.0f),
    diff(0.0f),
    vorticity(1.0f),
    timeStep(1.0f)
{
	vx.resize(totalSize);
	vy.resize(totalSize);
	vx0.resize(totalSize);
	vy0.resize(totalSize);
	d.resize(totalSize);
	d0.resize(totalSize);
	px.resize(totalSize);
	py.resize(totalSize);
	div.resize(totalSize);
	p.resize(totalSize);
	vort.resize(totalSize);
	absVort.resize(totalSize);
	gradVortX.resize(totalSize);
	gradVortY.resize(totalSize);
	lenGrad.resize(totalSize);
	vcfx.resize(totalSize);
	vcfy.resize(totalSize);

	for(size_t i = 0; i < rowSize; ++i)
    {
        for(size_t j=0; j < colSize; ++j)
        {
            px[Get1DIndex(i, j)] = (float)i + 0.5f;
            py[Get1DIndex(i, j)] = (float)j + 0.5f;
        }
    }
}

void StableFluidSimulator::Reset()
{
	for(size_t i = 0; i < rowSize; ++i)
    {
        vx[i] = 0.0f;
        vy[i] = 0.0f;
        d[i] = 0.0f;
    }
}

void StableFluidSimulator::SetBoundary(std::vector<float>& value, BoundaryType boundaryType)
{
	float rowMultiplier = 1.0f;
	float colMultiplier = 1.0f;

	switch (boundaryType)
	{
	case BoundaryType::Density :
		break;
	case BoundaryType::VelocityAlongX :
		colMultiplier = -1.0f;
		break;
	case BoundaryType::VelocityAlongY :
		rowMultiplier = -1.0f;
		break;
	default:
		break;
	}

	for(size_t i = 1; i <= rowSize-2; ++i)
	{
		value[Get1DIndex(i, 0)] = rowMultiplier * value[Get1DIndex(i, 1)];
		value[Get1DIndex(i, colSize-1)] = rowMultiplier * value[Get1DIndex(i, colSize-2)];
	}
	for(size_t j = 1; j <= colSize-2; ++j)
	{
		value[Get1DIndex(0, j)] = colMultiplier * value[Get1DIndex(1, j)];
		value[Get1DIndex(rowSize-1, j)] = colMultiplier * value[Get1DIndex(rowSize-2, j)];
	}

	value[Get1DIndex(0, 0)] = (value[Get1DIndex(0, 1)]+value[Get1DIndex(1, 0)])/2;
    value[Get1DIndex(rowSize-1, 0)] = (value[Get1DIndex(rowSize-2, 0)]+value[Get1DIndex(rowSize-1, 1)])/2;
    value[Get1DIndex(0, colSize-1)] = (value[Get1DIndex(0, colSize-2)]+value[Get1DIndex(1, colSize-1)])/2;
    value[Get1DIndex(rowSize-1, colSize-1)] = (value[Get1DIndex(rowSize-2, colSize-1)]+value[Get1DIndex(rowSize-1, colSize-2)])/2;
}

void StableFluidSimulator::Projection()
{
	for(size_t i=1; i<=rowSize-2; ++i)
    {
        for(size_t j=1; j<=colSize-2; ++j)
        {
            div[Get1DIndex(i, j)] = 0.5f * (vx[Get1DIndex(i+1, j)]-vx[Get1DIndex(i-1, j)]+vy[Get1DIndex(i, j+1)]-vy[Get1DIndex(i, j-1)]);
            p[Get1DIndex(i, j)] = 0.0f;;
        }
    }
    SetBoundary(div, BoundaryType::Density);
    SetBoundary(p, BoundaryType::Density);

    //projection iteration
    for(size_t k=0; k<20; k++)
    {
        for(size_t i=1; i<=rowSize-2; ++i)
        {
            for(size_t j=1; j<=colSize-2; ++j)
            {
                p[Get1DIndex(i, j)] = (p[Get1DIndex(i+1, j)]+p[Get1DIndex(i-1, j)]+p[Get1DIndex(i, j+1)]+p[Get1DIndex(i, j-1)]-div[Get1DIndex(i, j)])/4.0f;
            }
        }
        SetBoundary(p, BoundaryType::Density);
    }

    //velocity minus grad of Pressure
    for(size_t i=1; i<=rowSize-2; ++i)
    {
        for(size_t j=1; j<=colSize-2; ++j)
        {
            vx[Get1DIndex(i, j)] -= 0.5f*(p[Get1DIndex(i+1, j)]-p[Get1DIndex(i-1, j)]);
            vy[Get1DIndex(i, j)] -= 0.5f*(p[Get1DIndex(i, j+1)]-p[Get1DIndex(i, j-1)]);
        }
    }
    SetBoundary(vx, BoundaryType::VelocityAlongX);
    SetBoundary(vy, BoundaryType::VelocityAlongY);
}

void StableFluidSimulator::Advection(std::vector<float>& value, std::vector<float>& value0, 
		std::vector<float>& u, std::vector<float>& v, BoundaryType boundaryType)
{
	float oldX;
    float oldY;
    int i0;
    int i1;
    int j0;
    int j1;
    float wL;
    float wR;
    float wB;
    float wT;

    for(size_t i = 1; i <= rowSize-2; ++i)
    {
        for(size_t j = 1; j <= colSize-2; ++j)
        {
            oldX = px[Get1DIndex(i, j)] - u[Get1DIndex(i, j)]*timeStep;
            oldY = py[Get1DIndex(i, j)] - v[Get1DIndex(i, j)]*timeStep;

            if(oldX < minX) oldX = minX;
            if(oldX > maxX) oldX = maxX;
            if(oldY < minY) oldY = minY;
            if(oldY > maxY) oldY = maxY;

            i0 = (int)(oldX-0.5f);
            j0 = (int)(oldY-0.5f);
            i1 = i0+1;
            j1 = j0+1;
            
            wL = px[Get1DIndex(i1, j0)]-oldX;
            wR = 1.0f-wL;
            wB = py[Get1DIndex(i0, j1)]-oldY;
            wT = 1.0f-wB;

            value[Get1DIndex(i, j)] = wB*(wL*value0[Get1DIndex(i0, j0)]+wR*value0[Get1DIndex(i1, j0)])+
                                wT*(wL*value0[Get1DIndex(i0, j1)]+wR*value0[Get1DIndex(i1, j1)]);
        }
    }
    
    SetBoundary(value, boundaryType);
}

void StableFluidSimulator::Diffusion(std::vector<float>& value, std::vector<float>& value0, 
		float rate, BoundaryType boundaryType)
{
    for(size_t i = 0; i < totalSize; ++i)
	{
		value[i] = 0.0f;
	}
    float a = rate*timeStep;

    for(size_t k = 0; k < 20; ++k)
    {
        for(size_t i = 1; i <= rowSize-2; ++i)
        {
            for(size_t j = 1; j <= colSize-2; ++j)
            {
                value[Get1DIndex(i, j)] = (value0[Get1DIndex(i, j)]+a*(value[Get1DIndex(i+1, j)]+value[Get1DIndex(i-1, j)]+value[Get1DIndex(i, j+1)]+value[Get1DIndex(i, j-1)])) / (4.0f*a+1.0f);
            }
        }

        SetBoundary(value, boundaryType);
    }
}

void StableFluidSimulator::VortConfinement()
{
    for(size_t i = 1; i <= rowSize-2; ++i)
    {
        for(size_t j = 1; j <= colSize-2; ++j)
        {
            vort[Get1DIndex(i, j)] = 0.5f*(vy[Get1DIndex(i+1, j)]-vy[Get1DIndex(i-1, j)]-vx[Get1DIndex(i, j+1)]+vx[Get1DIndex(i, j-1)]);
            if(vort[Get1DIndex(i, j)] >= 0.0f) absVort[Get1DIndex(i, j)] = vort[Get1DIndex(i, j)];
            else absVort[Get1DIndex(i, j)] = -vort[Get1DIndex(i, j)];
        }
    }
    SetBoundary(vort, BoundaryType::Density);
    SetBoundary(absVort, BoundaryType::Density);

    for(size_t i = 1; i <= rowSize-2; ++i)
    {
        for(size_t j = 1; j <= colSize-2; ++j)
        {
            gradVortX[Get1DIndex(i, j)] = 0.5f*(absVort[Get1DIndex(i+1, j)]-absVort[Get1DIndex(i-1, j)]);
            gradVortY[Get1DIndex(i, j)] = 0.5f*(absVort[Get1DIndex(i, j+1)]-absVort[Get1DIndex(i, j-1)]);
            lenGrad[Get1DIndex(i, j)] = std::sqrt(gradVortX[Get1DIndex(i, j)]*gradVortX[Get1DIndex(i, j)]+gradVortY[Get1DIndex(i, j)]*gradVortY[Get1DIndex(i, j)]);
            if(lenGrad[Get1DIndex(i, j)] < 0.01f)
            {
                vcfx[Get1DIndex(i, j)] = 0.0f;
                vcfy[Get1DIndex(i, j)] = 0.0f;
            }
            else
            {
                vcfx[Get1DIndex(i, j)] = gradVortX[Get1DIndex(i, j)] / lenGrad[Get1DIndex(i, j)];
                vcfy[Get1DIndex(i, j)] = gradVortY[Get1DIndex(i, j)] / lenGrad[Get1DIndex(i, j)];
            }
        }
    }
    SetBoundary(vcfx, BoundaryType::Density);
    SetBoundary(vcfy, BoundaryType::Density);

    for(size_t i = 1; i <= rowSize-2; ++i)
    {
        for(size_t j = 1; j <= colSize-2; ++j)
        {
            vx[Get1DIndex(i, j)] += vorticity * (vcfy[Get1DIndex(i, j)] * vort[Get1DIndex(i, j)]);
            vy[Get1DIndex(i, j)] += vorticity * (-vcfx[Get1DIndex(i, j)] * vort[Get1DIndex(i, j)]);
        }
    }

    SetBoundary(vx, BoundaryType::VelocityAlongX);
    SetBoundary(vy, BoundaryType::VelocityAlongY);
}

void StableFluidSimulator::AddSource()
{
    for(size_t i = 1; i <= rowSize-2; ++i)
    {
        for(size_t j = 1; j <= colSize-2; ++j)
        {
            int index = Get1DIndex(i, j);
            vx[index] += vx0[index];
            vy[index] += vy0[index];
            d[index] += d0[index];
        }
    }

    SetBoundary(vx, BoundaryType::VelocityAlongX);
    SetBoundary(vy, BoundaryType::VelocityAlongY);
    SetBoundary(d, BoundaryType::Density);
}

void StableFluidSimulator::AnimateVel()
{
    if(diff > 0.0f)
    {
        std::swap(vx0, vx);
        std::swap(vy0, vy);
        Diffusion(vx, vx0, diff, BoundaryType::VelocityAlongX);
        Diffusion(vy, vy0, diff, BoundaryType::VelocityAlongY);
    }

    Projection();

    std::swap(vx0, vx);
    std::swap(vy0, vy);
    Advection(vx, vx0, vx0, vy0, BoundaryType::VelocityAlongX);
    Advection(vy, vy0, vx0, vy0, BoundaryType::VelocityAlongY);

    Projection();
}

void StableFluidSimulator::AnimateDen()
{
    if(visc > 0.0f)
    {
        std::swap(d0, d);
        Diffusion(d, d0, visc, BoundaryType::Density);
    }
    std::swap(d0, d);
    Advection(d, d0, vx, vy, BoundaryType::Density);
}
