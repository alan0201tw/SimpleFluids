#include "stable_fluids.hpp"

#include <cmath>
#include <iostream>

/**
 *
 *  Reference : https://github.com/finallyjustice/stablefluids/blob/master/GridStableFluid2D/GridStableSolver.cpp
 *
 */

StableFluidSimulator::StableFluidSimulator(size_t _rowSize, size_t _colSize) : 
	rowSize(_rowSize), 
	colSize(_colSize), 
	totalSize((rowSize+2)*(colSize+2)), 
	minX(1.0f),
    maxX(rowSize + 1.0f),
    minY(1.0f),
    maxY(colSize + 1.0f),
	visc(0.0f),
    diff(0.0f),
    vorticity(0.0f),
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
    tx0.resize(totalSize);
    ty0.resize(totalSize);
    tx.resize(totalSize);
    ty.resize(totalSize);

    for(size_t i = 0; i < totalSize; ++i)
    {
        px[i] = 0.0f;
        py[i] = 0.0f;
        vx[i] = 0.0f;
        vy[i] = 0.0f;
        vx0[i] = 0.0f;
        vy0[i] = 0.0f;
        d[i] = 0.0f;
        d0[i] = 0.0f;
        p[i] = 0.0f;
        div[i] = 0.0f;
        tx0[i] = 0.0f;
        ty0[i] = 0.0f;
        tx[i] = 0.0f;
        ty[i] = 0.0f;
    }

	for(size_t i = 0; i < rowSize+2; ++i)
    {
        for(size_t j=0; j < colSize+2; ++j)
        {
            px[Get1DIndex(i, j)] = (float)i + 0.5f;
            py[Get1DIndex(i, j)] = (float)j + 0.5f;

            tx[Get1DIndex(i, j)] = (float)i;
            ty[Get1DIndex(i, j)] = (float)j;
        }
    }
}

void StableFluidSimulator::Reset()
{
	for(size_t i = 0; i < totalSize; ++i)
    {
        vx[i] = 0.0f;
        vy[i] = 0.0f;
        d[i] = 0.0f;
    }
}

void StableFluidSimulator::SetBoundary(std::vector<float>& value, BoundaryType boundaryType)
{
	float xMultiplier = 1.0f;
	float yMultiplier = 1.0f;

	switch (boundaryType)
	{
	case BoundaryType::Density :
		break;
	case BoundaryType::VelocityAlongX :
		xMultiplier = -1.0f;
		break;
	case BoundaryType::VelocityAlongY :
		yMultiplier = -1.0f;
		break;
	default:
		break;
	}

    for(size_t i = 1; i <= rowSize; i++) 
    {
        value[Get1DIndex(0, i)]           = xMultiplier * value[Get1DIndex(1,i)];
        value[Get1DIndex(rowSize+1,i)]    = xMultiplier * value[Get1DIndex(rowSize,i)];
        value[Get1DIndex(i,0  )]          = yMultiplier * value[Get1DIndex(i,1)];
        value[Get1DIndex(i,rowSize+1)]    = yMultiplier * value[Get1DIndex(i,rowSize)];
    }
    value[Get1DIndex(0, 0)]                 = 0.5f*(value[Get1DIndex(1, 0)]+value[Get1DIndex(0, 1)]);
    value[Get1DIndex(0, rowSize+1)]         = 0.5f*(value[Get1DIndex(1, rowSize+1)]+value[Get1DIndex(0, rowSize)]);
    value[Get1DIndex(rowSize+1, 0)]         = 0.5f*(value[Get1DIndex(rowSize, 0)]+value[Get1DIndex(rowSize+1, 1)]);
    value[Get1DIndex(rowSize+1, rowSize+1)] = 0.5f*(value[Get1DIndex(rowSize, rowSize+1)]+value[Get1DIndex(rowSize+1, rowSize)]);
}

void StableFluidSimulator::LinearSolve(std::vector<float>& value, std::vector<float>& value0,
    float a, float c, BoundaryType boundaryType)
{
    for(size_t iteration=0; iteration < 20; ++iteration) 
    {
        for(size_t i = 1; i <= rowSize; ++i)
        { 
            for(size_t j = 1; j <= colSize; ++j)
            {
                value[Get1DIndex(i, j)] = (value0[Get1DIndex(i, j)] + a*(value[Get1DIndex(i-1, j)]+value[Get1DIndex(i+1, j)]+value[Get1DIndex(i, j-1)]+value[Get1DIndex(i, j+1)]))/c;
            } 
        }

        SetBoundary(value, boundaryType);
    }
}

void StableFluidSimulator::Projection()
{
    for(size_t i=1; i <= rowSize; ++i)
    { 
        for(size_t j=1; j <= colSize; ++j)
        {
            div[Get1DIndex(i, j)] = -0.5f*(vx[Get1DIndex(i+1,j)]-vx[Get1DIndex(i-1,j)]+vy[Get1DIndex(i,j+1)]-vy[Get1DIndex(i,j-1)]);
            p[Get1DIndex(i,j)] = 0;
        }
    }
    SetBoundary(div, BoundaryType::Density);
    SetBoundary(p, BoundaryType::Density);

    LinearSolve(p, div, 1.0, 4.0, BoundaryType::Density);

    //velocity minus grad of Pressure
    for(size_t i = 1; i <= rowSize; ++i)
    {
        for(size_t j = 1; j <= colSize; ++j)
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
    float iL;
    float iR;
    float iT;
    float iB;

    for(size_t i = 1; i <= rowSize; ++i)
    {
        for(size_t j = 1; j <= colSize; ++j)
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
            
            iR = oldX - px[Get1DIndex(i0, j0)];
            iT = oldY - py[Get1DIndex(i0, j0)];
            iL = 1.0f - iR;
            iB = 1.0f - iT;

            value[Get1DIndex(i, j)] = 
                iB * (iL*value0[Get1DIndex(i0, j0)] + iR*value0[Get1DIndex(i1, j0)]) +
                iT * (iL*value0[Get1DIndex(i0, j1)] + iR*value0[Get1DIndex(i1, j1)]);
        }
    }
    
    SetBoundary(value, boundaryType);
}

void StableFluidSimulator::Diffusion(std::vector<float>& value, std::vector<float>& value0, 
		float rate, BoundaryType boundaryType)
{
    float a = timeStep * diff; 
    LinearSolve(value, value0, a, 1 + 4*a, boundaryType);
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
    for(size_t i = 0; i <= totalSize; ++i)
    {
        vx[i] += vx0[i];
        vy[i] += vy0[i];
        d[i] += d0[i];
    }

    SetBoundary(vx, BoundaryType::VelocityAlongX);
    SetBoundary(vy, BoundaryType::VelocityAlongY);
    SetBoundary(d, BoundaryType::Density);
}

void StableFluidSimulator::AnimateVel()
{
    std::swap(vx0, vx); 
    std::swap(vy0, vy); 
    Diffusion(vx, vx0, visc, BoundaryType::VelocityAlongX);
    Diffusion(vy, vy0, visc, BoundaryType::VelocityAlongY);

    Projection();

    std::swap(vx0, vx);
    std::swap(vy0, vy);
    Advection(vx, vx0, vx0, vy0, BoundaryType::VelocityAlongX);
    Advection(vy, vy0, vx0, vy0, BoundaryType::VelocityAlongY);

    Projection();
}

void StableFluidSimulator::AnimateDen()
{
    std::swap(d0, d);
    Advection(d, d0, vx, vy, BoundaryType::Density);

    std::swap(d0, d);
    Diffusion(d, d0, diff, BoundaryType::Density);
}

void StableFluidSimulator::AnimateTex()
{
    std::swap(tx0, tx);
    std::swap(ty0, ty);
    Advection(tx, tx0, vx, vy, BoundaryType::Density);
    Advection(ty, ty0, vx, vy, BoundaryType::Density);

    std::swap(tx0, tx);
    std::swap(ty0, ty);
    Diffusion(tx, tx0, diff, BoundaryType::Density);
    Diffusion(ty, ty0, diff, BoundaryType::Density);
}