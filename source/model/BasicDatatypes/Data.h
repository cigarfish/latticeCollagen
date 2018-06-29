///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//  File Name:  Data.h                                                               //
//                                                                                   //
//     Author:  Andreas Buttenschoen <andreas@buttenschoen.ca>                       //
//    Created:  2016-05-10                                                           //
//                                                                                   //
//  This file is part of the CellSys7 code.                                          //
//                                                                                   //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig, and INRIA, Paris-Rocquencourt. //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////
#pragma once

#include "Matrix.h"
#include "Vector.h"
#include "Utils.h"

#include <iostream>
#include <cmath>
#include <utility>
#include <vector>
#include <string>
using namespace std;

struct Data
{
    Data() : mName("") {}
    Data(const string& name) : mName(name) {}
    Data(const Data& other) : mName(other.mName) {}

    virtual void Reset() = 0;
    virtual double operator()(void) { return 0.; };
    const string& getName() const { return mName; }

private:
    const string mName;
};

struct AverageData : public Data
{
    AverageData();
    AverageData(const string& name);
    AverageData(const AverageData& data);

    virtual void Reset();
    void Add(double value);
    void Add(double value, pair<int, int> _ids);
    double Average(void) const;
    virtual double operator()(void) const;

    double avg;
    std::size_t counter;
    // verification only remove me
    vector<pair<int, int>> ids;
};

struct StressData : public Data
{
    StressData(double& v);
    StressData(const StressData& data);

    virtual void Reset();
    void Add(const Vector3f& force, const Vector3f& distance);
    void Add(const Vector3f& force, const Vector3f& distance,
             pair<int, int> _ids);
    double Pressure(void) const;
    virtual double operator()(void) const;

    //
    // Keeps track of a cells stress data
    //
    // The virial stress is defined as
    //
    // s = 1 / V Sum f_ij x r_ij
    //
    // where V is taken to be the cell volume
    //
    BasicDatatypes::Tensor stress;
    bool empty;
    double& volume;
    vector<pair<int, int>> ids;
};
