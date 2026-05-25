#pragma once
#include <string>

struct Atom
{
    std::string name;

    double x;
    double y;

    double vx;
    double vy;

    double mass;
    double radius;

    int maxBonds;
    int currentBonds;
};

struct Bond
{
    int atomA;
    int atomB;
};

struct AngleConstraint
{
    int atomA;
    int atomB;
    int atomC;

    double idealAngleDegrees;
    double strength;
};
