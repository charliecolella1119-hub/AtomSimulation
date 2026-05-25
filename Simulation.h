#pragma once
#include <vector>
#include "Atom.h"
#include <string>

void createAtoms(std::vector<Atom>& atoms, int count);

void applyForces(std::vector<Atom>& atoms);

void applyBondForces(std::vector<Atom>& atoms,
                     const std::vector<Bond>& bonds);

void applyTemperature(std::vector<Atom>& atoms, double temperature);

void updateAtoms(std::vector<Atom>& atoms);

void updateBonds(std::vector<Atom>& atoms,
                 std::vector<Bond>& bonds);

void breakOverstretchedBonds(std::vector<Atom>& atoms,
                             std::vector<Bond>& bonds,
                             double temperature);
void resolveCollisions(std::vector<Atom>& atoms);
void handleReactions(std::vector<Atom>& atoms,
                     std::vector<Bond>& bonds);

bool bondExists(const std::vector<Bond>& bonds, int a, int b);

void spawnWater(std::vector<Atom>& atoms,
                std::vector<Bond>& bonds,
                std::vector<AngleConstraint>& angles);

void spawnCarbonDioxide(std::vector<Atom>& atoms,
                        std::vector<Bond>& bonds,
                        std::vector<AngleConstraint>& angles);

void spawnMethane(std::vector<Atom>& atoms,
                  std::vector<Bond>& bonds,
                  std::vector<AngleConstraint>& angles);


void applyAngleForces(std::vector<Atom>& atoms,
                      const std::vector<AngleConstraint>& angles);

std::string detectMolecules(const std::vector<Atom>& atoms,
                            const std::vector<Bond>& bonds);

void stabilizeBondLengths(std::vector<Atom>& atoms,
                          const std::vector<Bond>& bonds);

double calculateKineticEnergy(const std::vector<Atom>& atoms);

double calculateBondPotentialEnergy(const std::vector<Atom>& atoms,
                                    const std::vector<Bond>& bonds);

double calculateTotalEnergy(const std::vector<Atom>& atoms,
                            const std::vector<Bond>& bonds);




