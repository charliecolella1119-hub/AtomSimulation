#pragma once
#include <SFML/Graphics.hpp>
#include <vector>
#include "Atom.h"

sf::Color getAtomColor(const std::string& name);

sf::Vector2f atomToScreen(const Atom& atom);

void drawAtoms(sf::RenderWindow& window,
               const std::vector<Atom>& atoms,
               const sf::Font& font);

void drawBonds(sf::RenderWindow& window,
               const std::vector<Atom>& atoms,
               const std::vector<Bond>& bonds);


void drawUI(sf::RenderWindow& window,
            double temperature,
            bool paused,
            int atomCount,
            int bondCount,
            const std::string& moleculeLabel,
            const std::string& detectedMolecules,
            const std::string& simulationMode,
            double kineticEnergy,
            double bondEnergy,
            double totalEnergy,
            const sf::Font& font);

