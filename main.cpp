#include <SFML/Graphics.hpp>
#include <vector>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include "Atom.h"
#include "Simulation.h"
#include "Renderer.h"
#include <string>

int main()
{
    std::srand(static_cast<unsigned int>(std::time(nullptr)));

    std::vector<Atom> atoms;
    std::vector<Bond> bonds;
    std::vector<AngleConstraint> angles;
    std::string detectedMolecules = "";
    std::string simulationMode = "Free Simulation";


    double temperature = 0.0;
    bool paused = false;
    bool usingPreset = false;
    double kineticEnergy = 0.0;
    double bondEnergy = 0.0;
    double totalEnergy = 0.0;

    detectedMolecules = detectMolecules(atoms, bonds);


    std::string moleculeLabel = "Random Atom Mixture";
    createAtoms(atoms, 20);

    sf::RenderWindow window(sf::VideoMode({800, 600}), "Atom Simulation");
    window.setFramerateLimit(60);

    sf::Font font;

    if (!font.openFromFile("assets/Arial.ttf"))
    {
    std::cout << "Failed to load font." << std::endl;
    return 1;
    }

    while (window.isOpen())
{
    while (const std::optional event = window.pollEvent())
    {
        if (event->is<sf::Event::Closed>())
        {
            window.close();
        }

        if (const auto* keyPressed = event->getIf<sf::Event::KeyPressed>())
        {
            if (keyPressed->code == sf::Keyboard::Key::Escape)
            {
                window.close();
            }
            else if (keyPressed->code == sf::Keyboard::Key::Space)
            {
                paused = !paused;
            }
            else if (keyPressed->code == sf::Keyboard::Key::Up)
            {
                temperature += 0.2;
            }
            else if (keyPressed->code == sf::Keyboard::Key::Down)
            {
                temperature -= 0.2;

                if (temperature < 0.0)
                {
                    temperature = 0.0;
                }
            }
            else if (keyPressed->code == sf::Keyboard::Key::R)
            {
                atoms.clear();
                bonds.clear();
                createAtoms(atoms, 20);
            }
            else if (keyPressed->code == sf::Keyboard::Key::Num1)
            {
            spawnWater(atoms, bonds, angles);
            usingPreset = true;
            moleculeLabel = "H2O - Water";
            usingPreset = true;
            simulationMode = "Preset Molecule";
            }
            else if (keyPressed->code == sf::Keyboard::Key::Num2)
            {
            spawnCarbonDioxide(atoms, bonds, angles);
            usingPreset = true;
            moleculeLabel = "CO2 - Carbon Dioxide";
            usingPreset = true;
            simulationMode = "Preset Molecule";
            }
            else if (keyPressed->code == sf::Keyboard::Key::Num3)
            {
            spawnMethane(atoms, bonds, angles);
            usingPreset = true;
            moleculeLabel = "CH4 - Methane";
            usingPreset = true;
            simulationMode = "Preset Molecule";
            }
            else if (keyPressed->code == sf::Keyboard::Key::R)
            {
                atoms.clear();
                bonds.clear();
                angles.clear();

                createAtoms(atoms, 20);

                usingPreset = false;
                moleculeLabel = "Random Atom Mixture";
                
                usingPreset = false;
                simulationMode = "Free Simulation";
            }   
        }
    }

    std::string detectedMolecules = "";

    if (!paused)
    {
        applyForces(atoms);
        applyBondForces(atoms, bonds);
        //applyAngleForces(atoms, angles);
        applyTemperature(atoms, temperature);
    
        updateAtoms(atoms);

        for (int i = 0; i < 8; i++)
        {
        stabilizeBondLengths(atoms, bonds);
        resolveCollisions(atoms);
        handleReactions(atoms, bonds);
        }
        
        breakOverstretchedBonds(atoms, bonds, temperature);

        if (!usingPreset)
        {
        breakOverstretchedBonds(atoms, bonds, temperature);
        updateBonds(atoms, bonds);
        }

        kineticEnergy = calculateKineticEnergy(atoms);
        bondEnergy = calculateBondPotentialEnergy(atoms, bonds);
        totalEnergy = calculateTotalEnergy(atoms, bonds);

        detectedMolecules = detectMolecules(atoms, bonds);
    }

    window.clear(sf::Color::Black);

    drawBonds(window, atoms, bonds);
    drawAtoms(window, atoms, font);

    drawUI(window,
       temperature,
       paused,
       atoms.size(),
       bonds.size(),
       moleculeLabel,
       detectedMolecules,
       simulationMode,
       kineticEnergy,
       bondEnergy,
       totalEnergy,
       font);
       

    window.display();
}


    return 0;
}

//--- g++ -std=c++17 main.cpp Simulation.cpp Renderer.cpp -o main \
-I/opt/homebrew/include \
-L/opt/homebrew/lib \
-lsfml-graphics \
-lsfml-window \
-lsfml-system ---//

