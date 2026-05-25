#include "Renderer.h"
#include <array>

sf::Color getAtomColor(const std::string& name)
{
    if (name == "H") return sf::Color::White;
    if (name == "O") return sf::Color::Red;
    if (name == "C") return sf::Color(80, 80, 80);
    if (name == "Cl") return sf::Color::Green;

    return sf::Color::Magenta;
}


sf::Vector2f atomToScreen(const Atom& atom)
{
    const float windowWidth = 800.0f;
    const float windowHeight = 600.0f;

    float screenX =
        static_cast<float>((atom.x + 10.0) / 20.0 * windowWidth);

    float screenY =
        static_cast<float>((atom.y + 10.0) / 20.0 * windowHeight);

    return {screenX, screenY};
}

void drawAtoms(sf::RenderWindow& window,
               const std::vector<Atom>& atoms,
               const sf::Font& font)
{
   

    for (const Atom& atom : atoms)
    {
        sf::Vector2f screenPos = atomToScreen(atom);

        float renderRadius =
            static_cast<float>(atom.radius * 10.0);

        sf::CircleShape circle(renderRadius);

        circle.setOrigin(
        {
            renderRadius,
            renderRadius
        });

        circle.setPosition(screenPos);

        // atom colors
        if (atom.name == "H")
        {
            circle.setFillColor(sf::Color::White);
        }
        else if (atom.name == "O")
        {
            circle.setFillColor(sf::Color::Red);
        }
        else if (atom.name == "C")
        {
            circle.setFillColor(sf::Color(100, 100, 100));
        }
        else if (atom.name == "Cl")
        {
            circle.setFillColor(sf::Color::Green);
        }

        window.draw(circle);

        // LABELS
        sf::Text label(font, atom.name, 11);

        label.setFillColor(sf::Color::White);
        label.setOutlineColor(sf::Color::Black);
        label.setOutlineThickness(1.0f);

        label.setPosition(
        {
            screenPos.x - 7.0f,
            screenPos.y - 10.0f
        });

        window.draw(label);
    }
}


void drawBonds(sf::RenderWindow& window,
               const std::vector<Atom>& atoms,
               const std::vector<Bond>& bonds)
{
    for (const Bond& bond : bonds)
    {
        sf::Vector2f posA = atomToScreen(atoms[bond.atomA]);
        sf::Vector2f posB = atomToScreen(atoms[bond.atomB]);

        std::array line =
        {
            sf::Vertex{posA, sf::Color::Cyan},
            sf::Vertex{posB, sf::Color::Cyan}
        };

        window.draw(line.data(), line.size(), sf::PrimitiveType::Lines);
    }
}

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
            const sf::Font& font)
{
    std::string text =
        "Mode: " + simulationMode + "\n" +
        "Molecule: " + moleculeLabel + "\n" +
        "Detected: " + detectedMolecules + "\n" +
        "Temperature: " + std::to_string(temperature) + "\n" +
        "Kinetic Energy: " + std::to_string(kineticEnergy) + "\n" +
        "Bond Energy: " + std::to_string(bondEnergy) + "\n" +
        "Total Energy: " + std::to_string(totalEnergy) + "\n" +
        "Atoms: " + std::to_string(atomCount) + "\n" +
        "Bonds: " + std::to_string(bondCount) + "\n" +
        "Space: Pause | R: Reset | Up/Down: Temp\n"
        "1: H2O | 2: CO2 | 3: CH4";

    if (paused)
    {
        text += "\nPAUSED";
    }

    sf::Text uiText(font, text, 15);
    uiText.setFillColor(sf::Color::White);
    uiText.setPosition({10.0f, 10.0f});

    window.draw(uiText);
}


