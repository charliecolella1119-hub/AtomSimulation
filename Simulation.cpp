#include "Simulation.h"
#include <cstdlib>
#include <ctime>
#include <cmath>
#include "Atom.h"
#include <map>
#include <set>
#include <sstream>

void createAtoms(std::vector<Atom>& atoms, int count)
{
    for (int i = 0; i < count; i++)
    {
        Atom atom;

        int type = rand() % 4;

        if (type == 0)
        {
            atom.name = "H";
            atom.mass = 1.008;
            atom.radius = 0.35;
            atom.maxBonds = 1;
            atom.currentBonds = 0;
        }
        else if (type == 1)
        {
            atom.name = "C";
            atom.mass = 12.011;
            atom.radius = 0.70;
            atom.maxBonds = 4;
            atom.currentBonds = 0;
        }
        else if (type == 2)
        {
            atom.name = "O";
            atom.mass = 15.999;
            atom.radius = 0.60;
            atom.maxBonds = 2;
            atom.currentBonds = 0;
        }
        else
        {
            atom.name = "Cl";
            atom.mass = 35.45;
            atom.radius = 1.00;
            atom.maxBonds = 1;
            atom.currentBonds = 0;
        }

        atom.x = (rand() % 200 - 100) / 10.0;
        atom.y = (rand() % 200 - 100) / 10.0;

        atom.vx = ((rand() % 20 - 10) / 100.0) / atom.mass;
        atom.vy = ((rand() % 20 - 10) / 100.0) / atom.mass;

        atoms.push_back(atom);
    }
}

void applyForces(std::vector<Atom>& atoms)
{
    const double attractionStrength = 0.0005;
    const double repulsionStrength = 0.02;
    

    for (Atom& a : atoms)
    {
        for (const Atom& b : atoms)
        {
        if (&a == &b)
            {
            continue;
         }

        double minimumDistance = a.radius + b.radius;

        double dx = b.x - a.x;
        double dy = b.y - a.y;

        double distance = std::sqrt(dx * dx + dy * dy);

            if (distance < 0.01)
            {
                continue;
            }

            double force = 0.0;

            if (distance < minimumDistance)
            {
                force = -repulsionStrength / distance;
            }
            else
            {
                force = attractionStrength;
            }

            double directionX = dx / distance;
            double directionY = dy / distance;

            a.vx += (force * directionX) / a.mass;
            a.vy += (force * directionY) / a.mass;
        }
    }
}

double getIdealBondLength(const Atom& a, const Atom& b)
{
    std::string pair = a.name + "-" + b.name;

    if (pair == "H-H") return 0.75;
    if (pair == "O-H" || pair == "H-O") return 1.0;
    if (pair == "C-H" || pair == "H-C") return 1.1;
    if (pair == "C-O" || pair == "O-C") return 1.25;
    if (pair == "C-Cl" || pair == "Cl-C") return 1.7;
    if (pair == "O-O") return 1.2;

    return 1.4;
}

void applyBondForces(std::vector<Atom>& atoms,
                     const std::vector<Bond>& bonds)
{
    const double springStrength = 0.002;
    const double dampingStrength = 0.12;

    for (const Bond& bond : bonds)
    {
        Atom& a = atoms[bond.atomA];
        Atom& b = atoms[bond.atomB];

        double idealBondLength = getIdealBondLength(a, b);

        double dx = b.x - a.x;
        double dy = b.y - a.y;

        double distance = std::sqrt(dx * dx + dy * dy);

        if (distance < 0.001)
        {
            continue;
        }

        double directionX = dx / distance;
        double directionY = dy / distance;

        double displacement = distance - idealBondLength;

        double springForce = displacement * springStrength;

        double relativeVx = b.vx - a.vx;
        double relativeVy = b.vy - a.vy;

        double relativeVelocityAlongBond =
            relativeVx * directionX + relativeVy * directionY;

        double dampingForce =
            relativeVelocityAlongBond * dampingStrength;

        double totalForce = springForce + dampingForce;

        a.vx += totalForce * directionX / a.mass;
        a.vy += totalForce * directionY / a.mass;

        b.vx -= totalForce * directionX / b.mass;
        b.vy -= totalForce * directionY / b.mass;
    }
}

void applyAngleForces(std::vector<Atom>& atoms,
                      const std::vector<AngleConstraint>& angles)
{
    for (const AngleConstraint& angle : angles)
    {
        Atom& a = atoms[angle.atomA];
        Atom& b = atoms[angle.atomB];
        Atom& c = atoms[angle.atomC];

        double abx = a.x - b.x;
        double aby = a.y - b.y;
        double cbx = c.x - b.x;
        double cby = c.y - b.y;

        double angleA = std::atan2(aby, abx);
        double angleC = std::atan2(cby, cbx);

        double currentAngle = angleC - angleA;

        while (currentAngle > M_PI) currentAngle -= 2 * M_PI;
        while (currentAngle < -M_PI) currentAngle += 2 * M_PI;

        double idealRadians = angle.idealAngleDegrees * M_PI / 180.0;
        double error = currentAngle - idealRadians;

        a.vx += std::sin(angleA) * error * angle.strength;
        a.vy -= std::cos(angleA) * error * angle.strength;

        c.vx -= std::sin(angleC) * error * angle.strength;
        c.vy += std::cos(angleC) * error * angle.strength;
    }
}

void applyTemperature(std::vector<Atom>& atoms, double temperature)
{
    for (Atom& atom : atoms)
    {
        double randomX = (rand() % 200 - 100) / 10000.0;
        double randomY = (rand() % 200 - 100) / 10000.0;

        atom.vx += randomX * temperature / atom.mass;
        atom.vy += randomY * temperature / atom.mass;
    }
}

void updateAtoms(std::vector<Atom>& atoms)
{
    for (Atom& atom : atoms)
    {
        atom.x += atom.vx;
        atom.y += atom.vy;

        atom.vx *= 0.995;
        atom.vy *= 0.995;

        if (atom.x > 10)
        {
            atom.x = 10;
            atom.vx *= -1;
        }
        else if (atom.x < -10)
        {
            atom.x = -10;
            atom.vx *= -1;
        }

        if (atom.y > 10)
        {
            atom.y = 10;
            atom.vy *= -1;
        }
        else if (atom.y < -10)
        {
            atom.y = -10;
            atom.vy *= -1;
        }
    }
}

bool canBond(const Atom& a, const Atom& b)
{
    std::string pair = a.name + "-" + b.name;

    return pair == "H-H" ||
           pair == "H-O" || pair == "O-H" ||
           pair == "C-H" || pair == "H-C" ||
           pair == "C-O" || pair == "O-C" ||
           pair == "C-Cl" || pair == "Cl-C" ||
           pair == "O-O";
}

void updateBonds(std::vector<Atom>& atoms, std::vector<Bond>& bonds)
{
    const double bondDistance = 1.5;

    for (Atom& atom : atoms)
    {
        atom.currentBonds = 0;
    }

    for (const Bond& bond : bonds)
    {
        atoms[bond.atomA].currentBonds++;
        atoms[bond.atomB].currentBonds++;
    }

    for (int i = 0; i < atoms.size(); i++)
    {
        for (int j = i + 1; j < atoms.size(); j++)
        {
            if (!canBond(atoms[i], atoms[j]))
            {
                continue;
            }

            if (bondExists(bonds, i, j))
            {
                continue;
            }

            if (atoms[i].currentBonds >= atoms[i].maxBonds ||
                atoms[j].currentBonds >= atoms[j].maxBonds)
            {
                continue;
            }

            double dx = atoms[j].x - atoms[i].x;
            double dy = atoms[j].y - atoms[i].y;
            double distance = std::sqrt(dx * dx + dy * dy);

            if (distance < bondDistance)
            {
                bonds.push_back({i, j});
                atoms[i].currentBonds++;
                atoms[j].currentBonds++;
            }
        }
    }
}

bool bondExists(const std::vector<Bond>& bonds, int a, int b)
{
    for (const Bond& bond : bonds)
    {
        if ((bond.atomA == a && bond.atomB == b) ||
            (bond.atomA == b && bond.atomB == a))
        {
            return true;
        }
    }

    return false;
}

Atom makeAtom(std::string name, double x, double y)
{
    Atom atom;

    atom.name = name;
    atom.x = x;
    atom.y = y;
    atom.vx = 0.0;
    atom.vy = 0.0;
    atom.currentBonds = 0;

    if (name == "H")
    {
        atom.mass = 1.008;
        atom.radius = 0.35;
        atom.maxBonds = 1;
    }
    else if (name == "C")
    {
        atom.mass = 12.011;
        atom.radius = 0.70;
        atom.maxBonds = 4;
    }
    else if (name == "O")
    {
        atom.mass = 15.999;
        atom.radius = 0.60;
        atom.maxBonds = 2;
    }
    else
    {
        atom.mass = 35.45;
        atom.radius = 1.00;
        atom.maxBonds = 1;
    }

    return atom;
}

void spawnWater(std::vector<Atom>& atoms,
                std::vector<Bond>& bonds,
                std::vector<AngleConstraint>& angles)
{
    atoms.clear();
    bonds.clear();
    angles.clear();

    atoms.push_back(makeAtom("O", 0.0, 0.0));
    atoms.push_back(makeAtom("H", -1.0, 0.7));
    atoms.push_back(makeAtom("H", 1.0, 0.7));

    bonds.push_back({0, 1});
    bonds.push_back({0, 2});

    angles.push_back({1, 0, 2, 104.5, 0.002});
}


void spawnCarbonDioxide(std::vector<Atom>& atoms,
                        std::vector<Bond>& bonds,
                        std::vector<AngleConstraint>& angles)
{
    atoms.clear();
    bonds.clear();
    angles.clear();

    atoms.push_back(makeAtom("C", 0.0, 0.0));
    atoms.push_back(makeAtom("O", -1.4, 0.0));
    atoms.push_back(makeAtom("O", 1.4, 0.0));

    bonds.push_back({0, 1});
    bonds.push_back({0, 2});

    angles.push_back({1, 0, 2, 180.0, 0.002});
}

void spawnMethane(std::vector<Atom>& atoms,
                  std::vector<Bond>& bonds,
                  std::vector<AngleConstraint>& angles)
{
    atoms.clear();
    bonds.clear();
    angles.clear();

    atoms.push_back(makeAtom("C", 0.0, 0.0));
    atoms.push_back(makeAtom("H", 1.0, 0.0));
    atoms.push_back(makeAtom("H", -1.0, 0.0));
    atoms.push_back(makeAtom("H", 0.0, 1.0));
    atoms.push_back(makeAtom("H", 0.0, -1.0));

    bonds.push_back({0, 1});
    bonds.push_back({0, 2});
    bonds.push_back({0, 3});
    bonds.push_back({0, 4});

    angles.push_back({1, 0, 3, 109.5, 0.001});
    angles.push_back({1, 0, 4, 109.5, 0.001});
    angles.push_back({2, 0, 3, 109.5, 0.001});
    angles.push_back({2, 0, 4, 109.5, 0.001});
}

void breakOverstretchedBonds(std::vector<Atom>& atoms,
                             std::vector<Bond>& bonds,
                             double temperature)
{
    const double maxBondLength = 2.4;
    const double heatWeakening = 0.05;

    double breakDistance = maxBondLength - heatWeakening * temperature;

    if (breakDistance < 1.2)
    {
        breakDistance = 1.2;
    }

    for (int i = bonds.size() - 1; i >= 0; i--)
    {
        Atom& a = atoms[bonds[i].atomA];
        Atom& b = atoms[bonds[i].atomB];

        double dx = b.x - a.x;
        double dy = b.y - a.y;

        double distance = std::sqrt(dx * dx + dy * dy);

        if (distance > breakDistance)
        {
            bonds.erase(bonds.begin() + i);
        }
    }
}

void stabilizeBondLengths(std::vector<Atom>& atoms,
                          const std::vector<Bond>& bonds)
{
    const double correctionStrength = 0.08;

    for (const Bond& bond : bonds)
    {
        Atom& a = atoms[bond.atomA];
        Atom& b = atoms[bond.atomB];

        double idealLength = getIdealBondLength(a, b);

        double dx = b.x - a.x;
        double dy = b.y - a.y;

        double distance = std::sqrt(dx * dx + dy * dy);

        if (distance < 0.001)
        {
            continue;
        }

        double difference = distance - idealLength;

        double directionX = dx / distance;
        double directionY = dy / distance;

        double correction = difference * correctionStrength;

        double totalMass = a.mass + b.mass;

        double aShare = b.mass / totalMass;
        double bShare = a.mass / totalMass;

        a.x += directionX * correction * aShare;
        a.y += directionY * correction * aShare;

        b.x -= directionX * correction * bShare;
        b.y -= directionY * correction * bShare;
    }
}


void resolveCollisions(std::vector<Atom>& atoms)
{
    for (int i = 0; i < atoms.size(); i++)
    {
        for (int j = i + 1; j < atoms.size(); j++)
        {
            Atom& a = atoms[i];
            Atom& b = atoms[j];

            double dx = b.x - a.x;
            double dy = b.y - a.y;
            double distance = std::sqrt(dx * dx + dy * dy);

            double minDistance = a.radius + b.radius;

            if (distance > 0.001 && distance < minDistance)
            {
                double overlap = minDistance - distance;

                double directionX = dx / distance;
                double directionY = dy / distance;

                a.x -= directionX * overlap * 0.5;
                a.y -= directionY * overlap * 0.5;

                b.x += directionX * overlap * 0.5;
                b.y += directionY * overlap * 0.5;

                a.vx *= 0.9;
                a.vy *= 0.9;
                b.vx *= 0.9;
                b.vy *= 0.9;
            }
        }
    }
}

void visitMolecule(int atomIndex,
                   const std::vector<Bond>& bonds,
                   std::set<int>& visited,
                   std::vector<int>& moleculeAtoms)
{
    visited.insert(atomIndex);
    moleculeAtoms.push_back(atomIndex);

    for (const Bond& bond : bonds)
    {
        int neighbor = -1;

        if (bond.atomA == atomIndex)
        {
            neighbor = bond.atomB;
        }
        else if (bond.atomB == atomIndex)
        {
            neighbor = bond.atomA;
        }

        if (neighbor != -1 && visited.count(neighbor) == 0)
        {
            visitMolecule(neighbor, bonds, visited, moleculeAtoms);
        }
    }
}

std::string detectMolecules(const std::vector<Atom>& atoms,
                            const std::vector<Bond>& bonds)
{
    std::set<int> visited;
    std::vector<std::string> formulas;

    for (int i = 0; i < atoms.size(); i++)
    {
        if (visited.count(i) == 0)
        {
            std::vector<int> moleculeAtoms;

            visitMolecule(i, bonds, visited, moleculeAtoms);

            // Skip isolated atoms
            if (moleculeAtoms.size() < 2)
            {
                continue;
            }

            std::map<std::string, int> counts;

            for (int index : moleculeAtoms)
            {
                counts[atoms[index].name]++;
            }

            std::string formula = "";

            if (counts["C"] > 0)
            {
                formula += "C";
                if (counts["C"] > 1) formula += std::to_string(counts["C"]);
            }

            if (counts["H"] > 0)
            {
                formula += "H";
                if (counts["H"] > 1) formula += std::to_string(counts["H"]);
            }

            if (counts["O"] > 0)
            {
                formula += "O";
                if (counts["O"] > 1) formula += std::to_string(counts["O"]);
            }

            if (counts["Cl"] > 0)
            {
                formula += "Cl";
                if (counts["Cl"] > 1) formula += std::to_string(counts["Cl"]);
            }

            formulas.push_back(formula);
        }
    }

    if (formulas.empty())
    {
        return "None";
    }

    std::ostringstream output;

    const int maxDisplayed = 5;

    for (int i = 0; i < formulas.size() && i < maxDisplayed; i++)
    {
        output << formulas[i];

        if (i < formulas.size() - 1 && i < maxDisplayed - 1)
        {
            output << ", ";
        }
    }

    if (formulas.size() > maxDisplayed)
    {
        output << " ...";
    }

    return output.str();
}


double calculateKineticEnergy(const std::vector<Atom>& atoms)
{
    double totalKE = 0.0;

    for (const Atom& atom : atoms)
    {
        double speedSquared =
            atom.vx * atom.vx +
            atom.vy * atom.vy;

        totalKE += 0.5 * atom.mass * speedSquared;
    }

    return totalKE;
}

double calculateBondPotentialEnergy(const std::vector<Atom>& atoms,
                                    const std::vector<Bond>& bonds)
{
    const double springConstant = 0.002;
    double totalPE = 0.0;

    for (const Bond& bond : bonds)
    {
        const Atom& a = atoms[bond.atomA];
        const Atom& b = atoms[bond.atomB];

        double idealLength = getIdealBondLength(a, b);

        double dx = b.x - a.x;
        double dy = b.y - a.y;

        double distance = std::sqrt(dx * dx + dy * dy);

        double displacement = distance - idealLength;

        totalPE += 0.5 * springConstant * displacement * displacement;
    }

    return totalPE;
}

double calculateTotalEnergy(const std::vector<Atom>& atoms,
                            const std::vector<Bond>& bonds)
{
    return calculateKineticEnergy(atoms)
         + calculateBondPotentialEnergy(atoms, bonds);
}

int countBondsForAtom(const std::vector<Bond>& bonds, int atomIndex)
{
    int count = 0;

    for (const Bond& bond : bonds)
    {
        if (bond.atomA == atomIndex || bond.atomB == atomIndex)
        {
            count++;
        }
    }

    return count;
}

void handleReactions(std::vector<Atom>& atoms, std::vector<Bond>& bonds)
{
    for (int i = 0; i < atoms.size(); i++)
    {
        for (int j = i + 1; j < atoms.size(); j++)
        {
            Atom& a = atoms[i];
            Atom& b = atoms[j];

            if (!canBond(a, b))
            {
                continue;
            }

            if (bondExists(bonds, i, j))
            {
                continue;
            }

            int bondsA = countBondsForAtom(bonds, i);
            int bondsB = countBondsForAtom(bonds, j);

            if (bondsA >= a.maxBonds || bondsB >= b.maxBonds)
            {
                continue;
            }

            double dx = b.x - a.x;
            double dy = b.y - a.y;
            double distance = std::sqrt(dx * dx + dy * dy);

            if (distance > 1.0)
            {
                continue;
            }

            double relativeVX = b.vx - a.vx;
            double relativeVY = b.vy - a.vy;

            double relativeSpeed =
                std::sqrt(relativeVX * relativeVX +
                          relativeVY * relativeVY);

            if (relativeSpeed < 0.5)
            {
                bonds.push_back({i, j});
            }
        }
    }
}

