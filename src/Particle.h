/*
 * Particle.h
 *
 *  Created on: 23.02.2010
 *      Author: eckhardw
 */

#pragma once

#include <array>
#include <string>
#include <vector>

class Particle {

private:
    /**
     * Position of the particle
     */
    std::array<double, 3> x;

    /**
     * Velocity of the particle
     */
    std::array<double, 3> v;

    /**
     * Force effective on this particle
     */
    std::array<double, 3> f;

    /**
     * Force which was effective on this particle
     */
    std::array<double, 3> old_f;

    /**
     * Mass of this particle
     */
    double m;

    /**
     * Type of the particle. Use it for whatever you want (e.g. to separate
     * molecules belonging to different bodies, matters, and so on)
     */
    int type;

    double gGrav;

    double sigma;

    double epsilon;

    int xIndex;

    int yIndex;

public:
    explicit Particle(int type = 0);

    Particle(const Particle &other);

    Particle(
            std::array<double, 3> x_arg, std::array<double, 3> v_arg, double m_arg, double gGrav_type,
            double sigma_arg, double epsilon_arg, int type_arg);

    Particle(
            std::array<double, 3> x_arg, std::array<double, 3> v_arg, double m_arg, double gGrav_type,
            double sigma_arg, double epsilon_arg, int type_arg, int xIndex_arg, int yIndex_arg);

    Particle(std::array<double, 3> x_arg, std::array<double, 3> v_arg, std::array<double,
            3> f_arg, std::array<double, 3> oldf_arg, double m_arg, double gGrav_arg, double sigma_arg,
             double epsilon_arg, int type);

    Particle(std::array<double, 3> x_arg, std::array<double, 3> v_arg, std::array<double,
            3> f_arg, std::array<double, 3> oldf_arg, double m_arg, double gGrav_arg, double sigma_arg,
             double epsilon_arg, int type, int xIndex_arg, int yIndex_arg);

    virtual ~Particle();

    const std::array<double, 3> &getX() const;

    const std::array<double, 3> &getV() const;

    const std::array<double, 3> &getF() const;

    const std::array<double, 3> &getOldF() const;

    double getM() const;

    int getType() const;

    double getGGrav() const;

    double getSigma() const;

    double getEps() const;

    int getXIndex() const;

    int getYIndex() const;

    bool operator==(const Particle &other) const;

    std::string toString() const;

    void setX(const std::array<double, 3> &x_new);

    void setV(const std::array<double, 3> &v_new);

    void addF(const std::array<double, 3> &f_new);

    void resetF();
};

std::ostream &operator<<(std::ostream &stream, Particle &p);
