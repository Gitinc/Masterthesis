#ifndef PULSE_H
#define PULSE_H
#include "cmath"
#include "complex"
#include <random>

class Pulse
{
public:
    Pulse(double m_discretePoints, double m_height, double m_startPoint, double m_endPoint, double m_amplitudeNoiseStrength);
    std::vector<double> PulseVec;
    double height;
    void generateGaussPulse();
    void NoiseOnPulse();
    std::vector<double> addPulseWithPhaseshift(std::vector<double> Pulse_2, double phaseShift, double phaseNoiseStrength);

protected:

private:
    double discretePoints;
    double startPoint;
    double endPoint;
    double amplitudeNoiseStrength;
    double phaseNoiseStrength;
};

#endif // PULSE_H
