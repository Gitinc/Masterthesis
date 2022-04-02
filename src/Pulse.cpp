#include "Pulse.h"
#include "cmath"
#include "complex"
#include <random>
#include <time.h>

Pulse::Pulse(double m_discretePoints, double m_height, double m_startPoint, double m_endPoint, double m_amplitudeNoiseStrength):
        discretePoints(m_discretePoints), height(m_height), startPoint(m_startPoint), endPoint(m_endPoint), amplitudeNoiseStrength(m_amplitudeNoiseStrength), PulseVec(m_discretePoints, 0){};

void Pulse::generateGaussPulse(){
    const std::complex<double> I(0.0,1.0);
    double sigma = 1/(sqrt(M_PI_2)*height);
    double mu = 0;
    for(int i = 0; i < discretePoints; i++){
        double t = (startPoint + abs(startPoint-endPoint)/discretePoints*i);
        PulseVec[i] = height * exp(-(t-mu)*(t-mu)/(2*sigma*sigma));
    }
}

void Pulse::NoiseOnPulse(){
    std::default_random_engine generator_amplitude(time(NULL));
    for(int i = 0; i < discretePoints; i++){
        std::normal_distribution<double> dist_ampli(0,amplitudeNoiseStrength*PulseVec[i]);
        PulseVec[i] += dist_ampli(generator_amplitude);
    }
}

std::vector<double> Pulse::addPulseWithPhaseshift(std::vector<double> Pulse_2, double phaseShift, double phaseNoiseStrength){
    std::vector<double> PulseAbsAdded(discretePoints);
    std::default_random_engine generator_phase(time(NULL));
    const std::complex<double> I(0.0,1.0);
    for(int i = 0; i < discretePoints;i++){
        std::normal_distribution<double> dist_phase(0,phaseNoiseStrength);
        PulseAbsAdded[i] = abs(PulseVec[i]+Pulse_2[i]*exp(I*(phaseShift+dist_phase(generator_phase))));
        if(PulseAbsAdded[i]<0.0000000001){
            PulseAbsAdded[i] = 0;
        }
    }

    return PulseAbsAdded;
}
