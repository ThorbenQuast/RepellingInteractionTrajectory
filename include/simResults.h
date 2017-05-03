#ifndef SIMRESULTS_H
#define SIMRESULTS_H


#include <vector>
#include <map>
#include <string>
#include <fstream>
#include <iostream>


class simResult {
  public:
    simResult();
    void addEntry(std::string, double value);
    void printHeader(std::ofstream&);
    void printEntries(std::ofstream&);

    void set_energy(double energy) { _energy = energy; };
    void set_radialDistance(double radialDistance) { _radialDistance = radialDistance; };
    void set_coincidenceTime(double coincidenceTime) { _coincidenceTime = coincidenceTime; };
    void set_angleCut(double angleCut) { _angleCut = angleCut; };
    void set_electronBeamAngle(double electronBeamAngle) { _electronBeamAngle = electronBeamAngle; };
    void set_maxIteration(int maxIteration) { _maxIteration = maxIteration; };
    void set_T_max(double T_max) { _T_max = T_max; };
    void set_repulsion(bool repulsion) { _repulsion = repulsion; };

  private:
    std::map<std::string, std::vector<double> > _data;

    double _energy;
    double _radialDistance;
    double _coincidenceTime;
    double _electronBeamAngle;
    double _angleCut;
    int _maxIteration;
    double _T_max;
    bool _repulsion;

};



#endif
