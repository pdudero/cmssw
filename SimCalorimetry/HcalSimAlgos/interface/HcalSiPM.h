// -*- C++ -*-
#ifndef HcalSimAlgos_HcalSiPM_h
#define HcalSimAlgos_HcalSiPM_h

/**

  \class HcalSiPM

  \brief A general implementation for the response of a SiPM.

*/

#include <vector>
#include <algorithm>

#include "CLHEP/Random/RandGaussQ.h"
#include "CLHEP/Random/RandPoissonQ.h"
#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/RandBinomial.h"

class HcalSiPM {
 public:
  HcalSiPM(int nCells = 1, double tau = 15.);

  virtual ~HcalSiPM();

  void resetSiPM() { std::fill(theSiPM.begin(), theSiPM.end(), -999.); }
  //virtual int hitCells(unsigned int photons, unsigned int integral = 0) const;
  virtual double hitCells(unsigned int pes, double tempDiff = 0., 
			  double photonTime = 0.);


  virtual double totalCharge() const { return totalCharge(theLastHitTime); }
  virtual double totalCharge(double time) const;
  // virtual void recoverForTime(double time, double dt = 0.);

  int    getNCells()      const { return theCellCount; }
  double getTau()         const { return 1.0/theTauInv; }
  double getCrossTalk()   const { return theCrossTalk; }
  double getTempDep()     const { return theTempDep; }
  double getDarkCurrent() const { return darkCurrent_uA; }

  void setNCells(int nCells);
  void setTau(double tau) {theTauInv=1.0/tau;}
  void setCrossTalk(double xtalk);
  void setTemperatureDependence(double tempDep);
  void setDarkCurrent(double dc_uA) {darkCurrent_uA = dc_uA;}

  void initRandomEngine(CLHEP::HepRandomEngine& engine);


 protected:

  // void expRecover(double dt);

  double cellCharge(double deltaTime) const;
  double addCrossTalkCells_old(double in_pes);
  double addCrossTalkCells_binom(double in_pes);
  double addCrossTalkCells_borel(double in_pes);

  //numerical random generation from Borel-Tanner distribution
  double Borel(int n, double lambda, int k);
  std::vector<double> BorelCDF(double lambda, int k, int cdfsize=100);
  int    BorelRandomN(const std::vector<double>& cdf);

  unsigned int theCellCount;
  std::vector< double > theSiPM;
  double theTauInv;
  double theCrossTalk;
  double theTempDep;
  double theLastHitTime;
  double darkCurrent_uA;

  mutable CLHEP::RandGaussQ *theRndGauss;
  mutable CLHEP::RandPoissonQ *theRndPoisson;
  mutable CLHEP::RandFlat *theRndFlat;
  mutable CLHEP::RandBinomial *theRndBinom;

};

#endif //HcalSimAlgos_HcalSiPM_h
