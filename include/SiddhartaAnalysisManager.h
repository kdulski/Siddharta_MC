#ifndef SiddhartaAnalysisManager_h
#define SiddhartaAnalysisManager_h 1

//---------------------------------------------------------------------------
//
// ClassName:   SiddhartaAnalysisManager
//
// Description: Singleton class to hold analysis parameters and build histograms.
//              User cannot access to the constructor.
//              The pointer of the only existing object can be got via
//              SiddhartaAnalysisManager::GetInstance() static method.
//              The first invokation of this static method makes
//              the singleton object.
//
//----------------------------------------------------------------------------
//

#include "../include/SiddhartaEnergyDeposition.h"
#include "../include/SiddhartaCard.h"

#include <G4SystemOfUnits.hh>
#include <G4Event.hh>
#include <globals.hh>

#include <vector>

//?? -> maybe just attach the file?
class SiddhartaHisto;

class SiddhartaAnalysisManager
{

public:
  static SiddhartaAnalysisManager* getInstance();
  static void dispose();

  void bookHisto();

  void BeginOfRun();
  void EndOfRun();

  void BeginOfEvent();
  void EndOfEvent(const G4Event*);

  void AddParticle(G4double, G4double, G4double, G4double);
  void AddIsotope(G4double, G4double, G4double);
  void AddEnergy(G4double, G4double, G4double);

  void SetVerbose(G4int val) {verbose = val;};
  G4int GetVerbose() const {return verbose;};

  void SetFirstEventToDebug(G4int val) {nEvt1 = val;};
  G4int FirstEventToDebug() const {return nEvt1;};
  void SetLastEventToDebug(G4int val) {nEvt2 = val;};
  G4int LastEventToDebug() const {return nEvt2;};

  void SetMaxEnergyforHisto(G4double val) {histEMax = val;};
  G4double  GetMaxEnergyforHisto() const {return histEMax;};
  void SetMinEnergyforHisto(G4double val) {histEMin = val;};
  G4double  GetMinEnergyforHisto() const {return histEMin;};
  void SetNumBinforHisto(G4int val) {histNBin = val;};
  G4int  GeNumBinforHisto() const {return histNBin;};

  void SetThresholdEnergyforTarget(G4double val) {targetThresE = val;};
  G4double GetThresholdEnergyforTarget () const {return targetThresE;};
  void SetThresholdEnergyforDetector(G4double val) {detectorThresE = val;};
  G4double GetThresholdEnergyforDetector () const {return detectorThresE;};
  void SetPulseWidth(G4double val) {pulseWidth = val;};
  G4double GetPulseWidth () const {return pulseWidth;};

  bool CheckKaonIniMom(G4ThreeVector iniPos, G4ThreeVector momentum);

  SiddhartaHisto*  histo;
  G4bool YesHistos;
  G4double SDD_pos[384][3];//   G4double SDD_pos[144][3];
  G4double SDD_angle[384];//    G4double SDD_angle[144];

private:
//?? why private
  SiddhartaAnalysisManager();
  ~SiddhartaAnalysisManager();

  static SiddhartaAnalysisManager* fManager;
  SiddhartaCard* mycard;

  G4int verbose;
  G4int nEvt1;
  G4int nEvt2;

  G4double histEMax;
  G4double histEMin;
  G4int histNBin;

  G4double targetThresE;
  G4double detectorThresE;
  G4double pulseWidth;

  G4int iniKaonMomCond;

  std::vector <SiddhartaEnergyDeposition> Edepo;
};

#endif
