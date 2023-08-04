#ifndef SiddhartaCard_h
#define SiddhartaCard_h 1

#include <G4SystemOfUnits.hh>
#include <G4ThreeVector.hh>
#include <globals.hh>

#include <vector>
#include <map>

#include <iostream>
#include <fstream>

typedef CLHEP::Hep3Vector G4ThreeVector;

class SiddhartaCard
{

public:
  SiddhartaCard();
  ~SiddhartaCard();

  static SiddhartaCard* getInstance();
  int GetG4RandGaussSeed();
  void SetG4RandGaussSeed(int);
  std::map<G4String,G4double> variables;

  void ReadCard();
  void Initialization();

  void SetCardFileName(std::string filename) {cardName = filename;};

  void SetInteractiveTrue() {interactive = true;};
  void SetDoneStatus(double status) {done = status;};
  void SetNoDataStatus(bool status) {nodata = status;};
  void SetEventMult(double status) {eventmult = status;};

  void SetLumiBoostPosition(G4ThreeVector middle, G4ThreeVector sizes);
  void SetLumiAntiBoostPosition(G4ThreeVector middle, G4ThreeVector sizes);

  void AddDoneByOne() {done++;};
  
  std::string GetCardFileName() {
    std::string result = cardName;
    std::string dot = ".";
    std::size_t dotPlace = result.rfind(dot);
    result.replace(dotPlace, result.length() - dotPlace, "_");
    return result;
  };

  bool GetInteractive() {return interactive;};
  bool GetNoDataStatus() {return nodata;};
  int GetDoneStatus() {return done;};
  int GetEventMult() {return eventmult;};

  std::vector<double> GetCustomDensities() {return customDensitiesStep;};
  std::vector<G4ThreeVector> GetLumiBoostPosition() {return lumiBoostPlane;};
  std::vector<G4ThreeVector> GetLumiAntiBoostPosition() {return lumiAntiBoostPlane;};

  void OpenBeam(char *name){return inbeam.open(name);};
  void CloseBeam(){return inbeam.close();}
  bool BeamOpen(){return inbeam.is_open();}
  bool BeamEof(){return inbeam.eof();}
  double xx, xxp, yy, yyp, zz, de, rate, turn;

  void ResetCounter(){counter=0;};

  int GetBeam() {inbeam >> xx >> xxp >> yy >> yyp >> zz >> de >> rate >> turn; counter++; return counter;};


private:
  static SiddhartaCard* fManager;
  int G4RandGaussSeed;
  
  std::ifstream inbeam;
  int counter;
  std::string cardName = "CARD.DAT";

  bool interactive = false;
  bool nodata = true;
  int done = 0;
  int eventmult;

  std::vector<double> customDensitiesStep;
// 4 vertices of the plane of lumi
  std::vector<G4ThreeVector> lumiBoostPlane;
  std::vector<G4ThreeVector> lumiAntiBoostPlane;
};

#endif
