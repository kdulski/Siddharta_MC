#include "../include/SiddhartaCard.h"

#include <Randomize.hh>

#include <sstream>
#include <fstream>
#include <map>

SiddhartaCard* SiddhartaCard::fManager = nullptr;

SiddhartaCard::SiddhartaCard() {}

SiddhartaCard::~SiddhartaCard() {}

SiddhartaCard* SiddhartaCard::getInstance()
{
  if(!fManager) {
    fManager = new SiddhartaCard();
  }
  return fManager;
}

int SiddhartaCard::GetG4RandGaussSeed()
{
  return G4RandGaussSeed;
}

void SiddhartaCard::SetG4RandGaussSeed(int iseed)
{
  G4RandGaussSeed = iseed;
  G4RandGauss::setTheSeed(iseed);
}

void SiddhartaCard::Initialization()
{
  variables["SiddhartaSetupVersion"] = 1;
  variables["z_kmtop"] = 60.0*mm;
  variables["z_ghost"] = 130.0*mm;
  variables["dx_kmtop"] = 80.0*mm;
  variables["dy_kmtop"] = 80.0*mm;
  variables["z_kplus"] = -90.0*mm;
  variables["StoreTrajectory"] = 0;
  variables["lowlimit"] = 1000.; // in eV
  variables["defaultCutValue"] = 100.; // in um
  variables["OldTarget"] = 1;
  variables["Beam"] = 0;
  variables["externalBeam"] = 0;
  variables["backgroundType"] = 0;
  variables["nbOfSeconds"] = 1;
  variables["beamCurrent"] = 1000;
  variables["useDegFromCard"] = 0;
  variables["degFromCard"] = 550; // in um
  variables["useGasFromCard"] = 0;
  variables["atomicNumberOfTarget"] = 20;
  variables["targetDensityAsLiquidFraction"] = 0.3; // percent of liquid density
  variables["useGasGradient"] = 0;
  variables["useGasGradientCustomSteps"] = 0;
  variables["densityGradientNumberOfSteps"] = 1;
  variables["densityGradientSizeOfStep"] = 0; // in percent of liquid density
  variables["densityAddCustomGradientDensity"] = 0; // in percent of liquid density
  variables["kaonBottomBoostShift"] = 0; // in cm
  variables["kaonBottomVerticalShift"] = 0; // in cm
  variables["activeDetMode"]= 1; // 1 = SIDDHARTA, 2 = CZT, 3 = HPGE, 4 = All, 5 = SIDDHARTA&&CZT, 6 = SIDDHARTA&&HPGE, 7 = CZT&HPGe
  variables["initialKaonMomentumCheck"] = 0; // -1 Anti Boost, 0 No check, 1 Boost
  customDensitiesStep.clear();
}

void SiddhartaCard::ReadCard()
{
  G4String st;
  G4double db;

  G4String lineTemp;
  G4String firstChar;

  std::ifstream indata;
  indata.open(cardName.c_str());

  if(!indata.is_open()) {
    G4cout << "Problem opening file CARD.DAT" << G4endl;
  } else {
    G4cout << "Reading variables from CARD.DAT" << G4endl;
    std::map<G4String,G4double>::iterator it;

  //  indata >> st >> db;
  //  while(!indata.eof()) {
    while(std::getline(indata, lineTemp)) {
      std::istringstream is(lineTemp);
      firstChar = lineTemp[0];
      if (firstChar != "/" && firstChar != "\n" && firstChar != " " && !lineTemp.empty()) {
        is >> st >> db;
        if (st == "densityAddCustomGradientDensity") {
          std::cout << st << " " << db << std::endl;
          customDensitiesStep.push_back(db);
        } else if(variables.find(st) != variables.end()) {
          variables[st] = db;
          G4cout << st << " " << variables[st] << G4endl;
        } else {
          G4cout << st << " not defined" << G4endl;
        }
   //   indata >> st >> db;
      }
    }
  }
  indata.close();
}

void SiddhartaCard::SetLumiBoostPosition(G4ThreeVector middle, G4ThreeVector sizes)
{
//Ignoring Y axis
  G4ThreeVector temp;
  temp.setX(middle.x() - sizes.x());
  temp.setY(middle.y() - sizes.y());
  temp.setZ(middle.z() - sizes.z());
  lumiBoostPlane.push_back(temp);
  temp.setX(middle.x() - sizes.x());
  temp.setY(middle.y() + sizes.y());
  temp.setZ(middle.z() + sizes.z());
  lumiBoostPlane.push_back(temp);
}

void SiddhartaCard::SetLumiAntiBoostPosition(G4ThreeVector middle, G4ThreeVector sizes)
{
//Ignoring Y axis
  G4ThreeVector temp;
  temp.setX(middle.x() + sizes.x());
  temp.setY(middle.y() - sizes.y());
  temp.setZ(middle.z() - sizes.z());
  lumiAntiBoostPlane.push_back(temp);
  temp.setX(middle.x() + sizes.x());
  temp.setY(middle.y() + sizes.y());
  temp.setZ(middle.z() + sizes.z());
  lumiAntiBoostPlane.push_back(temp);
}

