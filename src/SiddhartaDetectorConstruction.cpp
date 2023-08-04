#include "../include/SiddhartaLumiDetectorAntiBoostSD.h"
#include "../include/SiddhartaKaonDetectorBottomSD.h"
#include "../include/SiddhartaDetectorConstruction.h"
#include "../include/SiddhartaLumiDetectorBoostSD.h"
#include "../include/SiddhartaKaonDetectorTopSD.h"
#include "../include/SiddhartaDetectorMessenger.h"
#include "../include/SiddhartaAnalysisManager.h"
#include "../include/SiddhartaMagneticField.h"
#include "../include/SiddhartaScintAntiSD.h"
#include "../include/SiddhartaTrackerSD.h"
#include "../include/SiddhartaKPlusSD.h"
#include "../include/SiddhartaGhostSD.h"
#include "../include/SiddhartaAnti1SD.h"
#include "../include/SiddhartaAnti2SD.h"
#include "../include/SiddhartaAnti3SD.h"
#include "../include/SiddhartaAnti4SD.h"
#include "../include/SiddhartaAnti5SD.h"
#include "../include/SiddhartaHisto.h"
#include "../include/SiddhartaCard.h"

#include "../include/KLIMAXBoxWindowAntiBoostSD.h"
#include "../include/KLIMAXTarget2AntiBoostSD.h"
#include "../include/KLIMAXTarget1AntiBoostSD.h"
#include "../include/KLIMAXDegAntiBoostSD.h"
#include "../include/KLIMAXCZTAntiBoostSD.h"
#include "../include/KLIMAXTargetBoostSD.h"
#include "../include/KLIMAXDegBoostSD.h"
#include "../include/KLIMAXHPGeBoostSD.h"


#include <G4PhysicalConstants.hh>
#include <G4IntersectionSolid.hh>
#include <G4GeometryTolerance.hh>
#include <G4SubtractionSolid.hh>
#include <G4PVParameterised.hh>
#include <G4GeometryManager.hh>
#include <G4RotationMatrix.hh>
#include <G4EllipticalTube.hh>
#include <G4VisAttributes.hh>
#include <G4LogicalVolume.hh>
#include <G4tgbVolumeMgr.hh>
#include <G4tgrMessenger.hh>
#include <G4PVPlacement.hh>
#include <G4UserLimits.hh>
#include <G4UnionSolid.hh>
#include <G4SDManager.hh>
#include <G4Polyhedra.hh>
#include <G4Material.hh>
#include <G4Colour.hh>
#include <G4Torus.hh>
#include <G4Tubs.hh>
#include <G4Trap.hh>
#include <G4Cons.hh>
#include <G4Trd.hh>
#include <G4Box.hh>
#include <G4ios.hh>

//Rewrite constructor to more clear form
SiddhartaDetectorConstruction::SiddhartaDetectorConstruction() : solidWorld(0),  logicWorld(0),  physiWorld(0),
                        solidTarget(0), logicTarget(0), physiTarget(0), solidBeamPipe(0), logicBeamPipe(0), physiBeamPipe(0),
                        solidKaonDetectorTop(0), logicKaonDetectorTop(0), physiKaonDetectorTop(0),
                        solidKaonDetectorBottom(0), logicKaonDetectorBottom(0), physiKaonDetectorBottom(0),
                        solidLumiDetectorBoost(0), logicLumiDetectorBoost(0), physiLumiDetectorBoost(0),
                        solidLumiDetectorAntiBoost(0), logicLumiDetectorAntiBoost(0), physiLumiDetectorAntiBoost(0),
                        solidKLIMAXDegAntiBoost(0),logicKLIMAXDegAntiBoost(0),physiKLIMAXDegAntiBoost(0),
                        solidKLIMAXTarget1AntiBoost(0),logicKLIMAXTarget1AntiBoost(0),physiKLIMAXTarget1AntiBoost(0),
                        solidKLIMAXTarget2AntiBoost(0),logicKLIMAXTarget2AntiBoost(0),physiKLIMAXTarget2AntiBoost(0),
                        solidKLIMAXBoxWindowAntiBoost(0),logicKLIMAXBoxWindowAntiBoost(0),physiKLIMAXBoxWindowAntiBoost(0),
                        solidKLIMAXCZTAntiBoost(0),logicKLIMAXCZTAntiBoost(0),physiKLIMAXCZTAntiBoost(0),
                        solidKLIMAXDegBoost(0),logicKLIMAXDegBoost(0),physiKLIMAXDegBoost(0),
                        solidKLIMAXTargetBoost(0),logicKLIMAXTargetBoost(0),physiKLIMAXTargetBoost(0),
                        solidKLIMAXHPGeBoost(0),logicKLIMAXHPGeBoost(0),physiKLIMAXHPGeBoost(0),
                        solidKPlusDetector(0), logicKPlusDetector(0), physiKPlusDetector(0),
                        solidDegrader(0), logicDegrader(0), physiDegrader(0), TargetMater(0), BeamPipeMater(0),
                        KaonDetectorTopMater(0), KaonDetectorBottomMater(0), KPlusDetectorMater(0), DegraderMater(0),
                        stepLimit(nullptr), fpMagField(0), fWorldLength(0.),  fTargetLength(0.)
{
  fpMagField = new SiddhartaMagneticField();
  detectorMessenger = new SiddhartaDetectorMessenger(this);
  fDegModif = new DegraderModificator();
}

SiddhartaDetectorConstruction::~SiddhartaDetectorConstruction()
{
  delete fpMagField;
  delete stepLimit;
  delete detectorMessenger;
}

G4VPhysicalVolume* SiddhartaDetectorConstruction::Construct()
{
  SiddhartaCard* mycard = SiddhartaCard::getInstance();
  int SiddhartaSetup = mycard->variables["SiddhartaSetupVersion"];
  int OldTarget = mycard->variables["OldTarget"];

  int degThickFromCardOption = mycard->variables["useDegFromCard"];
  G4double degThickFromCard = mycard->variables["degFromCard"]*um;

  G4cout << "---===--- Degrader read from CARD ---===---" << G4endl;
  G4cout << "---===--- Using card for degrader?: " << degThickFromCardOption << " ---===---" << G4endl;
  G4cout << "---===--- " << degThickFromCard << " [um] ---===---" << G4endl;
  fDegModif->SetThicknessFromHaLF(degThickFromCard);
  for (unsigned i=1; i<=fDegModif->GetLayerSize(); i++) {
    G4cout << "Layer nr " << i << " " << fDegModif->GetThickAtLayer(i) << G4endl;
  }
  G4cout << "---===--- ---===---" << G4endl;

  int targetMaterialFromCard = mycard->variables["useGasFromCard"];
  G4double atomicNumberOfTargetFromCard = mycard->variables["atomicNumberOfTarget"];
  G4double targetDensity = mycard->variables["targetDensityAsLiquidFraction"];

  G4cout << "+++===+++ Target read from CARD +++===+++" << G4endl;
  G4cout << "+++===+++ Using card for target?: " << targetMaterialFromCard << " +++===+++" << G4endl;
  G4cout << "+++===+++ Atomic number of target: " << atomicNumberOfTargetFromCard << " +++===+++" << G4endl;
  G4cout << "+++===+++ Liquid density fraction of target [%]: " << targetDensity << " +++===+++" << G4endl;
  G4cout << "+++===+++ +++   +++" << G4endl;

  int targetMaterialDensityGradient = mycard->variables["useGasGradient"];
  int targetMaterialCustomDensityGradient = mycard->variables["useGasGradientCustomSteps"];
  G4double densityGradientNumberOfSteps = mycard->variables["densityGradientNumberOfSteps"];
  G4double densityGradientSizeOfStep = mycard->variables["densityGradientSizeOfStep"];
  std::vector<double> customDensitiesContainer = mycard->GetCustomDensities();
  G4double liquidTargetDensityForPrinting = 1;

  if (targetMaterialCustomDensityGradient) {
    targetMaterialDensityGradient = 1; // To ensure that all of the options for density gradient holds. If checked properly can be removed
    densityGradientNumberOfSteps = customDensitiesContainer.size(); // Custom dominant over standard uniform distributed gradient
  }


  G4cout << "+++===+++ Using density gradient for target?: " << targetMaterialDensityGradient << " +++===+++" << G4endl;
  G4cout << "+++===+++ Using custom density gradient for target?: " << targetMaterialCustomDensityGradient << " +++===+++" << G4endl;
  G4cout << "+++===+++ Number of gradient steps: " << densityGradientNumberOfSteps << " +++===+++" << G4endl;
  G4cout << "+++===+++ Size of gradient step: " << densityGradientSizeOfStep << " +++===+++" << G4endl;
  G4cout << "+++===+++ Size of custom gradient step vector: " << customDensitiesContainer.size() << " +++===+++" << G4endl;
  G4cout << "+++===+++ +++===+++" << G4endl;

  SiddhartaAnalysisManager* analysis = SiddhartaAnalysisManager::getInstance();

//--------- Material definition ---------
  G4double a, z;
  G4double density;
  G4int nel;

  G4Element* H = new G4Element("Hydrogen", "H", z=1., a=1.00794*g/mole);
  G4Isotope* Deut = new G4Isotope("Deuteron", 1, 2, 2.0141018*g/mole);
  G4Element* D = new G4Element("Deuterium", "D", 1);
  D->AddIsotope(Deut, 1);
  G4Element* C = new G4Element("Carbon", "C", z=6., a=12.0107*g/mole);
  G4Element* N = new G4Element("Nitrogen", "N", z=7., a=14.01*g/mole);
  G4Element* O = new G4Element("Oxygen"  , "O", z=8., a=16.00*g/mole);
  G4Element* F = new G4Element("Fluor", "F", z=9., a=18.9984032*g/mole);
  G4Element* al = new G4Element("Aluminium","Al", z=13., a=26.981539*g/mole); //Why al not from capital?
  G4Element* Si = new G4Element("Silicon","Si", z=14., a=28.0855*g/mole);
  G4Element* Cl = new G4Element("Chlorine","Cl", z=17., a=35.453*g/mole);
  G4Element* Co = new G4Element("Cobalt","Co", z=27., a=58.933200*g/mole);
  G4Element* Zn = new G4Element("Zinc", "Zn", z=30., a= 65.38*g/mole);
  G4Element* Cd = new G4Element("Cadmium", "Cd", z=48., a= 112.411*g/mole);
  G4Element* Te = new G4Element("Tellurium", "Te", z=52., a= 127.60*g/mole);
  G4Element* Sm = new G4Element("Samarium","Sm", z=62., a=150.36*g/mole);
  G4Element* Ge = new G4Element("Germaniuìum","Ge", z=32., a=72.64*g/mole);

  G4Material* Al = new G4Material("Al", z=13., a=26.981539*g/mole, density=2.7*g/cm3);
  G4Material* Silicon = new G4Material("Silicon", z=14., a=28.0855*g/mole, density=2.3290*g/cm3);
  G4Material* Ti = new G4Material("Ti", z=22., a=47.867*g/mole, density=4.54*g/cm3);
  G4Material* Cu = new G4Material("Copper", z=29., a=63.546*g/mole, density=8.94*g/cm3);
  G4Material* Zr = new G4Material("Zr", z=40., a=91.224*g/mole, density=6.52*g/cm3);
  G4Material* Pb = new G4Material("Lead", z=82., a=207.19*g/mole, density=11.35*g/cm3);

  G4Material* Air = new G4Material("Air", density=1.29*mg/cm3, nel=2);
  Air->AddElement(N, 70*perCent);
  Air->AddElement(O, 30*perCent);

  density = universe_mean_density;
  G4double pressure = 3.e-18*pascal;
  G4double temperature = 2.73*kelvin;
  G4Material* vacuum = new G4Material("Galactic", z=1., a=1.01*g/mole, density, kStateGas, temperature, pressure);

  G4Material* Hydrogen;
  //G4double LHD = 0.0678*g/cm3;
  G4double LHD = 0.07*g/cm3;
  G4double LDD = 0.16*g/cm3; // 24.06.2021 Wikipedia value
  G4double LNeD = 1.207*g/cm3; // 16.02.2023 Wikipedia value

  if (SiddhartaSetup == 2) {
    Hydrogen = new G4Material("Hydrogen", z=1., a=1.00794*g/mole, density=3*0.01*LHD);
  } else if (SiddhartaSetup == 1) {
    Hydrogen = new G4Material("Hydrogen", z=1., a=1.00794*g/mole, density=1.5*0.01*LHD);
  } else if (SiddhartaSetup == 2020) {
    if (targetMaterialFromCard && atomicNumberOfTargetFromCard == 1)
      Hydrogen = new G4Material("Hydrogen", z=1., a=1.00794*g/mole, density=targetDensity*0.01*LHD);
    else
      Hydrogen = new G4Material("Hydrogen", z=1., a=1.00794*g/mole, density=10*0.01*LHD);
  } else {
    Hydrogen = new G4Material("Hydrogen", z=1., a=1.00794*g/mole, density=3*0.01*LHD);
  }

  G4Material* deuterium;
  if (targetMaterialFromCard && atomicNumberOfTargetFromCard == 2)
    deuterium = new G4Material("deuterium", density=targetDensity*0.01*LDD, 1);
  else
    deuterium = new G4Material("deuterium", density=1.04*0.01*LDD, 1);
  deuterium->AddElement(D, 1);

  G4Material* Helium4;
  G4double LHeD = 0.125*g/cm3;
  if (targetMaterialFromCard && atomicNumberOfTargetFromCard == 4)
    Helium4 = new G4Material("Helium4", z=2., a=4.00*g/mole, density=targetDensity*0.01*LHeD);
  else
    Helium4 = new G4Material("Helium4", z=2., a=4.00*g/mole, density=0.014*LHeD);

  G4Material* Neonium;
  if (targetMaterialFromCard && (atomicNumberOfTargetFromCard == 20 || atomicNumberOfTargetFromCard == 22))
    Neonium = new G4Material("Neonium", z=10., a=20.18*g/mole, density=targetDensity*0.01*LNeD);
  else
    Neonium = new G4Material("Neonium", z=10., a=20.18*g/mole, density=0.3*0.01*LNeD);

// Definition of the gradient material for target
  std::vector<G4Material*> gradMaterials;

  if (targetMaterialDensityGradient || targetMaterialCustomDensityGradient) {
    for (int i=0; i<densityGradientNumberOfSteps || i<customDensitiesContainer.size(); i++) {
      double tempDens = targetDensity - i*densityGradientSizeOfStep + (densityGradientNumberOfSteps - 1)*0.5*densityGradientSizeOfStep;

      if (targetMaterialCustomDensityGradient && i < customDensitiesContainer.size()) {
        tempDens = customDensitiesContainer.at(i);
        std::cout << "Cutom density step: " << i+1 << " with density " << tempDens << std::endl;
      } else if (targetMaterialCustomDensityGradient) { // Not enough steps defined in card/cfg
        std::cout << "!!POSSIBLE ERROR: More steps defined in config/card file than the custom densities given!!" << std::endl;
      }
      G4Material *tempMat;

      if (targetMaterialFromCard) {
        if (atomicNumberOfTargetFromCard == 1) {
          tempMat = new G4Material("HydrogenGrad" + std::to_string(i), z=1., a=1.00794*g/mole, density=tempDens*0.01*LHD);
          liquidTargetDensityForPrinting = LHD;
        } else if (atomicNumberOfTargetFromCard == 2) {
          tempMat = new G4Material("deuteriumGrad" + std::to_string(i), density=tempDens*0.01*LDD, 1);
          tempMat->AddElement(D, 1);
          liquidTargetDensityForPrinting = LDD;
        } else if (atomicNumberOfTargetFromCard == 4) {
          tempMat = new G4Material("Helium4Grad" + std::to_string(i), z=2., a=4.00*g/mole, density=tempDens*0.01*LHeD);
          liquidTargetDensityForPrinting = LHeD;
        } else if (atomicNumberOfTargetFromCard == 20 || atomicNumberOfTargetFromCard == 22) {
          tempMat = new G4Material("NeoniumGrad" + std::to_string(i), z=10., a=20.18*g/mole, density=tempDens*0.01*LNeD);
          liquidTargetDensityForPrinting = LNeD;
        } else {
          tempMat = new G4Material("HydrogenGrad" + std::to_string(i), z=1., a=1.00794*g/mole, density=tempDens*0.01*LHD);
          liquidTargetDensityForPrinting = LHD;
        }
      } else {
        tempMat = new G4Material("HydrogenGrad" + std::to_string(i), z=1., a=1.00794*g/mole, density=tempDens*0.01*LHD);
        liquidTargetDensityForPrinting = LHD;
      }
      gradMaterials.push_back(tempMat);
    }
  }
//

  G4Material* carbonFiber = new G4Material("fiber", density = 1.35*g/cm3, nel=2);  // estimation at Lab
  carbonFiber->AddElement(C, 7);
  carbonFiber->AddElement(C, 3);  // pure carbon presentation
//carbonFiber->AddElement(H, 3);

  G4Material* graphite = new G4Material("graphite", density = 2.2*g/cm3, nel=1);
  graphite->AddElement(C, 1);

  G4Material* Polyethylene = new G4Material("Polyethylene", density=0.930*g/cm3, nel=2); // DENSITYYYYY ////
  Polyethylene->AddElement(C, 2);
  Polyethylene->AddElement(H, 4);

  G4Material* SmCo = new G4Material("SmCo", density=8.4*g/cm3, nel=2);
  SmCo->AddElement(Sm, 1);
  SmCo->AddElement(Co, 5);

  G4Material* BC420 = new G4Material("BC420", density=1.032*g/cm3, nel=2);
  BC420->AddElement(C, 27);
  BC420->AddElement(H, 30);

  G4Material* teflon = new G4Material("teflon", density= 2.2*g/cm3, nel=2);
  teflon->AddElement(C, 1);
  teflon->AddElement(F, 4);

  G4Material* mylar = new G4Material("mylar", density= 1.370*g/cm3, nel=3); // DENSITY !!!!!! ////////////////
  mylar->AddElement(C, 10);
  mylar->AddElement(H, 8);
  mylar->AddElement(O, 4);

  G4Material* blacktape = new G4Material("blacktape", density= 1.370*g/cm3, nel=3); // DENSITY !!!!!! //////////////// EQUAL TO MYLAR ?? ////
  blacktape->AddElement(C, 10);
  blacktape->AddElement(H, 8);
  blacktape->AddElement(O, 4);

  G4Material* pvc = new G4Material("pvc", density = 1.370*g/cm3, nel=3);
  pvc->AddElement(C, 2);
  pvc->AddElement(H, 3);
  pvc->AddElement(Cl, 1);

  G4Material* polycarbonate = new G4Material("polycarbonate", density = 1.21*g/cm3, nel=3);
  polycarbonate->AddElement(C, 16);
  polycarbonate->AddElement(H, 14);
  polycarbonate->AddElement(O, 3);

  G4Material* paper = new G4Material("paper", density=1.370*g/cm3, nel=3); // DENSITYYYYY ////
  paper->AddElement(C, 10);
  paper->AddElement(H, 8);
  paper->AddElement(O, 4);

  G4Material* black_paper = new G4Material("black_paper", density=1.37*g/cm3, nel=3);
  black_paper->AddElement(C, 10);
  black_paper->AddElement(H, 8);
  black_paper->AddElement(O, 4);

  G4Material* Ceramic = new G4Material("Ceramic", density=2.7*g/cm3, nel=3);
  Ceramic->AddElement(al, 2);
  Ceramic->AddElement(Si, 1);
  Ceramic->AddElement(O, 3);

  G4Material* kapton = new G4Material("kapton", density=1.42*g/cm3, nel=4);
  kapton->AddElement(C, 22);
  kapton->AddElement(H, 10);
  kapton->AddElement(O, 5);
  kapton->AddElement(N, 2);

  G4Material* carbonPan = new G4Material("cpan", density = 1.84*g/cm3, nel=1);  // estimation at Lab
  carbonPan->AddElement(C, 1);
  
  G4Material* Epoxy = new G4Material("epox", density = 1.1*g/cm3, nel=3);  // estimation at Lab
  Epoxy->AddElement(C, 21);
  Epoxy->AddElement(H, 24);
  Epoxy->AddElement(O, 4);

  G4Material* CZT = new G4Material("CZT", density= 5.76*g/cm3,nel=3);
  CZT->AddElement(Cd, 45*perCent);
  CZT->AddElement(Zn, 5*perCent);
  CZT->AddElement(Te, 50*perCent);

  G4Material* HPGe = new G4Material("HPGe", density= 5.323*g/cm3,nel=1); // beware that density is at room temperature
  HPGe->AddElement(Ge, 1);
  
  G4MaterialTable matTable = *(G4Material::GetMaterialTable());
  analysis->histo->addMatIDFromString("G4_AIR");
  analysis->histo->addMatIDFromString("G4_Galactic");
  for (unsigned i=0; i<matTable.size(); i++) {
    analysis->histo->addMatIDFromString(matTable.at(i)->GetName());
  }

  G4cout << G4endl << "The materials defined are : " << G4endl << G4endl;
  G4cout << *(G4Material::GetMaterialTable()) << G4endl;

//--------- Rotation matrix  ---------
//There is alsow constructor for RotationMatrix with arguments like that: phi, theta, psi -> It could simplify the code
  G4double thetabp2 = 14.30*deg; //Used later
  G4RotationMatrix* rm00 = new G4RotationMatrix();
//Rotations by Pi/2
  G4double phi = 90.*deg;
  G4RotationMatrix* rm01 = new G4RotationMatrix();
  rm01->rotateX(phi);
  G4RotationMatrix* rmZ90 = new G4RotationMatrix();
  rmZ90->rotateZ(phi);
  G4RotationMatrix* rmY90 = new G4RotationMatrix();
  rmY90->rotateY(phi);
  G4RotationMatrix* rm20 = new G4RotationMatrix(); // Double definition of the same rotation
  rm20->rotateX(phi);
  G4RotationMatrix* rmXYZ90 = new G4RotationMatrix();
  rmXYZ90->rotateY(phi);
  rmXYZ90->rotateZ(phi);
  rmXYZ90->rotateX(phi);
//~Pi/2
//Rotations by Pi
  phi = 180.*deg;
  G4RotationMatrix* rm02 = new G4RotationMatrix();
  rm02->rotateX(phi);
  G4RotationMatrix* rm15 = new G4RotationMatrix();
  rm15->rotateX(phi);
//~Pi/2
//Rotations by predefined value 1
  G4double thetabp3 = 4.07*deg;
  G4RotationMatrix* rm16 = new G4RotationMatrix();
  rm16->rotateY(thetabp3);
  G4RotationMatrix* rm17 = new G4RotationMatrix();
  rm17->rotateY(-thetabp3);
//~1
//Rotations by -Pi
  phi = -90.*deg;
  G4RotationMatrix* rm18 = new G4RotationMatrix();
  rm18->rotateX(phi);
  G4RotationMatrix* rm22 = new G4RotationMatrix(); // Double definition
  rm22 -> rotateX(phi);
//~-Pi
//Rotation by -Pi/2 and Pi
  G4RotationMatrix* rm19 = new G4RotationMatrix();
  rm19->rotateX(phi);
  phi = 180.*deg;
  rm19->rotateY(phi);
//~-Pi/2 + Pi
//Rotation by Pi/2 + Pi
  G4RotationMatrix* rm21 = new G4RotationMatrix();
  phi = 90.*deg;
  rm21->rotateX(phi);
  rm21->rotateY(2*phi);
//~Pi/2 + Pi

  G4double anglecalor = twopi/12.0; // 30 degrees -> Pi/6
//Rotations by -Pi/2 + segmentation of 2*Pi onto 12 parts
  G4RotationMatrix* rm23 = new G4RotationMatrix();
  phi = -90.*deg;
  rm23->rotateX(phi);
  rm23->rotateY(anglecalor);
  G4RotationMatrix* rm24 = new G4RotationMatrix();
  rm24->rotateX(phi);
  rm24->rotateY(2.*anglecalor);
  G4RotationMatrix* rm25 = new G4RotationMatrix();
  rm25->rotateX(phi);
  rm25->rotateY(3.*anglecalor);
  G4RotationMatrix* rm26 = new G4RotationMatrix();
  rm26->rotateX(phi);
  rm26->rotateY(4.*anglecalor);
  G4RotationMatrix* rm27 = new G4RotationMatrix();
  rm27->rotateX(phi);
  rm27->rotateY(5.*anglecalor);
  G4RotationMatrix* rm28 = new G4RotationMatrix();
  rm28->rotateX(phi);
  rm28->rotateY(6.*anglecalor);
  G4RotationMatrix* rm29 = new G4RotationMatrix();
  rm29->rotateX(phi);
  rm29->rotateY(7.*anglecalor);
  G4RotationMatrix* rm30 = new G4RotationMatrix();
  rm30->rotateX(phi);
  rm30->rotateY(8.*anglecalor);
  G4RotationMatrix* rm31 = new G4RotationMatrix();
  rm31->rotateX(phi);
  rm31->rotateY(9.*anglecalor);
  G4RotationMatrix* rm32 = new G4RotationMatrix();
  rm32->rotateX(phi);
  rm32->rotateY(10.*anglecalor);
  G4RotationMatrix* rm33 = new G4RotationMatrix();
  rm33->rotateX(phi);
  rm33->rotateY(11.*anglecalor);
//~-Pi/2 + Pi/12-segmentation
//Rotations by Pi/2 + -Pi/2
  G4RotationMatrix* rm34 = new G4RotationMatrix();
  phi = 90.*deg;
  rm34->rotateX(phi);
  rm34->rotateY(-phi);
//~Pi/2 + -Pi/2
//Rotations by Pi/2 and Pi/2
  G4RotationMatrix* rm35 = new G4RotationMatrix();
  phi = 90.*deg;
  rm35->rotateX(phi);
  rm35->rotateY(phi);
  G4RotationMatrix* rm38 = new G4RotationMatrix();
  rm38 -> rotateX(phi);
  rm38 -> rotateZ(phi);
//~Pi/2 + Pi/2

//--------- Sizes of the principal geometrical components (solids)  ---------
    fTargetLength = 5.0*cm;                        // Full length of Target
    TargetMater = Pb;
    BeamPipeMater = Pb;
    KaonDetectorTopMater = Pb;
    KaonDetectorBottomMater = Pb;
    KPlusDetectorMater = Pb;
    DegraderMater = Pb;
    ShieldingMater = Pb;
    fWorldLength= 3.0*(fTargetLength + 560.*cm);

//--------- Definitions of Solids, Logical Volumes, Physical Volumes ---------

// ---------- Switches ------------------------
  G4bool USECALORIMETER;

  if (SiddhartaSetup == 1 ) {
    USECALORIMETER = 1;
  } else {
    USECALORIMETER = 0;
  }
  G4bool withMattoni = 0;
  G4bool withWall = 1;

//------------------------------------------------
// Define one or several text files containing the geometry description
//------------------------------------------------
  G4String filename = "geom.txt";
  G4tgbVolumeMgr* volmgr = G4tgbVolumeMgr::GetInstance();
  volmgr->AddTextFile(filename);

//------------------------------------------------
// Read the text files and construct the GEANT4 geometry
//------------------------------------------------
  G4VPhysicalVolume* physiWorld = volmgr->ReadAndConstructDetector();
  G4SDManager* SDman = G4SDManager::GetSDMpointer();

  G4String trackerSDD_SDname = "TextGeom/TrackerSDD_SD";
  SiddhartaTrackerSD* aTrackerSD = new SiddhartaTrackerSD(trackerSDD_SDname);
  SDman->AddNewDetector(aTrackerSD);
  G4LogicalVolume* sdd_Si = G4tgbVolumeMgr::GetInstance()->FindG4LogVol("SDD",0);

  G4String KaonDetectorTop_SDname = "TextGeom/KaonDetectorTop_SD";
  SiddhartaKaonDetectorTopSD* KDTop_SD = new SiddhartaKaonDetectorTopSD(KaonDetectorTop_SDname);
  SDman->AddNewDetector(KDTop_SD);
  G4LogicalVolume* sciTop = G4tgbVolumeMgr::GetInstance()->FindG4LogVol("KaonDetectorTop",0);

  G4String KaonDetectorBottom_SDname = "TextGeom/KaonDetectorBottom_SD";
  SiddhartaKaonDetectorBottomSD* KDBottom_SD = new SiddhartaKaonDetectorBottomSD(KaonDetectorBottom_SDname);
  SDman->AddNewDetector(KDBottom_SD);
  G4LogicalVolume* sciBottom = G4tgbVolumeMgr::GetInstance()->FindG4LogVol("KaonDetectorBottom",0);

  G4String LumiDetectorBoost_SDname = "TextGeom/LumiDetectorBoost_SD";
  SiddhartaLumiDetectorBoostSD* LDBoost_SD = new SiddhartaLumiDetectorBoostSD(LumiDetectorBoost_SDname);
  SDman->AddNewDetector(LDBoost_SD);
  G4LogicalVolume* sciBoost = G4tgbVolumeMgr::GetInstance()->FindG4LogVol("LumiDetectorBoost",0);

  G4String LumiDetectorAntiBoost_SDname = "TextGeom/LumiDetectorAntiBoost_SD";
  SiddhartaLumiDetectorAntiBoostSD* LDAntiBoost_SD = new SiddhartaLumiDetectorAntiBoostSD(LumiDetectorAntiBoost_SDname);
  SDman->AddNewDetector(LDAntiBoost_SD);
  G4LogicalVolume* sciAntiBoost = G4tgbVolumeMgr::GetInstance()->FindG4LogVol("LumiDetectorAntiBoost",0);

  G4String KLIMAXDegAntiBoost_SDname = "TextGeom/KLIMAXDegAntiBoost_SD";
  KLIMAXDegAntiBoostSD* KLDAntiBoost_SD = new KLIMAXDegAntiBoostSD( KLIMAXDegAntiBoost_SDname );
  SDman->AddNewDetector( KLDAntiBoost_SD );
  G4LogicalVolume* kldAntiBoost = G4tgbVolumeMgr::GetInstance()->FindG4LogVol("KLIMAXDegAntiBoost",0);

  G4String KLIMAXTarget1AntiBoost_SDname = "TextGeom/KLIMAXTarget1AntiBoost_SD";
  KLIMAXTarget1AntiBoostSD* KLT1AntiBoost_SD = new KLIMAXTarget1AntiBoostSD( KLIMAXTarget1AntiBoost_SDname );
  SDman->AddNewDetector( KLT1AntiBoost_SD );
  G4LogicalVolume* klt1AntiBoost = G4tgbVolumeMgr::GetInstance()->FindG4LogVol("KLIMAXTarget1AntiBoost",0);

  G4String KLIMAXTarget2AntiBoost_SDname = "TextGeom/KLIMAXTarget2AntiBoost_SD";
  KLIMAXTarget2AntiBoostSD* KLT2AntiBoost_SD = new KLIMAXTarget2AntiBoostSD( KLIMAXTarget2AntiBoost_SDname );
  SDman->AddNewDetector( KLT2AntiBoost_SD );
  G4LogicalVolume* klt2AntiBoost = G4tgbVolumeMgr::GetInstance()->FindG4LogVol("KLIMAXTarget2AntiBoost",0);

  G4String KLIMAXBoxWindowAntiBoost_SDname = "TextGeom/KLIMAXBoxWindowAntiBoost_SD";
  KLIMAXBoxWindowAntiBoostSD* KLBWAntiBoost_SD = new KLIMAXBoxWindowAntiBoostSD( KLIMAXBoxWindowAntiBoost_SDname );
  SDman->AddNewDetector( KLBWAntiBoost_SD );
  G4LogicalVolume* KLBWAntiBoost = G4tgbVolumeMgr::GetInstance()->FindG4LogVol("KLIMAXBoxWindowAntiBoost",0);

  G4String KLIMAXCZTAntiBoost_SDname = "TextGeom/KLIMAXCZTAntiBoost_SD";
  KLIMAXCZTAntiBoostSD* KLCZTAntiBoost_SD = new KLIMAXCZTAntiBoostSD( KLIMAXCZTAntiBoost_SDname );
  SDman->AddNewDetector( KLCZTAntiBoost_SD );
  G4LogicalVolume* klcztAntiBoost = G4tgbVolumeMgr::GetInstance()->FindG4LogVol("KLIMAXCZTAntiBoost",0);

  G4String KLIMAXDegBoost_SDname = "TextGeom/KLIMAXDegBoost_SD";
  KLIMAXDegBoostSD* KLDBoost_SD = new KLIMAXDegBoostSD( KLIMAXDegBoost_SDname );
  SDman->AddNewDetector( KLDBoost_SD );
  G4LogicalVolume* kldBoost = G4tgbVolumeMgr::GetInstance()->FindG4LogVol("KLIMAXDegBoost",0);

  G4String KLIMAXTargetBoost_SDname = "TextGeom/KLIMAXTargetBoost_SD";
  KLIMAXTargetBoostSD* KLTBoost_SD = new KLIMAXTargetBoostSD( KLIMAXTargetBoost_SDname );
  SDman->AddNewDetector( KLTBoost_SD );
  G4LogicalVolume* kltBoost = G4tgbVolumeMgr::GetInstance()->FindG4LogVol("KLIMAXTargetBoost",0);

  G4String KLIMAXHPGeBoost_SDname = "TextGeom/KLIMAXHPGeBoost_SD";
  KLIMAXHPGeBoostSD* KLHPGeBoost_SD = new KLIMAXHPGeBoostSD( KLIMAXHPGeBoost_SDname );
  SDman->AddNewDetector( KLHPGeBoost_SD );
  G4LogicalVolume* klHPGeBoost = G4tgbVolumeMgr::GetInstance()->FindG4LogVol("KLIMAXHPGeBoost",0);

  G4String ghost_SDname = "TextGeom/ghost_SD";
  SiddhartaGhostSD* ghostSD = new SiddhartaGhostSD(ghost_SDname);
  SDman->AddNewDetector(ghostSD);

  G4String KPlus_SDname = "TextGeom/KPlus_SD";
  SiddhartaKPlusSD* KPlus_SD = new SiddhartaKPlusSD(KPlus_SDname);
  SDman->AddNewDetector(KPlus_SD);
  G4LogicalVolume* kplusDet = G4tgbVolumeMgr::GetInstance()->FindG4LogVol("KPlus",0);

  G4String ScintAntiSD_SDname = "TextGeom/ScintAnti_SD";
  SiddhartaScintAntiSD* ScintAntiSD = new SiddhartaScintAntiSD(ScintAntiSD_SDname);
  SDman->AddNewDetector(ScintAntiSD);

  SiddhartaAnti1SD* Anti1_SD;
  SiddhartaAnti2SD* Anti2_SD;
  SiddhartaAnti3SD* Anti3_SD;
  SiddhartaAnti4SD* Anti4_SD;
  SiddhartaAnti5SD* Anti5_SD;

  G4String Anti1_SDname = "TextGeom/Anti1_SD";
  Anti1_SD = new SiddhartaAnti1SD(Anti1_SDname);
  SDman->AddNewDetector(Anti1_SD);
  G4String Anti2_SDname = "TextGeom/Anti2_SD";
  Anti2_SD = new SiddhartaAnti2SD(Anti2_SDname);
  SDman->AddNewDetector(Anti2_SD);
  G4String Anti3_SDname = "TextGeom/Anti3_SD";
  Anti3_SD = new SiddhartaAnti3SD(Anti3_SDname);
  SDman->AddNewDetector(Anti3_SD);
  G4String Anti4_SDname = "TextGeom/Anti4_SD";
  Anti4_SD = new SiddhartaAnti4SD(Anti4_SDname);
  SDman->AddNewDetector(Anti4_SD);
  G4String Anti5_SDname = "TextGeom/Anti5_SD";
  Anti5_SD = new SiddhartaAnti5SD(Anti5_SDname);
  SDman->AddNewDetector(Anti5_SD);

//--------- Magnetic Field Volume -------------------------------
  G4LogicalVolume* World = G4tgbVolumeMgr::GetInstance()->FindG4LogVol("World",0);
  World->SetVisAttributes(G4VisAttributes::GetInvisible());

//--------- Beam Pipe -------------------------------
  G4double bp_r3 = 0.5*59.0*mm;
  G4double bp_wd2 = 2.0*mm;

  G4double bp_wd1 = 350.*um;
  G4double bp_r1 = bp_r3-bp_wd2;
  G4double bp_r2 = bp_r1+bp_wd1;

///// 2020 configuration ///// ORIGINAL!!!

  G4double bp_fiber = 500.*um; //carbon fiber IP
  G4double bp_Al1 = 150.*um; //aluminum IP

  G4double bp_Al2 = 2.0*mm;  //aluminum other part thickness



// **** NEW value after Mihail's measurements of the BeamPipe in lab ***//
 /* G4double bp_fiber = 389.5*um;  // PAN carbon
  G4double bp_Al1 = 169.5*um; // aluminum IP
  G4double bp_epo = 141.7*um; // epoxy IP 
*/
  // ********************************************************************//



  G4double bp_r1_in = 0.5*58.3*mm;  // IP region
  G4double bp_r1_Al2Fiber = bp_r1_in + bp_Al1;    // interface between Carbon fiber and Al
  G4double bp_r1_out = bp_r1_in + bp_Al1 + bp_fiber;  //
  G4double bp_r2_in = 0.5*58.0*mm; //
  G4double bp_r2_out = bp_r2_in + bp_Al2; //
  G4double bp_L = 35.00*cm;          // IP region length
  G4double bp_L_p1 = (118 + 244) * mm;  // 2 mm thick straight part outside IP

// IP region, Aluminum in Fiber out - ORIGINAL!!!

  G4Tubs* solidBeamPipeIpFiber = new G4Tubs("BeamPipeIpFiber", bp_r1_Al2Fiber, bp_r1_out, bp_L, 0., 360.*deg);
  G4LogicalVolume* logicBeamPipeIpFiber = new G4LogicalVolume(solidBeamPipeIpFiber, carbonFiber, "World", 0, 0, 0);
  G4VPhysicalVolume* physiBeamPipe3 = new G4PVPlacement(0,                  // no rotation
            G4ThreeVector(0.,0.,0.), // at (0,0,0)
            logicBeamPipeIpFiber,      // its logical volume
            "BeamPipeIpFiber",         // its name
            World,               // its mother  volume
            false,           // no boolean operations
            0);              // copy number





// **** NEW part after Mihail's measurements of the BeamPipe in lab ***//

// New carbon Fiber beam pipe density 1.84g/cm3 instead of 1.35 g/cm3
/*  G4Tubs* solidBeamPipeIpFiber = new G4Tubs("BeamPipeIpFiber", bp_r1_Al2Fiber, bp_r1_out, bp_L, 0., 360.*deg);
  G4LogicalVolume* logicBeamPipeIpFiber = new G4LogicalVolume(solidBeamPipeIpFiber, carbonPan, "World", 0, 0, 0);
  G4VPhysicalVolume* physiBeamPipe3 = new G4PVPlacement(0,                  // no rotation
            G4ThreeVector(0.,0.,0.), // at (0,0,0)
            logicBeamPipeIpFiber,      // its logical volume
            "BeamPipeIpFiber",         // its name
            World,               // its mother  volume
            false,           // no boolean operations
            0);             // copy number
*/
//Epoxy beam pipe 
 /* G4double bp_r1_epo = bp_r1_out;
  G4double bp_r2_epo = bp_r1_out+bp_epo;
  G4Tubs* solidBeamPipeIpEpo = new G4Tubs("BeamPipeIpEpo", bp_r1_epo, bp_r2_epo, bp_L, 0., 360.*deg);
  G4LogicalVolume* logicBeamPipeIpEpo = new G4LogicalVolume(solidBeamPipeIpEpo, Epoxy, "World", 0, 0, 0);
  G4VPhysicalVolume* physiBeamPipe4 = new G4PVPlacement(0,                  // no rotation
            G4ThreeVector(0.,0.,0.), // at (0,0,0)
            logicBeamPipeIpEpo,      // its logical volume
            "BeamPipeIpEpo",         // its name
            World,               // its mother  volume
            false,           // no boolean operations
            0);              // copy number
*/
  // ********************************************************************//



  G4Tubs* solidBeamPipeIpAl = new G4Tubs("BeamPipeIpAl", bp_r1_in, bp_r1_Al2Fiber, bp_L, 0., 360.*deg);
  G4LogicalVolume* logicBeamPipeIpAl = new G4LogicalVolume(solidBeamPipeIpAl, Al, "World", 0, 0, 0);
  G4VPhysicalVolume* physiBeamPipeIpAl = new G4PVPlacement(0,               // no rotation
            G4ThreeVector(0.,0.,0.), // at (0,0,0)
            logicBeamPipeIpAl,      // its logical volume
            "BeamPipeIpAl",         // its name
            World,               // its mother  volume
            false,           // no boolean operations
            0);              // copy number

  G4Tubs* solidBeamPipeIpVacuum = new G4Tubs("BeamPipeIpVacuum", 0., bp_r1_in, bp_L, 0., 360.*deg);
  G4LogicalVolume* logicBeamPipeIpVacuum = new G4LogicalVolume( solidBeamPipeIpVacuum, vacuum, "World", 0, 0, 0);
  G4VPhysicalVolume* physiBeamPipeIpVacuum = new G4PVPlacement(0,         // no rotation
            G4ThreeVector(0.,0.,0.), // at (0,0,0)
            logicBeamPipeIpVacuum,      // its logical volume
            "BeamPipeIpVacuum",         // its name
            World,               // its mother  volume
            false,           // no boolean operations
            0);              // copy number

//// 2020 this part is obsolete?
////      have to define the Dafne luminometer
  G4double bpj_w;
  G4double bpj_L;
  bpj_w = 25*mm;
  bpj_L = 30*mm;
  G4Tubs* solidBeamPipeJunction;
  solidBeamPipeJunction = new G4Tubs("BeamPipeJunction", bp_r3, bp_r3+bpj_w, 0.5*bpj_L, 0., 360.*deg);
  G4LogicalVolume* logicBeamPipeJunction = new G4LogicalVolume( solidBeamPipeJunction, Al, "World", 0, 0, 0);
  if (SiddhartaSetup != 1) {
    G4VPhysicalVolume* physiBeamPipeJunction = new G4PVPlacement(0, G4ThreeVector(0., 0., 0.5*450.*mm+0.5*bpj_L),
                                                                 logicBeamPipeJunction, "BeamPipeJunction", World, false, 0);
    physiBeamPipeJunction = new G4PVPlacement(0, G4ThreeVector(0., 0., -0.5*450.*mm-0.5*bpj_L),
                                              logicBeamPipeJunction, "BeamPipeJunction", World, false, 1);
    logicBeamPipeJunction->SetVisAttributes(G4Colour(1.,0.,0.));
  }

//--------- Sputnik -------------------------------
  G4double sputnik_r1;
  G4double sputnik_r2;
  G4double sputnik_L;
  G4double sputnik_posz;

  if (SiddhartaSetup == 1) {
    sputnik_r1 = bp_r3 + 43.*mm;
    sputnik_r2 = bp_r3 + 13.*mm;
    sputnik_L = 164.*mm;
    sputnik_posz = 151.24*mm + 0.5*sputnik_L;

    const G4double zPlane[] = {0.5*sputnik_L, -0.5*sputnik_L};
    const G4double rInner[] = {0.5*109.*mm, 0.5*109.*mm};
    const G4double rOuter[] = {0.5*114.*mm, 0.5*200.*mm};

    const G4double zPlane1[] = {0.5*sputnik_L, -0.5*sputnik_L};
    const G4double rInner1[] = {0.5*109.*mm, 0.5*109.*mm};
    const G4double rOuter1[] = {0.5*114.*mm, 0.5*200.*mm};

    const G4double zPlane2[] = {0.5*sputnik_L-64.0*mm, 0.5*sputnik_L-64.0*mm-50.*mm};
    const G4double rInner2[] = {0.5*109.*mm, 0.5*109.*mm};
    const G4double rOuter2[] = {0.5*129.*mm, 0.5*129.*mm};

    const G4double zPlane3[] = {0.5*sputnik_L-64.0*mm-50.*mm, -0.5*sputnik_L};
    const G4double rInner3[] = {0.5*109.*mm, 0.5*109.*mm};
    const G4double rOuter3[] = {0.5*149.*mm, 0.5*149.*mm};

    G4Polyhedra* solidSputnik1 = new G4Polyhedra("Sputnik1", 0., 360.*deg, 8, 2, zPlane1, rInner1, rOuter1);
    G4Polyhedra* solidSputnik2 = new G4Polyhedra("Sputnik2", 0., 360.*deg, 8, 2, zPlane2, rInner2, rOuter2);
    G4Polyhedra* solidSputnik3 = new G4Polyhedra("Sputnik3", 0., 360.*deg, 8, 2, zPlane3, rInner3, rOuter3);

    G4Polyhedra* solidSputnik = new G4Polyhedra("Sputnik", 0., 360.*deg, 8, 2, zPlane, rInner, rOuter);
    G4LogicalVolume* logicSputnik = new G4LogicalVolume(solidSputnik, Pb, "World", 0, 0, 0);
    G4LogicalVolume* logicSputnik2 = new G4LogicalVolume(solidSputnik2, Air, "Sputnik", 0, 0, 0);
    G4LogicalVolume* logicSputnik3 = new G4LogicalVolume(solidSputnik3, Air, "Sputnik", 0, 0, 0);
    G4VPhysicalVolume* physiSputnik2 = new G4PVPlacement(0,
                G4ThreeVector(0.,0.,0.), // at (0,0,0)
                logicSputnik2,      // its logical volume
                "Sputnik2",         // its name
                logicSputnik,               // its mother  volume
                false,           // no boolean operations
                0);              // copy number
    G4VPhysicalVolume* physiSputnik3 = new G4PVPlacement(0,
                G4ThreeVector(0.,0.,0.), // at (0,0,0)
                logicSputnik3,      // its logical volume
                "Sputnik3",         // its name
                logicSputnik,               // its mother  volume
                false,           // no boolean operations
                0);              // copy number

    G4VPhysicalVolume* physiSputnik = new G4PVPlacement(0,
                G4ThreeVector(0.,0.,-sputnik_posz), // at (0,0,0)
                logicSputnik,      // its logical volume
                "Sputnik",         // its name
                World,               // its mother  volume
                false,           // no boolean operations
                0);              // copy number
    physiSputnik = new G4PVPlacement(rm15,
                G4ThreeVector(0.,0.,sputnik_posz), // at (0,0,0)
                logicSputnik,      // its logical volume
                "Sputnik",         // its name
                World,               // its mother  volume
                false,           // no boolean operations
                1);              // copy number

    G4Tubs* solidSputnikCyl = new G4Tubs("SputnikCyl", 0.5*60.0*mm, 0.5*108.0*mm, 0.5*64.0*mm,0., 360.*deg);
    G4LogicalVolume* logicSputnikCyl = new G4LogicalVolume(solidSputnikCyl,Pb, "World", 0, 0, 0);
    G4VPhysicalVolume* physiSputnikCyl = new G4PVPlacement(0, G4ThreeVector(0., 0., 151.24*mm + 0.5*64.0*mm), logicSputnikCyl,
                                                           "SputnikCyl", World, false, 1);
    physiSputnikCyl = new G4PVPlacement(0,G4ThreeVector(0., 0., -151.24*mm - 0.5*64.0*mm), logicSputnikCyl,
                                        "SputnikCyl", World, false, 2);
    logicSputnikCyl->SetVisAttributes(G4Colour(1.,0.70,0.));
    logicBeamPipeIpFiber->SetVisAttributes(G4Colour(0.0,1.0,0.0));
    logicBeamPipeIpAl->SetVisAttributes(G4Colour(1.0,0.0,0.0));
  } else {
    sputnik_r1 = bp_r3 + 43.*mm;
    sputnik_r2 = bp_r3 + 13.*mm;
    sputnik_L = 164.*mm;
    sputnik_posz = 151.24*mm + 0.5*sputnik_L;

    G4Cons* solidSputnik = new G4Cons("Sputnik", bp_r3, bp_r3 + 60.*mm, bp_r3, bp_r3 + 60.*mm, 0.5*50.*mm, 0., 360.*deg);
    G4LogicalVolume* logicSputnik = new G4LogicalVolume(solidSputnik, Pb, "World", 0, 0, 0);
  }

//---------Quadrupole 1 ---------------------//
  G4double quadrupole1_L;

  if (SiddhartaSetup != 1) {
    quadrupole1_L = 244.*mm;
  } else {
    quadrupole1_L = 15.0*cm;
  }
  G4double quadrupole1_W = 100.*mm-bp_r3;
  G4double pos_quad1;
  if (SiddhartaSetup != 1) {
    pos_quad1 = 0.5*450.*mm+68.*mm+0.5*quadrupole1_L;
  } else {
    pos_quad1 = -sputnik_posz-0.5*sputnik_L-0.5*quadrupole1_L;
  }
  G4Tubs* solidQuadrupole1 = new G4Tubs("Quadrupole1", bp_r3, bp_r3 + quadrupole1_W, 0.5*quadrupole1_L, 0., 360.*deg);
  G4LogicalVolume* logicQuadrupole1 = new G4LogicalVolume(solidQuadrupole1, SmCo, "World", 0, 0, 0);
  G4VPhysicalVolume* physiQuadrupole1 = new G4PVPlacement(rm00, G4ThreeVector(0.,0., pos_quad1), logicQuadrupole1,
                                                          "Quadrupole1", World, false, 0);
  physiQuadrupole1 = new G4PVPlacement(rm15, G4ThreeVector(0.,0.,-pos_quad1), logicQuadrupole1,
                                       "Quadrupole1", World, false, 1);
  logicQuadrupole1->SetVisAttributes(G4Colour(0.0,1.0,0.0));
//----------- Beam pipe junction -------------------
  G4double lj = 30.*cm;
  G4double rj = 1756.87*mm;//lj/thetabp2;
  G4Box* boxbpj = new G4Box("boxbpj", rj, rj, rj);
  G4Torus* torusbpj = new G4Torus("torusbpj", bp_r2, bp_r3, rj, 0, thetabp2);

  G4IntersectionSolid* interbpj = new G4IntersectionSolid("interbpj", boxbpj, torusbpj);
  G4LogicalVolume* logicinterbpj= new G4LogicalVolume(interbpj, Al, "World", 0, 0, 0);
  G4VPhysicalVolume* physiinterbpj = new G4PVPlacement(rm18,
            G4ThreeVector(-rj,0,bp_L),
            logicinterbpj,      // its logical volume
            "interbpj",         // its name
            World,               // its mother  volume
            false,           // no boolean operations
            0);              // copy number
  physiinterbpj = new G4PVPlacement(rm19,
            G4ThreeVector(rj,0,bp_L),
            logicinterbpj,      // its logical volume
            "interbpj",         // its name
            World,               // its mother  volume
            false,           // no boolean operations
            1);              // copy number
  physiinterbpj = new G4PVPlacement(rm20,
            G4ThreeVector(-rj,0,-bp_L),
            logicinterbpj,      // its logical volume
            "interbpj",         // its name
            World,               // its mother  volume
            false,           // no boolean operations
            2);              // copy number
  physiinterbpj = new G4PVPlacement(rm21,
            G4ThreeVector(rj,0,-bp_L),
            logicinterbpj,      // its logical volume
            "interbpj",         // its name
            World,               // its mother  volume
            false,           // no boolean operations
            3);              // copy number

  logicinterbpj->SetVisAttributes(G4Colour(0.0,0.0,1.0));

//----------- Beam pipe after junction -------------------
  G4double bp2_L = 150.*cm;
  G4double dx =  rj*(1 - cos(thetabp2)) + 0.5*bp2_L*sin(thetabp3);
  G4double dz =  rj*sin(thetabp2) + 0.5*bp2_L*cos(thetabp3);

  G4Tubs* solidBeamPipe21 = new G4Tubs("BeamPipe21", bp_r2, bp_r3, 0.5*bp2_L, 0., 360.*deg);
  G4LogicalVolume* logicBeamPipe21 = new G4LogicalVolume(solidBeamPipe21, Al, "World", 0, 0, 0);
  G4VPhysicalVolume* physiBeamPipe21 = new G4PVPlacement(rm16,
            G4ThreeVector(dx,0,-bp_L - dz),
            logicBeamPipe21,      // its logical volume
            "BeamPipe21",         // its name
            World,               // its mother  volume
            false,           // no boolean operations
            0);              // copy number

  physiBeamPipe21 = new G4PVPlacement(rm17,
            G4ThreeVector(-dx,0,-bp_L - dz),
            logicBeamPipe21,      // its logical volume
            "BeamPipe21",         // its name
            World,               // its mother  volume
            false,           // no boolean operations
            1);              // copy number

  physiBeamPipe21 = new G4PVPlacement(rm17,
            G4ThreeVector(dx,0,bp_L + dz),
            logicBeamPipe21,      // its logical volume
            "BeamPipe21",         // its name
            World,               // its mother  volume
            false,           // no boolean operations
            2);              // copy number
  physiBeamPipe21 = new G4PVPlacement(rm16,
            G4ThreeVector(-dx,0,bp_L + dz),
            logicBeamPipe21,      // its logical volume
            "BeamPipe21",         // its name
            World,               // its mother  volume
            false,           // no boolean operations
            3);              // copy number

  logicBeamPipe21 ->SetVisAttributes(G4Colour(1.0,0.0,0.0));

//---------Quadrupole 2 ---------------------//
  G4double quadrupole2_L = 25.0*cm;
  G4double quadrupole2_W = 1.5*cm;
  G4double dxq =  rj*(1 - cos(thetabp2)) + 0.5*quadrupole2_L*sin(thetabp3);
  G4double dzq =  rj*sin(thetabp2) + 0.5*quadrupole2_L*cos(thetabp3);

  G4Tubs* solidQuadrupole2 = new G4Tubs("Quadrupole2", bp_r3, bp_r3 + quadrupole2_W, 0.5*quadrupole2_L, 0., 360.*deg);
  G4LogicalVolume* logicQuadrupole2 = new G4LogicalVolume(solidQuadrupole2, SmCo, "World", 0, 0, 0);
  G4VPhysicalVolume* physiQuadrupole2 = new G4PVPlacement(rm16,
            G4ThreeVector(dxq,0,-bp_L - dzq),
            logicQuadrupole2,      // its logical volume
            "Quadrupole2",         // its name
            World,               // its mother  volume
            false,           // no boolean operations
            0);              // copy number

  physiQuadrupole2 = new G4PVPlacement(rm17,
            G4ThreeVector(-dxq,0,-bp_L - dzq),
            logicQuadrupole2,      // its logical volume
            "Quadrupole2",         // its name
            World,               // its mother  volume
            false,           // no boolean operations
            1);              // copy number

  physiQuadrupole2 = new G4PVPlacement(rm17,
            G4ThreeVector(dxq,0,bp_L + dzq),
            logicQuadrupole2,      // its logical volume
            "Quadrupole2",         // its name
            World,               // its mother  volume
            false,           // no boolean operations
            2);              // copy number
  physiQuadrupole2 = new G4PVPlacement(rm16,
            G4ThreeVector(-dxq,0,bp_L + dzq),
            logicQuadrupole2,      // its logical volume
            "Quadrupole2",         // its name
            World,               // its mother  volume
            false,           // no boolean operations
            3);              // copy number


  logicQuadrupole2->SetVisAttributes(G4Colour(0.0,1.0,0.0));

//------------- Calorimeter -----------------------
  if(USECALORIMETER) {
    G4double distQuadCal = 1.0*cm;
    G4double dx1 = 56.*mm;
    G4double dy1 = 1.0*cm;
    G4double dy2 = 1.0*cm;
    G4double dzL = 220.*mm;
    G4double calorR1 = 0.5*dx1/sin(0.5*anglecalor);
    G4double rcor1 = calorR1*cos(0.5*anglecalor);
    G4double rcor2 = rcor1 + dzL;
    G4double dx2 = 2.*rcor2*tan(0.5*anglecalor);
    G4double calorR2 = 0.5*dx2/sin(0.5*anglecalor);
    G4double rmean = 0.5*(rcor1 + rcor2);
    G4double calorLead5mm_W = 5.*mm;
    G4double calorLead10mm_W = 10.*mm;
    G4double calorLead5mm_D = 5.*mm;
    G4double calorLead10mm_D = 10.*mm;

    G4Trd* solidCalorLead5mm = new G4Trd("CalorimeterLead5mm", 0.5*dx1, 0.5*dx2, 0.5*calorLead5mm_W, 0.5*calorLead5mm_W, 0.5*dzL);
    G4LogicalVolume* logicCalorLead5mm = new G4LogicalVolume(solidCalorLead5mm, Pb, "World", 0, 0, 0);

    G4Trd* solidCalorLead10mm = new G4Trd("CalorimeterLead10mm", 0.5*dx1, 0.5*dx2, 0.5*calorLead10mm_W, 0.5*calorLead10mm_W, 0.5*dzL);
    G4LogicalVolume* logicCalorLead10mm = new G4LogicalVolume(solidCalorLead10mm, Pb, "World", 0, 0, 0);

    G4Trd* solidCalor = new G4Trd("Calorimeter", 0.5*dx1, 0.5*dx2, 0.5*dy1, 0.5*dy2, 0.5*dzL);
    G4LogicalVolume* logicCalor = new G4LogicalVolume(solidCalor, BC420, "World", 0, 0, 0); //CHANGE MATERIAALLLLLLL
    G4double calx1 = 0.;
    G4double caly1 = -rmean;
    G4double calzmean = (-sputnik_posz - 0.5*sputnik_L - 0.5*quadrupole1_L);
    G4double calz1 = calzmean;

    for (G4int j=0; j<12; j++) {
      G4double anglecalor = twopi/12.0;
      G4RotationMatrix* rmtmp = new G4RotationMatrix();
      phi = -90.*deg;
      rmtmp->rotateX(phi);
      rmtmp->rotateY(j*anglecalor);
      calx1 = -rmean*sin(j*anglecalor);
      caly1 = -rmean*cos(j*anglecalor);
      G4double caloffset = 0.0*cm;
      calz1 =  calzmean + caloffset;
      G4double distbCalor = 1.0*mm;
      G4double calor_L = 15.*cm;
      for (G4int j2=0; j2<12; j2++) {
        G4double postmp = 0;
        if (j2 < 4) {
          distbCalor = 1.0*cm;
          postmp = calz1 - calor_L + j2*(distbCalor + dy1);
        } else {
          distbCalor = 0.5*cm;
          postmp = calz1 - calor_L + 3.*(1.0*cm + dy1) + (j2 - 3)*(distbCalor + dy1);
        }
        G4VPhysicalVolume* physiCalor = new G4PVPlacement(rmtmp,
                G4ThreeVector(calx1, caly1, postmp),
                logicCalor,      // its logical volume
                "Calorimeter",         // its name
                World,               // its mother  volume
                false,           // no boolean operations
                j2*12 + j);              // copy number
        physiCalor = new G4PVPlacement(rmtmp,
                G4ThreeVector(calx1, caly1, -postmp),
                logicCalor,      // its logical volume
                "Calorimeter",         // its name
                World,               // its mother  volume
                false,           // no boolean operations
                136 + j2*12 + j);              // copy number
        if (j2 < 3) {
            G4VPhysicalVolume* physiCalorLead10mm = new G4PVPlacement(rmtmp,
                    G4ThreeVector(calx1, caly1, postmp + 0.5*dy1 + 0.5*calorLead10mm_D),
                    logicCalorLead10mm,      // its logical volume
                    "CalorimeterLead10mm",         // its name
                    World,               // its mother  volume
                    false,           // no boolean operations
                    j2*12 + j);              // copy number
            physiCalorLead10mm = new G4PVPlacement(rmtmp,
                    G4ThreeVector(calx1, caly1, -(postmp + 0.5*dy1 + 0.5*calorLead10mm_D)),
                    logicCalorLead10mm,      // its logical volume
                    "CalorimeterLead10mm",         // its name
                    World,               // its mother  volume
                    false,           // no boolean operations
                    36 + j2*12 + j);              // copy number
        } else if(j2 < 11) {
            G4VPhysicalVolume* physiCalorLead5mm = new G4PVPlacement(rmtmp,
                    G4ThreeVector(calx1, caly1, postmp + 0.5*dy1 + 0.5*calorLead5mm_D),
                    logicCalorLead5mm,      // its logical volume
                    "CalorimeterLead5mm",         // its name
                    World,               // its mother  volume
                    false,           // no boolean operations
                    j2*12 + j);              // copy number
            physiCalorLead5mm = new G4PVPlacement(rmtmp,
                    G4ThreeVector(calx1, caly1, -(postmp + 0.5*dy1 + 0.5*calorLead5mm_D)),
                    logicCalorLead5mm,      // its logical volume
                    "CalorimeterLead5mm",         // its name
                    World,               // its mother  volume
                    false,           // no boolean operations
                    96+j2*12+j);              // copy number
        }
      }
    }
    logicCalor->SetVisAttributes(G4Colour(0.0,0.0,1.0));
    logicCalorLead5mm->SetVisAttributes(G4Colour(1.0,0.7,0.0));
    logicCalorLead10mm->SetVisAttributes(G4Colour(1.0,0.7,0.0));
  }

//----------- Lead Wall -------------------
  G4double  zWall = 60.*cm;

  if (withMattoni) {
    G4double pbx = 20.0*cm;
    G4double pby = 4.9*cm;
    G4double pbz = 10.0*cm;

    G4Box* mattone = new G4Box("Mattone", 0.5*pbx, 0.5*pby, 0.5*pbz);
    G4LogicalVolume* logicMattone = new G4LogicalVolume(mattone, Pb, "World", 0, 0, 0);
    G4VPhysicalVolume* physiMattone;

    G4int ncopymatt = 0;
    G4int ncopymax = 6;
    for(G4int z=0; z<2; z++) {
      G4double zpos = 0.0;
      if (z == 0)
        zpos =  zWall;
      if (z == 1)
        zpos = -zWall;
      for(G4int j=0;j<5;j++) {
        for(G4int i=0; i<ncopymax-1; i+=2) {
          physiMattone = new G4PVPlacement(0, G4ThreeVector(-0.5*(i + 1)*pbx, -0.5*(3 + j*4)*pby, zpos), logicMattone,
                                           "Mattone", World, false, ++ncopymatt);
          physiMattone = new G4PVPlacement(0, G4ThreeVector(0.5*(i + 1)*pbx, -0.5*(3 + j*4)*pby, zpos), logicMattone,
                                           "Mattone", World, false, ++ncopymatt);
          physiMattone = new G4PVPlacement(0, G4ThreeVector(-0.5*(i + 1)*pbx, 0.5*(3 + j*4)*pby, zpos), logicMattone,
                                           "Mattone", World, false, ++ncopymatt);
          physiMattone = new G4PVPlacement(0, G4ThreeVector(0.5*(i + 1)*pbx, 0.5*(3 + j*4)*pby, zpos), logicMattone,
                                           "Mattone", World, false, ++ncopymatt);
          physiMattone = new G4PVPlacement(0, G4ThreeVector(0.5*(i + 2)*pbx, 0.5*(5 + j*4)*pby, zpos), logicMattone,
                                           "Mattone", World, false, ++ncopymatt);
          physiMattone = new G4PVPlacement(0, G4ThreeVector(-0.5*(i + 2)*pbx, 0.5*(5 + j*4)*pby, zpos), logicMattone,
                                           "Mattone", World, false, ++ncopymatt);
          physiMattone = new G4PVPlacement(0, G4ThreeVector(0.5*(i + 2)*pbx, -0.5*(5 + j*4)*pby, zpos), logicMattone,
                                        "Mattone", World, false, ++ncopymatt);
          physiMattone = new G4PVPlacement(0, G4ThreeVector(-0.5*(i + 2)*pbx, -0.5*(5 + j*4)*pby, zpos), logicMattone,
                                           "Mattone", World, false, ++ncopymatt);
        }
        physiMattone = new G4PVPlacement(0, G4ThreeVector(0.5*(ncopymax + 1)*pbx, 0.5*(3 + j*4)*pby, zpos), logicMattone,
                                         "Mattone", World, false, ++ncopymatt);
        physiMattone = new G4PVPlacement(0, G4ThreeVector(-0.5*(ncopymax + 1)*pbx, 0.5*(3 + j*4)*pby, zpos), logicMattone,
                                         "Mattone", World, false, ++ncopymatt);
        physiMattone = new G4PVPlacement(0, G4ThreeVector(0.5*(ncopymax + 1)*pbx, -0.5*(3 + j*4)*pby, zpos), logicMattone,
                                         "Mattone", World, false, ++ncopymatt);
        physiMattone = new G4PVPlacement(0, G4ThreeVector(-0.5*(ncopymax + 1)*pbx, -0.5*(3 + j*4)*pby, zpos), logicMattone,
                                         "Mattone", World, false, ++ncopymatt);
        physiMattone = new G4PVPlacement(0, G4ThreeVector(0, 0.5*(5 + j*4)*pby, zpos), logicMattone,
                                         "Mattone", World, false, ++ncopymatt);
        physiMattone = new G4PVPlacement(0, G4ThreeVector(0, -0.5*(5 + j*4)*pby, zpos), logicMattone,
                                         "Mattone", World, false, ++ncopymatt);
      }
      G4double xoffset = 7.0*cm;
      for(G4int i=0; i<ncopymax; i+=2) {
        physiMattone = new G4PVPlacement(0, G4ThreeVector(-0.5*(i + 1)*pbx - xoffset, -0.5*pby, zpos), logicMattone,
                                         "Mattone", World, false, ++ncopymatt);
        physiMattone = new G4PVPlacement(0, G4ThreeVector(0.5*(i + 1)*pbx + xoffset, -0.5*pby, zpos), logicMattone,
                                         "Mattone", World, false, ++ncopymatt);
        physiMattone = new G4PVPlacement(0, G4ThreeVector(-0.5*(i + 1)*pbx - xoffset, 0.5*pby, zpos), logicMattone,
                                         "Mattone", World, false, ++ncopymatt);
        physiMattone = new G4PVPlacement(0, G4ThreeVector(0.5*(i + 1)*pbx + xoffset, 0.5*pby, zpos), logicMattone,
                                         "Mattone", World, false, ++ncopymatt);
      }
    }
    G4VisAttributes* visMattone= new G4VisAttributes(0);
    visMattone->SetColor(1.0,1.0,1.0,0.1);
    logicMattone->SetVisAttributes(visMattone);
  } else if (withWall) {
    G4double wallz = 9.8*cm;
    G4double wally = 1.0*m;
    G4double wallx = 1.0*m;
    G4double wallz2 = wallz;
    G4double wally2 = 10.0*cm;
    G4double wallx2 = 10.0*cm;

    G4Box* solidWall1 = new G4Box("Wall1", 0.5*wallx, 0.5*wally, 0.5*wallz);
    G4Box* solidWall2 = new G4Box("Wall2", 0.5*wallx2, 0.5*wally2, 0.5*wallz2);
    G4SubtractionSolid* solidWall = new G4SubtractionSolid("Wall", solidWall1, solidWall2);
    G4LogicalVolume* logicWall = new G4LogicalVolume(solidWall, Pb, "World", 0, 0, 0);
    G4VPhysicalVolume* physiWall = new G4PVPlacement(0,
            G4ThreeVector(0.,0.,-zWall),
            logicWall,      // its logical volume
            "Wall",         // its name
            World,               // its mother  volume
            false,           // no boolean operations
            0);              // copy number
    physiWall = new G4PVPlacement(0,
            G4ThreeVector(0.,0.,zWall),
            logicWall,      // its logical volume
            "Wall",         // its name
            World,               // its mother  volume
            false,           // no boolean operations
            1);              // copy number

    G4VisAttributes* visWall = new G4VisAttributes(1);
    visWall->SetColor(1.0,1.0,1.0,0.);
    logicWall->SetVisAttributes(visWall);
  }

// Kaon Detector Top //
  G4double dx_kdtop;
  G4double dy_kdtop;
  G4double dz_kdtop;
  G4double pos_kdtop;

  if (SiddhartaSetup == 2) {
    dx_kdtop = mycard->variables["dx_kmtop"];
    dy_kdtop = mycard->variables["dy_kmtop"];
    dz_kdtop = 1450.*um;
    pos_kdtop = mycard->variables["z_kmtop"];
  } else if (SiddhartaSetup == 3) {
    dx_kdtop = 120.0*mm;
    dy_kdtop = 120.0*mm;
    dz_kdtop = 1450.*um;
    pos_kdtop = 102.*mm;
  } else if (SiddhartaSetup == 4) {
    dx_kdtop = 120.0*mm;
    dy_kdtop = 120.0*mm;
    dz_kdtop = 1450.*um;
    pos_kdtop = 115.*mm;
  } else if (SiddhartaSetup == 5) {
    dx_kdtop = 120.0*mm;
    dy_kdtop = 120.0*mm;
    dz_kdtop = 1450.*um;
    pos_kdtop = 115.*mm;
  } else if (SiddhartaSetup == 6 || SiddhartaSetup == 7) {
    dx_kdtop = 120.0*mm;
    dy_kdtop = 120.0*mm;
    dz_kdtop = 1450.*um;
    pos_kdtop = 115.*mm;
  } else if (SiddhartaSetup == 8) {
    dx_kdtop = 120.0*mm;
    dy_kdtop = 120.0*mm;
    dz_kdtop = 1450.*um;
    pos_kdtop = 115.*mm;
  } else if (SiddhartaSetup == 2020) {
    dx_kdtop = 100.0*mm;
    dy_kdtop = 100.0*mm;
    dz_kdtop = 1500.*um;
    pos_kdtop = 120.*mm;
  } else {
    dx_kdtop = 60.0*mm;
    dy_kdtop = 49.0*mm;
    dz_kdtop = 1450.*um;
    pos_kdtop = 60.*mm;
  }

// Scintillator
  G4Box* solidKaonDetectorTop = new G4Box("KaonDetectorTop", 0.5*dx_kdtop, 0.5*dy_kdtop, 0.5*dz_kdtop);
  G4LogicalVolume* logicKaonDetectorTop = new G4LogicalVolume(solidKaonDetectorTop, BC420, "World", 0, 0, 0);
  G4VPhysicalVolume* physiKaonDetectorTop = new G4PVPlacement(rm01, G4ThreeVector(0.,pos_kdtop,0.),
                                                              logicKaonDetectorTop, "KaonDetectorTop", World, false, 0);
  logicKaonDetectorTop->SetVisAttributes(G4Colour(0.,1.,0.));

// Aluminum potion in the Aluminized Mylar
  G4double kmtop_alt_x = dx_kdtop;
  G4double kmtop_alt_y = dy_kdtop;
  G4double km_alt_w = 2.*um;

  G4Box* solidKaonDetectorTopAl = new G4Box("KaonDetectorTopAl", 0.5*kmtop_alt_x, 0.5*kmtop_alt_y, 0.5*km_alt_w);
  G4LogicalVolume* logicKaonDetectorTopAl = new G4LogicalVolume(solidKaonDetectorTopAl, Al, "World", 0, 0, 0);
  G4VPhysicalVolume* physiKaonDetectorTopAl = new G4PVPlacement(rm01, G4ThreeVector(0., pos_kdtop + 0.5*dz_kdtop + 0.5*km_alt_w, 0.),
                                                                logicKaonDetectorTopAl, "KaonDetectorTopAl", World, false, 0);
  physiKaonDetectorTopAl = new G4PVPlacement(rm01,G4ThreeVector(0., pos_kdtop - 0.5*dz_kdtop - 0.5*km_alt_w, 0.),
                                             logicKaonDetectorTopAl, "KaonDetectorTopAl", World, false, 1);
  logicKaonDetectorTopAl->SetVisAttributes(G4Colour(1.,1.,0.));

// Mylar part in the Aluminized Mylar
  G4double km_my_w = 8.0*um;

  G4Box* solidKaonDetectorTopMy = new G4Box("KaonDetectorTopMy", 0.5*kmtop_alt_x, 0.5*kmtop_alt_y, 0.5*km_my_w);
  G4LogicalVolume* logicKaonDetectorTopMy = new G4LogicalVolume(solidKaonDetectorTopMy, mylar, "World", 0, 0, 0);
  G4VPhysicalVolume* physiKaonDetectorTopMy = new G4PVPlacement(rm01,G4ThreeVector(0., pos_kdtop + 0.5*dz_kdtop + km_alt_w + 0.5*km_my_w, 0.),
                                                                logicKaonDetectorTopMy, "KaonDetectorTopMy", World, false, 0);
  physiKaonDetectorTopMy = new G4PVPlacement(rm01, G4ThreeVector(0., pos_kdtop - 0.5*dz_kdtop - km_alt_w - 0.5*km_my_w, 0.)      ,
                                             logicKaonDetectorTopMy, "KaonDetectorTopMy", World, false, 1);
  logicKaonDetectorTopMy->SetVisAttributes(G4Colour(1.,1.,0.));

// Paper
  G4double km_paper_w = 50.0*um;

  G4Box* solidKaonDetectorTopPaper = new G4Box("KaonDetectorTopPaper",0.5*kmtop_alt_x,0.5*kmtop_alt_y,0.5*km_paper_w);
  G4LogicalVolume* logicKaonDetectorTopPaper = new G4LogicalVolume(solidKaonDetectorTopPaper, paper, "World", 0, 0, 0);
  G4VPhysicalVolume* physiKaonDetectorTopPaper =
                        new G4PVPlacement(rm01,G4ThreeVector(0., pos_kdtop + 0.5*dz_kdtop + km_alt_w + km_my_w + 0.5*km_paper_w, 0.),
                                          logicKaonDetectorTopPaper, "KaonDetectorTopPaper", World, false, 0);
  physiKaonDetectorTopPaper = new G4PVPlacement(rm01,
                                                G4ThreeVector(0., pos_kdtop - 0.5*dz_kdtop - km_alt_w - km_my_w - 0.5*km_paper_w, 0.),
                                                logicKaonDetectorTopPaper, "KaonDetectorTopPaper", World, false, 1);
  logicKaonDetectorTopPaper->SetVisAttributes(G4Colour(1.,1.,0.));

// Blacktape from PVC
  G4double km_bt_w = 250.0*um;

  G4Box* solidKaonDetectorTopBT = new G4Box("KaonDetectorTopBT", 0.5*kmtop_alt_x, 0.5*kmtop_alt_y, 0.5*km_bt_w);
  G4LogicalVolume* logicKaonDetectorTopBT = new G4LogicalVolume(solidKaonDetectorTopBT, pvc, "World", 0, 0, 0);
  G4VPhysicalVolume* physiKaonDetectorTopBT =
                        new G4PVPlacement(rm01, G4ThreeVector(0., pos_kdtop + 0.5*dz_kdtop + km_alt_w + km_my_w + km_paper_w + 0.5*km_bt_w, 0.),
                                          logicKaonDetectorTopBT, "KaonDetectorTopBT", World, false, 0);
  physiKaonDetectorTopBT = new G4PVPlacement(rm01, G4ThreeVector(0., pos_kdtop - 0.5*dz_kdtop - km_alt_w - km_my_w - km_paper_w - 0.5*km_bt_w, 0.),
                                             logicKaonDetectorTopBT, "KaonDetectorTopBT", World, false, 1);
  logicKaonDetectorTopBT->SetVisAttributes(G4Colour(1.,1.,0.));
  logicKaonDetectorTop->SetSensitiveDetector(KDTop_SD);

// Kaon Detector Bottom //
  G4double dx_kdbottom;
  G4double dy_kdbottom;
  G4double dz_kdbottom;
  G4double pos_kdbottom;
  G4double dx_kdbottom11;
  G4double dy_kdbottom11;
  G4double dx_kdbottom12;
  G4double dy_kdbottom12;
  G4double dz_kdbottom3;

  if ( SiddhartaSetup == 2020 ) {
    dx_kdbottom = 100.0*mm;
    dy_kdbottom = 100.0*mm;
    dz_kdbottom = 1500.*um;
    pos_kdbottom = -120.*mm;
    dx_kdbottom11 = dy_kdbottom;
    dy_kdbottom11 = dz_kdbottom;
    dx_kdbottom12 = 70.0*mm;
    dy_kdbottom12 = dz_kdbottom;
    dz_kdbottom3 = 70.0*mm;
  } else {
    if ( SiddhartaSetup != 1 ) {
      dx_kdbottom = 90.0*mm;
      dy_kdbottom = 144.0*mm;
      dz_kdbottom = 1450.*um;
      pos_kdbottom = -60.*mm;
      dx_kdbottom11 = dy_kdbottom;
      dy_kdbottom11 = dz_kdbottom;
      dx_kdbottom12 = 70.0*mm;
      dy_kdbottom12 = dz_kdbottom;
      dz_kdbottom3 = 70.0*mm;
    } else {
      dx_kdbottom = 72.0*mm;
      dy_kdbottom = 72.0*mm;
      dz_kdbottom = 1450.*um;
      pos_kdbottom = -60.*mm;
      dx_kdbottom11 = dy_kdbottom;
      dy_kdbottom11 = dz_kdbottom;
      dx_kdbottom12 = 49.0*mm;
      dy_kdbottom12 = dz_kdbottom;
      dz_kdbottom3 = 0.5*(152.0 - dx_kdbottom)*mm;
    }
  }

// Scintillator
// NEW MOVE THE KT BOTTOM 2 cm in BOOST DIRECTION --> on X: +2cm --> 20 mm
  G4double boostMove = mycard->variables["kaonBottomBoostShift"]*cm;
  G4double stepToVis = 0.1*cm;

  G4cout << "-------> <------" << G4endl;
  G4cout << "Kaon trigger bottom shift equal to " << boostMove << " mm" << G4endl;
  G4cout << "->  <- step with size " << stepToVis << " mm" << G4endl;
  for (double i=0; i<boostMove; i+=stepToVis) {
    G4cout << "->";
  }
  G4cout << " |=KTBottom=| ";
  for (double i=0; i>boostMove; i-=stepToVis) {
    G4cout << "<-";
  }
  G4cout << G4endl;
  G4cout << "-------> <------" << G4endl;
  G4cout << G4endl;

  G4double verticalMove = mycard->variables["kaonBottomVerticalShift"]*cm;
  pos_kdbottom += verticalMove;

  G4cout << "^^^^^^ vvvvvvvvv" << G4endl;
  G4cout << "Kaon trigger vertical shift equal to " << verticalMove << " mm" << G4endl;
  G4cout << "^  v step with size " << stepToVis << " mm" << G4endl;
  for (double i=0; i<verticalMove; i+=stepToVis) {
    G4cout << "v";
  }
  G4cout << " |=KTBottom=| ";
  for (double i=0; i>verticalMove; i-=stepToVis) {
    G4cout << "^";
  }
  G4cout << G4endl;
  G4cout << "^^^^^^ vvvvvvvvv" << G4endl;

  G4Box* solidKaonDetectorBottom0 = new G4Box("KaonDetectorBottom0", 0.5*dx_kdbottom, 0.5*dy_kdbottom, 0.5*dz_kdbottom);
  //G4Trd* solidKaonDetectorBottom2 = new G4Trd("KaonDetectorBottom2", 0.5*dx_kdbottom11, 0.5*dx_kdbottom12, 0.5*dy_kdbottom11,
   //                                           0.5*dy_kdbottom12, 0.5*dz_kdbottom3);
  //G4UnionSolid* solidKaonDetectorBottom3 = new G4UnionSolid("KaonDetectorBottom3", solidKaonDetectorBottom0, solidKaonDetectorBottom2, rm34,
  //                                                          G4ThreeVector(0.5*dx_kdbottom + 0.5*dz_kdbottom3,0*cm,0*cm));
 // G4UnionSolid* solidKaonDetectorBottom4 = new G4UnionSolid("KaonDetectorBottom4", solidKaonDetectorBottom3, solidKaonDetectorBottom2, rm35,
 //                                                           G4ThreeVector(-0.5*dx_kdbottom - 0.5*dz_kdbottom3,0*cm,0*cm));

  G4LogicalVolume* logicKaonDetectorBottom = new G4LogicalVolume(solidKaonDetectorBottom0, BC420, "World", 0, 0, 0);
  G4VPhysicalVolume* physiKaonDetectorBottom = new G4PVPlacement(rm01, G4ThreeVector(boostMove, pos_kdbottom,0.), logicKaonDetectorBottom,
                                                                "KaonDetectorBottom", World, false, 0); 
  logicKaonDetectorBottom->SetSensitiveDetector(KDBottom_SD);
  logicKaonDetectorBottom->SetVisAttributes(G4Colour(1.,1.,0.));

// Aluminum in Aluminized Mylar
  G4double km_alb_x = dx_kdbottom;
  G4double km_alb_y = dy_kdbottom;
  G4double km_alb_w = 2.*um;

  G4Box* solidKaonDetectorBottomAl = new G4Box("KaonDetectorBottomAl", 0.5*km_alb_x, 0.5*km_alb_y, 0.5*km_alb_w);
  G4LogicalVolume* logicKaonDetectorBottomAl = new G4LogicalVolume(solidKaonDetectorBottomAl, Al, "World", 0, 0, 0);
  G4VPhysicalVolume* physiKaonDetectorBottomAl = new G4PVPlacement(rm01, G4ThreeVector(boostMove, pos_kdbottom + 0.5*dz_kdbottom + 0.5*km_alb_w,0.),
                                                                   logicKaonDetectorBottomAl, "KaonDetectorBottomAl", World, false, 0);
  physiKaonDetectorBottomAl = new G4PVPlacement(rm01, G4ThreeVector(boostMove, pos_kdbottom - 0.5*dz_kdbottom - 0.5*km_alb_w,0.),
                                                logicKaonDetectorBottomAl, "KaonDetectorBottomAl", World, false, 1);
  logicKaonDetectorBottomAl->SetVisAttributes(G4Colour(1.,1.,0.));

// Mylar in Aluminized Mylar
  G4double km_myb_w = km_my_w;   // 2020 same thickness as top km mylar

  G4Box* solidKaonDetectorBottomMy = new G4Box("KaonDetectorBottomMy", 0.5*km_alb_x, 0.5*km_alb_y, 0.5*km_myb_w);
  G4LogicalVolume* logicKaonDetectorBottomMy = new G4LogicalVolume(solidKaonDetectorBottomMy, mylar, "World", 0, 0, 0);
  G4VPhysicalVolume* physiKaonDetectorBottomMy = new G4PVPlacement(rm01, G4ThreeVector(boostMove, pos_kdbottom + 0.5*dz_kdbottom + km_alb_w + 0.5*km_myb_w,0.),
                                                                   logicKaonDetectorBottomMy, "KaonDetectorBottomMy", World, false, 0);
  physiKaonDetectorBottomMy = new G4PVPlacement(rm01, G4ThreeVector(boostMove, pos_kdbottom - 0.5*dz_kdbottom - km_alb_w - 0.5*km_myb_w,0.),
                                                logicKaonDetectorBottomMy, "KaonDetectorBottomMy", World, false, 1);
  logicKaonDetectorBottomMy->SetVisAttributes(G4Colour(1.,1.,0.));

// Paper
  G4Box* solidKaonDetectorBotPaper = new G4Box("KaonDetectorBotPaper", 0.5*km_alb_x, 0.5*km_alb_y, 0.5*km_paper_w);
  G4LogicalVolume* logicKaonDetectorBotPaper = new G4LogicalVolume(solidKaonDetectorBotPaper, paper, "World", 0, 0, 0);
  G4VPhysicalVolume* physiKaonDetectorBotPaper =
                                new G4PVPlacement(rm01, G4ThreeVector(boostMove, pos_kdbottom + 0.5*dz_kdbottom + km_alb_w + km_myb_w + 0.5*km_paper_w,0.),
                                                  logicKaonDetectorBotPaper, "KaonDetectorBotPaper", World, false, 0);
  physiKaonDetectorBotPaper = new G4PVPlacement(rm01, G4ThreeVector(boostMove, pos_kdbottom - 0.5*dz_kdbottom - km_alb_w - km_myb_w - 0.5*km_paper_w,0.),
                                                logicKaonDetectorBotPaper, "KaonDetectorBotPaper", World, false, 1);
  logicKaonDetectorBotPaper->SetVisAttributes(G4Colour(1.,1.,0.));

// Blacktape from PVC
  G4double km_btb_w = km_bt_w;    // 2020 same as top KM

  G4Box* solidKaonDetectorBottomBT = new G4Box("KaonDetectorBottomBT", 0.5*km_alb_x, 0.5*km_alb_y, 0.5*km_btb_w);
  G4LogicalVolume* logicKaonDetectorBottomBT = new G4LogicalVolume(solidKaonDetectorBottomBT, pvc, "World", 0, 0, 0);
  G4VPhysicalVolume* physiKaonDetectorBottomBT =
                        new G4PVPlacement(rm01, G4ThreeVector(boostMove, pos_kdbottom + 0.5*dz_kdbottom + km_alb_w + km_myb_w + km_paper_w + 0.5*km_btb_w,0.),
                                          logicKaonDetectorBottomBT, "KaonDetectorBottomBT", World, false, 0);
  physiKaonDetectorBottomBT =
                        new G4PVPlacement(rm01, G4ThreeVector(boostMove, pos_kdbottom - 0.5*dz_kdbottom - km_alb_w - km_myb_w - km_paper_w - 0.5*km_btb_w,0.),
                                          logicKaonDetectorBottomBT, "KaonDetectorBottomBT", World, false, 1);
  logicKaonDetectorBottomBT->SetVisAttributes(G4Colour(1.,1.,0.));


/////////////////////////////////
///   Luminometer    2020      //
/////////////////////////////////
  G4double posx_lumi_anti  = 72*mm;
  G4double posx_lumi_boost = 102*mm;
  G4double dx_lumi = 80*mm;
  G4double dy_lumi = 40*mm;
  G4double dz_lumi = 2.*mm;

  G4double lumi_al_dx = dx_lumi;
  G4double lumi_al_dy = dy_lumi;
  G4double lumi_al_dz = 2.*um;

  G4double lumi_my_dx = dx_lumi;
  G4double lumi_my_dy = dy_lumi;
  G4double lumi_my_dz = 13.*um;

  G4double lumi_pc_dx = dx_lumi;
  G4double lumi_pc_dy = dy_lumi;
  G4double lumi_pc_dz = 25.*um;

////////////////////////////////
/// luminometer boost side
  G4Box* solidLumiDetectorBoost = new G4Box("LumiDetectorBoost", 0.5*dx_lumi, 0.5*dy_lumi, 0.5*dz_lumi);
  G4LogicalVolume* logicLumiDetectorBoost = new G4LogicalVolume(solidLumiDetectorBoost, BC420, "World", 0, 0, 0);
  G4VPhysicalVolume* physiLumiDetectorBoost = new G4PVPlacement(rmY90,G4ThreeVector(posx_lumi_boost,0.,0.), logicLumiDetectorBoost,
                                                                "LumiDetectorBoost", World, false, 0);
  logicLumiDetectorBoost->SetVisAttributes(G4Colour(0.,1.,0.));

// Alumi in Aluminized Mylar, overall 15 micron
  G4Box* solidLumiDetectorBoostAl = new G4Box("LumiDetectorBoostAl", 0.5*lumi_al_dx, 0.5*lumi_al_dy, 0.5*lumi_al_dz);
  G4LogicalVolume* logicLumiDetectorBoostAl = new G4LogicalVolume(solidLumiDetectorBoostAl, Al, "World", 0, 0, 0);
  G4VPhysicalVolume* physiLumiDetectorBoostAl = new G4PVPlacement(rmY90, G4ThreeVector(posx_lumi_boost + 0.5*dz_lumi + 0.5*lumi_al_dz,0.,0.),
                                                                  logicLumiDetectorBoostAl, "LumiDetectorBoostAl", World, false, 0);
  physiLumiDetectorBoostAl = new G4PVPlacement(rmY90,G4ThreeVector(posx_lumi_boost - 0.5*dz_lumi - 0.5*lumi_al_dz,0.,0.),
                                               logicLumiDetectorBoostAl, "LumiDetectorBoostAl", World, false, 1);
  logicLumiDetectorBoostAl->SetVisAttributes(G4Colour(1.,1.,0.));

// Mylar in Aluminized Mylar
  G4Box* solidLumiDetectorBoostMy = new G4Box("LumiDetectorBoostMy", 0.5*lumi_my_dx, 0.5*lumi_my_dy, 0.5*lumi_my_dz );
  G4LogicalVolume* logicLumiDetectorBoostMy = new G4LogicalVolume(solidLumiDetectorBoostMy, mylar, "World", 0, 0, 0);
  G4VPhysicalVolume* physiLumiDetectorBoostMy =
                                        new G4PVPlacement(rmY90, G4ThreeVector(posx_lumi_boost + 0.5*dz_lumi + lumi_al_dz + 0.5*lumi_my_dz,0.,0.),
                                                          logicLumiDetectorBoostMy, "LumiDetectorBoostMy", World, false, 0);
  physiLumiDetectorBoostMy = new G4PVPlacement(rmY90, G4ThreeVector(posx_lumi_boost - 0.5*dz_lumi - lumi_al_dz - 0.5*lumi_my_dz,0.,0.),
                                               logicLumiDetectorBoostMy, "LumiDetectorBoostMy", World, false, 1);
  logicLumiDetectorBoostMy->SetVisAttributes(G4Colour(1.,1.,0.));

// Pokalon - Polycarbonate?
  G4Box* solidLumiDetectorBoostPC = new G4Box("LumiDetectorBoostPC", 0.5*lumi_pc_dx, 0.5*lumi_pc_dy, 0.5*lumi_pc_dz );
  G4LogicalVolume* logicLumiDetectorBoostPC = new G4LogicalVolume(solidLumiDetectorBoostPC, kapton, "World", 0, 0, 0);
  G4VPhysicalVolume* physiLumiDetectorBoostPC =
                        new G4PVPlacement(rmY90, G4ThreeVector(posx_lumi_boost + 0.5*dz_lumi + lumi_al_dz + lumi_my_dz + 0.5*lumi_pc_dz,0.,0.),
                                          logicLumiDetectorBoostPC, "LumiDetectorBoostPC", World, false, 0);
  physiLumiDetectorBoostPC = new G4PVPlacement(rmY90, G4ThreeVector(posx_lumi_boost - 0.5*dz_lumi - lumi_al_dz - lumi_my_dz - 0.5*lumi_pc_dz,0.,0.),
                                               logicLumiDetectorBoostPC, "LumiDetectorBoostPC", World, false, 1);
  logicLumiDetectorBoostPC->SetVisAttributes(G4Colour(1.,1.,0.));
  logicLumiDetectorBoost->SetSensitiveDetector(LDBoost_SD);

///////////////////////////////////////////
/// luminometer anti-boost side
  G4Box* solidLumiDetectorAntiBoost = new G4Box("LumiDetectorAntiBoost", 0.5*dx_lumi, 0.5*dy_lumi, 0.5*dz_lumi);
  G4LogicalVolume* logicLumiDetectorAntiBoost = new G4LogicalVolume(solidLumiDetectorAntiBoost, BC420, "World", 0, 0, 0);
  G4VPhysicalVolume* physiLumiDetectorAntiBoost = new G4PVPlacement(rmY90, G4ThreeVector(-posx_lumi_anti,0.,0.),
                                                                    logicLumiDetectorAntiBoost, "LumiDetectorAntiBoost", World, false, 0);
  logicLumiDetectorAntiBoost->SetVisAttributes(G4Colour(0.,1.,0.));

// Alumi in Aluminized Mylar, overall 15 micron
  G4Box* solidLumiDetectorAntiBoostAl = new G4Box("LumiDetectorAntiBoostAl", 0.5*lumi_al_dx, 0.5*lumi_al_dy, 0.5*lumi_al_dz);
  G4LogicalVolume* logicLumiDetectorAntiBoostAl = new G4LogicalVolume(solidLumiDetectorAntiBoostAl, Al, "World", 0, 0, 0);
  G4VPhysicalVolume* physiLumiDetectorAntiBoostAl = new G4PVPlacement(rmY90, G4ThreeVector(-posx_lumi_anti + 0.5*dz_lumi + 0.5*lumi_al_dz,0.,0.),
                                                                      logicLumiDetectorAntiBoostAl, "LumiDetectorAntiBoostAl", World, false, 0);
  physiLumiDetectorAntiBoostAl = new G4PVPlacement(rmY90,G4ThreeVector(-posx_lumi_anti - 0.5*dz_lumi - 0.5*lumi_al_dz,0.,0.),
                                                   logicLumiDetectorAntiBoostAl, "LumiDetectorAntiBoostAl", World, false, 1);
  logicLumiDetectorAntiBoostAl->SetVisAttributes(G4Colour(1.,1.,0.));

// Mylar in Aluminized Mylar
  G4Box* solidLumiDetectorAntiBoostMy = new G4Box("LumiDetectorAntiBoostMy", 0.5*lumi_my_dx, 0.5*lumi_my_dy, 0.5*lumi_my_dz );
  G4LogicalVolume* logicLumiDetectorAntiBoostMy = new G4LogicalVolume(solidLumiDetectorAntiBoostMy, mylar, "World", 0, 0, 0);
  G4VPhysicalVolume* physiLumiDetectorAntiBoostMy =
                                        new G4PVPlacement(rmY90, G4ThreeVector(-posx_lumi_anti + 0.5*dz_lumi + lumi_al_dz + 0.5*lumi_my_dz,0.,0.),
                                                          logicLumiDetectorAntiBoostMy, "LumiDetectorAntiBoostMy", World, false, 0);
  physiLumiDetectorAntiBoostMy = new G4PVPlacement(rmY90, G4ThreeVector(-posx_lumi_anti - 0.5*dz_lumi - lumi_al_dz - 0.5*lumi_my_dz,0.,0.),
                                                   logicLumiDetectorAntiBoostMy, "LumiDetectorAntiBoostMy", World, false, 1);
  logicLumiDetectorAntiBoostMy->SetVisAttributes(G4Colour(1.,1.,0.));

// Pokalon - Polycarbonate?
  G4Box* solidLumiDetectorAntiBoostPC = new G4Box("LumiDetectorAntiBoostPC", 0.5*lumi_pc_dx, 0.5*lumi_pc_dy, 0.5*lumi_pc_dz );
  G4LogicalVolume* logicLumiDetectorAntiBoostPC = new G4LogicalVolume(solidLumiDetectorAntiBoostPC, kapton, "World", 0, 0, 0);
  G4VPhysicalVolume* physiLumiDetectorAntiBoostPC =
                            new G4PVPlacement(rmY90, G4ThreeVector(-posx_lumi_anti + 0.5*dz_lumi + lumi_al_dz + lumi_my_dz + 0.5*lumi_pc_dz,0.,0.),
                            logicLumiDetectorAntiBoostPC, "LumiDetectorAntiBoostPC", World, false, 0);
  physiLumiDetectorAntiBoostPC =
                                new G4PVPlacement(rmY90,G4ThreeVector(-posx_lumi_anti - 0.5*dz_lumi - lumi_al_dz - lumi_my_dz - 0.5*lumi_pc_dz,0.,0.),
                                                  logicLumiDetectorAntiBoostPC, "LumiDetectorAntiBoostPC", World, false, 1);
  logicLumiDetectorAntiBoostPC->SetVisAttributes(G4Colour(1.,1.,0.));
  logicLumiDetectorAntiBoost->SetSensitiveDetector(LDAntiBoost_SD);

// Rotation 90deg around Y axis -> Z = X, X = -Z
  mycard->SetLumiBoostPosition(G4ThreeVector(posx_lumi_boost, 0., 0), G4ThreeVector(0.5*dz_lumi, 0.5*dy_lumi, 0.5*dx_lumi));
  mycard->SetLumiAntiBoostPosition(G4ThreeVector(-posx_lumi_anti, 0., 0), G4ThreeVector(0.5*dz_lumi, 0.5*dy_lumi,0.5*dx_lumi ));

////////////////////////////////
/// KLIMAX AntiBoost PART
////////////////////////////////

/// KLIMAXDegAntiBoost
  G4double kldegAntiBoost_z = 0.*mm;  //toreplace_ab_deg//
  G4double posx_kldegAntiBoost = -posx_lumi_anti - 0.5*dz_lumi - lumi_al_dz - lumi_my_dz - lumi_pc_dz - 0.5*kldegAntiBoost_z;

/// KLIMAXTarget1AntiBoost
  G4double dx_target1_ab = 40*mm; //like in june-july 2023
  G4double dy_target1_ab = 43*mm; //like in june-july 2023
  G4double kltarget1AntiBoost_z = 0.07*mm;  //toreplace_ab_target1//
  G4double Lumi_Target1_distAntiBoost = 0*mm;
  G4double posx_kltarget1AntiBoost = posx_kldegAntiBoost - 0.5*kldegAntiBoost_z - Lumi_Target1_distAntiBoost- 0.5*kltarget1AntiBoost_z;

  dx_target1_ab = dx_lumi;
  dy_target1_ab = dy_lumi;
  //kltarget1AntiBoost_z += 0.2*mm; // Al box entrance window

  G4Box* solidKLIMAXTarget1AntiBoost = new G4Box("KLIMAXTarget1AntiBoost", 0.5*dx_target1_ab, 0.5*dy_target1_ab, 0.5*kltarget1AntiBoost_z);
  G4LogicalVolume* logicKLIMAXTarget1AntiBoost = new G4LogicalVolume(solidKLIMAXTarget1AntiBoost, Al, "World", 0, 0, 0);
  G4VPhysicalVolume* physiKLIMAXTarget1AntiBoost = new G4PVPlacement(rmY90,
                                                                     G4ThreeVector(posx_kltarget1AntiBoost,0.,0.),
                                                                     logicKLIMAXTarget1AntiBoost, "KLIMAXTarget1AntiBoost",
                                                                     World, false, 0, true);
  logicKLIMAXTarget1AntiBoost->SetVisAttributes(G4Colour(1.,1.,0.));
  logicKLIMAXTarget1AntiBoost->SetSensitiveDetector(KLT1AntiBoost_SD);

  G4cout << "TARGET 1 X: " << posx_kltarget1AntiBoost - 0.5*kltarget1AntiBoost_z << " " << posx_kltarget1AntiBoost << " ";
  G4cout << posx_kltarget1AntiBoost + 0.5*kltarget1AntiBoost_z << G4endl;
  G4cout << "TARGET 1 Y: " << -0.5*dy_target1_ab << " " << 0. << " " << 0.5*dy_target1_ab << G4endl;
  G4cout << "TARGET 1 Z: " << -0.5*dx_target1_ab << " " << 0. << " " << 0.5*dx_target1_ab << G4endl;

/// KLIMAXTarget2AntiBoost
  G4double dx_target2_ab = 40*mm; //like in june-july 2023
  G4double dy_target2_ab = 43*mm; //like in june-july 2023
  G4double kltarget2AntiBoost_z = 0.1*mm;  //toreplace_ab_target2//
  G4double Target1_Target2_distAntiBoost = 0*mm;
  G4double posx_kltarget2AntiBoost = posx_kltarget1AntiBoost - 0.5*kltarget1AntiBoost_z - Target1_Target2_distAntiBoost - 0.5*kltarget2AntiBoost_z;

  dx_target2_ab = dx_lumi;
  dy_target2_ab = dy_lumi;
  //kltarget2AntiBoost_z += 0.2*mm; // Al box entrance window

  G4Box* solidKLIMAXTarget2AntiBoost = new G4Box("KLIMAXTarget2AntiBoost", 0.5*dx_target2_ab, 0.5*dy_target2_ab, 0.5*kltarget2AntiBoost_z);
  G4LogicalVolume* logicKLIMAXTarget2AntiBoost = new G4LogicalVolume(solidKLIMAXTarget2AntiBoost, graphite, "World", 0, 0, 0);
  G4VPhysicalVolume* physiKLIMAXTarget2AntiBoost = new G4PVPlacement(rmY90, G4ThreeVector(posx_kltarget2AntiBoost, 0., 0.),
                                                                     logicKLIMAXTarget2AntiBoost, "KLIMAXTarget2AntiBoost",
                                                                     World, false, 0, true);
  logicKLIMAXTarget2AntiBoost->SetVisAttributes(G4Colour(1.,0.3,0.2));
  logicKLIMAXTarget2AntiBoost->SetSensitiveDetector(KLT2AntiBoost_SD);

  G4cout << "TARGET 2 X: " << posx_kltarget2AntiBoost - 0.5*kltarget2AntiBoost_z << " " << posx_kltarget2AntiBoost << " ";
  G4cout << posx_kltarget2AntiBoost+  0.5*kltarget2AntiBoost_z << G4endl;
  G4cout << "TARGET 2 Y: " << -0.5*dy_target2_ab << " " << 0. << " " << 0.5*dy_target2_ab << G4endl;
  G4cout << "TARGET 2 Z: " << -0.5*dx_target2_ab << " " << 0. << " " << 0.5*dx_target2_ab << G4endl;

/// KLIMAXBoxWindowAntiBoost
  G4double dx_boxwindow_ab = 40*mm; // like in june-july 2023
  G4double dy_boxwindow_ab = 43*mm;
  dx_boxwindow_ab = 80*mm;
  G4double klboxwindowAntiBoost_z = 0.27*mm;  //toreplace_ab_boxwindow//
  //klboxwindowAntiBoost_z += 0.2*mm; // Al box entrance window

  G4double posx_klboxwindowAntiBoost = -215*mm - 0.5*klboxwindowAntiBoost_z;

  G4Box* solidKLIMAXBoxWindowAntiBoost = new G4Box("KLIMAXBoxWindowAntiBoost", 0.5*dx_boxwindow_ab, 0.5*dy_boxwindow_ab, 0.5*klboxwindowAntiBoost_z);
  G4LogicalVolume* logicKLIMAXBoxWindowAntiBoost = new G4LogicalVolume(solidKLIMAXBoxWindowAntiBoost, Al, "World", 0, 0, 0);
  G4VPhysicalVolume* physiKLIMAXBoxWindowAntiBoost = new G4PVPlacement(rmY90, G4ThreeVector(posx_klboxwindowAntiBoost, 0., 0.),
                                                                       logicKLIMAXBoxWindowAntiBoost, "KLIMAXBoxWindowAntiBoost",
                                                                       World, false, 0, true);
  logicKLIMAXBoxWindowAntiBoost->SetVisAttributes(G4Colour(0.,0.,1.));
  logicKLIMAXBoxWindowAntiBoost->SetSensitiveDetector(KLT2AntiBoost_SD);

  G4cout << "Box Window X: " << posx_klboxwindowAntiBoost - 0.5*klboxwindowAntiBoost_z << " " << posx_klboxwindowAntiBoost << " ";
  G4cout << posx_klboxwindowAntiBoost + 0.5*klboxwindowAntiBoost_z << G4endl;
  G4cout << "Box Window Y: " << -0.5*dy_boxwindow_ab << " " << 0. << " " << 0.5*dy_boxwindow_ab << G4endl;
  G4cout << "Box Window Z: " << -0.5*dx_boxwindow_ab << " " << 0. << " " << 0.5*dx_boxwindow_ab << G4endl;

/// KLIMAXCZTAntiBoost
  G4double klcztAntiBoost_x = 26*mm; // like in june-july 2023
  klcztAntiBoost_x = 52*mm;
  G4double klcztAntiBoost_y = 30*mm;
  G4double klcztAntiBoost_z = 5.*mm;
  G4double CZT_BoxWindow_distAntiBoost = 4.*mm; // 3mm for the ABS frame + 1 mm for the screw

  G4double posx_klcztAntiBoost = posx_klboxwindowAntiBoost - 0.5*klboxwindowAntiBoost_z - CZT_BoxWindow_distAntiBoost- 0.5*klcztAntiBoost_z;

  G4Box* solidKLIMAXCZTAntiBoost = new G4Box("KLIMAXCZTAntiBoost", 0.5*klcztAntiBoost_x, 0.5*klcztAntiBoost_y, 0.5*klcztAntiBoost_z);
  G4LogicalVolume* logicKLIMAXCZTAntiBoost = new G4LogicalVolume(solidKLIMAXCZTAntiBoost, CZT, "World", 0, 0, 0);
  G4VPhysicalVolume* physiKLIMAXCZTAntiBoost = new G4PVPlacement(rmY90, G4ThreeVector(posx_klcztAntiBoost, 0., 0.),
                                                                 logicKLIMAXCZTAntiBoost, "KLIMAXCZTAntiBoost",
                                                                 World, false, 0, true);
  logicKLIMAXCZTAntiBoost->SetVisAttributes(G4Colour(0.,1.0,0.));
  logicKLIMAXCZTAntiBoost->SetSensitiveDetector(KLCZTAntiBoost_SD);

  G4cout << "CZT X: " << posx_klcztAntiBoost - 0.5*klcztAntiBoost_z << " " << posx_klcztAntiBoost << " " << posx_klcztAntiBoost + 0.5*klcztAntiBoost_z << G4endl;
  G4cout << "CZT Y: " << -0.5*klcztAntiBoost_y << " " << 0. << " " << 0.5*klcztAntiBoost_y << G4endl;
  G4cout << "CZT Z: " << -0.5*klcztAntiBoost_x << " " << 0. << " " << 0.5*klcztAntiBoost_x << G4endl;

/// KLIMAXCZT Shielding AntiBoost
  G4double klcztBoxShieldingDistance = 10.*mm; // 5mm from CZT to box wall + 5 mm of space (to be checked)
  G4double klcztShieldingAntiBoost_x = 50.*mm;
  G4double klcztShieldingAntiBoost_y = 100.*mm;
  G4double klcztShieldingAntiBoost_z = 200.*mm;
 
  G4double posy_KLCZTShieldAB12 = 0.*mm;
  G4double posx_KLCZTShieldAB12 = posx_klboxwindowAntiBoost - 0.5*klboxwindowAntiBoost_z - 0.5*klcztShieldingAntiBoost_z;
  G4double posz_KLCZTShieldAB12 = 0.5*klcztAntiBoost_x + 0.5*klcztShieldingAntiBoost_x + klcztBoxShieldingDistance;

  G4Box* solidKLIMAXCZTShieldingAntiBoost = new G4Box("KLIMAXCZTShieldingAntiBoost", 0.5*klcztShieldingAntiBoost_x,
                                                      0.5*klcztShieldingAntiBoost_y, 0.5*klcztShieldingAntiBoost_z);
  G4LogicalVolume* logicKLIMAXCZTShieldingAntiBoost = new G4LogicalVolume(solidKLIMAXCZTShieldingAntiBoost, Pb, "World", 0, 0, 0);
  G4VPhysicalVolume* physiKLIMAXCZTShieldingAntiBoost1 = new G4PVPlacement(rmY90, G4ThreeVector(posx_KLCZTShieldAB12, posy_KLCZTShieldAB12, posz_KLCZTShieldAB12),
                                                                           logicKLIMAXCZTShieldingAntiBoost, "KLIMAXCZTShieldingAntiBoost1",
                                                                           World, false, 0, true);
  G4VPhysicalVolume* physiKLIMAXCZTShieldingAntiBoost2 = new G4PVPlacement(rmY90, G4ThreeVector(posx_KLCZTShieldAB12, posy_KLCZTShieldAB12, -posz_KLCZTShieldAB12),
                                                                           logicKLIMAXCZTShieldingAntiBoost, "KLIMAXCZTShieldingAntiBoost2",
                                                                           World, false, 0, true);
  logicKLIMAXCZTShieldingAntiBoost->SetVisAttributes(G4Colour(0.5,0.5,0.5));

////////////////////////////////
/// KLIMAX Boost PART
////////////////////////////////

/// KLIMAXDegBoost
  G4double kldegBoost_z = 0.;  //toreplace_b_deg//
  G4double posx_kldegBoost = posx_lumi_boost + 0.5*dz_lumi + lumi_al_dz + lumi_my_dz + lumi_pc_dz + 0.5*kldegBoost_z;

  /// KLIMAXTargetBoost
  G4double dx_target_b = dx_lumi;
  G4double dy_target_b = dz_lumi;
  G4double kltargetBoost_z = 1.5*mm;  //toreplace_b_target//
  G4double Lumi_Target_distBoost = 0.;

  G4double posx_kltargetBoost = posx_kldegBoost + 0.5*kldegBoost_z + Lumi_Target_distBoost + 0.5*kltargetBoost_z;

  G4Box* solidKLIMAXTargetBoost = new G4Box("KLIMAXTargetBoost", 0.5*dx_lumi, 0.5*dy_lumi, 0.5*kltargetBoost_z);
  G4LogicalVolume* logicKLIMAXTargetBoost = new G4LogicalVolume(solidKLIMAXTargetBoost, Pb, "World", 0, 0, 0);
  G4VPhysicalVolume* physiKLIMAXTargetBoost = new G4PVPlacement(rmY90, G4ThreeVector(posx_kltargetBoost, 0., 0.),
                                                                logicKLIMAXTargetBoost, "KLIMAXTargetBoost", World, false, 0, true);
  logicKLIMAXTargetBoost->SetVisAttributes(G4Colour(0.,1.,1.));
  logicKLIMAXTargetBoost->SetSensitiveDetector(KLTBoost_SD);

  G4cout << "TARGET X: " << posx_kltargetBoost - 0.5*kltargetBoost_z << " " << posx_kltargetBoost << " " << posx_kltargetBoost + 0.5*kltargetBoost_z << G4endl;
  G4cout << "TARGET Y: " << -0.5*dy_lumi << " " << 0. << " " << 0.5*dy_lumi << G4endl;
  G4cout << "TARGET Z: " << -0.5*dx_lumi << " " << 0. << " " << 0.5*dx_lumi << G4endl;

/// KLIMAXHPGeBoost


  G4double klHPGeBoost_rad = 59.8*mm;
  G4double klHPGeBoost_z = 59.3*mm;
  G4double HPGe_IP_distBoost = 250.*mm + 43.*mm; // 250 mm for Main table + 43 mm before entrance window

  G4double posx_klHPGeBoost =  HPGe_IP_distBoost+ 0.5*klHPGeBoost_z;

  G4Tubs* solidKLIMAXHPGeBoost = new G4Tubs("KLIMAXHPGeBoost", 0.,0.5*klHPGeBoost_rad, 0.5*klHPGeBoost_z, 0., 360.*deg);
  G4LogicalVolume* logicKLIMAXHPGeBoost = new G4LogicalVolume(solidKLIMAXHPGeBoost, HPGe, "World", 0, 0, 0);
  G4VPhysicalVolume* physiKLIMAXHPGeBoost = new G4PVPlacement(rmY90, G4ThreeVector(posx_klHPGeBoost, 0., 0.),
                                                             logicKLIMAXHPGeBoost, "KLIMAXHPGeBoost", World, false, 0, true);
  logicKLIMAXHPGeBoost->SetVisAttributes(G4Colour(0.,1.,1.));
  logicKLIMAXHPGeBoost->SetSensitiveDetector(KLHPGeBoost_SD);

  G4cout << "HPGe Z: " << posx_klHPGeBoost - 0.5*klHPGeBoost_z << " " << posx_klHPGeBoost << " " << posx_klHPGeBoost + 0.5*klHPGeBoost_z << G4endl;

/// KLIMAX HPGe Shielding 
  
  G4double klHPGeShieldingBoost_x = 50.*mm;
  G4double klHPGeShieldingBoost_y = 100.*mm;
  G4double klHPGeShieldingBoost_z = 200.*mm;
 
  G4double posy_klHPGeShieldB12 = 0.*mm;
  G4double posx_klHPGeShieldB12 = posx_klHPGeBoost - 0.5*klHPGeBoost_z + 0.5*klHPGeShieldingBoost_z;
  G4double posz_klHPGeShieldB12 = 0.5*klHPGeShieldingBoost_x + 0.5*klHPGeShieldingBoost_x + 10.*mm;

  G4Box* solidKLIMAXHPGeShieldingBoost = new G4Box("KLIMAXHPGeShieldingBoost", 0.5*klHPGeShieldingBoost_x,
                                                      0.5*klHPGeShieldingBoost_y, 0.5*klHPGeShieldingBoost_z);
  G4LogicalVolume* logicKLIMAXHPGeShieldingBoost = new G4LogicalVolume(solidKLIMAXHPGeShieldingBoost, Pb, "World", 0, 0, 0);
  G4VPhysicalVolume* physiKLIMAXHPGeShieldingBoost1 = new G4PVPlacement(rmY90, G4ThreeVector(posx_klHPGeShieldB12, posy_klHPGeShieldB12, posz_klHPGeShieldB12),
                                                                           logicKLIMAXHPGeShieldingBoost, "KLIMAXHPGeShieldingBoost1",
                                                                           World, false, 0, true);
  G4VPhysicalVolume* physiKLIMAXHPGeShieldingBoost2 = new G4PVPlacement(rmY90, G4ThreeVector(posx_klHPGeShieldB12, posy_klHPGeShieldB12, -posz_klHPGeShieldB12),
                                                                           logicKLIMAXHPGeShieldingBoost, "KLIMAXHPGeShieldingBoost2",
                                                                           World, false, 0, true);

  G4double posy_klHPGeShieldB3 = 0.5*klHPGeBoost_rad + 0.5*klHPGeShieldingBoost_y;
  G4double posx_klHPGeShieldB3 = HPGe_IP_distBoost+ 0.5*klHPGeBoost_rad + 0.5*klHPGeShieldingBoost_x;
  G4double posz_klHPGeShieldB3 = 0.;

  G4VPhysicalVolume* physiKLIMAXHPGeShieldingBoost3 = new G4PVPlacement(rmXYZ90, G4ThreeVector(posx_klHPGeShieldB3, posy_klHPGeShieldB3, posz_klHPGeShieldB3),
                                                                           logicKLIMAXHPGeShieldingBoost, "KLIMAXHPGeShieldingBoost3",
                                                                           World, false, 0, true);

  G4double posy_klHPGeShieldB4 = 0.5*klHPGeBoost_rad + 0.5*klHPGeShieldingBoost_y;
  G4double posx_klHPGeShieldB4 = posx_klHPGeShieldB3 + 2.*klHPGeShieldingBoost_x;
  G4double posz_klHPGeShieldB4 = 0.;

  G4VPhysicalVolume* physiKLIMAXHPGeShieldingBoost4 = new G4PVPlacement(rmXYZ90, G4ThreeVector(posx_klHPGeShieldB4, posy_klHPGeShieldB4, posz_klHPGeShieldB4),
                                                                           logicKLIMAXHPGeShieldingBoost, "KLIMAXHPGeShieldingBoost4",
                                                                           World, false, 0, true);
  logicKLIMAXHPGeShieldingBoost->SetVisAttributes(G4Colour(0.5,0.5,0.5));

// Shielding //
  G4double pos_shield;
  G4double Z_SHIEL;
  G4double Y_SHIEL;
  G4double X_SHIEL;
  G4double r1_hole;
  G4double r2_hole;
  G4double r1_hole_polyethylene;
  G4double r2_hole_polyethylene;
  G4double r1_hole_Al;
  G4double r2_hole_Al;
  G4double r1_hole_Cu;
  G4double r2_hole_Cu;

  if (SiddhartaSetup != 2 && SiddhartaSetup != 4 && SiddhartaSetup != 5 && SiddhartaSetup != 6
        && SiddhartaSetup != 7 && SiddhartaSetup != 8 && SiddhartaSetup != 2020) {
    X_SHIEL = 500.*mm;
    Y_SHIEL = 310.*mm;
    Z_SHIEL = 52.0*mm;
    G4double z_shiel_Al = 5.0*mm;
    G4double SH_ANGLE = 120.*deg;
    G4double DY_SH = Z_SHIEL/tan(0.5*SH_ANGLE);
    G4double Y_SHIEL1 = Y_SHIEL;
    G4double Y_SHIEL2 = Y_SHIEL - 2.*DY_SH;
    G4double VX1 = 0;
    G4double VZ1 = 0;
    G4double VZ2 = 0.5*Y_SHIEL - 0.5*DY_SH;
    G4double dist_km_shield = 5*mm;

    if (SiddhartaSetup == 1) {
      pos_shield = 70*mm+z_shiel_Al+0.5*Z_SHIEL;
    } else if (SiddhartaSetup == 3) {
      pos_shield = bp_r3 + 0.5*Z_SHIEL + z_shiel_Al + 1*mm;
    } else {
      pos_shield = 70*mm + z_shiel_Al + 0.5*Z_SHIEL;
    }
    G4double VY1 = pos_shield;
    G4double VY2 = 0.;
    G4double VX2 = 0.;
    G4double VX3 = 0.;
    G4double L_SHIEL2 = 200.*mm;
    G4double L_SHIEL3 = 200.*mm;
    G4double angle0 = 60.*deg;
    G4double angle1 = -90.*deg;
    G4double angle2 = 0.*deg;
    G4double VY3 = 0.5*(L_SHIEL2 - 0.5*DY_SH)*sin(180.*deg - SH_ANGLE);
    G4double VZ3 = 0.5*(L_SHIEL2 - 0.5*DY_SH)*cos(180.*deg - SH_ANGLE);
    G4double VX4 = VX1 + VX2 + VX3;
    G4double VY4 = VY1 + VY2 + VY3;
    G4double VZ4 = VZ1 + VZ2 + VZ3;

    G4Trd* solidShielding1 = new G4Trd("Shielding1", 0.5*X_SHIEL, 0.5*X_SHIEL, 0.5*Y_SHIEL1, 0.5*Y_SHIEL2, 0.5*Z_SHIEL);
    G4LogicalVolume* logicShielding1 = new G4LogicalVolume(solidShielding1, Pb, "World", 0, 0, 0);
    G4VPhysicalVolume* physiShielding1 = new G4PVPlacement(rm01, G4ThreeVector(0.,pos_shield,0.), logicShielding1,
                                                           "Shielding1", World, false, 0);

    G4RotationMatrix* rm11 = new G4RotationMatrix();
    rm11->rotateX(180.*deg - SH_ANGLE);
    rm11->rotateY(angle1);
    rm11->rotateY(angle2);
    G4RotationMatrix* rm12 = new G4RotationMatrix();
    rm12->rotateX(180.*deg + SH_ANGLE);
    rm12->rotateY(angle1 + 180.*deg);
    rm12->rotateY(angle2);

    G4Trap* solidShielding2 = new G4Trap("Shielding2", X_SHIEL, Z_SHIEL, L_SHIEL2, L_SHIEL2 - DY_SH);
    G4LogicalVolume* logicShielding2 = new G4LogicalVolume(solidShielding2, Pb, "World", 0, 0, 0);
    G4VPhysicalVolume* physiShielding2 = new G4PVPlacement(rm11, G4ThreeVector(VX4,VY4,VZ4), logicShielding2,
                                                           "Shielding2", World, false, 1);
    physiShielding2 = new G4PVPlacement(rm12, G4ThreeVector(VX4,VY4,-VZ4), logicShielding2, "Shielding2", World, false, 2);

    r1_hole = 0.5*53.0*mm;
    r2_hole = 0.5*78.0*mm;
    r1_hole_polyethylene = 0.5*57.0*mm;
    r2_hole_polyethylene = 0.5*82.0*mm;
    r1_hole_Al = 0.5*67.0*mm;
    r2_hole_Al = 0.5*92.0*mm;
    r1_hole_Cu = 0.5*77.0*mm;
    r2_hole_Cu = 0.5*102.0*mm;

    G4Cons* solidWindowShield = new G4Cons("WindowShield", 0., r1_hole, 0., r2_hole, 0.5*Z_SHIEL, 0., 360.*deg);
    G4LogicalVolume* logicWindowShield = new G4LogicalVolume(solidWindowShield, Air, "Shielding1", 0, 0, 0);
    G4VPhysicalVolume* physiWindowShield = new G4PVPlacement(rm00, G4ThreeVector(0.,0.,0.), logicWindowShield,
                                                             "WindowShield", logicShielding1, false, 1);

    G4Cons* solidWindowShield2 = new G4Cons("WindowShield2", r1_hole, r1_hole_polyethylene, r2_hole, r2_hole_polyethylene,
                                            0.5*Z_SHIEL, 0., 360.*deg);
    G4LogicalVolume* logicWindowShield2 = new G4LogicalVolume(solidWindowShield2, Polyethylene, "Shielding1", 0, 0, 0);
    G4VPhysicalVolume* physiWindowShield2 = new G4PVPlacement(rm00, G4ThreeVector(0.,0.,0.), logicWindowShield2,
                                                              "WindowShield2", logicShielding1, false, 1);

    G4Cons* solidWindowShield3 = new G4Cons("WindowShield3", r1_hole_polyethylene, r1_hole_Al, r2_hole_polyethylene,
                                            r2_hole_Al, 0.5*Z_SHIEL, 0., 360.*deg);
    G4LogicalVolume* logicWindowShield3 = new G4LogicalVolume(solidWindowShield3, Al, "Shielding1", 0, 0, 0);
    G4VPhysicalVolume* physiWindowShield3 = new G4PVPlacement(rm00, G4ThreeVector(0.,0.,0.), logicWindowShield3,
                                                              "WindowShield3", logicShielding1, false, 1);

    G4Cons* solidWindowShield4 = new G4Cons("WindowShield4", r1_hole_Al, r1_hole_Cu, r2_hole_Al, r2_hole_Cu, 0.5*Z_SHIEL, 0., 360.*deg);
    G4LogicalVolume* logicWindowShield4 = new G4LogicalVolume(solidWindowShield4, Cu, "Shielding1", 0, 0, 0);
    G4VPhysicalVolume* physiWindowShield4 = new G4PVPlacement(rm00, G4ThreeVector(0.,0.,0.), logicWindowShield4,
                                                              "WindowShield4", logicShielding1, false, 1);

    G4double DY_SH2 = z_shiel_Al/tan(0.5*SH_ANGLE);
    G4double Y_SHIEL3 = Y_SHIEL1 + 2.*DY_SH2;

    G4Trd* solidShieldAl = new G4Trd("ShieldAl", 0.5*X_SHIEL, 0.5*X_SHIEL, 0.5*Y_SHIEL3, 0.5*Y_SHIEL1, 0.5*z_shiel_Al);
    G4LogicalVolume* logicShieldAl = new G4LogicalVolume(solidShieldAl, Al, "World", 0, 0, 0);
    G4VPhysicalVolume* physiShieldAl = new G4PVPlacement(rm01, G4ThreeVector(VX1,VY1 - 0.5*Z_SHIEL - 0.5*z_shiel_Al,VZ1), logicShieldAl,
                                                         "ShieldingAl", World, false, 0);

    G4Tubs* solidShield_Al_Hole = new G4Tubs("Shield_Al_Hole", 0., 70.*mm, 0.5*z_shiel_Al, 0., 360.*deg);
    G4LogicalVolume* logicShield_Al_Hole = new G4LogicalVolume(solidShield_Al_Hole, Air, "World", 0, 0, 0);
    G4VPhysicalVolume* physiShield_Al_Hole = new G4PVPlacement(rm00, G4ThreeVector(0.,0.,0.), logicShield_Al_Hole,
                                                               "Shield_Al_Hole", logicShieldAl, false, 0);

    G4double VY1P = pos_shield - 0.5*Z_SHIEL - 0.5*z_shiel_Al;
    G4double VY3P = 0.5*(L_SHIEL3 - 0.5*DY_SH2)*sin(180.*deg - SH_ANGLE);

    G4double VZ2P = 0.5*Y_SHIEL + 0.5*DY_SH2;
    G4double VZ3P = 0.5*(L_SHIEL3 - 0.5*DY_SH2)*cos(180.*deg - SH_ANGLE);

    G4double VX4Al = VX1 +VX2 +VX3;
    G4double VY4Al = VY1P+VY2 +VY3P;
    G4double VZ4Al = VZ1 +VZ2P+VZ3P;

    G4Trap* solidShieldAl2 = new G4Trap("ShieldAl2", X_SHIEL, z_shiel_Al, L_SHIEL3, L_SHIEL3 - DY_SH2);
    G4LogicalVolume* logicShieldAl2 = new G4LogicalVolume(solidShieldAl2, Al, "World", 0, 0, 0);
    G4VPhysicalVolume* physiShieldAl2 = new G4PVPlacement(rm11, G4ThreeVector(VX4Al,VY4Al,VZ4Al), logicShieldAl2,
                                                          "ShieldAl2", World, false, 1);
    physiShieldAl2 = new G4PVPlacement(rm12, G4ThreeVector(VX4Al,VY4Al,-VZ4Al), logicShieldAl2, "ShieldAl2", World, false, 2);

    } else {
      X_SHIEL = 600.*mm;
      if (SiddhartaSetup == 8) {
        Y_SHIEL = 450.*mm;
      } else {
        Y_SHIEL = 450.*mm;
      }
      Z_SHIEL = 65.0*mm;
      G4double z_shiel_Al = 5.0*mm;
      pos_shield = 40.*mm + 0.5*Z_SHIEL;  // height of the lead table center

      if (SiddhartaSetup == 2020) {
        X_SHIEL = 500.*mm;
        Y_SHIEL = 515.*mm;
        Z_SHIEL = 60.0*mm;
        pos_shield = 44.*mm + 0.5*Z_SHIEL;  // height of the lead table center

        r1_hole = 0.5*52.4*mm;
        r2_hole = 0.5*79.4*mm;
        r1_hole_polyethylene = 0.5*56.5*mm;
        r2_hole_polyethylene = 0.5*83.5*mm;
        r1_hole_Al = 0.5*66.7*mm;
        r2_hole_Al = 0.5*93.7*mm;
        r1_hole_Cu = 0.5*77.*mm;
        r2_hole_Cu = 0.5*104.*mm;
      }
// Lead table
      G4Box* solidShielding1 = new G4Box("Shielding1",0.5*X_SHIEL,0.5*Y_SHIEL,0.5*Z_SHIEL);
      G4LogicalVolume* logicShielding1 = new G4LogicalVolume( solidShielding1, Pb, "World", 0, 0, 0);
      G4VPhysicalVolume* physiShielding1 = new G4PVPlacement(rm01,G4ThreeVector(0,pos_shield,0),logicShielding1,"Shielding1",World,false,1);

      if (SiddhartaSetup == 6 || SiddhartaSetup == 7 || SiddhartaSetup == 8) {
        r1_hole = 0.5*60.*mm;
        r2_hole = 0.5*80.*mm;
        r1_hole_polyethylene = 0.5*70.*mm;
        r2_hole_polyethylene = 0.5*90.*mm;
        r1_hole_Al = 0.5*80.*mm;
        r2_hole_Al = 0.5*100.*mm;
        r1_hole_Cu = 0.5*90.*mm;
        r2_hole_Cu = 0.5*110.*mm;
      } else {
        r1_hole = 0.5*53.0*mm;
        r2_hole = 0.5*78.0*mm;
        r1_hole_polyethylene = 0.5*57.0*mm;
        r2_hole_polyethylene = 0.5*82.0*mm;
        r1_hole_Al = 0.5*67.0*mm;
        r2_hole_Al = 0.5*92.0*mm;
        r1_hole_Cu = 0.5*77.0*mm;
        r2_hole_Cu = 0.5*102.0*mm;
      }
      G4Cons* solidWindowShield = new G4Cons("WindowShield",0., r1_hole, 0., r2_hole, 0.5*Z_SHIEL, 0., 360.*deg);
      G4LogicalVolume* logicWindowShield = new G4LogicalVolume(solidWindowShield, Air, "Shielding1", 0, 0, 0);
      G4VPhysicalVolume* physiWindowShield = new G4PVPlacement(rm00, G4ThreeVector(0.,0.,0.), logicWindowShield,
                                                               "WindowShield", logicShielding1, false, 1);

      G4Cons* solidWindowShield2 = new G4Cons("WindowShield2",r1_hole, r1_hole_polyethylene, r2_hole, r2_hole_polyethylene,
                                              0.5*Z_SHIEL, 0., 360.*deg);
      G4LogicalVolume* logicWindowShield2 = new G4LogicalVolume(solidWindowShield2, teflon, "Shielding1", 0, 0, 0);
      G4VPhysicalVolume* physiWindowShield2 = new G4PVPlacement(rm00, G4ThreeVector(0.,0.,0.), logicWindowShield2,
                                                                "WindowShield2", logicShielding1, false, 1);

      G4Cons* solidWindowShield3 = new G4Cons("WindowShield3", r1_hole_polyethylene, r1_hole_Al, r2_hole_polyethylene,
                                              r2_hole_Al, 0.5*Z_SHIEL, 0., 360.*deg);
      G4LogicalVolume* logicWindowShield3 = new G4LogicalVolume(solidWindowShield3, Al, "Shielding1", 0, 0, 0);
      G4VPhysicalVolume* physiWindowShield3 = new G4PVPlacement(rm00, G4ThreeVector(0.,0.,0.), logicWindowShield3,
                                                                "WindowShield3", logicShielding1, false, 1);

      G4Cons* solidWindowShield4 = new G4Cons("WindowShield4", r1_hole_Al, r1_hole_Cu, r2_hole_Al, r2_hole_Cu, 0.5*Z_SHIEL, 0., 360.*deg);
      G4LogicalVolume* logicWindowShield4 = new G4LogicalVolume(solidWindowShield4, Cu, "Shielding1", 0, 0, 0);
      G4VPhysicalVolume* physiWindowShield4 = new G4PVPlacement(rm00, G4ThreeVector(0.,0.,0.), logicWindowShield4,
                                                                "WindowShield4", logicShielding1, false, 1);

      G4double Z_SHIEL4;
      G4double Y_SHIEL4;

      if ( SiddhartaSetup == 8 ) {
        Z_SHIEL4 = 200*mm + 100*mm;
        Y_SHIEL4 = 200*mm + 80*mm;
      } else {
        Z_SHIEL4 = 200*mm;
        Y_SHIEL4 = 200*mm;
      }
      G4double X_SHIEL4 = X_SHIEL;
      G4double posz_shiel4;
      G4double posy_shiel4;

      if ( SiddhartaSetup == 8 ) {
        posz_shiel4 = 0.5*450*mm + 20*mm + 0.5*Z_SHIEL4;
        posy_shiel4 = 100*mm + 0.5*Y_SHIEL4 + 2*mm;
      } else {
        posz_shiel4 = 0.5*450*mm + 33*mm + 0.5*Z_SHIEL4;
        posy_shiel4 = 100*mm + 0.5*Y_SHIEL4 + 2*mm;
      }
      G4Box* solidShielding4 = new G4Box("Shielding4",0.5*X_SHIEL4,0.5*Y_SHIEL4,0.5*Z_SHIEL4);
      G4LogicalVolume* logicShielding4 = new G4LogicalVolume( solidShielding4,Pb, "World", 0, 0, 0);
      G4VPhysicalVolume* physiShielding4 = new G4PVPlacement(0,G4ThreeVector(0.,posy_shiel4,posz_shiel4),logicShielding4,"Shielding4",World,false,0);
      physiShielding4 = new G4PVPlacement(0,G4ThreeVector(0.,posy_shiel4,-posz_shiel4),logicShielding4,"Shielding4",World,false,1);

      G4double Y_SHIEL5 = 100*mm + 2.*mm - bp_r3;
      G4double Z_SHIEL5 = 0.5*450.*mm + 68.*mm -(0.5*450.*mm + bpj_L);
      G4double X_SHIEL5 = X_SHIEL;
      G4double posz_shiel5 = 0.5*450*mm + 68.*mm - 0.5*Z_SHIEL5;
      G4double posy_shiel5 = bp_r3 + 0.5*Y_SHIEL5;

      G4Box* solidShielding5 = new G4Box("Shielding5", 0.5*X_SHIEL5, 0.5*Y_SHIEL5, 0.5*Z_SHIEL5);
      G4LogicalVolume* logicShielding5 = new G4LogicalVolume( solidShielding5,Pb, "World", 0, 0, 0);

      G4double Y_SHIEL6 = 2*(100*mm + 2.*mm);
      G4double Z_SHIEL6 = bpj_L;
      G4double X_SHIEL6 = X_SHIEL;
      G4double posz_shiel6 = 0.5*450*mm + 0.5*Z_SHIEL6;
      G4double posy_shiel6 = 0;

      G4Box* solidShielding6 = new G4Box("Shielding6", 0.5*X_SHIEL6, 0.5*Y_SHIEL6, 0.5*Z_SHIEL6);
      G4LogicalVolume* logicShielding6 = new G4LogicalVolume(solidShielding6,Pb, "World", 0, 0, 0);
      G4SubtractionSolid* solidShielding7 = new G4SubtractionSolid("Shielding7", solidShielding6, solidBeamPipeJunction,
                                                                   rm00, G4ThreeVector(0.,0.,0));
      G4IntersectionSolid* solidShielding8 = new G4IntersectionSolid("Shielding8", solidShielding7, solidShielding6,
                                                                     rm00, G4ThreeVector(0.,0.5*Y_SHIEL6 + bp_r3,0.));
      G4LogicalVolume* logicShielding8 = new G4LogicalVolume(solidShielding8,Pb, "World", 0, 0, 0);
    }

// KPlus Detector //
  G4double kplus_x;
  G4double kplus_y;
  G4double kplus_z;
  G4double kplus_teflon_z;
  G4double pos_kplus;
  G4double pos_kplus_teflon;

  if (SiddhartaSetup == 2) {
    kplus_x = 160.0*mm;
    kplus_y = 80.0*mm;
    kplus_z = 10.0*mm;
    pos_kplus = mycard->variables["z_kplus"];
  } else if (SiddhartaSetup == 4 || SiddhartaSetup == 5 || SiddhartaSetup == 6 || SiddhartaSetup == 7 || SiddhartaSetup == 8) {
    kplus_x = 160.0*mm;
    kplus_y = 160.0*mm;
    kplus_z = 10.0*mm;
    pos_kplus = -67.0*mm;
  } else if (SiddhartaSetup == 2020) {
    kplus_x = 200.0*mm;
    kplus_y = 200.0*mm;
    kplus_z = 5.0*mm;
    pos_kplus = -135.0*mm - 3.*mm - kplus_z/2;
    kplus_teflon_z = 3.*mm;
    pos_kplus_teflon = -135.0*mm - (kplus_teflon_z/2.)*mm;
  } else {
    kplus_x = 160.0*mm;
    kplus_y = 80.0*mm;
    kplus_z = 10.0*mm;
    pos_kplus = -67.0*mm;
  }

  G4Box *solidKPlus = new G4Box("KPlus",0.5*kplus_x,0.5*kplus_y,0.5*kplus_z);
  G4LogicalVolume* logicKPlus;

  if (SiddhartaSetup == 2 || SiddhartaSetup == 4 || SiddhartaSetup == 5 || SiddhartaSetup == 2020) {
    logicKPlus = new G4LogicalVolume(solidKPlus, BC420, "World", 0, 0, 0);
  } else {
    logicKPlus = new G4LogicalVolume(solidKPlus, Air, "World", 0, 0, 0);
  }
  G4VPhysicalVolume* physiKPlus  = new G4PVPlacement(rm01,G4ThreeVector(0,pos_kplus,0),logicKPlus,"KPlus",World,false,0);
  logicKPlus->SetSensitiveDetector(KPlus_SD);
  logicKPlus->SetVisAttributes(G4Colour(0.,0.,1.));

//teflon for KPlus
  G4Box *solidKPlusTeflon = new G4Box("KPlus",0.5*kplus_x,0.5*kplus_y,0.5*kplus_teflon_z);
  G4LogicalVolume* logicKPlusTeflon = new G4LogicalVolume(solidKPlusTeflon, teflon, "World", 0, 0, 0);
  G4VPhysicalVolume* physiKPlusTeflon  = new G4PVPlacement(rm01,G4ThreeVector(0,pos_kplus_teflon,0),logicKPlusTeflon,"KPlusTeflon",World,false,0);
  logicKPlusTeflon->SetVisAttributes(G4Colour(0.,0.,0.,1.));

// Anticoincidence 1 //
  G4double w_scint_thin = 10.0*mm;
  G4double w_scint_thick = 10.0*mm;
  G4double half_box = 0.5*450.*mm + 0.5*w_scint_thick; // 230 mm

  if (SiddhartaSetup == 2 || SiddhartaSetup == 6) {
    G4double anti1_x = 0.5*(2.*half_box - w_scint_thick - 120.*mm);
    G4double anti1_y = 120.*mm;
    G4double anti1_z = w_scint_thin;
    G4double pos_anti1 = pos_shield + 0.5*Z_SHIEL + 0.5*w_scint_thin + 1*mm;
    G4double posx_anti1 = half_box - 0.5*anti1_x - 0.5*w_scint_thick;

    G4Box *solidAnti1 = new G4Box("Anti1", 0.5*anti1_x, 0.5*anti1_y, 0.5*anti1_z);
    G4LogicalVolume* logicAnti1 = new G4LogicalVolume(solidAnti1, BC420, "World", 0, 0, 0);
    G4VPhysicalVolume* physiAnti1 = new G4PVPlacement(rm01, G4ThreeVector(posx_anti1,pos_anti1,0), logicAnti1,
                                                      "Anti1", World, false, 0);
    physiAnti1 = new G4PVPlacement(rm01, G4ThreeVector(-posx_anti1,pos_anti1,0), logicAnti1, "Anti1", World, false, 1);
    logicAnti1->SetSensitiveDetector(Anti1_SD);
    logicAnti1->SetVisAttributes(G4Colour(1.,1.,0.));

// Anticoincidence 2 //
    G4double anti2_x = 2.*half_box - w_scint_thick; // 2*230 - 10 = 450
    G4double anti2_y = 0.5*(2.*half_box - w_scint_thick - 120.*mm);
    G4double anti2_z = w_scint_thin;

    G4double pos_anti2 = pos_shield + 0.5*Z_SHIEL + 0.5*w_scint_thin + 1.*mm;
    G4double posz_anti2 = half_box - 0.5*anti2_y - 0.5*w_scint_thick;

    G4Box *solidAnti2 = new G4Box("Anti2", 0.5*anti2_x, 0.5*anti2_y, 0.5*anti2_z);
    G4LogicalVolume* logicAnti2 = new G4LogicalVolume(solidAnti2, BC420, "World", 0, 0, 0);
    G4VPhysicalVolume* physiAnti2 = new G4PVPlacement(rm01, G4ThreeVector(0,pos_anti2,posz_anti2), logicAnti2,
                                                      "Anti2", World, false, 0);
    physiAnti2 = new G4PVPlacement(rm01, G4ThreeVector(0,pos_anti2,-posz_anti2), logicAnti2, "Anti2", World, false, 1);
    logicAnti2->SetSensitiveDetector(Anti2_SD);
    logicAnti2->SetVisAttributes(G4Colour(0.,1.,0.));

// Anticoincidence 3 //
    G4double yoff = 21.*mm;
    G4double anti3_x = 2*half_box + w_scint_thick;
    G4double anti3_y = 360.0*mm - yoff;
    G4double anti3_z = w_scint_thick;
    G4double posx_anti3 = 0.;
    G4double posy_anti3 = pos_shield + 0.5*Z_SHIEL + 0.5*anti3_y + 1.*mm + yoff;
    G4double posz_anti3 = half_box;

    G4RotationMatrix* rm36 = new G4RotationMatrix();
    G4RotationMatrix* rm37 = new G4RotationMatrix();

    G4Box *solidAnti3 = new G4Box("Anti3", 0.5*anti3_x, 0.5*anti3_y, 0.5*anti3_z);
    G4LogicalVolume* logicAnti3 = new G4LogicalVolume(solidAnti3, BC420, "World", 0, 0, 0);
    G4VPhysicalVolume* physiAnti3 = new G4PVPlacement(rm36, G4ThreeVector(posx_anti3,posy_anti3,posz_anti3), logicAnti3,
                                                      "Anti3", World, false, 0);
    physiAnti3 = new G4PVPlacement(rm37, G4ThreeVector(posx_anti3,posy_anti3,-posz_anti3), logicAnti3, "Anti3", World, false, 1);
    logicAnti3->SetSensitiveDetector(Anti3_SD);
    logicAnti3->SetVisAttributes(G4Colour(1.,0.,0.));

// Anticoincidence 4 //
    G4double anti4_x1 = 2*half_box - w_scint_thick;
    G4double anti4_x2 = 2*half_box - w_scint_thick;
    G4double anti4_y1 = w_scint_thick;
    G4double anti4_y2 = w_scint_thick;
    G4double anti4_z = 360.0*mm;
    G4double posx_anti4 = half_box;
    G4double posy_anti4 = pos_shield + 0.5*Z_SHIEL + 0.5*anti4_z + 1.*mm;
    G4double posz_anti4 = 0;

    G4Trd* solidAnti4 = new G4Trd("Anti4", 0.5*anti4_x1, 0.5*anti4_x2, 0.5*anti4_y1, 0.5*anti4_y2, 0.5*anti4_z);
    G4LogicalVolume* logicAnti4 = new G4LogicalVolume(solidAnti4, BC420, "World", 0, 0, 0);
    G4VPhysicalVolume* physiAnti4 = new G4PVPlacement(rm38, G4ThreeVector(posx_anti4,posy_anti4,posz_anti4), logicAnti4,
                                                      "Anti4", World, false, 0);
    physiAnti4 = new G4PVPlacement(rm38, G4ThreeVector(-posx_anti4,posy_anti4,posz_anti4), logicAnti4, "Anti4", World, false, 1);
    logicAnti4->SetSensitiveDetector(Anti4_SD);
    logicAnti4->SetVisAttributes(G4Colour(0.,0.,1.));

// Anticoincidence 5 //
    G4double anti5_x = 2.*half_box + w_scint_thick;
    G4double anti5_y = w_scint_thick;
    G4double anti5_z = 2.*half_box + w_scint_thick;
    G4double posx_anti5 = 0.*mm;
    G4double posy_anti5 = pos_shield + 0.5*Z_SHIEL + 360.*mm + 0.5*w_scint_thick + 1.*mm - 90*mm;
    G4double posz_anti5 = 0;

    G4Box* solidAnti5 = new G4Box("Anti5", 0.5*anti5_x, 0.5*anti5_y, 0.5*anti5_z);
    G4LogicalVolume* logicAnti5 = new G4LogicalVolume(solidAnti5, vacuum, "World", 0, 0, 0);//BC420
    G4VPhysicalVolume* physiAnti5 = new G4PVPlacement(rm00, G4ThreeVector(posx_anti5,posy_anti5,posz_anti5), logicAnti5,
                                                      "Anti5", World, false, 0);
    logicAnti5->SetSensitiveDetector(Anti5_SD);
    logicAnti5->SetVisAttributes(G4Colour(0.,1.,0.));
  }
  G4double dist_shield_degrader = 1*mm;
  G4double dz_degrader = 0.84*mm;
  G4double pos_degrader = pos_shield + 0.5*Z_SHIEL + dist_shield_degrader + 0.5*dz_degrader;
  G4cout << "pos_degrader//////////////////////////////////////////////////////////" << pos_degrader << "\n";

  G4double dist_deg_vac;
  if (SiddhartaSetup == 2) {
    dist_deg_vac = 22*mm;
  } else {
    dist_deg_vac = 20*mm;
  }
  G4double dz_vac_chamb = 26.*cm;
  G4double thikness_vac_chamb = 10.0*mm;
  G4double dist_vac_target = 10.0*mm;
  G4double pos_vac_chamb;

// Vacuum Chamber //
  G4double r1_vac_chamb;
  G4double r2_vac_chamb;
  G4double r3_vac_chamb;
  G4double r4_vac_chamb;

  if (SiddhartaSetup == 2) {
    thikness_vac_chamb = 10.*mm;
    r1_vac_chamb = 0.5*450.*mm - thikness_vac_chamb - 2.*mm;
    r2_vac_chamb = r1_vac_chamb + thikness_vac_chamb ;
    r3_vac_chamb = 0.5*450.*mm - thikness_vac_chamb;
    r4_vac_chamb = r3_vac_chamb + thikness_vac_chamb;
    dz_vac_chamb = 260.*mm;
    pos_vac_chamb = 120.*mm + thikness_vac_chamb + 0.5*dz_vac_chamb;
  } else if (SiddhartaSetup == 3) {
    r1_vac_chamb = 110.*mm;
    thikness_vac_chamb = 10.*mm;
    r2_vac_chamb = r1_vac_chamb + thikness_vac_chamb;
    r3_vac_chamb = 190.*mm;
    r4_vac_chamb = r3_vac_chamb + thikness_vac_chamb;
    dz_vac_chamb = 260.*mm;
    pos_vac_chamb = 105.*mm + thikness_vac_chamb + 0.5*dz_vac_chamb;
  } else if (SiddhartaSetup == 4 || SiddhartaSetup == 5 || SiddhartaSetup == 6 || SiddhartaSetup == 7) {
    thikness_vac_chamb = 10.*mm;
    r1_vac_chamb = 0.5*450.*mm - thikness_vac_chamb;
    r2_vac_chamb = r1_vac_chamb + thikness_vac_chamb ;
    r3_vac_chamb = 0.5*450.*mm - thikness_vac_chamb;
    r4_vac_chamb = r3_vac_chamb + thikness_vac_chamb;
    dz_vac_chamb = 260.*mm;
    pos_vac_chamb = 120.*mm + thikness_vac_chamb + 0.5*dz_vac_chamb;
  } else if (SiddhartaSetup == 8) {
    thikness_vac_chamb = 10.*mm;
    r1_vac_chamb = 0.5*420.*mm - thikness_vac_chamb;
    r2_vac_chamb = r1_vac_chamb + thikness_vac_chamb ;
    r3_vac_chamb = 0.5*420.*mm - thikness_vac_chamb;
    r4_vac_chamb = r3_vac_chamb + thikness_vac_chamb;
    dz_vac_chamb = 260.*mm;
    pos_vac_chamb = 120.*mm + thikness_vac_chamb + 0.5*dz_vac_chamb;
  } else if (SiddhartaSetup == 2020) {
    thikness_vac_chamb = 10.*mm;
    r1_vac_chamb = 0.5*424.*mm;
    r2_vac_chamb = r1_vac_chamb + thikness_vac_chamb ;
    r3_vac_chamb = 0.5*424.*mm;
    r4_vac_chamb = r3_vac_chamb + thikness_vac_chamb;
    dz_vac_chamb = 260.*mm;
    pos_vac_chamb = 135.*mm + thikness_vac_chamb + 0.5*dz_vac_chamb;
   // pos_vac_chamb = 143.5*mm + thikness_vac_chamb + 0.5*dz_vac_chamb;
  } else {
    r1_vac_chamb = 110.*mm;
    thikness_vac_chamb = 10.*mm;
    r2_vac_chamb = r1_vac_chamb + thikness_vac_chamb;
    r3_vac_chamb = 190.*mm;
    r4_vac_chamb = r3_vac_chamb + thikness_vac_chamb;
    dz_vac_chamb = 260.*mm;
    pos_vac_chamb = 150.*mm + thikness_vac_chamb + 0.5*dz_vac_chamb;
  }
  G4Cons* solidVacCh = new G4Cons("VacCh", 0, r2_vac_chamb, 0, r4_vac_chamb, 0.5*dz_vac_chamb, 0., 360.*deg);
  G4LogicalVolume* logicVacCh = new G4LogicalVolume(solidVacCh,vacuum, "World", 0, 0, 0);
  G4VPhysicalVolume* physiVacCh = new G4PVPlacement(rm01, G4ThreeVector(0,pos_vac_chamb,0), logicVacCh, "VacCh", World, false, 0);
  logicVacCh->SetVisAttributes(G4Colour(0.,1.,0.,0.));

  G4Cons* solidVacChCon = new G4Cons("VacChCon", r1_vac_chamb, r2_vac_chamb, r3_vac_chamb, r4_vac_chamb, 0.5*dz_vac_chamb, 0., 360.*deg);
  G4LogicalVolume* logicVacChCon = new G4LogicalVolume(solidVacChCon, Al, "VacCh", 0, 0, 0);
  G4VPhysicalVolume* physiVacChCon  = new G4PVPlacement(rm00, G4ThreeVector(0,0,0), logicVacChCon, "VacChCon", logicVacCh, false, 0);
  logicVacChCon->SetVisAttributes(G4Colour(1.,0.,0.,0.));

  G4double r1_entrance;
  G4double r2_entrance;

  if (SiddhartaSetup == 2) {
    r1_entrance = 0.5*170.0*mm;
    r2_entrance = r1_entrance + 4.*mm;
  } else if (SiddhartaSetup == 4 || SiddhartaSetup == 5 || SiddhartaSetup == 6 || SiddhartaSetup == 7 || SiddhartaSetup == 8) {
    r1_entrance = 0.5*220.0*mm;
    r2_entrance = r1_entrance + 4.*mm;
  } else if ( SiddhartaSetup == 2020 ) {
    r1_entrance = 58.*mm;
    r2_entrance = r1_entrance + 4.*mm;
  } else {
    r1_entrance = 70.0*mm;
    r2_entrance = 80.0*mm;
  }
  G4double pos_VacChBase = pos_vac_chamb - 0.5*dz_vac_chamb - 0.5*thikness_vac_chamb;
  G4Tubs* solidVacChBase = new G4Tubs("VacChBase", r1_entrance, r2_vac_chamb, 0.5*thikness_vac_chamb, 0., 360.*deg);
  G4LogicalVolume* logicVacChBase = new G4LogicalVolume(solidVacChBase, Al, "World", 0, 0, 0);
  G4VPhysicalVolume* physiVacChBase = new G4PVPlacement(rm01, G4ThreeVector(0,pos_VacChBase,0), logicVacChBase, "VacChBase", World, false, 0);
  logicVacChBase->SetVisAttributes(G4Colour(1.,0.,0.,0.));

  G4double teflon_cone_w = 2.0*mm;
  G4double cryo_ent_kapton_w = 125.0*um;
  G4double teflon_cyl_w = 1.0*mm;
  G4double cryo_cyl_teflon_z = 0.5*thikness_vac_chamb - 0.5*cryo_ent_kapton_w;
  G4double thikness_entrance = 120.*um;
  G4double dz_cyl_entrance;

  if (SiddhartaSetup == 1) {
    dz_cyl_entrance = 10.*mm;
  } else {
    dz_cyl_entrance = 1.*mm;
  }

// Zr foil outside vacuum chamber entrance window
  G4double dz_zr_foil;
  dz_zr_foil = 0.1*um;

  G4Cons* solidVacChTef = new G4Cons("VacChTef", r1_vac_chamb - teflon_cone_w, r1_vac_chamb, r3_vac_chamb-teflon_cone_w,
                                     r3_vac_chamb, 0.5*dz_vac_chamb, 0., 360.*deg);
  G4LogicalVolume* logicVacChTef = new G4LogicalVolume(solidVacChTef, teflon, "VacCh", 0, 0, 0);
  G4VPhysicalVolume* physiVacChTef = new G4PVPlacement(rm00, G4ThreeVector(0,0,0), logicVacChTef, "VacChTef", logicVacCh, false, 0);
  logicVacChTef->SetVisAttributes(G4Colour(1.,0.,1.,1.));

  G4Tubs* solidVacChTefEntrance = new G4Tubs("VacChTefEntrance", r1_entrance, r1_vac_chamb - teflon_cone_w, 0.5*teflon_cone_w, 0., 360.*deg);
  G4LogicalVolume* logicVacChTefEntrance = new G4LogicalVolume(solidVacChTefEntrance, teflon, "VacCh", 0, 0, 0);
  G4VPhysicalVolume* physiVacChTefEntrance = new G4PVPlacement(rm00, G4ThreeVector(0,0,-0.5*dz_vac_chamb + 0.5*teflon_cone_w),
                                                               logicVacChTefEntrance, "VacChTefEntrance", logicVacCh, false, 0);
  logicVacChTefEntrance->SetVisAttributes(G4Colour(1.,0.,1.,1.));

  G4Tubs* solidCryoEntrance = new G4Tubs("CryoEntrance", 0, r1_entrance, 0.5*cryo_ent_kapton_w, 0., 360.*deg);
  G4LogicalVolume* logicCryoEntrance = new G4LogicalVolume(solidCryoEntrance, kapton, "World", 0, 0, 0);
  G4VPhysicalVolume* physiCryoEntrance = new G4PVPlacement(rm01, G4ThreeVector(0,pos_vac_chamb - 0.5*dz_vac_chamb - 0.5*thikness_vac_chamb,0),
                                                           logicCryoEntrance, "CryoEntrance", World, false, 0);
  logicCryoEntrance->SetVisAttributes(G4Colour(1.,1.,0.,1.));

  G4Tubs* solidCryoEntranceTeflon = new G4Tubs("CryoEntranceTeflon", r1_entrance - teflon_cyl_w, r1_entrance, 0.5*cryo_cyl_teflon_z, 0., 360.*deg);
  G4LogicalVolume* logicCryoEntranceTeflon = new G4LogicalVolume(solidCryoEntranceTeflon, teflon, "World", 0, 0, 0);
  G4VPhysicalVolume* physiCryoEntranceTeflon =
    new G4PVPlacement(rm01, G4ThreeVector(0,pos_vac_chamb - 0.5*dz_vac_chamb - 0.5*thikness_vac_chamb + 0.5*cryo_ent_kapton_w + 0.5*cryo_cyl_teflon_z,0),
                      logicCryoEntranceTeflon, "CryoEntranceTeflon", World, false, 0);
  logicCryoEntranceTeflon->SetVisAttributes(G4Colour(0.,1.,1.,1.));

  G4Tubs* solidCryoEntrance3 = new G4Tubs("CryoEntrance3", r1_entrance - thikness_entrance, r1_entrance, 0.5*dz_cyl_entrance, 0., 360.*deg);
  G4LogicalVolume* logicCryoEntrance3 = new G4LogicalVolume(solidCryoEntrance3,black_paper, "World", 0, 0, 0);
  G4VPhysicalVolume* physiCryoEntrance3 =
    new G4PVPlacement(rm01, G4ThreeVector(0,pos_vac_chamb - 0.5*dz_vac_chamb - 0.5*thikness_vac_chamb - 0.5*cryo_ent_kapton_w - 0.5*dz_cyl_entrance,0),
                      logicCryoEntrance3, "CryoEntrance3", World, false, 0);
  logicCryoEntrance3->SetVisAttributes(G4Colour(1.,0.,0.,1.));

  G4Tubs* solidCryoEntrance4 = new G4Tubs("CryoEntrance4", 0, r1_entrance, 0.5*thikness_entrance, 0., 360.*deg);
  G4LogicalVolume* logicCryoEntrance4 = new G4LogicalVolume(solidCryoEntrance4, black_paper, "World", 0, 0, 0);
  G4VPhysicalVolume* physiCryoEntrance4 =
    new G4PVPlacement(rm01, G4ThreeVector(0,
                                          pos_vac_chamb - 0.5*dz_vac_chamb - 0.5*thikness_vac_chamb -
                                          0.5*cryo_ent_kapton_w - 2*0.5*dz_cyl_entrance - 0.5*thikness_entrance,0),
                      logicCryoEntrance4, "CryoEntrance4", World, false, 0);
  logicCryoEntrance4->SetVisAttributes(G4Colour(1.,0.,0.,1.));

// Zr foil for test of the kaon timing PIXE lines
 /* G4Tubs* solidZrFoil = new G4Tubs("ZrFoilDeg", 0, r1_entrance, 0.5*dz_zr_foil, 0., 360.*deg);
  G4LogicalVolume* logicZrFoil = new G4LogicalVolume(solidZrFoil, Zr, "ZrFoil", 0, 0, 0);
  G4VPhysicalVolume* physiZrFoil = new G4PVPlacement(rm01, G4ThreeVector(0,pos_vac_chamb - 0.5*dz_vac_chamb - thikness_vac_chamb - 2*mm,0),
                                                     logicZrFoil, "ZrFoilDeg", World, false, 0);
  logicZrFoil->SetVisAttributes(G4Colour(1.,0.,0.,1.));
*/
// Target //
  G4double target_cyl_w = 75.0*um;
  G4double target_entrance_w = 50.0*um;
  G4double tgt_kap_ent_r  = 0.5 * 116 * mm;
  G4double r_target;
  G4double dz_target;

  if ((SiddhartaSetup == 2 || SiddhartaSetup == 5 || SiddhartaSetup == 6 || SiddhartaSetup == 7 || SiddhartaSetup == 8 ) && !OldTarget) {
    dz_target = 140.0*mm;
    r_target = 0.5*130.0*mm;
  } else if ( SiddhartaSetup == 2020 ) {
    dz_target = 125.0*mm;
    r_target = 0.5*144.0*mm;
    target_entrance_w = 125.*um;
    target_cyl_w = 150.0*um;
  } else {
    dz_target = 155.0*mm;
    r_target = 0.5*137.0*mm;
  }
  G4double pos_target;

  if (SiddhartaSetup == 2 || SiddhartaSetup == 4 || SiddhartaSetup == 5 || SiddhartaSetup == 6 || SiddhartaSetup == 7 || SiddhartaSetup == 8) {
    pos_target = -pos_vac_chamb + 150.*mm + target_entrance_w + 0.5*dz_target;
  } else if (SiddhartaSetup == 2020) {
    pos_target = -pos_vac_chamb + 182.*mm + target_entrance_w + 0.5*dz_target;
 //   pos_target = -pos_vac_chamb + 179.5*mm + target_entrance_w + 0.5*dz_target;
  } else if (SiddhartaSetup == 3) {
    pos_target = -pos_vac_chamb + 125.*mm + target_entrance_w + 0.5*dz_target;
  } else {
    pos_target = -pos_vac_chamb + 198.*mm + target_entrance_w + 0.5*dz_target;
  }

// 2020 Gas cell Ti top exit window
  G4double TiFoil_w = 100*um;
  G4double tgt_ti_exit_win_r = 0.5*110*mm;
// 2021 Ti Foil before target entrance for test
  G4double TiFoil_ent = 0.1*um;
  G4double tgt_ti_ent_win_r = 0.5*100*mm;
  G4Tubs* solidTarget = new G4Tubs("Target",0,r_target,0.5*(dz_target - TiFoil_w - 2*mm),0.,360.*deg);

// Definition of target density gradient elements
  G4double targetGradientStepSize = 0.5*(dz_target - TiFoil_w - 2*mm);
  targetGradientStepSize /= densityGradientNumberOfSteps;
  G4Tubs* solidTargetGradientStep = new G4Tubs("Target",0,r_target,targetGradientStepSize,0.,360.*deg);
// End definition

  if (SiddhartaSetup == 6 || SiddhartaSetup == 7 || SiddhartaSetup == 8) {
    G4LogicalVolume* logicTarget = new G4LogicalVolume(solidTarget, Hydrogen, "VacCh", 0, 0, 0);
    G4VPhysicalVolume* physiTarget = new G4PVPlacement(rm00, G4ThreeVector(0,0,pos_target - 0.5*TiFoil_w - 1*mm), logicTarget,
                                                        "Target", logicVacCh, false, 0);
    logicTarget->SetVisAttributes(G4Colour(1.,0.,0.));
  } else if (SiddhartaSetup == 2020) {
    G4LogicalVolume* logicTarget;
    if (targetMaterialFromCard) {
      G4cout << "------Target------" << G4endl;
      G4cout << "|                |" << G4endl;
      G4cout << "|                |" << G4endl;
      if (targetMaterialDensityGradient || targetMaterialCustomDensityGradient) {
        G4cout << "|     vacuum     |" << G4endl;
        logicTarget = new G4LogicalVolume(solidTarget, vacuum, "tgt_L", 0, 0, 0);
      } else if (atomicNumberOfTargetFromCard == 1) {
        G4cout << "|    Hydrogen    |" << G4endl;
        logicTarget = new G4LogicalVolume(solidTarget, Hydrogen, "tgt_L", 0, 0, 0);
      } else if (atomicNumberOfTargetFromCard == 2) {
        G4cout << "|    deuterium   |" << G4endl;
        logicTarget = new G4LogicalVolume(solidTarget, deuterium, "tgt_L", 0, 0, 0);
      } else if (atomicNumberOfTargetFromCard == 4) {
        G4cout << "|     Helium     |" << G4endl;
        logicTarget = new G4LogicalVolume(solidTarget, Helium4, "tgt_L", 0, 0, 0);
      } else if (atomicNumberOfTargetFromCard == 20 || atomicNumberOfTargetFromCard == 22) {
        G4cout << "|     Neonium    |" << G4endl;
        logicTarget = new G4LogicalVolume(solidTarget, Neonium, "tgt_L", 0, 0, 0);
      } else {
        G4cout << "|     vacuum     |" << G4endl;
        logicTarget = new G4LogicalVolume(solidTarget, vacuum, "tgt_L", 0, 0, 0);
      }
      G4cout << "|                |" << G4endl;
      G4cout << "|                |" << G4endl;
      G4cout << "------------------" << G4endl;
    } else {
      logicTarget = new G4LogicalVolume(solidTarget, Neonium, "tgt_L", 0, 0, 0);
    }
    G4VPhysicalVolume* physiTarget = new G4PVPlacement(rm00, G4ThreeVector(0,0,pos_target - 0.5*TiFoil_w - 1*mm), logicTarget,
                                                        "Target", logicVacCh, false, 0);
    logicTarget->SetVisAttributes(G4Colour(1.,0.,0.));

    if (targetMaterialDensityGradient || targetMaterialCustomDensityGradient) {
      std::vector<G4LogicalVolume*> logSteps;
// gradMaterials previously defined materials for step in the material section
      std::cout << std::endl;
      std::cout << "__--^^--__--^^--__" << std::endl;
      std::cout << "Creating density gradient in the target" << std::endl;
      G4cout << "------Target------" << G4endl;
      for (int i=0; i<densityGradientNumberOfSteps; i++) {
// gradient step defined inside the target logical so position of the target not needed
        G4double posZStep = 2*(- i*targetGradientStepSize + (densityGradientNumberOfSteps - 1)*0.5*targetGradientStepSize);
        if (i<10) {
          G4cout << "|  " << gradMaterials.at(i)->GetName() << "  | density / divided by liquid target density [%] ";
          G4cout << gradMaterials.at(i)->GetDensity()/liquidTargetDensityForPrinting*100;
          G4cout << " with PosZ in target " << posZStep << G4endl;
        } else if (i < 100) {
          G4cout << "|  " << gradMaterials.at(i)->GetName() << " | density / divided by liquid target density [%] ";
          G4cout << gradMaterials.at(i)->GetDensity()/liquidTargetDensityForPrinting*100;
          G4cout << " with PosZ in target " << posZStep << G4endl;
        } else {
          G4cout << "| " << gradMaterials.at(i)->GetName() << " | density / divided by liquid target density [%] ";
          G4cout << gradMaterials.at(i)->GetDensity()/liquidTargetDensityForPrinting*100;
          G4cout << " with PosZ in target " << posZStep << G4endl;
        }
        G4LogicalVolume* tempLog = new G4LogicalVolume(solidTargetGradientStep, gradMaterials.at(i), "tgt_L", 0, 0, 0);
	    new G4PVPlacement(rm00, G4ThreeVector(0,0,posZStep), tempLog, "Target", logicTarget, false, i);
	    tempLog->SetVisAttributes(G4Colour(1.,0.,0.));
      }
      G4cout << "------------------" << G4endl;
      std::cout << "__--^^--__--^^--__" << std::endl;
      std::cout << std::endl;
    }
  } else {
    G4LogicalVolume* logicTarget = new G4LogicalVolume(solidTarget, Hydrogen, "VacCh", 0, 0, 0);
    G4VPhysicalVolume* physiTarget = new G4PVPlacement(rm00, G4ThreeVector(0,0,pos_target), logicTarget, "Target", logicVacCh, false, 0);
    logicTarget->SetVisAttributes(G4Colour(1.,0.,0.));
  }
  G4Tubs* solidTargetMaterialBase = new G4Tubs("TargetMaterialBase", 0., tgt_kap_ent_r, 0.5*target_entrance_w, 0., 360.*deg);
  G4LogicalVolume* logicTargetMaterialBase = new G4LogicalVolume(solidTargetMaterialBase, kapton, "KapW_L", 0, 0, 0);
  G4VPhysicalVolume* physiTargetMaterialBase = new G4PVPlacement(rm00, G4ThreeVector(0,0,pos_target - 0.5*dz_target - 0.5*target_entrance_w),
                                                                 logicTargetMaterialBase, "TargetMaterialBase", logicVacCh, false, 0);
  logicTargetMaterialBase->SetVisAttributes(G4Colour(0.,1.,0.));

  G4Tubs* solidTargetMaterial = new G4Tubs("TargetMaterial", r_target, r_target + target_cyl_w, 0.5*dz_target - 1*mm, 0., 360.*deg);
  G4LogicalVolume* logicTargetMaterial = new G4LogicalVolume(solidTargetMaterial, kapton, "VacCh", 0, 0, 0);
  G4VPhysicalVolume* physiTargetMaterial = new G4PVPlacement(rm00, G4ThreeVector(0,0,pos_target - 1*mm), logicTargetMaterial,
                                                             "TargetMaterial", logicVacCh, false, 0);
  logicTargetMaterial->SetVisAttributes(G4Colour(0.,1.,0.,0.));

// additional foil on the Kapton walls
    G4double target_cyl_w_PETfoil = 12*um;

    G4Tubs* solidTargetMaterialPETfoil = new G4Tubs("TargetMaterialPETfoil", r_target + target_cyl_w, r_target + target_cyl_w + target_cyl_w_PETfoil, 
						    0.5*dz_target - 1*mm, 0., 360.*deg);
    G4LogicalVolume* logicTargetMaterialPETfoil = new G4LogicalVolume(solidTargetMaterialPETfoil, mylar, "VacCh", 0, 0, 0);
    G4VPhysicalVolume* physiTargetMaterialPETfoil = new G4PVPlacement(rm00, G4ThreeVector(0,0,pos_target - 1*mm), logicTargetMaterialPETfoil,
                                                             "TargetMaterial_PetFoil", logicVacCh, false, 0);
    logicTargetMaterialPETfoil->SetVisAttributes(G4Colour(1.,0.8,0.3,1.));

    G4double radius_al_bar = 2.5*mm;
    G4double dz_al_bar = 105.*mm;
    G4double radius_al_bar_PetFoil = 280.*um;

    G4double radius_distance_al_bar = 82*mm;
    G4int nmb_of_bars = 4;
    G4double angle_shift_al_bar = twopi/24.;
    G4double ang_increment = twopi/nmb_of_bars;

    G4Tubs* solidTargetMaterialAlBar = new G4Tubs("TargetMaterialAlBar", 0, radius_al_bar, 0.5*dz_al_bar, 0., 360.*deg);
    G4LogicalVolume* logicTargetMaterialAlBar = new G4LogicalVolume(solidTargetMaterialAlBar, Al, "AlBar", 0, 0, 0);
    logicTargetMaterialAlBar->SetVisAttributes(G4Colour(1.,0.8,0.3,1.));

    G4Tubs* solidTargetMaterialAlBarPetFoil = new G4Tubs("TargetMaterialPetFoil", radius_al_bar, radius_al_bar + radius_al_bar_PetFoil, 0.5*dz_al_bar, 0., 360.*deg);
    G4LogicalVolume* logicTargetMaterialAlBarPetFoil = new G4LogicalVolume(solidTargetMaterialAlBarPetFoil, mylar, "PetFoilAlBar", 0, 0, 0);
    logicTargetMaterialAlBarPetFoil->SetVisAttributes(G4Colour(1.,1.,0.,1.));

    for (unsigned i=0; i<nmb_of_bars; i++) {
      G4double angleBar = angle_shift_al_bar + ang_increment*i;
      G4double X_Al_Bar = sin(angleBar)*radius_distance_al_bar;
      G4double Y_Al_Bar = cos(angleBar)*radius_distance_al_bar;
      G4double Z_Al_Bar = pos_target;
      G4VPhysicalVolume* physiTargetMaterialAlBar = new G4PVPlacement(rm00, G4ThreeVector(X_Al_Bar,Y_Al_Bar,Z_Al_Bar), logicTargetMaterialAlBar,
                                                             "TargetMaterial_AlBar", logicVacCh, false, 0);
      G4VPhysicalVolume* physiTargetMaterialAlBarPetFoil = new G4PVPlacement(rm00, G4ThreeVector(X_Al_Bar,Y_Al_Bar,Z_Al_Bar), logicTargetMaterialAlBarPetFoil,
                                                             "TargetMaterial_AlBarPetFoil", logicVacCh, false, 0);
    }

// Ti Foil exit window
  G4Tubs* solidTitaniumFoil = new G4Tubs("TitaniumFoil", 0, tgt_ti_exit_win_r, 0.5*TiFoil_w, 0., 360.*deg);
  G4LogicalVolume* logicTitaniumFoil = new G4LogicalVolume(solidTitaniumFoil, Ti, "Titanium", 0, 0, 0);
  G4VPhysicalVolume* physiTitaniumFoil = new G4PVPlacement(rm00, G4ThreeVector(0.,0.,pos_target + 0.5*dz_target - 0.5*TiFoil_w),
                                                           logicTitaniumFoil, "TitaniumFoil", logicVacCh, false, 0);
  logicTitaniumFoil->SetVisAttributes(G4Colour(1.,0.8,0.3,1.));

////////////////////////////////////////
//// Alumminum holder structure 2020
////////////////////////////////////////
  G4double tgt_al_top_out_r_max = 0.5*175*mm;
  G4double tgt_al_top_out_r_min = tgt_ti_exit_win_r;
  G4double tgt_al_top_out_w = 12*mm ;

  G4double tgt_al_bot_out_r_max = 0.5*175*mm;
  G4double tgt_al_bot_out_r_min = 0.5*150*mm;
  G4double tgt_al_bot_out_w = 20*mm;

  G4double tgt_al_bot_in_r_max = tgt_al_bot_out_r_min;
  G4double tgt_al_bot_in_r_min = tgt_kap_ent_r;
  G4double tgt_al_bot_in_w = 10*mm ;

  G4double tgt_al_bot_out_shift = 0;
  G4double tgt_al_bot_in_shift = 0;

    tgt_al_bot_out_w = 8*mm;
    tgt_al_bot_in_w = 15*mm;
    tgt_al_bot_out_shift = 4*mm;
    tgt_al_bot_in_shift = 0.5*mm;

// Additional foil on the inner ring inside the target
    G4double tgt_bot_in_foil_w = 400*um;
    G4double tgt_bot_in_foil_out_w = 2*mm;
    G4double tgt_bot_in_foil_in_height = 12*mm;
    G4double tgt_bot_in_foil_out_height = 7*mm;

//Foils on the inner ring
    G4Tubs* solidTargetFoilOnBotRingIn = new G4Tubs("TargetFoilOnBotRingIn", tgt_kap_ent_r - tgt_bot_in_foil_w, tgt_kap_ent_r, 0.5*tgt_bot_in_foil_in_height, 0., 360.*deg);
    G4LogicalVolume* logicTargetFoilOnBotRingIn = new G4LogicalVolume(solidTargetFoilOnBotRingIn, Polyethylene, "VacCh", 0, 0, 0);
    G4VPhysicalVolume* physiTargetFoilOnBotRingIn = new G4PVPlacement(rm00, G4ThreeVector(0,0,pos_target - 0.5*dz_target + 0.5*tgt_bot_in_foil_in_height), 
								      logicTargetFoilOnBotRingIn, "TargetMaterial_FoilOnBotRingIn", logicVacCh, false, 0);
    logicTargetFoilOnBotRingIn->SetVisAttributes(G4Colour(1.,0.8,0.3,1.));

    G4Tubs* solidTargetFoilOnBotRingOut = new G4Tubs("TargetFoilOnBotRingOut", tgt_kap_ent_r - tgt_bot_in_foil_out_w, tgt_kap_ent_r, 0.5*tgt_bot_in_foil_out_height, 0., 360.*deg);
    G4LogicalVolume* logicTargetFoilOnBotRingOut = new G4LogicalVolume(solidTargetFoilOnBotRingOut, Polyethylene, "VacCh", 0, 0, 0);
    G4VPhysicalVolume* physiTargetFoilOnBotRingOut = new G4PVPlacement(rm00, G4ThreeVector(0,0,pos_target - 0.5*dz_target - target_entrance_w - 0.5*tgt_bot_in_foil_out_height), 
								      logicTargetFoilOnBotRingOut, "TargetMaterial_FoilOnBotRingOut", logicVacCh, false, 0);
    logicTargetFoilOnBotRingOut->SetVisAttributes(G4Colour(1.,0.8,0.3,1.));

//Isolation layers
    G4double r_isolation_out = 0.5*186*mm;
    G4double r_isolation_in = 0.5*155*mm;
    G4double w_isolation = 24*um;

    G4double shift_isolation_targetWindow = 8*mm;

    G4Tubs* solidIsolationBeforeTargetInside = new G4Tubs("IsolationBeforeTargetInside", 0, r_isolation_in, 0.5*w_isolation, 0., 360.*deg);
    G4LogicalVolume* logicIsolationBeforeTargetInside = new G4LogicalVolume(solidIsolationBeforeTargetInside, mylar, "VacCh", 0, 0, 0);
    G4VPhysicalVolume* physiIsolationBeforeTargetInside = new G4PVPlacement(rm00, G4ThreeVector(0,0,-0.5*dz_vac_chamb + teflon_cone_w + shift_isolation_targetWindow + 0.5*w_isolation), 
								            logicIsolationBeforeTargetInside, "IsolationBeforeTargetInside", logicVacCh, false, 0);
    logicIsolationBeforeTargetInside->SetVisAttributes(G4Colour(1.,0.8,0.3,1.));

  G4Tubs* solidTgtAlTopOutRing = new G4Tubs("TgtAlTopRing", tgt_al_top_out_r_min, tgt_al_top_out_r_max, 0.5*tgt_al_top_out_w, 0., 360.*deg);
  G4LogicalVolume* logicTgtAlTopOutRing = new G4LogicalVolume(solidTgtAlTopOutRing, Al, "VacCh", 0, 0, 0);
  G4VPhysicalVolume* physiTgtAlTopOutRing = new G4PVPlacement(rm00, G4ThreeVector(0,0,pos_target + 0.5*dz_target + 0.5*tgt_al_top_out_w),
                                                              logicTgtAlTopOutRing, "TgtAlTopRing", logicVacCh, false, 0);
  logicTgtAlTopOutRing->SetVisAttributes(G4Colour(0.,1.,0.,0.));

  G4Tubs* solidTgtAlBotOutRing = new G4Tubs("TgtAlBotOutRing", tgt_al_bot_out_r_min, tgt_al_bot_out_r_max, 0.5*tgt_al_bot_out_w, 0., 360.*deg);
  G4LogicalVolume* logicTgtAlBotOutRing = new G4LogicalVolume(solidTgtAlBotOutRing, Al, "AlBotOutRing", 0, 0, 0);
  G4VPhysicalVolume* physiTgtAlBotOutRing = new G4PVPlacement(rm00,G4ThreeVector(0,0, pos_target - 0.5*dz_target - tgt_al_bot_out_shift),
            logicTgtAlBotOutRing, "TgtAlBotOutRing",logicVacCh,false,0);
 // logicTgtAlBotOutRing->SetVisAttributes(G4Colour(0.,1.,0.,0.));

  G4Tubs* solidTgtAlBotInRing = new G4Tubs("TgtAlBotInRing", tgt_al_bot_in_r_min, tgt_al_bot_in_r_max, 0.5*tgt_al_bot_in_w, 0., 360.*deg);
  G4LogicalVolume* logicTgtAlBotInRing = new G4LogicalVolume(solidTgtAlBotInRing, Al, "AlBotInRing", 0, 0, 0);
  G4VPhysicalVolume* physiTgtAlBotInRing = new G4PVPlacement(rm00,G4ThreeVector(0,0, pos_target - 0.5*dz_target + tgt_al_bot_in_shift),
                                                             logicTgtAlBotInRing, "TgtAlBotInRing", logicVacCh, false, 0);

 /* G4Tubs* solidTgtAlBotOutRing = new G4Tubs("TgtAlBotOutRing", tgt_al_bot_out_r_min, tgt_al_bot_out_r_max, 0.5*tgt_al_bot_out_w, 0., 360.*deg);
  G4LogicalVolume* logicTgtAlBotOutRing = new G4LogicalVolume(solidTgtAlBotOutRing, Al, "VacCh", 0, 0, 0);
  G4VPhysicalVolume* physiTgtAlBotOutRing = new G4PVPlacement(rm00,G4ThreeVector(0,0, pos_target - 0.5*dz_target),
            logicTgtAlBotOutRing, "TgtAlBotOutRing",logicVacCh,false,0);
  logicTgtAlBotOutRing->SetVisAttributes(G4Colour(0.,1.,0.,0.));

  G4Tubs* solidTgtAlBotInRing = new G4Tubs("TgtAlBotInRing", tgt_al_bot_in_r_min, tgt_al_bot_in_r_max, 0.5*tgt_al_bot_in_w, 0., 360.*deg);
  G4LogicalVolume* logicTgtAlBotInRing = new G4LogicalVolume(solidTgtAlBotInRing, Al, "VacCh", 0, 0, 0);
  G4VPhysicalVolume* physiTgtAlBotInRing = new G4PVPlacement(rm00,G4ThreeVector(0,0, pos_target - 0.5*dz_target - 0.5*tgt_al_bot_in_w),
                                                             logicTgtAlBotInRing, "TgtAlBotInRing", logicVacCh, false, 0);
  logicTgtAlBotInRing->SetVisAttributes(G4Colour(0.,1.,0.,0.));*/

// Ti Foil before entrance window for test 15.06.2021 Shi
 /* G4Tubs* solidTitaniumFoilDeg = new G4Tubs("TitaniumFoilDeg", 0, tgt_ti_ent_win_r, 0.5*TiFoil_ent, 0., 360.*deg);
  G4LogicalVolume* logicTitaniumFoilDeg = new G4LogicalVolume(solidTitaniumFoilDeg, Ti, "TiDeg", 0, 0, 0);
  G4VPhysicalVolume* physiTitaniumFoilDeg  = new G4PVPlacement(rm00,G4ThreeVector(0.,0.,pos_target - 0.5*dz_target - tgt_al_bot_out_w - 2*mm),
                                                               logicTitaniumFoilDeg, "TitaniumFoilDeg", logicVacCh, false, 0);
  logicTitaniumFoilDeg->SetVisAttributes(G4Colour(1.,0.8,0.3,1.));*/

  G4cout << "***************************************************************** pos_target =  " << pos_target << G4endl;

// SDD's //
  G4double SDDBASE_X = 7.*mm; // Aluminium base. Previously 5*.mm
  G4double SDDSTRIPX = 1.*mm; // Ceramic base. Previously 2*.mm
  G4double SDDcaseWidth = 2.*mm;
  G4double scintillator_w = 5.*mm;  // Scintillator Anti. Previously 10*.mm
  G4double board_w = 0.;
  G4double dist_base_board = 0.;
  G4double dist_board_sci = 1.*mm;
  G4double xoffsetsdd;
  G4double xoffsetscin = 0;

  if ((SiddhartaSetup == 2 || SiddhartaSetup == 5 || SiddhartaSetup == 6 || SiddhartaSetup == 7 || SiddhartaSetup == 8)
      && !OldTarget) {
    xoffsetsdd = scintillator_w + board_w + dist_base_board + dist_board_sci;
  } else if (SiddhartaSetup == 2020) {
    xoffsetsdd = 0;
  } else {
    xoffsetsdd = 0;
  }
  G4double DX_SDDBOX = SDDBASE_X + SDDSTRIPX + SDDcaseWidth + xoffsetsdd;
  G4double DY_SDDBOX = 18.*mm;
  G4double DZ_SDDBOX = 35.*mm;
  G4double DY_SDD = 18.*mm;
  G4double DZ_SDD = 33.*mm;
  G4double DX_SDD = SDDBASE_X + SDDSTRIPX + SDDcaseWidth;
  G4double R_SDD;
  G4double SDDX = 450.*um;

  G4double DX_SDDScinBox;
  G4double DX_ScinBOX = 5.*mm;//10*mm;//5.*mm;
  G4double DY_ScinBOX = 2*13.*mm;//18*mm;;//2*13.*mm;
  G4double DZ_ScinBOX = 50.*mm;//35.*mm;//50.*mm;
  G4double R_Scin_INT = 114.5*mm + 0.5*DX_ScinBOX;
  G4double dist_SDDBox_ScinBox = 0;

  if ((SiddhartaSetup == 2 || SiddhartaSetup == 5 || SiddhartaSetup == 6 || SiddhartaSetup == 7  || SiddhartaSetup == 8)
      && !OldTarget) {
    R_SDD = 67.*mm + 0.5*DX_SDD;
  } else if (SiddhartaSetup == 2020) {
    R_SDD = 84.9*mm; //79.*mm;//84.9*mm;
  } else {
    R_SDD = 81.9*mm;
  }

  G4double R_SDD_INT = 67.*mm;

  if (SiddhartaSetup == 2020) {
 //   R_SDD_INT = R_SDD + DX_SDDBOX - SDDBASE_X;
    DX_SDDBOX = SDDBASE_X + SDDSTRIPX + SDDX; // Aluminium + Ceramic + SDD
    R_SDD_INT = R_SDD + DX_SDDBOX/2;

    SDDSTRIPX += SDDX;

//Added Sdd and scin box together
    dist_SDDBox_ScinBox = 1*mm;//21.15*mm;
    DX_SDDScinBox = SDDBASE_X + SDDSTRIPX + dist_SDDBox_ScinBox + DX_ScinBOX; // Aluminium + Ceramic + SDD + distance + scin
    R_SDD_INT = R_SDD + DX_SDDScinBox/2;
    xoffsetsdd = DX_SDDScinBox/2 - 0.5*DX_SDDBOX;
    xoffsetscin = DX_SDDScinBox/2 - 0.5*DX_ScinBOX;
  }

  G4cout << "**************************** SDD Scin Box dimensions (x, y, z) = " << DX_SDDScinBox << " " << DY_SDDBOX << " " << DZ_SDDBOX << G4endl;
  G4cout << "**************************** xOffset = " << xoffsetsdd << G4endl;

//Added Sdd and scin box together
  G4Box* solidSDDBOX = new G4Box("SDDBOX", 0.5*DX_SDDScinBox, 0.5*DY_ScinBOX, 0.5*DZ_ScinBOX);
  G4LogicalVolume* logicSDDBOX = new G4LogicalVolume(solidSDDBOX, vacuum, "VacCh", 0, 0, 0);
  G4VPhysicalVolume* physiSDDBOX;

  G4LogicalVolume* logicSDDBOX_top = new G4LogicalVolume(solidSDDBOX, vacuum, "VacCh", 0, 0, 0);
  G4VPhysicalVolume* physiSDDBOX_top;

  logicSDDBOX->SetVisAttributes(G4Colour(1.,1.,0.));
  logicSDDBOX_top->SetVisAttributes(G4Colour(1.,1.,0.));
  G4double XYZBOX[48][3];
  G4double THETA2_SDD = twopi/24.;

  G4Box* solidScinBOX = new G4Box("ScinBOX", 0.5*DX_ScinBOX, 0.5*DY_ScinBOX, 0.5*DZ_ScinBOX);
  G4LogicalVolume* logicScinBOX = new G4LogicalVolume(solidScinBOX, vacuum, "VacCh", 0, 0, 0);
  G4VPhysicalVolume* physiScinBOX;
  logicScinBOX->SetVisAttributes(G4Colour(1.,1.,0.));

  G4double sddoff = 1.*mm;
  if ((SiddhartaSetup == 2 || SiddhartaSetup == 5 || SiddhartaSetup == 6 || SiddhartaSetup == 7 || SiddhartaSetup == 8 || SiddhartaSetup == 2020)
      && !OldTarget) {
    THETA2_SDD = twopi/24.;
    G4double X2_00 = 0;
    G4double Y2_00 = -R_SDD_INT;
    G4double sddboxoff = 0.*mm;

    G4cout << "**************************** R_SDD_INT = " << R_SDD_INT << " R_SDD: " << R_SDD << G4endl;
    G4cout << "**************************** Z position of the SDD (top, bottom) = " << pos_target + 0.5*DZ_ScinBOX + 0.5*sddboxoff << " " << pos_target - 0.5*DZ_ScinBOX - 0.5*sddboxoff << G4endl;

    G4cout << "**************************** R_Scin_INT = " << R_Scin_INT << G4endl;
    G4cout << "**************************** Dimensions of the scin (x, y, z) = " << DX_ScinBOX << " " << DY_ScinBOX << " " << DZ_ScinBOX << G4endl;

    for (int j=0; j<24; j++) {
      G4double ang = (0.5+j)*THETA2_SDD;

      XYZBOX[2*j][0] = cos(ang)*X2_00 + sin(ang)*Y2_00;
      XYZBOX[2*j][1] = -sin(ang)*X2_00 + cos(ang)*Y2_00;
      XYZBOX[2*j][2] = pos_target - 0.5*DZ_ScinBOX - 0.5*sddboxoff;

      XYZBOX[2*j+1][0] = cos(ang)*X2_00 + sin(ang)*Y2_00;
      XYZBOX[2*j+1][1] = -sin(ang)*X2_00 + cos(ang)*Y2_00;
      XYZBOX[2*j+1][2] = pos_target + 0.5*DZ_ScinBOX + 0.5*sddboxoff;
    }
//Into loop will reduce the code
    G4RotationMatrix* rm39 = new G4RotationMatrix();
    phi = 270.*deg + 0.5*THETA2_SDD;
    rm39->rotateZ(phi);
    G4RotationMatrix* rm40 = new G4RotationMatrix();
    phi = 270.*deg + 1.5*THETA2_SDD;
    rm40->rotateZ(phi);
    G4RotationMatrix* rm41 = new G4RotationMatrix();
    phi = 270.*deg + 2.5*THETA2_SDD;
    rm41->rotateZ(phi);
    G4RotationMatrix* rm42 = new G4RotationMatrix();
    phi = 270.*deg + 3.5*THETA2_SDD;
    rm42->rotateZ(phi);
    G4RotationMatrix* rm43 = new G4RotationMatrix();
    phi = 270.*deg + 4.5*THETA2_SDD;
    rm43->rotateZ(phi);
    G4RotationMatrix* rm44 = new G4RotationMatrix();
    phi = 270.*deg + 5.5*THETA2_SDD;
    rm44->rotateZ(phi);
    G4RotationMatrix* rm45 = new G4RotationMatrix();
    phi = 270.*deg + 6.5*THETA2_SDD;
    rm45->rotateZ(phi);
    G4RotationMatrix* rm46 = new G4RotationMatrix();
    phi = 270.*deg + 7.5*THETA2_SDD;
    rm46->rotateZ(phi);
    G4RotationMatrix* rm47 = new G4RotationMatrix();
    phi = 270.*deg + 8.5*THETA2_SDD;
    rm47->rotateZ(phi);
    G4RotationMatrix* rm48 = new G4RotationMatrix();
    phi = 270.*deg + 9.5*THETA2_SDD;
    rm48->rotateZ(phi);
    G4RotationMatrix* rm49 = new G4RotationMatrix();
    phi = 270.*deg + 10.5*THETA2_SDD;
    rm49->rotateZ(phi);
    G4RotationMatrix* rm50 = new G4RotationMatrix();
    phi = 270.*deg + 11.5*THETA2_SDD;
    rm50->rotateZ(phi);
    G4RotationMatrix* rm51 = new G4RotationMatrix();
    phi = 270.*deg + 12.5*THETA2_SDD;
    rm51->rotateZ(phi);
    G4RotationMatrix* rm52 = new G4RotationMatrix();
    phi = 270.*deg + 13.5*THETA2_SDD;
    rm52->rotateZ(phi);
    G4RotationMatrix* rm53 = new G4RotationMatrix();
    phi = 270.*deg + 14.5*THETA2_SDD;
    rm53->rotateZ(phi);
    G4RotationMatrix* rm54 = new G4RotationMatrix();
    phi = 270.*deg + 15.5*THETA2_SDD;
    rm54->rotateZ(phi);
    G4RotationMatrix* rm55 = new G4RotationMatrix();
    phi = 270.*deg + 16.5*THETA2_SDD;
    rm55->rotateZ(phi);
    G4RotationMatrix* rm56 = new G4RotationMatrix();
    phi = 270.*deg + 17.5*THETA2_SDD;
    rm56->rotateZ(phi);
    G4RotationMatrix* rm57 = new G4RotationMatrix();
    phi = 270.*deg + 18.5*THETA2_SDD;
    rm57->rotateZ(phi);
    G4RotationMatrix* rm58 = new G4RotationMatrix();
    phi = 270.*deg + 19.5*THETA2_SDD;
    rm58->rotateZ(phi);
    G4RotationMatrix* rm59 = new G4RotationMatrix();
    phi = 270.*deg + 20.5*THETA2_SDD;
    rm59->rotateZ(phi);
    G4RotationMatrix* rm60 = new G4RotationMatrix();
    phi = 270.*deg + 21.5*THETA2_SDD;
    rm60->rotateZ(phi);
    G4RotationMatrix* rm61 = new G4RotationMatrix();
    phi = 270.*deg + 22.5*THETA2_SDD;
    rm61->rotateZ(phi);
    G4RotationMatrix* rm62 = new G4RotationMatrix();
    phi = 270.*deg + 23.5*THETA2_SDD;
    rm62->rotateZ(phi);

    physiSDDBOX = new G4PVPlacement(rm39, G4ThreeVector(XYZBOX[0][0],XYZBOX[0][1],XYZBOX[0][2]), logicSDDBOX,
                                    "SDDBOX", logicVacCh, false, 1);
    physiSDDBOX_top = new G4PVPlacement(rm39, G4ThreeVector(XYZBOX[1][0],XYZBOX[1][1],XYZBOX[1][2]), logicSDDBOX_top,
                                    "SDDBOX", logicVacCh, false, 2);
    physiSDDBOX = new G4PVPlacement(rm40, G4ThreeVector(XYZBOX[2][0],XYZBOX[2][1],XYZBOX[2][2]), logicSDDBOX,
                                    "SDDBOX", logicVacCh, false, 3);
    physiSDDBOX_top = new G4PVPlacement(rm40, G4ThreeVector(XYZBOX[3][0],XYZBOX[3][1],XYZBOX[3][2]), logicSDDBOX_top,
                                    "SDDBOX", logicVacCh, false, 4);
    physiSDDBOX = new G4PVPlacement(rm41, G4ThreeVector(XYZBOX[4][0],XYZBOX[4][1],XYZBOX[4][2]), logicSDDBOX,
                                    "SDDBOX", logicVacCh, false, 5);
    physiSDDBOX_top = new G4PVPlacement(rm41, G4ThreeVector(XYZBOX[5][0],XYZBOX[5][1],XYZBOX[5][2]), logicSDDBOX_top,
                                    "SDDBOX", logicVacCh, false, 6);
    physiSDDBOX = new G4PVPlacement(rm42, G4ThreeVector(XYZBOX[6][0],XYZBOX[6][1],XYZBOX[6][2]), logicSDDBOX,
                                    "SDDBOX", logicVacCh, false, 7);
    physiSDDBOX_top = new G4PVPlacement(rm42, G4ThreeVector(XYZBOX[7][0],XYZBOX[7][1],XYZBOX[7][2]), logicSDDBOX_top,
                                    "SDDBOX", logicVacCh, false, 8);
    physiSDDBOX = new G4PVPlacement(rm43, G4ThreeVector(XYZBOX[8][0],XYZBOX[8][1],XYZBOX[8][2]), logicSDDBOX,
                                    "SDDBOX", logicVacCh, false, 9);
    physiSDDBOX_top = new G4PVPlacement(rm43, G4ThreeVector(XYZBOX[9][0],XYZBOX[9][1],XYZBOX[9][2]), logicSDDBOX_top,
                                    "SDDBOX", logicVacCh, false, 10);
    physiSDDBOX = new G4PVPlacement(rm44, G4ThreeVector(XYZBOX[10][0],XYZBOX[10][1],XYZBOX[10][2]), logicSDDBOX,
                                    "SDDBOX", logicVacCh, false, 11);
    physiSDDBOX_top = new G4PVPlacement(rm44, G4ThreeVector(XYZBOX[11][0],XYZBOX[11][1],XYZBOX[11][2]), logicSDDBOX_top,
                                    "SDDBOX", logicVacCh, false, 12);
    physiSDDBOX = new G4PVPlacement(rm45, G4ThreeVector(XYZBOX[12][0],XYZBOX[12][1],XYZBOX[12][2]), logicSDDBOX,
                                    "SDDBOX", logicVacCh, false, 13);
    physiSDDBOX_top = new G4PVPlacement(rm45, G4ThreeVector(XYZBOX[13][0],XYZBOX[13][1],XYZBOX[13][2]), logicSDDBOX_top,
                                    "SDDBOX", logicVacCh, false, 14);
    physiSDDBOX = new G4PVPlacement(rm46, G4ThreeVector(XYZBOX[14][0],XYZBOX[14][1],XYZBOX[14][2]), logicSDDBOX,
                                    "SDDBOX", logicVacCh, false, 15);
    physiSDDBOX_top = new G4PVPlacement(rm46, G4ThreeVector(XYZBOX[15][0],XYZBOX[15][1],XYZBOX[15][2]), logicSDDBOX_top,
                                    "SDDBOX", logicVacCh, false, 16);
    physiSDDBOX = new G4PVPlacement(rm47, G4ThreeVector(XYZBOX[16][0],XYZBOX[16][1],XYZBOX[16][2]), logicSDDBOX,
                                    "SDDBOX", logicVacCh, false, 17);
    physiSDDBOX_top = new G4PVPlacement(rm47, G4ThreeVector(XYZBOX[17][0],XYZBOX[17][1],XYZBOX[17][2]), logicSDDBOX_top,
                                    "SDDBOX", logicVacCh, false, 18);
    physiSDDBOX = new G4PVPlacement(rm48, G4ThreeVector(XYZBOX[18][0],XYZBOX[18][1],XYZBOX[18][2]), logicSDDBOX,
                                    "SDDBOX", logicVacCh, false, 19);
    physiSDDBOX_top = new G4PVPlacement(rm48, G4ThreeVector(XYZBOX[19][0],XYZBOX[19][1],XYZBOX[19][2]), logicSDDBOX_top,
                                    "SDDBOX", logicVacCh, false, 20);
    physiSDDBOX = new G4PVPlacement(rm49, G4ThreeVector(XYZBOX[20][0],XYZBOX[20][1],XYZBOX[20][2]), logicSDDBOX,
                                    "SDDBOX", logicVacCh, false, 21);
    physiSDDBOX_top = new G4PVPlacement(rm49, G4ThreeVector(XYZBOX[21][0],XYZBOX[21][1],XYZBOX[21][2]), logicSDDBOX_top,
                                    "SDDBOX", logicVacCh, false, 22);
    physiSDDBOX = new G4PVPlacement(rm50, G4ThreeVector(XYZBOX[22][0],XYZBOX[22][1],XYZBOX[22][2]), logicSDDBOX,
                                    "SDDBOX", logicVacCh, false, 23);
    physiSDDBOX_top = new G4PVPlacement(rm50, G4ThreeVector(XYZBOX[23][0],XYZBOX[23][1],XYZBOX[23][2]), logicSDDBOX_top,
                                    "SDDBOX", logicVacCh, false, 24);
    physiSDDBOX = new G4PVPlacement(rm51, G4ThreeVector(XYZBOX[24][0],XYZBOX[24][1],XYZBOX[24][2]), logicSDDBOX,
                                    "SDDBOX", logicVacCh, false, 25);
    physiSDDBOX_top = new G4PVPlacement(rm51, G4ThreeVector(XYZBOX[25][0],XYZBOX[25][1],XYZBOX[25][2]), logicSDDBOX_top,
                                    "SDDBOX", logicVacCh, false, 26);
    physiSDDBOX = new G4PVPlacement(rm52, G4ThreeVector(XYZBOX[26][0],XYZBOX[26][1],XYZBOX[26][2]), logicSDDBOX,
                                    "SDDBOX", logicVacCh, false, 27);
    physiSDDBOX_top = new G4PVPlacement(rm52, G4ThreeVector(XYZBOX[27][0],XYZBOX[27][1],XYZBOX[27][2]), logicSDDBOX_top,
                                    "SDDBOX", logicVacCh, false, 28);
    physiSDDBOX = new G4PVPlacement(rm53, G4ThreeVector(XYZBOX[28][0],XYZBOX[28][1],XYZBOX[28][2]), logicSDDBOX,
                                    "SDDBOX", logicVacCh, false, 29);
    physiSDDBOX_top = new G4PVPlacement(rm53, G4ThreeVector(XYZBOX[29][0],XYZBOX[29][1],XYZBOX[29][2]), logicSDDBOX_top,
                                    "SDDBOX", logicVacCh, false, 30);
    physiSDDBOX = new G4PVPlacement(rm54, G4ThreeVector(XYZBOX[30][0],XYZBOX[30][1],XYZBOX[30][2]), logicSDDBOX,
                                    "SDDBOX", logicVacCh, false, 31);
    physiSDDBOX_top = new G4PVPlacement(rm54, G4ThreeVector(XYZBOX[31][0],XYZBOX[31][1],XYZBOX[31][2]), logicSDDBOX_top,
                                    "SDDBOX", logicVacCh, false, 32);
    physiSDDBOX = new G4PVPlacement(rm55, G4ThreeVector(XYZBOX[32][0],XYZBOX[32][1],XYZBOX[32][2]), logicSDDBOX,
                                    "SDDBOX", logicVacCh, false, 33);
    physiSDDBOX_top = new G4PVPlacement(rm55, G4ThreeVector(XYZBOX[33][0],XYZBOX[33][1],XYZBOX[33][2]), logicSDDBOX_top,
                                    "SDDBOX", logicVacCh, false, 34);
    physiSDDBOX = new G4PVPlacement(rm56, G4ThreeVector(XYZBOX[34][0],XYZBOX[34][1],XYZBOX[34][2]), logicSDDBOX,
                                    "SDDBOX", logicVacCh, false, 35);
    physiSDDBOX_top = new G4PVPlacement(rm56, G4ThreeVector(XYZBOX[35][0],XYZBOX[35][1],XYZBOX[35][2]), logicSDDBOX_top,
                                    "SDDBOX", logicVacCh, false, 36);
    physiSDDBOX = new G4PVPlacement(rm57, G4ThreeVector(XYZBOX[36][0],XYZBOX[36][1],XYZBOX[36][2]), logicSDDBOX,
                                    "SDDBOX", logicVacCh, false, 37);
    physiSDDBOX_top = new G4PVPlacement(rm57, G4ThreeVector(XYZBOX[37][0],XYZBOX[37][1],XYZBOX[37][2]), logicSDDBOX_top,
                                    "SDDBOX", logicVacCh, false, 38);
    physiSDDBOX = new G4PVPlacement(rm58, G4ThreeVector(XYZBOX[38][0],XYZBOX[38][1],XYZBOX[38][2]), logicSDDBOX,
                                    "SDDBOX", logicVacCh, false, 39);
    physiSDDBOX_top = new G4PVPlacement(rm58, G4ThreeVector(XYZBOX[39][0],XYZBOX[39][1],XYZBOX[39][2]), logicSDDBOX_top,
                                    "SDDBOX", logicVacCh, false, 40);
    physiSDDBOX = new G4PVPlacement(rm59, G4ThreeVector(XYZBOX[40][0],XYZBOX[40][1],XYZBOX[40][2]), logicSDDBOX,
                                    "SDDBOX", logicVacCh, false, 41);
    physiSDDBOX_top = new G4PVPlacement(rm59, G4ThreeVector(XYZBOX[41][0],XYZBOX[41][1],XYZBOX[41][2]), logicSDDBOX_top,
                                    "SDDBOX", logicVacCh, false, 42);
    physiSDDBOX = new G4PVPlacement(rm60, G4ThreeVector(XYZBOX[42][0],XYZBOX[42][1],XYZBOX[42][2]), logicSDDBOX,
                                    "SDDBOX", logicVacCh, false, 43);
    physiSDDBOX_top = new G4PVPlacement(rm60, G4ThreeVector(XYZBOX[43][0],XYZBOX[43][1],XYZBOX[43][2]), logicSDDBOX_top,
                                    "SDDBOX", logicVacCh, false, 44);
    physiSDDBOX = new G4PVPlacement(rm61, G4ThreeVector(XYZBOX[44][0],XYZBOX[44][1],XYZBOX[44][2]), logicSDDBOX,
                                    "SDDBOX", logicVacCh, false, 45);
    physiSDDBOX_top = new G4PVPlacement(rm61, G4ThreeVector(XYZBOX[45][0],XYZBOX[45][1],XYZBOX[45][2]), logicSDDBOX_top,
                                    "SDDBOX", logicVacCh, false, 46);
    physiSDDBOX = new G4PVPlacement(rm62, G4ThreeVector(XYZBOX[46][0],XYZBOX[46][1],XYZBOX[46][2]), logicSDDBOX,
                                    "SDDBOX", logicVacCh, false, 47);
    physiSDDBOX_top = new G4PVPlacement(rm62, G4ThreeVector(XYZBOX[47][0],XYZBOX[47][1],XYZBOX[47][2]), logicSDDBOX_top,
                                    "SDDBOX", logicVacCh, false, 48);
  }
//-------------------------
  G4double pos_sddbase = -0.5*DX_SDDBOX + 0.5*SDDBASE_X + xoffsetsdd;

  G4double pos_Al_Base = 0.;
  if (SiddhartaSetup == 2020) {
//    pos_Al_Base = -0.5*(DZ_SDDBOX - DZ_SDD);
    pos_Al_Base = 0.5*DZ_ScinBOX - 0.5*DZ_SDD - 0.5*sddoff;
  }

  G4Box* solidSDDBASE = new G4Box("SDDBASE", 0.5*SDDBASE_X, 0.5*DY_SDD, 0.5*DZ_SDD);
  G4LogicalVolume* logicSDDBASE = new G4LogicalVolume(solidSDDBASE, Al, "SDDBASE", 0, 0, 0);
  G4VPhysicalVolume* physiSDDBASE = new G4PVPlacement(rm00, G4ThreeVector(pos_sddbase, 0., pos_Al_Base), logicSDDBASE, "SDDBASE", logicSDDBOX, false, 1);

  physiSDDBASE = new G4PVPlacement(rm00, G4ThreeVector(pos_sddbase, 0., -pos_Al_Base), logicSDDBASE, "SDDBASE", logicSDDBOX_top, false, 1);
  logicSDDBASE->SetVisAttributes(G4Colour(0.,0.,0.));

  G4cout << "***************************************************************** pos_sddbase (x, z) =  " << pos_sddbase << " " << pos_Al_Base << G4endl;
  G4cout << "***************************************************************** pos_sddbase from target =  " << R_SDD_INT - pos_sddbase << G4endl;
  G4cout << "***************************************************************** SDD Aluminium Base dimensions (x,y,z) =  " << SDDBASE_X << " " << DY_SDD << " " << DZ_SDD << G4endl;

  G4double SDDSTRIPY = 1*mm;
  G4double sdd_cera_len = 39*mm; //Previously 34*mm
  G4double dx_film = 7*um;
  G4double SDDY = 0.8*cm;
  G4double SDDZ = 0.8*cm;
  G4double posx_sddstripm = -0.5*DX_SDDBOX + SDDBASE_X + 0.5*SDDSTRIPX + xoffsetsdd;
  G4double posx_sddstripb = -0.5*DX_SDDBOX + SDDBASE_X + 0.5*SDDSTRIPX + xoffsetsdd;
  G4double posz_sddstripb =  0.5*DZ_SDD - 0.5*SDDSTRIPY;
  G4double pos_sddprotection = -0.5*DX_SDDBOX + SDDBASE_X + SDDSTRIPX + 0.5*SDDcaseWidth + xoffsetsdd;  // 2020 obolete
  G4double pos_sddprothole = -0.5*(0.5*DY_SDD - 0.5*SDDSTRIPY) - 0.5*SDDSTRIPY;   // 2020 obsolete
  G4double posx_sddceramic = -0.5*DX_SDDBOX + SDDBASE_X + 0.5*SDDSTRIPX + xoffsetsdd;
  G4double posy_sddceramic = 0.5*SDDSTRIPY + 0.5*(0.5*DY_SDD - 0.5*SDDSTRIPY);
  G4double posz_sddceramic = 0.;

  if (SiddhartaSetup == 2020) {
    dx_film = 0;
    posy_sddceramic = 0.25*DY_SDDBOX;

    posz_sddceramic = 0.5*DZ_ScinBOX - 0.5*sdd_cera_len - 0.5*sddoff;
  }

  G4Box* solidSDDCeramic = new G4Box("SDDCeramic", 0.5*SDDSTRIPX, 0.5*0.5*DY_SDDBOX, 0.5*sdd_cera_len);
  G4LogicalVolume* logicSDDCeramic = new G4LogicalVolume(solidSDDCeramic, vacuum, "SDDCeramic", 0, 0, 0); // changed to ceramic 2023
  G4VPhysicalVolume* physiSDDCeramic = new G4PVPlacement(rm00, G4ThreeVector(posx_sddceramic,posy_sddceramic,posz_sddceramic),
                                                         logicSDDCeramic, "SDDCeramic", logicSDDBOX, false, 1);
  physiSDDCeramic = new G4PVPlacement(rm00, G4ThreeVector(posx_sddceramic,-1*posy_sddceramic,posz_sddceramic),
                                      logicSDDCeramic, "SDDCeramic", logicSDDBOX, false, 2);

  G4LogicalVolume* logicSDDCeramic_top = new G4LogicalVolume(solidSDDCeramic, vacuum, "SDDCeramic", 0, 0, 0); // changed to ceramic 2023
  G4VPhysicalVolume* physiSDDCeramic_top = new G4PVPlacement(rm00, G4ThreeVector(posx_sddceramic,posy_sddceramic,-posz_sddceramic),
                                                         logicSDDCeramic_top, "SDDCeramic", logicSDDBOX_top, false, 1);
  physiSDDCeramic_top = new G4PVPlacement(rm00, G4ThreeVector(posx_sddceramic,-1*posy_sddceramic,-posz_sddceramic),
                                      logicSDDCeramic_top, "SDDCeramic", logicSDDBOX_top, false, 2);

  logicSDDCeramic->SetVisAttributes(G4Colour(1.,0.,0.));
  logicSDDCeramic_top->SetVisAttributes(G4Colour(1.,0.,0.));

  G4cout << "***************************************************************** posxy_sddceramic =  " << posx_sddceramic
         << " " << posy_sddceramic << G4endl;
  G4cout << "***************************************************************** pos_sddceramic from target =  " << R_SDD_INT - posx_sddceramic << G4endl;

// SDD sensitive volume
/*  G4Box* solidSDDfilm = new G4Box("SDDfilm", 0.5*dx_film, 0.5*SDDY, 0.5*SDDZ);
  G4LogicalVolume* logicSDDfilm = new G4LogicalVolume(solidSDDfilm, vacuum, "SDDfilm", 0, 0, 0);
  G4VPhysicalVolume* physiSDDfilm = new G4PVPlacement(rm00, G4ThreeVector(0.5*SDDSTRIPX - 0.5*dx_film,0.,1.2*cm),
                                                      logicSDDfilm, "SDDfilm", logicSDDCeramic, false, 1);
  physiSDDfilm = new G4PVPlacement(rm00, G4ThreeVector(0.5*SDDSTRIPX - 0.5*dx_film,0.,0.4*cm),
                                   logicSDDfilm, "SDDfilm", logicSDDCeramic, false, 2);
  physiSDDfilm = new G4PVPlacement(rm00, G4ThreeVector(0.5*SDDSTRIPX - 0.5*dx_film,0.,-0.4*cm),
                                   logicSDDfilm, "SDDfilm", logicSDDCeramic, false, 3);
  physiSDDfilm = new G4PVPlacement(rm00, G4ThreeVector(0.5*SDDSTRIPX - 0.5*dx_film,0.,-1.2*cm),
                                   logicSDDfilm, "SDDfilm", logicSDDCeramic, false, 4);
  logicSDDfilm->SetVisAttributes(G4Colour(0.,0.,0.));*/

  G4Box* solidSDDCeramicProper = new G4Box("SDDCeramicProper", 0.5*(SDDSTRIPX - SDDX), 0.5*0.5*DY_SDDBOX, 0.5*sdd_cera_len);
  G4LogicalVolume* logicSDDCeramicProper = new G4LogicalVolume(solidSDDCeramicProper, Ceramic, "SDDCeramic", 0, 0, 0); // changed to ceramic 2023
  G4VPhysicalVolume* physiSDDCeramicProper = new G4PVPlacement(rm00, G4ThreeVector(-0.5*SDDX,0.0*mm,0.0*mm),
                                                         logicSDDCeramicProper, "SDDCeramicProper", logicSDDCeramic, false, 1);

  G4LogicalVolume* logicSDDCeramicProper_top = new G4LogicalVolume(solidSDDCeramicProper, Ceramic, "SDDCeramic", 0, 0, 0); // changed to ceramic 2023
  G4VPhysicalVolume* physiSDDCeramicProper_top = new G4PVPlacement(rm00, G4ThreeVector(-0.5*SDDX,0.0*mm,0.0*mm),
                                                         logicSDDCeramicProper_top, "SDDCeramicProper", logicSDDCeramic_top, false, 1);

  logicSDDCeramicProper->SetVisAttributes(G4Colour(1.,0.,0.));
  logicSDDCeramicProper_top->SetVisAttributes(G4Colour(1.,0.,0.));

//sdd and acin box added and sdd spearate with ceramic
 // G4double SDD_PosX = -0.5*DX_SDDBOX + SDDBASE_X + SDDSTRIPX + 0.5*SDDX + xoffsetsdd;

  G4double SDD_PosX = 0.5*SDDSTRIPX - 0.5*SDDX;

  G4Box* solidSDD = new G4Box("SDD", 0.5*SDDX, 0.5*SDDY, 0.5*SDDZ);
  G4LogicalVolume* logicSDD = new G4LogicalVolume(solidSDD, Silicon, "SDD", 0, 0, 0);
 
  G4VPhysicalVolume* physiSDD = new G4PVPlacement(rm00, G4ThreeVector(SDD_PosX, 0, 0.5*sdd_cera_len - 1*mm - 0.5*SDDZ),
                                                  logicSDD, "SDD", logicSDDCeramic, false, 1);
  physiSDD = new G4PVPlacement(rm00, G4ThreeVector(SDD_PosX, 0, 0.5*sdd_cera_len - 1*mm - 1.5*SDDZ),
                                                  logicSDD, "SDD", logicSDDCeramic, false, 2);
  physiSDD = new G4PVPlacement(rm00, G4ThreeVector(SDD_PosX, 0, 0.5*sdd_cera_len - 1*mm - 2.5*SDDZ),
                                                  logicSDD, "SDD", logicSDDCeramic, false, 3);
  physiSDD = new G4PVPlacement(rm00, G4ThreeVector(SDD_PosX, 0, 0.5*sdd_cera_len - 1*mm - 3.5*SDDZ),
                                                  logicSDD, "SDD", logicSDDCeramic, false, 4);

  G4VPhysicalVolume* physiSDD_top = new G4PVPlacement(rm00, G4ThreeVector(SDD_PosX, 0, -0.5*sdd_cera_len + 1*mm + 3.5*SDDZ),
                                                  logicSDD, "SDD", logicSDDCeramic_top, false, 1);
  physiSDD_top = new G4PVPlacement(rm00, G4ThreeVector(SDD_PosX, 0, -0.5*sdd_cera_len + 1*mm + 2.5*SDDZ),
                                                  logicSDD, "SDD", logicSDDCeramic_top, false, 2);
  physiSDD_top = new G4PVPlacement(rm00, G4ThreeVector(SDD_PosX, 0, -0.5*sdd_cera_len + 1*mm + 1.5*SDDZ),
                                                  logicSDD, "SDD", logicSDDCeramic_top, false, 3);
  physiSDD_top = new G4PVPlacement(rm00, G4ThreeVector(SDD_PosX, 0, -0.5*sdd_cera_len + 1*mm + 0.5*SDDZ),
                                                  logicSDD, "SDD", logicSDDCeramic_top, false, 4);

  G4cout << "***************************************************************** Pos of SDD from ceramic =  " << SDD_PosX << G4endl;
  G4cout << "***************************************************************** Pos of SDD from target =  " << R_SDD_INT - posx_sddceramic - SDD_PosX << G4endl;

//
//
//Need correcting positions after correcting SDDBox > Is it used in the analysis ?
//
//
  if ((SiddhartaSetup == 2 || SiddhartaSetup == 5 || SiddhartaSetup == 6 || SiddhartaSetup == 7 || SiddhartaSetup == 8 || SiddhartaSetup == 2020)
      && !OldTarget) {
    G4double Xb;
    G4double Yb;
    G4double Zb;
    G4double Xd;
    G4double Yd;
    G4double Zd;

    for (int j=0; j<24; j++) {
      Xb = XYZBOX[2*j][0];
      Yb = pos_vac_chamb + XYZBOX[2*j][2];
      Zb = -XYZBOX[2*j][1];

      analysis->SDD_angle[16*j+0]  = 90.*deg - (0.5 + j)*twopi/24.;
      analysis->SDD_angle[16*j+1]  = 90.*deg - (0.5 + j)*twopi/24.;
      analysis->SDD_angle[16*j+2]  = 90.*deg - (0.5 + j)*twopi/24.;
      analysis->SDD_angle[16*j+3]  = 90.*deg - (0.5 + j)*twopi/24.;
      analysis->SDD_angle[16*j+4]  = 90.*deg - (0.5 + j)*twopi/24.;
      analysis->SDD_angle[16*j+5]  = 90.*deg - (0.5 + j)*twopi/24.;
      analysis->SDD_angle[16*j+6]  = 90.*deg - (0.5 + j)*twopi/24.;
      analysis->SDD_angle[16*j+7]  = 90.*deg - (0.5 + j)*twopi/24.;
      analysis->SDD_angle[16*j+8]  = 90.*deg - (0.5 + j)*twopi/24.;
      analysis->SDD_angle[16*j+9]  = 90.*deg - (0.5 + j)*twopi/24.;
      analysis->SDD_angle[16*j+10] = 90.*deg - (0.5 + j)*twopi/24.;
      analysis->SDD_angle[16*j+11] = 90.*deg - (0.5 + j)*twopi/24.;
      analysis->SDD_angle[16*j+12] = 90.*deg - (0.5 + j)*twopi/24.;
      analysis->SDD_angle[16*j+13] = 90.*deg - (0.5 + j)*twopi/24.;
      analysis->SDD_angle[16*j+14] = 90.*deg - (0.5 + j)*twopi/24.;
      analysis->SDD_angle[16*j+15] = 90.*deg - (0.5 + j)*twopi/24.;

// Y in G4Placement is Z; Z in G4Placement is Y
      Xd = XYZBOX[2*j][0] + posx_sddceramic + SDD_PosX;
      Yd = pos_vac_chamb + XYZBOX[2*j][2] + 0.5*sdd_cera_len - 1*mm - 0.5*SDDZ;
//      Xd = XYZBOX[2*j][0] + posx_sddceramic + 0.5*SDDSTRIPX - dx_film - 0.5*SDDX;
//      Yd = pos_vac_chamb + XYZBOX[2*j][2] + 1.2*cm;
      Zd = -XYZBOX[2*j][1] - posy_sddceramic;
      analysis->SDD_pos[16*j+0][0] = (Xd - Xb)*cos(analysis->SDD_angle[16*j+0]) + (Zd - Zb)*sin(analysis->SDD_angle[16*j+0]) + Xb;
      analysis->SDD_pos[16*j+0][1] = Yd;
      analysis->SDD_pos[16*j+0][2] = -(Xd - Xb)*sin(analysis->SDD_angle[16*j+0]) + (Zd - Zb)*cos(analysis->SDD_angle[16*j+0]) + Zb;

      Yd = pos_vac_chamb + XYZBOX[2*j][2] + 0.5*sdd_cera_len - 1*mm - 1.5*SDDZ;
//      Yd = pos_vac_chamb + XYZBOX[2*j][2] + 0.4*cm;
      Zd = -XYZBOX[2*j][1] - posy_sddceramic;
      analysis->SDD_pos[16*j+1][0] = (Xd - Xb)*cos(analysis->SDD_angle[16*j+1]) + (Zd - Zb)*sin(analysis->SDD_angle[16*j+1]) + Xb;
      analysis->SDD_pos[16*j+1][1] = Yd;
      analysis->SDD_pos[16*j+1][2] = -(Xd - Xb)*sin(analysis->SDD_angle[16*j+1]) + (Zd - Zb)*cos(analysis->SDD_angle[16*j+1]) + Zb;

      Yd = pos_vac_chamb + XYZBOX[2*j][2] + 0.5*sdd_cera_len - 1*mm - 2.5*SDDZ;
//      Yd = pos_vac_chamb + XYZBOX[2*j][2] - 0.4*cm;
      Zd = -XYZBOX[2*j][1] - posy_sddceramic;
      analysis->SDD_pos[16*j+2][0] = (Xd - Xb)*cos(analysis->SDD_angle[16*j+2]) + (Zd - Zb)*sin(analysis->SDD_angle[16*j+2]) + Xb;
      analysis->SDD_pos[16*j+2][1] = Yd;
      analysis->SDD_pos[16*j+2][2] = -(Xd - Xb)*sin(analysis->SDD_angle[16*j+2]) + (Zd - Zb)*cos(analysis->SDD_angle[16*j+2]) + Zb;

      Yd = pos_vac_chamb + XYZBOX[2*j][2] + 0.5*sdd_cera_len - 1*mm - 3.5*SDDZ;
//      Yd = pos_vac_chamb + XYZBOX[2*j][2] - 1.2*cm;
      Zd = -XYZBOX[2*j][1] + posy_sddceramic;
      analysis->SDD_pos[16*j+3][0] = (Xd - Xb)*cos(analysis->SDD_angle[16*j+3]) + (Zd - Zb)*sin(analysis->SDD_angle[16*j+3]) + Xb;
      analysis->SDD_pos[16*j+3][1] = Yd;
      analysis->SDD_pos[16*j+3][2] = -(Xd - Xb)*sin(analysis->SDD_angle[16*j+3]) + (Zd - Zb)*cos(analysis->SDD_angle[16*j+3]) + Zb;

      Yd = pos_vac_chamb + XYZBOX[2*j][2] + 0.5*sdd_cera_len - 1*mm - 0.5*SDDZ;
//      Yd = pos_vac_chamb + XYZBOX[2*j][2] + 1.2*cm;
      Zd = -XYZBOX[2*j][1] + posy_sddceramic;
      analysis->SDD_pos[16*j+4][0] = (Xd - Xb)*cos(analysis->SDD_angle[16*j+4]) + (Zd - Zb)*sin(analysis->SDD_angle[16*j+4]) + Xb;
      analysis->SDD_pos[16*j+4][1] = Yd;
      analysis->SDD_pos[16*j+4][2] = -(Xd - Xb)*sin(analysis->SDD_angle[16*j+4]) + (Zd - Zb)*cos(analysis->SDD_angle[16*j+4]) + Zb;

      Yd = pos_vac_chamb + XYZBOX[2*j][2] + 0.5*sdd_cera_len - 1*mm - 1.5*SDDZ;
//      Yd = pos_vac_chamb + XYZBOX[2*j][2] + 0.4*cm;
      Zd = -XYZBOX[2*j][1] + posy_sddceramic;
      analysis->SDD_pos[16*j+5][0] = (Xd - Xb)*cos(analysis->SDD_angle[16*j+5]) + (Zd - Zb)*sin(analysis->SDD_angle[16*j+5]) + Xb;
      analysis->SDD_pos[16*j+5][1] = Yd;
      analysis->SDD_pos[16*j+5][2] = -(Xd - Xb)*sin(analysis->SDD_angle[16*j+5]) + (Zd - Zb)*cos(analysis->SDD_angle[16*j+5]) + Zb;

      Yd = pos_vac_chamb + XYZBOX[2*j][2] + 0.5*sdd_cera_len - 1*mm - 2.5*SDDZ;
//      Yd = pos_vac_chamb + XYZBOX[2*j][2] - 0.4*cm;
      Zd = -XYZBOX[2*j][1] - posy_sddceramic;
      analysis->SDD_pos[16*j+6][0] = (Xd - Xb)*cos(analysis->SDD_angle[16*j+6]) + (Zd - Zb)*sin(analysis->SDD_angle[16*j+6]) + Xb;
      analysis->SDD_pos[16*j+6][1] = Yd;
      analysis->SDD_pos[16*j+6][2] = -(Xd - Xb)*sin(analysis->SDD_angle[16*j+6]) + (Zd - Zb)*cos(analysis->SDD_angle[16*j+6]) + Zb;

      Yd = pos_vac_chamb + XYZBOX[2*j][2] + 0.5*sdd_cera_len - 1*mm - 3.5*SDDZ;
//      Yd = pos_vac_chamb + XYZBOX[2*j][2] - 1.2*cm;
      Zd = -XYZBOX[2*j][1] - posy_sddceramic;
      analysis->SDD_pos[16*j+7][0] = (Xd - Xb)*cos(analysis->SDD_angle[16*j+7]) + (Zd - Zb)*sin(analysis->SDD_angle[16*j+7]) + Xb;
      analysis->SDD_pos[16*j+7][1] = Yd;
      analysis->SDD_pos[16*j+7][2] = -(Xd - Xb)*sin(analysis->SDD_angle[16*j+7]) + (Zd - Zb)*cos(analysis->SDD_angle[16*j+7]) + Zb;

      Xd = XYZBOX[2*j+1][0] + posx_sddceramic + SDD_PosX;
      Yd = pos_vac_chamb + XYZBOX[2*j][2] + 0.5*sdd_cera_len - 1*mm - 1.5*SDDZ;
//      Xd = XYZBOX[2*j+1][0] + posx_sddceramic + 0.5*SDDSTRIPX - dx_film - 0.5*SDDX;
//      Yd = pos_vac_chamb + XYZBOX[2*j+1][2] + 1.2*cm;
      Zd = -XYZBOX[2*j+1][1] - posy_sddceramic;
      analysis->SDD_pos[16*j+8][0] = (Xd - Xb)*cos(analysis->SDD_angle[16*j+8]) + (Zd - Zb)*sin(analysis->SDD_angle[16*j+8]) + Xb;
      analysis->SDD_pos[16*j+8][1] = Yd;
      analysis->SDD_pos[16*j+8][2] = -(Xd - Xb)*sin(analysis->SDD_angle[16*j+8]) + (Zd - Zb)*cos(analysis->SDD_angle[16*j+8]) + Zb;

      Yd = pos_vac_chamb + XYZBOX[2*j][2] + 0.5*sdd_cera_len - 1*mm - 1.5*SDDZ;
//      Yd = pos_vac_chamb + XYZBOX[2*j+1][2] + 0.4*cm;
      Zd = -XYZBOX[2*j+1][1] + posy_sddceramic;
      analysis->SDD_pos[16*j+9][0] = (Xd - Xb)*cos(analysis->SDD_angle[16*j+9]) + (Zd - Zb)*sin(analysis->SDD_angle[16*j+9]) + Xb;
      analysis->SDD_pos[16*j+9][1] = Yd;
      analysis->SDD_pos[16*j+9][2] = -(Xd - Xb)*sin(analysis->SDD_angle[16*j+9]) + (Zd - Zb)*cos(analysis->SDD_angle[16*j+9]) + Zb;

      Yd = pos_vac_chamb + XYZBOX[2*j][2] + 0.5*sdd_cera_len - 1*mm - 2.5*SDDZ;
//      Yd = pos_vac_chamb + XYZBOX[2*j][2] - 0.4*cm;
      Zd = -XYZBOX[2*j+1][1] + posy_sddceramic;
      analysis->SDD_pos[16*j+10][0] = (Xd - Xb)*cos(analysis->SDD_angle[16*j+10]) + (Zd - Zb)*sin(analysis->SDD_angle[16*j+10]) + Xb;
      analysis->SDD_pos[16*j+10][1] = Yd;
      analysis->SDD_pos[16*j+10][2] = -(Xd - Xb)*sin(analysis->SDD_angle[16*j+10]) + (Zd - Zb)*cos(analysis->SDD_angle[16*j+10]) + Zb;

      Yd = pos_vac_chamb + XYZBOX[2*j][2] + 0.5*sdd_cera_len - 1*mm - 3.5*SDDZ;
//      Yd = pos_vac_chamb + XYZBOX[2*j][2] - 1.2*cm;
      Zd = -XYZBOX[2*j+1][1] + posy_sddceramic;
      analysis->SDD_pos[16*j+11][0] = (Xd - Xb)*cos(analysis->SDD_angle[16*j+11]) + (Zd - Zb)*sin(analysis->SDD_angle[16*j+11]) + Xb;
      analysis->SDD_pos[16*j+11][1] = Yd;
      analysis->SDD_pos[16*j+11][2] = -(Xd - Xb)*sin(analysis->SDD_angle[16*j+11]) + (Zd - Zb)*cos(analysis->SDD_angle[16*j+11]) + Zb;

      Yd = pos_vac_chamb + XYZBOX[2*j][2] + 0.5*sdd_cera_len - 1*mm - 0.5*SDDZ;
//      Yd = pos_vac_chamb + XYZBOX[2*j][2] + 1.2*cm;
      Zd = -XYZBOX[2*j+1][1] - posy_sddceramic;
      analysis->SDD_pos[16*j+12][0] =  (Xd-Xb)*cos(analysis->SDD_angle[16*j+12]) + (Zd-Zb)*sin(analysis->SDD_angle[16*j+12]) + Xb;
      analysis->SDD_pos[16*j+12][1] = Yd;
      analysis->SDD_pos[16*j+12][2] = -(Xd-Xb)*sin(analysis->SDD_angle[16*j+12]) + (Zd-Zb)*cos(analysis->SDD_angle[16*j+12]) + Zb;

      Yd = pos_vac_chamb + XYZBOX[2*j][2] + 0.5*sdd_cera_len - 1*mm - 1.5*SDDZ;
//      Yd = pos_vac_chamb + XYZBOX[2*j][2] + 0.4*cm;
      Zd = -XYZBOX[2*j+1][1] - posy_sddceramic;
      analysis->SDD_pos[16*j+13][0] = (Xd - Xb)*cos(analysis->SDD_angle[16*j+13]) + (Zd - Zb)*sin(analysis->SDD_angle[16*j+13]) + Xb;
      analysis->SDD_pos[16*j+13][1] = Yd;
      analysis->SDD_pos[16*j+13][2] = -(Xd - Xb)*sin(analysis->SDD_angle[16*j+13]) + (Zd - Zb)*cos(analysis->SDD_angle[16*j+13]) + Zb;

      Yd = pos_vac_chamb + XYZBOX[2*j][2] + 0.5*sdd_cera_len - 1*mm - 2.5*SDDZ;
//      Yd = pos_vac_chamb + XYZBOX[2*j][2] - 0.4*cm;
      Zd = -XYZBOX[2*j+1][1] - posy_sddceramic;
      analysis->SDD_pos[16*j+14][0] = (Xd - Xb)*cos(analysis->SDD_angle[16*j+14]) + (Zd - Zb)*sin(analysis->SDD_angle[16*j+14]) + Xb;
      analysis->SDD_pos[16*j+14][1] = Yd;
      analysis->SDD_pos[16*j+14][2] = -(Xd - Xb)*sin(analysis->SDD_angle[16*j+14]) + (Zd - Zb)*cos(analysis->SDD_angle[16*j+14]) + Zb;

      Yd = pos_vac_chamb + XYZBOX[2*j][2] + 0.5*sdd_cera_len - 1*mm - 3.5*SDDZ;
//      Yd = pos_vac_chamb + XYZBOX[2*j][2] - 1.2*cm;
      Zd = -XYZBOX[2*j+1][1] + posy_sddceramic;
      analysis->SDD_pos[16*j+15][0] = (Xd - Xb)*cos(analysis->SDD_angle[16*j+15]) + (Zd - Zb)*sin(analysis->SDD_angle[16*j+15]) + Xb;
      analysis->SDD_pos[16*j+15][1] = Yd;
      analysis->SDD_pos[16*j+15][2] = -(Xd - Xb)*sin(analysis->SDD_angle[16*j+15]) + (Zd - Zb)*cos(analysis->SDD_angle[16*j+15]) + Zb;
    }
  }
  logicSDD->SetSensitiveDetector(aTrackerSD);
  logicSDD->SetVisAttributes(G4Colour(0.,1.,0.));

// Scintillator //
  if (SiddhartaSetup == 2 || SiddhartaSetup == 6 || SiddhartaSetup == 7 || SiddhartaSetup == 8 || SiddhartaSetup == 2020) {
    G4double pos_scint = -0.5*DX_SDDBOX + 0.5*scintillator_w;

    if (SiddhartaSetup == 2020) {
   //   pos_scint = 0;

//Sdd and scin box together
      pos_scint = -xoffsetscin;
    }

    G4Box* solidScintAnti = new G4Box("ScintAnti", 0.5*DX_ScinBOX, 0.5*DY_ScinBOX, 0.5*DZ_ScinBOX);
    G4LogicalVolume* logicScintAnti = new G4LogicalVolume(solidScintAnti, BC420, "ScintAnti", 0, 0, 0);
  //  G4VPhysicalVolume* physiScintAnti = new G4PVPlacement(rm00, G4ThreeVector(pos_scint,0.,0.), logicScintAnti,
  //                                                        "ScintAnti", logicScinBOX, false, 1);

//Sdd and scin box together
    G4VPhysicalVolume* physiScintAnti = new G4PVPlacement(rm00, G4ThreeVector(pos_scint,0.,0.), logicScintAnti,
                                                          "ScintAnti", logicSDDBOX, false, 1);

    physiScintAnti = new G4PVPlacement(rm00, G4ThreeVector(pos_scint,0.,0.), logicScintAnti,
                                                          "ScintAnti", logicSDDBOX_top, false, 1);
    logicScintAnti->SetSensitiveDetector(ScintAntiSD);
  }

  if (SiddhartaSetup == 7 || SiddhartaSetup == 8 || SiddhartaSetup == 2020) {
// Anticoincidence 1 //
    G4double anti1_x = 0.5*(2.*half_box - w_scint_thick - 120.*mm);
    G4double anti1_y = 120.*mm;
    G4double anti1_z = w_scint_thin;
    G4double pos_anti1 = pos_shield + 0.5*Z_SHIEL + 0.5*w_scint_thin + 1*mm;
    G4double posx_anti1 = half_box - 0.5*anti1_x - 0.5*w_scint_thick;

    if (SiddhartaSetup != 2020) {
      G4Box *solidAnti1 = new G4Box("Anti1", 0.5*anti1_x, 0.5*anti1_y, 0.5*anti1_z);
      G4LogicalVolume* logicAnti1 = new G4LogicalVolume(solidAnti1, BC420, "World", 0, 0, 0);
      G4VPhysicalVolume* physiAnti1 = new G4PVPlacement(rm01, G4ThreeVector(posx_anti1,pos_anti1,0), logicAnti1,
                                                      "Anti1", World, false, 0);
      physiAnti1 = new G4PVPlacement(rm01, G4ThreeVector(-posx_anti1,pos_anti1,0), logicAnti1, "Anti1", World, false, 1);
      logicAnti1->SetSensitiveDetector(Anti1_SD);
      logicAnti1->SetVisAttributes(G4Colour(1.,1.,0.));
    }
// Anticoincidence 2 //  To correct the dimensions and positioning
  //  G4double anti2_x = 2.*half_box - w_scint_thick;
    G4double xadd;

    if (SiddhartaSetup == 8 || SiddhartaSetup == 2020  ) {
      xadd = 5.*mm;
    } else {
      xadd = 30.*mm;
    }
   // G4double anti2_y = 0.5*(2.*half_box - w_scint_thick - 120.*mm) + xadd;
   // G4double anti2_z = w_scint_thin;
  //  G4double pos_anti2 = pos_shield + 0.5*Z_SHIEL + 0.5*w_scint_thin + 1.*mm;
  //  G4double posz_anti2 = half_box - 0.5*anti2_y - 0.5*w_scint_thick + xadd;

// 2023 setup
    G4double anti2_x = 550.*mm;
    G4double anti2_y = 110.*mm;
    G4double anti2_z = 10.*mm;

    G4double pos_anti2 = pos_kdtop;
    G4double posz_anti2 = 0.5*dy_kdtop + 10.*mm + 0.5*anti2_y;
// 2023 setup

    G4Box *solidAnti2 = new G4Box("Anti2", 0.5*anti2_x, 0.5*anti2_y, 0.5*anti2_z);
    G4LogicalVolume* logicAnti2 = new G4LogicalVolume(solidAnti2, BC420, "World", 0, 0, 0);
    G4VPhysicalVolume* physiAnti2 = new G4PVPlacement(rm01, G4ThreeVector(0,pos_anti2,posz_anti2), logicAnti2,
                                                      "Anti2", World, false, 0);
    physiAnti2 = new G4PVPlacement(rm01, G4ThreeVector(0,pos_anti2,-posz_anti2), logicAnti2, "Anti2", World, false, 1);
    logicAnti2->SetSensitiveDetector(Anti2_SD);
    logicAnti2->SetVisAttributes(G4Colour(0.,1.,0.));

// Anticoincidence 3 //
    G4double anti3_y = 250.*mm;
    G4double anti3_z = w_scint_thick;
    G4double pos_anti3 = pos_vac_chamb - 15.*mm;
    G4double posz_anti3;

    if (SiddhartaSetup == 8) {
      posz_anti3 = 220.*mm;
    } else if (SiddhartaSetup == 2020) {
      anti3_y = 260.*mm;
      posz_anti3 = r1_vac_chamb + 8*mm + 10*mm;
      pos_anti3 = pos_vac_chamb;
    } else  {
      posz_anti3 = 235.*mm;
    }
    G4double Length = twopi*(posz_anti3 - 0.5*w_scint_thick);
    G4double anti3_x;

    if (SiddhartaSetup == 8) {
      anti3_x = Length/12. + 2.*mm;
    } else if (SiddhartaSetup == 2020) {
      anti3_x = Length/12. + 2.*mm;
    } else {
      anti3_x = 123.*mm;
    }
    G4double THETA_ANTI = -twopi/12.;

    G4double X2_00 = 0;
    G4double Y2_00 = posz_anti3;
    G4double sddoff = 0.*mm;
    G4double X2_100 = cos(0.5*THETA_ANTI)*X2_00 + sin(0.5*THETA_ANTI)*Y2_00;
    G4double Z2_100 = -sin(0.5*THETA_ANTI)*X2_00 + cos(0.5*THETA_ANTI)*Y2_00;
    G4double Y2_100 = pos_anti3;
    G4double X2_101 = cos(1.5*THETA_ANTI)*X2_00 + sin(1.5*THETA_ANTI)*Y2_00;
    G4double Z2_101 = -sin(1.5*THETA_ANTI)*X2_00 + cos(1.5*THETA_ANTI)*Y2_00;
    G4double Y2_101 = pos_anti3;
    G4double X2_102 = cos(2.5*THETA_ANTI)*X2_00 + sin(2.5*THETA_ANTI)*Y2_00;
    G4double Z2_102 = -sin(2.5*THETA_ANTI)*X2_00 + cos(2.5*THETA_ANTI)*Y2_00;
    G4double Y2_102 = pos_anti3;
    G4double X2_103 = cos(3.5*THETA_ANTI)*X2_00 + sin(3.5*THETA_ANTI)*Y2_00;
    G4double Z2_103 = -sin(3.5*THETA_ANTI)*X2_00 + cos(3.5*THETA_ANTI)*Y2_00;
    G4double Y2_103 = pos_anti3;
    G4double X2_104 = cos(4.5*THETA_ANTI)*X2_00 + sin(4.5*THETA_ANTI)*Y2_00;
    G4double Z2_104 = -sin(4.5*THETA_ANTI)*X2_00 + cos(4.5*THETA_ANTI)*Y2_00;
    G4double Y2_104 = pos_anti3;
    G4double X2_105 = cos(5.5*THETA_ANTI)*X2_00 + sin(5.5*THETA_ANTI)*Y2_00;
    G4double Z2_105 = -sin(5.5*THETA_ANTI)*X2_00 + cos(5.5*THETA_ANTI)*Y2_00;
    G4double Y2_105 = pos_anti3;
    G4double X2_106 = cos(6.5*THETA_ANTI)*X2_00 + sin(6.5*THETA_ANTI)*Y2_00;
    G4double Z2_106 = -sin(6.5*THETA_ANTI)*X2_00 + cos(6.5*THETA_ANTI)*Y2_00;
    G4double Y2_106 = pos_anti3;
    G4double X2_107 = cos(7.5*THETA_ANTI)*X2_00 + sin(7.5*THETA_ANTI)*Y2_00;
    G4double Z2_107 = -sin(7.5*THETA_ANTI)*X2_00 + cos(7.5*THETA_ANTI)*Y2_00;
    G4double Y2_107 = pos_anti3;
    G4double X2_108 = cos(8.5*THETA_ANTI)*X2_00 + sin(8.5*THETA_ANTI)*Y2_00;
    G4double Z2_108 = -sin(8.5*THETA_ANTI)*X2_00 + cos(8.5*THETA_ANTI)*Y2_00;
    G4double Y2_108 = pos_anti3;
    G4double X2_109 = cos(9.5*THETA_ANTI)*X2_00 + sin(9.5*THETA_ANTI)*Y2_00;
    G4double Z2_109 = -sin(9.5*THETA_ANTI)*X2_00 + cos(9.5*THETA_ANTI)*Y2_00;
    G4double Y2_109 = pos_anti3;
    G4double X2_110 = cos(10.5*THETA_ANTI)*X2_00 + sin(10.5*THETA_ANTI)*Y2_00;
    G4double Z2_110 = -sin(10.5*THETA_ANTI)*X2_00 + cos(10.5*THETA_ANTI)*Y2_00;
    G4double Y2_110 = pos_anti3;
    G4double X2_111 = cos(11.5*THETA_ANTI)*X2_00 + sin(11.5*THETA_ANTI)*Y2_00;
    G4double Z2_111 = -sin(11.5*THETA_ANTI)*X2_00 + cos(11.5*THETA_ANTI)*Y2_00;
    G4double Y2_111 = pos_anti3;

    G4RotationMatrix* rm39 = new G4RotationMatrix();
    phi = -0.5*THETA_ANTI;
    rm39->rotateY(phi);
    G4RotationMatrix* rm40 = new G4RotationMatrix();
    phi = -1.5*THETA_ANTI;
    rm40->rotateY(phi);
    G4RotationMatrix* rm41 = new G4RotationMatrix();
    phi = -2.5*THETA_ANTI;
    rm41->rotateY(phi);
    G4RotationMatrix* rm42 = new G4RotationMatrix();
    phi = -3.5*THETA_ANTI;
    rm42->rotateY(phi);
    G4RotationMatrix* rm43 = new G4RotationMatrix();
    phi = -4.5*THETA_ANTI;
    rm43->rotateY(phi);
    G4RotationMatrix* rm44 = new G4RotationMatrix();
    phi = -5.5*THETA_ANTI;
    rm44->rotateY(phi);
    G4RotationMatrix* rm45 = new G4RotationMatrix();
    phi = -6.5*THETA_ANTI;
    rm45->rotateY(phi);
    G4RotationMatrix* rm46 = new G4RotationMatrix();
    phi = -7.5*THETA_ANTI;
    rm46->rotateY(phi);
    G4RotationMatrix* rm47 = new G4RotationMatrix();
    phi = -8.5*THETA_ANTI;
    rm47->rotateY(phi);
    G4RotationMatrix* rm48 = new G4RotationMatrix();
    phi = -9.5*THETA_ANTI;
    rm48->rotateY(phi);
    G4RotationMatrix* rm49 = new G4RotationMatrix();
    phi = -10.5*THETA_ANTI;
    rm49->rotateY(phi);
    G4RotationMatrix* rm50 = new G4RotationMatrix();
    phi = -11.5*THETA_ANTI;
    rm50->rotateY(phi);

    G4Box *solidAnti3 = new G4Box("Anti3", 0.5*anti3_x, 0.5*anti3_y, 0.5*anti3_z);
    G4LogicalVolume* logicAnti3 = new G4LogicalVolume(solidAnti3, BC420, "World", 0, 0, 0);
    G4VPhysicalVolume* physiAnti3 = new G4PVPlacement(rm39, G4ThreeVector(X2_100,Y2_100,Z2_100), logicAnti3,
                                                      "Anti3", World, false, 0);
    physiAnti3 = new G4PVPlacement(rm40, G4ThreeVector(X2_101,Y2_101,Z2_101), logicAnti3, "Anti3", World, false, 1);
    physiAnti3 = new G4PVPlacement(rm41, G4ThreeVector(X2_102,Y2_102,Z2_102), logicAnti3, "Anti3", World, false, 2);
    physiAnti3 = new G4PVPlacement(rm42, G4ThreeVector(X2_103,Y2_103,Z2_103), logicAnti3, "Anti3", World, false, 3);
    physiAnti3 = new G4PVPlacement(rm43, G4ThreeVector(X2_104,Y2_104,Z2_104), logicAnti3, "Anti3", World, false, 4);
    physiAnti3 = new G4PVPlacement(rm44, G4ThreeVector(X2_105,Y2_105,Z2_105), logicAnti3, "Anti3", World, false, 5);
    physiAnti3 = new G4PVPlacement(rm45, G4ThreeVector(X2_106,Y2_106,Z2_106), logicAnti3, "Anti3", World, false, 6);
    physiAnti3 = new G4PVPlacement(rm46, G4ThreeVector(X2_107,Y2_107,Z2_107), logicAnti3, "Anti3", World, false, 7);
    physiAnti3 = new G4PVPlacement(rm47, G4ThreeVector(X2_108,Y2_108,Z2_108), logicAnti3, "Anti3", World, false, 8);
    physiAnti3 = new G4PVPlacement(rm48, G4ThreeVector(X2_109,Y2_109,Z2_109), logicAnti3, "Anti3", World, false, 9);
    physiAnti3 = new G4PVPlacement(rm49, G4ThreeVector(X2_110,Y2_110,Z2_110), logicAnti3, "Anti3", World, false, 10);
    physiAnti3 = new G4PVPlacement(rm50, G4ThreeVector(X2_111,Y2_111,Z2_111), logicAnti3, "Anti3", World, false, 11);
    logicAnti3->SetSensitiveDetector(Anti3_SD);
    logicAnti3->SetVisAttributes(G4Colour(0.,0.,1.));

// Anticoincidence 5 //
    G4double anti5_x = 2.*half_box + w_scint_thick;
    G4double anti5_y = w_scint_thick;
    G4double anti5_z = 2.*half_box + w_scint_thick + 50.*mm;
    G4double posx_anti5 = 0.*mm;
    G4double posy_anti5 = pos_shield + 0.5*Z_SHIEL + 360.*mm + 0.5*w_scint_thick + 1.*mm - 80.*mm;
    G4double posz_anti5 = 0;

    if (SiddhartaSetup == 2020) {
      posy_anti5 = pos_shield + 0.5*Z_SHIEL + 380.*mm;
    }
    G4Box* solidAnti5 = new G4Box("Anti5", 0.5*anti5_x, 0.5*anti5_y, 0.5*anti5_z);
    G4LogicalVolume* logicAnti5 = new G4LogicalVolume(solidAnti5, vacuum, "World", 0, 0, 0);
    G4VPhysicalVolume* physiAnti5 = new G4PVPlacement(rm00, G4ThreeVector(posx_anti5,posy_anti5,posz_anti5), logicAnti5,
                                                      "Anti5", World, false, 0);
    logicAnti5->SetSensitiveDetector(Anti5_SD);
    logicAnti5->SetVisAttributes(G4Colour(0.,1.,0.));
  }
// Ghost Detector //
  G4double ghost_x;
  G4double ghost_y;
  G4double ghost_z;
  G4double pos_ghost;

  if (SiddhartaSetup == 2) {
    ghost_x = 200.0*mm;
    ghost_y = 200.0*mm;
    ghost_z = 1.0*mm;
    pos_ghost = -120.*mm;
  } else  if (SiddhartaSetup == 1) {
    ghost_x = 200.0*mm;
    ghost_y = 200.0*mm;
    ghost_z = 1.0*mm;
    pos_ghost = -120.*mm;
  } else if (SiddhartaSetup == 4 || SiddhartaSetup == 5) {
    ghost_x = 200.0*mm;
    ghost_y = 200.0*mm;
    ghost_z = 1.0*mm;
    pos_ghost = -115.*mm;
  } else if (SiddhartaSetup == 6 || SiddhartaSetup == 7 || SiddhartaSetup == 8) {
    ghost_x = 200.0*mm;
    ghost_y = 200.0*mm;
    ghost_z = 10.*um;
    pos_ghost = -130.*mm;
  } else {
    ghost_x = 200.0*mm;
    ghost_y = 200.0*mm;
    ghost_z = 1.0*mm;
    pos_ghost = -125.*mm;
  }
  G4Box *solidGhost = new G4Box("Ghost", 0.5*ghost_x, 0.5*ghost_y, 0.5*ghost_z);
  G4LogicalVolume* logicGhost;

  if (SiddhartaSetup == 1) {
    logicGhost = new G4LogicalVolume(solidGhost, vacuum, "logicVacCh", 0, 0, 0);
    G4VPhysicalVolume* physiGhost  = new G4PVPlacement(rm00, G4ThreeVector(0,0,pos_ghost), logicGhost, "Ghost", logicVacCh, false, 0);
  } else if (SiddhartaSetup == 6 || SiddhartaSetup == 7 || SiddhartaSetup == 8) {
    logicGhost = new G4LogicalVolume(solidGhost, vacuum, "logicVacCh", 0, 0, 0);
    G4VPhysicalVolume* physiGhost  = new G4PVPlacement(rm00, G4ThreeVector(0,0,pos_ghost), logicGhost, "Ghost", logicVacCh, false, 0);
  } else {
    logicGhost = new G4LogicalVolume(solidGhost, vacuum, "logicVacCh", 0, 0, 0);
    G4VPhysicalVolume* physiGhost  = new G4PVPlacement(rm00, G4ThreeVector(0,0,pos_ghost), logicGhost, "Ghost", logicVacCh, false, 0);
  }
  logicGhost->SetSensitiveDetector( ghostSD );
  logicGhost->SetVisAttributes(G4Colour(0.,0.,1.,0.));

// Degrader //
  if (SiddhartaSetup != 6 && SiddhartaSetup != 7 && SiddhartaSetup != 8 && SiddhartaSetup != 2020) {
    G4double degrader_step = 20.0*mm;
    G4double dx_degrader = 6*degrader_step;
    G4double dy_degrader = 6*degrader_step;
    G4double dx1_deg = dx_degrader;
    G4double dy1_deg = dy_degrader;
    G4double dz1_deg;

    if (SiddhartaSetup == 1) {
      dz1_deg = 575.*um;
    } else {
      dz1_deg = 525.*um;//
    }
    G4double pos1_deg;
    pos1_deg = pos_shield + 0.5*Z_SHIEL + dist_shield_degrader + 0.5*dz1_deg;
    G4double dz2_deg;
    G4double dz3_deg;
    G4double dz4_deg;
    G4double dz5_deg;
    G4double dz6_deg;
    dz2_deg = 100.*um;
    dz3_deg = 100.*um;
    dz4_deg = 200.*um;
    dz5_deg = 100.*um;
    dz6_deg = 100.*um;

    G4double dx2_deg = dx_degrader - degrader_step;
    G4double dy2_deg = dy_degrader;
    G4double pos2_deg = pos1_deg + 0.5*dz1_deg + 0.5*dz2_deg;
    G4double posx_deg2 = 0.5*degrader_step;
    G4double dx3_deg = dx_degrader - 2.*degrader_step;
    G4double dy3_deg = dy_degrader;
    G4double pos3_deg = pos2_deg + 0.5*dz2_deg + 0.5*dz3_deg;
    G4double posx_deg3 = 0.5*2.*degrader_step;
    G4double dx4_deg = dx_degrader - 3.*degrader_step;
    G4double dy4_deg = dy_degrader;
    G4double pos4_deg = pos3_deg + 0.5*dz3_deg + 0.5*dz4_deg;
    G4double posx_deg4 = 0.5*3.*degrader_step;
    G4double dx5_deg = dx_degrader - 4.*degrader_step;
    G4double dy5_deg = dy_degrader;
    G4double pos5_deg = pos4_deg + 0.5*dz4_deg + 0.5*dz5_deg;
    G4double posx_deg5 = 0.5*4.*degrader_step;
    G4double dx6_deg = dx_degrader - 5.*degrader_step;
    G4double dy6_deg = dy_degrader;
    G4double pos6_deg = pos5_deg + 0.5*dz5_deg + 0.5*dz6_deg;
    G4double posx_deg6 = 0.5*5.*degrader_step;

    G4Box* solidDegraderBase = new G4Box("DegraderBase", 0.5*dx1_deg, 0.5*dy1_deg, 0.5*dz1_deg);
    G4LogicalVolume* logicDegraderBase = new G4LogicalVolume(solidDegraderBase, mylar, "DegraderBase", 0, 0, 0);
    G4VPhysicalVolume* physiDegraderBase = new G4PVPlacement(rm01, G4ThreeVector(0.,pos1_deg,0.), logicDegraderBase,
                                                             "DegraderBase", World, false, 1);
    logicDegraderBase->SetVisAttributes(G4Colour(0.,1.,0.,0.));
    G4Box* solidDegrader2 = new G4Box("Degrader2", 0.5*dx2_deg, 0.5*dy2_deg, 0.5*dz2_deg);
    G4LogicalVolume* logicDegrader2 = new G4LogicalVolume(solidDegrader2, mylar, "Degrader2", 0, 0, 0);
    G4VPhysicalVolume* physiDegrader2 = new G4PVPlacement(rm01, G4ThreeVector(posx_deg2,pos2_deg,0.), logicDegrader2,
                                                          "Degrader2", World, false, 1);
    logicDegrader2->SetVisAttributes(G4Colour(1.,0.,0.,0.));
    G4Box* solidDegrader3 = new G4Box("Degrader3", 0.5*dx3_deg, 0.5*dy3_deg, 0.5*dz3_deg);
    G4LogicalVolume* logicDegrader3 = new G4LogicalVolume(solidDegrader3, mylar, "Degrader3", 0, 0, 0);
    G4VPhysicalVolume* physiDegrader3 = new G4PVPlacement(rm01, G4ThreeVector(posx_deg3,pos3_deg,0.), logicDegrader3,
                                                          "Degrader3", World, false, 1);
    logicDegrader3->SetVisAttributes(G4Colour(0.,0.,1.,0.));
    G4Box* solidDegrader4 = new G4Box("Degrader4", 0.5*dx4_deg, 0.5*dy4_deg, 0.5*dz4_deg);
    G4LogicalVolume* logicDegrader4 = new G4LogicalVolume(solidDegrader4, mylar, "Degrader4", 0, 0, 0);
    G4VPhysicalVolume* physiDegrader4 = new G4PVPlacement(rm01, G4ThreeVector(posx_deg4,pos4_deg,0.), logicDegrader4,
                                                          "Degrader4", World, false, 1);
    logicDegrader4->SetVisAttributes(G4Colour(1.,1.,0.,0.));
    G4Box* solidDegrader5 = new G4Box("Degrader5", 0.5*dx5_deg, 0.5*dy5_deg, 0.5*dz5_deg);
    G4LogicalVolume* logicDegrader5 = new G4LogicalVolume(solidDegrader5, mylar, "Degrader5", 0, 0, 0);
    G4VPhysicalVolume* physiDegrader5 = new G4PVPlacement(rm01, G4ThreeVector(posx_deg5,pos5_deg,0.), logicDegrader5,
                                                          "Degrader5", World, false, 1);
    logicDegrader5->SetVisAttributes(G4Colour(1.,0.,0.,0.));
    G4Box* solidDegrader6 = new G4Box("Degrader6", 0.5*dx6_deg, 0.5*dy6_deg, 0.5*dz6_deg);
    G4LogicalVolume* logicDegrader6 = new G4LogicalVolume(solidDegrader6, mylar, "Degrader6", 0, 0, 0);
    G4VPhysicalVolume* physiDegrader6 = new G4PVPlacement(rm01, G4ThreeVector(posx_deg6,pos6_deg,0.), logicDegrader6,
                                                          "Degrader6", World, false, 1);
    logicDegrader6->SetVisAttributes(G4Colour(0.,0.,1.,0.));
  } else {
    G4double degrader_step = 10.0*mm;
    G4double dx_degrader = 12*degrader_step;
    G4double dy_degrader = 12*degrader_step;
    G4double dx1_deg = dx_degrader;
    G4double dy1_deg = dy_degrader;

    G4double dz1_deg;
    G4double dz2_deg;
    G4double dz3_deg;
    G4double dz4_deg;
    G4double dz5_deg;
    G4double dz6_deg;
    G4double dz7_deg;
    G4double dz8_deg;
    G4double dz9_deg;
    G4double dz10_deg;
    G4double dz11_deg;
    G4double dz12_deg;

    if (SiddhartaSetup == 2020) {
      dz1_deg = 0.1*um;
      dz2_deg = 0.1*um;
      dz3_deg = 200.*um;
      dz4_deg = 100.*um;
      dz5_deg = 100.*um;
      dz6_deg = 150.*um;
      dz7_deg = 100.*um;
      dz8_deg = 100.*um;
      dz9_deg = 100.*um;
      dz10_deg = 100.*um;
      dz11_deg = 0.1*um;
      dz12_deg = 0.1*um;
    } else {
      dz1_deg = 210.*um;
      dz2_deg = 125.*um;
      dz3_deg = 125.*um;
      dz4_deg = 110.*um;
      dz5_deg = 105.*um;
      dz6_deg = 95.*um;
      dz7_deg = 140.*um;
      dz8_deg = 80.*um;
      dz9_deg = 75.*um;
      dz10_deg = 65.*um;
      dz11_deg = 55.*um;
      dz12_deg = 45.*um;
    }

    if (degThickFromCardOption) {
      dz1_deg = fDegModif->GetThickAtLayer(1);
      dz2_deg = fDegModif->GetThickAtLayer(2);
      dz3_deg = fDegModif->GetThickAtLayer(3);
      dz4_deg = fDegModif->GetThickAtLayer(4);
      dz5_deg = fDegModif->GetThickAtLayer(5);
      dz6_deg = fDegModif->GetThickAtLayer(6);
      dz7_deg = fDegModif->GetThickAtLayer(7);
      dz8_deg = fDegModif->GetThickAtLayer(8);
      dz9_deg = fDegModif->GetThickAtLayer(9);
      dz10_deg = fDegModif->GetThickAtLayer(10);
      dz11_deg = fDegModif->GetThickAtLayer(11);
      dz12_deg = fDegModif->GetThickAtLayer(12);
    }

    G4double pos1_deg;
    pos1_deg = pos_shield + 0.5*Z_SHIEL + dist_shield_degrader + 0.5*dz1_deg;
    G4double dx2_deg = dx_degrader - degrader_step;
    G4double dy2_deg = dy_degrader;
    G4double pos2_deg = pos1_deg + 0.5*dz1_deg + 0.5*dz2_deg;
    G4double posx_deg2 = 0.5*degrader_step;
    G4double dx3_deg = dx_degrader - 2.*degrader_step;
    G4double dy3_deg = dy_degrader;
    G4double pos3_deg = pos2_deg + 0.5*dz2_deg + 0.5*dz3_deg;
    G4double posx_deg3 = 0.5*2.*degrader_step;
    G4double dx4_deg = dx_degrader - 3.*degrader_step;
    G4double dy4_deg = dy_degrader;
    G4double pos4_deg = pos3_deg + 0.5*dz3_deg + 0.5*dz4_deg;
    G4double posx_deg4 = 0.5*3.*degrader_step;
    G4double dx5_deg = dx_degrader - 4.*degrader_step;
    G4double dy5_deg = dy_degrader;
    G4double pos5_deg = pos4_deg + 0.5*dz4_deg + 0.5*dz5_deg;
    G4double posx_deg5 = 0.5*4.*degrader_step;
    G4double dx6_deg = dx_degrader - 5.*degrader_step;
    G4double dy6_deg = dy_degrader;
    G4double pos6_deg = pos5_deg + 0.5*dz5_deg + 0.5*dz6_deg;
    G4double posx_deg6 = 0.5*5.*degrader_step;
    G4double dx7_deg = dx_degrader - 6.*degrader_step;
    G4double dy7_deg = dy_degrader;
    G4double pos7_deg = pos6_deg + 0.5*dz6_deg + 0.5*dz7_deg;
    G4double posx_deg7 = 0.5*6.*degrader_step;
    G4double dx8_deg = dx_degrader - 7.*degrader_step;
    G4double dy8_deg = dy_degrader;
    G4double pos8_deg = pos7_deg + 0.5*dz7_deg + 0.5*dz8_deg;
    G4double posx_deg8 = 0.5*7.*degrader_step;
    G4double dx9_deg = dx_degrader - 8.*degrader_step;
    G4double dy9_deg = dy_degrader;
    G4double pos9_deg = pos8_deg + 0.5*dz8_deg + 0.5*dz9_deg;
    G4double posx_deg9 = 0.5*8.*degrader_step;
    G4double dx10_deg = dx_degrader - 9.*degrader_step;
    G4double dy10_deg = dy_degrader;
    G4double pos10_deg = pos9_deg + 0.5*dz9_deg + 0.5*dz10_deg;
    G4double posx_deg10 = 0.5*9.*degrader_step;
    G4double dx11_deg = dx_degrader - 10.*degrader_step;
    G4double dy11_deg = dy_degrader;
    G4double pos11_deg = pos10_deg + 0.5*dz10_deg + 0.5*dz11_deg;
    G4double posx_deg11 = 0.5*10.*degrader_step;
    G4double dx12_deg = dx_degrader - 11.*degrader_step;
    G4double dy12_deg = dy_degrader;
    G4double pos12_deg = pos11_deg + 0.5*dz11_deg + 0.5*dz12_deg;
    G4double posx_deg12 = 0.5*11.*degrader_step;

    G4Box* solidDegraderBase = new G4Box("DegraderBase", 0.5*dx1_deg, 0.5*dy1_deg, 0.5*dz1_deg);
    G4LogicalVolume* logicDegraderBase = new G4LogicalVolume(solidDegraderBase, mylar, "DegraderBase", 0, 0, 0);
    G4VPhysicalVolume* physiDegraderBase = new G4PVPlacement(rm01, G4ThreeVector(0.,pos1_deg,0.), logicDegraderBase,
                                                             "DegraderBase", World, false, 1);
    logicDegraderBase->SetVisAttributes(G4Colour(1.,0.8,0.3,1.));
    G4Box* solidDegrader2 = new G4Box("Degrader2", 0.5*dx2_deg, 0.5*dy2_deg, 0.5*dz2_deg);
    G4LogicalVolume* logicDegrader2 = new G4LogicalVolume(solidDegrader2, mylar, "Degrader2", 0, 0, 0);
    G4VPhysicalVolume* physiDegrader2 = new G4PVPlacement(rm01, G4ThreeVector(posx_deg2,pos2_deg,0.), logicDegrader2,
                                                          "Degrader2", World, false, 1);
    logicDegrader2->SetVisAttributes(G4Colour(1.,0.8,0.3,1.));
    G4Box* solidDegrader3 = new G4Box("Degrader3", 0.5*dx3_deg, 0.5*dy3_deg, 0.5*dz3_deg);
    G4LogicalVolume* logicDegrader3 = new G4LogicalVolume(solidDegrader3, mylar, "Degrader3", 0, 0, 0);
    G4VPhysicalVolume* physiDegrader3 = new G4PVPlacement(rm01, G4ThreeVector(posx_deg3,pos3_deg,0.), logicDegrader3,
                                                          "Degrader3", World, false, 1);
    logicDegrader3->SetVisAttributes(G4Colour(1.,0.8,0.3,1.));
    G4Box* solidDegrader4 = new G4Box("Degrader4", 0.5*dx4_deg, 0.5*dy4_deg, 0.5*dz4_deg);
    G4LogicalVolume* logicDegrader4 = new G4LogicalVolume(solidDegrader4, mylar, "Degrader4", 0, 0, 0);
    G4VPhysicalVolume* physiDegrader4 = new G4PVPlacement(rm01, G4ThreeVector(posx_deg4,pos4_deg,0.), logicDegrader4,
                                                          "Degrader4", World, false, 1);
    logicDegrader4->SetVisAttributes(G4Colour(1.,0.8,0.3,1.));
    G4Box* solidDegrader5 = new G4Box("Degrader5", 0.5*dx5_deg, 0.5*dy5_deg, 0.5*dz5_deg);
    G4LogicalVolume* logicDegrader5 = new G4LogicalVolume(solidDegrader5, mylar, "Degrader5", 0, 0, 0);
    G4VPhysicalVolume* physiDegrader5 = new G4PVPlacement(rm01, G4ThreeVector(posx_deg5,pos5_deg,0.), logicDegrader5,
                                                          "Degrader5", World, false, 1);
    logicDegrader5->SetVisAttributes(G4Colour(1.,0.8,0.3,1.));
    G4Box* solidDegrader6 = new G4Box("Degrader6", 0.5*dx6_deg, 0.5*dy6_deg, 0.5*dz6_deg);
    G4LogicalVolume* logicDegrader6 = new G4LogicalVolume(solidDegrader6, mylar, "Degrader6", 0, 0, 0);
    G4VPhysicalVolume* physiDegrader6 = new G4PVPlacement(rm01, G4ThreeVector(posx_deg6,pos6_deg,0.), logicDegrader6,
                                                          "Degrader6", World, false, 1);
    logicDegrader6->SetVisAttributes(G4Colour(1.,0.8,0.3,1.));
    G4Box* solidDegrader7 = new G4Box("Degrader7", 0.5*dx7_deg, 0.5*dy7_deg, 0.5*dz7_deg);
    G4LogicalVolume* logicDegrader7 = new G4LogicalVolume(solidDegrader7, mylar, "Degrader7", 0, 0, 0);
    G4VPhysicalVolume* physiDegrader7 = new G4PVPlacement(rm01, G4ThreeVector(posx_deg7,pos7_deg,0.), logicDegrader7,
                                                          "Degrader7", World, false, 1);
    logicDegrader7->SetVisAttributes(G4Colour(1.,0.8,0.3,1.));
    G4Box* solidDegrader8 = new G4Box("Degrader8", 0.5*dx8_deg, 0.5*dy8_deg, 0.5*dz8_deg);
    G4LogicalVolume* logicDegrader8 = new G4LogicalVolume(solidDegrader8, mylar, "Degrader8", 0, 0, 0);
    G4VPhysicalVolume* physiDegrader8 = new G4PVPlacement(rm01, G4ThreeVector(posx_deg8,pos8_deg,0.), logicDegrader8,
                                                          "Degrader8", World, false, 1);
    logicDegrader8->SetVisAttributes(G4Colour(1.,0.8,0.3,1.));
    G4Box* solidDegrader9 = new G4Box("Degrader9", 0.5*dx9_deg, 0.5*dy9_deg, 0.5*dz9_deg);
    G4LogicalVolume* logicDegrader9 = new G4LogicalVolume(solidDegrader9, mylar, "Degrader9", 0, 0, 0);
    G4VPhysicalVolume* physiDegrader9 = new G4PVPlacement(rm01, G4ThreeVector(posx_deg9,pos9_deg,0.), logicDegrader9,
                                                          "Degrader9", World, false, 1);
    logicDegrader9->SetVisAttributes(G4Colour(1.,0.8,0.3,1.));
    G4Box* solidDegrader10 = new G4Box("Degrader10",0.5*dx10_deg,0.5*dy10_deg,0.5*dz10_deg);
    G4LogicalVolume* logicDegrader10 = new G4LogicalVolume(solidDegrader10, mylar, "Degrader10", 0, 0, 0);
    G4VPhysicalVolume* physiDegrader10 = new G4PVPlacement(rm01, G4ThreeVector(posx_deg10,pos10_deg,0.), logicDegrader10,
                                                           "Degrader10", World, false, 1);
    logicDegrader10->SetVisAttributes(G4Colour(1.,0.8,0.3,1.));
    G4Box* solidDegrader11 = new G4Box("Degrader11",0.5*dx11_deg,0.5*dy11_deg,0.5*dz11_deg);
    G4LogicalVolume* logicDegrader11 = new G4LogicalVolume(solidDegrader11, mylar, "Degrader11", 0, 0, 0);
    G4VPhysicalVolume* physiDegrader11 = new G4PVPlacement(rm01, G4ThreeVector(posx_deg11,pos11_deg,0.), logicDegrader11,
                                                           "Degrader11", World, false, 1);
    logicDegrader11->SetVisAttributes(G4Colour(1.,0.8,0.3,1.));
    G4Box* solidDegrader12 = new G4Box("Degrader12",0.5*dx12_deg,0.5*dy12_deg,0.5*dz12_deg);
    G4LogicalVolume* logicDegrader12 = new G4LogicalVolume(solidDegrader12, mylar, "Degrader12", 0, 0, 0);
    G4VPhysicalVolume* physiDegrader12 = new G4PVPlacement(rm01, G4ThreeVector(posx_deg12,pos12_deg,0.), logicDegrader12,
                                                           "Degrader12", World, false, 1);
    logicDegrader12->SetVisAttributes(G4Colour(1.,0.8,0.3,1.));
  }

  if (SiddhartaSetup == 1) {
    G4Box* solidCriogenics = new G4Box("Criogenics", 0.5*50*mm, 0.5*20*mm, 0.5*50*mm);
    G4LogicalVolume* logicCriogenics = new G4LogicalVolume(solidCriogenics, Al, "Criogenics", 0, 0, 0);
    G4VPhysicalVolume* physiCriogenics = new G4PVPlacement(rm00, G4ThreeVector(0,650*mm,0), logicCriogenics,
                                                           "Criogenics", World, false, 0);
    logicCriogenics->SetVisAttributes(G4Colour(0.,0.,1.,0.));
  } else if (SiddhartaSetup == 7) {
    G4Box* solidCriogenics = new G4Box("Criogenics", 0.5*800*mm, 0.5*200*mm, 0.5*800*mm);
    G4LogicalVolume* logicCriogenics = new G4LogicalVolume(solidCriogenics, Al, "Criogenics", 0, 0, 0);
    G4VPhysicalVolume* physiCriogenics = new G4PVPlacement(rm00, G4ThreeVector(0,650*mm,0), logicCriogenics,
                                                           "Criogenics", World, false, 0);
    logicCriogenics->SetVisAttributes(G4Colour(0.,0.,1.,0.));
  } else if (SiddhartaSetup == 8) {
    G4Box* solidCriogenics = new G4Box("Criogenics", 0.5*80*mm, 0.5*20*mm, 0.5*80*mm);
    G4LogicalVolume* logicCriogenics = new G4LogicalVolume(solidCriogenics, Al, "Criogenics", 0, 0, 0);
    G4VPhysicalVolume* physiCriogenics = new G4PVPlacement(rm00, G4ThreeVector(0,650*mm,0), logicCriogenics,
                                                           "Criogenics", World, false, 0);
    logicCriogenics->SetVisAttributes(G4Colour(0.,0.,1.,0.));
  }
  G4cout << "Scintillators X " << 3 << G4endl;
  return physiWorld;
}

void SiddhartaDetectorConstruction::setTargetMaterial(G4String materialName)
{
  G4Material* pttoMaterial = G4Material::GetMaterial(materialName);
  if (pttoMaterial) {
    TargetMater = pttoMaterial;
    logicTarget->SetMaterial(pttoMaterial);
    G4cout << "\n----> The target is " << fTargetLength/cm << " cm of " << materialName << G4endl;
  }
}

void SiddhartaDetectorConstruction::SetMagField(G4double fieldValue)
{
  fpMagField->SetMagFieldValue(fieldValue);
}

void SiddhartaDetectorConstruction::SetMaxStep(G4double maxStep)
{
  if ((stepLimit) && (maxStep > 0.))
    stepLimit->SetMaxAllowedStep(maxStep);
}
