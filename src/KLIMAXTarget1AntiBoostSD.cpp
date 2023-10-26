#include "../include/KLIMAXTarget1AntiBoostSD.h"
#include "../include/SiddhartaAnalysisManager.h"
#include "../include/SiddhartaHisto.h"

#include <G4ThreeVector.hh>
#include <G4SDManager.hh>
#include <G4VSolid.hh>
#include <G4Step.hh>
#include <G4ios.hh>

KLIMAXTarget1AntiBoostSD::KLIMAXTarget1AntiBoostSD(G4String name) : G4VSensitiveDetector(name)
{
  G4String HCname;
  collectionName.insert(HCname = "trackerCollectionKLT1AntiBoost");
}

KLIMAXTarget1AntiBoostSD::~KLIMAXTarget1AntiBoostSD(){ }

void KLIMAXTarget1AntiBoostSD::Initialize(G4HCofThisEvent* HCE)
{
  trackerCollection = new SiddhartaTrackerHitsCollection(SensitiveDetectorName,collectionName[0]);
  static G4int HCID = -1;

  if(HCID<0) {
    HCID = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]);
  }
  HCE->AddHitsCollection(HCID, trackerCollection);

  SiddhartaAnalysisManager* analysis = SiddhartaAnalysisManager::getInstance();

  if (analysis->YesHistos) {
    sciEnergy = 0.0;
    sciEnergyMax = 0.0;
    TimeKLT1AntiBoost = -9999.;
    X = -99999.;
    Y = -99999.;
    Z = -99999.;

    analysis->histo->ntuData.KLT1AntiBooststop[0] = -1000000 / mm;
    analysis->histo->ntuData.KLT1AntiBooststop[1] = -1000000 /  mm;
    analysis->histo->ntuData.KLT1AntiBooststop[2] = -1000000 / mm;
    analysis->histo->ntuData.EnergyDepKLT1AntiBoost = -1000000 / eV;
    analysis->histo->ntuData.TimeKLT1AntiBoost = -1000000 / ns;
    analysis->histo->ntuData.XYZKLT1AntiBoost[0] = -1000000 / mm;
    analysis->histo->ntuData.XYZKLT1AntiBoost[1] = -1000000 / mm;
    analysis->histo->ntuData.XYZKLT1AntiBoost[2] = -1000000 / mm;
    analysis->histo->ntuData.KLT1AntiBoostKaonstop[0] = -1000000 / mm; 
    analysis->histo->ntuData.KLT1AntiBoostKaonstop[1] = -1000000 / mm;
    analysis->histo->ntuData.KLT1AntiBoostKaonstop[2] = -1000000 / mm;
    analysis->histo->ntuData.kaonKinEKLT1AntiBoost = -1000000 / eV;
    analysis->histo->ntuData.gammaKinEKLT1AntiBoost = -1000000 / eV;
    analysis->histo->ntuData.lastkaonKinEKLT1AntiBoost = -1000000 / eV;
    analysis->histo->ntuData.pdgcodeKLT1AntiBoost = -1000000;
  }
}

G4bool KLIMAXTarget1AntiBoostSD::ProcessHits(G4Step* aStep,G4TouchableHistory*)
{
  G4double edep = aStep->GetTotalEnergyDeposit();

  if (edep == 0.)
    return false;

  SiddhartaTrackerHit* newHit = new SiddhartaTrackerHit();
  newHit->SetTrackID(aStep->GetTrack()->GetTrackID());
  newHit->SetEdep(edep);
  newHit->SetPos(aStep->GetPostStepPoint()->GetPosition());
  trackerCollection->insert(newHit);

  G4String pname = "";
  pname = aStep->GetTrack()->GetDefinition()->GetParticleName();
  G4String prematname = "";
  G4String prevolname = "";
  G4String postmatname = "";
  G4String postvolname = "";

  G4StepStatus prestepState = aStep->GetPreStepPoint()->GetStepStatus();
  G4StepStatus poststepState = aStep->GetPostStepPoint()->GetStepStatus();
  prematname = aStep->GetPreStepPoint()->GetMaterial()->GetName();
  prevolname = aStep->GetPreStepPoint()->GetPhysicalVolume()->GetName();
  postmatname = aStep->GetPostStepPoint()->GetMaterial()->GetName();
  postvolname = aStep->GetPostStepPoint()->GetPhysicalVolume()->GetName();
  G4double Kinener = aStep->GetTrack()->GetMomentum().mag();
  G4ThreeVector pos = aStep->GetTrack()->GetPosition();
  G4String preproc = "undefined";
  G4String postproc = "undefined";

  SiddhartaAnalysisManager* analysis = SiddhartaAnalysisManager::getInstance();

  sciEnergy += edep;

  if (analysis->YesHistos) {
    analysis->histo->ntuData.pdgcodeKLT1AntiBoost = aStep->GetTrack()->GetDynamicParticle()->GetPDGcode();
    X = (aStep->GetTrack()->GetPosition())[0] ;
    Y = (aStep->GetTrack()->GetPosition())[1] ;
    Z = (aStep->GetTrack()->GetPosition())[2] ;

    if (edep > sciEnergyMax) {
      sciEnergyMax = edep;
      TimeKLT1AntiBoost = aStep->GetTrack()->GetGlobalTime();
      pname = aStep->GetTrack()->GetDefinition()->GetParticleName();

      analysis->histo->ntuData.TimeKLT1AntiBoost = TimeKLT1AntiBoost / ns;
      analysis->histo->ntuData.XYZKLT1AntiBoost[0] = X / mm;
      analysis->histo->ntuData.XYZKLT1AntiBoost[1] = Y / mm;
      analysis->histo->ntuData.XYZKLT1AntiBoost[2] = Z / mm;
    }
    analysis->histo->ntuData.KLT1AntiBooststop[0] = X / mm;
    analysis->histo->ntuData.KLT1AntiBooststop[1] = Y / mm;
    analysis->histo->ntuData.KLT1AntiBooststop[2] = Z / mm;

    if(analysis->histo->ntuData.pdgcodeKLT1AntiBoost == -321) { //K-//
      G4cout << "K- in Target 1 " << (aStep->GetTrack()->GetKineticEnergy()) / eV << G4endl;
      kaonCounter ++;

      if( kaonCounter == 1 )
        analysis->histo->ntuData.kaonKinEKLT1AntiBoost = (aStep->GetTrack()->GetKineticEnergy()) / eV; // get the Energy of the first kaon hit

      analysis->histo->ntuData.lastkaonKinEKLT1AntiBoost = (aStep->GetTrack()->GetKineticEnergy()) / eV; // get the Energy of the last kaon hit
      analysis->histo->ntuData.KLT1AntiBoostKaonstop[0] = X / mm ;
      analysis->histo->ntuData.KLT1AntiBoostKaonstop[1] = Y / mm ;
      analysis->histo->ntuData.KLT1AntiBoostKaonstop[2] = Z / mm ;
    }
    
    if (analysis->histo->ntuData.pdgcodeKLT1AntiBoost == 22)
      analysis->histo->ntuData.gammaKinEKLT1AntiBoost = (aStep->GetTrack()->GetKineticEnergy()) / eV; // get the Energy of the last photon hit
  }

  return true;
}

void KLIMAXTarget1AntiBoostSD::EndOfEvent(G4HCofThisEvent*)
{
  kaonCounter = 0;
  kaonKinE = 0;
}
