#include "../include/KLIMAXTargetBoostSD.h"
#include "../include/SiddhartaAnalysisManager.h"
#include "../include/SiddhartaHisto.h"

#include <G4ThreeVector.hh>
#include <G4SDManager.hh>
#include <G4VSolid.hh>
#include <G4Step.hh>
#include <G4ios.hh>

KLIMAXTargetBoostSD::KLIMAXTargetBoostSD(G4String name) : G4VSensitiveDetector(name)
{
  G4String HCname;
  collectionName.insert(HCname = "trackerCollectionKLTBoost");
}

KLIMAXTargetBoostSD::~KLIMAXTargetBoostSD(){ }

void KLIMAXTargetBoostSD::Initialize(G4HCofThisEvent* HCE)
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
    TimeKLTBoost = -9999.;
    X = -99999.;
    Y = -99999.;
    Z = -99999.;

    analysis->histo->ntuData.KLTBooststop[0] = -1000000 / mm;
    analysis->histo->ntuData.KLTBooststop[1] = -1000000 /  mm;
    analysis->histo->ntuData.KLTBooststop[2] = -1000000 / mm;
    analysis->histo->ntuData.EnergyDepKLTBoost = -1000000 / eV;
    analysis->histo->ntuData.TimeKLTBoost = -1000000 / ns;
    analysis->histo->ntuData.XYZKLTBoost[0] = -1000000 / mm;
    analysis->histo->ntuData.XYZKLTBoost[1] = -1000000 / mm;
    analysis->histo->ntuData.XYZKLTBoost[2] = -1000000 / mm;
    analysis->histo->ntuData.KLTBoostKaonstop[0] = -1000000 / mm; 
    analysis->histo->ntuData.KLTBoostKaonstop[1] = -1000000 / mm;
    analysis->histo->ntuData.KLTBoostKaonstop[2] = -1000000 / mm;
    analysis->histo->ntuData.kaonKinEKLTBoost = -1000000 / eV;
    analysis->histo->ntuData.gammaKinEKLTBoost = -1000000 / eV;
    analysis->histo->ntuData.lastkaonKinEKLTBoost = -1000000 / eV;
    analysis->histo->ntuData.pdgcodeKLTBoost = -1000000;
  }
}

G4bool KLIMAXTargetBoostSD::ProcessHits(G4Step* aStep,G4TouchableHistory*)
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
    analysis->histo->ntuData.pdgcodeKLTBoost = aStep->GetTrack()->GetDynamicParticle()->GetPDGcode();
    X = (aStep->GetTrack()->GetPosition())[0] ;
    Y = (aStep->GetTrack()->GetPosition())[1] ;
    Z = (aStep->GetTrack()->GetPosition())[2] ;

    if (edep > sciEnergyMax) {
      sciEnergyMax = edep;
      TimeKLTBoost = aStep->GetTrack()->GetGlobalTime();
      pname = aStep->GetTrack()->GetDefinition()->GetParticleName();

      analysis->histo->ntuData.TimeKLTBoost = TimeKLTBoost / ns;
      analysis->histo->ntuData.XYZKLTBoost[0] = X / mm;
      analysis->histo->ntuData.XYZKLTBoost[1] = Y / mm;
      analysis->histo->ntuData.XYZKLTBoost[2] = Z / mm;
    }
    analysis->histo->ntuData.KLTBooststop[0] = X / mm;
    analysis->histo->ntuData.KLTBooststop[1] = Y / mm;
    analysis->histo->ntuData.KLTBooststop[2] = Z / mm;

    if(analysis->histo->ntuData.pdgcodeKLTBoost == -321) { //K-//
      G4cout << "K- in Target 1 " << (aStep->GetTrack()->GetKineticEnergy()) / eV << G4endl;
      kaonCounter ++;

      if( kaonCounter == 1 )
        analysis->histo->ntuData.kaonKinEKLTBoost = (aStep->GetTrack()->GetKineticEnergy()) / eV; // get the Energy of the first kaon hit

      analysis->histo->ntuData.lastkaonKinEKLTBoost = (aStep->GetTrack()->GetKineticEnergy()) / eV; // get the Energy of the last kaon hit
      analysis->histo->ntuData.KLTBoostKaonstop[0] = X / mm ;
      analysis->histo->ntuData.KLTBoostKaonstop[1] = Y / mm ;
      analysis->histo->ntuData.KLTBoostKaonstop[2] = Z / mm ;
    }
    
    if (analysis->histo->ntuData.pdgcodeKLTBoost == 22)
      analysis->histo->ntuData.gammaKinEKLTBoost = (aStep->GetTrack()->GetKineticEnergy()) / eV; // get the Energy of the last gamma hit
  }

  return true;
}

void KLIMAXTargetBoostSD::EndOfEvent(G4HCofThisEvent*)
{
  kaonCounter = 0;
  kaonKinE = 0;
}

