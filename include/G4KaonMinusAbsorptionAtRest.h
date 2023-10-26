//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************

#ifndef G4KaonMinusAbsorptionAtRest_h
#define G4KaonMinusAbsorptionAtRest_h 1

// Class Description:
//
// Process for nuclear absorption of K- at rest.
// To be used in your physics list in case you need this physics.

//To remove after applying BertiniAbsorption
#include <G4Fragment.hh>
#include <G4ThreeVector.hh>
#include <G4LorentzVector.hh>
#include <G4ReactionProduct.hh>
#include <G4ExcitationHandler.hh>
#include <G4ReactionProductVector.hh>
//To remove after applying BertiniAbsorption

//To remove after applying G4PhaseSpaceDecayChannel
#include <G4PhysicalConstants.hh>
#include <G4SystemOfUnits.hh>
//To remove after applying G4PhaseSpaceDecayChannel

#include <G4DynamicParticleVector.hh>
#include <G4HadronicProcessType.hh>
#include <G4NucleiProperties.hh>
#include <G4DynamicParticle.hh>
#include <G4ParticleTypes.hh>
#include <G4VRestProcess.hh>
#include <Randomize.hh>
#include <G4Nucleus.hh>
#include <globals.hh>

// *********************************************************
class G4KaonMinusAbsorptionAtRest : public G4VRestProcess
// *********************************************************
{
public:
  G4KaonMinusAbsorptionAtRest(const G4String& processName = "KaonMinusAbsorptionAtRest", G4ProcessType aType = fHadronic);
  ~G4KaonMinusAbsorptionAtRest();

  G4bool IsApplicable(const G4ParticleDefinition& particle) {
    return (particle == *(G4KaonMinus::KaonMinus()));
  }

  void PreparePhysicsTable(const G4ParticleDefinition&);
  void BuildPhysicsTable(const G4ParticleDefinition&);
  G4VParticleChange* AtRestDoIt(const G4Track& aTrack, const G4Step& aStep);

  void GenerateKaonNucleusAbsLines();
  
private:
  G4KaonMinusAbsorptionAtRest& operator=(const G4KaonMinusAbsorptionAtRest &right);
  G4KaonMinusAbsorptionAtRest(const G4KaonMinusAbsorptionAtRest& );

  G4double GetMeanLifeTime(const G4Track& aTrack, G4ForceCondition*) {return 0;}
// returns proton or neutron with fermi-momentum
  G4DynamicParticle GetAbsorbingNucleon();
// returns proton or neutron particle definition;
  G4ParticleDefinition* SelectAbsorbingNucleon();
// provides the neutron halo factor for absorption on nucleus surface.
  G4double NeutronHaloFactor(G4double Z, G4double N);
//  creates the reaction products
  G4DynamicParticleVector* KaonNucleonReaction();
// secondary pion absorption in parent nucleus
// if TRUE, then add excitation energy to the Nucleus
  G4bool AbsorbPionByNucleus(G4DynamicParticle* aPion);
//  secondary Sigma-Lambda conversion
// if conversion Done, then add excitation energy to the Nucleus
  G4DynamicParticle *SigmaLambdaConversion(G4DynamicParticle* aSigma);

  G4Material* currentMaterial;
  const G4DynamicParticle *stoppedHadron;
  G4Nucleus* nucleus;

  G4double pionAbsorptionRate;
  G4double rateLambdaZeroPiZero;
  G4double rateLambdaZeroPiMinus;
  G4double rateSigmaMinusPiPlus;
  G4double rateSigmaMinusPiZero;
  G4double rateSigmaPlusPiMinus;
  G4double rateSigmaZeroPiZero;
  G4double rateSigmaZeroPiMinus;

// Sigma Lambda Conversion rates
// for sigma- p -> lambda n
//     sigma+ n -> lambda p
//     sigma- n -> lambda

  G4double sigmaPlusLambdaConversionRate;
  G4double sigmaMinusLambdaConversionRate;
  G4double sigmaZeroLambdaConversionRate;

// Atomic cascade with Xrays
  G4double yieldKa,probKa,yieldKb,probKb,yieldKg,probKg,yieldKd,probKd,yieldKe,probKe,yieldKz,probKz,yieldKi,probKi;
  G4double yieldXaOther,probXaOther;

  std::map<G4int, unsigned> atomicNumberVSnucleusRef;
  std::vector<KaonNucleusAbsLines> nucleiImplemented;
  
  bool fPlotProducedXrayLies = true;
};

//------------------Implementing structure for nuclei absorption lines
enum AbsLineSimType {absSimulateAll, absSimulateByIntensity};

class KaonNucleusAbsLines
{
public:
  KaonNucleusAbsLines() {nmbOfIsotopes = 0; atomicNumber = 0; cleanMaterialsOnly = false;};
  ~KaonNucleusAbsLines() {nmbOfIsotopes = 0; atomicNumber = 0; nucleusNumber.clear(); photonEnergy.clear(); intensityOfGivenEnergy.clear();};

  void AddNucleusNumber(G4int newA) {nucleusNumber.push_back(newA);};
  void AddPhotonEnergy(G4double energy, G4double intensity) {
    photonEnergy.push_back(energy);
    intensityOfGivenEnergy.push_back(intensity);
  };
  
  void SetAtomicNumber(G4int newZ) {atomicNumber = newZ;};
  void SetSimType(AbsLineSimType type) {typeOfSimulation = type;};
  void SetCleanMaterialOnly() {cleanMaterialsOnly = true;};
  
  bool CheckNucleusNumber(G4int givenZ);
  
  unsigned GetNumberOfIsotopes() {return nmbOfIsotopes;};
  G4int GetAtomicNumber() {return atomicNumber;};
  std::vector<G4double> GetPhotonEnergies(G4int noOfElements);

private:
  unsigned nmbOfIsotopes;
  G4int atomicNumber; //number of protons
  std::vector<G4int> nucleusNumber;
  std::vector<G4double> photonEnergy;
  std::vector<G4double> intensityOfGivenEnergy;
  AbsLineSimType typeOfSimulation;
  bool cleanMaterialsOnly;
}

bool KaonNucleusAbsLines::CheckNucleusNumber(G4int givenZ)
{
  for (unsigned i=0; i<nucleusNumber.size(); i++) {
    if (givenZ == nucleusNumber.at(i))
      return true;
  }
  return false;
}

std::vector<G4double> KaonNucleusAbsLines::GetPhotonEnergies(G4int noOfElements)
{
  std::vector<G4double> photonEnergyToReturn;
  
  if (cleanMaterialsOnly && noOfElements > 1)
    return photonEnergyToReturn;
  
  if (typeOfSimulation == AbsLineSimType::absSimulateAll) {
    photonEnergyToReturn = photonEnergy;
  } else {
    for unsigned i=0; i<intensityOfGivenEnergy.size() && i<photonEnergyToReturn.size(); i++) {
      if (intensityOfGivenEnergy.at(i) > G4UniformRand())
        photonEnergyToReturn.push_back(photonEnergy.at(i));
    }
  }
  return photonEnergyToReturn;
}
//--------------------------------------------------------------------

//To remove after implementing G4PhaseSpaceDecayChannel
//--
class G4ReactionKinematics {
public:
  void TwoBodyScattering(const G4DynamicParticle* pIn1, const G4DynamicParticle* pIn2,
                         G4DynamicParticle* pOut1, G4DynamicParticle* pOut2)
  {
    G4LorentzVector sumIn(pIn1->Get4Momentum() + pIn2->Get4Momentum());
    G4double invariantMass = sumIn.mag();
    G4ThreeVector betaCMS = sumIn.boostVector();
    G4double massOut1 = pOut1->GetMass();
    G4double massOut2 = pOut2->GetMass();
    G4double breakupMomentum = BreakupMomentum(invariantMass, massOut1, massOut2);
    G4double costheta = 2.0*G4UniformRand() - 1.0;
    G4double sintheta = std::sqrt(1.0 - costheta*costheta);
    G4double phi = 2.0*pi*G4UniformRand();
    G4double pz = costheta*breakupMomentum;
    G4double px = sintheta*std::cos(phi)*breakupMomentum;
    G4double py = sintheta*std::sin(phi)*breakupMomentum;
    G4double breakupMomentumSquared = breakupMomentum*breakupMomentum;
    G4double energy1 = std::sqrt(breakupMomentumSquared + massOut1*massOut1);
    G4double energy2 = std::sqrt(breakupMomentumSquared + massOut2*massOut2);
    G4LorentzVector lorentz1(px, py, pz, energy1);
    G4LorentzVector lorentz2(px, py, pz, energy2);
    lorentz1.boost(betaCMS);
    lorentz2.boost(betaCMS);
    pOut1->Set4Momentum(lorentz1);
    pOut2->Set4Momentum(lorentz2);
    return;
  }
  inline G4double BreakupMomentum(G4double totalMass, G4double m1, G4double m2);
};

inline G4double G4ReactionKinematics::BreakupMomentum(G4double totalMass, G4double massA, G4double massB)
{
  G4double m0squared = totalMass*totalMass;
  G4double breakupMomentumSquared = (m0squared - (massA + massB)*(massA + massB))*(m0squared - (massA - massB)*(massA - massB))/(4*m0squared);
  if (breakupMomentumSquared > 0)
    return std::sqrt(breakupMomentumSquared);
  else
    return -1.;
}
//--
//To remove after implementing G4PhaseSpaceDecayChannel



//To remove after applying BertiniAbsorption
//--
class G4StopDeexcitationAlgorithm
{
public:
  G4StopDeexcitationAlgorithm() {};
  virtual ~G4StopDeexcitationAlgorithm() {};
  virtual G4ReactionProductVector* BreakUp(G4int A, G4int Z, G4double excitation, const G4ThreeVector& p) = 0;
private:
  G4StopDeexcitationAlgorithm& operator=(const G4StopDeexcitationAlgorithm &right);
  G4StopDeexcitationAlgorithm(const G4StopDeexcitationAlgorithm&);
};


class G4StopTheoDeexcitation: public G4StopDeexcitationAlgorithm
{
public:
  G4StopTheoDeexcitation() {};
  virtual ~G4StopTheoDeexcitation() {};
  virtual G4ReactionProductVector* BreakUp(G4int A, G4int Z, G4double excitation, const G4ThreeVector& p)
  {
    G4ExcitationHandler theHandler;
    theHandler.SetMinEForMultiFrag(300*GeV);
    G4double atomicMass = G4NucleiProperties::GetNuclearMass(static_cast<G4int>(A),static_cast<G4int>(Z));
    G4double mass = atomicMass + excitation;
    G4double pMag = p.mag();
    G4LorentzVector initialMomentum(p.x(),p.y(),p.z(),std::sqrt(pMag*pMag + mass*mass));
    G4Fragment theExcitedNucleus(A,Z,initialMomentum);
    return theHandler.BreakItUp(theExcitedNucleus);
  };
private:
  G4StopTheoDeexcitation& operator=(const G4StopTheoDeexcitation &right);
  G4StopTheoDeexcitation(const G4StopTheoDeexcitation&);
};

class G4StopDeexcitation
{
public:
  G4StopDeexcitation(G4StopDeexcitationAlgorithm* algorithm) {fAlgorithm = algorithm;};
  ~G4StopDeexcitation() {delete fAlgorithm;};
  G4ReactionProductVector* DoBreakUp(G4int A, G4int Z, G4double excitation, const G4ThreeVector& p) const
  {
    G4ReactionProductVector* v = 0;
    if (fAlgorithm != 0) {
      v = fAlgorithm->BreakUp(A,Z,excitation,p);
    }
    return v;
  };
private:
  G4StopDeexcitation& operator=(const G4StopDeexcitation &right);
  G4StopDeexcitation(const G4StopDeexcitation&);
  G4StopDeexcitationAlgorithm* fAlgorithm;
};
//--
//To remove after applying BertiniAbsorption



#endif
