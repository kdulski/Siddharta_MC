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

#include "../include/G4KaonMinusAbsorptionAtRest.h"
#include "../include/SiddhartaAnalysisManager.h"
#include "../include/SiddhartaHisto.h"

//#include <G4StopDeexcitationAlgorithm.hh>
//#include <G4StopTheoDeexcitation.hh>
//#include <G4StopDeexcitation.hh>
//
//equivalent to G4PhaseSpaceDecayChannel
//#include <G4ReactionKinematics.hh>

#include <G4HadronicProcessStore.hh>

#define G4RandBreitWigner CLHEP::RandBreitWigner

G4KaonMinusAbsorptionAtRest::G4KaonMinusAbsorptionAtRest(const G4String& processName, G4ProcessType aType) : G4VRestProcess (processName, aType)
{
  if (verboseLevel > 0) {
    G4cout << GetProcessName() << " is created "<< G4endl;
  }
  SetProcessSubType(fHadronAtRest);

  // see Cohn et al, PLB27(1968) 527;
  //     Davis et al, PLB1(1967) 434;
  pionAbsorptionRate = 0.07;

  // see  VanderVelde-Wilquet et al, Nuov.Cim.39A(1978)538;
  // see  VanderVelde-Wilquet et al, Nuov.Cim.38A(1977)178;
  // see  VanderVelde-Wilquet et al, Nucl.Phys.A241(1975)511;
  // primary production rates ( for absorption on Carbon)
  // .. other elements are extrapolated by the halo factor.
  rateLambdaZeroPiZero = 0.052;
  rateSigmaMinusPiPlus = 0.199;
  rateSigmaPlusPiMinus = 0.446;
  rateSigmaZeroPiZero  = 0.303;
  rateLambdaZeroPiMinus = 0.568;
  rateSigmaZeroPiMinus  = 0.216;
  rateSigmaMinusPiZero  = 0.216;

  // for sigma- p -> lambda n
  //     sigma+ n -> lambda p
  //     sigma- n -> lambda
  // all values compatible with 0.55 same literature as above.
  sigmaPlusLambdaConversionRate = 0.55;
  sigmaMinusLambdaConversionRate = 0.55;
  sigmaZeroLambdaConversionRate = 0.55;

  //atomic cascade transitions
  //see yields from DEAR-SIDDHARTA bibliography
  //M.Iliescu, Feb. 2010
  yieldKa=1.;probKa=.3;
  yieldKb=1.;probKb=.2;
  yieldKg=1.;probKg=.15;
  yieldKd=1.;probKd=.1;
  yieldKe=1.;probKe=.1;
  yieldKz=1.;probKz=.1;
  yieldKi=1.;probKi=.05;

  yieldXaOther=1.;probXaOther=.75;
  G4HadronicProcessStore::Instance()->RegisterExtraProcess(this);
  
  GenerateKaonNucleusAbsLines();
}


G4KaonMinusAbsorptionAtRest::~G4KaonMinusAbsorptionAtRest()
{
  G4HadronicProcessStore::Instance()->DeRegisterExtraProcess(this);
}

void G4KaonMinusAbsorptionAtRest::PreparePhysicsTable(const G4ParticleDefinition& p)
{
  G4HadronicProcessStore::Instance()->RegisterParticleForExtraProcess(this, &p);
}

void G4KaonMinusAbsorptionAtRest::BuildPhysicsTable(const G4ParticleDefinition& p)
{
  G4HadronicProcessStore::Instance()->PrintInfo(&p);
}

G4VParticleChange* G4KaonMinusAbsorptionAtRest::AtRestDoIt(const G4Track& track, const G4Step& )
{
  stoppedHadron = track.GetDynamicParticle();

  if (!IsApplicable(*(stoppedHadron->GetDefinition()))) {
    G4cerr << "G4KaonMinusAbsorptionAtRest:ERROR, particle must be a Kaon!" << G4endl;
    return nullptr;
  }
  currentMaterial = track.GetMaterial();
  G4double currentTime = track.GetGlobalTime();

  nucleus = nullptr;
  do {
    nucleus = new G4Nucleus(currentMaterial);

// cascade modified, original value 1.5
    if (nucleus->GetA_asInt() < 0.5) {
      delete nucleus;
      nucleus = nullptr;
    }

  } while (nucleus == nullptr);

  G4int Z = nucleus->GetZ_asInt();
  G4int A = nucleus->GetA_asInt();
  G4DynamicParticleVector* absorptionProducts = KaonNucleonReaction();

// For A=3 N=0 fragmentation case in He target 12.2020 H. Shi
  G4int fragZ = nucleus->GetZ_asInt();
  G4int fragA = nucleus->GetA_asInt();
  if (A >= 1.5 && !(fragA == 3 && fragZ == 0) )     // not removing the A=3 Z=0 error case..
  {
//Secondary interactions
    G4DynamicParticle* thePion;
    for (int i=0; i<absorptionProducts->size(); i++) {
      thePion = (*absorptionProducts)[i];
      if (thePion->GetDefinition() == G4PionMinus::PionMinus()
            || thePion->GetDefinition() == G4PionPlus::PionPlus()
            || thePion->GetDefinition() == G4PionZero::PionZero())
      {
        if (AbsorbPionByNucleus(thePion)) {
          absorptionProducts->erase(absorptionProducts->begin()+i);
          i--;
          delete thePion;
          if (verboseLevel > 1)
            G4cout << "G4KaonMinusAbsorption::AtRestDoIt: Pion absorbed in Nucleus" << G4endl;
        }
      }
    }

    G4DynamicParticle* theSigma;
    G4DynamicParticle* theLambda;
    for (int i=0; i<absorptionProducts->size(); i++) {
      theSigma = (*absorptionProducts)[i];
      if (theSigma->GetDefinition() == G4SigmaMinus::SigmaMinus()
            || theSigma->GetDefinition() == G4SigmaPlus::SigmaPlus()
            || theSigma->GetDefinition() == G4SigmaZero::SigmaZero())
      {
        theLambda = SigmaLambdaConversion(theSigma);
        if (theLambda != 0) {
          absorptionProducts->erase(absorptionProducts->begin()+i);
          i--;
          delete theSigma;
          absorptionProducts->push_back(theLambda);

          if (verboseLevel > 1)
            G4cout << "G4KaonMinusAbsorption::AtRestDoIt: SigmaLambdaConversion Done" << G4endl;
        }
      }
    }

// Nucleus deexcitation
    G4double productEnergy = 0.;
    G4ThreeVector pProducts(0.,0.,0.);

    unsigned int nAbsorptionProducts = 0;
    if (absorptionProducts != nullptr)
      nAbsorptionProducts = absorptionProducts->size();

    for (int i=0; i<nAbsorptionProducts; i++) {
      pProducts += (*absorptionProducts)[i]->GetMomentum();
      productEnergy += (*absorptionProducts)[i]->GetKineticEnergy();
    }

    G4int newZ = nucleus->GetZ_asInt();
    G4int newA = nucleus->GetA_asInt();

    G4double bDiff = G4NucleiProperties::GetBindingEnergy(A,Z) -
                            G4NucleiProperties::GetBindingEnergy(newA, newZ);

    G4StopDeexcitationAlgorithm* nucleusAlgorithm = new G4StopTheoDeexcitation();
    G4StopDeexcitation stopDeexcitation(nucleusAlgorithm);

    nucleus->AddExcitationEnergy(bDiff);

    G4double energyDeposit = nucleus->GetEnergyDeposit();
    if (verboseLevel > 0) {
      G4cout << " -- KaonAtRest -- excitation = " << energyDeposit
             << ", pNucleus = " << pProducts
             << ", A: " << A
             << ", " << newA
             << ", Z: " << Z
             << ", " << newZ << G4endl;
    }

    if (energyDeposit < 0.) {
      G4Exception("G4KaonMinusAbsorptionAtRest", "007", FatalException, "AtRestDoIt -- excitation energy < 0");
    }
    delete nucleus;

    if (newA == 3 && newZ == 0) {
      G4cout << "---- Abnormal fragmentation products A = 3  Z = 0 -----" << G4endl;
      G4cout << "     Old A = " << A << "  old Z = " << Z << G4endl;
      G4cout << "     NewA = "  << newA << "  newZ = " << newZ << G4endl;
      G4cout << "     Skip stopDeexcitation break up     " << G4endl;
      G4cout << "     After KaonNucleonReaction(): A = "  << fragA << "  Z = " << fragZ << G4endl;
    } else {
      G4ReactionProductVector* fragmentationProducts = stopDeexcitation.DoBreakUp(newA, newZ, energyDeposit, pProducts);

      unsigned nFragmentationProducts = 0;
      if (fragmentationProducts != 0)
        nFragmentationProducts = fragmentationProducts->size();

//Initialize ParticleChange -> Internal variable in G4ParticleChange <- Note in Geant4 notes as possible to remove!!!
      aParticleChange.Initialize(track);
      aParticleChange.SetNumberOfSecondaries(G4int(nAbsorptionProducts + nFragmentationProducts));

//Update List of alive particles. put energy deposit at the right place ...
      for (int i=0; i<nAbsorptionProducts; i++) {
        aParticleChange.AddSecondary((*absorptionProducts)[i], currentTime);
      }
      if (absorptionProducts != nullptr)
        delete absorptionProducts;

      for(int i=0; i<nFragmentationProducts; i++) {
        G4DynamicParticle * aNew = new G4DynamicParticle((*fragmentationProducts)[i]->GetDefinition(),
                                                            (*fragmentationProducts)[i]->GetTotalEnergy(),
                                                            (*fragmentationProducts)[i]->GetMomentum());
        G4double newTime = aParticleChange.GetGlobalTime((*fragmentationProducts)[i]->GetFormationTime());
        aParticleChange.AddSecondary(aNew, newTime);
        delete (*fragmentationProducts)[i];
      }
      if (fragmentationProducts != 0)
        delete fragmentationProducts;
    }
  } else { //else works for -> A<1.5 || fragA==3 && fragZ==0
    unsigned nAbsorptionProducts = 0;

    delete nucleus;
    if (absorptionProducts != nullptr)
      nAbsorptionProducts = absorptionProducts->size();

//Initialize ParticleChange -> Internal variable in G4ParticleChange <- Note in Geant4 notes as possible to remove!!!
    aParticleChange.Initialize(track);
    aParticleChange.SetNumberOfSecondaries(G4int(nAbsorptionProducts) );

//Update List of alive particles.
    for (int i=0; i<nAbsorptionProducts; i++) {
      aParticleChange.AddSecondary((*absorptionProducts)[i],currentTime);
    }
    if (absorptionProducts != nullptr)
      delete absorptionProducts;
  }
  aParticleChange.ProposeTrackStatus(fStopAndKill); //fStopAndKill -> internal state (enum) in G4TrackStatus
  return &aParticleChange;
}


G4DynamicParticle G4KaonMinusAbsorptionAtRest::GetAbsorbingNucleon()
{
  G4DynamicParticle aNucleon;
//Atomic cascade modified
  if (nucleus->GetA_asInt() >= 1.5) {
//Get nucleon definition, based on Z,N of current Nucleus
  	aNucleon.SetDefinition(SelectAbsorbingNucleon());

//Fermi momentum distribution in three dimensions
  	G4ThreeVector pFermi = nucleus->GetFermiMomentum();
  	aNucleon.SetMomentum(pFermi);
  } else {
    aNucleon.SetDefinition(G4Proton::Proton());  //kaonic hydrogen case
    aNucleon.SetMomentum(G4ThreeVector(0.,0.,0.));
  }
  return aNucleon;
}

G4ParticleDefinition* G4KaonMinusAbsorptionAtRest::SelectAbsorbingNucleon()
{
// (Ch. Voelcker) extended from ReturnTargetParticle():
// Choose a proton or a neutron as the absorbing particle,
// taking weight into account!
// Update nucleon's atomic numbers.
  G4ParticleDefinition* absorbingParticleDef;
  G4double ranflat = G4UniformRand();

  G4int myZ = nucleus->GetZ_asInt();   // number of protons
  G4int myN = nucleus->GetA_asInt();   // number of nucleons (not neutrons!!)

// See  VanderVelde-Wilquet et al, Nuov.Cim.39A(1978)538;
  G4double carbonRatioNP = 0.18;  // (Rn/Rp)c, see page 544

  G4double neutronProtonRatio = NeutronHaloFactor(myZ,myN)*carbonRatioNP*(double)(myN-myZ)/(double)myZ;
  G4double protonProbability = 1./(1. + neutronProtonRatio);

  if (ranflat < protonProbability) {
    absorbingParticleDef = G4Proton::Proton();
    myZ -= 1.;
  } else {
    absorbingParticleDef = G4Neutron::Neutron();
  }
//It may call an exception. Is it inteded?
  myN -= 1.;
// remove the interacting nucleon from the current nucleus
  if (myN < 1) {
    myN = 1;
  }
  if (myZ < 0) {
    myZ = 0;
  }
  nucleus->SetParameters(myN,myZ);
  return absorbingParticleDef;
}


G4double G4KaonMinusAbsorptionAtRest::NeutronHaloFactor(G4double Z, G4double N)
{
// this function should take care of the probability for absorption
// on neutrons, depending on number of protons Z and number of neutrons N-Z
// parametrisation from fit to
// VanderVelde-Wilquet et al, Nuov.Cim.39A(1978)538;

  if (Z == 1.)
    return 1.389; // deuterium
  else if (Z == 2.)
    return 1.78; // helium
  else if (Z == 10.)
    return 0.66; // neon
  else
    return 0.6742+(N-Z)*0.06524;
}


G4DynamicParticleVector* G4KaonMinusAbsorptionAtRest::KaonNucleonReaction()
{
  G4DynamicParticleVector* products = new G4DynamicParticleVector();

  G4double ranflat = G4UniformRand();
  G4double prob = 0;

  G4ParticleDefinition* producedBaryonDef;
  G4ParticleDefinition* producedMesonDef;
  G4ParticleDefinition* producedBosonDef;  // atomic cascade photon

  G4int iniZ = nucleus->GetZ_asInt();
  G4int iniA = nucleus->GetA_asInt();

  G4DynamicParticle aNucleon = GetAbsorbingNucleon();
  G4double nucleonMass;

  if (aNucleon.GetDefinition() == G4Proton::Proton()) {
    nucleonMass = proton_mass_c2 + electron_mass_c2;
    if ((prob += rateLambdaZeroPiZero) > ranflat) {
      producedBaryonDef = G4Lambda::Lambda();
      producedMesonDef  = G4PionZero::PionZero();
    } else if ((prob += rateSigmaPlusPiMinus) > ranflat) {
      producedBaryonDef = G4SigmaPlus::SigmaPlus();
      producedMesonDef  = G4PionMinus::PionMinus();
    } else if ((prob += rateSigmaMinusPiPlus) > ranflat) {
      producedBaryonDef = G4SigmaMinus::SigmaMinus();
      producedMesonDef  = G4PionPlus::PionPlus();
    } else {
      producedBaryonDef = G4SigmaZero::SigmaZero();
      producedMesonDef  = G4PionZero::PionZero();
    }
  } else if (aNucleon.GetDefinition() == G4Neutron::Neutron()) {
    nucleonMass = neutron_mass_c2;
    if ((prob += rateLambdaZeroPiMinus) > ranflat) {
      producedBaryonDef = G4Lambda::Lambda();
      producedMesonDef  = G4PionMinus::PionMinus();
    } else if ((prob += rateSigmaZeroPiMinus) > ranflat) {
      producedBaryonDef = G4SigmaZero::SigmaZero();
      producedMesonDef = G4PionMinus::PionMinus();
    } else {
      producedBaryonDef = G4SigmaMinus::SigmaMinus();
      producedMesonDef  = G4PionZero::PionZero();
    }
  } else {
    if (verboseLevel > 0) {
      G4cout << "G4KaonMinusAbsorption::KaonNucleonReaction: " << aNucleon.GetDefinition()->GetParticleName()
             << " is not a good nucleon - check G4Nucleus::ReturnTargetParticle()!" << G4endl;
    }
    return nullptr;
  }
  G4DynamicParticle modifiedHadron = (*stoppedHadron);

  if (iniA >= 1.5) {            // kaonic hydrogen cascade
    G4int newZ = nucleus->GetZ_asInt();
    G4int newA = nucleus->GetA_asInt();

// Modify the Kaon mass to take nuclear binding energy into account
// .. using mass formula ..
// .. using mass table ..
// equivalent to '-initialBindingEnergy+nucleus.GetBindingEnergy' !

    G4double nucleonBindingEnergy = G4NucleiProperties::GetBindingEnergy(newA, newZ) - G4NucleiProperties::GetBindingEnergy(iniA, iniZ);
    modifiedHadron.SetMass(stoppedHadron->GetMass() + nucleonBindingEnergy);
  }
//Is it needed still?
//else {							//implicit
//	modifiedHadron.SetMass(stoppedHadron->GetMass());
//}

// Setup outgoing dynamic particles
  G4ThreeVector dummy(0.,0.,0.);
  G4DynamicParticle* producedBaryon = new G4DynamicParticle(producedBaryonDef, dummy);
  G4DynamicParticle* producedMeson = new G4DynamicParticle(producedMesonDef, dummy);

// Produce the secondary particles in a twobody process:
  G4ReactionKinematics theReactionKinematics;
  theReactionKinematics.TwoBodyScattering(&modifiedHadron, &aNucleon, producedBaryon, producedMeson);

//====kaonic atom cascade photons=================================================

  producedBosonDef = G4Gamma::Gamma();
  G4double photonEnergy;
  std::vector<G4double> photonEnergy1;
    
  std::vector<G4double> photonEnergies;
  if (iniA == 2 && iniZ == 1) {
    photonEnergies = nucleiImplemented.at(atomicNumberVSnucleusRef.at(-iniZ)).GetPhotonEnergies(currentMaterial->GetNumberOfElements());
  } else {
    photonEnergies = nucleiImplemented.at(atomicNumberVSnucleusRef.at(iniZ)).GetPhotonEnergies(currentMaterial->GetNumberOfElements());
  }
    
  for (unsigned i=0; i<photonEnergies.size(); i++) {
    G4DynamicParticle* producedBoson = new G4DynamicParticle(producedBosonDef, dummy);

    G4double costheta = 2.0*G4UniformRand() - 1.0;
    G4double sintheta = std::sqrt(1.0 - costheta*costheta);
    G4double phi = 2.0*pi*G4UniformRand();
    G4double photonMomentum = photonEnergies.at(i);

    G4double pz=costheta*photonMomentum;
    G4double px=sintheta*std::cos(phi)*photonMomentum;
    G4double py=sintheta*std::sin(phi)*photonMomentum;

    G4ThreeVector photMomentum(px,py,pz);
    producedBoson->SetMomentum(photMomentum);
    products->push_back(producedBoson);
  }

//==========================end cascade==========================================
  if (fPlotProducedXrayLies) {
    SiddhartaAnalysisManager* analysis = SiddhartaAnalysisManager::getInstance();
    for (unsigned it=0; it<products->size(); it++) {
      analysis->histo->fillHistogram("kaonAbsorptionEnergies_vs_Nucleus", products->at(it)->GetTotalMomentum()/eV, doubleCheck(iniA), doubleCheck(iniA));
    }
  }

  products->push_back(producedBaryon);
  products->push_back(producedMeson);

  if (verboseLevel > 1) {
    G4cout << "G4KaonMinusAbsorption::KaonNucleonReaction: Number of primaries = " << products->size()
           << ": " << producedMesonDef->GetParticleName() << ", " << producedBaryonDef->GetParticleName()
           << ", " << producedBosonDef->GetParticleName() << G4endl;
  }
  return products;
}


G4bool G4KaonMinusAbsorptionAtRest::AbsorbPionByNucleus(G4DynamicParticle* aPion)
{
// Needs some more investigation!
  G4double ranflat = G4UniformRand();

  if (ranflat < pionAbsorptionRate) {
// Add pion energy to ExcitationEnergy and NucleusMomentum
    nucleus->AddExcitationEnergy(aPion->GetTotalEnergy());
    nucleus->AddMomentum(aPion->GetMomentum());
  }

  return (ranflat < pionAbsorptionRate);
}


G4DynamicParticle* G4KaonMinusAbsorptionAtRest::SigmaLambdaConversion(G4DynamicParticle* aSigma)
{
  G4double  ranflat = G4UniformRand();
  G4double  sigmaLambdaConversionRate;

  G4int A = nucleus->GetA_asInt();
  G4int Z = nucleus->GetZ_asInt();

  G4int newZ = Z;
  G4double nucleonMassDifference = 0;

  G4ParticleDefinition* inNucleonDef = nullptr;
  G4ParticleDefinition* outNucleonDef = nullptr;

  // Decide which sigma
  switch((int)aSigma->GetDefinition()->GetPDGCharge()) {
    case 1:
      sigmaLambdaConversionRate = sigmaPlusLambdaConversionRate;
      inNucleonDef   = G4Neutron::Neutron();
      outNucleonDef  = G4Proton::Proton();
      newZ = Z+1;
      nucleonMassDifference = neutron_mass_c2 - (proton_mass_c2 + electron_mass_c2);
      break;
    case -1:
      sigmaLambdaConversionRate = sigmaMinusLambdaConversionRate;
      inNucleonDef   = G4Proton::Proton();
      outNucleonDef  = G4Neutron::Neutron();
      newZ = Z-1;
      nucleonMassDifference = (proton_mass_c2 + electron_mass_c2) - neutron_mass_c2;
      break;
    case 0:
      sigmaLambdaConversionRate = sigmaZeroLambdaConversionRate;
// The 'outgoing' nucleon is just virtual, to keep the energy-momentum
// balance and will not appear in the ParticleChange. Therefore no need
// choose between neutron and proton here!
      inNucleonDef   = G4Neutron::Neutron();
      outNucleonDef  = G4Neutron::Neutron();
      break;
    default:
      sigmaLambdaConversionRate = 0.;
  }
  if (ranflat >= sigmaLambdaConversionRate)
    return nullptr;

  G4ThreeVector dummy(0.,0.,0.);

// Fermi momentum distribution in three dimensions
  G4ThreeVector momentum = nucleus->GetFermiMomentum();
  G4ParticleDefinition* lambdaDef  = G4Lambda::Lambda();

  G4DynamicParticle inNucleon(inNucleonDef, momentum);
  G4DynamicParticle outNucleon(outNucleonDef, dummy);
  G4DynamicParticle* outLambda = new G4DynamicParticle(lambdaDef, dummy);
  G4ReactionKinematics theReactionKinematics;
// Now do the twobody scattering
  theReactionKinematics.TwoBodyScattering(aSigma, &inNucleon, &outNucleon, outLambda);

// Binding energy of nucleus has changed. This will change the
// ExcitationEnergy.
// .. using mass formula ..
// .. using mass table ..
// equivalent to -'initialBindingEnergy+nucleus.GetBindingEnergy' !

// Add energy and momentum to nucleus, change Z,A
  nucleus->AddExcitationEnergy(outNucleon.GetKineticEnergy());
  nucleus->AddMomentum(outNucleon.GetMomentum());

  if (newZ<0) //Problem for deuterium -> what is a neutron electron pair?!
    newZ = 0;
  nucleus->SetParameters(A,newZ);
// The calling routine is responsible to delete the sigma!!
  return outLambda;
}

void G4KaonMinusAbsorptionAtRest::GenerateKaonNucleusAbsLines()
{
  std::map<G4int, unsigned> atomicNumberVSnucleusRef;
  std::vector<KaonNucleusAbsLines> nucleiImplemented;
  
//Hydrogen
  KaonNucleusAbsLines hydrogen;
  hydrogen.SetAtomicNumber(1);
  hydrogen.AddNucleusNumber(1);
  hydrogen.AddPhotonEnergy(6191.8*eV, 1.0);
  hydrogen.AddPhotonEnergy(7388.6*eV, 1.0);
  hydrogen.SetSimType(AbsLineSimType::absSimulateAll);
  hydrogen.SetCleanMaterialOnly();
  atomicNumberVSnucleusRef.insert(std::pair{1,nucleiImplemented.size()});
  nucleiImplemented.push_back(hydrogen);
//Deuterium -> giving minus atomic number to differentiate between H and D
  KaonNucleusAbsLines deuterium;
  deuterium.SetAtomicNumber(1);
  deuterium.AddNucleusNumber(2);
  G4double shiftForD = -800.*eV;
  deuterium.AddPhotonEnergy(7834*eV + shiftForD, 1.0);
  deuterium.AddPhotonEnergy(9280.2*eV + shiftForD, 1.0);
  deuterium.SetSimType(AbsLineSimType::absSimulateAll);
  deuterium.SetCleanMaterialOnly();
  atomicNumberVSnucleusRef.insert(std::pair{-1,nucleiImplemented.size()});
  nucleiImplemented.push_back(deuterium);
//Helium
  KaonNucleusAbsLines helium;
  helium.SetAtomicNumber(2);
  helium.AddNucleusNumber(4);
  G4double shiftForHe = 0.*eV;
  helium.AddPhotonEnergy(6463.6*eV + shiftForHe, 1.0);
  helium.SetSimType(AbsLineSimType::absSimulateAll);
  helium.SetCleanMaterialOnly();
  atomicNumberVSnucleusRef.insert(std::pair{2,nucleiImplemented.size()});
  nucleiImplemented.push_back(helium);
//Carbon
  KaonNucleusAbsLines carbon;
  carbon.SetAtomicNumber(6);
  carbon.AddNucleusNumber(12);
  carbon.AddPhotonEnergy(62881.1*eV, 1.0); // 3-->2
  carbon.AddPhotonEnergy(22008.4*eV, 1.0); // 4-->3
  carbon.AddPhotonEnergy(32195.1*eV, 1.0); // 5-->3
  carbon.AddPhotonEnergy(10216.5*eV, 1.0); // 5-->4
  carbon.AddPhotonEnergy(15809.0*eV, 1.0); // 6-->4
  carbon.AddPhotonEnergy(5544.9*eV, 1.0); // 6-->5
  carbon.AddPhotonEnergy(8885.8*eV, 1.0); // 7-->5
  carbon.SetSimType(AbsLineSimType::absSimulateAll);
  atomicNumberVSnucleusRef.insert(std::pair{6,nucleiImplemented.size()});
  nucleiImplemented.push_back(carbon);
//Nitrogen
  KaonNucleusAbsLines nitrogen;
  nitrogen.SetAtomicNumber(7);
  nitrogen.AddNucleusNumber(14);
  nitrogen.AddPhotonEnergy(13995.4*eV, 1.0);
  nitrogen.AddPhotonEnergy(21588.2*eV, 1.0);
  nitrogen.AddPhotonEnergy(7595.4*eV, 1.0);
  nitrogen.AddPhotonEnergy(12223.0*eV, 1.0);
  nitrogen.AddPhotonEnergy(4577.1*eV, 1.0);
  nitrogen.AddPhotonEnergy(7578.0*eV, 1.0);
  nitrogen.SetSimType(AbsLineSimType::absSimulateAll);
  atomicNumberVSnucleusRef.insert(std::pair{7,nucleiImplemented.size()});
  nucleiImplemented.push_back(nitrogen);
//Oxygen
  KaonNucleusAbsLines oxygen;
  oxygen.SetAtomicNumber(8);
  oxygen.AddNucleusNumber(16);
  oxygen.AddPhotonEnergy(9968.7*eV, 1.0);
  oxygen.AddPhotonEnergy(16062.0*eV, 1.0);
  oxygen.AddPhotonEnergy(6006.8*eV, 1.0);
  oxygen.AddPhotonEnergy(9958.0*eV, 1.0);
  oxygen.SetSimType(AbsLineSimType::absSimulateAll);
  atomicNumberVSnucleusRef.insert(std::pair{8,nucleiImplemented.size()});
  nucleiImplemented.push_back(oxygen);
//Fluorine
  KaonNucleusAbsLines fluorine;
  fluorine.SetAtomicNumber(9);
  fluorine.AddNucleusNumber(19);
  fluorine.AddPhotonEnergy(7642.7*eV, 1.0);
  fluorine.AddPhotonEnergy(12599.0*eV, 1.0);
  fluorine.SetSimType(AbsLineSimType::absSimulateAll);
  atomicNumberVSnucleusRef.insert(std::pair{9,nucleiImplemented.size()});
  nucleiImplemented.push_back(fluorine);
//Neon
  KaonNucleusAbsLines neon;
  neon.SetAtomicNumber(10);
  neon.AddNucleusNumber(20);
  neon.AddNucleusNumber(22);
  neon.AddPhotonEnergy(6118.91*eV, 1.0);
  neon.AddPhotonEnergy(9427.65*eV, 1.0);
  neon.AddPhotonEnergy(15635.4*eV, 1.0);
  neon.SetSimType(AbsLineSimType::absSimulateAll);
  neon.SetCleanMaterialOnly();
  atomicNumberVSnucleusRef.insert(std::pair{10,nucleiImplemented.size()});
  nucleiImplemented.push_back(neon);
//Aluminium
  KaonNucleusAbsLines aluminium;
  aluminium.SetAtomicNumber(13);
  aluminium.AddNucleusNumber(27);
  aluminium.AddPhotonEnergy(302293.0*eV, 1.0); // 3-->2
  aluminium.AddPhotonEnergy(105803.0*eV, 1.0); // 4-->3
  aluminium.AddPhotonEnergy(48972.0*eV, 1.0); // 5-->4
  aluminium.AddPhotonEnergy(75573.0*eV, 1.0); // 6-->4
  aluminium.AddPhotonEnergy(154774.0*eV, 1.0); // 5-->3
  aluminium.AddPhotonEnergy(26602.0*eV, 1.0); // 6-->5
  aluminium.AddPhotonEnergy(42642.0*eV, 1.0); // 7-->5
  aluminium.SetSimType(AbsLineSimType::absSimulateAll);
  atomicNumberVSnucleusRef.insert(std::pair{13,nucleiImplemented.size()});
  nucleiImplemented.push_back(aluminium);
//Silicon
  KaonNucleusAbsLines silicon;
  silicon.SetAtomicNumber(14);
  silicon.AddNucleusNumber(28);
  silicon.AddPhotonEnergy(8300.1*eV, 1.0);
  silicon.AddPhotonEnergy(14233.1*eV, 1.0);
  silicon.AddPhotonEnergy(5935.0*eV, 1.0);
  silicon.AddPhotonEnergy(10324.0*eV, 1.0);
  silicon.AddPhotonEnergy(4390.1*eV, 1.0);
  silicon.AddPhotonEnergy(7728.0*eV, 1.0);
  silicon.SetSimType(AbsLineSimType::absSimulateAll);
  atomicNumberVSnucleusRef.insert(std::pair{14,nucleiImplemented.size()});
  nucleiImplemented.push_back(silicon);
//Titanium
  KaonNucleusAbsLines titanium;
  titanium.SetAtomicNumber(22);
  titanium.AddNucleusNumber(48);
  titanium.AddPhotonEnergy(8312.0*eV, 1.0);
  titanium.AddPhotonEnergy(14770.0*eV, 1.0);
  titanium.AddPhotonEnergy(6467.0*eV, 1.0);
  titanium.AddPhotonEnergy(11590.0*eV, 1.0);
  titanium.SetSimType(AbsLineSimType::absSimulateAll);
  atomicNumberVSnucleusRef.insert(std::pair{22,nucleiImplemented.size()});
  nucleiImplemented.push_back(titanium);
//Molibden
  KaonNucleusAbsLines molibden;
  molibden.SetAtomicNumber(42);
  molibden.AddNucleusNumber(96);
  molibden.AddNucleusNumber(98);
  molibden.AddPhotonEnergy(284000.0*eV, 1.0); // 6-->5
  molibden.AddPhotonEnergy(514000.0*eV, 1.0); // 5-->4
 // molibden.AddPhotonEnergy(0.0*eV, 1.0); // 4-->3
 // molibden.AddPhotonEnergy(0.0*eV, 1.0); // 7-->6
  molibden.SetSimType(AbsLineSimType::absSimulateAll);
  atomicNumberVSnucleusRef.insert(std::pair{22,nucleiImplemented.size()});
  nucleiImplemented.push_back(molibden);
//Lead
  KaonNucleusAbsLines lead;
  lead.SetAtomicNumber(82);
  lead.AddNucleusNumber(206);
  lead.AddNucleusNumber(207);
  lead.AddNucleusNumber(208);
  lead.AddPhotonEnergy(288822.0*eV, 1.0); // 3-->2
  lead.AddPhotonEnergy(206593.0*eV, 1.0); // 4-->3
  lead.AddPhotonEnergy(152855.0*eV, 1.0); // 5-->4
  lead.AddPhotonEnergy(269114.0*eV, 1.0); // 6-->4
  lead.AddPhotonEnergy(359448.0*eV, 1.0); // 5-->3
  lead.AddPhotonEnergy(116259.0*eV, 1.0); // 6-->5
  lead.AddPhotonEnergy(206736.0*eV, 1.0); // 7-->5
  lead.SetSimType(AbsLineSimType::absSimulateAll);
  atomicNumberVSnucleusRef.insert(std::pair{82,nucleiImplemented.size()});
  nucleiImplemented.push_back(lead);     
}
