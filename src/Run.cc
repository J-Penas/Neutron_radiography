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
//
/// \file Run.cc
/// \brief Implementation of the Run class
//
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "Run.hh"
#include "DetectorConstruction.hh"
#include "PrimaryGeneratorAction.hh"
#include "HistoManager.hh"

#include "G4SDManager.hh"
#include "G4Event.hh"
#include "G4THitsMap.hh"

#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Run::Run(DetectorConstruction* det)
: G4Run(),
  fDetector(det), fParticle(0), fEkin(0.),
  fNbStep1(0), fNbStep2(0),
  fTrackLen1(0.), fTrackLen2(0.),
  fTime1(0.),fTime2(0.),
  reID(0),
  nEvent(0),idx(0),
  Edata(),Ndata(),Hitsdata()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Run::~Run()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::SetPrimary(G4ParticleDefinition* particle, G4double energy)
{ 
  fParticle = particle;
  fEkin = energy;
} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::CountProcesses(const G4VProcess* process) 
{
  if (process == nullptr) return;
  G4String procName = process->GetProcessName();
  std::map<G4String,G4int>::iterator it = fProcCounter.find(procName);
  if ( it == fProcCounter.end()) {
    fProcCounter[procName] = 1;
  }
  else {
    fProcCounter[procName]++; 
  }
}                 
                  
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::ParticleCount(G4String name, G4double Ekin)
{
  std::map<G4String, ParticleData>::iterator it = fParticleDataMap.find(name);
  if ( it == fParticleDataMap.end()) {
    fParticleDataMap[name] = ParticleData(1, Ekin, Ekin, Ekin);
  }
  else {
    ParticleData& data = it->second;
    data.fCount++;
    data.fEmean += Ekin;
    //update min max
    G4double emin = data.fEmin;
    if (Ekin < emin) data.fEmin = Ekin;
    G4double emax = data.fEmax;
    if (Ekin > emax) data.fEmax = Ekin; 
  }   
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::SumTrackLength(G4int nstep1, G4int nstep2, 
                         G4double trackl1, G4double trackl2,
                         G4double time1, G4double time2)
{
  fNbStep1   += nstep1;  fNbStep2   += nstep2;
  fTrackLen1 += trackl1; fTrackLen2 += trackl2;
  fTime1 += time1; fTime2 += time2;  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::RecordEvent(const G4Event* evt)
{
  
  idx = 0;

  for (G4int i=-8; i<1; i++) {
    for (G4int j=1; j<10; j++) {

    //generate name
    G4String Sname = "score_r/re" + std::to_string(i+8) + "_" + std::to_string(j);

    //indentify and get hitCollection
    reID = G4SDManager::GetSDMpointer()->GetCollectionID(Sname);

    G4HCofThisEvent* HCE = evt->GetHCofThisEvent();
    G4THitsMap<G4double>* evt_re = (G4THitsMap<G4double>*)(HCE->GetHC(reID));

    //fill the score with the number of events
    Hitsdata[idx] += *evt_re;

    //G4cout << "Este evento se ha registrado" << idx << G4endl;
    idx++;

    }
  }

  //identify each hitcollection
  //re0ID = G4SDManager::GetSDMpointer()->GetCollectionID("score_r/re0");
  //re1ID = G4SDManager::GetSDMpointer()->GetCollectionID("score_r/re1");
  //re2ID = G4SDManager::GetSDMpointer()->GetCollectionID("score_r/re2");
  //re3ID = G4SDManager::GetSDMpointer()->GetCollectionID("score_r/re3");
  //re4ID = G4SDManager::GetSDMpointer()->GetCollectionID("score_r/re4");
  //re5ID = G4SDManager::GetSDMpointer()->GetCollectionID("score_r/re5");
  //reTID = G4SDManager::GetSDMpointer()->GetCollectionID("score_r/reT");

  //get hitcollection of each score
  //G4HCofThisEvent* HCE = evt->GetHCofThisEvent();
  //G4THitsMap<G4double>* evt_re0 = (G4THitsMap<G4double>*)(HCE->GetHC(re0ID));
  //G4THitsMap<G4double>* evt_re1 = (G4THitsMap<G4double>*)(HCE->GetHC(re1ID));
  //G4THitsMap<G4double>* evt_re2 = (G4THitsMap<G4double>*)(HCE->GetHC(re2ID));
  //G4THitsMap<G4double>* evt_re3 = (G4THitsMap<G4double>*)(HCE->GetHC(re3ID));
  //G4THitsMap<G4double>* evt_re4 = (G4THitsMap<G4double>*)(HCE->GetHC(re4ID));
  //G4THitsMap<G4double>* evt_re5 = (G4THitsMap<G4double>*)(HCE->GetHC(re5ID));
  //G4THitsMap<G4double>* evt_reT = (G4THitsMap<G4double>*)(HCE->GetHC(reTID));

  //fill the scores with the number of events
  //re0 += *evt_re0;
  //re1 += *evt_re1;
  //re2 += *evt_re2;
  //re3 += *evt_re3;
  //re4 += *evt_re4;
  //re5 += *evt_re5;
  //reT += *evt_reT;

  G4Run::RecordEvent(evt);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::Merge(const G4Run* run)
{
  const Run* localRun = static_cast<const Run*>(run);
  
  //primary particle info
  //
  fParticle = localRun->fParticle;
  fEkin     = localRun->fEkin;
  
  // accumulate sums
  //
  fNbStep1   += localRun->fNbStep1;
  fNbStep2   += localRun->fNbStep2;   
  fTrackLen1 += localRun->fTrackLen1;  
  fTrackLen2 += localRun->fTrackLen2;
  fTime1     += localRun->fTime1;  
  fTime2     += localRun->fTime2;
  
  //map: processes count
  std::map<G4String,G4int>::const_iterator itp;
  for ( itp = localRun->fProcCounter.begin();
        itp != localRun->fProcCounter.end(); ++itp ) {

    G4String procName = itp->first;
    G4int localCount = itp->second;
    if ( fProcCounter.find(procName) == fProcCounter.end()) {
      fProcCounter[procName] = localCount;
    }
    else {
      fProcCounter[procName] += localCount;
    }  
  }
   
  //map: created particles count         
  std::map<G4String,ParticleData>::const_iterator itn;
  for (itn = localRun->fParticleDataMap.begin(); 
       itn != localRun->fParticleDataMap.end(); ++itn) {
    
    G4String name = itn->first;
    const ParticleData& localData = itn->second;   
    if ( fParticleDataMap.find(name) == fParticleDataMap.end()) {
      fParticleDataMap[name]
       = ParticleData(localData.fCount, 
                      localData.fEmean, 
                      localData.fEmin, 
                      localData.fEmax);
    }
    else {
      ParticleData& data = fParticleDataMap[name];   
      data.fCount += localData.fCount;
      data.fEmean += localData.fEmean;
      G4double emin = localData.fEmin;
      if (emin < data.fEmin) data.fEmin = emin;
      G4double emax = localData.fEmax;
      if (emax > data.fEmax) data.fEmax = emax; 
    }   
  }

  idx = 0;
  
  for (G4int i=-8; i<1; i++) {
    for (G4int j=1; j<10; j++) {
      
      //accumulate score counts
      Hitsdata[idx] += (localRun->Hitsdata[idx]);

      //G4cout << "Esta run se ha fusionado" << idx << G4endl;
      idx++;
    }
  }
  
  //accumulate score counts
  //re0 += (localRun->re0);
  //re1 += (localRun->re1);
  //re2 += (localRun->re2);
  //re3 += (localRun->re3);
  //re4 += (localRun->re4);
  //re5 += (localRun->re5);
  //reT += (localRun->reT);

  G4Run::Merge(run); 
} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::EndOfRun() 
{
  G4int prec = 5, wid = prec + 2;  
  G4int dfprec = G4cout.precision(prec);
  
  //run condition
  //
  G4Material* material = fDetector->GetMaterial();
  G4double density = material->GetDensity();
   
  G4String Particle = fParticle->GetParticleName();    
  G4cout << "\n The run is " << numberOfEvent << " "<< Particle << " of "
         << G4BestUnit(fEkin,"Energy") << " through " 
         << G4BestUnit((fDetector->GetSize()),"Length") << " of "
         << material->GetName() << " (density: " 
         << G4BestUnit(density,"Volumic Mass") << ")" << G4endl;

  if (numberOfEvent == 0) { G4cout.precision(dfprec);   return;}
             
  //frequency of processes
  //
  G4cout << "\n Process calls frequency :" << G4endl;  
  G4int survive = 0;
  std::map<G4String,G4int>::iterator it;    
  for (it = fProcCounter.begin(); it != fProcCounter.end(); it++) {
     G4String procName = it->first;
     G4int    count    = it->second;
     G4cout << "\t" << procName << "= " << count;
     if (procName == "Transportation") survive = count;
  }
  G4cout << G4endl;
      
  if (survive > 0) {
    G4cout << "\n Nb of incident particles surviving after "
           << G4BestUnit(0.5*(fDetector->GetSize()),"Length") << " of "
           << fDetector->GetMaterial()->GetName() << " : " << survive << G4endl;
  }

 // total track length of incident neutron
 //
 G4cout << "\n Parcours of incident neutron:";
  
 G4double meanCollision1  = (G4double)fNbStep1/numberOfEvent;
 G4double meanCollision2  = (G4double)fNbStep2/numberOfEvent;
 G4double meanCollisTota  = meanCollision1 + meanCollision2;

 G4cout << "\n   nb of collisions    E>1*eV= " << meanCollision1
        << "      E<1*eV= " << meanCollision2
        << "       total= " << meanCollisTota;        
        
 G4double meanTrackLen1  = fTrackLen1/numberOfEvent;
 G4double meanTrackLen2  = fTrackLen2/numberOfEvent;
 G4double meanTrackLtot  =  meanTrackLen1 + meanTrackLen2;  

 G4cout 
   << "\n   track length        E>1*eV= " << G4BestUnit(meanTrackLen1,"Length")
   << "  E<1*eV= " << G4BestUnit(meanTrackLen2, "Length")
   << "   total= " << G4BestUnit(meanTrackLtot, "Length");   
   
 G4double meanTime1  = fTime1/numberOfEvent;
 G4double meanTime2  = fTime2/numberOfEvent;
 G4double meanTimeTo = meanTime1 + meanTime2;  

 G4cout 
   << "\n   time of flight      E>1*eV= " << G4BestUnit(meanTime1,"Time")
   << "  E<1*eV= " << G4BestUnit(meanTime2, "Time")
   << "   total= " << G4BestUnit(meanTimeTo, "Time") << G4endl;   
             
 //particles count
 //
 G4cout << "\n List of generated particles:" << G4endl;
     
 std::map<G4String,ParticleData>::iterator itn;               
 for (itn = fParticleDataMap.begin(); itn != fParticleDataMap.end(); itn++) { 
    G4String name = itn->first;
    ParticleData data = itn->second;
    G4int count = data.fCount;
    G4double eMean = data.fEmean/count;
    G4double eMin = data.fEmin;
    G4double eMax = data.fEmax;    
         
    G4cout << "  " << std::setw(13) << name << ": " << std::setw(7) << count
           << "  Emean = " << std::setw(wid) << G4BestUnit(eMean, "Energy")
           << "\t( "  << G4BestUnit(eMin, "Energy")
           << " --> " << G4BestUnit(eMax, "Energy") 
           << ")" << G4endl;           
 }

 //print score results

 idx = 0;
 //std::map<G4String,G4int>::iterator it;

  for (G4int i=-8; i<1; i++) {
    for (G4int j=1; j<10; j++) {

      for (auto it : *Hitsdata[idx].GetMap()) {
        Ndata[idx] += *(it.second);;
      }

      Edata.push_back(j*pow(10,i)*MeV);

      idx++;
    }
  }

  G4double Size = fDetector->GetSize()/10;
  G4String filename = "Sp_" + material->GetName() + "_" + std::to_string(Size) + ".csv";
  
  static std::ofstream results(filename);
  results << material->GetName() << std::endl;
  results << fDetector->GetSize()/10 << std::endl;
  results << numberOfEvent << std::endl;
  results << "#,E (MeV),N (part.)" << std::endl;

  idx = 0;
  for (G4int i=-8; i<1; i++) {
    for (G4int j=1; j<10; j++) {
      results << "," << Edata[idx] << "," << Ndata[idx] << std::endl;
      idx++;
    }
  }

  //fill histogram
  //G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  //for (G4int i=-8; i<1; i++) {
  //  for (G4int j=1; j<10; j++) {
  //    analysisManager->FillH1(8,Ndata[idx]);
  //  }
  //}

  //G4cout << Edata << G4endl;
  //G4cout << Ndata << G4endl;

 //auto m0 = re0.GetMap();
 //G4int N0 = *(m0.second);
 //G4double N0;
 //G4double N1;
 //G4double N2;
 //G4double N3;
 //G4double N4;
 //G4double N5;
 //G4double NT;

  //for ( auto it0 : *re0.GetMap() ) {
  //  N0 += *(it0.second);
  //}
  //for ( auto it1 : *re1.GetMap() ) {
  //  N1 += *(it1.second);
  //}
  //for ( auto it2 : *re2.GetMap() ) {
  //  N2 += *(it2.second);
  //}
  //for ( auto it3 : *re3.GetMap() ) {
  //  N3 += *(it3.second);
  //}
  //for ( auto it4 : *re4.GetMap() ) {
  //  N4 += *(it4.second);
  //}
  //for ( auto it5 : *re5.GetMap() ) {
  //  N5 += *(it5.second);
  //}
  //for ( auto itT : *reT.GetMap() ) {
  //  NT += *(itT.second);
  //}
  // hitsMap->GetMap() returns the map of std::map<G4int, G4double*>


 //G4cout << "\n Number of events: " << numberOfEvent << G4endl;
 //G4cout << " Rear side ------------------" << G4endl;
 //G4cout << "\n E:\t    0 eV - 1 eV     : " << N0 << G4endl;
 //G4cout << " E:\t    1 eV - 500 eV   : " << N1 << G4endl;
 //G4cout << " E:\t 0.5 keV - 1 keV    : " << N2 << G4endl;
 //G4cout << " E:\t   1 keV - 500 keV  : " << N3 << G4endl;
 //G4cout << " E:\t 0.5 MeV - 1 MeV    : " << N4 << G4endl;
 //G4cout << " E:\t   1 MeV - 2 MeV    : " << N5 << G4endl;
 //G4cout << " Total number: " << NT << G4endl;
 //G4cout << "\n ----------------------------" << G4endl;
 
  //normalize histograms      
  ////G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  ////G4double factor = 1./numberOfEvent;
  ////analysisManager->ScaleH1(3,factor);
           
  //remove all contents in fProcCounter, fCount 
  fProcCounter.clear();
  fParticleDataMap.clear();
                          
  //restore default format         
  G4cout.precision(dfprec);   
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
