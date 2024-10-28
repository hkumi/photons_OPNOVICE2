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
/// \file optical/OpNovice2/include/DetectorConstruction.hh
/// \brief Definition of the DetectorConstruction class
//
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1

#include "globals.hh"
#include "G4OpticalSurface.hh"
#include "G4RunManager.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4SubtractionSolid.hh"

class DetectorMessenger;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class DetectorConstruction : public G4VUserDetectorConstruction
{
 public:
  DetectorConstruction();
  virtual ~DetectorConstruction();

  G4VPhysicalVolume* GetTank() { return fTank; }
  G4double GetTankXSize() { return fTank_x; }

  G4OpticalSurface* GetSurface(void) { return fSurface; }

  void SetSurfaceFinish(const G4OpticalSurfaceFinish finish)
  {
    fSurface->SetFinish(finish);
    G4RunManager::GetRunManager()->GeometryHasBeenModified();
  }
  G4OpticalSurfaceFinish GetSurfaceFinish(void)
  {
    return fSurface->GetFinish();
  }

  void SetSurfaceType(const G4SurfaceType type)
  {
    fSurface->SetType(type);
    G4RunManager::GetRunManager()->GeometryHasBeenModified();
  }

  void SetSurfaceModel(const G4OpticalSurfaceModel model)
  {
    fSurface->SetModel(model);
    G4RunManager::GetRunManager()->GeometryHasBeenModified();
  }
  G4OpticalSurfaceModel GetSurfaceModel(void) { return fSurface->GetModel(); }

  void SetSurfaceSigmaAlpha(G4double v);
  void SetSurfacePolish(G4double v);

  void AddTankMPV(const G4String& prop, G4MaterialPropertyVector* mpv);
  void AddTankMPC(const G4String& prop, G4double v);
  G4MaterialPropertiesTable* GetTankMaterialPropertiesTable()
  {
    return fTankMPT;
  }

  void AddWorldMPV(const G4String& prop, G4MaterialPropertyVector* mpv);
  void AddWorldMPC(const G4String& prop, G4double v);
  G4MaterialPropertiesTable* GetWorldMaterialPropertiesTable()
  {
    return fWorldMPT;
  }

  void AddSurfaceMPV(const G4String& prop, G4MaterialPropertyVector* mpv);
  void AddSurfaceMPC(const G4String& prop, G4double v);
  G4MaterialPropertiesTable* GetSurfaceMaterialPropertiesTable()
  {
    return fSurfaceMPT;
  }

  void SetWorldMaterial(const G4String&);
  G4Material* GetWorldMaterial() const { return fWorldMaterial; }
  void SetTankMaterial(const G4String&);
  G4Material* GetTankMaterial() const { return fTankMaterial; }

  virtual G4VPhysicalVolume* Construct();

 private:
  G4double fExpHall_x;
  G4double fExpHall_y;
  G4double fExpHall_z;

  G4VPhysicalVolume* fTank;

  G4double fTank_x;
  G4double fTank_y;
  G4double fTank_z;

  G4LogicalVolume* fWorld_LV;
  G4LogicalVolume* fTank_LV;

  G4Material* fWorldMaterial;
  G4Material* fTankMaterial;

  G4OpticalSurface* fSurface;

  DetectorMessenger* fDetectorMessenger;
//  void DefineMaterials();


  G4MaterialPropertiesTable* fTankMPT;
  G4MaterialPropertiesTable* fWorldMPT;
  G4MaterialPropertiesTable* fSurfaceMPT;
  G4Box  *HDPE_Box1,*Lead_Box,*HDPE_Box2, *sBox,*HDPE_Box3,*HDPE_Box4,*HDPE_Box5,*HDPE_Box6,*HDPE_Box7,*HDPE_Box8,*HDPE_Box9,*HDPE_Box10,*Lead_Box4,*Lead_Box5,*Graphite_Box;
  G4Box  *HDPE_Box11,*HDPE_Box12,*HDPE_Box13,*HDPE_Box14,*HDPE_Box15,*HDPE_Box16,*HDPE_Box17,*Borated_Box1,*solidScintillator,*Lead_Box1,*Lead_Box2,*Lead_Box3,*Graphite_Box2,*Graphite_Box3;
  G4Box  *HDPE_Box18,*HDPE_Box19,*HDPE_Box20,*HDPE_Box21,*HDPE_Box22,*HDPE_Box23,*HDPE_Box24,*HDPE_Box25;
  G4VPhysicalVolume *Lead_PV,*HDPE_PV1,*HDPE_PV2, *HDPE_PV3,*HDPE_PV4,*HDPE_PV5,*HDPE_PV6,*HDPE_PV7,*HDPE_PV8,*HDPE_PV9,*HDPE_PV10,*Lead_PV4,*Lead_PV5,*Graphite_PV,*Graphite_PV2,*Graphite_PV3;
  G4VPhysicalVolume *HDPE_PV11,*HDPE_PV12,*HDPE_PV13,*HDPE_PV14,*HDPE_PV15,*HDPE_PV16,*HDPE_PV17,*Borated_PV1,*physScintillator,*Lead_PV1,*Lead_PV2,*Lead_PV3,*Hole_PV3;
  G4VPhysicalVolume *HDPE_PV18,*HDPE_PV19,*HDPE_PV20,*HDPE_PV21,*HDPE_PV22,*HDPE_PV23,*HDPE_PV24,*HDPE_PV25;
  G4LogicalVolume   *fLBox,*Lead_LV,*HDPE_LV1,*HDPE_LV2,*HDPE_LV3,*HDPE_LV4,*HDPE_LV5,*HDPE_LV6,*HDPE_LV7,*HDPE_LV8,*HDPE_LV9,*HDPE_LV10,*Lead_LV4,*Lead_LV5,*Graphite_LV,*Graphite_LV2,*Graphite_LV3;
  G4LogicalVolume   *HDPE_LV11,*HDPE_LV12,*HDPE_LV13,*HDPE_LV14,*HDPE_LV15,*HDPE_LV16,*HDPE_LV17,*Borated_LV1,*logicScintillator,*Lead_LV1,*Lead_LV2,*Lead_LV3;
  G4LogicalVolume   *HDPE_LV18,*HDPE_LV19,*HDPE_LV20,*HDPE_LV21,*HDPE_LV22,*HDPE_LV23,*HDPE_LV24,*HDPE_LV25;
  G4LogicalVolume   *fScoringVolume_1, *Hole_LV,*Hole_LV3;
      // Store scoring volumes in a vector
  std::vector<G4LogicalVolume*> fScoringVolumes;
  G4SubtractionSolid    *collimator;
  G4VSolid *Hole, *Hole3;

  G4VPhysicalVolume *fPBox, *Hole_PV,*Hole_PV2;

  G4double           fBoxSize, fblockSize,fLeadSize,BoratedSize,Borated_thickness,LeadSize,fGraphiteSize;
  G4Material*        fMaterial;
  G4Material *Air, *b_polyethylene,  *polyethylene, *NaI,*mat_graphite;
  G4Material  *leadMaterial,*Aluminium,*PP;
  G4Element  *Na, *I, *C,*N,*O,*F,*Al;
  G4Material *Vacc;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif /*DetectorConstruction_h*/
