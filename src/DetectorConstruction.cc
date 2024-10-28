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
/// \file optical/OpNovice2/src/DetectorConstruction.cc
/// \brief Implementation of the DetectorConstruction class
//
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "DetectorConstruction.hh"

#include "DetectorMessenger.hh"

#include "G4NistManager.hh"
#include "G4Material.hh"
#include "G4Element.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4OpticalSurface.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4Colour.hh"
#include "G4Tubs.hh"
#include "G4VisAttributes.hh"
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction()
  : G4VUserDetectorConstruction()
  , fDetectorMessenger(nullptr)
{
  fExpHall_x = fExpHall_y = fExpHall_z = 5.0 * m;
  fTank_x = fTank_y =70*cm; 
  fTank_z = 0.3*cm; // 3mm for CF4 gas

  fTank = nullptr;

  fTankMPT    = new G4MaterialPropertiesTable();
  fWorldMPT   = new G4MaterialPropertiesTable();
  fSurfaceMPT = new G4MaterialPropertiesTable();

  fSurface = new G4OpticalSurface("Surface");
  fSurface->SetType(dielectric_metal);
  fSurface->SetFinish(polished);
  fSurface->SetModel(unified);
  fSurface->SetMaterialPropertiesTable(fSurfaceMPT);

  fTank_LV  = nullptr;
  fWorld_LV = nullptr;
  G4double pressure = 0.046*atmosphere; //35 torr
  G4double temperature = 293.15*kelvin; // 
  G4Material* CF4 = new G4Material("CF4", 0.1223*mg/cm3,2,kStateGas,temperature,pressure);
 
  C = G4NistManager::Instance()->FindOrBuildElement("C");
  F  = new G4Element("Fluorine","F",9.,18.998*g/mole);
  CF4->AddElement(C,1);
  CF4->AddElement(F,4);
  const G4int iNbEntries = 300;

  G4double CF4PhotonMomentum[iNbEntries] = {6.2*eV,6.138613861*eV,6.078431373*eV,6.019417476*eV,5.961538462*eV,5.904761905*eV,5.849056604*eV,5.794392523*eV,5.740740741*eV,5.688073394*eV,5.636363636*eV,5.585585586*eV,
                                                5.535714286*eV,5.486725664*eV,5.438596491*eV,5.391304348*eV,5.344827586*eV,5.299145299*eV,5.254237288*eV,5.210084034*eV,5.166666667*eV,5.123966942*eV,5.081967213*eV,5.040650407*eV,5*eV,4.96*eV,4.920634921*eV,
                                                4.881889764*eV,4.84375*eV,4.80620155*eV,4.769230769*eV,4.732824427*eV,4.696969697*eV,4.661654135*eV,4.626865672*eV,4.592592593*eV,4.558823529*eV,4.525547445*eV,4.492753623*eV,4.460431655*eV,
                                                4.428571429*eV,4.397163121*eV,4.366197183*eV,4.335664336*eV,4.305555556*eV,4.275862069*eV,4.246575342*eV,4.217687075*eV,4.189189189*eV,4.161073826*eV,4.133333333*eV,4.105960265*eV,4.078947368*eV,4.052287582*eV,
                                                4.025974026*eV,4*eV,3.974358974*eV,3.949044586*eV,3.924050633*eV,3.899371069*eV,3.875*eV,3.850931677*eV,3.827160494*eV,3.803680982*eV,3.780487805*eV,3.757575758*eV,3.734939759*eV,
                                                3.71257485*eV,3.69047619*eV,3.668639053*eV,3.647058824*eV,3.625730994*eV,3.604651163*eV,3.583815029*eV,3.563218391*eV,3.542857143*eV,3.522727273*eV,3.502824859*eV,3.483146067*eV,
                                                3.463687151*eV,3.444444444*eV,3.425414365*eV,3.406593407*eV,3.387978142*eV,3.369565217*eV,3.351351351*eV,3.333333333*eV,3.315508021*eV,3.29787234*eV,3.28042328*eV,3.263157895*eV,
                                                3.246073298*eV,3.229166667*eV,3.212435233*eV,3.195876289*eV,3.179487179*eV,3.163265306*eV,3.147208122*eV,3.131313131*eV,3.115577889*eV,3.1*eV,3.084577114*eV,3.069306931*eV,3.054187192*eV,
                                                3.039215686*eV,3.024390244*eV,3.009708738*eV,2.995169082*eV,2.980769231*eV,2.966507177*eV,2.952380952*eV,2.938388626*eV,2.924528302*eV,2.910798122*eV,2.897196262*eV,2.88372093*eV,
                                                2.87037037*eV,2.857142857*eV,2.844036697*eV,2.831050228*eV,2.818181818*eV,2.805429864*eV,2.792792793*eV,2.780269058*eV,2.767857143*eV,2.755555556*eV,2.743362832*eV,
                                                2.731277533*eV,2.719298246*eV,2.707423581*eV,2.695652174*eV,2.683982684*eV,2.672413793*eV,2.660944206*eV,2.64957265*eV,2.638297872*eV,2.627118644*eV,2.616033755*eV,2.605042017*eV,
                                                2.594142259*eV,2.583333333*eV,2.572614108*eV,2.561983471*eV,2.551440329*eV,2.540983607*eV,2.530612245*eV,2.520325203*eV,2.510121457*eV,2.5*eV,2.489959839*eV,2.48*eV,2.470119522*eV,2.46031746*eV,
                                                2.450592885*eV,2.440944882*eV,2.431372549*eV,2.421875*eV,2.412451362*eV,2.403100775*eV,2.393822394*eV,2.384615385*eV,2.375478927*eV,2.366412214*eV,2.357414449*eV,
                                                2.348484848*eV,2.339622642*eV,2.330827068*eV,2.322097378*eV,2.313432836*eV,2.304832714*eV,2.296296296*eV,2.287822878*eV,2.279411765*eV,2.271062271*eV,2.262773723*eV,2.254545455*eV,2.246376812*eV,2.238267148*eV,
                                                2.230215827*eV,2.222222222*eV,2.214285714*eV,2.206405694*eV,2.19858156*eV,2.190812721*eV,2.183098592*eV,2.175438596*eV,2.167832168*eV,2.160278746*eV,2.152777778*eV,2.14532872*eV,2.137931034*eV,2.130584192*eV,
                                                2.123287671*eV,2.116040956*eV,2.108843537*eV,2.101694915*eV,2.094594595*eV,2.087542088*eV,2.080536913*eV,2.073578595*eV,2.066666667*eV,2.059800664*eV,2.052980132*eV,2.04620462*eV,2.039473684*eV,2.032786885*eV,
                                                2.026143791*eV,2.019543974*eV,2.012987013*eV,2.006472492*eV,2*eV,1.993569132*eV,1.987179487*eV,1.980830671*eV,1.974522293*eV,1.968253968*eV,1.962025316*eV,1.955835962*eV,
                                                1.949685535*eV,1.943573668*eV,1.9375*eV,1.931464174*eV,1.925465839*eV,1.919504644*eV,1.913580247*eV,1.907692308*eV,1.901840491*eV,1.896024465*eV,1.890243902*eV,1.88449848*eV,1.878787879*eV,1.873111782*eV,
                                                1.86746988*eV,1.861861862*eV,1.856287425*eV,1.850746269*eV,1.845238095*eV,1.839762611*eV,1.834319527*eV,1.828908555*eV,1.823529412*eV,1.818181818*eV,
                                                1.812865497*eV,1.807580175*eV,1.802325581*eV,1.797101449*eV,1.791907514*eV,1.786743516*eV,1.781609195*eV,1.776504298*eV,1.771428571*eV,1.766381766*eV,1.761363636*eV,1.756373938*eV,1.751412429*eV,
                                                1.746478873*eV,1.741573034*eV,1.736694678*eV,1.731843575*eV,1.727019499*eV,1.722222222*eV,1.717451524*eV,1.712707182*eV,1.707988981*eV,1.703296703*eV,1.698630137*eV,1.693989071*eV,1.689373297*eV,1.684782609*eV,
                                                1.680216802*eV,1.675675676*eV,1.67115903*eV,1.666666667*eV,1.662198391*eV,1.657754011*eV,1.653333333*eV,1.64893617*eV,1.644562334*eV,1.64021164*eV,1.635883905*eV,1.631578947*eV,1.627296588*eV,1.623036649*eV,
                                                1.618798956*eV,1.614583333*eV,1.61038961*eV,1.606217617*eV,1.602067183*eV,1.597938144*eV,1.593830334*eV,1.58974359*eV,1.585677749*eV,1.581632653*eV,1.577608142*eV,1.573604061*eV,1.569620253*eV,
                                                1.565656566*eV,1.561712846*eV,1.557788945*eV,1.553884712*eV};

  // Sort the array in ascending order
  std::sort(CF4PhotonMomentum, CF4PhotonMomentum + iNbEntries);
   
 G4double CF4Scintillation_Fast[iNbEntries]    = {0.0029,0.0029,0.0017,0.0024,0.0018,0.0011,0.0027,0.0009,0.0003,0.0019,0.0030,0.0024,0.0023,0.0036,0.0039,0.0056,
                                                0.0049,0.0061,0.0053,0.0052,0.0056,0.0064,0.0072,0.0064,0.0080,0.0071,0.0056,0.0069,0.0053,0.0070,0.0060,0.0057,0.0071,0.0066,0.0066,
                                                0.0055,0.0082,0.0076,0.0093,0.0089,0.0106,0.0109,0.0105,0.0102,0.0120,0.0121,0.0102,0.0097,0.0120,0.0126,0.0097,0.0103,0.0097,0.0084,
                                                0.0119,0.0112,0.0096,0.0171,0.0235,0.0078,0.0089,0.0071,0.0065,0.0074,0.0073,0.0074,0.0074,0.0080,0.0143,0.0522,0.0069,0.0076,0.0042,
                                                0.0059,0.0039,0.0053,0.0054,0.0185,0.0077,0.0599,0.0048,0.0034,0.0041,0.0041,0.0047,0.0059,0.0046,0.0065,0.0128,0.0037,0.0167,0.0053,
                                                0.0038,0.0042,0.0046,0.0032,0.0037,0.0073,0.0049,0.0067,0.0116,0.0054,0.0077,0.0111,0.0042,0.0043,0.0037,0.0046,0.0041,0.0028,0.0055,
                                                0.0031,0.0048,0.0057,0.0056,0.0035,0.0039,0.0068,0.0051,0.0037,0.0054,0.0048,0.0061,0.0033,0.0050,0.0052,0.0047,0.0014,0.0043,0.0041,
                                                0.0023,0.0062,0.0036,0.0038,0.0039,0.0043,0.0049,0.0049,0.0036,0.0048,0.0039,0.0023,0.0035,0.0025,0.0036,0.0010,0.0044,0.0013,0.0041,
                                                0.0021,0.0016,0.0046,0.0040,0.0034,0.0027,0.0026,0.0034,0.0004,0.0037,0.0004,0.0036,0.0029,0.0029,0.0036,0.0055,0.0034,0.0034,0.0025,
                                                0.0028,0.0055,0.0064,0.0037,0.0029,0.0047,0.0058,0.0040,0.0062,0.0055,0.0029,0.0067,0.0070,0.0080,0.0060,0.0094,0.0082,0.0072,0.0089,
                                                0.0117,0.0102,0.0134,0.0131,0.0131,0.0120,0.0135,0.0096,0.0107,0.0179,0.0210,0.0172,0.0165,0.0167,0.0176,0.0137,0.0196,0.0217,0.0175,
                                                0.0223,0.0192,0.0222,0.0188,0.0184,0.0183,0.0156,0.0098,0.0198,0.0268,0.0188,0.0236,0.0208,0.0171,0.0229,0.0228,0.0227,0.0204,0.0184,
                                                0.0190,0.0185,0.0145,0.0138,0.0122,0.0180,0.0132,0.0146,0.0087,0.0039,0.0147,0.0000,0.0000,0.0137,0.0084,0.0094,0.0114,0.0078,0.0100,
                                                0.0069,0.0055,0.0164,0.0113,0.0148,0.0053,0.0054,0.0065,0.0092,0.0000,0.0047,0.0000,0.0071,0.0000,0.0057,0.0063,0.0064,0.0050,0.0077,
                                                0.0034,0.0025,0.0000,0.0041,0.0025,0.0019,0.0042,0.0030,0.0000,0.0030,0.0000,0.0000,0.0000,0.0027,0.0000,0.0000,0.0000,0.0000,0.0006,
                                                0.0051,0.0083,0.0000,0.0000,0.0064,0.0003,0.0002,0.0074,0.0038,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000};

  std::sort(CF4Scintillation_Fast, CF4Scintillation_Fast + iNbEntries);  
  G4double CF4Scintillation_Slow[iNbEntries]    = {0.0029,0.0029,0.0017,0.0024,0.0018,0.0011,0.0027,0.0009,0.0003,0.0019,0.0030,0.0024,0.0023,0.0036,0.0039,0.0056,
                                                0.0049,0.0061,0.0053,0.0052,0.0056,0.0064,0.0072,0.0064,0.0080,0.0071,0.0056,0.0069,0.0053,0.0070,0.0060,0.0057,0.0071,0.0066,0.0066,
                                                0.0055,0.0082,0.0076,0.0093,0.0089,0.0106,0.0109,0.0105,0.0102,0.0120,0.0121,0.0102,0.0097,0.0120,0.0126,0.0097,0.0103,0.0097,0.0084,
                                                0.0119,0.0112,0.0096,0.0171,0.0235,0.0078,0.0089,0.0071,0.0065,0.0074,0.0073,0.0074,0.0074,0.0080,0.0143,0.0522,0.0069,0.0076,0.0042,
                                                0.0059,0.0039,0.0053,0.0054,0.0185,0.0077,0.0599,0.0048,0.0034,0.0041,0.0041,0.0047,0.0059,0.0046,0.0065,0.0128,0.0037,0.0167,0.0053,
                                                0.0038,0.0042,0.0046,0.0032,0.0037,0.0073,0.0049,0.0067,0.0116,0.0054,0.0077,0.0111,0.0042,0.0043,0.0037,0.0046,0.0041,0.0028,0.0055,
                                                0.0031,0.0048,0.0057,0.0056,0.0035,0.0039,0.0068,0.0051,0.0037,0.0054,0.0048,0.0061,0.0033,0.0050,0.0052,0.0047,0.0014,0.0043,0.0041,
                                                0.0023,0.0062,0.0036,0.0038,0.0039,0.0043,0.0049,0.0049,0.0036,0.0048,0.0039,0.0023,0.0035,0.0025,0.0036,0.0010,0.0044,0.0013,0.0041,
                                                0.0021,0.0016,0.0046,0.0040,0.0034,0.0027,0.0026,0.0034,0.0004,0.0037,0.0004,0.0036,0.0029,0.0029,0.0036,0.0055,0.0034,0.0034,0.0025,
                                                0.0028,0.0055,0.0064,0.0037,0.0029,0.0047,0.0058,0.0040,0.0062,0.0055,0.0029,0.0067,0.0070,0.0080,0.0060,0.0094,0.0082,0.0072,0.0089,
                                                0.0117,0.0102,0.0134,0.0131,0.0131,0.0120,0.0135,0.0096,0.0107,0.0179,0.0210,0.0172,0.0165,0.0167,0.0176,0.0137,0.0196,0.0217,0.0175,
                                                0.0223,0.0192,0.0222,0.0188,0.0184,0.0183,0.0156,0.0098,0.0198,0.0268,0.0188,0.0236,0.0208,0.0171,0.0229,0.0228,0.0227,0.0204,0.0184,
                                                0.0190,0.0185,0.0145,0.0138,0.0122,0.0180,0.0132,0.0146,0.0087,0.0039,0.0147,0.0000,0.0000,0.0137,0.0084,0.0094,0.0114,0.0078,0.0100,
                                                0.0069,0.0055,0.0164,0.0113,0.0148,0.0053,0.0054,0.0065,0.0092,0.0000,0.0047,0.0000,0.0071,0.0000,0.0057,0.0063,0.0064,0.0050,0.0077,
                                                0.0034,0.0025,0.0000,0.0041,0.0025,0.0019,0.0042,0.0030,0.0000,0.0030,0.0000,0.0000,0.0000,0.0027,0.0000,0.0000,0.0000,0.0000,0.0006,
                                                0.0051,0.0083,0.0000,0.0000,0.0064,0.0003,0.0002,0.0074,0.0038,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000};

   std::sort(CF4Scintillation_Slow, CF4Scintillation_Slow + iNbEntries);

  const G4int iNbEntries_1 = 3;
  G4double CF4PhotonMomentum_1[iNbEntries_1] = {200*eV,500*eV,700*eV};
  G4double CF4RefractiveIndex[iNbEntries_1]  = {1.004,1.004,1.004};
  G4double CF4AbsorbtionLength[iNbEntries_1] = {100.*cm, 100.*cm, 100.*cm};
  G4double CF4ScatteringLength[iNbEntries_1] = {30.*cm,  30.*cm,  30.*cm};
  G4MaterialPropertiesTable *CF4PropertiesTable = new G4MaterialPropertiesTable();
  CF4PropertiesTable->AddProperty("FASTCOMPONENT", CF4PhotonMomentum, CF4Scintillation_Fast, iNbEntries,true);
  CF4PropertiesTable->AddProperty("SLOWCOMPONENT", CF4PhotonMomentum, CF4Scintillation_Slow, iNbEntries,true);
  CF4PropertiesTable->AddProperty("RINDEX", CF4PhotonMomentum_1, CF4RefractiveIndex, iNbEntries_1);
  CF4PropertiesTable->AddProperty("ABSLENGTH", CF4PhotonMomentum_1, CF4AbsorbtionLength, iNbEntries_1);
  CF4PropertiesTable->AddProperty("RAYLEIGH", CF4PhotonMomentum_1, CF4ScatteringLength, iNbEntries_1);
  CF4PropertiesTable->AddConstProperty("SCINTILLATIONYIELD", 2500./keV,true);  // for electron recoil
  CF4PropertiesTable->AddConstProperty("RESOLUTIONSCALE", 1.0);
  CF4PropertiesTable->AddConstProperty("FASTTIMECONSTANT", 3.*ns,true);
  CF4PropertiesTable->AddConstProperty("SLOWTIMECONSTANT", 10.*ns,true);
  CF4PropertiesTable->AddConstProperty("YIELDRATIO", 1.0,true);
  CF4->SetMaterialPropertiesTable(CF4PropertiesTable);



  fTankMaterial = CF4;
  fWorldMaterial = G4NistManager::Instance()->FindOrBuildMaterial("G4_Galactic");

  fDetectorMessenger = new DetectorMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction() {
  delete fTankMPT;
  delete fWorldMPT;
  delete fSurfaceMPT;
  delete fSurface;
  delete fDetectorMessenger;
}


G4VPhysicalVolume* DetectorConstruction::Construct()
{

  G4NistManager *nist = G4NistManager::Instance();
  G4int ncomponents, natoms;
  G4double massfraction;

  G4double Vdens = 1.e-25*g/cm3;
  G4double Vpres = 1.e-19*pascal;
  G4double Vtemp = 0.1*kelvin;

  G4double a, z;
  C = nist->FindOrBuildElement("C");
  N  = new G4Element("Nitrogen","N",7.,14.007*g/mole);
  O  = new G4Element("Oxygen","O",8.,15.999*g/mole);
  F  = new G4Element("Fluorine","F",9.,18.998*g/mole);

  // boron
  G4Isotope* B10 = new G4Isotope("B10", 5, 10);
  G4Isotope* B11 = new G4Isotope("B11", 5, 11);
  G4Element* B = new G4Element("Boron", "B", ncomponents=2);
  B->AddIsotope(B10, 19.9*perCent);
  B->AddIsotope(B11, 80.1*perCent);
  G4Material* boron = new G4Material("boron", 2.46*g/cm3, ncomponents=1, kStateSolid,293*kelvin, 1*atmosphere);
  boron->AddElement(B, natoms=1);

   // pressurized water
  G4Element* H  = new G4Element("TS_H_of_Water" ,"H" , 1., 1.0079*g/mole);
  G4Material* H2O =
  new G4Material("Water_ts", 1.000*g/cm3, ncomponents=2,
                         kStateLiquid, 593*kelvin, 150*bar);
  H2O->AddElement(H, natoms=2);
  H2O->AddElement(O, natoms=1);
  H2O->GetIonisation()->SetMeanExcitationEnergy(78.0*eV);
  // vacuum
  Vacc = new G4Material("Galactic", z=1, a=1.01*g/mole, Vdens, kStateGas, Vtemp, Vpres);

  //Graphite 
  mat_graphite = nist->FindOrBuildMaterial("G4_GRAPHITE");
  // air
  G4Element* N = new G4Element("Nitrogen", "N", 7., 14.01*g/mole);
  Air = new G4Material("air", 1.290*mg/cm3, ncomponents=2, kStateGas, 293*kelvin, 1*atmosphere);
  Air->AddElement(N, massfraction=70.*perCent);
  Air->AddElement(O, massfraction=30.*perCent);

  // polyethilene
  G4Element* Hpe = new G4Element("TS_H_of_Polyethylene", "H", 1, 1.0079*g/mole);
  G4Element* Cpe = new G4Element("Carbon", "C", 6, 12.01*g/mole);
  polyethylene = new G4Material("polyethylene", 0.93*g/cm3, ncomponents=2, kStateSolid, 293*kelvin, 1*atmosphere);
  polyethylene->AddElement(Hpe, natoms=4);
  polyethylene->AddElement(Cpe, natoms=2);

  // borated polyethilene
  b_polyethylene = new G4Material("b_polyethylene",0.94*g/cm3,ncomponents=4,kStateSolid,293*kelvin,1*atmosphere);
  b_polyethylene->AddElement(Hpe, 11.6*perCent);
  b_polyethylene->AddElement(Cpe, 61.2*perCent);
  b_polyethylene->AddElement(B, 5*perCent);
  b_polyethylene->AddElement(O, 22.2*perCent);
  

  // Define the lead material
  leadMaterial = new G4Material("Lead", 82, 207.2 * g/mole, 11.35 * g/cm3);


    // ------------------------------------ Polypropilene ------------------------------------

  nist->FindOrBuildMaterial("G4_POLYPROPYLENE");
  PP = G4Material::GetMaterial("G4_POLYPROPYLENE");

  
//----------------------------------- Aluminium ------------------------------------
 
   G4double density;
   Aluminium = new G4Material("Al", z = 13., a = 26.98 * g / mole,
                        density = 2.7 * g / cm3);
 


//....................End of scintillator material........................................


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

  fTankMaterial->SetMaterialPropertiesTable(fTankMPT);
  fTankMaterial->GetIonisation()->SetBirksConstant(0.126 * mm / MeV);

  fWorldMaterial->SetMaterialPropertiesTable(fWorldMPT);

  // ------------- Volumes --------------
  // The experimental Hall
  G4Box* world_box = new G4Box("World", fExpHall_x, fExpHall_y, fExpHall_z);

  fWorld_LV = new G4LogicalVolume(world_box, fWorldMaterial, "World");

  G4VPhysicalVolume* world_PV =
    new G4PVPlacement(0, G4ThreeVector(), fWorld_LV, "World", 0, false, 0);

  //The HDPE_block1

  fblockSize = 10*cm;


  HDPE_Box1 = new G4Box("HDPE1",                             //its name
                   10*cm/2,10*cm/2,5*cm/2);   //its dimensions

  HDPE_LV1 = new G4LogicalVolume(HDPE_Box1,                     //its shape
                              polyethylene,                      //its material
                             "HDPE1");                  //its name

  HDPE_PV1 = new G4PVPlacement(0,                          //no rotation
                            G4ThreeVector(0,0,2.5*cm),            //at (0,0,0)
                             HDPE_LV1,                      //its logical volume
                            "HDPE1",                    //its name
                            fWorld_LV,                          //its mother  volume
                            false,                      //no boolean operation
                            0,true);                         //copy number

//The HDPE_block2


  HDPE_Box2 = new G4Box("HDPE2",                             //its name
                   fblockSize/2,fblockSize/2,fblockSize/2);   //its dimensions

  HDPE_LV2 = new G4LogicalVolume(HDPE_Box2,                     //its shape
                              polyethylene,                      //its material
                             "HDPE2");                  //its name

  HDPE_PV2 = new G4PVPlacement(0,                          //no rotation
                            G4ThreeVector(0,10*cm,5*cm),            //at (0,0,0)
                             HDPE_LV2,                      //its logical volume
                            "HDPE2",                    //its name
                            fWorld_LV,                          //its mother  volume
                            false,                      //no boolean operation
                            0,true);                         //copy number

 //The HDPE_block3


  HDPE_Box3 = new G4Box("HDPE3",                             //its name
                   fblockSize/2,fblockSize/2,fblockSize/2);   //its dimensions

  HDPE_LV3 = new G4LogicalVolume(HDPE_Box3,                     //its shape
                              polyethylene,                      //its material
                             "HDPE3");                  //its name

  HDPE_PV3 = new G4PVPlacement(0,                          //no rotation
                            G4ThreeVector(0,-10*cm,5*cm),            //at (0,0,0)
                             HDPE_LV3,                      //its logical volume
                            "HDPE3",                    //its name
                            fWorld_LV,                          //its mother  volume
                            false,                      //no boolean operation
                            0,true);                         //copy number

   //The HDPE_block4


  HDPE_Box4 = new G4Box("HDPE4",                             //its name
                   fblockSize/2,fblockSize/2,fblockSize/2);   //its dimensions

  HDPE_LV4 = new G4LogicalVolume(HDPE_Box4,                     //its shape
                              polyethylene,                      //its material
                             "HDPE4");                  //its name

  HDPE_PV4 = new G4PVPlacement(0,                          //no rotation
                            G4ThreeVector(0,10*cm,15*cm),            //at (0,0,0)
                             HDPE_LV4,                      //its logical volume
                            "HDPE4",                    //its name
                            fWorld_LV,                          //its mother  volume
                            false,                      //no boolean operation
                            0,true);                         //copy number

//The HDPE_block5


  HDPE_Box5 = new G4Box("HDPE5",                             //its name
                   fblockSize/2,fblockSize/2,fblockSize/2);   //its dimensions

  HDPE_LV5 = new G4LogicalVolume(HDPE_Box5,                     //its shape
                              polyethylene,                      //its material
                             "HDPE5");                  //its name

  HDPE_PV5 = new G4PVPlacement(0,                          //no rotation
                            G4ThreeVector(0,-10*cm,15*cm),            //at (0,0,0)
                             HDPE_LV5,                      //its logical volume
                            "HDPE5",                    //its name
                             fWorld_LV,                          //its mother  volume
                            false,                      //no boolean operation
                            0,true);                         //copy number
 //The HDPE_block6


  HDPE_Box6 = new G4Box("HDPE6",                             //its name
                   fblockSize/2,fblockSize/2,fblockSize/2);   //its dimensions

  HDPE_LV6 = new G4LogicalVolume(HDPE_Box6,                     //its shape
                              polyethylene,                      //its material
                             "HDPE6");                  //its name

  HDPE_PV6 = new G4PVPlacement(0,                          //no rotation
                            G4ThreeVector(10*cm,0,5*cm),            //at (0,0,0)
                             HDPE_LV6,                      //its logical volume
                            "HDPE6",                    //its name
                             fWorld_LV,                          //its mother  volume
                            false,                      //no boolean operation
                            0,true);                         //copy number

  //The HDPE_block7


  HDPE_Box7 = new G4Box("HDPE7",                             //its name
                   fblockSize/2,fblockSize/2,fblockSize/2);   //its dimensions

  HDPE_LV7 = new G4LogicalVolume(HDPE_Box7,                     //its shape
                              polyethylene,                      //its material
                             "HDPE7");                  //its name


  HDPE_PV7 = new G4PVPlacement(0,                          //no rotation
                            G4ThreeVector(-10*cm,0,5*cm),            //at (0,0,0)
                             HDPE_LV7,                      //its logical volume
                            "HDPE7",                    //its name
                             fWorld_LV,                          //its mother  volume
                            false,                      //no boolean operation
                            0,true);                         //copy number


  //The HDPE_block8


  HDPE_Box8 = new G4Box("HDPE8",                             //its name
                   fblockSize/2,fblockSize/2,fblockSize/2);   //its dimensions

  HDPE_LV8 = new G4LogicalVolume(HDPE_Box8,                     //its shape
                              polyethylene,                      //its material
                             "HDPE8");                  //its name

  HDPE_PV8 = new G4PVPlacement(0,                          //no rotation
                            G4ThreeVector(10*cm,10*cm,5*cm),            //at (0,0,0)
                             HDPE_LV8,                      //its logical volume
                            "HDPE8",                    //its name
                             fWorld_LV,                          //its mother  volume
                            false,                      //no boolean operation
                            0,true);                         //copy number


     //The HDPE_block9


  HDPE_Box9 = new G4Box("HDPE9",                             //its name
                   fblockSize/2,fblockSize/2,fblockSize/2);   //its dimensions

  HDPE_LV9 = new G4LogicalVolume(HDPE_Box9,                     //its shape
                              polyethylene,                      //its material
                             "HDPE9");                  //its name

  HDPE_PV9 = new G4PVPlacement(0,                          //no rotation
                            G4ThreeVector(10*cm,-10*cm,5*cm),            //at (0,0,0)
                             HDPE_LV9,                      //its logical volume
                            "HDPE9",                    //its name
                             fWorld_LV,                          //its mother  volume
                            false,                      //no boolean operation
                            0,true);                         //copy number

//The HDPE_block10

  HDPE_Box10 = new G4Box("HDPE10",                             //its name
                   fblockSize/2,fblockSize/2,fblockSize/2);   //its dimensions

  HDPE_LV10 = new G4LogicalVolume(HDPE_Box10,                     //its shape
                              polyethylene,                      //its material
                             "HDPE10");                  //its name

  HDPE_PV10 = new G4PVPlacement(0,                          //no rotation
                            G4ThreeVector(-10*cm,10*cm,5*cm),            //at (0,0,0)
                             HDPE_LV10,                      //its logical volume
                            "HDPE10",                    //its name
                             fWorld_LV,                          //its mother  volume
                            false,                      //no boolean operation
                            0,true);                         //copy number

   //The HDPE_block11

  HDPE_Box11 = new G4Box("HDPE11",                             //its name
                   fblockSize/2,fblockSize/2,fblockSize/2);   //its dimensions

  HDPE_LV11 = new G4LogicalVolume(HDPE_Box11,                     //its shape
                              polyethylene,                      //its material
                             "HDPE11");                  //its name

  HDPE_PV11 = new G4PVPlacement(0,                          //no rotation
                            G4ThreeVector(-10*cm,-10*cm,5*cm),            //at (0,0,0)
                             HDPE_LV11,                      //its logical volume
                            "HDPE11",                    //its name
                             fWorld_LV,                          //its mother  volume
                            false,                      //no boolean operation
                            0,true);                         //copy number


//The HDPE_block12


  HDPE_Box12 = new G4Box("HDPE12",                             //its name
                   fblockSize/2,fblockSize/2,fblockSize/2);   //its dimensions

  HDPE_LV12 = new G4LogicalVolume(HDPE_Box12,                     //its shape
                              polyethylene,                      //its material
                             "HDPE12");                  //its name

  HDPE_PV12 = new G4PVPlacement(0,                          //no rotation
                            G4ThreeVector(-10*cm,-10*cm,15*cm),            //at (0,0,0)
                             HDPE_LV12,                      //its logical volume
                            "HDPE12",                    //its name
                             fWorld_LV,                          //its mother  volume
                            false,                      //no boolean operation
                            0,true);                         //copy number
    
   //The HDPE_block13


  HDPE_Box13 = new G4Box("HDPE13",                             //its name
                   fblockSize/2,fblockSize/2,fblockSize/2);   //its dimensions

  HDPE_LV13 = new G4LogicalVolume(HDPE_Box13,                     //its shape
                              polyethylene,                      //its material
                             "HDPE13");                  //its name

  HDPE_PV13 = new G4PVPlacement(0,                          //no rotation
                            G4ThreeVector(10*cm,10*cm,15*cm),            //at (0,0,0)
                             HDPE_LV13,                      //its logical volume
                            "HDPE13",                    //its name
                            fWorld_LV,                          //its mother  volume
                            false,                      //no boolean operation
                            0,true);                         //copy number

//The HDPE_block14


  HDPE_Box14 = new G4Box("HDPE14",                             //its name
                   fblockSize/2,fblockSize/2,fblockSize/2);   //its dimensions

  HDPE_LV14 = new G4LogicalVolume(HDPE_Box14,                     //its shape
                              polyethylene,                      //its material
                             "HDPE14");                  //its name

  HDPE_PV14 = new G4PVPlacement(0,                          //no rotation
                            G4ThreeVector(10*cm,-10*cm,15*cm),            //at (0,0,0)
                             HDPE_LV14,                      //its logical volume
                            "HDPE14",                    //its name
                            fWorld_LV,                          //its mother  volume
                            false,                      //no boolean operation
                            0,true);                         //copy number



  //The HDPE_block15

  HDPE_Box15 = new G4Box("HDPE15",                             //its name
                   fblockSize/2,fblockSize/2,fblockSize/2);   //its dimensions

  HDPE_LV15 = new G4LogicalVolume(HDPE_Box15,                     //its shape
                              polyethylene,                      //its material
                             "HDPE15");                  //its name

  HDPE_PV15 = new G4PVPlacement(0,                          //no rotation
                            G4ThreeVector(-10*cm,10*cm,15*cm),            //at (0,0,0)
                             HDPE_LV15,                      //its logical volume
                            "HDPE15",                    //its name
                            fWorld_LV,                          //its mother  volume
                            false,                      //no boolean operation
                            0,true);                         //copy number

//The HDPE_block16


  HDPE_Box16 = new G4Box("HDPE16",                             //its name
                   fblockSize/2,fblockSize/2,fblockSize/2);   //its dimensions

  HDPE_LV16 = new G4LogicalVolume(HDPE_Box16,                     //its shape
                              polyethylene,                      //its material
                             "HDPE16");                  //its name

  HDPE_PV16 = new G4PVPlacement(0,                          //no rotation
                            G4ThreeVector(-10*cm,0*cm,15*cm),            //at (0,0,0)
                             HDPE_LV16,                      //its logical volume
                            "HDPE16",                    //its name
                            fWorld_LV,                          //its mother  volume
                            false,                      //no boolean operation
                            0,true);                         //copy number

  //The HDPE_block17

  HDPE_Box17 = new G4Box("HDPE17",                             //its name
                   fblockSize/2,fblockSize/2,fblockSize/2);   //its dimensions

  HDPE_LV17 = new G4LogicalVolume(HDPE_Box17,                     //its shape
                              polyethylene,                      //its material
                             "HDPE17");                  //its name

  HDPE_PV17 = new G4PVPlacement(0,                          //no rotation
                            G4ThreeVector(10*cm,0*cm,15*cm),            //at (0,0,0)
                             HDPE_LV17,                      //its logical volume
                            "HDPE17",                    //its name
                             fWorld_LV,                          //its mother  volume
                            false,                      //no boolean operation
                            0,true);                         //copy number

    //The HDPE_block18


  HDPE_Box18 = new G4Box("HDPE18",                             //its name
                   fblockSize/2,fblockSize/2,10*cm/2);   //its dimensions

  HDPE_LV18 = new G4LogicalVolume(HDPE_Box18,                     //its shape
                              polyethylene,                      //its material
                             "HDPE18");                  //its name

  HDPE_PV18 = new G4PVPlacement(0,                          //no rotation
                            G4ThreeVector(10*cm,-10*cm,25*cm),            //at (0,0,0)
                             HDPE_LV18,                      //its logical volume
                            "HDPE18",                    //its name
                             fWorld_LV,                          //its mother  volume
                            false,                      //no boolean operation
                            0,true);                         //copy number

    //The HDPE_block19


  HDPE_Box19 = new G4Box("HDPE19",                             //its name
                   fblockSize/2,fblockSize/2,10*cm/2);   //its dimensions

  HDPE_LV19 = new G4LogicalVolume(HDPE_Box19,                     //its shape
                              polyethylene,                      //its material
                             "HDPE19");                  //its name

  HDPE_PV19 = new G4PVPlacement(0,                          //no rotation
                            G4ThreeVector(-10*cm,-10*cm,25*cm),            //at (0,0,0)
                             HDPE_LV19,                      //its logical volume
                            "HDPE19",                    //its name
                             fWorld_LV,                          //its mother  volume
                            false,                      //no boolean operation
                            0,true);

     //The HDPE_block20


  HDPE_Box20 = new G4Box("HDPE20",                             //its name
                   fblockSize/2,fblockSize/2,10*cm/2);   //its dimensions

  HDPE_LV20 = new G4LogicalVolume(HDPE_Box20,                     //its shape
                              polyethylene,                      //its material
                             "HDPE20");                  //its name

  HDPE_PV20 = new G4PVPlacement(0,                          //no rotation
                            G4ThreeVector(10*cm,10*cm,25*cm),            //at (0,0,0)
                             HDPE_LV20,                      //its logical volume
                            "HDPE20",                    //its name
                             fWorld_LV,                          //its mother  volume
                            false,                      //no boolean operation
                            0,true);                         //copy number

 //The HDPE_block21

  HDPE_Box21 = new G4Box("HDPE21",                             //its name
                   fblockSize/2,fblockSize/2,10*cm/2);   //its dimensions

  HDPE_LV21 = new G4LogicalVolume(HDPE_Box21,                     //its shape
                              polyethylene,                      //its material
                             "HDPE21");                  //its name

  HDPE_PV21 = new G4PVPlacement(0,                          //no rotation
                            G4ThreeVector(-10*cm,10*cm,25*cm),            //at (0,0,0)
                             HDPE_LV21,                      //its logical volume
                            "HDPE21",                    //its name
                             fWorld_LV,                          //its mother  volume
                            false,                      //no boolean operation
                            0,true);                         //copy number

  //The HDPE_block22


  HDPE_Box22 = new G4Box("HDPE22",                             //its name
                   fblockSize/2,fblockSize/2,10*cm/2);   //its dimensions

  HDPE_LV22 = new G4LogicalVolume(HDPE_Box22,                     //its shape
                              polyethylene,                      //its material
                             "HDPE22");                  //its name

  HDPE_PV22 = new G4PVPlacement(0,                          //no rotation
                            G4ThreeVector(-10*cm,0*cm,25*cm),            //at (0,0,0)
                             HDPE_LV22,                      //its logical volume
                            "HDPE22",                    //its name
                            fWorld_LV,                          //its mother  volume
                            false,                      //no boolean operation
                            0,true);                         //copy number

//The HDPE_block23

  HDPE_Box23 = new G4Box("HDPE23",                             //its name
                   fblockSize/2,fblockSize/2,10*cm/2);   //its dimensions

  HDPE_LV23 = new G4LogicalVolume(HDPE_Box23,                     //its shape
                              polyethylene,                      //its material
                             "HDPE23");                  //its name

  HDPE_PV23 = new G4PVPlacement(0,                          //no rotation
                            G4ThreeVector(10*cm,0*cm,25*cm),            //at (0,0,0)
                             HDPE_LV23,                      //its logical volume
                            "HDPE23",                    //its name
                            fWorld_LV,                          //its mother  volume
                            false,                      //no boolean operation
                            0,true);                         //copy number


 //The HDPE_block24
  

  HDPE_Box24 = new G4Box("HDPE24",                             //its name
                   fblockSize/2,fblockSize/2,10*cm/2);   //its dimensions

  HDPE_LV24 = new G4LogicalVolume(HDPE_Box24,                     //its shape
                              polyethylene,                      //its material
                             "HDPE24");                  //its name

  HDPE_PV24 = new G4PVPlacement(0,                          //no rotation
                            G4ThreeVector(0,10*cm,25*cm),            //at (0,0,0)
                             HDPE_LV24,                      //its logical volume
                            "HDPE24",                    //its name
                            fWorld_LV,                          //its mother  volume
                            false,                      //no boolean operation
                            0,true);                         //copy number

//The HDPE_block25


  HDPE_Box25 = new G4Box("HDPE25",                             //its name
                   fblockSize/2,fblockSize/2,10*cm/2);   //its dimensions

  HDPE_LV25 = new G4LogicalVolume(HDPE_Box25,                     //its shape
                              polyethylene,                      //its material
                             "HDPE25");                  //its name

  HDPE_PV25 = new G4PVPlacement(0,                          //no rotation
                            G4ThreeVector(0,-10*cm,25*cm),            //at (0,0,0)
                             HDPE_LV25,                      //its logical volume
                            "HDPE25",                    //its name
                            fWorld_LV,                          //its mother  volume
                            false,                      //no boolean operation
                            0,true);                         //copy number

  //The lead1
  fLeadSize = 10*cm;


  Lead_Box = new G4Box("Lead",                             //its name
                   fLeadSize/2,fLeadSize/2, 5*cm/2);   //its dimensions

  Lead_LV = new G4LogicalVolume(Lead_Box,                     //its shape
                              leadMaterial,                      //its material
                             "Lead");                  //its name

  Lead_PV = new G4PVPlacement(0,                          //no rotation
                            G4ThreeVector(0,0,7.5*cm),            //at (0,0,0)
                             Lead_LV,                      //its logical volume
                            "Lead",                    //its name
                            fWorld_LV,                          //its mother  volume
                            false,                      //no boolean operation
                            0,true);                         //copy number


   G4VisAttributes* red = new G4VisAttributes(G4Colour::Red());

   red->SetVisibility(true);
   red->SetForceAuxEdgeVisible(true);

   Lead_LV->SetVisAttributes(red);


  //.........................................................................................................................................
  //The Borated polythylene_block1 with pinhole

  BoratedSize = 30*cm;
  Borated_thickness = 3*cm;
  Borated_Box1 = new G4Box("Borated1",                             //its name
                   BoratedSize/2,  BoratedSize/2,Borated_thickness/2);   //its dimensions



  Hole = new G4Tubs("BoxHole", 0.0*cm, 0.75*cm, 1.5*cm, 0*deg, 360*deg);  // the diameter of the exit(pinhole) is 3cm. In G4 we use halfsize of the radius. 

  Hole_LV = new G4LogicalVolume(Hole,                     //its shape
                              Vacc,                      //its material
                             "H1");                  //its name

   Borated_LV1 = new G4LogicalVolume(Borated_Box1,                     //its shape
                              b_polyethylene,                      //its material
                             "Borated1", 0,0,0);                  //its name



  Borated_PV1 = new G4PVPlacement(0,                          //no rotation
                            G4ThreeVector(0*cm,0*cm,31.5*cm),            //at (0,0,0)
                             Borated_LV1,                      //its logical volume
                            "Borated1",                    //its name
                            fWorld_LV,                          //its mother  volume
                            false,                      //no boolean operation
                            0,true);                         //copy number

  Hole_PV = new G4PVPlacement(0,                          //no rotation
                            G4ThreeVector(0*cm,0*cm,0*cm),            //at (0,0,0)
                             Hole_LV,                      //its logical volume
                            "H1",                    //its name
                            Borated_LV1,                          //its mother  volume
                            false,                      //no boolean operation
                            0,true);                         //copy number


  G4VisAttributes* green = new G4VisAttributes(G4Colour::Green());

   green->SetVisibility(true);
   green->SetForceAuxEdgeVisible(true);

   Borated_LV1->SetVisAttributes(green);

     //The lead2
  

  LeadSize = 3*cm;
  Lead_Box2 = new G4Box("Lead2",                             //its name
                   30*cm/2,LeadSize/2,33*cm/2);   //its dimensions

  Lead_LV2 = new G4LogicalVolume(Lead_Box2,                     //its shape
                              leadMaterial,                      //its material
                             "Lead2");                  //its name

  Lead_PV2 = new G4PVPlacement(0,                          //no rotation
                            G4ThreeVector(0*cm,16.5*cm,16.5*cm),            //at (0,0,0)
                             Lead_LV2,                      //its logical volume
                            "Lead2",                    //its name
                            fWorld_LV,                          //its mother  volume
                            false,                      //no boolean operation
                            0,true);                         //copy number


  G4VisAttributes* yellow= new G4VisAttributes(G4Colour::Yellow());

  yellow->SetVisibility(true);
  yellow->SetForceAuxEdgeVisible(true);

  Lead_LV2->SetVisAttributes(red);



     //The lead3
 
  Lead_Box3 = new G4Box("Lead3",                             //its name
                   3*cm/2,30*cm/2, 33*cm/2);   //its dimensions

  Lead_LV3 = new G4LogicalVolume(Lead_Box3,                     //its shape
                              leadMaterial,                      //its material
                             "Lead3");                  //its name

  Lead_PV3 = new G4PVPlacement(0,                          //no rotation
                            G4ThreeVector(16.5*cm,0*cm,16.5*cm),            //at (0,0,0)
                             Lead_LV3,                      //its logical volume
                            "Lead3",                    //its name
                            fWorld_LV,                          //its mother  volume
                            false,                      //no boolean operation
                            0,true);                         //copy number


   
  Lead_LV3->SetVisAttributes(red);

  //The lead4

  Lead_Box4 = new G4Box("Lead4",                             //its name
                   3*cm/2,30*cm/2, 33*cm/2);   //its dimensions

  Lead_LV4 = new G4LogicalVolume(Lead_Box4,                     //its shape
                              leadMaterial,                      //its material
                             "Lead4");                  //its name

  Lead_PV4 = new G4PVPlacement(0,                          //no rotation
                            G4ThreeVector(-16.5*cm,0*cm,16.5*cm),            //at (0,0,0)
                             Lead_LV4,                      //its logical volume
                            "Lead4",                    //its name
                            fWorld_LV,                          //its mother  volume
                            false,                      //no boolean operation
                            0,true);                         //copy number


   
  Lead_LV4->SetVisAttributes(red);


  //The lead5

  Lead_Box5 = new G4Box("Lead5",                             //its name
                   30*cm/2,30*cm/2, 3*cm/2);   //its dimensions

  Lead_LV5 = new G4LogicalVolume(Lead_Box5,                     //its shape
                              leadMaterial,                      //its material
                             "Lead5");                  //its name

  Lead_PV5 = new G4PVPlacement(0,                          //no rotation
                            G4ThreeVector(0*cm,0*cm,34.5*cm),            //at (0,0,0)
                             Lead_LV5,                      //its logical volume
                            "Lead5",                    //its name
                             fWorld_LV,                          //its mother  volume
                            false,                      //no boolean operation
                            0,true);                         //copy number
  Hole_PV2 = new G4PVPlacement(0,                          //no rotation
                            G4ThreeVector(0*cm,0*cm,0*cm),            //at (0,0,0)
                             Hole_LV,                      //its logical volume
                            "H2",                    //its name
                             Lead_LV5,                          //its mother  volume
                            false,                      //no boolean operation
                            0,true);                         //copy number




   
  Lead_LV5->SetVisAttributes(red);


  // Define dimensions and materials
  G4double shieldThickness = 0.001*cm; // 10 microns
  G4double ppThickness = 0.0001*cm;    // 1 micrometer
  G4double coatingThickness = 0.00001*cm; // 0.1 micrometer
  G4double gasThickness = 0.3*cm; // 3mm for CF4 gas
  G4double size = 70*cm;

  // Loop for 100 iterations

  //for (int i = 0; i < 1; ++i) {
    // Calculate z position for this iteration
    G4double zStart = 112.0*cm  ; // Adjust z position for each iteration

    // Create and place shield (polyethylene)
    auto shield = new G4Box("shield", size/2, size/2, shieldThickness/2);
    auto lShield = new G4LogicalVolume(shield, polyethylene, "Shield");

    auto pShield = new G4PVPlacement(0,
                                      G4ThreeVector(0.*cm, 0.*cm, zStart),
                                      lShield,
                                      "Shield",
                                      fWorld_LV,
                                      false,
                                      0, true);

    // Create and place polypropylene foil
    G4Box* foilSolid = new G4Box("foilSolid", size/2, size/2, ppThickness/2);
    G4LogicalVolume* PPCoating = new G4LogicalVolume(foilSolid, PP, "PPCoating");

    G4double zFoil = zStart + shieldThickness + ppThickness/2; // Place foil right after shield
    G4PVPlacement* foilCoating = new G4PVPlacement(0,
                                                   G4ThreeVector(0.*cm, 0.*cm, zFoil),
                                                   PPCoating,
                                                   "PPCoating",
                                                   fWorld_LV,
                                                   false,
                                                   0, true);

    PPCoating->SetVisAttributes(green);

    // Create and place aluminum coating
    G4Box* coatingSolid = new G4Box("coatingSolid", size/2, size/2, coatingThickness/2);
    G4LogicalVolume* lCoating = new G4LogicalVolume(coatingSolid, Aluminium, "AluminiumCoating");

    G4double zCoating = zFoil + ppThickness/2 + coatingThickness/2; // Place aluminum coating after foil
    auto pCoating = new G4PVPlacement(0,
                                      G4ThreeVector(0.*cm, 0.*cm, zCoating),
                                      lCoating,
                                      "AluminiumCoating",
                                      fWorld_LV,
                                      false,
                                      0, true);

    lCoating->SetVisAttributes(red);

    // Create and place CF4 gas


    G4Box* tank_box = new G4Box("Tank", fTank_x/2, fTank_y/2, fTank_z/2);
    fTank_LV = new G4LogicalVolume(tank_box,  fTankMaterial  , "Tank");

    G4double zGas = zCoating + coatingThickness/2 + gasThickness/2; // Place CF4 gas after coating
    fTank = new G4PVPlacement(0,
                                     G4ThreeVector(0.*cm, 0.*cm, zGas),
                                     fTank_LV,
                                     "fPScore",
                                     fWorld_LV,
                                     false,
                                     0, true);

    //fScoringVolume_1 = fLScore;
    // After creating each scoring volume, push it to the vector
    //fScoringVolumes.push_back(fLScore);


    // Place a second aluminum coating
    G4Box* coatingSolid2 = new G4Box("coatingSolid2", size/2, size/2, coatingThickness/2);
    G4LogicalVolume* lCoating2 = new G4LogicalVolume(coatingSolid2, Aluminium, "AluminiumCoating2");

    G4double zCoating2 = zGas + gasThickness/2 + coatingThickness/2; // Place second aluminum coating
    auto pCoating2 = new G4PVPlacement(0,
                                       G4ThreeVector(0.*cm, 0.*cm, zCoating2),
                                       lCoating2,
                                       "AluminiumCoating2",
                                       fWorld_LV,
                                       false,
                                       0, true);

    lCoating2->SetVisAttributes(red);

    // Place a second polypropylene foil
    G4Box* foilSolid2 = new G4Box("foilSolid2", size/2, size/2, ppThickness/2);
    G4LogicalVolume* PPCoating2 = new G4LogicalVolume(foilSolid2, PP, "PPCoating2");

    G4double zFoil2 = zCoating2 + coatingThickness/2 + ppThickness/2; // Place second polypropylene foil
    G4PVPlacement* foilCoating2 = new G4PVPlacement(0,
                                                    G4ThreeVector(0.*cm, 0.*cm, zFoil2),
                                                    PPCoating2,
                                                    "PPCoating2",
                                                    fWorld_LV,
                                                    false,
                                                    0, true);

    PPCoating2->SetVisAttributes(green);


    //surfaces

    // surface reflecting
    G4OpticalSurface* oppac_Al_gas = new G4OpticalSurface("Reflecting");
    oppac_Al_gas->SetModel(unified);
    oppac_Al_gas->SetType(dielectric_metal);
    oppac_Al_gas->SetFinish(polished);
    G4LogicalBorderSurface* oppac_1 = new G4LogicalBorderSurface("oppac_1", fTank,pCoating2,oppac_Al_gas);
    G4LogicalBorderSurface* oppac_2 = new G4LogicalBorderSurface("oppac_2", fTank,pCoating,oppac_Al_gas);










/* 


  // ------------- Surface --------------

  G4LogicalBorderSurface* surface =
    new G4LogicalBorderSurface("Surface", fTank, pCoating2, fSurface);

  G4LogicalBorderSurface* surface_2 = new G4LogicalBorderSurface("surface_2", fTank,pCoating,fSurface);

  G4OpticalSurface* opticalSurface = dynamic_cast<G4OpticalSurface*>(
    surface->GetSurface(fTank, pCoating2)->GetSurfaceProperty());
  G4cout << "******  opticalSurface->DumpInfo:" << G4endl;
  if(opticalSurface)
  {
    opticalSurface->DumpInfo();
  }
  G4cout << "******  end of opticalSurface->DumpInfo" << G4endl;
*/
  return world_PV;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstruction::SetSurfaceSigmaAlpha(G4double v)
{
  fSurface->SetSigmaAlpha(v);
  G4RunManager::GetRunManager()->GeometryHasBeenModified();

  G4cout << "Surface sigma alpha set to: " << fSurface->GetSigmaAlpha()
         << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstruction::SetSurfacePolish(G4double v)
{
  fSurface->SetPolish(v);
  G4RunManager::GetRunManager()->GeometryHasBeenModified();

  G4cout << "Surface polish set to: " << fSurface->GetPolish() << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstruction::AddTankMPV(const G4String& prop,
                                      G4MaterialPropertyVector* mpv)
{
  fTankMPT->AddProperty(prop, mpv);
  G4cout << "The MPT for the box is now: " << G4endl;
  fTankMPT->DumpTable();
  G4cout << "............." << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstruction::AddWorldMPV(const G4String& prop,
                                       G4MaterialPropertyVector* mpv)
{
  fWorldMPT->AddProperty(prop, mpv);
  G4cout << "The MPT for the world is now: " << G4endl;
  fWorldMPT->DumpTable();
  G4cout << "............." << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstruction::AddSurfaceMPV(const G4String& prop,
                                         G4MaterialPropertyVector* mpv)
{
  fSurfaceMPT->AddProperty(prop, mpv);
  G4cout << "The MPT for the surface is now: " << G4endl;
  fSurfaceMPT->DumpTable();
  G4cout << "............." << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstruction::AddTankMPC(const G4String& prop, G4double v)
{
  fTankMPT->AddConstProperty(prop, v);
  G4cout << "The MPT for the box is now: " << G4endl;
  fTankMPT->DumpTable();
  G4cout << "............." << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstruction::AddWorldMPC(const G4String& prop, G4double v)
{
  fWorldMPT->AddConstProperty(prop, v);
  G4cout << "The MPT for the world is now: " << G4endl;
  fWorldMPT->DumpTable();
  G4cout << "............." << G4endl;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstruction::AddSurfaceMPC(const G4String& prop, G4double v)
{
  fSurfaceMPT->AddConstProperty(prop, v);
  G4cout << "The MPT for the surface is now: " << G4endl;
  fSurfaceMPT->DumpTable();
  G4cout << "............." << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstruction::SetWorldMaterial(const G4String& mat)
{
  G4Material* pmat = G4NistManager::Instance()->FindOrBuildMaterial(mat);
  if(pmat && fWorldMaterial != pmat)
  {
    fWorldMaterial = pmat;
    if(fWorld_LV)
    {
      fWorld_LV->SetMaterial(fWorldMaterial);
      fWorldMaterial->SetMaterialPropertiesTable(fWorldMPT);
    }
    G4RunManager::GetRunManager()->PhysicsHasBeenModified();
    G4cout << "World material set to " << fWorldMaterial->GetName() << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstruction::SetTankMaterial(const G4String& mat)
{
  G4Material* pmat = G4NistManager::Instance()->FindOrBuildMaterial(mat);
  if(pmat && fTankMaterial != pmat)
  {
    fTankMaterial = pmat;
    if(fTank_LV)
    {
      fTank_LV->SetMaterial(fTankMaterial);
      fTankMaterial->SetMaterialPropertiesTable(fTankMPT);
      fTankMaterial->GetIonisation()->SetBirksConstant(0.126 * mm / MeV);
    }
    G4RunManager::GetRunManager()->PhysicsHasBeenModified();
    G4cout << "Tank material set to " << fTankMaterial->GetName() << G4endl;
  }
}
