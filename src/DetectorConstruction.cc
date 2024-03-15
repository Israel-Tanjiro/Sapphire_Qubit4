//This is the Class where I build the parts of the Sapphire QUBIT
//		G4CMPPhononElectrode to demonstrate KaplanQP.

#include "DetectorConstruction.hh"
//#include "PhononSensitivity.hh"
#include "G4CMPLogicalBorderSurface.hh"
#include "G4CMPPhononElectrode.hh"
#include "G4CMPSurfaceProperty.hh"
#include "G4Box.hh"
#include "G4Para.hh"
#include "G4Trd.hh"
#include "G4Colour.hh"
#include "G4GeometryManager.hh"
#include "G4LatticeLogical.hh"
#include "G4LatticeManager.hh"
#include "G4LatticePhysical.hh"
#include "G4LogicalVolume.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4Material.hh"
#include "G4MaterialPropertiesTable.hh"
#include "G4NistManager.hh"
#include "G4PVPlacement.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4RunManager.hh"
#include "G4SDManager.hh"
#include "G4SolidStore.hh"
#include "G4SystemOfUnits.hh"
#include "G4TransportationManager.hh"
#include "G4Tubs.hh"
#include "G4UserLimits.hh"
#include "G4VisAttributes.hh"
#include "CLHEP/Units/SystemOfUnits.h"
#include "G4CMPPhononElectrode.hh"
#include <cmath>
#include "G4UnionSolid.hh"
#include "G4VisAttributes.hh"
#include "G4SubtractionSolid.hh"
#include "G4Tet.hh"
#include "G4IntersectionSolid.hh"
#include "G4MultiUnion.hh"
////////---------Note I need to add the Nb Interface and Copper Interface, and the Lattice to the Sapphire
////////-------------Note the Big Trick is to define the Superconductor Ground Plane as Mother Volume

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

DetectorConstruction::DetectorConstruction()
  : fLiquidHelium(0), fSapphire(0),fAluminum(0),fNiobium (0),copper_mat(0), fTungsten(0),fOxigen(0),fAl2O3(0),
    fWorldPhys(0),fConstructed(false) {;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

DetectorConstruction::~DetectorConstruction() {
  // delete topSurfProp;
  // delete botSurfProp;
  // delete wallSurfProp;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4VPhysicalVolume* DetectorConstruction::Construct()
{
  if (fConstructed) {
    if (!G4RunManager::IfGeometryHasBeenDestroyed()) {
      // Run manager hasn't cleaned volume stores. This code shouldn't execute
      G4GeometryManager::GetInstance()->OpenGeometry();
      G4PhysicalVolumeStore::GetInstance()->Clean();
      G4LogicalVolumeStore::GetInstance()->Clean();
      G4SolidStore::GetInstance()->Clean();
    }
    // Have to completely remove all lattices to avoid warning on reconstruction
    // G4LatticeManager::GetLatticeManager()->Reset();
    // // Clear all LogicalSurfaces
    // // NOTE: No need to redefine the G4CMPSurfaceProperties
    // G4CMPLogicalBorderSurface::CleanSurfaceTable();
  }

  DefineMaterials();
  SetupGeometry();
  fConstructed = true;

  return fWorldPhys;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void DetectorConstruction::DefineMaterials()
{
  G4NistManager* nistManager = G4NistManager::Instance();

  fLiquidHelium = nistManager->FindOrBuildMaterial("G4_AIR"); // to be corrected
  //fGermanium = nistManager->FindOrBuildMaterial("G4_Ge");
  //fGermanium = nistManager->FindOrBuildMaterial("G4_Si");
  fAluminum = nistManager->FindOrBuildMaterial("G4_Al");
  fOxigen = nistManager->FindOrBuildMaterial("G4_O");
  fTungsten = nistManager->FindOrBuildMaterial("G4_W");
    fNiobium = nistManager->FindOrBuildMaterial("G4_Nb");
    copper_mat = nistManager->FindOrBuildMaterial("G4_Cu");
 fSapphire = new G4Material("fAl2O3", 3.98*g/cm3, 2);
    fSapphire->AddElement(nistManager->FindOrBuildElement("Al"), 2);
 fSapphire->AddElement(nistManager->FindOrBuildElement("O"), 3);


}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void DetectorConstruction::SetupGeometry()
{



G4VSolid* worldSolid = new G4Box("World",16.*cm,16.*cm,16.*cm);
G4LogicalVolume* worldLogical =new G4LogicalVolume(worldSolid,fLiquidHelium,"world");

// G4VSolid* solid_world = new G4Box("World",55.*cm,55.*cm,55.*cm);
// G4LogicalVolume* log_world = new G4LogicalVolume(solid_world,fLiquidHelium,"World");
// fWorldPhys = new G4PVPlacement(0,
//        G4ThreeVector(),
//        log_world,
//        "World",
//        0,
//                                false,
//        0);
constexpr double Position_Superconductor=430*CLHEP::um;
worldLogical->SetUserLimits(new G4UserLimits(10*mm, DBL_MAX, DBL_MAX, 0, 0));
fWorldPhys = new G4PVPlacement(0,G4ThreeVector(),worldLogical,"World",0,false,0);
// Set visual attributes
G4VisAttributes* WorldA = new G4VisAttributes(G4VisAttributes::Invisible); // Green
worldLogical->SetVisAttributes(WorldA);
//Setting the Substrate
G4VSolid* fSapphireSolid = new G4Box("fSapphireSolid",0.7*cm,0.7*cm,430*CLHEP::um);
G4LogicalVolume* fSapphireLogical =new G4LogicalVolume(fSapphireSolid,fSapphire,"fSapphireLogical");
G4VPhysicalVolume* SapphirePhys =  new G4PVPlacement(0,G4ThreeVector(0.0,0.0,0.0*cm),fSapphireLogical,"fSapphirePhysical",worldLogical,false,0);
G4VisAttributes* SuperColorSub = new G4VisAttributes(G4Colour(0.0,0.0,0.0)); // Green
fSapphireLogical->SetVisAttributes(SuperColorSub);
//Loading the Lattice
//Set up the G4CMP silicon lattice information using the G4LatticeManager
// G4LatticeManager gives physics processes access to lattices by volume
G4LatticeManager* LM = G4LatticeManager::GetLatticeManager();
G4LatticeLogical* log_sapphireLattice = LM->LoadLattice(fSapphire, "Al2O3");// Modifying for the Momento to Sapphire

// G4LatticePhysical assigns G4LatticeLogical a physical orientation
G4LatticePhysical* phys_sapphireLattice = new G4LatticePhysical(log_sapphireLattice);
phys_sapphireLattice->SetMillerOrientation(0,0,1); //No idea what this should be yet... Going to keep it unmotivated for now.
// Previous Value phys_siliconLattice->SetMillerOrientation(1,0,0); //No idea what this should be yet... Going to keep it unmotivated for now.

LM->RegisterLattice(SapphirePhys,phys_sapphireLattice);




//Setting the GroundPlane This will be the Mother Volume
G4VSolid* GroundPlane= new G4Box("GroundPlane",0.7*cm,0.7*cm,0.01*cm);
G4LogicalVolume* GroundPlaneLogical =new G4LogicalVolume(GroundPlane,fSapphire,"GroundPlaneLogical");
G4VPhysicalVolume* GroundPlanePhys =  new G4PVPlacement(0,G4ThreeVector(0.0,0.0,Position_Superconductor),GroundPlaneLogical,"GroundPlanePhysical",  worldLogical,false,0);
G4VisAttributes* SuperColor = new G4VisAttributes(G4Colour(0.0,1.0,1.0,0.5)); // Green
G4VisAttributes* SuperColor1 = new G4VisAttributes(G4Colour(0.0,0.0,1.0,0.5)); // Green
GroundPlaneLogical->SetVisAttributes(SuperColor);
/////////Defining the SupeSolidAdditions


constexpr double Vacu_GroundPlane= 10.0* CLHEP::um;
constexpr double Z_Draw= 1000.0* CLHEP::um;///Variable to draw the componente separate from the ground plane

/////--------------------------------------------
/////////////////////////////////////// List of All logical Volumes on the Detector.
//////--------------------
const G4double GHz = 1e9 * hertz;

//  the following coefficients and cutoff values are not well-motivated
//  the code below is used only to demonstrate how to set these values.
// const std::vector<G4double> anhCoeffs = {0,0,0,0,0,0};//{0, 0, 0, 0, 0, 1.51e-14}; //Turn this off temporarily
// const std::vector<G4double> diffCoeffs = {1,0,0,0,0,0};//{5.88e-2, 7.83e-4, -2.47e-6, 1.71e-8, -2.98e-11}; //Explicitly make this 1 so that there's no mistake, because the way these are accessed is confusing.
// const std::vector<G4double> specCoeffs = {0,0,0,0,0,0};//{0,928, -2.03e-4, -3.21e-6, 3.1e-9, 2.9e-13}; //Turn this off temporarily

//////////////////Mikes Valesu

const std::vector<G4double> anhCoeffs = {0, 0, 0, 0, 0, 1.51e-14};
const std::vector<G4double> diffCoeffs =
  {5.88e-2, 7.83e-4, -2.47e-6, 1.71e-8, -2.98e-11};
const std::vector<G4double> specCoeffs =
  {0,928, -2.03e-4, -3.21e-6, 3.1e-9, 2.9e-13};

 const G4double anhCutoff = 520., reflCutoff = 350.;   // Units external


//These are currently not motivated at all. Need to understand what the parameters are and how to define for these various surfaces.
//Only the surfaces facing the silicon are needed, I think, since G4CMP doesn't know how to propagate things in supeconductors.
//These are just the definitions of the interface TYPES, not the interfaces themselves. These must be called in a set of loops
//below, and invoke these surface definitions.

  //filmAbsorptiones=0.3445;
  //G4cout<<"**************This to test the parameter of the filmAbsorption Nb ************************ \t"<<filmAbsorptiones<<G4endl;

   G4CMPSurfaceProperty *fSiNbInterface = new G4CMPSurfaceProperty("SiNbInterface",
              1., 0., 0., 0.,
              1.0, 1.0, 0., 0.);

  G4CMPSurfaceProperty *fSiCopperInterface = new G4CMPSurfaceProperty("SiCopperInterface",
            1.0, 0.0, 0.0, 0.0,
            0.0, 1.0, 0.0, 0.0 );
  G4CMPSurfaceProperty *fSiVacuumInterface = new G4CMPSurfaceProperty("SiVacuumInterface",
            1.0, 0.0, 0.0, 0.0,
            0.0, 1.0, 0.0, 0.0 );

  fSiNbInterface->AddScatteringProperties(anhCutoff, reflCutoff, anhCoeffs,
            diffCoeffs, specCoeffs, GHz, GHz, GHz);
  fSiCopperInterface->AddScatteringProperties(anhCutoff, reflCutoff, anhCoeffs,
          diffCoeffs, specCoeffs, GHz, GHz, GHz);
  fSiVacuumInterface->AddScatteringProperties(anhCutoff, reflCutoff, anhCoeffs,
          diffCoeffs, specCoeffs, GHz, GHz, GHz);









//Building the properties of the Resonator
///////////-----------------------
// 7 Rectangles and 8 Half Circles This is the General Class to include all resonator on the Superconductor Part.
// MARCH 2024 08 I ADD A NEW VARIABLE TO CREATE THE GROUND PLANE HOLES FOR THIS I NEED TO REDEFINE ALLL THE SOLID BY USING GROUND PLANE AND ADDING MORE WITH
constexpr double dp_resonatorAssemblyBaseNbDimX_1 = 644.0 * CLHEP::um;
constexpr double dp_resonatorAssemblyBaseNbDimY_1 = 10.0* CLHEP::um;
constexpr double dp_resonatorAssemblyBaseNbDimZ_1 = 0.01*cm;////Thicnkess Superconductor


G4VSolid* Resonator_Square_1 = new G4Box("Resonator_Square_1",dp_resonatorAssemblyBaseNbDimX_1,dp_resonatorAssemblyBaseNbDimY_1,dp_resonatorAssemblyBaseNbDimZ_1);
G4LogicalVolume* Resonator_SquareLogical_1 = new G4LogicalVolume( Resonator_Square_1,fNiobium,"Resonator_SquareLogical_1");
G4VSolid* Resonator_Square_1_GP = new G4Box("Resonator_Square_1_GP",dp_resonatorAssemblyBaseNbDimX_1,dp_resonatorAssemblyBaseNbDimY_1+Vacu_GroundPlane,dp_resonatorAssemblyBaseNbDimZ_1);
G4LogicalVolume* Resonator_SquareLogical_1_GP = new G4LogicalVolume( Resonator_Square_1_GP,fLiquidHelium,"Resonator_SquareLogical_1_GP");
G4double tube_dPhi = 0.5*M_PI*rad;
Resonator_SquareLogical_1->SetVisAttributes(SuperColor);
Resonator_SquareLogical_1_GP->SetVisAttributes(WorldA);
/////
G4VSolid* Half_Cirlce_R_1= new G4Tubs("Half_Cirlce_R_1",  30.0* CLHEP::um, 50.0* CLHEP::um,dp_resonatorAssemblyBaseNbDimZ_1,-1.0*tube_dPhi,2*tube_dPhi ); // segment angle
G4LogicalVolume* Half_Cirlce_RLogical_1 = new G4LogicalVolume( Half_Cirlce_R_1,fNiobium,"Half_Cirlce_RLogical_1");
G4VSolid* Half_Cirlce_R_1_GP= new G4Tubs("Half_Cirlce_R_1_GP",  30.0* CLHEP::um-Vacu_GroundPlane, 50.0* CLHEP::um+Vacu_GroundPlane,dp_resonatorAssemblyBaseNbDimZ_1,-1.0*tube_dPhi,2*tube_dPhi ); // segment angle
G4LogicalVolume* Half_Cirlce_RLogical_1_GP = new G4LogicalVolume( Half_Cirlce_R_1_GP,fLiquidHelium,"Half_Cirlce_RLogical_1_GP");

G4VSolid* Half_Cirlce_L_1= new G4Tubs("Half_Cirlce_L_1",  30.0* CLHEP::um, 50.0* CLHEP::um,dp_resonatorAssemblyBaseNbDimZ_1,tube_dPhi,2*tube_dPhi ); // segment angle
G4LogicalVolume* Half_Cirlce_LLogical_1 = new G4LogicalVolume( Half_Cirlce_L_1,fNiobium,"Half_Cirlce_LLogical_1");
G4VSolid* Half_Cirlce_L_1_GP= new G4Tubs("Half_Cirlce_L_1_GP",  30.0* CLHEP::um-Vacu_GroundPlane, 50.0* CLHEP::um+Vacu_GroundPlane,dp_resonatorAssemblyBaseNbDimZ_1,tube_dPhi,2*tube_dPhi ); // segment angle
G4LogicalVolume* Half_Cirlce_LLogical_1_GP = new G4LogicalVolume( Half_Cirlce_L_1_GP,fLiquidHelium,"Half_Cirlce_LLogical_1_GP");//////////--------------VARiables for ground Plane  RESONATOR 1-----------------------








//////////--------------VARiables for ground Plane  RESONATOR 1-----------------------






// for(G4int i=0; i<8;i++){
//   if(i<7){
//   G4VPhysicalVolume* Resonator_SquarePhys_1 =new G4PVPlacement(0,G4ThreeVector(X_Center_1,i*80.0*CLHEP::um+Y_Center_1,0.37*cm),Resonator_SquareLogical_1,"Resonator_Square_1",
//                       worldLogical,false,i,true);
// }
//
//                       if (pow(-1,i)>0) {
//                         G4VPhysicalVolume* Half_Cirlce_RPhys_1 =new G4PVPlacement(0,G4ThreeVector(X_Center_1+pow(-1,i)*dp_resonatorAssemblyBaseNbDimX_1,40.0*CLHEP::um+i*80.0*CLHEP::um+Y_Center_1,0.37*cm),Half_Cirlce_RLogical_1,"Half_Cirlce_R_1",
//                                                                 worldLogical,false,i,true);
//                       }
//                       else{G4VPhysicalVolume* Half_Cirlce_LPhys_1 =new G4PVPlacement(0,G4ThreeVector(X_Center_1+pow(-1,i)*dp_resonatorAssemblyBaseNbDimX_1,-120.0*CLHEP::um+i*80.0*CLHEP::um+Y_Center_1,0.37*cm),Half_Cirlce_LLogical_1,"Half_Cirlce_L_1",
//                                                               worldLogical,false,i,true);}
//
//
// }
//-------------------------------~~~~~~~~~~~~~~~~~~~``

//Creating the Cross_Qubit
constexpr double dp_CrossQubitNbDimX = 320.0 * CLHEP::um;
constexpr double dp_CrossQubitNbDimY = 25.0* CLHEP::um;
constexpr double dp_CrossQubitNbDimZ = dp_resonatorAssemblyBaseNbDimZ_1;

G4VSolid* Cross_1 = new G4Box("Cross_1",dp_CrossQubitNbDimX,dp_CrossQubitNbDimY,dp_CrossQubitNbDimZ);
G4VSolid* Cross_2 = new G4Box("Cross_2",dp_CrossQubitNbDimY,dp_CrossQubitNbDimX,dp_CrossQubitNbDimZ);
G4VSolid* Cross_Qubit= new G4UnionSolid("Cross_Qubit", Cross_1, Cross_2);

G4VSolid* Cross_1_GP = new G4Box("Cross_1_GP",dp_CrossQubitNbDimX+3*Vacu_GroundPlane,dp_CrossQubitNbDimY+2*Vacu_GroundPlane,dp_CrossQubitNbDimZ);
G4VSolid* Cross_2_GP = new G4Box("Cross_2_GP",dp_CrossQubitNbDimY+2*Vacu_GroundPlane,dp_CrossQubitNbDimX+3*Vacu_GroundPlane,dp_CrossQubitNbDimZ);
G4VSolid* Cross_Qubit_GP= new G4UnionSolid("Cross_Qubit_GP", Cross_1_GP, Cross_2_GP);

G4LogicalVolume* Cross_QubitLogical = new G4LogicalVolume( Cross_Qubit,fNiobium,"Cross_QubitLogical");
G4LogicalVolume* Cross_QubitLogical_GP = new G4LogicalVolume( Cross_Qubit_GP,fLiquidHelium,"Cross_QubitLogical_GP");
////// DOing the Part where I Fill with Vacuum





//The Cross Qubits is Working well------------------
/////////////-------------------- I need 6 Panels





//Note Remember is the half lenght for the x and y and z coordiantes of every box

//February 29 2024 ---------_Creating the Second Resonators  -- Checking the measurments
//////////////---------------------RESONATOR 2--------------------
constexpr double dp_resonatorAssemblyBaseNbDimX_2 = 694.0 * CLHEP::um;
constexpr double dp_resonatorAssemblyBaseNbDimY_2 = 10.0* CLHEP::um;
constexpr double dp_resonatorAssemblyBaseNbDimZ_2 = 0.01*cm;


G4VSolid* Resonator_Square_2 = new G4Box("Resonator_Square_2",dp_resonatorAssemblyBaseNbDimX_2,dp_resonatorAssemblyBaseNbDimY_2,dp_resonatorAssemblyBaseNbDimZ_2);
G4LogicalVolume* Resonator_SquareLogical_2 = new G4LogicalVolume( Resonator_Square_2,fNiobium,"Resonator_SquareLogical_2");
G4VSolid* Resonator_Square_2_GP = new G4Box("Resonator_Square_2_GP",dp_resonatorAssemblyBaseNbDimX_2,dp_resonatorAssemblyBaseNbDimY_2+Vacu_GroundPlane,dp_resonatorAssemblyBaseNbDimZ_2);
G4LogicalVolume* Resonator_SquareLogical_2_GP = new G4LogicalVolume( Resonator_Square_2_GP,fLiquidHelium,"Resonator_SquareLogical_2_GP");
G4double tube_dPhi_2 = 0.5*M_PI*rad;
/////
G4VSolid* Half_Cirlce_R_2= new G4Tubs("Half_Cirlce_R_2",  30.0* CLHEP::um, 50.0* CLHEP::um,10* CLHEP::um,-1.0*tube_dPhi_2,2*tube_dPhi_2 ); // segment angle
G4LogicalVolume* Half_Cirlce_RLogical_2 = new G4LogicalVolume( Half_Cirlce_R_2,fNiobium,"Half_Cirlce_RLogical_2_GP");
G4VSolid* Half_Cirlce_R_2_GP= new G4Tubs("Half_Cirlce_R_2_GP",  30.0* CLHEP::um-Vacu_GroundPlane, 50.0* CLHEP::um+Vacu_GroundPlane,10* CLHEP::um,-1.0*tube_dPhi_2,2*tube_dPhi_2 ); // segment angle
G4LogicalVolume* Half_Cirlce_RLogical_2_GP = new G4LogicalVolume( Half_Cirlce_R_2_GP,fLiquidHelium,"Half_Cirlce_RLogical_2_GP");

G4VSolid* Half_Cirlce_L_2= new G4Tubs("Half_Cirlce_L_2",  30.0* CLHEP::um, 50.0* CLHEP::um,dp_resonatorAssemblyBaseNbDimZ_1,tube_dPhi_2,2*tube_dPhi_2 ); // segment angle
G4LogicalVolume* Half_Cirlce_LLogical_2 = new G4LogicalVolume( Half_Cirlce_L_2,fNiobium,"Half_Cirlce_LLogical_2");
G4VSolid* Half_Cirlce_L_2_GP= new G4Tubs("Half_Cirlce_L_2_GP",  30.0* CLHEP::um-Vacu_GroundPlane, 50.0* CLHEP::um+Vacu_GroundPlane,dp_resonatorAssemblyBaseNbDimZ_1,tube_dPhi_2,2*tube_dPhi_2 ); // segment angle
G4LogicalVolume* Half_Cirlce_LLogical_2_GP = new G4LogicalVolume( Half_Cirlce_L_2_GP,fLiquidHelium,"Half_Cirlce_LLogical_2_GP");
// for(G4int i=0; i<8;i++){
//   if(i<7){
//   G4VPhysicalVolume* Resonator_SquarePhys_2 =new G4PVPlacement(0,G4ThreeVector(X_Center_2,i*80.0*CLHEP::um+Y_Center_2,Position_Superconductor),Resonator_SquareLogical_2,"Resonator_Square_2",
//                       worldLogical,false,i,true);
// }
//
//                       if (pow(-1,i)>0) {
//                         G4VPhysicalVolume* Half_Cirlce_RPhys_2 =new G4PVPlacement(0,G4ThreeVector(X_Center_2+pow(-1,i)*dp_resonatorAssemblyBaseNbDimX_2,40.0*CLHEP::um+i*80.0*CLHEP::um+Y_Center_2,Position_Superconductor),Half_Cirlce_RLogical_2,"Half_Cirlce_R_2",
//                                                                 worldLogical,false,i,true);
//                       }
//                       else{G4VPhysicalVolume* Half_Cirlce_LPhys_2 =new G4PVPlacement(0,G4ThreeVector(X_Center_2+pow(-1,i)*dp_resonatorAssemblyBaseNbDimX_2,-120.0*CLHEP::um+i*80.0*CLHEP::um+Y_Center_2,Position_Superconductor),Half_Cirlce_LLogical_2,"Half_Cirlce_L_2",
//                                                               worldLogical,false,i,true);}
//
//
// }

/////////////-------------------------RESONATOR 3---------------------------
constexpr double dp_resonatorAssemblyBaseNbDimX_3 = 743.0 * CLHEP::um;
constexpr double dp_resonatorAssemblyBaseNbDimY_3 = 10.0* CLHEP::um;
constexpr double dp_resonatorAssemblyBaseNbDimZ_3 = 0.01*cm;


G4VSolid* Resonator_Square_3 = new G4Box("Resonator_Square_3",dp_resonatorAssemblyBaseNbDimX_3,dp_resonatorAssemblyBaseNbDimY_3,dp_resonatorAssemblyBaseNbDimZ_3);
G4LogicalVolume* Resonator_SquareLogical_3 = new G4LogicalVolume( Resonator_Square_3,fNiobium,"Resonator_SquareLogical_3");
G4VSolid* Resonator_Square_3_GP = new G4Box("Resonator_Square_3_GP",dp_resonatorAssemblyBaseNbDimX_3,dp_resonatorAssemblyBaseNbDimY_3+Vacu_GroundPlane,dp_resonatorAssemblyBaseNbDimZ_3);
G4LogicalVolume* Resonator_SquareLogical_3_GP = new G4LogicalVolume( Resonator_Square_3_GP,fLiquidHelium,"Resonator_SquareLogical_3_GP");
G4double tube_dPhi_3 = 0.5*M_PI*rad;
/////
G4VSolid* Half_Cirlce_R_3= new G4Tubs("Half_Cirlce_R_3",  30.0* CLHEP::um, 50.0* CLHEP::um,dp_resonatorAssemblyBaseNbDimZ_1,-1.0*tube_dPhi_3,2*tube_dPhi_3 ); // segment angle
G4LogicalVolume* Half_Cirlce_RLogical_3 = new G4LogicalVolume( Half_Cirlce_R_3,fNiobium,"Half_Cirlce_RLogical_3");
G4VSolid* Half_Cirlce_R_3_GP= new G4Tubs("Half_Cirlce_R_3_GP",  30.0* CLHEP::um-Vacu_GroundPlane, 50.0* CLHEP::um+Vacu_GroundPlane,dp_resonatorAssemblyBaseNbDimZ_1,-1.0*tube_dPhi_3,2*tube_dPhi_3 ); // segment angle
G4LogicalVolume* Half_Cirlce_RLogical_3_GP = new G4LogicalVolume( Half_Cirlce_R_3_GP,fLiquidHelium,"Half_Cirlce_RLogical_3_GP");

G4VSolid* Half_Cirlce_L_3= new G4Tubs("Half_Cirlce_L_3",  30.0* CLHEP::um, 50.0* CLHEP::um,dp_resonatorAssemblyBaseNbDimZ_1,tube_dPhi_3,2*tube_dPhi_3 ); // segment angle
G4LogicalVolume* Half_Cirlce_LLogical_3 = new G4LogicalVolume( Half_Cirlce_L_3,fNiobium,"Half_Cirlce_LLogical_3");
G4VSolid* Half_Cirlce_L_3_GP= new G4Tubs("Half_Cirlce_L_3_GP",  30.0* CLHEP::um-Vacu_GroundPlane, 50.0* CLHEP::um+Vacu_GroundPlane,dp_resonatorAssemblyBaseNbDimZ_1,tube_dPhi_3,2*tube_dPhi_3 ); // segment angle
G4LogicalVolume* Half_Cirlce_LLogical_3_GP = new G4LogicalVolume( Half_Cirlce_L_3_GP,fLiquidHelium,"Half_Cirlce_LLogical_3_GP");

// for(G4int i=0; i<8;i++){
//   if(i<7){
//   G4VPhysicalVolume* Resonator_SquarePhys_3 =new G4PVPlacement(0,G4ThreeVector(X_Center_3,i*80.0*CLHEP::um+Y_Center_3,Position_Superconductor),Resonator_SquareLogical_3,"Resonator_Square_3",
//                       worldLogical,false,i,true);
// }
//
//                       if (pow(-1,i)>0) {
//                         G4VPhysicalVolume* Half_Cirlce_RPhys_3 =new G4PVPlacement(0,G4ThreeVector(X_Center_3+pow(-1,i)*dp_resonatorAssemblyBaseNbDimX_3,40.0*CLHEP::um+i*80.0*CLHEP::um+Y_Center_3,Position_Superconductor),Half_Cirlce_RLogical_3,"Half_Cirlce_R_3",
//                                                                 worldLogical,false,i,true);
//                       }
//                       else{G4VPhysicalVolume* Half_Cirlce_LPhys_3 =new G4PVPlacement(0,G4ThreeVector(X_Center_3+pow(-1,i)*dp_resonatorAssemblyBaseNbDimX_3,-120.0*CLHEP::um+i*80.0*CLHEP::um+Y_Center_3,Position_Superconductor),Half_Cirlce_LLogical_3,"Half_Cirlce_L_3",
//                                                               worldLogical,false,i,true);}
//
//
// }
/////////////////////////-------------------RESONATOR4-----------------/////////////////
constexpr double dp_resonatorAssemblyBaseNbDimX_4 = 743.0 * CLHEP::um;
constexpr double dp_resonatorAssemblyBaseNbDimY_4 = 10.0* CLHEP::um;
//constexpr double dp_resonatorAssemblyBaseNbDimZ_1 = 10*CLHEP::um;


G4VSolid* Resonator_Square_4 = new G4Box("Resonator_Square_4",dp_resonatorAssemblyBaseNbDimX_4,dp_resonatorAssemblyBaseNbDimY_4,dp_resonatorAssemblyBaseNbDimZ_1);
G4LogicalVolume* Resonator_SquareLogical_4 = new G4LogicalVolume( Resonator_Square_4,fNiobium,"Resonator_SquareLogical_4");
G4VSolid* Resonator_Square_4_GP = new G4Box("Resonator_Square_4_GP",dp_resonatorAssemblyBaseNbDimX_4,dp_resonatorAssemblyBaseNbDimY_4+Vacu_GroundPlane,dp_resonatorAssemblyBaseNbDimZ_1);
G4LogicalVolume* Resonator_SquareLogical_4_GP = new G4LogicalVolume( Resonator_Square_4_GP,fLiquidHelium,"Resonator_SquareLogical_4_GP");
G4double tube_dPhi_4 = 0.5*M_PI*rad;
/////
G4VSolid* Half_Cirlce_R_4= new G4Tubs("Half_Cirlce_R_4",  30.0* CLHEP::um, 50.0* CLHEP::um,dp_resonatorAssemblyBaseNbDimZ_1,-1.0*tube_dPhi_4,2*tube_dPhi_4 ); // segment angle
G4VSolid* Half_Cirlce_R_4_GP= new G4Tubs("Half_Cirlce_R_4_GP",  30.0* CLHEP::um-Vacu_GroundPlane, 50.0* CLHEP::um+Vacu_GroundPlane,dp_resonatorAssemblyBaseNbDimZ_1,-1.0*tube_dPhi_4,2*tube_dPhi_4 ); // segment angle
G4LogicalVolume* Half_Cirlce_RLogical_4 = new G4LogicalVolume( Half_Cirlce_R_4,fNiobium,"Half_Cirlce_RLogical_3");
G4LogicalVolume* Half_Cirlce_RLogical_4_GP = new G4LogicalVolume( Half_Cirlce_R_4_GP,fLiquidHelium,"Half_Cirlce_RLogical_3_GP");

G4VSolid* Half_Cirlce_L_4= new G4Tubs("Half_Cirlce_L_4",  30.0* CLHEP::um, 50.0* CLHEP::um,dp_resonatorAssemblyBaseNbDimZ_1,tube_dPhi_4,2*tube_dPhi_4 ); // segment angle
G4LogicalVolume* Half_Cirlce_LLogical_4 = new G4LogicalVolume( Half_Cirlce_L_4,fNiobium,"Half_Cirlce_LLogical_4");
G4VSolid* Half_Cirlce_L_4_GP= new G4Tubs("Half_Cirlce_L_4_GP",  30.0* CLHEP::um-Vacu_GroundPlane, 50.0* CLHEP::um+Vacu_GroundPlane,dp_resonatorAssemblyBaseNbDimZ_1,tube_dPhi_4,2*tube_dPhi_4 ); // segment angle
G4LogicalVolume* Half_Cirlce_LLogical_4_GP = new G4LogicalVolume( Half_Cirlce_L_4_GP,fLiquidHelium,"Half_Cirlce_LLogical_4_GP");
// for(G4int i=0; i<8;i++){
//   if(i<7){
//   G4VPhysicalVolume* Resonator_SquarePhys_4 =new G4PVPlacement(0,G4ThreeVector(X_Center_4,i*80.0*CLHEP::um+Y_Center_4,Position_Superconductor),Resonator_SquareLogical_4,"Resonator_Square_4",
//                       worldLogical,false,i,true);
// }
//
//                       if (pow(-1,i)>0 ) {
//                         G4VPhysicalVolume* Half_Cirlce_RPhys_4 =new G4PVPlacement(0,G4ThreeVector(X_Center_4+pow(-1,i)*dp_resonatorAssemblyBaseNbDimX_4,40.0*CLHEP::um+i*80.0*CLHEP::um+Y_Center_4-80.0*CLHEP::um,Position_Superconductor),Half_Cirlce_RLogical_4,"Half_Cirlce_R_4",
//                                                                 worldLogical,false,i,true);
//                       }
//                       else{G4VPhysicalVolume* Half_Cirlce_LPhys_4 =new G4PVPlacement(0,G4ThreeVector(X_Center_4+pow(-1,i)*dp_resonatorAssemblyBaseNbDimX_4,-120.0*CLHEP::um+i*80.0*CLHEP::um+Y_Center_4+80.0*CLHEP::um,Position_Superconductor),Half_Cirlce_LLogical_4,"Half_Cirlce_L_4",
//                                                               worldLogical,false,i,true);}
//
//
// }


/////////////////////////////-------------------------_Creating the Top Transsmission line ----------------------
////We have 2 half circles, 2 quarter circles, one big rectangle and two small Rectangles

////////////////////
//Defining the Circles
constexpr double dp_TransmissionLineBaseNbDimX = 4350.0 * CLHEP::um;
constexpr double dp_TransmissionLineBaseNbDimY = 20.0* CLHEP::um;
constexpr double dp_TransmissionLineBaseNbDimZ = 10* CLHEP::um;
constexpr double X_Center_TransmissionLine = 0* CLHEP::um;
constexpr double Y_Center_TransmissionLine = 3905* CLHEP::um;///Parameter to Move the tranmission line

constexpr double dp_TransmissionLineBaseNbDimX_1 = 150.0 * CLHEP::um;
constexpr double dp_TransmissionLineBaseNbDimY_1 = 20.0* CLHEP::um;
constexpr double dp_TransmissionLineBaseNbDimZ_1 = 10* CLHEP::um;
////Staritng Here
G4VSolid* TransmissionLine_Resonator_Square = new G4Box("TransmissionLine_Resonator_Square",dp_TransmissionLineBaseNbDimX,dp_TransmissionLineBaseNbDimY,dp_TransmissionLineBaseNbDimZ);
G4LogicalVolume* TransmissionLine_Resonator_SquareLogical = new G4LogicalVolume( TransmissionLine_Resonator_Square,fNiobium,"TransmissionLine_Resonator_SquareLogical");
G4VSolid* TransmissionLine_Resonator_Square_GP = new G4Box("TransmissionLine_Resonator_Square_GP",dp_TransmissionLineBaseNbDimX,dp_TransmissionLineBaseNbDimY+2*Vacu_GroundPlane,dp_TransmissionLineBaseNbDimZ);
G4LogicalVolume* TransmissionLine_Resonator_Square_GPLogical = new G4LogicalVolume( TransmissionLine_Resonator_Square_GP,fLiquidHelium,"TransmissionLine_Resonator_SquareGPLogical");







G4VSolid* Half_Cirlce_R_Transmission= new G4Tubs("Half_Cirlce_R__Transmission",  70.0* CLHEP::um, 110.0* CLHEP::um,dp_resonatorAssemblyBaseNbDimZ_1,-1.0*tube_dPhi_4,2*tube_dPhi_4 ); // segment angle
G4LogicalVolume* Half_Cirlce_RLogical__Transmission = new G4LogicalVolume( Half_Cirlce_R_Transmission,fNiobium,"Half_Cirlce_RLogical__Transmission");
G4VSolid* Half_Cirlce_R_Transmission_GP= new G4Tubs("Half_Cirlce_R__Transmission_GP",  70.0* CLHEP::um-2*Vacu_GroundPlane, 110.0* CLHEP::um+2*Vacu_GroundPlane,dp_resonatorAssemblyBaseNbDimZ_1,-1.0*tube_dPhi_4,2*tube_dPhi_4 ); // segment angle
G4LogicalVolume* Half_Cirlce_RLogical__Transmission_GP = new G4LogicalVolume( Half_Cirlce_R_Transmission_GP,fLiquidHelium,"Half_Cirlce_RLogical__Transmission_GP");


G4VSolid* Half_Cirlce_L__Transmission= new G4Tubs("Half_Cirlce_L__Transmission",  70.0* CLHEP::um, 110.0* CLHEP::um,dp_resonatorAssemblyBaseNbDimZ_1,tube_dPhi_4,2*tube_dPhi_4 ); // segment angle
G4LogicalVolume* Half_Cirlce_LLogical__Transmission = new G4LogicalVolume( Half_Cirlce_L__Transmission,fNiobium,"Half_Cirlce_LLogical__Transmission");
G4VSolid* Half_Cirlce_L__Transmission_GP= new G4Tubs("Half_Cirlce_L__Transmission_GP",  70.0* CLHEP::um-2*Vacu_GroundPlane, 110.0* CLHEP::um+2*Vacu_GroundPlane,dp_resonatorAssemblyBaseNbDimZ_1,tube_dPhi_4,2*tube_dPhi_4 ); // segment angle
G4LogicalVolume* Half_Cirlce_LLogical__Transmission_GP = new G4LogicalVolume( Half_Cirlce_L__Transmission_GP,fLiquidHelium,"Half_Cirlce_LLogical__Transmission_GP");


G4VSolid* Quarter_Cirlce_R_Transmission= new G4Tubs("Quarter_Cirlce_R__Transmission",  70.0* CLHEP::um, 110.0* CLHEP::um,dp_resonatorAssemblyBaseNbDimZ_1,-1.0*tube_dPhi_4,1*tube_dPhi_4 ); // segment angle
G4LogicalVolume* Quarter_Cirlce_RLogical__Transmission = new G4LogicalVolume(Quarter_Cirlce_R_Transmission,fNiobium,"Quarter_Cirlce_RLogical__Transmission");
G4VSolid* Quarter_Cirlce_R_Transmission_GP= new G4Tubs("Quarter_Cirlce_R__Transmission_GP",  70.0* CLHEP::um-2*Vacu_GroundPlane, 110.0* CLHEP::um+2*Vacu_GroundPlane,dp_resonatorAssemblyBaseNbDimZ_1,-1.0*tube_dPhi_4,1*tube_dPhi_4 ); // segment angle
G4LogicalVolume* Quarter_Cirlce_RLogical__Transmission_GP = new G4LogicalVolume(Quarter_Cirlce_R_Transmission_GP,fLiquidHelium,"Quarter_Cirlce_RLogical__Transmission_GP");


G4VSolid* Quarter_Cirlce_L__Transmission= new G4Tubs("Quarter_Cirlce_L__Transmission",  70.0* CLHEP::um, 110.0* CLHEP::um,dp_resonatorAssemblyBaseNbDimZ_1,-2*tube_dPhi_4,1*tube_dPhi_4 ); // segment angle
G4LogicalVolume* Quarter_Cirlce_LLogical__Transmission = new G4LogicalVolume(Quarter_Cirlce_L__Transmission,fNiobium,"Quarter_Cirlce_LLogical__Transmission");
G4VSolid* Quarter_Cirlce_L__Transmission_GP= new G4Tubs("Quarter_Cirlce_L__Transmission_GP",  70.0* CLHEP::um-2*Vacu_GroundPlane, 110.0* CLHEP::um+2*Vacu_GroundPlane,dp_resonatorAssemblyBaseNbDimZ_1,-2*tube_dPhi_4,1*tube_dPhi_4 ); // segment angle
G4LogicalVolume* Quarter_Cirlce_LLogical__Transmission_GP = new G4LogicalVolume(Quarter_Cirlce_L__Transmission_GP,fLiquidHelium,"Quarter_Cirlce_LLogical__Transmission_GP");



G4VSolid* TransmissionLine_Resonator_Square_1 = new G4Box("TransmissionLine_Resonator_Square_1",dp_TransmissionLineBaseNbDimX_1,dp_TransmissionLineBaseNbDimY_1,dp_TransmissionLineBaseNbDimZ_1);
G4LogicalVolume* TransmissionLine_Resonator_SquareLogical_1 = new G4LogicalVolume( TransmissionLine_Resonator_Square_1,fNiobium,"TransmissionLine_Resonator_SquareLogical_1");
G4VSolid* TransmissionLine_Resonator_Square_1_GP = new G4Box("TransmissionLine_Resonator_Square_1_GP",dp_TransmissionLineBaseNbDimX_1,dp_TransmissionLineBaseNbDimY_1+2*Vacu_GroundPlane,dp_TransmissionLineBaseNbDimZ_1);
G4LogicalVolume* TransmissionLine_Resonator_SquareLogical_1_GP = new G4LogicalVolume(TransmissionLine_Resonator_Square_1_GP,fLiquidHelium,"TransmissionLine_Resonator_SquareLogical_1_GP");

G4VPhysicalVolume* Resonator_TransmissionLine_GP=new G4PVPlacement(0,G4ThreeVector(X_Center_TransmissionLine,Y_Center_TransmissionLine,0),TransmissionLine_Resonator_Square_GPLogical,"Resonator_TransmissionLine_GP",
                                        GroundPlaneLogical,false,0,true);

// G4VPhysicalVolume* Resonator_TransmissionLine=new G4PVPlacement(0,G4ThreeVector(X_Center_TransmissionLine,Y_Center_TransmissionLine,0.0*cm),TransmissionLine_Resonator_SquareLogical,"Resonator_TransmissionLine",
//                     TransmissionLine_Resonator_Square_GPLogical,false,0,true);
G4VPhysicalVolume* Resonator_TransmissionLine=new G4PVPlacement(0,G4ThreeVector(0.0*cm,0.0*cm,0.0*cm),TransmissionLine_Resonator_SquareLogical,"Resonator_TransmissionLine",
                                        TransmissionLine_Resonator_Square_GPLogical,false,0,true);


G4VPhysicalVolume* Half_Cirlce_RPhys_Transmission_GP =new G4PVPlacement(0,G4ThreeVector(dp_TransmissionLineBaseNbDimX,Y_Center_TransmissionLine+90.0* CLHEP::um,0),Half_Cirlce_RLogical__Transmission_GP,"Half_Cirlce_RPhys_Transmission_GP",GroundPlaneLogical,false,0,true);
G4VPhysicalVolume* Half_Cirlce_RPhys_Transmission =new G4PVPlacement(0,G4ThreeVector(),Half_Cirlce_RLogical__Transmission,"Half_Cirlce_RPhys_Transmission",Half_Cirlce_RLogical__Transmission_GP,false,0,true);




G4VPhysicalVolume* Half_Cirlce_LPhys_Transmission_GP =new G4PVPlacement(0,G4ThreeVector(-1.0*dp_TransmissionLineBaseNbDimX,Y_Center_TransmissionLine+90.0* CLHEP::um,0),Half_Cirlce_LLogical__Transmission_GP,"Half_Cirlce_LPhys_Transmission_GP",
                                                                                                    GroundPlaneLogical,false,0,true);
G4VPhysicalVolume* Half_Cirlce_LPhys_Transmission =new G4PVPlacement(0,G4ThreeVector(),Half_Cirlce_LLogical__Transmission,"Half_Cirlce_LPhys_Transmission",
                                                                                                  Half_Cirlce_LLogical__Transmission_GP,false,0,true);




G4VPhysicalVolume* Quarter_Cirlce_RPhys_Transmission_GP =new G4PVPlacement(0,G4ThreeVector(-1.0*dp_TransmissionLineBaseNbDimX+2*dp_TransmissionLineBaseNbDimX_1,Y_Center_TransmissionLine+3*90.0* CLHEP::um,0),Quarter_Cirlce_RLogical__Transmission_GP,"Quarter_Cirlce_RPhys_Transmission_GP",
                                                                                                    GroundPlaneLogical,false,0,true);
G4VPhysicalVolume* Quarter_Cirlce_RPhys_Transmission =new G4PVPlacement(0,G4ThreeVector(),Quarter_Cirlce_RLogical__Transmission,"Quarter_Cirlce_RPhys_Transmission",
                                                                                                      Quarter_Cirlce_RLogical__Transmission_GP,false,0,true);


G4VPhysicalVolume* Quarter_Cirlce_LPhys_Transmission_GP =new G4PVPlacement(0,G4ThreeVector(1.0*dp_TransmissionLineBaseNbDimX-2*dp_TransmissionLineBaseNbDimX_1,Y_Center_TransmissionLine+3*90.0* CLHEP::um,0),Quarter_Cirlce_LLogical__Transmission_GP,"Quarter_Cirlce_LPhys_Transmission_GP",GroundPlaneLogical,false,0,true);
G4VPhysicalVolume* Quarter_Cirlce_LPhys_Transmission =new G4PVPlacement(0,G4ThreeVector(),Quarter_Cirlce_LLogical__Transmission,"Quarter_Cirlce_LPhys_Transmission",
                                                                                                        Quarter_Cirlce_LLogical__Transmission_GP,false,0,true);


G4VPhysicalVolume* Resonator_TransmissionLine_1_GP=new G4PVPlacement(0,G4ThreeVector(-1.0*dp_TransmissionLineBaseNbDimX+1*dp_TransmissionLineBaseNbDimX_1,Y_Center_TransmissionLine+2*90.0* CLHEP::um,0),TransmissionLine_Resonator_SquareLogical_1_GP,"Resonator_TransmissionLine_1_GP",GroundPlaneLogical,false,0,true);
G4VPhysicalVolume* Resonator_TransmissionLine_2_GP=new G4PVPlacement(0,G4ThreeVector(1.0*dp_TransmissionLineBaseNbDimX-1*dp_TransmissionLineBaseNbDimX_1,Y_Center_TransmissionLine+2*90.0* CLHEP::um,0),TransmissionLine_Resonator_SquareLogical_1_GP,"Resonator_TransmissionLine_2_GP",GroundPlaneLogical,false,0,true);
G4VPhysicalVolume* Resonator_TransmissionLine_1=new G4PVPlacement(0,G4ThreeVector(),TransmissionLine_Resonator_SquareLogical_1,"Resonator_TransmissionLine_1",TransmissionLine_Resonator_SquareLogical_1_GP,false,0,true);
G4VPhysicalVolume* Resonator_TransmissionLine_2=new G4PVPlacement(0,G4ThreeVector(),TransmissionLine_Resonator_SquareLogical_1,"Resonator_TransmissionLine_2",TransmissionLine_Resonator_SquareLogical_1_GP,false,0,true);




///Creating the Panels of the Connecting Trannsmitions Lines. ///Changin to create the last part as a Triangle
constexpr double dp_PanelNbDimX1 = dp_resonatorAssemblyBaseNbDimZ_1;// Here must be the Depp og 1000
constexpr double dp_PanelNbDimY1 = 150.0* CLHEP::um;
constexpr double dp_PanelNbDimX2 = dp_resonatorAssemblyBaseNbDimZ_1;
constexpr double dp_PanelNbDimY2 = 10.0* CLHEP::um;
constexpr double dp_PanelNbDimZ = 150* CLHEP::um;
// G4Trd(const G4String& pName,
//             G4double  dx1,
//             G4double  dx2,
//             G4double  dy1,
//             G4double  dy2,
// G4double dz);
// G4RotationMatrix* yRot = new G4RotationMatrix;
// G4RotationMatrix* XYRot = new G4RotationMatrix;
// G4RotationMatrix* XYRotP = new G4RotationMatrix;
// yRot->rotateY(-1*M_PI/2.*rad);
// XYRot->rotateZ(-1*M_PI/2.*rad);
// XYRotP->rotateZ(1*M_PI/2.*rad);
// G4VSolid* Panel= new G4Trd("Panel",dp_PanelNbDimX1,dp_PanelNbDimX2,dp_PanelNbDimY1,dp_PanelNbDimY2,dp_PanelNbDimZ);
//G4LogicalVolume* PanelLogical = new G4LogicalVolume( Panel,fNiobium,"PanelLogical");


//Creating the Box atttached to the previuos
constexpr double dp_Square_panelNbDimX1 = 175.0 * CLHEP::um;
///Not Neccesaru deleting later
// G4VSolid* Square_panel = new G4Box("Square_panel", dp_Square_panelNbDimX1,dp_PanelNbDimY1,dp_resonatorAssemblyBaseNbDimZ_1);
// G4LogicalVolume* Square_panelLogical = new G4LogicalVolume( Square_panel,fNiobium,"Square_panelLogical");
// G4ThreeVector zTrans(2*160.0* CLHEP::um, 0, 0);
//
// G4UnionSolid* Transmtion_Panel =new G4UnionSolid("Transmtion_Panel", Square_panel, Panel, yRot,zTrans);
// G4LogicalVolume* Transmtion_PanelLogical = new G4LogicalVolume( Transmtion_Panel,fNiobium,"Transmtion_PanelLogical");
// ///GODDD
// for(G4int i=0; i<6;i++){
// if(i<1){
// G4VPhysicalVolume* Transmtion_PanelPhys =
//   new G4PVPlacement(0,G4ThreeVector(-1000.0*CLHEP::um,1000.0*CLHEP::um,0.37*cm),Transmtion_PanelLogical,"Transmtion_PanelPhys",
//                     worldLogical,false,i,true);}
//   else if( i>=1 && i<3){
// G4VPhysicalVolume* Transmtion_PanelPhys =
//                       new G4PVPlacement(XYRotP,G4ThreeVector(pow(-1,i)*(dp_TransmissionLineBaseNbDimX-2*dp_TransmissionLineBaseNbDimX_1-90.0* CLHEP::um),Y_Center_TransmissionLine+8*90.0* CLHEP::um+20*CLHEP::um,Position_Superconductor),Transmtion_PanelLogical,"Transmtion_PanelPhys",
//                                         worldLogical,false,i,true);}
//
// else{G4VPhysicalVolume* Transmtion_PanelPhys =
//                       new G4PVPlacement(XYRot,G4ThreeVector(-1000.0*CLHEP::um+800*i*CLHEP::um,-2000.0*CLHEP::um+i*10*CLHEP::um,0.37*cm),Transmtion_PanelLogical,"Transmtion_PanelPhys",
//                                         worldLogical,false,i,true);  }
//
//
// }
//



/////////////////Capacitor Connected to the Transmission Line 4 in total
///Definin the Rectangle
constexpr double dp_Square_Capacitor_to_TTLDimX1 = 480.0 * CLHEP::um;
constexpr double dp_Square_Capacitor_to_TTLDimY1  = 20.0 * CLHEP::um;


constexpr double dp_Square_connection_from_TTLD_to_ResonatorX1 = 10.0 * CLHEP::um;
constexpr double dp_Square_connection_from_TTLD_to_ResonatorY1  = 455.0 * CLHEP::um;

constexpr double dp_Square_connection_from_TTLD_to_ResonatorX2 = 10.0 * CLHEP::um;
constexpr double dp_Square_connection_from_TTLD_to_ResonatorY2  = 1400.0 * CLHEP::um;

G4VSolid* Square_Capacitor_to_TTLD = new G4Box("Square_Capacitor_to_TTLD",dp_Square_Capacitor_to_TTLDimX1 ,dp_Square_Capacitor_to_TTLDimY1 ,dp_resonatorAssemblyBaseNbDimZ_1);
G4LogicalVolume* Square_Capacitor_to_TTLDLogical = new G4LogicalVolume( Square_Capacitor_to_TTLD,fNiobium,"Square_Capacitor_to_TTLDLogical");
G4VSolid* Square_Capacitor_to_TTLD_GP = new G4Box("Square_Capacitor_to_TTLD_GP",dp_Square_Capacitor_to_TTLDimX1 ,dp_Square_Capacitor_to_TTLDimY1+Vacu_GroundPlane,dp_resonatorAssemblyBaseNbDimZ_1);
G4LogicalVolume* Square_Capacitor_to_TTLDLogical_GP = new G4LogicalVolume( Square_Capacitor_to_TTLD_GP,fLiquidHelium,"Square_Capacitor_to_TTLDLogical_GP");



///Adding the two Circles
G4VSolid* Half_Cirlce_R_Transmission_Square_Capacitor_to_TTDL= new G4Tubs("Half_Cirlce_R_Transmission_Square_Capacitor_to_TTDL",  0.0* CLHEP::um, 20.0* CLHEP::um,dp_resonatorAssemblyBaseNbDimZ_1,-1.0*tube_dPhi_4,2*tube_dPhi_4 ); // segment angle
G4LogicalVolume* Half_Cirlce_R_Transmission_Square_Capacitor_to_TTDLLogical = new G4LogicalVolume(Half_Cirlce_R_Transmission_Square_Capacitor_to_TTDL,fNiobium,"Half_Cirlce_R_Transmission_Square_Capacitor_to_TTDLLogical");
G4VSolid* Half_Cirlce_R_Transmission_Square_Capacitor_to_TTDL_GP= new G4Tubs("Half_Cirlce_R_Transmission_Square_Capacitor_to_TTDL_GP",  0.0* CLHEP::um, 20.0* CLHEP::um+Vacu_GroundPlane,dp_resonatorAssemblyBaseNbDimZ_1,-1.0*tube_dPhi_4,2*tube_dPhi_4 ); // segment angle
G4LogicalVolume* Half_Cirlce_R_Transmission_Square_Capacitor_to_TTDLLogical_GP = new G4LogicalVolume(Half_Cirlce_R_Transmission_Square_Capacitor_to_TTDL_GP,fLiquidHelium,"Half_Cirlce_R_Transmission_Square_Capacitor_to_TTDLLogical_GP");

G4VSolid* Half_Cirlce_L_Transmission_Square_Capacitor_to_TTDL= new G4Tubs("Half_Cirlce_L_Transmission_Square_Capacitor_to_TTDL",  0.0* CLHEP::um, 20.0* CLHEP::um,dp_resonatorAssemblyBaseNbDimZ_1,1.0*tube_dPhi_4,2*tube_dPhi_4 ); // segment angle
G4LogicalVolume* Half_Cirlce_L_Transmission_Square_Capacitor_to_TTDLLogical = new G4LogicalVolume(Half_Cirlce_L_Transmission_Square_Capacitor_to_TTDL,fNiobium,"Half_Cirlce_L_Transmission_Square_Capacitor_to_TTDLLogical");
G4VSolid* Half_Cirlce_L_Transmission_Square_Capacitor_to_TTDL_GP= new G4Tubs("Half_Cirlce_L_Transmission_Square_Capacitor_to_TTDL_GP",  0.0* CLHEP::um, 20.0* CLHEP::um+Vacu_GroundPlane,dp_resonatorAssemblyBaseNbDimZ_1,1.0*tube_dPhi_4,2*tube_dPhi_4 ); // segment angle
G4LogicalVolume* Half_Cirlce_L_Transmission_Square_Capacitor_to_TTDLLogical_GP = new G4LogicalVolume(Half_Cirlce_L_Transmission_Square_Capacitor_to_TTDL_GP,fLiquidHelium,"Half_Cirlce_L_Transmission_Square_Capacitor_to_TTDLLogical_GP");





///////Adding the Rectangles to the connection with the Resonantors
G4VSolid* Square_connection_from_TTLD_to_Resonator = new G4Box("Square_connection_from_TTLD_to_Resonator",dp_Square_connection_from_TTLD_to_ResonatorX1,dp_Square_connection_from_TTLD_to_ResonatorY1 ,dp_resonatorAssemblyBaseNbDimZ_1);
G4LogicalVolume* Square_connection_from_TTLD_to_ResonatorLogical = new G4LogicalVolume(Square_connection_from_TTLD_to_Resonator,fNiobium,"Square_connection_from_TTLD_to_ResonatorLogical");
G4VSolid* Square_connection_from_TTLD_to_Resonator_GP = new G4Box("Square_connection_from_TTLD_to_Resonator_GP",dp_Square_connection_from_TTLD_to_ResonatorX1+Vacu_GroundPlane,dp_Square_connection_from_TTLD_to_ResonatorY1 ,dp_resonatorAssemblyBaseNbDimZ_1);
G4LogicalVolume* Square_connection_from_TTLD_to_ResonatorLogical_GP = new G4LogicalVolume(Square_connection_from_TTLD_to_Resonator_GP,fLiquidHelium,"Square_connection_from_TTLD_to_ResonatorLogical_GP");


G4VSolid* Square_connection_from_TTLD_to_Resonator_B = new G4Box("Square_connection_from_TTLD_to_Resonator_B",dp_Square_connection_from_TTLD_to_ResonatorX2,dp_Square_connection_from_TTLD_to_ResonatorY2 ,dp_resonatorAssemblyBaseNbDimZ_1);
G4LogicalVolume* Square_connection_from_TTLD_to_ResonatorLogical_B = new G4LogicalVolume(Square_connection_from_TTLD_to_Resonator_B,fNiobium,"Square_connection_from_TTLD_to_ResonatorLogical_B");
G4VSolid* Square_connection_from_TTLD_to_Resonator_B_GP = new G4Box("Square_connection_from_TTLD_to_Resonator_B_GP",dp_Square_connection_from_TTLD_to_ResonatorX2+Vacu_GroundPlane,dp_Square_connection_from_TTLD_to_ResonatorY2 ,dp_resonatorAssemblyBaseNbDimZ_1);
G4LogicalVolume* Square_connection_from_TTLD_to_ResonatorLogical_B_GP = new G4LogicalVolume(Square_connection_from_TTLD_to_Resonator_B_GP,fLiquidHelium,"Square_connection_from_TTLD_to_ResonatorLogical_B_GP");

//G4double tube_dPhi_3 = 0.5*M_PI*rad;
///////////////////////////////Left quarter Circles conneted to the Transmissions Lines
G4VSolid* Half_Cirlce_L_TTDL= new G4Tubs("Half_Cirlce_L_TTDL",  30.0* CLHEP::um,50.0* CLHEP::um,dp_resonatorAssemblyBaseNbDimZ_1,M_PI*rad,0.5*M_PI*rad ); // segment angle
G4LogicalVolume* Half_Cirlce_L_TTDLogical = new G4LogicalVolume( Half_Cirlce_L_TTDL,fNiobium,"Half_Cirlce_L_TTDLogical");
G4VSolid* Half_Cirlce_L_TTDL_GP= new G4Tubs("Half_Cirlce_L_TTDL_GP",  30.0* CLHEP::um-Vacu_GroundPlane,50.0* CLHEP::um+Vacu_GroundPlane,dp_resonatorAssemblyBaseNbDimZ_1,M_PI*rad,0.5*M_PI*rad ); // segment angle
G4LogicalVolume* Half_Cirlce_L_TTDLogical_GP = new G4LogicalVolume( Half_Cirlce_L_TTDL_GP,fLiquidHelium,"Half_Cirlce_L_TTDLogical_GP");


G4VSolid* Half_Cirlce_R_TTDL= new G4Tubs("Half_Cirlce_R_TTDL",  30.0* CLHEP::um,50.0* CLHEP::um,dp_resonatorAssemblyBaseNbDimZ_1,-0.5*M_PI*rad,0.5*M_PI*rad ); // segment angle
G4LogicalVolume* Half_Cirlce_R_TTDLogical = new G4LogicalVolume( Half_Cirlce_R_TTDL,fNiobium,"Half_Cirlce_R_TTDLogical");
G4VSolid* Half_Cirlce_R_TTDL_GP= new G4Tubs("Half_Cirlce_R_TTDL_GP",  30.0* CLHEP::um-Vacu_GroundPlane,50.0* CLHEP::um+Vacu_GroundPlane,dp_resonatorAssemblyBaseNbDimZ_1,-0.5*M_PI*rad,0.5*M_PI*rad ); // segment angle
G4LogicalVolume* Half_Cirlce_R_TTDLogical_GP = new G4LogicalVolume( Half_Cirlce_R_TTDL_GP,fLiquidHelium,"Half_Cirlce_R_TTDLogical_GP");





for(G4int i=0; i<4;i++){
G4VPhysicalVolume* Square_Capacitor_to_TTLDPhys_GP=new G4PVPlacement(0,G4ThreeVector(X_Center_TransmissionLine-3*1430*CLHEP::um+2*265.0*CLHEP::um+i*(2*620*CLHEP::um+2*dp_Square_Capacitor_to_TTLDimX1),Y_Center_TransmissionLine-2*20.0*CLHEP::um-2*dp_TransmissionLineBaseNbDimY ,0),Square_Capacitor_to_TTLDLogical_GP,"Square_Capacitor_to_TTLDPhys_GP",GroundPlaneLogical,false,i,true);
G4CMPLogicalBorderSurface *border_sapphire_Vacuum_Square_Capacitor_to_TTLDPhys_GP=new G4CMPLogicalBorderSurface("border_sapphire_Vacuum_Square_Capacitor_to_TTLDPhys_GP",SapphirePhys,Square_Capacitor_to_TTLDPhys_GP,fSiVacuumInterface);
G4VPhysicalVolume* Square_Capacitor_to_TTLDPhys=new G4PVPlacement(0,G4ThreeVector(),Square_Capacitor_to_TTLDLogical,"Square_Capacitor_to_TTLDPhys",Square_Capacitor_to_TTLDLogical_GP,false,i,true);
G4CMPLogicalBorderSurface * border_sapphire_Vacuum_Square_Capacitor_to_TTLDPhys=new G4CMPLogicalBorderSurface("border_sapphire_Vacuum_Square_Capacitor_to_TTLDPhys",SapphirePhys,Square_Capacitor_to_TTLDPhys,fSiNbInterface);


G4VPhysicalVolume* Half_Cirlce_R_Transmission_Square_Capacitor_to_TTDLPhys_GP=new G4PVPlacement(0,G4ThreeVector(X_Center_TransmissionLine+dp_Square_Capacitor_to_TTLDimX1-3*1430*CLHEP::um+2*265.0*CLHEP::um+i*(2*620*CLHEP::um+2*dp_Square_Capacitor_to_TTLDimX1),Y_Center_TransmissionLine-2*20.0*CLHEP::um-2*dp_TransmissionLineBaseNbDimY ,0),Half_Cirlce_R_Transmission_Square_Capacitor_to_TTDLLogical_GP,"Half_Cirlce_R_Transmission_Square_Capacitor_to_TTDLPhys_GP",GroundPlaneLogical,false,i,true);
G4CMPLogicalBorderSurface * border_sapphire_Vacuum_Half_Cirlce_R_Transmission_Square_Capacitor_to_TTDLPhys_GP=new G4CMPLogicalBorderSurface("border_sapphire_Vacuum_Half_Cirlce_R_Transmission_Square_Capacitor_to_TTDLPhys_GP",SapphirePhys,Half_Cirlce_R_Transmission_Square_Capacitor_to_TTDLPhys_GP,fSiVacuumInterface);
G4VPhysicalVolume* Half_Cirlce_R_Transmission_Square_Capacitor_to_TTDLPhys=new G4PVPlacement(0,G4ThreeVector(),Half_Cirlce_R_Transmission_Square_Capacitor_to_TTDLLogical,"Half_Cirlce_R_Transmission_Square_Capacitor_to_TTDLPhys",Half_Cirlce_R_Transmission_Square_Capacitor_to_TTDLLogical_GP,false,i,true);
G4CMPLogicalBorderSurface * border_sapphire_Vacuum_Half_Cirlce_R_Transmission_Square_Capacitor_to_TTDLPhys=new G4CMPLogicalBorderSurface("border_sapphire_Vacuum_Half_Cirlce_R_Transmission_Square_Capacitor_to_TTDLPhys",SapphirePhys,Half_Cirlce_R_Transmission_Square_Capacitor_to_TTDLPhys,fSiNbInterface);


G4VPhysicalVolume* Half_Cirlce_L_Transmission_Square_Capacitor_to_TTDLPhys_GP=new G4PVPlacement(0,G4ThreeVector(X_Center_TransmissionLine-dp_Square_Capacitor_to_TTLDimX1-3*1430*CLHEP::um+2*265.0*CLHEP::um+i*(2*620*CLHEP::um+2*dp_Square_Capacitor_to_TTLDimX1),Y_Center_TransmissionLine-2*20.0*CLHEP::um-2*dp_TransmissionLineBaseNbDimY ,0),Half_Cirlce_L_Transmission_Square_Capacitor_to_TTDLLogical_GP,"Half_Cirlce_L_Transmission_Square_Capacitor_to_TTDLPhys_GP",GroundPlaneLogical,false,i,true);
G4CMPLogicalBorderSurface * border_sapphire_Vacuum_Half_Cirlce_L_Transmission_Square_Capacitor_to_TTDLPhys_GP=new G4CMPLogicalBorderSurface("border_sapphire_Vacuum_Half_Cirlce_L_Transmission_Square_Capacitor_to_TTDLPhys_GP",SapphirePhys,Half_Cirlce_L_Transmission_Square_Capacitor_to_TTDLPhys_GP,fSiVacuumInterface);
G4VPhysicalVolume* Half_Cirlce_L_Transmission_Square_Capacitor_to_TTDLPhys=new G4PVPlacement(0,G4ThreeVector(),Half_Cirlce_L_Transmission_Square_Capacitor_to_TTDLLogical,"Half_Cirlce_L_Transmission_Square_Capacitor_to_TTDLPhys",Half_Cirlce_L_Transmission_Square_Capacitor_to_TTDLLogical_GP,false,i,true);
G4CMPLogicalBorderSurface * border_sapphire_Vacuum_Half_Cirlce_L_Transmission_Square_Capacitor_to_TTDLPhys=new G4CMPLogicalBorderSurface("border_sapphire_Vacuum_Half_Cirlce_L_Transmission_Square_Capacitor_to_TTDLPhys",SapphirePhys,Half_Cirlce_L_Transmission_Square_Capacitor_to_TTDLPhys,fSiNbInterface);


if( i==0 or i==2 ){
G4VPhysicalVolume* Square_connection_from_TTLD_to_ResonatorPhys_GP=new G4PVPlacement(0,G4ThreeVector(X_Center_TransmissionLine-3*1430*CLHEP::um+2*265.0*CLHEP::um+i*(2*620*CLHEP::um+2*dp_Square_Capacitor_to_TTLDimX1),Y_Center_TransmissionLine-70.0*CLHEP::um-dp_TransmissionLineBaseNbDimY-dp_Square_connection_from_TTLD_to_ResonatorY1-dp_Square_Capacitor_to_TTLDimY1,0),Square_connection_from_TTLD_to_ResonatorLogical_GP,"Square_connection_from_TTLD_to_ResonatorPhys_GP",GroundPlaneLogical,false,i,true);
G4CMPLogicalBorderSurface * border_sapphire_Vacuum_Square_connection_from_TTLD_to_ResonatorPhys_GP=new G4CMPLogicalBorderSurface("border_sapphire_Vacuum_Square_connection_from_TTLD_to_ResonatorPhys_GP",SapphirePhys,Square_connection_from_TTLD_to_ResonatorPhys_GP,fSiVacuumInterface);

G4VPhysicalVolume* Square_connection_from_TTLD_to_ResonatorPhys=new G4PVPlacement(0,G4ThreeVector(),Square_connection_from_TTLD_to_ResonatorLogical,"Square_connection_from_TTLD_to_ResonatorPhys",Square_connection_from_TTLD_to_ResonatorLogical_GP,false,i,true);
G4CMPLogicalBorderSurface * border_sapphire_Vacuum_Square_connection_from_TTLD_to_ResonatorPhys=new G4CMPLogicalBorderSurface("border_sapphire_Vacuum_Square_connection_from_TTLD_to_ResonatorPhys",SapphirePhys,Square_connection_from_TTLD_to_ResonatorPhys,fSiNbInterface);

                                      }
else{
G4VPhysicalVolume* Square_connection_from_TTLD_to_ResonatorPhys_B_GP=new G4PVPlacement(0,G4ThreeVector(X_Center_TransmissionLine-3*1430*CLHEP::um+2*265.0*CLHEP::um+i*(2*620*CLHEP::um+2*dp_Square_Capacitor_to_TTLDimX1),Y_Center_TransmissionLine-70.0*CLHEP::um-dp_TransmissionLineBaseNbDimY-dp_Square_connection_from_TTLD_to_ResonatorY2-dp_Square_Capacitor_to_TTLDimY1,0),Square_connection_from_TTLD_to_ResonatorLogical_B_GP,"Square_connection_from_TTLD_to_ResonatorPhys_B_GP",GroundPlaneLogical,false,i,true);
G4CMPLogicalBorderSurface * border_sapphire_Vacuum_Square_connection_from_TTLD_to_ResonatorPhys_B_GP=new G4CMPLogicalBorderSurface("border_sapphire_Vacuum_Square_connection_from_TTLD_to_ResonatorPhys_B_GP",SapphirePhys,Square_connection_from_TTLD_to_ResonatorPhys_B_GP,fSiVacuumInterface);


G4VPhysicalVolume* Square_connection_from_TTLD_to_ResonatorPhys_B=new G4PVPlacement(0,G4ThreeVector(),Square_connection_from_TTLD_to_ResonatorLogical_B,"Square_connection_from_TTLD_to_ResonatorPhys_B",Square_connection_from_TTLD_to_ResonatorLogical_B_GP,false,i,true);
G4CMPLogicalBorderSurface * border_sapphire_Vacuum_Square_connection_from_TTLD_to_ResonatorPhys_B=new G4CMPLogicalBorderSurface("border_sapphire_Vacuum_Square_connection_from_TTLD_to_ResonatorPhys_B",SapphirePhys,Square_connection_from_TTLD_to_ResonatorPhys_B,fSiNbInterface);


  }


if( i==0 or i==2 ){

G4VPhysicalVolume* Half_Cirlce_L_TTDPhys_GP=new G4PVPlacement(0,G4ThreeVector(X_Center_TransmissionLine-3*1430*CLHEP::um+2*265.0*CLHEP::um+i*(2*620*CLHEP::um+2*dp_Square_Capacitor_to_TTLDimX1)+40.0*CLHEP::um,Y_Center_TransmissionLine-70.0*CLHEP::um-dp_TransmissionLineBaseNbDimY-2*dp_Square_connection_from_TTLD_to_ResonatorY1-dp_Square_Capacitor_to_TTLDimY1,0),Half_Cirlce_L_TTDLogical_GP,"Half_Cirlce_L_TTDPhys_GP",GroundPlaneLogical,false,i,true);
G4CMPLogicalBorderSurface * border_sapphire_Vacuum_Half_Cirlce_L_TTDPhys_GP=new G4CMPLogicalBorderSurface("border_sapphire_Vacuum_Half_Cirlce_L_TTDPhys_GP",SapphirePhys,Half_Cirlce_L_TTDPhys_GP,fSiVacuumInterface);

G4VPhysicalVolume* Half_Cirlce_L_TTDPhys=new G4PVPlacement(0,G4ThreeVector(),Half_Cirlce_L_TTDLogical,"Half_Cirlce_L_TTDPhys",Half_Cirlce_L_TTDLogical_GP,false,i,true);
G4CMPLogicalBorderSurface * border_sapphire_Vacuum_Half_Cirlce_L_TTDPhys=new G4CMPLogicalBorderSurface("border_sapphire_Vacuum_Half_Cirlce_L_TTDPhys",SapphirePhys,Half_Cirlce_L_TTDPhys,fSiNbInterface);

                  }
else if(i==1){
G4VPhysicalVolume* Half_Cirlce_L_TTDPhys_GP=new G4PVPlacement(0,G4ThreeVector(X_Center_TransmissionLine-3*1430*CLHEP::um+2*265.0*CLHEP::um+i*(2*620*CLHEP::um+2*dp_Square_Capacitor_to_TTLDimX1)+40.0*CLHEP::um,Y_Center_TransmissionLine-70.0*CLHEP::um-dp_TransmissionLineBaseNbDimY-2*dp_Square_connection_from_TTLD_to_ResonatorY2-dp_Square_Capacitor_to_TTLDimY1,0),Half_Cirlce_L_TTDLogical_GP,"Half_Cirlce_L_TTDPhys_GP",GroundPlaneLogical,false,i,true);
G4CMPLogicalBorderSurface * border_sapphire_Vacuum_Half_Cirlce_L_TTDPhys_GP=new G4CMPLogicalBorderSurface("border_sapphire_Vacuum_Half_Cirlce_L_TTDPhys_GP",SapphirePhys,Half_Cirlce_L_TTDPhys_GP,fSiVacuumInterface);

G4VPhysicalVolume* Half_Cirlce_L_TTDPhys=new G4PVPlacement(0,G4ThreeVector(),Half_Cirlce_L_TTDLogical,"Half_Cirlce_L_TTDPhys",Half_Cirlce_L_TTDLogical_GP,false,i,true);
G4CMPLogicalBorderSurface * border_sapphire_Vacuum_Half_Cirlce_L_TTDPhys=new G4CMPLogicalBorderSurface("border_sapphire_Vacuum_Half_Cirlce_L_TTDPhys",SapphirePhys,Half_Cirlce_L_TTDPhys,fSiNbInterface);
              }

else{G4VPhysicalVolume* Half_Cirlce_R_TTDPhys_GP=new G4PVPlacement(0,G4ThreeVector(X_Center_TransmissionLine-3*1430*CLHEP::um+2*265.0*CLHEP::um+i*(2*620*CLHEP::um+2*dp_Square_Capacitor_to_TTLDimX1)-40.0*CLHEP::um,Y_Center_TransmissionLine-70.0*CLHEP::um-dp_TransmissionLineBaseNbDimY-2*dp_Square_connection_from_TTLD_to_ResonatorY2-dp_Square_Capacitor_to_TTLDimY1,0),Half_Cirlce_R_TTDLogical_GP,"Half_Cirlce_R_TTDPhys_GP",GroundPlaneLogical,false,i,true);
G4CMPLogicalBorderSurface * border_sapphire_Vacuum_Half_Cirlce_R_TTDPhys_GP=new G4CMPLogicalBorderSurface("border_sapphire_Vacuum_Half_Cirlce_R_TTDPhys_GP",SapphirePhys,Half_Cirlce_R_TTDPhys_GP,fSiVacuumInterface);


G4VPhysicalVolume* Half_Cirlce_R_TTDPhys=new G4PVPlacement(0,G4ThreeVector(),Half_Cirlce_R_TTDLogical,"Half_Cirlce_R_TTDPhys",Half_Cirlce_R_TTDLogical_GP,false,i,true);
G4CMPLogicalBorderSurface * border_sapphire_Vacuum_Half_Cirlce_R_TTDPhys=new G4CMPLogicalBorderSurface("border_sapphire_Vacuum_Half_Cirlce_R_TTDPhys",SapphirePhys,Half_Cirlce_R_TTDPhys,fSiNbInterface);

                                      }



}
/////////////////////////--------------Creating the Square to connect with the  Resonators--------///////////
//Creating an array for the dimentions
constexpr double dp_Square_to_ResonatorX[4] = {280.0 * CLHEP::um,305.0 * CLHEP::um,332.0 * CLHEP::um,330.0*CLHEP::um};
constexpr double dp_Square_to_ResonatorY= 10.0*CLHEP::um;

G4VSolid* Square_To_Resonator1 = new G4Box("Square_To_Resonator1 ",dp_Square_to_ResonatorX[0] ,dp_Square_to_ResonatorY ,dp_resonatorAssemblyBaseNbDimZ_1);
G4VSolid* Square_To_Resonator2 = new G4Box("Square_To_Resonator2 ",dp_Square_to_ResonatorX[1] ,dp_Square_to_ResonatorY ,dp_resonatorAssemblyBaseNbDimZ_1);
G4VSolid* Square_To_Resonator3 = new G4Box("Square_To_Resonator3 ",dp_Square_to_ResonatorX[2] ,dp_Square_to_ResonatorY ,dp_resonatorAssemblyBaseNbDimZ_1);
G4VSolid* Square_To_Resonator4 = new G4Box("Square_To_Resonator4 ",dp_Square_to_ResonatorX[3] ,dp_Square_to_ResonatorY ,dp_resonatorAssemblyBaseNbDimZ_1);
G4VSolid* Square_To_Resonator1_GP = new G4Box("Square_To_Resonator1_GP ",dp_Square_to_ResonatorX[0] ,dp_Square_to_ResonatorY+Vacu_GroundPlane ,dp_resonatorAssemblyBaseNbDimZ_1);
G4VSolid* Square_To_Resonator2_GP = new G4Box("Square_To_Resonator2_GP ",dp_Square_to_ResonatorX[1] ,dp_Square_to_ResonatorY +Vacu_GroundPlane,dp_resonatorAssemblyBaseNbDimZ_1);
G4VSolid* Square_To_Resonator3_GP = new G4Box("Square_To_Resonator3_GP ",dp_Square_to_ResonatorX[2] ,dp_Square_to_ResonatorY +Vacu_GroundPlane,dp_resonatorAssemblyBaseNbDimZ_1);
G4VSolid* Square_To_Resonator4_GP = new G4Box("Square_To_Resonator4_GP ",dp_Square_to_ResonatorX[3] ,dp_Square_to_ResonatorY+Vacu_GroundPlane ,dp_resonatorAssemblyBaseNbDimZ_1);

G4LogicalVolume* Square_To_ResonatorLogical1 = new G4LogicalVolume( Square_To_Resonator1,fNiobium,"Square_To_ResonatorLogical1");
G4LogicalVolume* Square_To_ResonatorLogical2 = new G4LogicalVolume( Square_To_Resonator2,fNiobium,"Square_To_ResonatorLogical2");
G4LogicalVolume* Square_To_ResonatorLogical3 = new G4LogicalVolume( Square_To_Resonator3,fNiobium,"Square_To_ResonatorLogical3");
G4LogicalVolume* Square_To_ResonatorLogical4 = new G4LogicalVolume( Square_To_Resonator4,fNiobium,"Square_To_ResonatorLogical4");
G4LogicalVolume* Square_To_ResonatorLogical1_GP = new G4LogicalVolume( Square_To_Resonator1_GP,fLiquidHelium,"Square_To_ResonatorLogical1_GP");
G4LogicalVolume* Square_To_ResonatorLogical2_GP = new G4LogicalVolume( Square_To_Resonator2_GP,fLiquidHelium,"Square_To_ResonatorLogical2_GP");
G4LogicalVolume* Square_To_ResonatorLogical3_GP = new G4LogicalVolume( Square_To_Resonator3_GP,fLiquidHelium,"Square_To_ResonatorLogical3_GP");
G4LogicalVolume* Square_To_ResonatorLogical4_GP = new G4LogicalVolume( Square_To_Resonator4_GP,fLiquidHelium,"Square_To_ResonatorLogical4_GP");

G4VPhysicalVolume* Square_To_ResonatorLogicalPhys1_GP=new G4PVPlacement(0,G4ThreeVector(X_Center_TransmissionLine+dp_Square_to_ResonatorX[0]-3*1430*CLHEP::um+2*265.0*CLHEP::um+0*(2*620*CLHEP::um+2*dp_Square_Capacitor_to_TTLDimX1)+40.0*CLHEP::um,Y_Center_TransmissionLine-70.0*CLHEP::um-dp_TransmissionLineBaseNbDimY-2*dp_Square_connection_from_TTLD_to_ResonatorY1-dp_Square_Capacitor_to_TTLDimY1-4*dp_Square_to_ResonatorY,0),Square_To_ResonatorLogical1_GP,"Square_To_ResonatorLogicalPhys1_GP",
                                                                                GroundPlaneLogical,false,0,true);
G4VPhysicalVolume* Square_To_ResonatorLogicalPhys1=new G4PVPlacement(0,G4ThreeVector(),Square_To_ResonatorLogical1,"Square_To_ResonatorLogicalPhys1",
                                        Square_To_ResonatorLogical1_GP,false,0,true);


G4VPhysicalVolume* Square_To_ResonatorLogicalPhys2_GP=new G4PVPlacement(0,G4ThreeVector(X_Center_TransmissionLine+dp_Square_to_ResonatorX[1]-3*1430*CLHEP::um+2*265.0*CLHEP::um+1*(2*620*CLHEP::um+2*dp_Square_Capacitor_to_TTLDimX1)+40.0*CLHEP::um,Y_Center_TransmissionLine-70.0*CLHEP::um-dp_TransmissionLineBaseNbDimY-2*dp_Square_connection_from_TTLD_to_ResonatorY2-dp_Square_Capacitor_to_TTLDimY1-4*dp_Square_to_ResonatorY,0),Square_To_ResonatorLogical2_GP,"Square_To_ResonatorLogicalPhys2_GP",GroundPlaneLogical,false,0,true);
G4VPhysicalVolume* Square_To_ResonatorLogicalPhys2=new G4PVPlacement(0,G4ThreeVector(),Square_To_ResonatorLogical2,"Square_To_ResonatorLogicalPhys2",
                                                                              Square_To_ResonatorLogical2_GP,false,0,true);


G4VPhysicalVolume* Square_To_ResonatorLogicalPhys3_GP=new G4PVPlacement(0,G4ThreeVector(X_Center_TransmissionLine+dp_Square_to_ResonatorX[2]-3*1430*CLHEP::um+2*265.0*CLHEP::um+2*(2*620*CLHEP::um+2*dp_Square_Capacitor_to_TTLDimX1)+40.0*CLHEP::um,Y_Center_TransmissionLine-70.0*CLHEP::um-dp_TransmissionLineBaseNbDimY-2*dp_Square_connection_from_TTLD_to_ResonatorY1-dp_Square_Capacitor_to_TTLDimY1-4*dp_Square_to_ResonatorY,0),Square_To_ResonatorLogical3_GP,"Square_To_ResonatorLogicalPhys3_GP",
                                                                                                GroundPlaneLogical,false,0,true);
G4VPhysicalVolume* Square_To_ResonatorLogicalPhys3=new G4PVPlacement(0,G4ThreeVector(),Square_To_ResonatorLogical3,"Square_To_ResonatorLogicalPhys3",
                  Square_To_ResonatorLogical3_GP,false,0,true);


G4VPhysicalVolume* Square_To_ResonatorLogicalPhys4_GP=new G4PVPlacement(0,G4ThreeVector(X_Center_TransmissionLine-dp_Square_to_ResonatorX[3]-3*1430*CLHEP::um+2*265.0*CLHEP::um+3*(2*620*CLHEP::um+2*dp_Square_Capacitor_to_TTLDimX1)-40.0*CLHEP::um,Y_Center_TransmissionLine-70.0*CLHEP::um-dp_TransmissionLineBaseNbDimY-2*dp_Square_connection_from_TTLD_to_ResonatorY2-dp_Square_Capacitor_to_TTLDimY1-4*dp_Square_to_ResonatorY,0),Square_To_ResonatorLogical4_GP,"Square_To_ResonatorLogicalPhys4_GP",
                                  GroundPlaneLogical,false,0,true);
G4VPhysicalVolume* Square_To_ResonatorLogicalPhys4=new G4PVPlacement(0,G4ThreeVector(),Square_To_ResonatorLogical4,"Square_To_ResonatorLogicalPhys4",Square_To_ResonatorLogical4_GP,false,0,true);


////////////////////////////------------------------Here Finish the Transmission for the connection with the Resonators.

////////////////////Implementing the Rsonators
constexpr double X_Center_1 =X_Center_TransmissionLine+dp_Square_to_ResonatorX[0]-3*1430*CLHEP::um+2*265.0*CLHEP::um+0*(2*620*CLHEP::um+2*dp_Square_Capacitor_to_TTLDimX1)+40.0*CLHEP::um-0.5*dp_resonatorAssemblyBaseNbDimX_1-1*42.0*CLHEP::um;
constexpr double Y_Center_1 =Y_Center_TransmissionLine-70.0*CLHEP::um-dp_TransmissionLineBaseNbDimY-2*dp_Square_connection_from_TTLD_to_ResonatorY1-dp_Square_Capacitor_to_TTLDimY1-4*dp_Square_to_ResonatorY-14*40.0*CLHEP::um;
for(G4int i=0; i<8;i++){
  if(i<7){
G4VPhysicalVolume* Resonator_SquarePhys_1_GP =new G4PVPlacement(0,G4ThreeVector(X_Center_1,i*80.0*CLHEP::um+Y_Center_1,0),Resonator_SquareLogical_1_GP,"Resonator_Square_1_GP",GroundPlaneLogical,false,i,true);
G4CMPLogicalBorderSurface * border_sapphire_Vacuum_Resonator_SquarePhys_1_GP=new G4CMPLogicalBorderSurface("border_sapphire_Vacuum_Resonator_SquarePhys_1_GP",SapphirePhys,Resonator_SquarePhys_1_GP,fSiVacuumInterface);

G4VPhysicalVolume* Resonator_SquarePhys_1 =new G4PVPlacement(0,G4ThreeVector(),Resonator_SquareLogical_1,"Resonator_Square_1",Resonator_SquareLogical_1_GP,false,i,true);
G4CMPLogicalBorderSurface * border_sapphire_Vacuum_Resonator_SquarePhys_1=new G4CMPLogicalBorderSurface("border_sapphire_Vacuum_Resonator_SquarePhys_1",SapphirePhys,Resonator_SquarePhys_1,fSiNbInterface);


}

      if (pow(-1,i)>0) {
          G4VPhysicalVolume* Half_Cirlce_RPhys_1_GP =new G4PVPlacement(0,G4ThreeVector(X_Center_1+pow(-1,i)*dp_resonatorAssemblyBaseNbDimX_1,40.0*CLHEP::um+i*80.0*CLHEP::um+Y_Center_1,0),Half_Cirlce_RLogical_1_GP,"Half_Cirlce_R_1_GP",GroundPlaneLogical,false,i,true);
          G4CMPLogicalBorderSurface * border_sapphire_Vacuum_Half_Cirlce_RPhys_1_GP=new G4CMPLogicalBorderSurface("border_sapphire_Vacuum_Half_Cirlce_RPhys_1_GP",SapphirePhys,Half_Cirlce_RPhys_1_GP,fSiVacuumInterface);

          G4VPhysicalVolume* Half_Cirlce_RPhys_1 =new G4PVPlacement(0,G4ThreeVector(),Half_Cirlce_RLogical_1,"Half_Cirlce_R_1",Half_Cirlce_RLogical_1_GP,false,i,true);
          G4CMPLogicalBorderSurface * border_sapphire_Vacuum_Half_Cirlce_RPhys_1=new G4CMPLogicalBorderSurface("border_sapphire_Vacuum_Half_Cirlce_RPhys_1",SapphirePhys,Half_Cirlce_RPhys_1,fSiNbInterface);


                      }
        else{
    G4VPhysicalVolume* Half_Cirlce_LPhys_1_GP =new G4PVPlacement(0,G4ThreeVector(X_Center_1+pow(-1,i)*dp_resonatorAssemblyBaseNbDimX_1,-120.0*CLHEP::um+i*80.0*CLHEP::um+Y_Center_1,0),Half_Cirlce_LLogical_1_GP,"Half_Cirlce_L_1_GP",GroundPlaneLogical,false,i,true);
    G4CMPLogicalBorderSurface * border_sapphire_Vacuum_Half_Cirlce_LPhys_1_GP=new G4CMPLogicalBorderSurface("border_sapphire_Vacuum_Half_Cirlce_LPhys_1_GP",SapphirePhys,Half_Cirlce_LPhys_1_GP,fSiVacuumInterface);

  G4VPhysicalVolume* Half_Cirlce_LPhys_1 =new G4PVPlacement(0,G4ThreeVector(),Half_Cirlce_LLogical_1,"Half_Cirlce_L_1",Half_Cirlce_LLogical_1_GP,false,i,true);
  G4CMPLogicalBorderSurface * border_sapphire_Vacuum_Half_Cirlce_LPhys_1=new G4CMPLogicalBorderSurface("border_sapphire_Vacuum_Half_Cirlce_LPhys_1",SapphirePhys,Half_Cirlce_LPhys_1,fSiNbInterface);




                                                                                                    }
}//////////For for Resonator 1
///////////////----------Resonator2
constexpr double X_Center_2 =X_Center_TransmissionLine+dp_Square_to_ResonatorX[1]-3*1430*CLHEP::um+2*265.0*CLHEP::um+1*(2*620*CLHEP::um+2*dp_Square_Capacitor_to_TTLDimX1)+40.0*CLHEP::um-0.5*dp_resonatorAssemblyBaseNbDimX_1-1*67.0*CLHEP::um;
constexpr double Y_Center_2 =Y_Center_TransmissionLine-70.0*CLHEP::um-dp_TransmissionLineBaseNbDimY-2*dp_Square_connection_from_TTLD_to_ResonatorY2-dp_Square_Capacitor_to_TTLDimY1-4*dp_Square_to_ResonatorY-14*40.0*CLHEP::um;
for(G4int i=0; i<8;i++){
  if(i<7){
  G4VPhysicalVolume* Resonator_SquarePhys_2_GP =new G4PVPlacement(0,G4ThreeVector(X_Center_2,i*80.0*CLHEP::um+Y_Center_2,0),Resonator_SquareLogical_2_GP,"Resonator_Square_2_GP",GroundPlaneLogical,false,i,true);
  G4CMPLogicalBorderSurface * border_sapphire_Vacuum_Resonator_SquarePhys_2_GP=new G4CMPLogicalBorderSurface("border_sapphire_Vacuum_Resonator_SquarePhys_2_GP",SapphirePhys,Resonator_SquarePhys_2_GP,fSiVacuumInterface);

  G4VPhysicalVolume* Resonator_SquarePhys_2 =new G4PVPlacement(0,G4ThreeVector(),Resonator_SquareLogical_2,"Resonator_Square_2",Resonator_SquareLogical_2_GP,false,i,true);
  G4CMPLogicalBorderSurface * border_sapphire_Vacuum_Resonator_SquarePhys_2=new G4CMPLogicalBorderSurface("border_sapphire_Vacuum_Resonator_SquarePhys_2",SapphirePhys,Resonator_SquarePhys_2,fSiNbInterface);


}

    if (pow(-1,i)>0) {
      G4VPhysicalVolume* Half_Cirlce_RPhys_2_GP =new G4PVPlacement(0,G4ThreeVector(X_Center_2+pow(-1,i)*dp_resonatorAssemblyBaseNbDimX_2,40.0*CLHEP::um+i*80.0*CLHEP::um+Y_Center_2,0),Half_Cirlce_RLogical_2_GP,"Half_Cirlce_R_2_GP",GroundPlaneLogical,false,i,true);
      G4CMPLogicalBorderSurface * border_sapphire_Vacuum_Half_Cirlce_RPhys_2_GP=new G4CMPLogicalBorderSurface("border_sapphire_Vacuum_Half_Cirlce_RPhys_2_GP",SapphirePhys,Half_Cirlce_RPhys_2_GP,fSiVacuumInterface);

      G4VPhysicalVolume* Half_Cirlce_RPhys_2 =new G4PVPlacement(0,G4ThreeVector(),Half_Cirlce_RLogical_2,"Half_Cirlce_R_2",Half_Cirlce_RLogical_2_GP,false,i,true);
      G4CMPLogicalBorderSurface * border_sapphire_Vacuum_Half_Cirlce_RPhys_2=new G4CMPLogicalBorderSurface("border_sapphire_Vacuum_Half_Cirlce_RPhys_2",SapphirePhys,Half_Cirlce_RPhys_2,fSiNbInterface);
                    }
  else{
    G4VPhysicalVolume* Half_Cirlce_LPhys_2_GP =new G4PVPlacement(0,G4ThreeVector(X_Center_2+pow(-1,i)*dp_resonatorAssemblyBaseNbDimX_2,-120.0*CLHEP::um+i*80.0*CLHEP::um+Y_Center_2,0),Half_Cirlce_LLogical_2_GP,"Half_Cirlce_L_2_GP",GroundPlaneLogical,false,i,true);
    G4CMPLogicalBorderSurface * border_sapphire_Vacuum_Half_Cirlce_LPhys_2_GP=new G4CMPLogicalBorderSurface("border_sapphire_Vacuum_Half_Cirlce_LPhys_2_GP",SapphirePhys,Half_Cirlce_LPhys_2_GP,fSiVacuumInterface);

    G4VPhysicalVolume* Half_Cirlce_LPhys_2 =new G4PVPlacement(0,G4ThreeVector(),Half_Cirlce_LLogical_2,"Half_Cirlce_L_2",Half_Cirlce_LLogical_2_GP,false,i,true);
    G4CMPLogicalBorderSurface * border_sapphire_Vacuum_Half_Cirlce_LPhys_2=new G4CMPLogicalBorderSurface("border_sapphire_Vacuum_Half_Cirlce_LPhys_2",SapphirePhys,Half_Cirlce_LPhys_2,fSiNbInterface);


      }


}
////////////////////RESONATOR 3 -------------//////////
constexpr double X_Center_3 =X_Center_TransmissionLine+dp_Square_to_ResonatorX[2]-3*1430*CLHEP::um+2*265.0*CLHEP::um+2*(2*620*CLHEP::um+2*dp_Square_Capacitor_to_TTLDimX1)+40.0*CLHEP::um-0.5*dp_resonatorAssemblyBaseNbDimX_1-1*89.0*CLHEP::um;
constexpr double Y_Center_3 =Y_Center_TransmissionLine-70.0*CLHEP::um-dp_TransmissionLineBaseNbDimY-2*dp_Square_connection_from_TTLD_to_ResonatorY1-dp_Square_Capacitor_to_TTLDimY1-4*dp_Square_to_ResonatorY-14*40.0*CLHEP::um;

for(G4int i=0; i<8;i++){
  if(i<7){

G4VPhysicalVolume* Resonator_SquarePhys_3_GP =new G4PVPlacement(0,G4ThreeVector(X_Center_3,i*80.0*CLHEP::um+Y_Center_3,0),Resonator_SquareLogical_3_GP,"Resonator_Square_3_GP",GroundPlaneLogical,false,i,true);
G4CMPLogicalBorderSurface * border_sapphire_Vacuum_Resonator_SquarePhys_3_GP=new G4CMPLogicalBorderSurface("border_sapphire_Vacuum_Resonator_SquarePhys_3_GP",SapphirePhys,Resonator_SquarePhys_3_GP,fSiVacuumInterface);

G4VPhysicalVolume* Resonator_SquarePhys_3 =new G4PVPlacement(0,G4ThreeVector(),Resonator_SquareLogical_3,"Resonator_Square_3",Resonator_SquareLogical_3_GP,false,i,true);
G4CMPLogicalBorderSurface * border_sapphire_Vacuum_Resonator_SquarePhys_3=new G4CMPLogicalBorderSurface("border_sapphire_Vacuum_Resonator_SquarePhys_3",SapphirePhys,Resonator_SquarePhys_3,fSiNbInterface);

}

  if (pow(-1,i)>0) {
      G4VPhysicalVolume* Half_Cirlce_RPhys_3_GP =new G4PVPlacement(0,G4ThreeVector(X_Center_3+pow(-1,i)*dp_resonatorAssemblyBaseNbDimX_3,40.0*CLHEP::um+i*80.0*CLHEP::um+Y_Center_3,0),Half_Cirlce_RLogical_3_GP,"Half_Cirlce_R_3_GP",GroundPlaneLogical,false,i,true);
      G4CMPLogicalBorderSurface * border_sapphire_Vacuum_Half_Cirlce_RPhys_3_GP=new G4CMPLogicalBorderSurface("border_sapphire_Vacuum_Half_Cirlce_RPhys_3_GP",SapphirePhys,Half_Cirlce_RPhys_3_GP,fSiVacuumInterface);

      G4VPhysicalVolume* Half_Cirlce_RPhys_3 =new G4PVPlacement(0,G4ThreeVector(),Half_Cirlce_RLogical_3,"Half_Cirlce_R_3",Half_Cirlce_RLogical_3_GP,false,i,true);
      G4CMPLogicalBorderSurface * border_sapphire_Vacuum_Half_Cirlce_RPhys_3=new G4CMPLogicalBorderSurface("border_sapphire_Vacuum_Half_Cirlce_RPhys_3",SapphirePhys,Half_Cirlce_RPhys_3,fSiNbInterface);


                  }
  else{
        G4VPhysicalVolume* Half_Cirlce_LPhys_3_GP =new G4PVPlacement(0,G4ThreeVector(X_Center_3+pow(-1,i)*dp_resonatorAssemblyBaseNbDimX_3,-120.0*CLHEP::um+i*80.0*CLHEP::um+Y_Center_3,0),Half_Cirlce_LLogical_3_GP,"Half_Cirlce_L_3_GP",GroundPlaneLogical,false,i,true);
        G4CMPLogicalBorderSurface * border_sapphire_Vacuum_Half_Cirlce_LPhys_3_GP=new G4CMPLogicalBorderSurface("border_sapphire_Vacuum_Half_Cirlce_LPhys_3_GP",SapphirePhys,Half_Cirlce_LPhys_3_GP,fSiVacuumInterface);

        G4VPhysicalVolume* Half_Cirlce_LPhys_3 =new G4PVPlacement(0,G4ThreeVector(),Half_Cirlce_LLogical_3,"Half_Cirlce_L_3",Half_Cirlce_LLogical_3_GP,false,i,true);
        G4CMPLogicalBorderSurface * border_sapphire_Vacuum_Half_Cirlce_LPhys_3=new G4CMPLogicalBorderSurface("border_sapphire_Vacuum_Half_Cirlce_LPhys_3",SapphirePhys,Half_Cirlce_LPhys_3,fSiNbInterface);

      }


}
//////////////////Resonator 4
constexpr double X_Center_4 =X_Center_TransmissionLine-dp_Square_to_ResonatorX[3]-3*1430*CLHEP::um+2*265.0*CLHEP::um+3*(2*620*CLHEP::um+2*dp_Square_Capacitor_to_TTLDimX1)-40.0*CLHEP::um+0.5*dp_resonatorAssemblyBaseNbDimX_1+1*91.0*CLHEP::um;
constexpr double Y_Center_4 =Y_Center_TransmissionLine-70.0*CLHEP::um-dp_TransmissionLineBaseNbDimY-2*dp_Square_connection_from_TTLD_to_ResonatorY2-dp_Square_Capacitor_to_TTLDimY1-4*dp_Square_to_ResonatorY-14*40.0*CLHEP::um;
for(G4int i=0; i<8;i++){
  if(i<7){
  G4VPhysicalVolume* Resonator_SquarePhys_4_GP =new G4PVPlacement(0,G4ThreeVector(X_Center_4,i*80.0*CLHEP::um+Y_Center_4,0),Resonator_SquareLogical_4_GP,"Resonator_Square_4_GP",GroundPlaneLogical,false,i,true);
  G4CMPLogicalBorderSurface * border_sapphire_Vacuum_Resonator_SquarePhys_4_GP=new G4CMPLogicalBorderSurface("border_sapphire_Vacuum_Resonator_SquarePhys_4_GP",SapphirePhys,Resonator_SquarePhys_4_GP,fSiVacuumInterface);

  G4VPhysicalVolume* Resonator_SquarePhys_4 =new G4PVPlacement(0,G4ThreeVector(),Resonator_SquareLogical_4,"Resonator_Square_4",Resonator_SquareLogical_4_GP,false,i,true);
  G4CMPLogicalBorderSurface * border_sapphire_Vacuum_Resonator_SquarePhys_4=new G4CMPLogicalBorderSurface("border_sapphire_Vacuum_Resonator_SquarePhys_4",SapphirePhys,Resonator_SquarePhys_4,fSiNbInterface);


}

  if (pow(-1,i)>0 ) {
    G4VPhysicalVolume* Half_Cirlce_RPhys_4_GP =new G4PVPlacement(0,G4ThreeVector(X_Center_4+pow(-1,i)*dp_resonatorAssemblyBaseNbDimX_4,40.0*CLHEP::um+i*80.0*CLHEP::um+Y_Center_4-80.0*CLHEP::um,0),Half_Cirlce_RLogical_4_GP,"Half_Cirlce_R_4_GP",GroundPlaneLogical,false,i,true);
    G4CMPLogicalBorderSurface * border_sapphire_Vacuum_Half_Cirlce_RPhys_4_GP=new G4CMPLogicalBorderSurface("border_sapphire_Vacuum_Half_Cirlce_RPhys_4_GP",SapphirePhys,Half_Cirlce_RPhys_4_GP,fSiVacuumInterface);

    G4VPhysicalVolume* Half_Cirlce_RPhys_4 =new G4PVPlacement(0,G4ThreeVector(),Half_Cirlce_RLogical_4,"Half_Cirlce_R_4",Half_Cirlce_RLogical_4_GP,false,i,true);
    G4CMPLogicalBorderSurface * border_sapphire_Vacuum_Half_Cirlce_RPhys_4=new G4CMPLogicalBorderSurface("border_sapphire_Vacuum_Half_Cirlce_RPhys_4",SapphirePhys,Half_Cirlce_RPhys_4,fSiNbInterface);


                      }
  else{G4VPhysicalVolume* Half_Cirlce_LPhys_4_GP =new G4PVPlacement(0,G4ThreeVector(X_Center_4+pow(-1,i)*dp_resonatorAssemblyBaseNbDimX_4,-120.0*CLHEP::um+i*80.0*CLHEP::um+Y_Center_4+80.0*CLHEP::um,0),Half_Cirlce_LLogical_4_GP,"Half_Cirlce_L_4_GP",GroundPlaneLogical,false,i,true);
    G4CMPLogicalBorderSurface * border_sapphire_Vacuum_Half_Cirlce_LPhys_4_GP=new G4CMPLogicalBorderSurface("border_sapphire_Vacuum_Half_Cirlce_LPhys_4_GP",SapphirePhys,Half_Cirlce_LPhys_4_GP,fSiVacuumInterface);

      G4VPhysicalVolume* Half_Cirlce_LPhys_4 =new G4PVPlacement(0,G4ThreeVector(),Half_Cirlce_LLogical_4,"Half_Cirlce_L_4",Half_Cirlce_LLogical_4_GP,false,i,true);

      G4CMPLogicalBorderSurface * border_sapphire_Vacuum_Half_Cirlce_LPhys_4=new G4CMPLogicalBorderSurface("border_sapphire_Vacuum_Half_Cirlce_LPhys_4",SapphirePhys,Half_Cirlce_LPhys_4,fSiNbInterface);

                                                            }


}

////////////////////////////Creating the Connections to the resonator to the squares and Qubit------
constexpr double dp_Square_to_CapcitorQX1[4] = {280.0*CLHEP::um,300.0*CLHEP::um,330.0*CLHEP::um,330.0*CLHEP::um};
constexpr double dp_Square_to_CapcitorQY1= 10.0*CLHEP::um;
constexpr double dp_Square_to_CapcitorQX2[4] = {670.0 * CLHEP::um,350.0 * CLHEP::um,30.0 * CLHEP::um,130.0*CLHEP::um};
constexpr double dp_Square_to_CapcitorQY2= 10.0*CLHEP::um;
///Vertical
constexpr double dp_Square_to_CapcitorQXV2 = 10.0*CLHEP::um;
constexpr double dp_Square_to_CapcitorQYV2[2]= {1243.0*CLHEP::um,465.0*CLHEP::um};
constexpr double dp_Square_to_CapcitorQBX2 = 10.0*CLHEP::um;
constexpr double dp_Square_to_CapcitorQBY2[4]= {315.0*CLHEP::um,610.0*CLHEP::um,1092.0*CLHEP::um,610.0*CLHEP::um};

G4VSolid* Square_to_CapcitorH1 = new G4Box("Square_to_CapcitorH1",dp_Square_to_CapcitorQX1[0] ,dp_Square_to_CapcitorQY1 ,dp_resonatorAssemblyBaseNbDimZ_1);
G4LogicalVolume* Square_to_CapcitorH1Logical = new G4LogicalVolume( Square_to_CapcitorH1,fNiobium,"Square_to_CapcitorH1Logical");
G4VSolid* Square_to_CapcitorH1_GP = new G4Box("Square_to_CapcitorH1_GP",dp_Square_to_CapcitorQX1[0] ,dp_Square_to_CapcitorQY1+Vacu_GroundPlane ,dp_resonatorAssemblyBaseNbDimZ_1);
G4LogicalVolume* Square_to_CapcitorH1Logical_GP = new G4LogicalVolume( Square_to_CapcitorH1_GP,fLiquidHelium,"Square_to_CapcitorH1Logical_GP");

G4VSolid* Square_to_CapcitorH2 = new G4Box("Square_to_CapcitorH2",dp_Square_to_CapcitorQX1[1] ,dp_Square_to_CapcitorQY1 ,dp_resonatorAssemblyBaseNbDimZ_1);
G4LogicalVolume* Square_to_CapcitorH2Logical = new G4LogicalVolume( Square_to_CapcitorH2,fNiobium,"Square_to_CapcitorH2Logical");
G4VSolid* Square_to_CapcitorH2_GP = new G4Box("Square_to_CapcitorH2_GP",dp_Square_to_CapcitorQX1[1] ,dp_Square_to_CapcitorQY1+Vacu_GroundPlane ,dp_resonatorAssemblyBaseNbDimZ_1);
G4LogicalVolume* Square_to_CapcitorH2Logical_GP = new G4LogicalVolume( Square_to_CapcitorH2_GP,fLiquidHelium,"Square_to_CapcitorH2Logical_GP");

G4VSolid* Square_to_CapcitorH3 = new G4Box("Square_to_CapcitorH3 ",dp_Square_to_CapcitorQX1[2] ,dp_Square_to_CapcitorQY1 ,dp_resonatorAssemblyBaseNbDimZ_1);
G4LogicalVolume* Square_to_CapcitorH3Logical = new G4LogicalVolume( Square_to_CapcitorH3,fNiobium,"Square_to_CapcitorH3Logical");
G4VSolid* Square_to_CapcitorH3_GP = new G4Box("Square_to_CapcitorH3_GP ",dp_Square_to_CapcitorQX1[2] ,dp_Square_to_CapcitorQY1+Vacu_GroundPlane ,dp_resonatorAssemblyBaseNbDimZ_1);
G4LogicalVolume* Square_to_CapcitorH3Logical_GP = new G4LogicalVolume( Square_to_CapcitorH3_GP,fLiquidHelium,"Square_to_CapcitorH3Logical_GP");

G4VSolid* Square_to_CapcitorH4 = new G4Box("Square_to_CapcitorH4",dp_Square_to_CapcitorQX1[3] ,dp_Square_to_CapcitorQY1 ,dp_resonatorAssemblyBaseNbDimZ_1);
G4LogicalVolume* Square_to_CapcitorH4Logical = new G4LogicalVolume( Square_to_CapcitorH4,fNiobium,"Square_to_CapcitorH4Logical");
G4VSolid* Square_to_CapcitorH4_GP = new G4Box("Square_to_CapcitorH4_GP",dp_Square_to_CapcitorQX1[3] ,dp_Square_to_CapcitorQY1+Vacu_GroundPlane ,dp_resonatorAssemblyBaseNbDimZ_1);
G4LogicalVolume* Square_to_CapcitorH4Logical_GP = new G4LogicalVolume( Square_to_CapcitorH4_GP,fLiquidHelium,"Square_to_CapcitorH4Logical_GP");



G4VSolid* Half_Cirlce_L_SquareCapacitor= new G4Tubs("Half_Cirlce_L_SquareCapacitor",  30.0* CLHEP::um,50.0* CLHEP::um,dp_resonatorAssemblyBaseNbDimZ_1,0.5*M_PI*rad,0.5*M_PI*rad ); // segment angle
G4LogicalVolume* Half_Cirlce_L_SquareCapacitorLogical = new G4LogicalVolume( Half_Cirlce_L_SquareCapacitor,fNiobium,"Half_Cirlce_L_SquareCapacitorLogical");
G4VSolid* Half_Cirlce_L_SquareCapacitor_GP= new G4Tubs("Half_Cirlce_L_SquareCapacitor_GP",  30.0* CLHEP::um-Vacu_GroundPlane,50.0* CLHEP::um+Vacu_GroundPlane,dp_resonatorAssemblyBaseNbDimZ_1,0.5*M_PI*rad,0.5*M_PI*rad ); // segment angle
G4LogicalVolume* Half_Cirlce_L_SquareCapacitorLogical_GP = new G4LogicalVolume( Half_Cirlce_L_SquareCapacitor_GP,fLiquidHelium,"Half_Cirlce_L_SquareCapacitorLogical_GP");

G4VSolid* Half_Cirlce_R_SquareCapacitor= new G4Tubs("Half_Cirlce_R_SquareCapacitor",  30.0* CLHEP::um,50.0* CLHEP::um,dp_resonatorAssemblyBaseNbDimZ_1,0.0*M_PI*rad,0.5*M_PI*rad ); // segment angle
G4LogicalVolume* Half_Cirlce_R_SquareCapacitorLogical = new G4LogicalVolume( Half_Cirlce_R_SquareCapacitor,fNiobium,"Half_Cirlce_R_SquareCapacitorLogical");
G4VSolid* Half_Cirlce_R_SquareCapacitor_GP= new G4Tubs("Half_Cirlce_R_SquareCapacitor_GP",  30.0* CLHEP::um-Vacu_GroundPlane,50.0* CLHEP::um+Vacu_GroundPlane,dp_resonatorAssemblyBaseNbDimZ_1,0.0*M_PI*rad,0.5*M_PI*rad ); // segment angle
G4LogicalVolume* Half_Cirlce_R_SquareCapacitorLogical_GP = new G4LogicalVolume( Half_Cirlce_R_SquareCapacitor_GP,fLiquidHelium,"Half_Cirlce_R_SquareCapacitorLogical_GP");


//////////////////The vertical connector for Qubit 1 and 3
G4VSolid* Square_to_CapcitorV1 = new G4Box("Square_to_CapcitorV1 ",dp_Square_to_CapcitorQXV2 ,dp_Square_to_CapcitorQYV2[0] ,dp_resonatorAssemblyBaseNbDimZ_1);
G4LogicalVolume* Square_to_CapcitorV1Logical = new G4LogicalVolume(Square_to_CapcitorV1,fNiobium,"Square_to_CapcitorV1Logical");
G4VSolid* Square_to_CapcitorV1_GP = new G4Box("Square_to_CapcitorV1_GP",dp_Square_to_CapcitorQXV2+Vacu_GroundPlane ,dp_Square_to_CapcitorQYV2[0] ,dp_resonatorAssemblyBaseNbDimZ_1);
G4LogicalVolume* Square_to_CapcitorV1Logical_GP = new G4LogicalVolume(Square_to_CapcitorV1_GP,fLiquidHelium,"Square_to_CapcitorV1Logical_GP");

G4VSolid* Square_to_CapcitorV2 = new G4Box("Square_to_CapcitorV2 ",dp_Square_to_CapcitorQXV2 ,dp_Square_to_CapcitorQYV2[1] ,dp_resonatorAssemblyBaseNbDimZ_1);
G4LogicalVolume* Square_to_CapcitorV2Logical = new G4LogicalVolume(Square_to_CapcitorV2,fNiobium,"Square_to_CapcitorV2Logical");
G4VSolid* Square_to_CapcitorV2_GP = new G4Box("Square_to_CapcitorV2_GP ",dp_Square_to_CapcitorQXV2+Vacu_GroundPlane ,dp_Square_to_CapcitorQYV2[1] ,dp_resonatorAssemblyBaseNbDimZ_1);
G4LogicalVolume* Square_to_CapcitorV2Logical_GP = new G4LogicalVolume(Square_to_CapcitorV2_GP,fLiquidHelium,"Square_to_CapcitorV2Logical_GP");
/////////////////////////////


////////////////////
G4VSolid* Square_to_CapcitorVB1 = new G4Box("Square_to_CapcitorVB1 ",dp_Square_to_CapcitorQBX2 ,dp_Square_to_CapcitorQBY2[0] ,dp_resonatorAssemblyBaseNbDimZ_1);
G4LogicalVolume* Square_to_CapcitorVB1Logical = new G4LogicalVolume(Square_to_CapcitorVB1,fNiobium,"Square_to_CapcitorVB1Logical");
G4VSolid* Square_to_CapcitorVB1_GP = new G4Box("Square_to_CapcitorVB1_GP ",dp_Square_to_CapcitorQBX2+Vacu_GroundPlane ,dp_Square_to_CapcitorQBY2[0] ,dp_resonatorAssemblyBaseNbDimZ_1);
G4LogicalVolume* Square_to_CapcitorVB1Logical_GP = new G4LogicalVolume(Square_to_CapcitorVB1_GP,fLiquidHelium,"Square_to_CapcitorVB1Logical_GP");

G4VSolid* Square_to_CapcitorVB2 = new G4Box("Square_to_CapcitorVB2",dp_Square_to_CapcitorQBX2 ,dp_Square_to_CapcitorQBY2[1] ,dp_resonatorAssemblyBaseNbDimZ_1);
G4LogicalVolume* Square_to_CapcitorVB2Logical = new G4LogicalVolume(Square_to_CapcitorVB2,fNiobium,"Square_to_CapcitorVB2Logical");
G4VSolid* Square_to_CapcitorVB2_GP = new G4Box("Square_to_CapcitorVB2_GP",dp_Square_to_CapcitorQBX2+Vacu_GroundPlane ,dp_Square_to_CapcitorQBY2[1] ,dp_resonatorAssemblyBaseNbDimZ_1);
G4LogicalVolume* Square_to_CapcitorVB2Logical_GP = new G4LogicalVolume(Square_to_CapcitorVB2_GP,fLiquidHelium,"Square_to_CapcitorVB2Logical_GP");

G4VSolid* Square_to_CapcitorVB3 = new G4Box("Square_to_CapcitorVB3",dp_Square_to_CapcitorQBX2 ,dp_Square_to_CapcitorQBY2[2] ,dp_resonatorAssemblyBaseNbDimZ_1);
G4LogicalVolume* Square_to_CapcitorVB3Logical = new G4LogicalVolume(Square_to_CapcitorVB3,fNiobium,"Square_to_CapcitorVB3Logical");
G4VSolid* Square_to_CapcitorVB3_GP = new G4Box("Square_to_CapcitorVB3_GP",dp_Square_to_CapcitorQBX2+Vacu_GroundPlane ,dp_Square_to_CapcitorQBY2[2] ,dp_resonatorAssemblyBaseNbDimZ_1);
G4LogicalVolume* Square_to_CapcitorVB3Logical_GP = new G4LogicalVolume(Square_to_CapcitorVB3_GP,fLiquidHelium,"Square_to_CapcitorVB3Logical_GP");

G4VSolid* Square_to_CapcitorVB4 = new G4Box("Square_to_CapcitorVB4 ",dp_Square_to_CapcitorQBX2 ,dp_Square_to_CapcitorQBY2[3] ,dp_resonatorAssemblyBaseNbDimZ_1);
G4LogicalVolume* Square_to_CapcitorVB4Logical = new G4LogicalVolume(Square_to_CapcitorVB4,fNiobium,"Square_to_CapcitorVB4Logical");
G4VSolid* Square_to_CapcitorVB4_GP = new G4Box("Square_to_CapcitorVB4_GP_GP ",dp_Square_to_CapcitorQBX2+Vacu_GroundPlane ,dp_Square_to_CapcitorQBY2[3] ,dp_resonatorAssemblyBaseNbDimZ_1);
G4LogicalVolume* Square_to_CapcitorVB4Logical_GP = new G4LogicalVolume(Square_to_CapcitorVB4_GP,fLiquidHelium,"Square_to_CapcitorVB4Logical_GP");
//////////////////////
/// Coordiantes Variables

///////The 4 firts Vertiacals lines
G4VPhysicalVolume* Square_to_CapcitorH1Phys_GP=new G4PVPlacement(0,G4ThreeVector(X_Center_1-dp_Square_to_CapcitorQX1[0]-84*CLHEP::um,Y_Center_1-80*CLHEP::um,0),Square_to_CapcitorH1Logical_GP,"Square_to_CapcitorH1Phys_GP",GroundPlaneLogical,false,0,true);
G4VPhysicalVolume* Square_to_CapcitorH1Phys=new G4PVPlacement(0,G4ThreeVector(),Square_to_CapcitorH1Logical,"Square_to_CapcitorH1Phys",Square_to_CapcitorH1Logical_GP,false,0,true);

G4VPhysicalVolume* Square_to_CapcitorH2Phys_GP=new G4PVPlacement(0,G4ThreeVector(X_Center_2-dp_Square_to_CapcitorQX1[1]-94*CLHEP::um,Y_Center_2-80*CLHEP::um,0),Square_to_CapcitorH2Logical_GP,"Square_to_CapcitorH2Phys_GP",
GroundPlaneLogical,false,0,true);
G4VPhysicalVolume* Square_to_CapcitorH2Phys=new G4PVPlacement(0,G4ThreeVector(),Square_to_CapcitorH2Logical,"Square_to_CapcitorH2Phys",
Square_to_CapcitorH2Logical_GP,false,0,true);

G4VPhysicalVolume* Square_to_CapcitorH3Phys_GP=new G4PVPlacement(0,G4ThreeVector(X_Center_3-dp_Square_to_CapcitorQX1[2]-83*CLHEP::um,Y_Center_3-80*CLHEP::um,0),Square_to_CapcitorH3Logical_GP,"Square_to_CapcitorH3Phys_GP",GroundPlaneLogical,false,0,true);
G4VPhysicalVolume* Square_to_CapcitorH3Phys=new G4PVPlacement(0,G4ThreeVector(),Square_to_CapcitorH3Logical,"Square_to_CapcitorH3Phys",Square_to_CapcitorH3Logical_GP,false,0,true);

G4VPhysicalVolume* Square_to_CapcitorH4Phys_GP=new G4PVPlacement(0,G4ThreeVector(X_Center_4+dp_Square_to_CapcitorQX1[3]+83*CLHEP::um,Y_Center_4-80*CLHEP::um,0),Square_to_CapcitorH4Logical_GP,"Square_to_CapcitorH4Phys_GP",GroundPlaneLogical,false,0,true);
G4VPhysicalVolume* Square_to_CapcitorH4Phys=new G4PVPlacement(0,G4ThreeVector(),Square_to_CapcitorH4Logical,"Square_to_CapcitorH4Phys",Square_to_CapcitorH4Logical_GP,false,0,true);

constexpr double X_Center_5[4]={X_Center_1-0.5*dp_Square_to_CapcitorQX1[0]+84*CLHEP::um,X_Center_2-0.5*dp_Square_to_CapcitorQX1[1]+94*CLHEP::um,X_Center_3-0.5*dp_Square_to_CapcitorQX1[2]+83*CLHEP::um,X_Center_4+0.5*dp_Square_to_CapcitorQX1[3]-83*CLHEP::um};
constexpr double Y_Center_5[4]={Y_Center_1-120*CLHEP::um,Y_Center_2-120*CLHEP::um,Y_Center_3-120*CLHEP::um,Y_Center_4-120*CLHEP::um};
///Connecting the Half Circles Rights
G4VPhysicalVolume* Half_Cirlce_R_SquareCapacitorPhys1_GP=new G4PVPlacement(0,G4ThreeVector(X_Center_5[0]-28*CLHEP::um,Y_Center_5[0],0),Half_Cirlce_R_SquareCapacitorLogical_GP,"Half_Cirlce_R_SquareCapacitorPhys1_GP",GroundPlaneLogical,false,0,true);
G4VPhysicalVolume* Half_Cirlce_R_SquareCapacitorPhys2_GP=new G4PVPlacement(0,G4ThreeVector(X_Center_5[1]-38*CLHEP::um,Y_Center_5[1],0),Half_Cirlce_R_SquareCapacitorLogical_GP,"Half_Cirlce_R_SquareCapacitorPhys2_GP",GroundPlaneLogical,false,0,true);
G4VPhysicalVolume* Half_Cirlce_R_SquareCapacitorPhys3_GP=new G4PVPlacement(0,G4ThreeVector(X_Center_5[2]-1*CLHEP::um,Y_Center_5[2],0),Half_Cirlce_R_SquareCapacitorLogical_GP,"Half_Cirlce_R_SquareCapacitorPhys3_GP",GroundPlaneLogical,false,0,true);
G4VPhysicalVolume* Half_Cirlce_L_SquareCapacitorPhys4_GP=new G4PVPlacement(0,G4ThreeVector(X_Center_5[3]+1*CLHEP::um,Y_Center_5[3],0),Half_Cirlce_L_SquareCapacitorLogical_GP,"Half_Cirlce_L_SquareCapacitorPhys4_GP",GroundPlaneLogical,false,0,true);

G4VPhysicalVolume* Half_Cirlce_R_SquareCapacitorPhys1=new G4PVPlacement(0,G4ThreeVector(),Half_Cirlce_R_SquareCapacitorLogical,"Half_Cirlce_R_SquareCapacitorPhys1",Half_Cirlce_R_SquareCapacitorLogical_GP,false,0,true);
G4VPhysicalVolume* Half_Cirlce_R_SquareCapacitorPhys2=new G4PVPlacement(0,G4ThreeVector(),Half_Cirlce_R_SquareCapacitorLogical,"Half_Cirlce_R_SquareCapacitorPhys2",Half_Cirlce_R_SquareCapacitorLogical_GP,false,0,true);
G4VPhysicalVolume* Half_Cirlce_R_SquareCapacitorPhys3=new G4PVPlacement(0,G4ThreeVector(),Half_Cirlce_R_SquareCapacitorLogical,"Half_Cirlce_R_SquareCapacitorPhys3",Half_Cirlce_R_SquareCapacitorLogical_GP,false,0,true);
G4VPhysicalVolume* Half_Cirlce_L_SquareCapacitorPhys4=new G4PVPlacement(0,G4ThreeVector(),Half_Cirlce_L_SquareCapacitorLogical,"Half_Cirlce_L_SquareCapacitorPhys4",Half_Cirlce_L_SquareCapacitorLogical_GP,false,0,true);


//////////////////////----------------------------------////////////////////
///// New Vertical Connectors and Quarters Circles ------------//////////////
G4VPhysicalVolume* Square_to_CapcitorV1Phys1_GP=new G4PVPlacement(0,G4ThreeVector(X_Center_5[0]-28*CLHEP::um+40*CLHEP::um,Y_Center_5[0]-1*dp_Square_to_CapcitorQYV2[0] ,0),Square_to_CapcitorV1Logical_GP,"Square_to_CapcitorV1Phys1_GP",GroundPlaneLogical,false,0,true);
G4VPhysicalVolume* Half_Cirlce_L_TTDPhys_1_GP=new G4PVPlacement(0,G4ThreeVector(X_Center_5[1]-28*CLHEP::um+70*CLHEP::um,Y_Center_5[1],0),Half_Cirlce_L_TTDLogical,"Half_Cirlce_L_TTDPhys_1_GP",GroundPlaneLogical,false,0,true);
G4VPhysicalVolume* Square_to_CapcitorV1Phys2_GP=new G4PVPlacement(0,G4ThreeVector(X_Center_5[2]-1*CLHEP::um+40*CLHEP::um,Y_Center_5[2]-1*dp_Square_to_CapcitorQYV2[1] ,0),Square_to_CapcitorV2Logical_GP,"Square_to_CapcitorV1Phys2_GP",GroundPlaneLogical,false,0,true);

G4VPhysicalVolume* Square_to_CapcitorV1Phys1=new G4PVPlacement(0,G4ThreeVector(),Square_to_CapcitorV1Logical,"Square_to_CapcitorV1Phys1",Square_to_CapcitorV1Logical_GP,false,0,true);
G4VPhysicalVolume* Half_Cirlce_L_TTDPhys_1=new G4PVPlacement(0,G4ThreeVector(),Half_Cirlce_L_TTDLogical,"Half_Cirlce_L_TTDPhys_1",Half_Cirlce_L_TTDLogical_GP,false,0,true);
G4VPhysicalVolume* Square_to_CapcitorV1Phys2=new G4PVPlacement(0,G4ThreeVector(),Square_to_CapcitorV2Logical,"Square_to_CapcitorV1Phys2",Square_to_CapcitorV2Logical_GP,false,0,true);


//Those are the circles connectors for the final line with qubit

G4VPhysicalVolume* Half_Cirlce_R_TTDPhys_1_GP=new G4PVPlacement(0,G4ThreeVector(X_Center_5[0]-28*CLHEP::um+80*CLHEP::um,Y_Center_5[0]-2*dp_Square_to_CapcitorQYV2[0] ,0),Half_Cirlce_L_TTDLogical_GP,"Half_Cirlce_R_TTDPhys_1_GP",GroundPlaneLogical,false,0,true);
G4VPhysicalVolume* Half_Cirlce_R_TTDPhys_2_GP=new G4PVPlacement(0,G4ThreeVector(X_Center_5[2]-28*CLHEP::um+107*CLHEP::um,Y_Center_5[2]-2*dp_Square_to_CapcitorQYV2[1],0),Half_Cirlce_L_TTDLogical_GP,"Half_Cirlce_R_TTDPhys_2_GP",GroundPlaneLogical,false,0,true);
G4VPhysicalVolume* Half_Cirlce_R_TTDPhys_3_GP=new G4PVPlacement(0,G4ThreeVector(X_Center_5[3]-28*CLHEP::um-50*CLHEP::um,Y_Center_5[3],0),Half_Cirlce_R_TTDLogical_GP,"Half_Cirlce_R_TTDPhys_3_GP",GroundPlaneLogical,false,0,true);

G4VPhysicalVolume* Half_Cirlce_R_TTDPhys_1=new G4PVPlacement(0,G4ThreeVector(),Half_Cirlce_L_TTDLogical,"Half_Cirlce_R_TTDPhys_1",Half_Cirlce_L_TTDLogical_GP,false,0,true);
G4VPhysicalVolume* Half_Cirlce_R_TTDPhys_2=new G4PVPlacement(0,G4ThreeVector(),Half_Cirlce_L_TTDLogical,"Half_Cirlce_R_TTDPhys_2",Half_Cirlce_L_TTDLogical_GP,false,0,true);
G4VPhysicalVolume* Half_Cirlce_R_TTDPhys_3=new G4PVPlacement(0,G4ThreeVector(),Half_Cirlce_R_TTDLogical,"Half_Cirlce_R_TTDPhys_3",Half_Cirlce_R_TTDLogical_GP,false,0,true);

/////orizontal Lines ///////////////////

/////////////////////////////-------------------------------CONTINUIM HERE------------------

/////////Horizontal Line for Qubits 1 ////////////
G4VSolid* Square_to_Qubit1 = new G4Box("Square_to_Qubit1",dp_Square_to_CapcitorQX2[0] ,dp_Square_to_CapcitorQY1 ,dp_resonatorAssemblyBaseNbDimZ_1);
G4LogicalVolume* Square_to_Qubit1Logical = new G4LogicalVolume( Square_to_Qubit1,fNiobium,"Square_to_Qubit1Logical");
G4VSolid* Square_to_Qubit1_GP = new G4Box("Square_to_Qubit1_GP",dp_Square_to_CapcitorQX2[0] ,dp_Square_to_CapcitorQY1+Vacu_GroundPlane ,dp_resonatorAssemblyBaseNbDimZ_1);
G4LogicalVolume* Square_to_Qubit1Logical_GP = new G4LogicalVolume( Square_to_Qubit1_GP,fLiquidHelium,"Square_to_Qubit1Logical_GP");
constexpr double X_Center_6=X_Center_5[0]+52*CLHEP::um;
constexpr double Y_Center_6=Y_Center_5[0]-2*dp_Square_to_CapcitorQYV2[0]-40*CLHEP::um;

G4VPhysicalVolume* Square_to_Qubit1Phys_GP=new G4PVPlacement(0,G4ThreeVector(X_Center_6+1.0*dp_Square_to_CapcitorQX2[0],Y_Center_6,0),Square_to_Qubit1Logical_GP,"Square_to_Capcitor1Phys_GP",GroundPlaneLogical,false,0,true);
G4VPhysicalVolume* Square_to_Qubit1Phys=new G4PVPlacement(0,G4ThreeVector(),Square_to_Qubit1Logical,"Square_to_Capcitor1Phys",Square_to_Qubit1Logical_GP,false,0,true);

///////Horizontal line to qubit 2
G4VSolid* Square_to_Qubit2 = new G4Box("Square_to_Qubit2",dp_Square_to_CapcitorQX2[1] ,dp_Square_to_CapcitorQY1 ,dp_resonatorAssemblyBaseNbDimZ_1);
G4LogicalVolume* Square_to_Qubit2Logical = new G4LogicalVolume( Square_to_Qubit2,fNiobium,"Square_to_Qubit2Logical");
G4VSolid* Square_to_Qubit2_GP = new G4Box("Square_to_Qubit2_GP",dp_Square_to_CapcitorQX2[1] ,dp_Square_to_CapcitorQY1+Vacu_GroundPlane ,dp_resonatorAssemblyBaseNbDimZ_1);
G4LogicalVolume* Square_to_Qubit2Logical_GP = new G4LogicalVolume( Square_to_Qubit2_GP,fLiquidHelium,"Square_to_Qubit2Logical_GP");
constexpr double X_Center_7=X_Center_5[1]+42*CLHEP::um;
constexpr double Y_Center_7=Y_Center_5[1]-40*CLHEP::um;
G4VPhysicalVolume* Square_to_Qubit2Phys_GP=new G4PVPlacement(0,G4ThreeVector(X_Center_7+1.0*dp_Square_to_CapcitorQX2[1],Y_Center_7,0),Square_to_Qubit2Logical_GP,"Square_to_Capcitor2Phys_GP",GroundPlaneLogical,false,0,true);
G4VPhysicalVolume* Square_to_Qubit2Phys=new G4PVPlacement(0,G4ThreeVector(),Square_to_Qubit2Logical,"Square_to_Capcitor2Phys",Square_to_Qubit2Logical_GP,false,0,true);

///////Horizontal line to qubit 3
G4VSolid* Square_to_Qubit3 = new G4Box("Square_to_Qubit3",dp_Square_to_CapcitorQX2[2] ,dp_Square_to_CapcitorQY1 ,dp_resonatorAssemblyBaseNbDimZ_1);
G4LogicalVolume* Square_to_Qubit3Logical = new G4LogicalVolume( Square_to_Qubit3,fNiobium,"Square_to_Qubit3Logical");
G4VSolid* Square_to_Qubit3_GP = new G4Box("Square_to_Qubit3_GP",dp_Square_to_CapcitorQX2[2] ,dp_Square_to_CapcitorQY1+Vacu_GroundPlane ,dp_resonatorAssemblyBaseNbDimZ_1);
G4LogicalVolume* Square_to_Qubit3Logical_GP = new G4LogicalVolume( Square_to_Qubit3_GP,fLiquidHelium,"Square_to_Qubit3Logical_GP");
constexpr double X_Center_8=X_Center_5[2]+79*CLHEP::um;
constexpr double Y_Center_8=Y_Center_5[2]-40*CLHEP::um-2*dp_Square_to_CapcitorQYV2[1];
G4VPhysicalVolume* Square_to_Qubit3Phys_GP=new G4PVPlacement(0,G4ThreeVector(X_Center_8+1.0*dp_Square_to_CapcitorQX2[2],Y_Center_8,0),Square_to_Qubit3Logical_GP,"Square_to_CapcitorH3Phys_GP",GroundPlaneLogical,false,0,true);
G4VPhysicalVolume* Square_to_Qubit3Phys=new G4PVPlacement(0,G4ThreeVector(),Square_to_Qubit3Logical,"Square_to_Capcitor3Phys",Square_to_Qubit3Logical_GP,false,0,true);

///////Horizontal line to qubit 4
G4VSolid* Square_to_Qubit4 = new G4Box("Square_to_Qubit3",dp_Square_to_CapcitorQX2[3] ,dp_Square_to_CapcitorQY1 ,dp_resonatorAssemblyBaseNbDimZ_1);
G4LogicalVolume* Square_to_Qubit4Logical = new G4LogicalVolume( Square_to_Qubit4,fNiobium,"Square_to_Qubit4Logical");
G4VSolid* Square_to_Qubit4_GP = new G4Box("Square_to_Qubit3_GP",dp_Square_to_CapcitorQX2[3] ,dp_Square_to_CapcitorQY1+Vacu_GroundPlane ,dp_resonatorAssemblyBaseNbDimZ_1);
G4LogicalVolume* Square_to_Qubit4Logical_GP = new G4LogicalVolume( Square_to_Qubit4_GP,fLiquidHelium,"Square_to_Qubit4Logical_GP");
constexpr double X_Center_9=X_Center_5[3]-79*CLHEP::um;
constexpr double Y_Center_9=Y_Center_5[3]-40*CLHEP::um;
G4VPhysicalVolume* Square_to_Qubit4Phys_GP=new G4PVPlacement(0,G4ThreeVector(X_Center_9-1.0*dp_Square_to_CapcitorQX2[3],Y_Center_9,0),Square_to_Qubit4Logical_GP,"Square_to_Capcitor4Phys_GP",GroundPlaneLogical,false,0,true);
G4VPhysicalVolume* Square_to_Qubit4Phys=new G4PVPlacement(0,G4ThreeVector(),Square_to_Qubit4Logical,"Square_to_Capcitor4Phys",Square_to_Qubit4Logical_GP,false,0,true);

///////////////////////////


//////////////////////////////////---------------Conencting the Quaters Circles 3 Right and one Left//////////////

G4VPhysicalVolume* Half_Cirlce_R_TTDPhys_Q_1_GP=new G4PVPlacement(0,G4ThreeVector(X_Center_6+2.0*dp_Square_to_CapcitorQX2[0],Y_Center_6-40*CLHEP::um,0),Half_Cirlce_R_SquareCapacitorLogical_GP,"Half_Cirlce_R_TTDPhys_Q_1_GP",GroundPlaneLogical,false,0,true);
G4VPhysicalVolume* Half_Cirlce_R_TTDPhys_Q_2_GP=new G4PVPlacement(0,G4ThreeVector(X_Center_7+2.0*dp_Square_to_CapcitorQX2[1],Y_Center_7-40*CLHEP::um,0),Half_Cirlce_R_SquareCapacitorLogical_GP,"Half_Cirlce_R_TTDPhys_Q_2_GP",GroundPlaneLogical,false,0,true);
G4VPhysicalVolume* Half_Cirlce_R_TTDPhys_Q_3_GP=new G4PVPlacement(0,G4ThreeVector(X_Center_8+2.0*dp_Square_to_CapcitorQX2[2],Y_Center_8-40*CLHEP::um,0),Half_Cirlce_R_SquareCapacitorLogical_GP,"Half_Cirlce_R_TTDPhys_Q_3_GP",GroundPlaneLogical,false,0,true);
G4VPhysicalVolume* Half_Cirlce_R_TTDPhys_Q_4_GP=new G4PVPlacement(0,G4ThreeVector(X_Center_9-2.0*dp_Square_to_CapcitorQX2[3],Y_Center_9-40*CLHEP::um,0),Half_Cirlce_L_SquareCapacitorLogical_GP,"Half_Cirlce_R_TTDPhys_Q_4_GP",GroundPlaneLogical,false,0,true);


G4VPhysicalVolume* Half_Cirlce_R_TTDPhys_Q_1=new G4PVPlacement(0,G4ThreeVector(),Half_Cirlce_R_SquareCapacitorLogical,"Half_Cirlce_R_TTDPhys_Q_1",Half_Cirlce_R_SquareCapacitorLogical_GP,false,0,true);
G4VPhysicalVolume* Half_Cirlce_R_TTDPhys_Q_2=new G4PVPlacement(0,G4ThreeVector(),Half_Cirlce_R_SquareCapacitorLogical,"Half_Cirlce_R_TTDPhys_Q_2",Half_Cirlce_R_SquareCapacitorLogical_GP,false,0,true);
G4VPhysicalVolume* Half_Cirlce_R_TTDPhys_Q_3=new G4PVPlacement(0,G4ThreeVector(),Half_Cirlce_R_SquareCapacitorLogical,"Half_Cirlce_R_TTDPhys_Q_3",Half_Cirlce_R_SquareCapacitorLogical_GP,false,0,true);
G4VPhysicalVolume* Half_Cirlce_R_TTDPhys_Q_4=new G4PVPlacement(0,G4ThreeVector(),Half_Cirlce_L_SquareCapacitorLogical,"Half_Cirlce_R_TTDPhys_Q_4",Half_Cirlce_L_SquareCapacitorLogical_GP,false,0,true);


//////////////////////----------------------------

///////////////--------___Finals Vertical Lines Connected to Qubits
//constexpr double dp_Square_to_CapcitorQBY2[4]= {315.0*CLHEP::um,610.0*CLHEP::um,1092.0*CLHEP::um,610.0*CLHEP::um};
G4VSolid* SquareV_to_Qubit1 = new G4Box("SquareV_to_Qubit1",dp_Square_to_CapcitorQY1,dp_Square_to_CapcitorQBY2[0] ,dp_resonatorAssemblyBaseNbDimZ_1);
G4LogicalVolume* SquareV_to_Qubit1Logical = new G4LogicalVolume( SquareV_to_Qubit1,fNiobium,"SquareV_to_Qubit1Logical");
G4VSolid* SquareV_to_Qubit1_GP = new G4Box("SquareV_to_Qubit1_GP",dp_Square_to_CapcitorQY1+Vacu_GroundPlane,dp_Square_to_CapcitorQBY2[0] ,dp_resonatorAssemblyBaseNbDimZ_1);
G4LogicalVolume* SquareV_to_Qubit1Logical_GP = new G4LogicalVolume( SquareV_to_Qubit1_GP,fLiquidHelium,"SquareV_to_Qubit1Logical_GP");
constexpr double X_Center_10=X_Center_6+2.0*dp_Square_to_CapcitorQX2[0]+40*CLHEP::um;
constexpr double Y_Center_10=Y_Center_6-40*CLHEP::um-dp_Square_to_CapcitorQBY2[0];
G4VPhysicalVolume* SquareV_to_Qubit1Phys_GP=new G4PVPlacement(0,G4ThreeVector(X_Center_10,Y_Center_10,0),SquareV_to_Qubit1Logical_GP,"SquareV_to_Qubit1Phys_GP",GroundPlaneLogical,false,0,true);
G4VPhysicalVolume* SquareV_to_Qubit1Phys=new G4PVPlacement(0,G4ThreeVector(),SquareV_to_Qubit1Logical,"SquareV_to_Qubit1Phys",SquareV_to_Qubit1Logical_GP,false,0,true);

///////////////////Verticsl Line to Qubit 2
G4VSolid* SquareV_to_Qubit2 = new G4Box("SquareV_to_Qubit2",dp_Square_to_CapcitorQY1,dp_Square_to_CapcitorQBY2[1] ,dp_resonatorAssemblyBaseNbDimZ_1);
G4LogicalVolume* SquareV_to_Qubit2Logical = new G4LogicalVolume( SquareV_to_Qubit2,fNiobium,"SquareV_to_Qubit2Logical");
G4VSolid* SquareV_to_Qubit2_GP = new G4Box("SquareV_to_Qubit2_GP",dp_Square_to_CapcitorQY1+Vacu_GroundPlane,dp_Square_to_CapcitorQBY2[1] ,dp_resonatorAssemblyBaseNbDimZ_1);
G4LogicalVolume* SquareV_to_Qubit2Logical_GP = new G4LogicalVolume( SquareV_to_Qubit2_GP,fLiquidHelium,"SquareV_to_Qubit2Logical_GP");
constexpr double X_Center_11=X_Center_7+2.0*dp_Square_to_CapcitorQX2[1]+40*CLHEP::um;
constexpr double Y_Center_11=Y_Center_7-40*CLHEP::um-dp_Square_to_CapcitorQBY2[1];
G4VPhysicalVolume* SquareV_to_Qubit2Phys_GP=new G4PVPlacement(0,G4ThreeVector(X_Center_11,Y_Center_11,0),SquareV_to_Qubit2Logical_GP,"SquareV_to_Qubit2Phys_GP",GroundPlaneLogical,false,0,true);
G4VPhysicalVolume* SquareV_to_Qubit2Phys=new G4PVPlacement(0,G4ThreeVector(),SquareV_to_Qubit2Logical,"SquareV_to_Qubit2Phys",SquareV_to_Qubit2Logical_GP,false,0,true);

///////////////////Verticsl Line to Qubit 3
G4VSolid* SquareV_to_Qubit3= new G4Box("SquareV_to_Qubit3",dp_Square_to_CapcitorQY1,dp_Square_to_CapcitorQBY2[2] ,dp_resonatorAssemblyBaseNbDimZ_1);
G4LogicalVolume* SquareV_to_Qubit3Logical = new G4LogicalVolume( SquareV_to_Qubit3,fNiobium,"SquareV_to_Qubit3Logical");
G4VSolid* SquareV_to_Qubit3_GP= new G4Box("SquareV_to_Qubit3_GP",dp_Square_to_CapcitorQY1+Vacu_GroundPlane,dp_Square_to_CapcitorQBY2[2] ,dp_resonatorAssemblyBaseNbDimZ_1);
G4LogicalVolume* SquareV_to_Qubit3Logical_GP = new G4LogicalVolume( SquareV_to_Qubit3_GP,fLiquidHelium,"SquareV_to_Qubit3Logical_GP");
constexpr double X_Center_12=X_Center_8+2.0*dp_Square_to_CapcitorQX2[2]+40*CLHEP::um;
constexpr double Y_Center_12=Y_Center_8-40*CLHEP::um-dp_Square_to_CapcitorQBY2[2];
G4VPhysicalVolume* SquareV_to_Qubit3Phys_GP=new G4PVPlacement(0,G4ThreeVector(X_Center_12,Y_Center_12,0),SquareV_to_Qubit3Logical_GP,"SquareV_to_Qubit3Phys_GP",GroundPlaneLogical,false,0,true);
G4VPhysicalVolume* SquareV_to_Qubit3Phys=new G4PVPlacement(0,G4ThreeVector(),SquareV_to_Qubit3Logical,"SquareV_to_Qubit3Phys",SquareV_to_Qubit3Logical_GP,false,0,true);

///////////////////Vertical Line to Qubit 4
G4VSolid* SquareV_to_Qubit4= new G4Box("SquareV_to_Qubit4",dp_Square_to_CapcitorQY1,dp_Square_to_CapcitorQBY2[3] ,dp_resonatorAssemblyBaseNbDimZ_1);
G4LogicalVolume* SquareV_to_Qubit4Logical = new G4LogicalVolume( SquareV_to_Qubit4,fNiobium,"SquareV_to_Qubit4Logical");
G4VSolid* SquareV_to_Qubit4_GP= new G4Box("SquareV_to_Qubit4_GP",dp_Square_to_CapcitorQY1+Vacu_GroundPlane,dp_Square_to_CapcitorQBY2[3] ,dp_resonatorAssemblyBaseNbDimZ_1);
G4LogicalVolume* SquareV_to_Qubit4Logical_GP = new G4LogicalVolume( SquareV_to_Qubit4_GP,fLiquidHelium,"SquareV_to_Qubit4Logical_GP");
constexpr double X_Center_13=X_Center_9-2.0*dp_Square_to_CapcitorQX2[3]-40*CLHEP::um;
constexpr double Y_Center_13=Y_Center_9-40*CLHEP::um-dp_Square_to_CapcitorQBY2[3];
G4VPhysicalVolume* SquareV_to_Qubit4Phys_GP=new G4PVPlacement(0,G4ThreeVector(X_Center_13,Y_Center_13,0),SquareV_to_Qubit4Logical_GP,"SquareV_to_Qubit4Phys_GP",GroundPlaneLogical,false,0,true);
G4VPhysicalVolume* SquareV_to_Qubit4Phys=new G4PVPlacement(0,G4ThreeVector(),SquareV_to_Qubit4Logical,"SquareV_to_Qubit4Phys",SquareV_to_Qubit4Logical_GP,false,0,true);

/////////////////////////------------


//###############################New Section Qubits Square Connections  Cutted with Small Square
//120 microns on X and 145 on Y
//Distance 20 micorn s
G4VSolid* Square_Capacitor_B= new G4Box("Square_Capacitor_B",120*CLHEP::um,145*CLHEP::um,dp_resonatorAssemblyBaseNbDimZ_1);
G4VSolid* Square_Capacitor_S= new G4Box("Square_Capacitor_S",80*CLHEP::um,145*CLHEP::um,dp_resonatorAssemblyBaseNbDimZ_1);
G4VSolid* Square_Capacitor_B_GP= new G4Box("Square_Capacitor_B_GP",120*CLHEP::um+Vacu_GroundPlane,145*CLHEP::um+0*Vacu_GroundPlane,dp_resonatorAssemblyBaseNbDimZ_1);
G4VSolid* Square_Capacitor_S_GP= new G4Box("Square_Capacitor_S_GP",80*CLHEP::um,145*CLHEP::um,dp_resonatorAssemblyBaseNbDimZ_1);

G4VSolid* Square_Capacitor_Qubit= new G4SubtractionSolid("Square_Capacitor_B-Square_Capacitor_S",Square_Capacitor_B, Square_Capacitor_S,0, G4ThreeVector(0,-40*CLHEP::um,0.));
G4VSolid* Square_Capacitor_Qubit_GP= new G4SubtractionSolid("Square_Capacitor_B_GP-Square_Capacitor_S_GP",Square_Capacitor_B_GP, Square_Capacitor_S_GP,0, G4ThreeVector(0,-40*CLHEP::um,0.));

//////New 3 Squares to Fill the space bewtten the qubits
G4VSolid* Square_Capacitor_Fill_1= new G4Box("Square_Capacitor_Fill_1",70*CLHEP::um,3*CLHEP::um,dp_resonatorAssemblyBaseNbDimZ_1);
G4VSolid* Square_Capacitor_Fill_2= new G4Box("Square_Capacitor_Fill_2",3*CLHEP::um,123.5*CLHEP::um,dp_resonatorAssemblyBaseNbDimZ_1);
G4LogicalVolume* Square_Capacitor_QubitLogical = new G4LogicalVolume( Square_Capacitor_Qubit,fNiobium,"Square_Capacitor_QubitLogical");

G4LogicalVolume* Square_Capacitor_QubitLogical_GP1 = new G4LogicalVolume( Square_Capacitor_Qubit_GP,fLiquidHelium,"Square_Capacitor_QubitLogical_GP");
G4LogicalVolume* Square_Capacitor_QubitLogical_GP2 = new G4LogicalVolume( Square_Capacitor_Qubit_GP,fLiquidHelium,"Square_Capacitor_QubitLogical_GP2");
G4LogicalVolume* Square_Capacitor_QubitLogical_GP3 = new G4LogicalVolume( Square_Capacitor_Qubit_GP,fLiquidHelium,"Square_Capacitor_QubitLogical_GP3");
G4LogicalVolume* Square_Capacitor_QubitLogical_GP4 = new G4LogicalVolume( Square_Capacitor_Qubit_GP,fLiquidHelium,"Square_Capacitor_QubitLogical_GP4");

G4LogicalVolume* Square_Capacitor_Fill_1Logical = new G4LogicalVolume( Square_Capacitor_Fill_1,fNiobium,"Square_Capacitor_Fill_1Logical");
G4LogicalVolume* Square_Capacitor_Fill_2Logical = new G4LogicalVolume( Square_Capacitor_Fill_2,fNiobium,"Square_Capacitor_Fill_1Logical");
constexpr double X_Center_14=X_Center_10;
constexpr double Y_Center_14=Y_Center_10-3*154*CLHEP::um;
G4VPhysicalVolume* Square_Capacitor_QubitPhys1_GP=new G4PVPlacement(0,G4ThreeVector(X_Center_14,Y_Center_14+2*CLHEP::um,0),Square_Capacitor_QubitLogical_GP1,"Square_Capacitor_QubitPhys1_GP",GroundPlaneLogical,false,0,true);
G4VPhysicalVolume* Square_Capacitor_QubitPhys1=new G4PVPlacement(0,G4ThreeVector(),Square_Capacitor_QubitLogical,"Square_Capacitor_QubitPhys1",Square_Capacitor_QubitLogical_GP1,false,0,true);
// G4VPhysicalVolume* Square_Capacitor_Fill_1Phys=new G4PVPlacement(0,G4ThreeVector(0,99*CLHEP::um,0),Square_Capacitor_Fill_1Logical,"Square_Capacitor_Fill_1Phys",Square_Capacitor_QubitLogical_GP1,false,0,true);
// G4VPhysicalVolume* Square_Capacitor_Fill_2Phys=new G4PVPlacement(0,G4ThreeVector(-73*CLHEP::um,+99*CLHEP::um-120.5*CLHEP::um,0),Square_Capacitor_Fill_2Logical,"Square_Capacitor_Fill_2Phys",Square_Capacitor_QubitLogical_GP1,false,0,true);
// G4VPhysicalVolume* Square_Capacitor_Fill_3Phys=new G4PVPlacement(0,G4ThreeVector(+73*CLHEP::um,99*CLHEP::um-120.5*CLHEP::um,0),Square_Capacitor_Fill_2Logical,"Square_Capacitor_Fill_3Phys",Square_Capacitor_QubitLogical_GP1,false,0,true);

G4VPhysicalVolume* Square_Capacitor_QubitPhys2_GP=new G4PVPlacement(0,G4ThreeVector(X_Center_11,Y_Center_14+8*CLHEP::um,0),Square_Capacitor_QubitLogical_GP2,"Square_Capacitor_QubitPhys2_GP",GroundPlaneLogical,false,0,true);
G4VPhysicalVolume* Square_Capacitor_QubitPhys2=new G4PVPlacement(0,G4ThreeVector(),Square_Capacitor_QubitLogical,"Square_Capacitor_QubitPhys2",Square_Capacitor_QubitLogical_GP2,false,0,true);
// G4VPhysicalVolume* Square_Capacitor_Fill_1Phys2=new G4PVPlacement(0,G4ThreeVector(0,8*CLHEP::um+99*CLHEP::um,0),Square_Capacitor_Fill_1Logical,"Square_Capacitor_Fill_1Phys2",Square_Capacitor_QubitLogical_GP2,false,0,true);
// G4VPhysicalVolume* Square_Capacitor_Fill_2Phys2=new G4PVPlacement(0,G4ThreeVector(-73*CLHEP::um,8*CLHEP::um+99*CLHEP::um-120.5*CLHEP::um,0),Square_Capacitor_Fill_2Logical,"Square_Capacitor_Fill_2Phys2",Square_Capacitor_QubitLogical_GP2,false,0,true);
// G4VPhysicalVolume* Square_Capacitor_Fill_3Phys2=new G4PVPlacement(0,G4ThreeVector(+73*CLHEP::um,8*CLHEP::um+99*CLHEP::um-120.5*CLHEP::um,0),Square_Capacitor_Fill_2Logical,"Square_Capacitor_Fill_3Phys2",Square_Capacitor_QubitLogical_GP2,false,0,true);

G4VPhysicalVolume* Square_Capacitor_QubitPhys3_GP=new G4PVPlacement(0,G4ThreeVector(X_Center_12,Y_Center_14+4*CLHEP::um,0),Square_Capacitor_QubitLogical_GP3,"Square_Capacitor_QubitPhys3_GP",GroundPlaneLogical,false,0,true);
G4VPhysicalVolume* Square_Capacitor_QubitPhys3=new G4PVPlacement(0,G4ThreeVector(),Square_Capacitor_QubitLogical,"Square_Capacitor_QubitPhys3",Square_Capacitor_QubitLogical_GP3,false,0,true);
// G4VPhysicalVolume* Square_Capacitor_Fill_1Phys3=new G4PVPlacement(0,G4ThreeVector(0,4*CLHEP::um+99*CLHEP::um,0),Square_Capacitor_Fill_1Logical,"Square_Capacitor_Fill_1Phys3",Square_Capacitor_QubitLogical_GP3,false,0,true);
// G4VPhysicalVolume* Square_Capacitor_Fill_2Phys3=new G4PVPlacement(0,G4ThreeVector(-73*CLHEP::um,4*CLHEP::um+99*CLHEP::um-120.5*CLHEP::um,0),Square_Capacitor_Fill_2Logical,"Square_Capacitor_Fill_2Phys3",Square_Capacitor_QubitLogical_GP3,false,0,true);
// G4VPhysicalVolume* Square_Capacitor_Fill_3Phys3=new G4PVPlacement(0,G4ThreeVector(+73*CLHEP::um,4*CLHEP::um+99*CLHEP::um-120.5*CLHEP::um,0),Square_Capacitor_Fill_2Logical,"Square_Capacitor_Fill_3Phys3",Square_Capacitor_QubitLogical_GP3,false,0,true);

G4VPhysicalVolume* Square_Capacitor_QubitPhys4_GP=new G4PVPlacement(0,G4ThreeVector(X_Center_13,Y_Center_14+8*CLHEP::um,0),Square_Capacitor_QubitLogical_GP4,"Square_Capacitor_QubitPhys4_GP",GroundPlaneLogical,false,0,true);
G4VPhysicalVolume* Square_Capacitor_QubitPhys4=new G4PVPlacement(0,G4ThreeVector(),Square_Capacitor_QubitLogical,"Square_Capacitor_QubitPhys4",Square_Capacitor_QubitLogical_GP4,false,0,true);
// G4VPhysicalVolume* Square_Capacitor_Fill_1Phys4=new G4PVPlacement(0,G4ThreeVector(0,8*CLHEP::um+99*CLHEP::um,0),Square_Capacitor_Fill_1Logical,"Square_Capacitor_Fill_1Phys4",Square_Capacitor_QubitLogical_GP4,false,0,true);
// G4VPhysicalVolume* Square_Capacitor_Fill_2Phys4=new G4PVPlacement(0,G4ThreeVector(-73*CLHEP::um,8*CLHEP::um+99*CLHEP::um-120.5*CLHEP::um,0),Square_Capacitor_Fill_2Logical,"Square_Capacitor_Fill_2Phys4",Square_Capacitor_QubitLogical_GP4,false,0,true);
// G4VPhysicalVolume* Square_Capacitor_Fill_3Phys4=new G4PVPlacement(0,G4ThreeVector(73*CLHEP::um,8*CLHEP::um+99*CLHEP::um-120.5*CLHEP::um,0),Square_Capacitor_Fill_2Logical,"Square_Capacitor_Fill_3Phys4",Square_Capacitor_QubitLogical_GP4,false,0,true);

//////////////////----------------------Inseritn the More important part of the Qubits the Cross Qubits
////// I need to check and rename for individual qubits
constexpr double X_Center_15[4]={X_Center_14,X_Center_11,X_Center_12,X_Center_13};
constexpr double Y_Center_15[4]={Y_Center_14+2*CLHEP::um,Y_Center_14+8*CLHEP::um,Y_Center_14+4*CLHEP::um,Y_Center_14+8*CLHEP::um};
for(G4int i=0;i<4;i++){
  G4VPhysicalVolume* Cross_QubitPhys_GP =new G4PVPlacement(0,G4ThreeVector(X_Center_15[i],Y_Center_15[i]-0.5*dp_CrossQubitNbDimX-2*58*CLHEP::um,0),Cross_QubitLogical_GP,"Cross_Qubit_GP",GroundPlaneLogical,false,i,true);
  G4CMPLogicalBorderSurface * border_sapphire_Vacuum_Cross_QubitPhys_GP=new G4CMPLogicalBorderSurface("border_sapphire_Vacuum_Cross_QubitPhys_GP",SapphirePhys,Cross_QubitPhys_GP,fSiVacuumInterface);

  G4VPhysicalVolume* Cross_QubitPhys =new G4PVPlacement(0,G4ThreeVector(),Cross_QubitLogical,"Cross_Qubit",Cross_QubitLogical_GP,false,i,true);
G4CMPLogicalBorderSurface * border_sapphire_Vacuum_Cross_QubitPhys=new G4CMPLogicalBorderSurface("border_sapphire_Vacuum_Cross_QubitPhys",SapphirePhys,Cross_QubitPhys,fSiNbInterface);


}




//////////////////////////////////////////Finals Steps The Transmission Lines to the Panels-----------**********************
/////////////Creating the First Vertical lines
G4VSolid* Square_Capacitor_Q_To_P_VS1= new G4Box("Square_Capacitor_Q_To_P_VS1",10.0*CLHEP::um,113*CLHEP::um,dp_resonatorAssemblyBaseNbDimZ_1);
G4VSolid* Square_Capacitor_Q_To_P_VS2= new G4Box("Square_Capacitor_Q_To_P_VS2",10*CLHEP::um,495*CLHEP::um,dp_resonatorAssemblyBaseNbDimZ_1);
G4VSolid* Square_Capacitor_Q_To_P_VS3= new G4Box("Square_Capacitor_Q_To_P_VS3",10*CLHEP::um,620*CLHEP::um,dp_resonatorAssemblyBaseNbDimZ_1);
G4VSolid* Square_Capacitor_Q_To_P_VS4= new G4Box("Square_Capacitor_Q_To_P_VS4",10*CLHEP::um,620*CLHEP::um,dp_resonatorAssemblyBaseNbDimZ_1);

G4VSolid* Square_Capacitor_Q_To_P_VS1_GP= new G4Box("Square_Capacitor_Q_To_P_VS1_GP",10.0*CLHEP::um+Vacu_GroundPlane,113*CLHEP::um,dp_resonatorAssemblyBaseNbDimZ_1);
G4VSolid* Square_Capacitor_Q_To_P_VS2_GP= new G4Box("Square_Capacitor_Q_To_P_VS2_GP",10*CLHEP::um+Vacu_GroundPlane,495*CLHEP::um,dp_resonatorAssemblyBaseNbDimZ_1);
G4VSolid* Square_Capacitor_Q_To_P_VS3_GP= new G4Box("Square_Capacitor_Q_To_P_VS3_GP",10*CLHEP::um+Vacu_GroundPlane,620*CLHEP::um,dp_resonatorAssemblyBaseNbDimZ_1);
G4VSolid* Square_Capacitor_Q_To_P_VS4_GP= new G4Box("Square_Capacitor_Q_To_P_VS4_GP",10*CLHEP::um+Vacu_GroundPlane,620*CLHEP::um,dp_resonatorAssemblyBaseNbDimZ_1);

G4LogicalVolume* Square_Capacitor_Q_To_P_VS1Logical = new G4LogicalVolume( Square_Capacitor_Q_To_P_VS1,fNiobium,"Square_Capacitor_Q_To_P_VS1Logical");
G4LogicalVolume* Square_Capacitor_Q_To_P_VS2Logical = new G4LogicalVolume( Square_Capacitor_Q_To_P_VS2,fNiobium,"Square_Capacitor_Q_To_P_VS2Logical");
G4LogicalVolume* Square_Capacitor_Q_To_P_VS3Logical = new G4LogicalVolume( Square_Capacitor_Q_To_P_VS3,fNiobium,"Square_Capacitor_Q_To_P_VS3Logical");
G4LogicalVolume* Square_Capacitor_Q_To_P_VS4Logical = new G4LogicalVolume( Square_Capacitor_Q_To_P_VS4,fNiobium,"Square_Capacitor_Q_To_P_VS4Logical");

G4LogicalVolume* Square_Capacitor_Q_To_P_VS1Logical_GP = new G4LogicalVolume( Square_Capacitor_Q_To_P_VS1_GP,fLiquidHelium,"Square_Capacitor_Q_To_P_VS1Logical_GP");
G4LogicalVolume* Square_Capacitor_Q_To_P_VS2Logical_GP = new G4LogicalVolume( Square_Capacitor_Q_To_P_VS2_GP,fLiquidHelium,"Square_Capacitor_Q_To_P_VS2Logical_GP");
G4LogicalVolume* Square_Capacitor_Q_To_P_VS3Logical_GP = new G4LogicalVolume( Square_Capacitor_Q_To_P_VS3_GP,fLiquidHelium,"Square_Capacitor_Q_To_P_VS3Logical_GP");
G4LogicalVolume* Square_Capacitor_Q_To_P_VS4Logical_GP = new G4LogicalVolume( Square_Capacitor_Q_To_P_VS4_GP,fLiquidHelium,"Square_Capacitor_Q_To_P_VS4Logical_GP");




G4VPhysicalVolume* Square_Capacitor_Q_To_P_VS1Phys_GP=new G4PVPlacement(0,G4ThreeVector(X_Center_15[0],Y_Center_15[0]-2*113*CLHEP::um-1.0*dp_CrossQubitNbDimX-2*82*CLHEP::um+1*CLHEP::um-40*CLHEP::um,0),Square_Capacitor_Q_To_P_VS1Logical_GP,"Square_Capacitor_Q_To_P_VS1Phys_GP",GroundPlaneLogical,false,0,true);
G4VPhysicalVolume* Square_Capacitor_Q_To_P_VS2Phys_GP=new G4PVPlacement(0,G4ThreeVector(X_Center_15[1],Y_Center_15[1]-1*495*CLHEP::um-1.0*dp_CrossQubitNbDimX-2*140*CLHEP::um+4*CLHEP::um-40*CLHEP::um,0),Square_Capacitor_Q_To_P_VS2Logical_GP,"Square_Capacitor_Q_To_P_VS2Phys_GP",GroundPlaneLogical,false,0,true);
G4VPhysicalVolume* Square_Capacitor_Q_To_P_VS3Phys_GP=new G4PVPlacement(0,G4ThreeVector(X_Center_15[2],Y_Center_15[2]-1*620*CLHEP::um-1.0*dp_CrossQubitNbDimX-2*140*CLHEP::um+4*CLHEP::um-40*CLHEP::um,0),Square_Capacitor_Q_To_P_VS3Logical_GP,"Square_Capacitor_Q_To_P_VS3Phys_GP",GroundPlaneLogical,false,0,true);
G4VPhysicalVolume* Square_Capacitor_Q_To_P_VS4Phys_GP=new G4PVPlacement(0,G4ThreeVector(X_Center_15[3],Y_Center_15[3]-1*620*CLHEP::um-1.0*dp_CrossQubitNbDimX-2*140*CLHEP::um+4*CLHEP::um-40*CLHEP::um,0),Square_Capacitor_Q_To_P_VS4Logical_GP,"Square_Capacitor_Q_To_P_VS4Phys_GP",GroundPlaneLogical,false,0,true);



G4VPhysicalVolume* Square_Capacitor_Q_To_P_VS1Phys=new G4PVPlacement(0,G4ThreeVector(),Square_Capacitor_Q_To_P_VS1Logical,"Square_Capacitor_Q_To_P_VS1Phys",Square_Capacitor_Q_To_P_VS1Logical_GP,false,0,true);
G4VPhysicalVolume* Square_Capacitor_Q_To_P_VS2Phys=new G4PVPlacement(0,G4ThreeVector(),Square_Capacitor_Q_To_P_VS2Logical,"Square_Capacitor_Q_To_P_VS2Phys",Square_Capacitor_Q_To_P_VS2Logical_GP,false,0,true);
G4VPhysicalVolume* Square_Capacitor_Q_To_P_VS3Phys=new G4PVPlacement(0,G4ThreeVector(),Square_Capacitor_Q_To_P_VS3Logical,"Square_Capacitor_Q_To_P_VS3Phys",Square_Capacitor_Q_To_P_VS3Logical_GP,false,0,true);
G4VPhysicalVolume* Square_Capacitor_Q_To_P_VS4Phys=new G4PVPlacement(0,G4ThreeVector(),Square_Capacitor_Q_To_P_VS4Logical,"Square_Capacitor_Q_To_P_VS4Phys",Square_Capacitor_Q_To_P_VS4Logical_GP,false,0,true);
//
//



/////////////////---------_Defining the Right and Left Circles for the conenctions----------
G4VSolid* Half_Cirlce_C_To_P= new G4Tubs("Half_Cirlce_C_To_P",  90.0* CLHEP::um,110.0* CLHEP::um,dp_resonatorAssemblyBaseNbDimZ_1,0*M_PI*rad,0.5*M_PI*rad ); // segment angle
G4LogicalVolume* Half_Cirlce_C_To_PLogical = new G4LogicalVolume( Half_Cirlce_C_To_P,fNiobium,"Half_Cirlce_C_To_PLogical");
G4VSolid* Half_Cirlce_C_To_P_GP= new G4Tubs("Half_Cirlce_C_To_P_GP",  90.0* CLHEP::um-Vacu_GroundPlane,110.0* CLHEP::um+Vacu_GroundPlane,dp_resonatorAssemblyBaseNbDimZ_1,0*M_PI*rad,0.5*M_PI*rad ); // segment angle
G4LogicalVolume* Half_Cirlce_C_To_PLogical_GP = new G4LogicalVolume(Half_Cirlce_C_To_P_GP ,fLiquidHelium,"Half_Cirlce_C_To_PLogical_GP");

G4RotationMatrix* XYRotC = new G4RotationMatrix;
G4RotationMatrix* XYRotC1 = new G4RotationMatrix;
XYRotC->rotateZ(1*M_PI/2.*rad);
XYRotC1->rotateZ(1.0*M_PI*rad);

G4VPhysicalVolume* Half_Cirlce_C_To_PPhys_GP=new G4PVPlacement(XYRotC,G4ThreeVector(X_Center_15[0]-100.0*CLHEP::um,Y_Center_15[0]-2*113*CLHEP::um-1.0*dp_CrossQubitNbDimX-2*82*CLHEP::um+3*CLHEP::um-155*CLHEP::um,0),Half_Cirlce_C_To_PLogical_GP,"Half_Cirlce_C_To_PPhys_GP",GroundPlaneLogical,false,0,true);
G4VPhysicalVolume* Half_Cirlce_C_To_PPhys2_GP=new G4PVPlacement(XYRotC,G4ThreeVector(X_Center_15[1]-100.0*CLHEP::um,Y_Center_15[1]-2*495*CLHEP::um-1.0*dp_CrossQubitNbDimX-2*82*CLHEP::um+3*CLHEP::um-155*CLHEP::um,0),Half_Cirlce_C_To_PLogical_GP,"Half_Cirlce_C_To_PPhys2_GP",GroundPlaneLogical,false,0,true);
G4VPhysicalVolume* Half_Cirlce_C_To_PPhys3_GP=new G4PVPlacement(XYRotC,G4ThreeVector(X_Center_15[2]-100.0*CLHEP::um,Y_Center_15[2]-2*620*CLHEP::um-1.0*dp_CrossQubitNbDimX-2*82*CLHEP::um+3*CLHEP::um-155*CLHEP::um,0),Half_Cirlce_C_To_PLogical_GP,"Half_Cirlce_C_To_PPhys3_GP",GroundPlaneLogical,false,0,true);
G4VPhysicalVolume* Half_Cirlce_C_To_PPhys4_GP=new G4PVPlacement(XYRotC1,G4ThreeVector(X_Center_15[3]+100.0*CLHEP::um,Y_Center_15[3]-2*620*CLHEP::um-1.0*dp_CrossQubitNbDimX-2*82*CLHEP::um+3*CLHEP::um-155*CLHEP::um,0),Half_Cirlce_C_To_PLogical_GP,"Half_Cirlce_C_To_PPhys4_GP",GroundPlaneLogical,false,0,true);


G4VPhysicalVolume* Half_Cirlce_C_To_PPhys=new G4PVPlacement(0,G4ThreeVector(),Half_Cirlce_C_To_PLogical,"Half_Cirlce_C_To_PPhys",Half_Cirlce_C_To_PLogical_GP,false,0,true);
G4VPhysicalVolume* Half_Cirlce_C_To_PPhys2=new G4PVPlacement(0,G4ThreeVector(),Half_Cirlce_C_To_PLogical,"Half_Cirlce_C_To_PPhys2",Half_Cirlce_C_To_PLogical_GP,false,0,true);
G4VPhysicalVolume* Half_Cirlce_C_To_PPhys3=new G4PVPlacement(0,G4ThreeVector(),Half_Cirlce_C_To_PLogical,"Half_Cirlce_C_To_PPhys3",Half_Cirlce_C_To_PLogical_GP,false,0,true);
G4VPhysicalVolume* Half_Cirlce_C_To_PPhys4=new G4PVPlacement(0,G4ThreeVector(),Half_Cirlce_C_To_PLogical,"Half_Cirlce_C_To_PPhys4",Half_Cirlce_C_To_PLogical_GP,false,0,true);




////////////////////////////////////////Horizontal  Lines  ////////////////////
G4VSolid* Square_Capacitor_Q_To_P_HS1= new G4Box("Square_Capacitor_Q_To_P_HS1",770.0*CLHEP::um,10.0*CLHEP::um,dp_resonatorAssemblyBaseNbDimZ_1);
G4VSolid* Square_Capacitor_Q_To_P_HS2= new G4Box("Square_Capacitor_Q_To_P_HS2",1350.0*CLHEP::um,10*CLHEP::um,dp_resonatorAssemblyBaseNbDimZ_1);
G4VSolid* Square_Capacitor_Q_To_P_HS3= new G4Box("Square_Capacitor_Q_To_P_HS3",230.0*CLHEP::um,10*CLHEP::um,dp_resonatorAssemblyBaseNbDimZ_1);
G4VSolid* Square_Capacitor_Q_To_P_HS4= new G4Box("Square_Capacitor_Q_To_P_HS4",490.0*CLHEP::um,10*CLHEP::um,dp_resonatorAssemblyBaseNbDimZ_1);

G4VSolid* Square_Capacitor_Q_To_P_HS1_GP= new G4Box("Square_Capacitor_Q_To_P_HS1_GP",770.0*CLHEP::um,10.0*CLHEP::um+Vacu_GroundPlane,dp_resonatorAssemblyBaseNbDimZ_1);
G4VSolid* Square_Capacitor_Q_To_P_HS2_GP= new G4Box("Square_Capacitor_Q_To_P_HS2_GP",1350.0*CLHEP::um,10*CLHEP::um+Vacu_GroundPlane,dp_resonatorAssemblyBaseNbDimZ_1);
G4VSolid* Square_Capacitor_Q_To_P_HS3_GP= new G4Box("Square_Capacitor_Q_To_P_HS3_GP",230.0*CLHEP::um,10*CLHEP::um+Vacu_GroundPlane,dp_resonatorAssemblyBaseNbDimZ_1);
G4VSolid* Square_Capacitor_Q_To_P_HS4_GP= new G4Box("Square_Capacitor_Q_To_P_HS4_GP",490.0*CLHEP::um,10*CLHEP::um+Vacu_GroundPlane,dp_resonatorAssemblyBaseNbDimZ_1);

G4LogicalVolume* Square_Capacitor_Q_To_P_HS1Logical = new G4LogicalVolume( Square_Capacitor_Q_To_P_HS1,fNiobium,"Square_Capacitor_Q_To_P_HS1Logical");
G4LogicalVolume* Square_Capacitor_Q_To_P_HS2Logical = new G4LogicalVolume( Square_Capacitor_Q_To_P_HS2,fNiobium,"Square_Capacitor_Q_To_P_HS2Logical");
G4LogicalVolume* Square_Capacitor_Q_To_P_HS3Logical = new G4LogicalVolume( Square_Capacitor_Q_To_P_HS3,fNiobium,"Square_Capacitor_Q_To_P_HS3Logical");
G4LogicalVolume* Square_Capacitor_Q_To_P_HS4Logical = new G4LogicalVolume( Square_Capacitor_Q_To_P_HS4,fNiobium,"Square_Capacitor_Q_To_P_HS4Logical");

G4LogicalVolume* Square_Capacitor_Q_To_P_HS1Logical_GP = new G4LogicalVolume( Square_Capacitor_Q_To_P_HS1_GP,fLiquidHelium,"Square_Capacitor_Q_To_P_HS1Logical_GP");
G4LogicalVolume* Square_Capacitor_Q_To_P_HS2Logical_GP = new G4LogicalVolume( Square_Capacitor_Q_To_P_HS2_GP,fLiquidHelium,"Square_Capacitor_Q_To_P_HS2Logical_GP");
G4LogicalVolume* Square_Capacitor_Q_To_P_HS3Logical_GP = new G4LogicalVolume( Square_Capacitor_Q_To_P_HS3_GP,fLiquidHelium,"Square_Capacitor_Q_To_P_HS3Logical_GP");
G4LogicalVolume* Square_Capacitor_Q_To_P_HS4Logical_GP = new G4LogicalVolume( Square_Capacitor_Q_To_P_HS4_GP,fLiquidHelium,"Square_Capacitor_Q_To_P_HS4Logical_GP");



G4VPhysicalVolume* Square_Capacitor_Q_To_P_HS1Phys_GP=new G4PVPlacement(0,G4ThreeVector(X_Center_15[0]-1*770*CLHEP::um-100*CLHEP::um,Y_Center_15[0]-2*113*CLHEP::um-1.0*dp_CrossQubitNbDimX-2*113*CLHEP::um+1*CLHEP::um-191*CLHEP::um,0),Square_Capacitor_Q_To_P_HS1Logical_GP,"Square_Capacitor_Q_To_P_HS1Phys_GP",GroundPlaneLogical,false,0,true);
G4VPhysicalVolume* Square_Capacitor_Q_To_P_HS2Phys_GP=new G4PVPlacement(0,G4ThreeVector(X_Center_15[1]-1*1350*CLHEP::um-100*CLHEP::um,Y_Center_15[1]-1*495*CLHEP::um-1.0*dp_CrossQubitNbDimX-2*495*CLHEP::um+4*CLHEP::um+75*CLHEP::um,0),Square_Capacitor_Q_To_P_HS2Logical_GP,"Square_Capacitor_Q_To_P_HS2Phys_GP",GroundPlaneLogical,false,0,true);
G4VPhysicalVolume* Square_Capacitor_Q_To_P_HS3Phys_GP=new G4PVPlacement(0,G4ThreeVector(X_Center_15[2]-1*230*CLHEP::um-100*CLHEP::um,Y_Center_15[2]-1*620*CLHEP::um-1.0*dp_CrossQubitNbDimX-2*620*CLHEP::um+4*CLHEP::um+200*CLHEP::um,0),Square_Capacitor_Q_To_P_HS3Logical_GP,"Square_Capacitor_Q_To_P_HS3Phys_GP",GroundPlaneLogical,false,0,true);
G4VPhysicalVolume* Square_Capacitor_Q_To_P_HS4Phys_GP=new G4PVPlacement(0,G4ThreeVector(X_Center_15[3]+1*490*CLHEP::um+100*CLHEP::um,Y_Center_15[3]-1*620*CLHEP::um-1.0*dp_CrossQubitNbDimX-2*620*CLHEP::um+4*CLHEP::um+200*CLHEP::um,0),Square_Capacitor_Q_To_P_HS4Logical_GP,"Square_Capacitor_Q_To_P_HS4Phys_GP",GroundPlaneLogical,false,0,true);

G4VPhysicalVolume* Square_Capacitor_Q_To_P_HS1Phys=new G4PVPlacement(0,G4ThreeVector(),Square_Capacitor_Q_To_P_HS1Logical,"Square_Capacitor_Q_To_P_HS1Phys",Square_Capacitor_Q_To_P_HS1Logical_GP,false,0,true);
G4VPhysicalVolume* Square_Capacitor_Q_To_P_HS2Phys=new G4PVPlacement(0,G4ThreeVector(),Square_Capacitor_Q_To_P_HS2Logical,"Square_Capacitor_Q_To_P_HS2Phys",Square_Capacitor_Q_To_P_HS2Logical_GP,false,0,true);
G4VPhysicalVolume* Square_Capacitor_Q_To_P_HS3Phys=new G4PVPlacement(0,G4ThreeVector(),Square_Capacitor_Q_To_P_HS3Logical,"Square_Capacitor_Q_To_P_HS3Phys",Square_Capacitor_Q_To_P_HS3Logical_GP,false,0,true);
G4VPhysicalVolume* Square_Capacitor_Q_To_P_HS4Phys=new G4PVPlacement(0,G4ThreeVector(),Square_Capacitor_Q_To_P_HS4Logical,"Square_Capacitor_Q_To_P_HS4Phys",Square_Capacitor_Q_To_P_HS4Logical_GP,false,0,true);



//////////Connecting the Circles
G4RotationMatrix* XYRotC2 = new G4RotationMatrix;
G4RotationMatrix* XYRotC3 = new G4RotationMatrix;
XYRotC2->rotateZ(1.0*M_PI*rad);
XYRotC3->rotateZ(1.5*M_PI*rad);
G4VPhysicalVolume* Half_Cirlce_C_To_PPhys_P_1_GP=new G4PVPlacement(XYRotC2,G4ThreeVector(X_Center_15[0]-2*770*CLHEP::um-100*CLHEP::um,Y_Center_15[0]-2*113*CLHEP::um-1.0*dp_CrossQubitNbDimX-2*113*CLHEP::um+1*CLHEP::um-91*CLHEP::um,0),Half_Cirlce_C_To_PLogical_GP,"Half_Cirlce_C_To_PPhys_P_1_GP",GroundPlaneLogical,false,0,true);
G4VPhysicalVolume* Half_Cirlce_C_To_PPhys_P_2_GP=new G4PVPlacement(XYRotC3,G4ThreeVector(X_Center_15[1]-2*1350*CLHEP::um-100*CLHEP::um,Y_Center_15[1]-1*495*CLHEP::um-1.0*dp_CrossQubitNbDimX-2*495*CLHEP::um+4*CLHEP::um-25*CLHEP::um,0),Half_Cirlce_C_To_PLogical_GP,"Half_Cirlce_C_To_PPhys_P_2_GP",GroundPlaneLogical,false,0,true);
G4VPhysicalVolume* Half_Cirlce_C_To_PPhys_P_3_GP=new G4PVPlacement(XYRotC3,G4ThreeVector(X_Center_15[2]-2*230*CLHEP::um-100*CLHEP::um,Y_Center_15[2]-1*620*CLHEP::um-1.0*dp_CrossQubitNbDimX-2*620*CLHEP::um+4*CLHEP::um+100*CLHEP::um,0),Half_Cirlce_C_To_PLogical_GP,"Half_Cirlce_C_To_PPhys_P_3_GP",GroundPlaneLogical,false,0,true);
G4VPhysicalVolume* Half_Cirlce_C_To_PPhys_P_4_GP=new G4PVPlacement(0,G4ThreeVector(X_Center_15[3]+2*490*CLHEP::um+100*CLHEP::um,Y_Center_15[3]-1*620*CLHEP::um-1.0*dp_CrossQubitNbDimX-2*620*CLHEP::um+4*CLHEP::um+100*CLHEP::um,0),Half_Cirlce_C_To_PLogical_GP,"Half_Cirlce_C_To_PPhys_P_4_GP",GroundPlaneLogical,false,0,true);

G4VPhysicalVolume* Half_Cirlce_C_To_PPhys_P_1=new G4PVPlacement(0,G4ThreeVector(),Half_Cirlce_C_To_PLogical,"Half_Cirlce_C_To_PPhys_P_1",Half_Cirlce_C_To_PLogical_GP,false,0,true);
G4VPhysicalVolume* Half_Cirlce_C_To_PPhys_P_2=new G4PVPlacement(0,G4ThreeVector(),Half_Cirlce_C_To_PLogical,"Half_Cirlce_C_To_PPhys_P_2",Half_Cirlce_C_To_PLogical_GP,false,0,true);
G4VPhysicalVolume* Half_Cirlce_C_To_PPhys_P_3=new G4PVPlacement(0,G4ThreeVector(),Half_Cirlce_C_To_PLogical,"Half_Cirlce_C_To_PPhys_P_3",Half_Cirlce_C_To_PLogical_GP,false,0,true);
G4VPhysicalVolume* Half_Cirlce_C_To_PPhys_P_4=new G4PVPlacement(0,G4ThreeVector(),Half_Cirlce_C_To_PLogical,"Half_Cirlce_C_To_PPhys_P_4",Half_Cirlce_C_To_PLogical_GP,false,0,true);




/////////////-----------
/////////Final Connectionss for The panels for Qubit 1 and Qubit 2
//One Verrical Line and One Cirlce For Qubit 1
G4VSolid* Square_Capacitor_Q_To_P_V_P1= new G4Box("Square_Capacitor_Q_To_P_VS1",10.0*CLHEP::um,190*CLHEP::um,dp_resonatorAssemblyBaseNbDimZ_1);
G4LogicalVolume* Square_Capacitor_Q_To_P_V_P1Logical = new G4LogicalVolume( Square_Capacitor_Q_To_P_V_P1,fNiobium,"Square_Capacitor_Q_To_P_V_P1Logical");
G4VSolid* Square_Capacitor_Q_To_P_V_P1_GP= new G4Box("Square_Capacitor_Q_To_P_VS1_GP",10.0*CLHEP::um+Vacu_GroundPlane,190*CLHEP::um,dp_resonatorAssemblyBaseNbDimZ_1);
G4LogicalVolume* Square_Capacitor_Q_To_P_V_P1Logical_GP = new G4LogicalVolume( Square_Capacitor_Q_To_P_V_P1_GP,fLiquidHelium,"Square_Capacitor_Q_To_P_V_P1Logical_GP");

constexpr double X_Center_16=X_Center_15[0]-2*770*CLHEP::um-200*CLHEP::um;
constexpr double Y_Center_16=Y_Center_15[0]-2*113*CLHEP::um-1.0*dp_CrossQubitNbDimX-2*113*CLHEP::um+1*CLHEP::um+99*CLHEP::um;

G4VPhysicalVolume* Square_Capacitor_Q_To_P_V_P1Phys_GP=new G4PVPlacement(0,G4ThreeVector(X_Center_16,Y_Center_16,0),Square_Capacitor_Q_To_P_V_P1Logical_GP,"Square_Capacitor_Q_To_P_V_P1Phys_GP",GroundPlaneLogical,false,0,true);
G4VPhysicalVolume* Half_Cirlce_C_To_PPhys_P_5_GP=new G4PVPlacement(0,G4ThreeVector(X_Center_16-100*CLHEP::um,Y_Center_16+190*CLHEP::um,0),Half_Cirlce_C_To_PLogical_GP,"Half_Cirlce_C_To_PPhys_P_5_GP",GroundPlaneLogical,false,0,true);
G4VPhysicalVolume* Square_Capacitor_Q_To_P_V_P1Phys=new G4PVPlacement(0,G4ThreeVector(),Square_Capacitor_Q_To_P_V_P1Logical,"Square_Capacitor_Q_To_P_V_P1Phys",Square_Capacitor_Q_To_P_V_P1Logical_GP,false,0,true);
G4VPhysicalVolume* Half_Cirlce_C_To_PPhys_P_5=new G4PVPlacement(0,G4ThreeVector(),Half_Cirlce_C_To_PLogical,"Half_Cirlce_C_To_PPhys_P_5",Half_Cirlce_C_To_PLogical_GP,false,0,true);


////////////////////////??Connecting the Pannelsss///////////////////////
//////Defining The coordinates for the pannels

constexpr double X_Center_17[3]={X_Center_15[1]-2*1350*CLHEP::um-100*CLHEP::um,
X_Center_15[2]-2*230*CLHEP::um-100*CLHEP::um,X_Center_15[3]+2*490*CLHEP::um+100*CLHEP::um};
constexpr double Y_Center_17[3]={Y_Center_15[1]-1*495*CLHEP::um-1.0*dp_CrossQubitNbDimX-2*495*CLHEP::um+4*CLHEP::um-25*CLHEP::um,
Y_Center_15[2]-1*620*CLHEP::um-1.0*dp_CrossQubitNbDimX-2*620*CLHEP::um+4*CLHEP::um+100*CLHEP::um,
Y_Center_15[3]-1*620*CLHEP::um-1.0*dp_CrossQubitNbDimX-2*620*CLHEP::um+4*CLHEP::um+100*CLHEP::um};
////////////////


////////////---------------RENAME ALL THE QUBITS FOR ANNALYSIS ON THE DETECTOR ///-----------------------
G4VSolid*Calis=new G4Tet("Calis",G4ThreeVector(0,0,30*CLHEP::um),G4ThreeVector(0*CLHEP::um,150*CLHEP::um,0),
G4ThreeVector(0*CLHEP::um,-150*CLHEP::um,0),G4ThreeVector(310*CLHEP::um,0*CLHEP::um,0),0);
G4VSolid*Calis_GP=new G4Tet("Calis_GP",G4ThreeVector(0,0,30*CLHEP::um),G4ThreeVector(0*CLHEP::um,150*CLHEP::um+2*75*CLHEP::um,0),
G4ThreeVector(0*CLHEP::um,-150*CLHEP::um-2*75*CLHEP::um,0),G4ThreeVector(310*CLHEP::um,0*CLHEP::um,0),0);

G4VSolid* Square_Calis= new G4Box("Square_Calis",10.0*CLHEP::um,20*CLHEP::um,20*CLHEP::um);
G4VSolid* Big_Trian_Panel_0= new G4SubtractionSolid("Calis-Square_Calis",Calis, Square_Calis,0, G4ThreeVector(300*CLHEP::um,0*CLHEP::um,0.));
G4VSolid* Big_Trian_Panel_0_GP= new G4SubtractionSolid("Calis_GP-Square_Calis",Calis_GP, Square_Calis,0, G4ThreeVector(300*CLHEP::um,0*CLHEP::um,0.));
G4VSolid* Square_Panel= new G4Box("Square_Panel",150.0*CLHEP::um,175*CLHEP::um,10*CLHEP::um);
G4VSolid* Square_Panel_GP= new G4Box("Square_Panel_GP",150.0*CLHEP::um+2*75*CLHEP::um,175*CLHEP::um,10*CLHEP::um);
G4VSolid* Square_Panel_GP1= new G4Box("Square_Panel_GP1",300*CLHEP::um,75*CLHEP::um,10*CLHEP::um);




// G4VSolid* Big_Trian_Panel= new G4UnionSolid("Square_Panel+Big_Trian_Panel_0", Square_Panel, Big_Trian_Panel_0,nullptr,G4ThreeVector(0.0*CLHEP::um, 0., 0.));
G4RotationMatrix* XYRotC4 = new G4RotationMatrix;
XYRotC4->rotateZ(-0.5*M_PI*rad);
G4RotationMatrix* XYRotC41 = new G4RotationMatrix;
XYRotC41->rotateZ(0.5*M_PI*rad);
 G4LogicalVolume* Big_Trian_PanelLogical = new G4LogicalVolume( Big_Trian_Panel_0,fNiobium,"Big_Trian_PanelLogical");
 G4LogicalVolume* Big_Trian_PanelLogical_GP = new G4LogicalVolume( Big_Trian_Panel_0_GP,fLiquidHelium,"Big_Trian_PanelLogical_GP");

  G4LogicalVolume* Square_PanelLogical = new G4LogicalVolume( Square_Panel,fNiobium,"Square_PanelLogical");
  G4LogicalVolume* Square_PanelLogical_GP = new G4LogicalVolume( Square_Panel_GP,fLiquidHelium,"Square_PanelLogical_GP");
  G4LogicalVolume* Square_PanelLogical_GP1 = new G4LogicalVolume( Square_Panel_GP1,fLiquidHelium,"Square_PanelLogical_GP1");
 // G4VPhysicalVolume* Big_Trian_PanelPhys=new G4PVPlacement(XYRotC4,G4ThreeVector(X_Center_17[2]+100*CLHEP::um,Y_Center_17[2]-2*dp_Square_panelNbDimX1+60*CLHEP::um,0.26*cm),Big_Trian_PanelLogical,"Big_Trian_PanelPhys",worldLogical,false,0,true);

//
// G4ThreeVector(X_Center_16-570*CLHEP::um,Y_Center_16+290*CLHEP::um,Position_Superconductor)

for(G4int i=0; i<6;i++){
if(i<1){
  G4VPhysicalVolume* Big_Trian_PanelPhys_GP=new G4PVPlacement(0,G4ThreeVector(X_Center_16-390*CLHEP::um,Y_Center_16+290*CLHEP::um,0),Big_Trian_PanelLogical_GP,"Big_Trian_PanelPhys_GP",GroundPlaneLogical,false,i,true);
  G4VPhysicalVolume* Square_PanelPhys_GP=new G4PVPlacement(XYRotC41,G4ThreeVector(X_Center_16-390*CLHEP::um-175*CLHEP::um,Y_Center_16+290*CLHEP::um,0),Square_PanelLogical_GP,"Square_PanellPhys_GP",GroundPlaneLogical,false,i,true);
  G4VPhysicalVolume* Square_PanelPhys_GP1=new G4PVPlacement(XYRotC41,G4ThreeVector(X_Center_16-365*CLHEP::um-3*150*CLHEP::um,Y_Center_16+290*CLHEP::um,0),Square_PanelLogical_GP1,"Square_PanellPhys_GP1",GroundPlaneLogical,false,i,true);
  G4CMPLogicalBorderSurface * border_sapphire_Vacuum_Big_Trian_PanelPhys_GP=new G4CMPLogicalBorderSurface("border_sapphire_Vacuum_Big_Trian_PanelPhys_GP",SapphirePhys,Big_Trian_PanelPhys_GP,fSiVacuumInterface);
  G4CMPLogicalBorderSurface * border_sapphire_Vacuum_Square_PanelPhys_GP=new G4CMPLogicalBorderSurface("border_sapphire_Vacuum_Square_PanelPhys_GP",SapphirePhys,Square_PanelPhys_GP,fSiVacuumInterface);
  G4CMPLogicalBorderSurface * border_sapphire_Vacuum_Square_PanelPhys_GP1=new G4CMPLogicalBorderSurface("border_sapphire_Vacuum_Square_PanelPhys_GP1",SapphirePhys,Square_PanelPhys_GP1,fSiVacuumInterface);

  G4VPhysicalVolume* Big_Trian_PanelPhys=new G4PVPlacement(0,G4ThreeVector(),Big_Trian_PanelLogical,"Big_Trian_PanelPhys",Big_Trian_PanelLogical_GP,false,i,true);
  G4VPhysicalVolume* Square_PanelPhys=new G4PVPlacement(0,G4ThreeVector(),Square_PanelLogical,"Square_PanellPhys",Square_PanelLogical_GP,false,i,true);
  G4CMPLogicalBorderSurface * border_sapphire_Vacuum_Big_Trian_PanelPhys=new G4CMPLogicalBorderSurface("border_sapphire_Vacuum_Big_Trian_PanelPhys",SapphirePhys,Big_Trian_PanelPhys,fSiNbInterface);
  G4CMPLogicalBorderSurface * border_sapphire_Vacuum_Square_PanelPhys=new G4CMPLogicalBorderSurface("border_sapphire_Vacuum_Square_PanelPhys",SapphirePhys,Square_PanelPhys,fSiNbInterface);



}



  else if( i>=1 && i<3){
  G4VPhysicalVolume* Big_Trian_PanelPhys_GP=new G4PVPlacement(XYRotC41,G4ThreeVector(pow(-1,i)*(dp_TransmissionLineBaseNbDimX-2*dp_TransmissionLineBaseNbDimX_1-90.0* CLHEP::um),Y_Center_TransmissionLine+6.5*90.0*CLHEP::um-25.0*CLHEP::um,0),Big_Trian_PanelLogical_GP,"Big_Trian_PanelPhys_GP",GroundPlaneLogical,false,i,true);
  G4VPhysicalVolume* Square_PanelPhys_GP=new G4PVPlacement(0,G4ThreeVector(pow(-1,i)*(dp_TransmissionLineBaseNbDimX-2*dp_TransmissionLineBaseNbDimX_1-90.0* CLHEP::um),Y_Center_TransmissionLine+6.5*90.0*CLHEP::um+150.0*CLHEP::um,0),Square_PanelLogical_GP,"Square_PanellPhys_GP",GroundPlaneLogical,false,i,true);
G4VPhysicalVolume* Square_PanelPhys_GP1=new G4PVPlacement(0,G4ThreeVector(pow(-1,i)*(dp_TransmissionLineBaseNbDimX-2*dp_TransmissionLineBaseNbDimX_1-90.0* CLHEP::um),Y_Center_TransmissionLine+6.5*90.0*CLHEP::um+3*134*CLHEP::um-2*CLHEP::um,0),Square_PanelLogical_GP1,"Square_PanellPhys_GP1",GroundPlaneLogical,false,i,true);
G4CMPLogicalBorderSurface * border_sapphire_Vacuum_Big_Trian_PanelPhys_GP=new G4CMPLogicalBorderSurface("border_sapphire_Vacuum_Big_Trian_PanelPhys_GP",SapphirePhys,Big_Trian_PanelPhys_GP,fSiVacuumInterface);
G4CMPLogicalBorderSurface * border_sapphire_Vacuum_Square_PanelPhys_GP=new G4CMPLogicalBorderSurface("border_sapphire_Vacuum_Square_PanelPhys_GP",SapphirePhys,Square_PanelPhys_GP,fSiVacuumInterface);
G4CMPLogicalBorderSurface * border_sapphire_Vacuum_Square_PanelPhys_GP1=new G4CMPLogicalBorderSurface("border_sapphire_Vacuum_Square_PanelPhys_GP1",SapphirePhys,Square_PanelPhys_GP1,fSiVacuumInterface);


G4VPhysicalVolume* Big_Trian_PanelPhys=new G4PVPlacement(0,G4ThreeVector(),Big_Trian_PanelLogical,"Big_Trian_PanelPhys",Big_Trian_PanelLogical_GP,false,i,true);
G4VPhysicalVolume*  Square_PanelPhys=new G4PVPlacement(0,G4ThreeVector(), Square_PanelLogical," Square_PanelPhys",Square_PanelLogical_GP,false,i,true);
G4CMPLogicalBorderSurface * border_sapphire_Vacuum_Big_Trian_PanelPhys=new G4CMPLogicalBorderSurface("border_sapphire_Vacuum_Big_Trian_PanelPhys",SapphirePhys,Big_Trian_PanelPhys,fSiNbInterface);
G4CMPLogicalBorderSurface * border_sapphire_Vacuum_Square_PanelPhys=new G4CMPLogicalBorderSurface("border_sapphire_Vacuum_Square_PanelPhys",SapphirePhys,Square_PanelPhys,fSiNbInterface);



}



//G4ThreeVector(pow(-1,i)*(dp_TransmissionLineBaseNbDimX-2*dp_TransmissionLineBaseNbDimX_1-90.0* CLHEP::um),Y_Center_TransmissionLine+8*90.0* CLHEP::um+20*CLHEP::um,Position_Superconductor)
else if(i>=3 && i<5) {
  G4VPhysicalVolume* Big_Trian_PanelPhys_GP=new G4PVPlacement(XYRotC4,G4ThreeVector(X_Center_17[i-3]-100*CLHEP::um,Y_Center_17[i-3]-2*dp_Square_panelNbDimX1+60*CLHEP::um,0),Big_Trian_PanelLogical_GP,"Big_Trian_PanelPhys_GP",GroundPlaneLogical,false,i,true);
  G4VPhysicalVolume* Square_PanelPhys_GP=new G4PVPlacement(0,G4ThreeVector(X_Center_17[i-3]-100*CLHEP::um,Y_Center_17[i-3]-2*dp_Square_panelNbDimX1+60*CLHEP::um-175*CLHEP::um,0),Square_PanelLogical_GP,"Square_PanellPhys_GP",GroundPlaneLogical,false,i,true);
  G4VPhysicalVolume* Square_PanelPhys_GP1=new G4PVPlacement(0,G4ThreeVector(X_Center_17[i-3]-100*CLHEP::um,Y_Center_17[i-3]-2*dp_Square_panelNbDimX1+60*CLHEP::um-3*140*CLHEP::um-5*CLHEP::um,0),Square_PanelLogical_GP1,"Square_PanellPhys_GP1",GroundPlaneLogical,false,i,true);
  G4CMPLogicalBorderSurface * border_sapphire_Vacuum_Big_Trian_PanelPhys_GP=new G4CMPLogicalBorderSurface("border_sapphire_Vacuum_Big_Trian_PanelPhys_GP",SapphirePhys,Big_Trian_PanelPhys_GP,fSiVacuumInterface);
  G4CMPLogicalBorderSurface * border_sapphire_Vacuum_Square_PanelPhys_GP=new G4CMPLogicalBorderSurface("border_sapphire_Vacuum_Square_PanelPhys_GP",SapphirePhys,Square_PanelPhys_GP,fSiVacuumInterface);
  G4CMPLogicalBorderSurface * border_sapphire_Vacuum_Square_PanelPhys_GP1=new G4CMPLogicalBorderSurface("border_sapphire_Vacuum_Square_PanelPhys_GP1",SapphirePhys,Square_PanelPhys_GP1,fSiVacuumInterface);



  G4VPhysicalVolume* Big_Trian_PanelPhys=new G4PVPlacement(0,G4ThreeVector(),Big_Trian_PanelLogical,"Big_Trian_PanelPhys",Big_Trian_PanelLogical_GP,false,i,true);
G4VPhysicalVolume* Square_PanelPhys=new G4PVPlacement(0,G4ThreeVector(),Square_PanelLogical,"Square_PanellPhys",Square_PanelLogical_GP,false,i,true);
G4CMPLogicalBorderSurface * border_sapphire_Vacuum_Big_Trian_PanelPhys=new G4CMPLogicalBorderSurface("border_sapphire_Vacuum_Big_Trian_PanelPhys",SapphirePhys,Big_Trian_PanelPhys,fSiNbInterface);
G4CMPLogicalBorderSurface * border_sapphire_Vacuum_Square_PanelPhys=new G4CMPLogicalBorderSurface("border_sapphire_Vacuum_Square_PanelPhys",SapphirePhys,Square_PanelPhys,fSiNbInterface);




}

else{
  G4VPhysicalVolume* Big_Trian_PanelPhys_GP=new G4PVPlacement(XYRotC4,G4ThreeVector(X_Center_17[i-3]+100*CLHEP::um,Y_Center_17[i-3]-2*dp_Square_panelNbDimX1+60*CLHEP::um,0),Big_Trian_PanelLogical_GP,"Big_Trian_PanelPhys_GP",GroundPlaneLogical,false,i,true);
  G4VPhysicalVolume* Square_PanelPhys_GP=new G4PVPlacement(0,G4ThreeVector(X_Center_17[i-3]+100*CLHEP::um,Y_Center_17[i-3]-2*dp_Square_panelNbDimX1+60*CLHEP::um-175*CLHEP::um,0),Square_PanelLogical_GP,"Square_PanellPhys_GP",GroundPlaneLogical,false,i,true);
  G4VPhysicalVolume* Square_PanelPhys_GP1=new G4PVPlacement(0,G4ThreeVector(X_Center_17[i-3]+100*CLHEP::um,Y_Center_17[i-3]-2*dp_Square_panelNbDimX1+60*CLHEP::um-3*140*CLHEP::um-5*CLHEP::um,0),Square_PanelLogical_GP1,"Square_PanellPhys_GP1",GroundPlaneLogical,false,i,true);
  G4CMPLogicalBorderSurface * border_sapphire_Vacuum_Big_Trian_PanelPhys_GP=new G4CMPLogicalBorderSurface("border_sapphire_Vacuum_Big_Trian_PanelPhys_GP",SapphirePhys,Big_Trian_PanelPhys_GP,fSiVacuumInterface);
  G4CMPLogicalBorderSurface * border_sapphire_Vacuum_Square_PanelPhys_GP=new G4CMPLogicalBorderSurface("border_sapphire_Vacuum_Square_PanelPhys_GP",SapphirePhys,Square_PanelPhys_GP,fSiVacuumInterface);
  G4CMPLogicalBorderSurface * border_sapphire_Vacuum_Square_PanelPhys_GP1=new G4CMPLogicalBorderSurface("border_sapphire_Vacuum_Square_PanelPhys_GP1",SapphirePhys,Square_PanelPhys_GP1,fSiVacuumInterface);



G4VPhysicalVolume* Big_Trian_PanelPhys=new G4PVPlacement(0,G4ThreeVector(),Big_Trian_PanelLogical,"Big_Trian_PanelPhys",Big_Trian_PanelLogical_GP,false,i,true);
G4VPhysicalVolume* Square_PanelPhys=new G4PVPlacement(0,G4ThreeVector(),Square_PanelLogical,"Square_PanellPhys",Square_PanelLogical_GP,false,i,true);
G4CMPLogicalBorderSurface * border_sapphire_Vacuum_Big_Trian_PanelPhys=new G4CMPLogicalBorderSurface("border_sapphire_Vacuum_Big_Trian_PanelPhys",SapphirePhys,Big_Trian_PanelPhys,fSiNbInterface);
G4CMPLogicalBorderSurface * border_sapphire_Vacuum_Square_PanelPhys=new G4CMPLogicalBorderSurface("border_sapphire_Vacuum_Square_PanelPhys",SapphirePhys,Square_PanelPhys,fSiNbInterface);




}

}
/////////////////////////-------------------------CREATING THE COOPPER BOX -------------------/////////////////////
///It is will be a Cylinder  and Inside of  the Cylinder the Sapphire, the sapphire only tocuh the corners
// I need to Substrate one Square minus 4 corners circlles
///The area that touch the Sapphgire with tyeh copper is 0.02 inc
///

G4VSolid* Copper_Cylinder= new G4Tubs("Copper_Cylinder",0.0, 1.905*cm,0.635*cm,0,4*tube_dPhi ); // segment angle
G4VSolid* Box_Copper_Air= new G4Box("Box_Copper",0.7*cm,0.7*cm,0.3*cm);
// G4VSolid* Box_Copper_Air= new G4Box("Box_Copper",(0.7-0.0508)*cm,(0.7-0.0508)*cm,0.3*cm);
// G4LogicalVolume* Copper_Box_AirLogical = new G4LogicalVolume(Box_Copper_Air,fLiquidHelium,"Copper_Box_AirLogical");
//Invicible
//Copper_Box_AirLogical->SetVisAttributes(WorldA);
// Copper_Box_AirLogical->SetVisAttributes(SuperColor);
G4VSolid* Box_Copper= new G4Box("Box_Copper",0.7*cm,0.7*cm,0.3*cm);
G4VSolid* Copper_Box = new G4SubtractionSolid("Copper_Cylinder-Box_Copper", Copper_Cylinder,Box_Copper, 0, G4ThreeVector(0.*cm,0.,(0.6356-0.3)*cm));
G4VSolid* Sonten_Copper1= new G4Box("Sonten_Copper",0.7*cm,0.0508*cm,0.3*cm-430*CLHEP::um);
G4VSolid* Sonten_Copper2= new G4Box("Sonten_Copper2",0.0508*cm,0.7*cm-2*0.0508*cm,0.3*cm-430*CLHEP::um);

G4LogicalVolume* Copper_BoxLogical = new G4LogicalVolume(Copper_Box,copper_mat,"Copper_BoxLogical");
G4LogicalVolume* Sonten_Copper1Logical = new G4LogicalVolume(Sonten_Copper1,copper_mat,"Sonten_Copper1Logical");
G4LogicalVolume* Sonten_Copper2Logical = new G4LogicalVolume(Sonten_Copper2,copper_mat,"Sonten_Copper2Logical");
constexpr double Z_Copper=-0.5*0.635*cm-7*430*CLHEP::um+259*CLHEP::um;

G4VPhysicalVolume* Copper_BoxPhys=new G4PVPlacement(0,G4ThreeVector(0,0,Z_Copper),Copper_BoxLogical,"Copper_BoxPhys",worldLogical,false,0,true);
G4VPhysicalVolume* Sonten_Copper1Phys1=new G4PVPlacement(0,G4ThreeVector(0,1.0*0.7*cm-0.5*0.0508*cm-254*CLHEP::um,Z_Copper+1*0.3*cm-74*CLHEP::um),Sonten_Copper1Logical,"Sonten_Copper1xPhys1",worldLogical,false,0,true);
G4VPhysicalVolume* Sonten_Copper1Phys2=new G4PVPlacement(0,G4ThreeVector(0,-1.0*0.7*cm+0.5*0.0508*cm+254*CLHEP::um,Z_Copper+1*0.3*cm-74*CLHEP::um),Sonten_Copper1Logical,"Sonten_Copper1xPhys2",worldLogical,false,0,true);
G4VPhysicalVolume* Sonten_Copper1Phys3=new G4PVPlacement(0,G4ThreeVector(1.0*0.7*cm-0.5*0.0508*cm-254*CLHEP::um,0,Z_Copper+1*0.3*cm-74*CLHEP::um),Sonten_Copper2Logical,"Sonten_Copper1xPhys3",worldLogical,false,0,true);
G4VPhysicalVolume* Sonten_Copper1Phys4=new G4PVPlacement(0,G4ThreeVector(-1.0*0.7*cm+0.5*0.0508*cm+254*CLHEP::um,0,Z_Copper+1*0.3*cm-74*CLHEP::um),Sonten_Copper2Logical,"Sonten_Copper1xPhys4",worldLogical,false,0,true);


G4VisAttributes* S= new G4VisAttributes(G4Colour(1.0,0.4,0.1,0.5)); // Green
Copper_BoxLogical->SetVisAttributes(S);

//Add a phonon sensor to the interface properties here.
  G4CMPLogicalBorderSurface * border_sapphire_Vacuum_Resonator_TransmissionLine=new G4CMPLogicalBorderSurface("border_sapphire_Vacuum_Resonator_TransmissionLine",SapphirePhys,Resonator_TransmissionLine,fSiNbInterface);
   G4CMPLogicalBorderSurface * border_sapphire_Vacuum_Half_Cirlce_RPhys_Transmission=new G4CMPLogicalBorderSurface("border_sapphire_Vacuum_Half_Cirlce_RPhys_Transmission",SapphirePhys,Half_Cirlce_RPhys_Transmission,fSiNbInterface);
 G4CMPLogicalBorderSurface * border_sapphire_Vacuum_Half_Cirlce_LPhys_Transmission=new G4CMPLogicalBorderSurface("border_sapphire_Vacuum_Half_Cirlce_LPhys_Transmission",SapphirePhys,Half_Cirlce_LPhys_Transmission,fSiNbInterface);
  G4CMPLogicalBorderSurface * border_sapphire_Vacuum_Quarter_Cirlce_RPhys_Transmission=new G4CMPLogicalBorderSurface("border_sapphire_Vacuum_Quarter_Cirlce_RPhys_Transmission",SapphirePhys,Quarter_Cirlce_RPhys_Transmission,fSiNbInterface);
   G4CMPLogicalBorderSurface * border_sapphire_Vacuum_Quarter_Cirlce_LPhys_Transmission=new G4CMPLogicalBorderSurface("border_sapphire_Vacuum_Quarter_Cirlce_LPhys_Transmission",SapphirePhys,Quarter_Cirlce_LPhys_Transmission,fSiNbInterface);
   G4CMPLogicalBorderSurface * border_sapphire_Vacuum_Resonator_TransmissionLine_1=new G4CMPLogicalBorderSurface("border_sapphire_Vacuum_Resonator_TransmissionLine_1",SapphirePhys,Resonator_TransmissionLine_1,fSiNbInterface);
   G4CMPLogicalBorderSurface * border_sapphire_Vacuum_Resonator_TransmissionLine_2=new G4CMPLogicalBorderSurface("border_sapphire_Vacuum_Resonator_TransmissionLine_2",SapphirePhys,Resonator_TransmissionLine_2,fSiNbInterface);
   // G4CMPLogicalBorderSurface * border_sapphire_Vacuum_Half_Cirlce_R_TTDPhys=new G4CMPLogicalBorderSurface("border_sapphire_Vacuum_Half_Cirlce_R_TTDPhys",SapphirePhys,Half_Cirlce_R_TTDPhys,fSiNbInterface);
    G4CMPLogicalBorderSurface * border_sapphire_Vacuum_Square_To_ResonatorLogicalPhys1=new G4CMPLogicalBorderSurface("border_sapphire_Vacuum_Square_To_ResonatorLogicalPhys1",SapphirePhys,Square_To_ResonatorLogicalPhys1,fSiNbInterface);
    G4CMPLogicalBorderSurface * border_sapphire_Vacuum_Square_To_ResonatorLogicalPhys2=new G4CMPLogicalBorderSurface("border_sapphire_Vacuum_Square_To_ResonatorLogicalPhys2",SapphirePhys,Square_To_ResonatorLogicalPhys2,fSiNbInterface);
    G4CMPLogicalBorderSurface * border_sapphire_Vacuum_Square_To_ResonatorLogicalPhys3=new G4CMPLogicalBorderSurface("border_sapphire_Vacuum_Square_To_ResonatorLogicalPhys3",SapphirePhys,Square_To_ResonatorLogicalPhys3,fSiNbInterface);
    G4CMPLogicalBorderSurface * border_sapphire_Vacuum_Square_To_ResonatorLogicalPhys4=new G4CMPLogicalBorderSurface("border_sapphire_Vacuum_Square_To_ResonatorLogicalPhys4",SapphirePhys,Square_To_ResonatorLogicalPhys4,fSiNbInterface);
   G4CMPLogicalBorderSurface * border_sapphire_Vacuum_Square_to_CapcitorH1Phys=new G4CMPLogicalBorderSurface("border_sapphire_Vacuum_Square_to_CapcitorH1Phys",SapphirePhys,Square_to_CapcitorH1Phys,fSiNbInterface);
   G4CMPLogicalBorderSurface * border_sapphire_Vacuum_Square_to_CapcitorH2Phys=new G4CMPLogicalBorderSurface("border_sapphire_Vacuum_Square_to_CapcitorH2Phys",SapphirePhys,Square_to_CapcitorH2Phys,fSiNbInterface);
   G4CMPLogicalBorderSurface * border_sapphire_Vacuum_Square_to_CapcitorH3Phys=new G4CMPLogicalBorderSurface("border_sapphire_Vacuum_Square_to_CapcitorH3Phys",SapphirePhys,Square_to_CapcitorH3Phys,fSiNbInterface);
G4CMPLogicalBorderSurface * border_sapphire_Vacuum_Square_to_CapcitorH4Phys=new G4CMPLogicalBorderSurface("border_sapphire_Vacuum_Square_to_CapcitorH4Phys",SapphirePhys,Square_to_CapcitorH4Phys,fSiNbInterface);
G4CMPLogicalBorderSurface * border_sapphire_Vacuum_Half_Cirlce_R_SquareCapacitorPhys2=new G4CMPLogicalBorderSurface("border_sapphire_Vacuum_Half_Cirlce_R_SquareCapacitorPhys2",SapphirePhys,Half_Cirlce_R_SquareCapacitorPhys2,fSiNbInterface);
   G4CMPLogicalBorderSurface * border_sapphire_Vacuum_Half_Cirlce_R_SquareCapacitorPhys3=new G4CMPLogicalBorderSurface("border_sapphire_Vacuum_Half_Cirlce_R_SquareCapacitorPhys3",SapphirePhys,Half_Cirlce_R_SquareCapacitorPhys3,fSiNbInterface);
   G4CMPLogicalBorderSurface * border_sapphire_Vacuum_Half_Cirlce_L_SquareCapacitorPhys4=new G4CMPLogicalBorderSurface("border_sapphire_Vacuum_Half_Cirlce_L_SquareCapacitorPhys4",SapphirePhys,Half_Cirlce_L_SquareCapacitorPhys4,fSiNbInterface);
   G4CMPLogicalBorderSurface * border_sapphire_Vacuum_Half_Cirlce_R_SquareCapacitorPhys1=new G4CMPLogicalBorderSurface("border_sapphire_Vacuum_Half_Cirlce_R_SquareCapacitorPhys1",SapphirePhys,Half_Cirlce_R_SquareCapacitorPhys1,fSiNbInterface);
   G4CMPLogicalBorderSurface * border_sapphire_Vacuum_Square_to_CapcitorV1Phys1=new G4CMPLogicalBorderSurface("border_sapphire_Vacuum_Square_to_CapcitorV1Phys1",SapphirePhys,Square_to_CapcitorV1Phys1,fSiNbInterface);
  G4CMPLogicalBorderSurface * border_sapphire_Vacuum_Half_Cirlce_L_TTDPhys_1=new G4CMPLogicalBorderSurface("border_sapphire_Vacuum_Half_Cirlce_L_TTDPhys_1",SapphirePhys,Half_Cirlce_L_TTDPhys_1,fSiNbInterface);
  G4CMPLogicalBorderSurface * border_sapphire_Vacuum_Square_to_CapcitorV1Phys2=new G4CMPLogicalBorderSurface("border_sapphire_Vacuum_Square_to_CapcitorV1Phys2",SapphirePhys,Square_to_CapcitorV1Phys2,fSiNbInterface);
   G4CMPLogicalBorderSurface * border_sapphire_Vacuum_Half_Cirlce_R_TTDPhys_1=new G4CMPLogicalBorderSurface("border_sapphire_Vacuum_Half_Cirlce_R_TTDPhys_1",SapphirePhys,Half_Cirlce_R_TTDPhys_1,fSiNbInterface);
  G4CMPLogicalBorderSurface * border_sapphire_Vacuum_Half_Cirlce_R_TTDPhys_2=new G4CMPLogicalBorderSurface("border_sapphire_Vacuum_Half_Cirlce_R_TTDPhys_2",SapphirePhys,Half_Cirlce_R_TTDPhys_2,fSiNbInterface);
   G4CMPLogicalBorderSurface * border_sapphire_Vacuum_Half_Cirlce_R_TTDPhys_3=new G4CMPLogicalBorderSurface("border_sapphire_Vacuum_Half_Cirlce_R_TTDPhys_3",SapphirePhys,Half_Cirlce_R_TTDPhys_3,fSiNbInterface);
 G4CMPLogicalBorderSurface * border_sapphire_Vacuum_Square_to_Qubit1Phys=new G4CMPLogicalBorderSurface("border_sapphire_Vacuum_Square_to_Qubit1Phys",SapphirePhys,Square_to_Qubit1Phys,fSiNbInterface);
   G4CMPLogicalBorderSurface * border_sapphire_Vacuum_Square_to_Qubit2Phys=new G4CMPLogicalBorderSurface("border_sapphire_Vacuum_Square_to_Qubit2Phys",SapphirePhys,Square_to_Qubit2Phys,fSiNbInterface);
   G4CMPLogicalBorderSurface * border_sapphire_Vacuum_Square_to_Qubit3Phys=new G4CMPLogicalBorderSurface("border_sapphire_Vacuum_Square_to_Qubit3Phys",SapphirePhys,Square_to_Qubit3Phys,fSiNbInterface);
   G4CMPLogicalBorderSurface * border_sapphire_Vacuum_Square_to_Qubit4Phys=new G4CMPLogicalBorderSurface("border_sapphire_Vacuum_Square_to_Qubit4Phys",SapphirePhys,Square_to_Qubit4Phys,fSiNbInterface);
   G4CMPLogicalBorderSurface * border_sapphire_Vacuum_Half_Cirlce_R_TTDPhys_Q_1=new G4CMPLogicalBorderSurface("border_sapphire_Vacuum_Half_Cirlce_R_TTDPhys_Q_1",SapphirePhys,Half_Cirlce_R_TTDPhys_Q_1,fSiNbInterface);
   G4CMPLogicalBorderSurface * border_sapphire_Vacuum_Half_Cirlce_R_TTDPhys_Q_2=new G4CMPLogicalBorderSurface("border_sapphire_Vacuum_Half_Cirlce_R_TTDPhys_Q_2",SapphirePhys,Half_Cirlce_R_TTDPhys_Q_2,fSiNbInterface);
   G4CMPLogicalBorderSurface * border_sapphire_Vacuum_Half_Cirlce_R_TTDPhys_Q_3=new G4CMPLogicalBorderSurface("border_sapphire_Vacuum_Half_Cirlce_R_TTDPhys_Q_3",SapphirePhys,Half_Cirlce_R_TTDPhys_Q_3,fSiNbInterface);
   G4CMPLogicalBorderSurface * border_sapphire_Vacuum_Half_Cirlce_R_TTDPhys_Q_4=new G4CMPLogicalBorderSurface("border_sapphire_Vacuum_Half_Cirlce_R_TTDPhys_Q_4",SapphirePhys,Half_Cirlce_R_TTDPhys_Q_4,fSiNbInterface);
   G4CMPLogicalBorderSurface * border_sapphire_Vacuum_SquareV_to_Qubit1Phys=new G4CMPLogicalBorderSurface("border_sapphire_Vacuum_SquareV_to_Qubit1Phys",SapphirePhys,SquareV_to_Qubit1Phys,fSiNbInterface);
   G4CMPLogicalBorderSurface * border_sapphire_Vacuum_SquareV_to_Qubit2Phys=new G4CMPLogicalBorderSurface("border_sapphire_Vacuum_SquareV_to_Qubit2Phys",SapphirePhys,SquareV_to_Qubit2Phys,fSiNbInterface);
   G4CMPLogicalBorderSurface * border_sapphire_Vacuum_SquareV_to_Qubit3Phys=new G4CMPLogicalBorderSurface("border_sapphire_Vacuum_SquareV_to_Qubit3Phys",SapphirePhys,SquareV_to_Qubit3Phys,fSiNbInterface);
   G4CMPLogicalBorderSurface * border_sapphire_Vacuum_SquareV_to_Qubit4Phys=new G4CMPLogicalBorderSurface("border_sapphire_Vacuum_SquareV_to_Qubit4Phys",SapphirePhys,SquareV_to_Qubit4Phys,fSiNbInterface);
G4CMPLogicalBorderSurface * border_sapphire_Vacuum_Square_Capacitor_QubitPhys1=new G4CMPLogicalBorderSurface("border_sapphire_Vacuum_Square_Capacitor_QubitPhys1",SapphirePhys,Square_Capacitor_QubitPhys1,fSiNbInterface);
   G4CMPLogicalBorderSurface * border_sapphire_Vacuum_Square_Capacitor_QubitPhys2=new G4CMPLogicalBorderSurface("border_sapphire_Vacuum_Square_Capacitor_QubitPhys2",SapphirePhys,Square_Capacitor_QubitPhys2,fSiNbInterface);
G4CMPLogicalBorderSurface * border_sapphire_Vacuum_Square_Capacitor_QubitPhys3=new G4CMPLogicalBorderSurface("border_sapphire_Vacuum_Square_Capacitor_QubitPhys3",SapphirePhys,Square_Capacitor_QubitPhys3,fSiNbInterface);
   G4CMPLogicalBorderSurface * border_sapphire_Vacuum_Square_Capacitor_QubitPhys4=new G4CMPLogicalBorderSurface("border_sapphire_Vacuum_Square_Capacitor_QubitPhys4",SapphirePhys,Square_Capacitor_QubitPhys4,fSiNbInterface);
   G4CMPLogicalBorderSurface * border_sapphire_Vacuum_Square_Capacitor_Q_To_P_VS1Phys=new G4CMPLogicalBorderSurface("border_sapphire_Vacuum_Square_Capacitor_Q_To_P_VS1Phys",SapphirePhys,Square_Capacitor_Q_To_P_VS1Phys,fSiNbInterface);
   G4CMPLogicalBorderSurface * border_sapphire_Vacuum_Square_Capacitor_Q_To_P_VS2Phys=new G4CMPLogicalBorderSurface("border_sapphire_Vacuum_Square_Capacitor_Q_To_P_VS2Phys",SapphirePhys,Square_Capacitor_Q_To_P_VS2Phys,fSiNbInterface);
   G4CMPLogicalBorderSurface * border_sapphire_Vacuum_Square_Capacitor_Q_To_P_VS3Phys=new G4CMPLogicalBorderSurface("border_sapphire_Vacuum_Square_Capacitor_Q_To_P_VS3Phys",SapphirePhys,Square_Capacitor_Q_To_P_VS3Phys,fSiNbInterface);
   G4CMPLogicalBorderSurface * border_sapphire_Vacuum_Square_Capacitor_Q_To_P_VS4Phys=new G4CMPLogicalBorderSurface("border_sapphire_Vacuum_Square_Capacitor_Q_To_P_VS4Phys",SapphirePhys,Square_Capacitor_Q_To_P_VS4Phys,fSiNbInterface);
   G4CMPLogicalBorderSurface * border_sapphire_Vacuum_Half_Cirlce_C_To_PPhys=new G4CMPLogicalBorderSurface("border_sapphire_Vacuum_Half_Cirlce_C_To_PPhys",SapphirePhys,Half_Cirlce_C_To_PPhys,fSiNbInterface);
   G4CMPLogicalBorderSurface * border_sapphire_Vacuum_Half_Cirlce_C_To_PPhys2=new G4CMPLogicalBorderSurface("border_sapphire_Vacuum_Half_Cirlce_C_To_PPhys2",SapphirePhys,Half_Cirlce_C_To_PPhys2,fSiNbInterface);
 G4CMPLogicalBorderSurface * border_sapphire_Vacuum_Half_Cirlce_C_To_PPhys3=new G4CMPLogicalBorderSurface("border_sapphire_Vacuum_Half_Cirlce_C_To_PPhys3",SapphirePhys,Half_Cirlce_C_To_PPhys3,fSiNbInterface);
   G4CMPLogicalBorderSurface * border_sapphire_Vacuum_Half_Cirlce_C_To_PPhys4=new G4CMPLogicalBorderSurface("border_sapphire_Vacuum_Half_Cirlce_C_To_PPhys4",SapphirePhys,Half_Cirlce_C_To_PPhys4,fSiNbInterface);
 G4CMPLogicalBorderSurface * border_sapphire_Vacuum_Half_Cirlce_C_To_PPhys_P_1=new G4CMPLogicalBorderSurface("border_sapphire_Vacuum_Half_Cirlce_C_To_PPhys_P_1",SapphirePhys,Half_Cirlce_C_To_PPhys_P_1,fSiNbInterface);
G4CMPLogicalBorderSurface * border_sapphire_Vacuum_Half_Cirlce_C_To_PPhys_P_2=new G4CMPLogicalBorderSurface("border_sapphire_Vacuum_Half_Cirlce_C_To_PPhys_P_2",SapphirePhys,Half_Cirlce_C_To_PPhys_P_2,fSiNbInterface);
   G4CMPLogicalBorderSurface * border_sapphire_Vacuum_Half_Cirlce_C_To_PPhys_P_3=new G4CMPLogicalBorderSurface("border_sapphire_Vacuum_Half_Cirlce_C_To_PPhys_P_3",SapphirePhys,Half_Cirlce_C_To_PPhys_P_3,fSiNbInterface);
 G4CMPLogicalBorderSurface * border_sapphire_Vacuum_Half_Cirlce_C_To_PPhys_P_4=new G4CMPLogicalBorderSurface("border_sapphire_Vacuum_Half_Cirlce_C_To_PPhys_P_4",SapphirePhys,Half_Cirlce_C_To_PPhys_P_4,fSiNbInterface);
   G4CMPLogicalBorderSurface * border_sapphire_Vacuum_Square_Capacitor_Q_To_P_HS1Phys=new G4CMPLogicalBorderSurface("border_sapphire_Vacuum_Square_Capacitor_Q_To_P_HS1Phys",SapphirePhys,Square_Capacitor_Q_To_P_HS1Phys,fSiNbInterface);
   G4CMPLogicalBorderSurface * border_sapphire_Vacuum_Square_Capacitor_Q_To_P_HS2Phys=new G4CMPLogicalBorderSurface("border_sapphire_Vacuum_Square_Capacitor_Q_To_P_HS2Phys",SapphirePhys,Square_Capacitor_Q_To_P_HS2Phys,fSiNbInterface);
   G4CMPLogicalBorderSurface * border_sapphire_Vacuum_Square_Capacitor_Q_To_P_HS3Phys=new G4CMPLogicalBorderSurface("border_sapphire_Vacuum_Square_Capacitor_Q_To_P_HS3Phys",SapphirePhys,Square_Capacitor_Q_To_P_HS3Phys,fSiNbInterface);
   G4CMPLogicalBorderSurface * border_sapphire_Vacuum_Square_Capacitor_Q_To_P_HS4Phys=new G4CMPLogicalBorderSurface("border_sapphire_Vacuum_Square_Capacitor_Q_To_P_HS4Phys",SapphirePhys,Square_Capacitor_Q_To_P_HS4Phys,fSiNbInterface);
   G4CMPLogicalBorderSurface * border_sapphire_Vacuum_Square_Capacitor_Q_To_P_V_P1Phys=new G4CMPLogicalBorderSurface("border_sapphire_Vacuum_Square_Capacitor_Q_To_P_V_P1Phys",SapphirePhys,Square_Capacitor_Q_To_P_V_P1Phys,fSiNbInterface);
   G4CMPLogicalBorderSurface * border_sapphire_Vacuum_Half_Cirlce_C_To_PPhys_P_5=new G4CMPLogicalBorderSurface("border_sapphire_Vacuum_Half_Cirlce_C_To_PPhys_P_5",SapphirePhys,Half_Cirlce_C_To_PPhys_P_5,fSiNbInterface);

G4CMPLogicalBorderSurface * border_sapphire_Vacuum_Copper_BoxPhys=new G4CMPLogicalBorderSurface("border_sapphire_Vacuum_Copper_BoxPhys",SapphirePhys,Copper_BoxPhys,fSiNbInterface);///
G4CMPLogicalBorderSurface * border_sapphire_Vacuum_Sonten_Copper1Phys1=new G4CMPLogicalBorderSurface("border_sapphire_Vacuum_Sonten_Copper1Phys1",SapphirePhys,Sonten_Copper1Phys1,fSiNbInterface);
G4CMPLogicalBorderSurface * border_sapphire_Vacuum_Sonten_Copper1Phys2=new G4CMPLogicalBorderSurface("border_sapphire_Vacuum_Sonten_Copper1Phys2",SapphirePhys,Sonten_Copper1Phys2,fSiNbInterface);
G4CMPLogicalBorderSurface * border_sapphire_Vacuum_Sonten_Copper1Phys3=new G4CMPLogicalBorderSurface("border_sapphire_Vacuum_Sonten_Copper1Phys3",SapphirePhys,Sonten_Copper1Phys3,fSiNbInterface);
G4CMPLogicalBorderSurface * border_sapphire_Vacuum_Sonten_Copper1Phys4=new G4CMPLogicalBorderSurface("border_sapphire_Vacuum_Sonten_Copper1Phys4",SapphirePhys,Sonten_Copper1Phys4,fSiNbInterface);





///////////-------------_Adiing the Vaccum Interface
G4CMPLogicalBorderSurface * border_sapphire_Vacuum_Resonator_TransmissionLine_GP=new G4CMPLogicalBorderSurface("border_sapphire_Vacuum_Resonator_TransmissionLine_GP",SapphirePhys,Resonator_TransmissionLine_GP,fSiVacuumInterface);
G4CMPLogicalBorderSurface * border_sapphire_Vacuum_Half_Cirlce_RPhys_Transmission_GP=new G4CMPLogicalBorderSurface("border_sapphire_Vacuum_Half_Cirlce_RPhys_Transmission_GP",SapphirePhys,Half_Cirlce_RPhys_Transmission_GP,fSiVacuumInterface);
//G4CMPLogicalBorderSurface * border_sapphire_Vacuum_Half_Cirlce_LPhys_Transmission_GP=new G4CMPLogicalBorderSurface("border_sapphire_Vacuum_Half_Cirlce_LPhys_Transmission_GP",SapphirePhys,Half_Cirlce_LPhys_Transmission_GP,fSiVacuumInterface);
G4CMPLogicalBorderSurface * border_sapphire_Vacuum_Quarter_Cirlce_RPhys_Transmission_GP=new G4CMPLogicalBorderSurface("border_sapphire_Vacuum_Quarter_Cirlce_RPhys_Transmission_GP",SapphirePhys,Quarter_Cirlce_RPhys_Transmission_GP,fSiVacuumInterface);
G4CMPLogicalBorderSurface * border_sapphire_Vacuum_Quarter_Cirlce_LPhys_Transmission_GP=new G4CMPLogicalBorderSurface("border_sapphire_Vacuum_Quarter_Cirlce_LPhys_Transmission_GP",SapphirePhys,Quarter_Cirlce_LPhys_Transmission_GP,fSiVacuumInterface);
G4CMPLogicalBorderSurface * border_sapphire_Vacuum_Resonator_TransmissionLine_1_GP=new G4CMPLogicalBorderSurface("border_sapphire_Vacuum_Resonator_TransmissionLine_1_GP",SapphirePhys,Resonator_TransmissionLine_1_GP,fSiVacuumInterface);
G4CMPLogicalBorderSurface * border_sapphire_Vacuum_Resonator_TransmissionLine_2_GP=new G4CMPLogicalBorderSurface("border_sapphire_Vacuum_Resonator_TransmissionLine_2_GP",SapphirePhys,Resonator_TransmissionLine_2_GP,fSiVacuumInterface);
G4CMPLogicalBorderSurface * border_sapphire_Vacuum_Square_To_ResonatorLogicalPhys1_GP=new G4CMPLogicalBorderSurface("border_sapphire_Vacuum_Square_To_ResonatorLogicalPhys1_GP",SapphirePhys,Square_To_ResonatorLogicalPhys1_GP,fSiVacuumInterface);
G4CMPLogicalBorderSurface * border_sapphire_Vacuum_Square_To_ResonatorLogicalPhys2_GP=new G4CMPLogicalBorderSurface("border_sapphire_Vacuum_Square_To_ResonatorLogicalPhys2_GP",SapphirePhys,Square_To_ResonatorLogicalPhys2_GP,fSiVacuumInterface);
G4CMPLogicalBorderSurface * border_sapphire_Vacuum_Square_To_ResonatorLogicalPhys3_GP=new G4CMPLogicalBorderSurface("border_sapphire_Vacuum_Square_To_ResonatorLogicalPhys3_GP",SapphirePhys,Square_To_ResonatorLogicalPhys3_GP,fSiVacuumInterface);
G4CMPLogicalBorderSurface * border_sapphire_Vacuum_Square_To_ResonatorLogicalPhys4_GP=new G4CMPLogicalBorderSurface("border_sapphire_Vacuum_Square_To_ResonatorLogicalPhys4_GP",SapphirePhys,Square_To_ResonatorLogicalPhys4_GP,fSiVacuumInterface);
G4CMPLogicalBorderSurface * border_sapphire_Vacuum_Square_to_CapcitorH1Phys_GP=new G4CMPLogicalBorderSurface("border_sapphire_Vacuum_Square_to_CapcitorH1Phys_GP",SapphirePhys,Square_to_CapcitorH1Phys_GP,fSiVacuumInterface);
G4CMPLogicalBorderSurface * border_sapphire_Vacuum_Square_to_CapcitorH2Phys_GP=new G4CMPLogicalBorderSurface("border_sapphire_Vacuum_Square_to_CapcitorH2Phys_GP",SapphirePhys,Square_to_CapcitorH2Phys_GP,fSiVacuumInterface);
G4CMPLogicalBorderSurface * border_sapphire_Vacuum_Square_to_CapcitorH3Phys_GP=new G4CMPLogicalBorderSurface("border_sapphire_Vacuum_Square_to_CapcitorH3Phys_GP",SapphirePhys,Square_to_CapcitorH3Phys_GP,fSiVacuumInterface);
G4CMPLogicalBorderSurface * border_sapphire_Vacuum_Square_to_CapcitorH4Phys_GP=new G4CMPLogicalBorderSurface("border_sapphire_Vacuum_Square_to_CapcitorH4Phys_GP",SapphirePhys,Square_to_CapcitorH4Phys_GP,fSiVacuumInterface);
G4CMPLogicalBorderSurface * border_sapphire_Vacuum_Half_Cirlce_R_SquareCapacitorPhys2_GP=new G4CMPLogicalBorderSurface("border_sapphire_Vacuum_Half_Cirlce_R_SquareCapacitorPhys2_GP",SapphirePhys,Half_Cirlce_R_SquareCapacitorPhys2_GP,fSiVacuumInterface);
G4CMPLogicalBorderSurface * border_sapphire_Vacuum_Half_Cirlce_R_SquareCapacitorPhys3_GP=new G4CMPLogicalBorderSurface("border_sapphire_Vacuum_Half_Cirlce_R_SquareCapacitorPhys3_GP",SapphirePhys,Half_Cirlce_R_SquareCapacitorPhys3_GP,fSiVacuumInterface);
G4CMPLogicalBorderSurface * border_sapphire_Vacuum_Half_Cirlce_L_SquareCapacitorPhys4_GP=new G4CMPLogicalBorderSurface("border_sapphire_Vacuum_Half_Cirlce_L_SquareCapacitorPhys4_GP",SapphirePhys,Half_Cirlce_L_SquareCapacitorPhys4_GP,fSiVacuumInterface);
G4CMPLogicalBorderSurface * border_sapphire_Vacuum_Half_Cirlce_R_SquareCapacitorPhys1_GP=new G4CMPLogicalBorderSurface("border_sapphire_Vacuum_Half_Cirlce_R_SquareCapacitorPhys1_GP",SapphirePhys,Half_Cirlce_R_SquareCapacitorPhys1_GP,fSiVacuumInterface);
G4CMPLogicalBorderSurface * border_sapphire_Vacuum_Square_to_CapcitorV1Phys1_GP=new G4CMPLogicalBorderSurface("border_sapphire_Vacuum_Square_to_CapcitorV1Phys1_GP",SapphirePhys,Square_to_CapcitorV1Phys1_GP,fSiVacuumInterface);
G4CMPLogicalBorderSurface * border_sapphire_Vacuum_Half_Cirlce_L_TTDPhys_1_GP=new G4CMPLogicalBorderSurface("border_sapphire_Vacuum_Half_Cirlce_L_TTDPhys_1_GP",SapphirePhys,Half_Cirlce_L_TTDPhys_1_GP,fSiVacuumInterface);
G4CMPLogicalBorderSurface * border_sapphire_Vacuum_Square_to_CapcitorV1Phys2_GP=new G4CMPLogicalBorderSurface("border_sapphire_Vacuum_Square_to_CapcitorV1Phys2_GP",SapphirePhys,Square_to_CapcitorV1Phys2_GP,fSiVacuumInterface);
G4CMPLogicalBorderSurface * border_sapphire_Vacuum_Half_Cirlce_R_TTDPhys_1_GP=new G4CMPLogicalBorderSurface("border_sapphire_Vacuum_Half_Cirlce_R_TTDPhys_1_GP",SapphirePhys,Half_Cirlce_R_TTDPhys_1_GP,fSiVacuumInterface);
G4CMPLogicalBorderSurface * border_sapphire_Vacuum_Half_Cirlce_R_TTDPhys_2_GP=new G4CMPLogicalBorderSurface("border_sapphire_Vacuum_Half_Cirlce_R_TTDPhys_2_GP",SapphirePhys,Half_Cirlce_R_TTDPhys_2_GP,fSiVacuumInterface);
G4CMPLogicalBorderSurface * border_sapphire_Vacuum_Half_Cirlce_R_TTDPhys_3_GP=new G4CMPLogicalBorderSurface("border_sapphire_Vacuum_Half_Cirlce_R_TTDPhys_3_GP",SapphirePhys,Half_Cirlce_R_TTDPhys_3_GP,fSiVacuumInterface);
G4CMPLogicalBorderSurface * border_sapphire_Vacuum_Square_to_Qubit1Phys_GP=new G4CMPLogicalBorderSurface("border_sapphire_Vacuum_Square_to_Qubit1Phys_GP",SapphirePhys,Square_to_Qubit1Phys_GP,fSiVacuumInterface);
G4CMPLogicalBorderSurface * border_sapphire_Vacuum_Square_to_Qubit2Phys_GP=new G4CMPLogicalBorderSurface("border_sapphire_Vacuum_Square_to_Qubit2Phys_GP",SapphirePhys,Square_to_Qubit2Phys_GP,fSiVacuumInterface);
G4CMPLogicalBorderSurface * border_sapphire_Vacuum_Square_to_Qubit3Phys_GP=new G4CMPLogicalBorderSurface("border_sapphire_Vacuum_Square_to_Qubit3Phys_GP",SapphirePhys,Square_to_Qubit3Phys_GP,fSiVacuumInterface);
G4CMPLogicalBorderSurface * border_sapphire_Vacuum_Square_to_Qubit4Phys_GP=new G4CMPLogicalBorderSurface("border_sapphire_Vacuum_Square_to_Qubit4Phys_GP",SapphirePhys,Square_to_Qubit4Phys_GP,fSiVacuumInterface);
G4CMPLogicalBorderSurface * border_sapphire_Vacuum_Half_Cirlce_R_TTDPhys_Q_1_GP=new G4CMPLogicalBorderSurface("border_sapphire_Vacuum_Half_Cirlce_R_TTDPhys_Q_1_GP",SapphirePhys,Half_Cirlce_R_TTDPhys_Q_1_GP,fSiVacuumInterface);
G4CMPLogicalBorderSurface * border_sapphire_Vacuum_Half_Cirlce_R_TTDPhys_Q_2_GP=new G4CMPLogicalBorderSurface("border_sapphire_Vacuum_Half_Cirlce_R_TTDPhys_Q_2_GP",SapphirePhys,Half_Cirlce_R_TTDPhys_Q_2_GP,fSiVacuumInterface);
G4CMPLogicalBorderSurface * border_sapphire_Vacuum_Half_Cirlce_R_TTDPhys_Q_3_GP=new G4CMPLogicalBorderSurface("border_sapphire_Vacuum_Half_Cirlce_R_TTDPhys_Q_3_GP",SapphirePhys,Half_Cirlce_R_TTDPhys_Q_3_GP,fSiVacuumInterface);
G4CMPLogicalBorderSurface * border_sapphire_Vacuum_Half_Cirlce_R_TTDPhys_Q_4_GP=new G4CMPLogicalBorderSurface("border_sapphire_Vacuum_Half_Cirlce_R_TTDPhys_Q_4_GP",SapphirePhys,Half_Cirlce_R_TTDPhys_Q_4_GP,fSiVacuumInterface);
G4CMPLogicalBorderSurface * border_sapphire_Vacuum_SquareV_to_Qubit1Phys_GP=new G4CMPLogicalBorderSurface("border_sapphire_Vacuum_SquareV_to_Qubit1Phys_GP",SapphirePhys,SquareV_to_Qubit1Phys_GP,fSiVacuumInterface);
G4CMPLogicalBorderSurface * border_sapphire_Vacuum_SquareV_to_Qubit2Phys_GP=new G4CMPLogicalBorderSurface("border_sapphire_Vacuum_SquareV_to_Qubit2Phys_GP",SapphirePhys,SquareV_to_Qubit2Phys_GP,fSiVacuumInterface);
G4CMPLogicalBorderSurface * border_sapphire_Vacuum_SquareV_to_Qubit3Phys_GP=new G4CMPLogicalBorderSurface("border_sapphire_Vacuum_SquareV_to_Qubit3Phys_GP",SapphirePhys,SquareV_to_Qubit3Phys_GP,fSiVacuumInterface);
G4CMPLogicalBorderSurface * border_sapphire_Vacuum_SquareV_to_Qubit4Phys_GP=new G4CMPLogicalBorderSurface("border_sapphire_Vacuum_SquareV_to_Qubit4Phys_GP",SapphirePhys,SquareV_to_Qubit4Phys_GP,fSiVacuumInterface);
G4CMPLogicalBorderSurface * border_sapphire_Vacuum_Square_Capacitor_QubitPhys1_GP=new G4CMPLogicalBorderSurface("border_sapphire_Vacuum_Square_Capacitor_QubitPhys1_GP",SapphirePhys,Square_Capacitor_QubitPhys1_GP,fSiVacuumInterface);
G4CMPLogicalBorderSurface * border_sapphire_Vacuum_Square_Capacitor_QubitPhys2_GP=new G4CMPLogicalBorderSurface("border_sapphire_Vacuum_Square_Capacitor_QubitPhys2_GP",SapphirePhys,Square_Capacitor_QubitPhys2_GP,fSiVacuumInterface);
G4CMPLogicalBorderSurface * border_sapphire_Vacuum_Square_Capacitor_QubitPhys3_GP=new G4CMPLogicalBorderSurface("border_sapphire_Vacuum_Square_Capacitor_QubitPhys3_GP",SapphirePhys,Square_Capacitor_QubitPhys3_GP,fSiVacuumInterface);
G4CMPLogicalBorderSurface * border_sapphire_Vacuum_Square_Capacitor_QubitPhys4_GP=new G4CMPLogicalBorderSurface("border_sapphire_Vacuum_Square_Capacitor_QubitPhys4_GP",SapphirePhys,Square_Capacitor_QubitPhys4_GP,fSiVacuumInterface);
G4CMPLogicalBorderSurface * border_sapphire_Vacuum_Square_Capacitor_Q_To_P_VS1Phys_GP=new G4CMPLogicalBorderSurface("border_sapphire_Vacuum_Square_Capacitor_Q_To_P_VS1Phys_GP",SapphirePhys,Square_Capacitor_Q_To_P_VS1Phys_GP,fSiVacuumInterface);
G4CMPLogicalBorderSurface * border_sapphire_Vacuum_Square_Capacitor_Q_To_P_VS2Phys_GP=new G4CMPLogicalBorderSurface("border_sapphire_Vacuum_Square_Capacitor_Q_To_P_VS2Phys_GP",SapphirePhys,Square_Capacitor_Q_To_P_VS2Phys_GP,fSiVacuumInterface);
G4CMPLogicalBorderSurface * border_sapphire_Vacuum_Square_Capacitor_Q_To_P_VS3Phys_GP=new G4CMPLogicalBorderSurface("border_sapphire_Vacuum_Square_Capacitor_Q_To_P_VS3Phys_GP",SapphirePhys,Square_Capacitor_Q_To_P_VS3Phys_GP,fSiVacuumInterface);
G4CMPLogicalBorderSurface * border_sapphire_Vacuum_Square_Capacitor_Q_To_P_VS4Phys_GP=new G4CMPLogicalBorderSurface("border_sapphire_Vacuum_Square_Capacitor_Q_To_P_VS4Phys_GP",SapphirePhys,Square_Capacitor_Q_To_P_VS4Phys_GP,fSiVacuumInterface);
G4CMPLogicalBorderSurface * border_sapphire_Vacuum_Half_Cirlce_C_To_PPhys_GP=new G4CMPLogicalBorderSurface("border_sapphire_Vacuum_Half_Cirlce_C_To_PPhys_GP",SapphirePhys,Half_Cirlce_C_To_PPhys_GP,fSiVacuumInterface);
G4CMPLogicalBorderSurface * border_sapphire_Vacuum_Half_Cirlce_C_To_PPhys2_GP=new G4CMPLogicalBorderSurface("border_sapphire_Vacuum_Half_Cirlce_C_To_PPhys2_GP",SapphirePhys,Half_Cirlce_C_To_PPhys2_GP,fSiVacuumInterface);
G4CMPLogicalBorderSurface * border_sapphire_Vacuum_Half_Cirlce_C_To_PPhys3_GP=new G4CMPLogicalBorderSurface("border_sapphire_Vacuum_Half_Cirlce_C_To_PPhys3_GP",SapphirePhys,Half_Cirlce_C_To_PPhys3_GP,fSiVacuumInterface);
G4CMPLogicalBorderSurface * border_sapphire_Vacuum_Half_Cirlce_C_To_PPhys4_GP=new G4CMPLogicalBorderSurface("border_sapphire_Vacuum_Half_Cirlce_C_To_PPhys4_GP",SapphirePhys,Half_Cirlce_C_To_PPhys4_GP,fSiVacuumInterface);
G4CMPLogicalBorderSurface * border_sapphire_Vacuum_Half_Cirlce_C_To_PPhys_P_1_GP=new G4CMPLogicalBorderSurface("border_sapphire_Vacuum_Half_Cirlce_C_To_PPhys_P_1_GP",SapphirePhys,Half_Cirlce_C_To_PPhys_P_1_GP,fSiVacuumInterface);
G4CMPLogicalBorderSurface * border_sapphire_Vacuum_Half_Cirlce_C_To_PPhys_P_2_GP=new G4CMPLogicalBorderSurface("border_sapphire_Vacuum_Half_Cirlce_C_To_PPhys_P_2_GP",SapphirePhys,Half_Cirlce_C_To_PPhys_P_2_GP,fSiVacuumInterface);
G4CMPLogicalBorderSurface * border_sapphire_Vacuum_Half_Cirlce_C_To_PPhys_P_3_GP=new G4CMPLogicalBorderSurface("border_sapphire_Vacuum_Half_Cirlce_C_To_PPhys_P_3_GP",SapphirePhys,Half_Cirlce_C_To_PPhys_P_3_GP,fSiVacuumInterface);
G4CMPLogicalBorderSurface * border_sapphire_Vacuum_Half_Cirlce_C_To_PPhys_P_4_GP=new G4CMPLogicalBorderSurface("border_sapphire_Vacuum_Half_Cirlce_C_To_PPhys_P_4_GP",SapphirePhys,Half_Cirlce_C_To_PPhys_P_4_GP,fSiVacuumInterface);
G4CMPLogicalBorderSurface * border_sapphire_Vacuum_Square_Capacitor_Q_To_P_HS1Phys_GP=new G4CMPLogicalBorderSurface("border_sapphire_Vacuum_Square_Capacitor_Q_To_P_HS1Phys_GP",SapphirePhys,Square_Capacitor_Q_To_P_HS1Phys_GP,fSiVacuumInterface);
G4CMPLogicalBorderSurface * border_sapphire_Vacuum_Square_Capacitor_Q_To_P_HS2Phys_GP=new G4CMPLogicalBorderSurface("border_sapphire_Vacuum_Square_Capacitor_Q_To_P_HS2Phys_GP",SapphirePhys,Square_Capacitor_Q_To_P_HS2Phys_GP,fSiVacuumInterface);
G4CMPLogicalBorderSurface * border_sapphire_Vacuum_Square_Capacitor_Q_To_P_HS3Phys_GP=new G4CMPLogicalBorderSurface("border_sapphire_Vacuum_Square_Capacitor_Q_To_P_HS3Phys_GP",SapphirePhys,Square_Capacitor_Q_To_P_HS3Phys_GP,fSiVacuumInterface);
G4CMPLogicalBorderSurface * border_sapphire_Vacuum_Square_Capacitor_Q_To_P_HS4Phys_GP=new G4CMPLogicalBorderSurface("border_sapphire_Vacuum_Square_Capacitor_Q_To_P_HS4Phys_GP",SapphirePhys,Square_Capacitor_Q_To_P_HS4Phys_GP,fSiVacuumInterface);
G4CMPLogicalBorderSurface * border_sapphire_Vacuum_Square_Capacitor_Q_To_P_V_P1Phys_GP=new G4CMPLogicalBorderSurface("border_sapphire_Vacuum_Square_Capacitor_Q_To_P_V_P1Phys_GP",SapphirePhys,Square_Capacitor_Q_To_P_V_P1Phys_GP,fSiVacuumInterface);
G4CMPLogicalBorderSurface * border_sapphire_Vacuum_Half_Cirlce_C_To_PPhys_P_5_GP=new G4CMPLogicalBorderSurface("border_sapphire_Vacuum_Half_Cirlce_C_To_PPhys_P_5_GP",SapphirePhys,Half_Cirlce_C_To_PPhys_P_5_GP,fSiVacuumInterface);
//G4CMPLogicalBorderSurface * border_sapphire_Vacuum_Copper_BoxPhys_GP=new G4CMPLogicalBorderSurface("border_sapphire_Vacuum_Copper_BoxPhys_GP",SapphirePhys,Copper_BoxPhys_GP,fSiVacuumInterface);
//G4CMPLogicalBorderSurface * border_sapphire_Vacuum_Sonten_Copper1Phys1_GP=new G4CMPLogicalBorderSurface("border_sapphire_Vacuum_Sonten_Copper1Phys1_GP",SapphirePhys,Sonten_Copper1Phys1_GP,fSiVacuumInterface);
//G4CMPLogicalBorderSurface * border_sapphire_Vacuum_Sonten_Copper1Phys2_GP=new G4CMPLogicalBorderSurface("border_sapphire_Vacuum_Sonten_Copper1Phys2_GP",SapphirePhys,Sonten_Copper1Phys2_GP,fSiVacuumInterface);
//G4CMPLogicalBorderSurface * border_sapphire_Vacuum_Sonten_Copper1Phys3_GP=new G4CMPLogicalBorderSurface("border_sapphire_Vacuum_Sonten_Copper1Phys3_GP",SapphirePhys,Sonten_Copper1Phys3_GP,fSiVacuumInterface);
//G4CMPLogicalBorderSurface * border_sapphire_Vacuum_Sonten_Copper1Phys4_GP=new G4CMPLogicalBorderSurface("border_sapphire_Vacuum_Sonten_Copper1Phys4_GP",SapphirePhys,Sonten_Copper1Phys4_GP,fSiVacuumInterface);
////////////?Physical Volumes
//G4VPhysicalVolume* SapphirePhys

//G4CMPLogicalBorderSurface * border_sapphire_Vacuum_Resonator_TransmissionLine= new G4CMPLogicalBorderSurface("border_sapphire_Vacuum_Resonator_TransmissionLine", SapphirePhys,Resonator_TransmissionLine ,fSiVacuumInterface);
//
// G4VPhysicalVolume* Resonator_TransmissionLine
// G4VPhysicalVolume* Half_Cirlce_RPhys_Transmission
// G4VPhysicalVolume* Half_Cirlce_RPhys_Transmission
// G4VPhysicalVolume* Half_Cirlce_LPhys_Transmission
// G4VPhysicalVolume* Quarter_Cirlce_RPhys_Transmission
// G4VPhysicalVolume* Quarter_Cirlce_LPhys_Transmission
// G4VPhysicalVolume* Resonator_TransmissionLine_1
// G4VPhysicalVolume* Resonator_TransmissionLine_2
// G4VPhysicalVolume* Square_Capacitor_to_TTLDPhys
// G4VPhysicalVolume* Half_Cirlce_R_Transmission_Square_Capacitor_to_TTDLPhys
// G4VPhysicalVolume* Square_connection_from_TTLD_to_ResonatorPhys
// G4VPhysicalVolume* Half_Cirlce_L_Transmission_Square_Capacitor_to_TTDLPhys
// G4VPhysicalVolume* Half_Cirlce_L_TTDPhys
// G4VPhysicalVolume* Half_Cirlce_L_TTDPhys
// G4VPhysicalVolume* Half_Cirlce_R_TTDPhys
// G4VPhysicalVolume* Square_To_ResonatorLogicalPhys1
// G4VPhysicalVolume* Square_To_ResonatorLogicalPhys2
// G4VPhysicalVolume* Square_To_ResonatorLogicalPhys3
// G4VPhysicalVolume* Square_To_ResonatorLogicalPhys4
// G4VPhysicalVolume* Resonator_SquarePhys_1
// G4VPhysicalVolume* Half_Cirlce_RPhys_1
// G4VPhysicalVolume* Half_Cirlce_LPhys_1
// G4VPhysicalVolume* Resonator_SquarePhys_2
// G4VPhysicalVolume* Half_Cirlce_RPhys_2
// G4VPhysicalVolume* Half_Cirlce_LPhys_2
// G4VPhysicalVolume* Resonator_SquarePhys_3
// G4VPhysicalVolume* Half_Cirlce_RPhys_3
// G4VPhysicalVolume* Half_Cirlce_LPhys_3
// G4VPhysicalVolume* Resonator_SquarePhys_4
// G4VPhysicalVolume* Half_Cirlce_RPhys_4
// G4VPhysicalVolume* Half_Cirlce_LPhys_4
// G4VPhysicalVolume* Square_to_CapcitorH1Phys
// G4VPhysicalVolume* Square_to_CapcitorH2Phys
// G4VPhysicalVolume* Square_to_CapcitorH3Phys
// G4VPhysicalVolume* Square_to_CapcitorH4Phys
// G4VPhysicalVolume* Half_Cirlce_R_SquareCapacitorPhys2
// G4VPhysicalVolume* Half_Cirlce_R_SquareCapacitorPhys3
// G4VPhysicalVolume* Half_Cirlce_L_SquareCapacitorPhys4
// G4VPhysicalVolume* Half_Cirlce_R_SquareCapacitorPhys1
// G4VPhysicalVolume* Square_to_CapcitorV1Phys1
// G4VPhysicalVolume* Half_Cirlce_L_TTDPhys_1
// G4VPhysicalVolume* Square_to_CapcitorV1Phys2
// G4VPhysicalVolume* Half_Cirlce_R_TTDPhys_1
// G4VPhysicalVolume* Half_Cirlce_R_TTDPhys_2
// G4VPhysicalVolume* Half_Cirlce_R_TTDPhys_3
// G4VPhysicalVolume* Square_to_Qubit1Phys
// G4VPhysicalVolume* Square_to_Qubit2Phys
// G4VPhysicalVolume* Square_to_Qubit3Phys
// G4VPhysicalVolume* Square_to_Qubit4Phys
// G4VPhysicalVolume* Half_Cirlce_R_TTDPhys_Q_1
// G4VPhysicalVolume* Half_Cirlce_R_TTDPhys_Q_2
// G4VPhysicalVolume* Half_Cirlce_R_TTDPhys_Q_3
// G4VPhysicalVolume* Half_Cirlce_R_TTDPhys_Q_4
// G4VPhysicalVolume* SquareV_to_Qubit1Phys
// G4VPhysicalVolume* SquareV_to_Qubit2Phys
// G4VPhysicalVolume* SquareV_to_Qubit3Phys
// G4VPhysicalVolume* SquareV_to_Qubit4Phys
// G4VPhysicalVolume* Square_Capacitor_QubitPhys1
// G4VPhysicalVolume* Square_Capacitor_QubitPhys2
// G4VPhysicalVolume* Square_Capacitor_QubitPhys3
// G4VPhysicalVolume* Cross_QubitPhys
// G4VPhysicalVolume* Square_Capacitor_QubitPhys4
// G4VPhysicalVolume* Square_Capacitor_Q_To_P_VS1Phys
// G4VPhysicalVolume* Square_Capacitor_Q_To_P_VS2Phys
// G4VPhysicalVolume* Square_Capacitor_Q_To_P_VS3Phys
// G4VPhysicalVolume* Square_Capacitor_Q_To_P_VS4Phys
// G4VPhysicalVolume* Half_Cirlce_C_To_PPhys
// G4VPhysicalVolume* Half_Cirlce_C_To_PPhys2
// G4VPhysicalVolume* Half_Cirlce_C_To_PPhys3
// G4VPhysicalVolume* Half_Cirlce_C_To_PPhys4
// G4VPhysicalVolume* Half_Cirlce_C_To_PPhys_P_1
// G4VPhysicalVolume* Half_Cirlce_C_To_PPhys_P_2
// G4VPhysicalVolume* Half_Cirlce_C_To_PPhys_P_3
// G4VPhysicalVolume* Half_Cirlce_C_To_PPhys_P_4
// G4VPhysicalVolume* Square_Capacitor_Q_To_P_HS1Phys
// G4VPhysicalVolume* Square_Capacitor_Q_To_P_HS2Phys
// G4VPhysicalVolume* Square_Capacitor_Q_To_P_HS3Phys
// G4VPhysicalVolume* Square_Capacitor_Q_To_P_HS4Phys
// G4VPhysicalVolume* Square_Capacitor_Q_To_P_V_P1Phys
// G4VPhysicalVolume* Half_Cirlce_C_To_PPhys_P_5
// G4VPhysicalVolume* Big_Trian_PanelPhys
// G4VPhysicalVolume* Square_PanelPhys
// G4VPhysicalVolume* Big_Trian_PanelPhys
// G4VPhysicalVolume*  Square_PanelPhys
// G4VPhysicalVolume* Big_Trian_PanelPhys
// G4VPhysicalVolume* Square_PanelPhys
// G4VPhysicalVolume* Big_Trian_PanelPhys
// G4VPhysicalVolume* Square_PanelPhys
// G4VPhysicalVolume* Copper_BoxPhys
// G4VPhysicalVolume* Sonten_Copper1Phys1
// G4VPhysicalVolume* Sonten_Copper1Phys2
// G4VPhysicalVolume* Sonten_Copper1Phys3
// G4VPhysicalVolume* Sonten_Copper1Phys4
//
// /////Vacumm Interface
//
//
// G4VPhysicalVolume* Resonator_TransmissionLine_GP
// G4VPhysicalVolume* Half_Cirlce_RPhys_Transmission_GP
// G4VPhysicalVolume* Half_Cirlce_LPhys_Transmission_GP
// G4VPhysicalVolume* Quarter_Cirlce_RPhys_Transmission_GP
// G4VPhysicalVolume* Quarter_Cirlce_LPhys_Transmission_GP
// G4VPhysicalVolume* Resonator_TransmissionLine_1_GP
// G4VPhysicalVolume* Resonator_TransmissionLine_2_GP
// G4VPhysicalVolume* Square_Capacitor_to_TTLDPhys_GP
// G4VPhysicalVolume* Half_Cirlce_R_Transmission_Square_Capacitor_to_TTDLPhys_GP
// G4VPhysicalVolume* Half_Cirlce_L_Transmission_Square_Capacitor_to_TTDLPhys_GP
// G4VPhysicalVolume* Square_connection_from_TTLD_to_ResonatorPhys_GP
// G4VPhysicalVolume* Square_connection_from_TTLD_to_ResonatorPhys_B_GP
// G4VPhysicalVolume* Square_connection_from_TTLD_to_ResonatorPhys_B
// G4VPhysicalVolume* Half_Cirlce_L_TTDPhys_GP
// G4VPhysicalVolume* Half_Cirlce_L_TTDPhys_GP
// G4VPhysicalVolume* Half_Cirlce_R_TTDPhys_GP
// G4VPhysicalVolume* Square_To_ResonatorLogicalPhys1_GP
// G4VPhysicalVolume* Square_To_ResonatorLogicalPhys2_GP
// G4VPhysicalVolume* Square_To_ResonatorLogicalPhys3_GP
// G4VPhysicalVolume* Square_To_ResonatorLogicalPhys4_GP
// G4VPhysicalVolume* Resonator_SquarePhys_1_GP
// G4VPhysicalVolume* Half_Cirlce_RPhys_1_GP
// G4VPhysicalVolume* Half_Cirlce_LPhys_1_GP
// G4VPhysicalVolume* Resonator_SquarePhys_2_GP
// G4VPhysicalVolume* Half_Cirlce_RPhys_2_GP
// G4VPhysicalVolume* Half_Cirlce_LPhys_2_GP
// G4VPhysicalVolume* Resonator_SquarePhys_3_GP
// G4VPhysicalVolume* Half_Cirlce_RPhys_3_GP
// G4VPhysicalVolume* Half_Cirlce_LPhys_3_GP
// G4VPhysicalVolume* Resonator_SquarePhys_4_GP
// G4VPhysicalVolume* Half_Cirlce_RPhys_4_GP
// G4VPhysicalVolume* Half_Cirlce_LPhys_4_GP
// G4VPhysicalVolume* Square_to_CapcitorH1Phys_GP
// G4VPhysicalVolume* Square_to_CapcitorH2Phys_GP
// G4VPhysicalVolume* Square_to_CapcitorH3Phys_GP
// G4VPhysicalVolume* Square_to_CapcitorH4Phys_GP
// G4VPhysicalVolume* Half_Cirlce_R_SquareCapacitorPhys1_GP
// G4VPhysicalVolume* Half_Cirlce_R_SquareCapacitorPhys2_GP
// G4VPhysicalVolume* Half_Cirlce_R_SquareCapacitorPhys3_GP
// G4VPhysicalVolume* Half_Cirlce_L_SquareCapacitorPhys4_GP
// G4VPhysicalVolume* Square_to_CapcitorV1Phys1_GP
// G4VPhysicalVolume* Half_Cirlce_L_TTDPhys_1_GP
// G4VPhysicalVolume* Square_to_CapcitorV1Phys2_GP
// G4VPhysicalVolume* Half_Cirlce_R_TTDPhys_1_GP
// G4VPhysicalVolume* Half_Cirlce_R_TTDPhys_2_GP
// G4VPhysicalVolume* Half_Cirlce_R_TTDPhys_3_GP
// G4VPhysicalVolume* Square_to_Qubit1Phys_GP
// G4VPhysicalVolume* Square_to_Qubit2Phys_GP
// G4VPhysicalVolume* Square_to_Qubit3Phys_GP
// G4VPhysicalVolume* Square_to_Qubit4Phys_GP
// G4VPhysicalVolume* Half_Cirlce_R_TTDPhys_Q_1_GP
// G4VPhysicalVolume* Half_Cirlce_R_TTDPhys_Q_2_GP
// G4VPhysicalVolume* Half_Cirlce_R_TTDPhys_Q_3_GP
// G4VPhysicalVolume* Half_Cirlce_R_TTDPhys_Q_4_GP
// G4VPhysicalVolume* SquareV_to_Qubit1Phys_GP
// G4VPhysicalVolume* SquareV_to_Qubit2Phys_GP
// G4VPhysicalVolume* SquareV_to_Qubit3Phys_GP
// G4VPhysicalVolume* SquareV_to_Qubit4Phys_GP
// G4VPhysicalVolume* Square_Capacitor_QubitPhys1_GP
// G4VPhysicalVolume* Square_Capacitor_QubitPhys2_GP
// G4VPhysicalVolume* Square_Capacitor_QubitPhys3_GP
// G4VPhysicalVolume* Square_Capacitor_QubitPhys4_GP
// G4VPhysicalVolume* Cross_QubitPhys_GP
// G4VPhysicalVolume* Square_Capacitor_Q_To_P_VS1Phys_GP
// G4VPhysicalVolume* Square_Capacitor_Q_To_P_VS2Phys_GP
// G4VPhysicalVolume* Square_Capacitor_Q_To_P_VS3Phys_GP
// G4VPhysicalVolume* Square_Capacitor_Q_To_P_VS4Phys_GP
// G4VPhysicalVolume* Half_Cirlce_C_To_PPhys_GP
// G4VPhysicalVolume* Half_Cirlce_C_To_PPhys2_GP
// G4VPhysicalVolume* Half_Cirlce_C_To_PPhys3_GP
// G4VPhysicalVolume* Half_Cirlce_C_To_PPhys4_GP
// G4VPhysicalVolume* Square_Capacitor_Q_To_P_HS1Phys_GP
// G4VPhysicalVolume* Square_Capacitor_Q_To_P_HS2Phys_GP
// G4VPhysicalVolume* Square_Capacitor_Q_To_P_HS3Phys_GP
// G4VPhysicalVolume* Square_Capacitor_Q_To_P_HS4Phys_GP
// G4VPhysicalVolume* Half_Cirlce_C_To_PPhys_P_1_GP
// G4VPhysicalVolume* Half_Cirlce_C_To_PPhys_P_2_GP
// G4VPhysicalVolume* Half_Cirlce_C_To_PPhys_P_3_GP
// G4VPhysicalVolume* Half_Cirlce_C_To_PPhys_P_4_GP
// G4VPhysicalVolume* Square_Capacitor_Q_To_P_V_P1Phys_GP
// G4VPhysicalVolume* Half_Cirlce_C_To_PPhys_P_5_GP
// G4VPhysicalVolume* Big_Trian_PanelPhys_GP
// G4VPhysicalVolume* Square_PanelPhys_GP
// G4VPhysicalVolume* Square_PanelPhys_GP1
// G4VPhysicalVolume* Big_Trian_PanelPhys_GP
// G4VPhysicalVolume* Square_PanelPhys_GP
// G4VPhysicalVolume* Square_PanelPhys_GP1
// G4VPhysicalVolume* Big_Trian_PanelPhys_GP
// G4VPhysicalVolume* Square_PanelPhys_GP
// G4VPhysicalVolume* Square_PanelPhys_GP1
// G4VPhysicalVolume* Big_Trian_PanelPhys_GP
// G4VPhysicalVolume* Square_PanelPhys_GP
// G4VPhysicalVolume* Square_PanelPhys_GP1


//Miust be the Physical Property.
//So I need to separate of





// G4LogicalVolume* fGermaniumLogical
// G4LogicalVolume* Resonator_SquareLogical_1
// G4LogicalVolume* Half_Cirlce_RLogical_1
// G4LogicalVolume* Half_Cirlce_LLogical_1
// G4LogicalVolume* Cross_QubitLogical
// G4LogicalVolume* Half_Cirlce_RLogical_2
// G4LogicalVolume* Resonator_SquareLogical_2
// G4LogicalVolume* Half_Cirlce_LLogical_2
// G4LogicalVolume* Resonator_SquareLogical_3
// G4LogicalVolume* Half_Cirlce_RLogical_3
// G4LogicalVolume* Half_Cirlce_LLogical_3
// G4LogicalVolume* Cross_QubitLogical_GP
// G4LogicalVolume* Half_Cirlce_RLogical_4
// G4LogicalVolume* Resonator_SquareLogical_4
// G4LogicalVolume* Half_Cirlce_LLogical_4
// G4LogicalVolume* TransmissionLine_Resonator_SquareLogical
// G4LogicalVolume* TransmissionLine_Resonator_Square_GPLogical
// G4LogicalVolume* Half_Cirlce_RLogical__Transmission
// G4LogicalVolume* Square_Capacitor_to_TTLDLogical
// G4LogicalVolume* Half_Cirlce_R_Transmission_Square_Capacitor_to_TTDLLogical
// G4LogicalVolume* Half_Cirlce_L_Transmission_Square_Capacitor_to_TTDLLogical
// G4LogicalVolume* Square_connection_from_TTLD_to_ResonatorLogical
// G4LogicalVolume* Square_connection_from_TTLD_to_ResonatorLogical_B
// G4LogicalVolume* Half_Cirlce_L_TTDLogical
// G4LogicalVolume* Half_Cirlce_R_TTDLogica
// G4LogicalVolume* Square_To_ResonatorLogical1
// G4LogicalVolume* Square_To_ResonatorLogical2
// G4LogicalVolume* Square_To_ResonatorLogical3
// G4LogicalVolume* Square_To_ResonatorLogical4
// G4LogicalVolume* Square_to_CapcitorH1Logical
// G4LogicalVolume* Square_to_CapcitorH2Logical
// G4LogicalVolume* Square_to_CapcitorH3Logical
// G4LogicalVolume* Square_to_CapcitorH4Logical
// G4LogicalVolume* Half_Cirlce_L_SquareCapacitorLogical
// G4LogicalVolume* Half_Cirlce_R_SquareCapacitorLogical
// G4LogicalVolume* Square_to_CapcitorV2Logical
// G4LogicalVolume* Square_to_CapcitorV1Logical
// G4LogicalVolume* Square_to_CapcitorVB1Logical
// G4LogicalVolume* Square_to_CapcitorVB2Logical
// G4LogicalVolume* Square_to_CapcitorVB3Logical
// G4LogicalVolume* Square_to_Qubit1Logical
// G4LogicalVolume* Square_to_Qubit2Logical
// G4LogicalVolume* Square_to_Qubit3Logical
// G4LogicalVolume* Square_to_Qubit4Logical
// G4LogicalVolume* SquareV_to_Qubit1Logical
// G4LogicalVolume* SquareV_to_Qubit2Logical
// G4LogicalVolume* SquareV_to_Qubit3Logical
// G4LogicalVolume* SquareV_to_Qubit4Logical
// G4LogicalVolume* Square_Capacitor_QubitLogical
// G4LogicalVolume* Square_Capacitor_Fill_1Logical
// G4LogicalVolume* Square_Capacitor_Fill_2Logical
// G4LogicalVolume* Square_Capacitor_Q_To_P_VS1Logical
// G4LogicalVolume* Square_Capacitor_Q_To_P_VS2Logical
// G4LogicalVolume* Square_Capacitor_Q_To_P_VS3Logical
// G4LogicalVolume* Square_Capacitor_Q_To_P_VS4Logical
// G4LogicalVolume* Half_Cirlce_C_To_PLogical
// G4LogicalVolume* Square_Capacitor_Q_To_P_HS1Logical
// G4LogicalVolume* Square_Capacitor_Q_To_P_HS2Logical
// G4LogicalVolume* Square_Capacitor_Q_To_P_HS3Logical
// G4LogicalVolume* Square_Capacitor_Q_To_P_HS4Logical
// G4LogicalVolume* Square_Capacitor_Q_To_P_V_P1Logical
// G4LogicalVolume* Copper_BoxLogical
// G4LogicalVolume* Sonten_Copper1Logical
// G4LogicalVolume* Sonten_Copper2Logical
//
//
//
// G4LogicalVolume* Resonator_SquareLogical_1_GP
// G4LogicalVolume* Half_Cirlce_LLogical_1_GP
// G4LogicalVolume* Half_Cirlce_RLogical_1_GP
// G4LogicalVolume* Resonator_SquareLogical_2_GP
// G4LogicalVolume* Half_Cirlce_RLogical_2_GP
// G4LogicalVolume* Half_Cirlce_LLogical_2_GP
// G4LogicalVolume* Half_Cirlce_LLogical__Transmission
// G4LogicalVolume* Quarter_Cirlce_RLogical__Transmission
// G4LogicalVolume* Quarter_Cirlce_LLogical__Transmission
// G4LogicalVolume* TransmissionLine_Resonator_SquareLogical_1
// G4LogicalVolume* Resonator_SquareLogical_3_GP
// G4LogicalVolume* Half_Cirlce_RLogical_3_GP
// G4LogicalVolume* Half_Cirlce_LLogical_3_GP
// G4LogicalVolume* Resonator_SquareLogical_4_GP
// G4LogicalVolume* Half_Cirlce_RLogical_4_GP
// G4LogicalVolume* Half_Cirlce_LLogical_4_GP
// G4LogicalVolume* Half_Cirlce_RLogical__Transmission_GP
// G4LogicalVolume* Half_Cirlce_LLogical__Transmission_GP
// G4LogicalVolume* Quarter_Cirlce_RLogical__Transmission_GP
// G4LogicalVolume* Quarter_Cirlce_LLogical__Transmission_GP
// G4LogicalVolume* TransmissionLine_Resonator_SquareLogical_1_GP
// G4LogicalVolume* Square_Capacitor_to_TTLDLogical_GP
// G4LogicalVolume* Half_Cirlce_R_Transmission_Square_Capacitor_to_TTDLLogical_GP
// G4LogicalVolume* Half_Cirlce_L_Transmission_Square_Capacitor_to_TTDLLogical_GP
// G4LogicalVolume* Square_connection_from_TTLD_to_ResonatorLogical_GP
// G4LogicalVolume* Square_connection_from_TTLD_to_ResonatorLogical_B_GP
// G4LogicalVolume* Half_Cirlce_L_TTDLogical_GP
// G4LogicalVolume* Half_Cirlce_R_TTDLogical_GP
// G4LogicalVolume* Big_Trian_PanelLogical
// G4LogicalVolume* Square_PanelLogical
// ////////////////////649 Code Line
// G4LogicalVolume* Square_To_ResonatorLogical1_GP
// G4LogicalVolume* Square_To_ResonatorLogical2_GP
// G4LogicalVolume* Square_To_ResonatorLogical3_GP
// G4LogicalVolume* Square_To_ResonatorLogical4_GP
// G4LogicalVolume* Square_to_CapcitorH1Logical_GP
// G4LogicalVolume* Square_to_CapcitorH2Logical_GP
// G4LogicalVolume* Square_to_CapcitorH3Logical_GP
// G4LogicalVolume* Square_to_CapcitorH4Logical_GP
// G4LogicalVolume* Half_Cirlce_L_SquareCapacitorLogical_GP
// G4LogicalVolume* Half_Cirlce_R_SquareCapacitorLogical_GP
// G4LogicalVolume* Square_to_CapcitorV1Logical_GP
// G4LogicalVolume* Square_to_CapcitorV2Logical_GP
// G4LogicalVolume* Square_to_CapcitorVB1Logical_GP
// G4LogicalVolume* Square_to_CapcitorVB2Logical_GP
// G4LogicalVolume* Square_to_CapcitorVB3Logical_GP
// G4LogicalVolume* Square_to_CapcitorVB4Logical
// G4LogicalVolume* Square_to_CapcitorVB4Logical_GP
// G4LogicalVolume* Square_to_Qubit1Logical_GP
// G4LogicalVolume* Square_to_Qubit2Logical_GP
// G4LogicalVolume* Square_to_Qubit3Logical_GP
// G4LogicalVolume* Square_to_Qubit4Logical_GP
// G4LogicalVolume* SquareV_to_Qubit1Logical_GP
// G4LogicalVolume* SquareV_to_Qubit2Logical_GP
// G4LogicalVolume* SquareV_to_Qubit3Logical_GP
// G4LogicalVolume* SquareV_to_Qubit4Logical_GP
// G4LogicalVolume* Square_Capacitor_QubitLogical_GP1
// G4LogicalVolume* Square_Capacitor_QubitLogical_GP2
// G4LogicalVolume* Square_Capacitor_QubitLogical_GP3
// G4LogicalVolume* Square_Capacitor_QubitLogical_GP4
// G4LogicalVolume* Square_Capacitor_Q_To_P_VS1Logical_GP
// G4LogicalVolume* Square_Capacitor_Q_To_P_VS2Logical_GP
// G4LogicalVolume* Square_Capacitor_Q_To_P_VS3Logical_GP
// G4LogicalVolume* Square_Capacitor_Q_To_P_VS4Logical_GP
// G4LogicalVolume* Half_Cirlce_C_To_PLogical_GP
// G4LogicalVolume* Square_Capacitor_Q_To_P_HS1Logical_GP
// G4LogicalVolume* Square_Capacitor_Q_To_P_HS2Logical_GP
// G4LogicalVolume* Square_Capacitor_Q_To_P_HS3Logical_GP
// G4LogicalVolume* Square_Capacitor_Q_To_P_HS4Logical_GP
// G4LogicalVolume* Square_Capacitor_Q_To_P_V_P1Logical_GP
// G4LogicalVolume* Square_PanelLogical_GP
// G4LogicalVolume* Big_Trian_PanelLogical_GP
// G4LogicalVolume* Square_PanelLogical_GP1



}

void DetectorConstruction::AttachPhononSensor(G4CMPSurfaceProperty * surfProp)
{
  //If no surface, don't do anything
  if(!surfProp) return;

  //Specify properties of the niobium sensors
  auto sensorProp = surfProp->GetPhononMaterialPropertiesTablePointer();
  //G4cout<<"#############This to test the parameter of the filmAbsorption ************************ \t"<<filmAbsorptiones<<G4endl;

  sensorProp->AddConstProperty("filmAbsorption",1.0);              //NOT WELL MOTIVATED - probably parametrize and put on slider?
  sensorProp->AddConstProperty("filmThickness",90.*CLHEP::nm);     //Accurate for our thin film.
  sensorProp->AddConstProperty("gapEnergy",1.6e-3*CLHEP::eV);       //Reasonably motivated. Actually, looks like Novotny and Meincke are quoting 2Delta, and this is delta. Nuss and Goossen mention that Nb has a delta value closer to this.
  sensorProp->AddConstProperty("lowQPLimit",3.);                   //NOT WELL MOTIVATED YET -- Dunno how to inform this...
  sensorProp->AddConstProperty("phononLifetime",4.17*CLHEP::ps);   //Kaplan paper says 242ps for Al, same table says 4.17ps for characteristic time for Nb.
  sensorProp->AddConstProperty("phononLifetimeSlope",0.29);        //Based on guessing from Kaplan paper, I think this is material-agnostic?
  sensorProp->AddConstProperty("vSound",3.480*CLHEP::km/CLHEP::s); //True for room temperature, probably good to 10%ish - should follow up
  sensorProp->AddConstProperty("subgapAbsorption",1.0);            //Assuming that since we're mostly sensitive to quasiparticle density, phonon "heat" here isn't something that we're sensitive to? Unsure how to select this.

  //  sensorProp->AddConstProperty("gapEnergy",3.0e-3*CLHEP::eV);      //Reasonably motivated. Novotny and Meincke, 1975 (2.8-3.14 meV)
  //  sensorProp->AddConstProperty("phononLifetime",242.*ps);      //Kaplan paper says 242ps for Al, same table says 4.17ps for characteristic time for Nb.

  surfProp->SetPhononElectrode(new G4CMPPhononElectrode);

}



  //
