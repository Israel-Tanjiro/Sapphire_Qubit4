/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file exoticphysics/phonon/include/DetectorConstructiontruction.hh
/// \brief Definition of the DetectorConstructiontruction class
//
// $Id: 4c06153e9ea08f2a90b22c53e5c39bde4b847c07 $
//
// 20221006  M. Kelsey -- Remove "IsField" flag, unnecessary with phonons.
//		Add material properties for aluminum phonon sensors

#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"

class G4Material;
class G4VPhysicalVolume;
class G4CMPSurfaceProperty;
class G4CMPElectrodeSensitivity;


class DetectorConstruction : public G4VUserDetectorConstruction {
public:
  DetectorConstruction();
  virtual ~DetectorConstruction();

public:
  virtual G4VPhysicalVolume* Construct();
  void AttachPhononSensor(G4CMPSurfaceProperty * surfProp);
private:
  void DefineMaterials();
  void SetupGeometry();
  //void AttachPhononSensor(G4CMPSurfaceProperty* surfProp);

private:
  G4Material* fLiquidHelium;
  G4Material* fSapphire;
  G4Material* fAluminum;
  G4Material* fNiobium;
  G4Material* copper_mat;
  G4Material* fTungsten;
  G4Material* fOxigen;
  G4Material* fAl2O3;
  G4VPhysicalVolume* fWorldPhys;
  // G4CMPSurfaceProperty* topSurfProp;
  // G4CMPSurfaceProperty* botSurfProp;
  // G4CMPSurfaceProperty* wallSurfProp;
  // G4CMPElectrodeSensitivity* electrodeSensitivity;

  G4bool fConstructed;		// Flag to not re-recreate surface properties
};

#endif
