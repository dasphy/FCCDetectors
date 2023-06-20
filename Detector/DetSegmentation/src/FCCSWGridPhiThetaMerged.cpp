#include "DetSegmentation/FCCSWGridPhiThetaMerged.h"

#include <iostream>

namespace dd4hep {
namespace DDSegmentation {

/// default constructor using an encoding string
FCCSWGridPhiThetaMerged::FCCSWGridPhiThetaMerged(const std::string& cellEncoding) : GridTheta(cellEncoding) {
  // define type and description
  _type = "FCCSWGridPhiThetaMerged";
  _description = "Phi-theta segmentation in the global coordinates";

  // register all necessary parameters (additional to those registered in GridTheta)
  registerParameter("phi_bins", "Number of bins phi", m_phiBins, 1);
  registerParameter("offset_phi", "Angular offset in phi", m_offsetPhi, 0., SegmentationParameter::AngleUnit, true);
  registerParameter("offset_phi_layer", "Angular offset in phi for each layer", m_offsetPhi_Layer, std::vector<double>(), SegmentationParameter::AngleUnit, true);
  registerIdentifier("identifier_phi", "Cell ID identifier for phi", m_phiID, "phi");
  registerIdentifier("identifier_layer", "Cell ID identifier for layer", m_layerID, "layer");
  registerParameter("mergedCells_Theta", "Numbers of merged cells in theta per layer", m_mergedCellsTheta, std::vector<int>());
  registerParameter("mergedCells_Phi", "Numbers of merged cells in phi per layer", m_mergedCellsPhi, std::vector<int>());
}

FCCSWGridPhiThetaMerged::FCCSWGridPhiThetaMerged(const BitFieldCoder* decoder) : GridTheta(decoder) {
  // define type and description
  _type = "FCCSWGridPhiThetaMerged";
  _description = "Phi-theta segmentation in the global coordinates";

  // register all necessary parameters (additional to those registered in GridTheta)
  registerParameter("phi_bins", "Number of bins phi", m_phiBins, 1);
  registerParameter("offset_phi", "Angular offset in phi", m_offsetPhi, 0., SegmentationParameter::AngleUnit, true);
  registerParameter("offset_phi_layer", "Angular offset in phi for each layer", m_offsetPhi_Layer, std::vector<double>(), SegmentationParameter::AngleUnit, true);
  registerIdentifier("identifier_phi", "Cell ID identifier for phi", m_phiID, "phi");
  registerIdentifier("identifier_layer", "Cell ID identifier for layer", m_layerID, "layer");
  registerParameter("mergedCells_Theta", "Numbers of merged cells in theta per layer", m_mergedCellsTheta, std::vector<int>());
  registerParameter("mergedCells_Phi", "Numbers of merged cells in phi per layer", m_mergedCellsPhi, std::vector<int>());
}

/// determine the local position based on the cell ID
/// no need to change anything
Vector3D FCCSWGridPhiThetaMerged::position(const CellID& cID) const {
  return positionFromRThetaPhi(1.0, theta(cID), phi(cID));
}

/// determine the cell ID based on the global position
CellID FCCSWGridPhiThetaMerged::cellID(const Vector3D& /* localPosition */, const Vector3D& globalPosition,
                          const VolumeID& vID) const {
  CellID cID = vID;

  // debug
  // std::cout << "x: " << globalPosition.X << " y: " << globalPosition.Y << " z: " << globalPosition.Z << std::endl;

  int layer = _decoder->get(cID, m_layerID);
  // debug
  // std::cout << "Layer " << layer << " from volume ID: " << vID << "\n" << std::endl;

  // double lTheta = thetaFromXYZ(globalPosition);
  // lTheta is always pi/2 because z=0 after CreateCaloCells (could this be fixed easily or not)?
  // temporary hack below - but to make it work I had to modify code in the caller such that vID contains
  // already theta in the bitfield corresponding to the new readout encoding..
  int thetaBin = _decoder->get(cID, m_thetaID);
  // std::cout << "Theta bin: " << thetaBin << std::endl;
  // std::cout << "lTheta: " << lTheta << " m_gridSizeTheta: " << m_gridSizeTheta << " m_offsetTheta: " << m_offsetTheta << std::endl;
  // adjust theta bin if cells are merged along theta in this layer
  if (m_mergedCellsTheta[layer]>1) {

    // debug
    // std::cout << "Number of cells to merge along theta: " << m_mergedCellsTheta[layer] << std::endl;

    if ((thetaBin % m_mergedCellsTheta[layer])>0)
      thetaBin -= (thetaBin % m_mergedCellsTheta[layer]);

    // debug
    // std::cout << "New theta bin: " << thetaBin << std::endl;
  }
  _decoder->set(cID, m_thetaID, thetaBin);


  double lPhi = phiFromXYZ(globalPosition);
  int phiBin = positionToBin(lPhi, 2 * M_PI / (double)m_phiBins, m_offsetPhi);
  // debug
  // std::cout << "Phi bin: " << phiBin << std::endl;
  // std::cout << "lPhi: " << lPhi << " m_gridSizePhi: " << 2 * M_PI / (double)m_phiBins << " m_offsetPhi: " << m_offsetPhi << std::endl;
  // adjust phi bin if cells are merged along phi in this layer
  if (m_mergedCellsPhi[layer]>1) {

    // debug
    // std::cout << "Number of cells to merge along phi: " << m_mergedCellsPhi[layer] << std::endl;

    if ((phiBin % m_mergedCellsPhi[layer])>0)
      phiBin -= (phiBin % m_mergedCellsPhi[layer]);
  }

  // debug
  // std::cout << "New phi bin: " << phiBin << std::endl;

  _decoder->set(cID, m_phiID, phiBin);

  return cID;
}

/// determine the azimuthal angle phi based on the cell ID
double FCCSWGridPhiThetaMerged::phi(const CellID& cID) const {

  // retrieve layer
  int layer = _decoder->get(cID, m_layerID);

  // retrieve phi
  CellID phiValue = _decoder->get(cID, m_phiID);
  // double _phi = binToPosition(phiValue, 2. * M_PI / (double)m_phiBins, m_offsetPhi_Layer[layer]);
  double _phi = binToPosition(phiValue, 2. * M_PI / (double)m_phiBins, m_offsetPhi);

  // debug
  // std::cout << "Layer: " << layer << " phi bin: " << phiValue << " initial phi: " << _phi << std::endl;

  // adjust return value if cells are merged along phi in this layer
  // shift by (N-1)*half phi grid size
  if (m_mergedCellsPhi[layer]>1) {
    // debug
    // std::cout << "Number of cells to merge along phi: " << m_mergedCellsPhi[layer] << std::endl;
    _phi += (m_mergedCellsPhi[layer]-1)* M_PI / (double)m_phiBins ;
  }

  // debug
  // std::cout << "Final phi: " << _phi << std::endl;

  return _phi;
}

double FCCSWGridPhiThetaMerged::theta(const CellID& cID) const {

  // retrieve layer
  int layer = _decoder->get(cID, m_layerID);

  // retrieve theta
  CellID thetaValue = _decoder->get(cID, m_thetaID);
  double _theta = binToPosition(thetaValue, m_gridSizeTheta, m_offsetTheta);

  // debug
  // std::cout << "Layer: " << layer << " theta bin: " << thetaValue << " initial theta: " << _theta << std::endl;

  // adjust return value if cells are merged along theta in this layer
  // shift by (N-1)*half theta grid size
  if (m_mergedCellsTheta[layer]>1) {
    // debug
    // std::cout << "Number of cells to merge along theta: " << m_mergedCellsTheta[layer] << std::endl;
    _theta += (m_mergedCellsTheta[layer]-1) * m_gridSizeTheta / 2.0 ;
  }

  // debug
  // std::cout << "Final theta: " << _theta << std::endl;

  return _theta;
}

}
}
