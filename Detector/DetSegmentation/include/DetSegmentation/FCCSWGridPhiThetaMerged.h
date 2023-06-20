#ifndef DETSEGMENTATION_FCCSWGRIDPHITHETAMERGED_H
#define DETSEGMENTATION_FCCSWGRIDPHITHETAMERGED_H

// FCCSW
#include "DetSegmentation/GridTheta.h"

/** FCCSWGridPhiThetaMerged Detector/DetSegmentation/DetSegmentation/FCCSWGridPhiThetaMerged.h FCCSWGridPhiThetaMerged.h
 *
 *  Segmentation in theta and phi.
 *  Based on GridTheta, addition of azimuthal angle coordinate and merging of cells in eta/phi
 *
 */

namespace dd4hep {
namespace DDSegmentation {
class FCCSWGridPhiThetaMerged : public GridTheta {
public:
  /// default constructor using an arbitrary type
  FCCSWGridPhiThetaMerged(const std::string& aCellEncoding);
  /// Default constructor used by derived classes passing an existing decoder
  FCCSWGridPhiThetaMerged(const BitFieldCoder* decoder);

  /// destructor
  virtual ~FCCSWGridPhiThetaMerged() = default;

  /**  Determine the global position based on the cell ID.
   *   @warning This segmentation has no knowledge of radius, so radius = 1 is taken into calculations.
   *   @param[in] aCellId ID of a cell.
   *   return Position (radius = 1).
   */
  virtual Vector3D position(const CellID& aCellID) const;
  /**  Determine the cell ID based on the position.
   *   @param[in] aLocalPosition (not used).
   *   @param[in] aGlobalPosition position in the global coordinates.
   *   @param[in] aVolumeId ID of a volume.
   *   return Cell ID.
   */
  virtual CellID cellID(const Vector3D& aLocalPosition, const Vector3D& aGlobalPosition,
                        const VolumeID& aVolumeID) const;
  /**  Determine the azimuthal angle based on the cell ID.
   *   @param[in] aCellId ID of a cell.
   *   return Phi.
   */
  double phi(const CellID& aCellID) const;
  /**  Get the grid size in phi.
   *   return Grid size in phi.
   */
  inline double gridSizePhi() const { return 2 * M_PI / static_cast<double>(m_phiBins); }
  /**  Get the number of bins in azimuthal angle.
   *   return Number of bins in phi.
   */
  inline int phiBins() const { return m_phiBins; }
  /**  Get the coordinate offset in azimuthal angle.
   *   return The offset in phi.
   */
  inline double offsetPhi() const { return m_offsetPhi; }
  /**  Get the field name for azimuthal angle.
   *   return The field name for phi.
   */
  inline const std::string& fieldNamePhi() const { return m_phiID; }
  /**  Get the number of merged cells in phi for given layer
   *   @param[in] layer
   *   return The number of merged cells in phi
   */
  inline int mergedThetaCells(unsigned int layer) const { return m_mergedCellsTheta[layer]; }
  /**  Get the number of merged cells in theta for given layer
   *   @param[in] layer
   *   return The number of merged cells in theta
   */
  inline int mergedPhiCells(unsigned int layer) const { return m_mergedCellsPhi[layer]; }
  /**  Set the number of bins in azimuthal angle.
   *   @param[in] aNumberBins Number of bins in phi.
   */
  inline void setPhiBins(int bins) { m_phiBins = bins; }
  /**  Set the coordinate offset in azimuthal angle.
   *   @param[in] aOffset Offset in phi.
   */
  inline void setOffsetPhi(double offset) { m_offsetPhi = offset; }
  /**  Set the field name used for azimuthal angle.
   *   @param[in] aFieldName Field name for phi.
   */
  inline void setFieldNamePhi(const std::string& fieldName) { m_phiID = fieldName; }
  /**  Determine the polar angle based on the cell ID.
   *   @param[in] aCellId ID of a cell.
   *   return Theta.
   */
  double theta(const CellID& aCellID) const;

protected:
  /// determine the azimuthal angle phi based on the current cell ID
  double phi() const;
  /// the number of bins in phi
  int m_phiBins;
  /// the coordinate offset in phi
  double m_offsetPhi;
  /// vector of offsets in phi for each layer
  std::vector<double> m_offsetPhi_Layer;
  /// the field name used for phi
  std::string m_phiID;
  /// the field name used for layer
  std::string m_layerID;
  /// vector of number of cells to be merged along theta for each layer
  std::vector<int> m_mergedCellsTheta;
  /// vector of number of cells to be merged along phi for each layer
  std::vector<int> m_mergedCellsPhi;
};
}
}
#endif /* DETSEGMENTATION_FCCSWGRIDPHITHETAMERGED_H */
