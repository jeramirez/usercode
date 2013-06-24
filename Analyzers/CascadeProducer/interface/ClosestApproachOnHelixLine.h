#ifndef _ClosestApproachOnHelixLine_H_
#define _ClosestApproachOnHelixLine_H_

#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/GeometryVector/interface/GlobalVector.h"
#include "DataFormats/TrajectoryState/interface/TrackCharge.h"

class MagneticField;

class ClosestApproachOnHelixLine {
public:
   ClosestApproachOnHelixLine() {status_      = false;
                                 dtransversal = -1;
                                 dlongitudinal= -999;
                                }
  ~ClosestApproachOnHelixLine() {}

   
   bool status() const {return status_;}

  /**
   * Returns the two PCA on the trajectories.
   */
  std::pair<GlobalPoint, GlobalPoint> points() const;

   /** arithmetic mean of the two points of closest approach */
  GlobalPoint crossingPoint() const;
   /** distance between the two points of closest approach in 3D */
  double distance() const;

  /**
   *  Clone method
   */
  ClosestApproachOnHelixLine* clone() const {
    return new ClosestApproachOnHelixLine(* this);
  }

  bool calculate(const TrackCharge & chargeCircle, 
		 const GlobalVector & momentumCircle, 
		 const GlobalPoint & positionCircle, 
		 const GlobalVector & momentumLine, 
		 const GlobalPoint & positionLine,
		 const MagneticField& magField);


  // Computes center coordinates and unsigned radius of circle;
  void circleParameters(const TrackCharge& charge, 
		        const GlobalVector& momentum, 
		        const GlobalPoint& position, 
		        double& xc, double& yc, double& r,
		        const MagneticField& magField) const;

  // Computes crossing points of circle with centre (cx, cy)
  // and unsigned radii r with line passing by a global point
  // and pointing in a given direction
  // Two cases: - circle do not cross the line
  //                take 1 solution
  //              return value = 1;
  //            - circle  have two intersection points;
  //                take 2 solutions
  //              return value = 2;
  // if the calculation fails (e.g. line towards circle), return value = 0;
  int transverseCoord(double xca, double yca, double ra,
                      const GlobalVector& momentum,
                      const GlobalPoint& position,
                      double & xg1, double & yg1, double & zg1,
                      double & xg2, double & yg2, double & zg2
                     );

  // Computes z-coordinate on helix at given transverse coordinates
  double zCoord(const GlobalVector& mom,
                const GlobalPoint& pos,
                double r, double xc, double yc,
                double xg, double yg) const;

  double d_t() const {return dtransversal;}
  double d_l() const {return dlongitudinal;}
  double radius() const {return radii;}
private:
  bool status_;
  GlobalPoint posA, posB;
  double dtransversal,dlongitudinal,radii;

};//end class

#endif

