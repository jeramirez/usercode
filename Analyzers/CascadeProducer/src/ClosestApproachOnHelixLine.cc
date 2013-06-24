// -*- C++ -*-
//
// Package:    CascadeProducer
// Class:      ClosestApproachOnHelixLine
// 
/**\class CascadeFitter ClosestApproachOnHelixLine.cc Analyzers/CascadeProducer/src/ClosestApproachOnHelixLine.cc

 Description: compute the points of closest  approach of an helix and a straight line.

 Implementation:
     Adapted from ClosestApproachInRPhi
*/
//
// Original Author:  Juan Eduardo Ramirez Vargas
//         Created:  Tue Jun 23 17:52:47 CDT 2009
// $Id: ClosestApproachOnHelixLine.cc,v 1.1 2009/07/13 23:02:48 jramirez Exp $
//
//
//
#include "Analyzers/CascadeProducer/interface/ClosestApproachOnHelixLine.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "FWCore/Utilities/interface/Exception.h"


std::pair<GlobalPoint, GlobalPoint> ClosestApproachOnHelixLine::points() const
{
  if (!status_)
    throw cms::Exception("Analyzers/CascadeProducer","ClosestApproachOnHelixLine::could not compute track crossing. Check status before calling this method!");
  return  std::pair<GlobalPoint, GlobalPoint> (posA, posB);
}

GlobalPoint 
ClosestApproachOnHelixLine::crossingPoint() const
{
  if (!status_)
    throw cms::Exception("Analyzers/CascadeProducer","ClosestApproachOnHelixLine::could not compute track crossing. Check status before calling this method!");
  return GlobalPoint((posA.x() + posB.x())/2., 
                     (posA.y() + posB.y())/2., 
                     (posA.z() + posB.z())/2.);
}
 
double ClosestApproachOnHelixLine::distance() const
{
  if (!status_)
    throw cms::Exception("Analyzers/CascadeProducer","ClosestApproachOnHelixLine::could not compute track crossing. Check status before calling this method!");
  return (posB - posA).mag();
}

bool ClosestApproachOnHelixLine::calculate(const TrackCharge & chargeCircle, 
                                            const GlobalVector & momentumCircle, 
                                            const GlobalPoint & positionCircle, 
                                            const GlobalVector & momentumLine, 
                                            const GlobalPoint & positionLine,
                                            const MagneticField& magField) 
{
   // centre and radii of track circle
  double xca, yca, ra;
  circleParameters(chargeCircle, momentumCircle, positionCircle, xca, yca, ra, magField);

   // points of closest approach
  double xg1,yg1,zg1,xg2,yg2,zg2;
  int flag = transverseCoord(xca, yca, ra, momentumLine, positionLine, xg1,yg1,zg1,xg2,yg2,zg2);
  if (flag == 0) {
    status_ = false;
    return false;
  }

  double zgtrack1 = zCoord(momentumCircle,positionCircle,ra,xca,yca,xg1,yg1);
  double zgtrack2 = zCoord(momentumCircle,positionCircle,ra,xca,yca,xg2,yg2);

  double xcir,ycir,zcir,xlin,ylin,zlin;
  if (flag ==1){
     //circle don't cross line 
     xcir=xg1; ycir=yg1; zcir=zgtrack1;
     xlin=xg2; ylin=yg2; zlin=zg2; 
  }else{
     //circle cross line, must choose the closest in z
     if (std::abs(zg1-zgtrack1) < std::abs(zg2-zgtrack2)){
        xcir=xg1; ycir=yg1; zcir=zgtrack1;
        xlin=xg1; ylin=yg1; zlin=zg1;
     }else{
        xcir=xg2; ycir=yg2; zcir=zgtrack2;
        xlin=xg2; ylin=yg2; zlin=zg2;
     }//break tie using z
  }//deal with solutions

  posA = GlobalPoint(xcir,ycir,zcir);
  posB = GlobalPoint(xlin,ylin,zlin);
  status_ = true;
  return true;

}
// Handful code to simplify copied from ClosestApproachInRPhi 
void ClosestApproachOnHelixLine::circleParameters(const TrackCharge& charge, 
					const GlobalVector& momentum, 
					const GlobalPoint& position, 
					double& xc, double& yc, double& r,
					const MagneticField& magField) 
const
{

  // compute radius of circle
  /** temporary code, to be replaced by call to curvature() when bug 
   *  is fixed. 
   */
//   double bz = MagneticField::inInverseGeV(position).z();
  double bz = magField.inTesla(position).z() * 2.99792458e-3;

  // signed_r directed towards circle center, along F_Lorentz = q*v X B
  double signed_r = charge*momentum.transverse() / bz;
  r = std::abs(signed_r);
  /** end of temporary code
   */

  // compute centre of circle
  double phi = momentum.phi();
  xc =  signed_r*sin(phi) + position.x();
  yc = -signed_r*cos(phi) + position.y();

}

int 
ClosestApproachOnHelixLine::transverseCoord(
                          double xca, double yca, double ra, 
                          const GlobalVector& momentum,
			  const GlobalPoint& position, 
			  double & xg1, double & yg1, double & zg1, 
                          double & xg2, double & yg2, double & zg2
			  )
{
   int flag = 0;
   double x1, y1, z1, x2, y2, z2;
 
   //in the tranversal plane (2D) compute 
   //distance components (xv,yv)-(xc,yc) in the XY plane
   //dtransversal is perpendicular to pt direction in the 2D plane
   //dlongitudinal is along the pt direction in the 2D plane
   dtransversal = std::abs(
	                          (position.x()-xca)*momentum.y()-
				  (position.y()-yca)*momentum.x()
	                         )/momentum.transverse();

   dlongitudinal =((position.x()-xca)*momentum.x() +
	                  (position.y()-yca)*momentum.y()
			 )/momentum.transverse();		
   radii = ra;			 
   double rlongitudinal2 = ra*ra -  dtransversal*dtransversal;
   double rl=rlongitudinal2>0?sqrt(rlongitudinal2):0;
   x1 = position.x()-(dlongitudinal-rl)*momentum.x()/momentum.transverse();
   y1 = position.y()-(dlongitudinal-rl)*momentum.y()/momentum.transverse();
   z1 = position.z()-(dlongitudinal-rl)*momentum.z()/momentum.transverse();

   x2 = position.x()-(dlongitudinal+rl)*momentum.x()/momentum.transverse();
   y2 = position.y()-(dlongitudinal+rl)*momentum.y()/momentum.transverse();
   z2 = position.z()-(dlongitudinal+rl)*momentum.z()/momentum.transverse();
 
   if ( rl == 0 ){
     // circle is external to line and don't cross
     // one solution
     flag = 1;
     xg1 = (1-ra/dtransversal)*xca + ra/dtransversal*x1;
     yg1 = (1-ra/dtransversal)*yca + ra/dtransversal*y1;
     zg1 = z1;
     xg2 = x2; yg2=y2; zg2=z2; //save to compute dca to line
   } 
   else {  
     // circle cross line two solutions
     flag = 2;
     xg1 = x1; yg1=y1; zg1=z1;
     xg2 = x2; yg2=y2; zg2=z2; 
   }
 
  return flag;
}

double ClosestApproachOnHelixLine::zCoord(
                             const GlobalVector& mom, 
                             const GlobalPoint & pos, 
			     double r, double xc, double yc, 
			     double xg, double yg) const
{

  // starting point
  double x = pos.x(); double y = pos.y(); double z = pos.z();

  double px = mom.x(); double py = mom.y(); double pz = mom.z();

  // rotation angle phi from starting point to crossing point (absolute value)
  // -- compute sin(phi/2) if phi smaller than pi/4, 
  // -- cos(phi) if phi larger than pi/4
  double phi = 0.;
  double sinHalfPhi = sqrt((x-xg)*(x-xg) + (y-yg)*(y-yg))/(2*r);
  if (sinHalfPhi < 0.383) { // sin(pi/8)
    phi = 2*asin(sinHalfPhi);
  }
  else {
    double cosPhi = ((x-xc)*(xg-xc) + (y-yc)*(yg-yc))/(r*r);
    phi = std::abs(acos(cosPhi));
  }
  // -- sign of phi
  double signPhi = ((x - xc)*(yg - yc) - (xg - xc)*(y - yc) > 0) ? 1. : -1.;

  // sign of track angular momentum
  // if rotation is along angular momentum, delta z is along pz
  double signOmega = ((x - xc)*py - (y - yc)*px > 0) ? 1. : -1.;

  // delta z
  // -- |dz| = |cos(theta) * path along helix|
  //         = |cos(theta) * arc length along circle / sin(theta)|
  double dz = signPhi*signOmega*(pz/mom.transverse())*phi*r;

  return z + dz;
 }


