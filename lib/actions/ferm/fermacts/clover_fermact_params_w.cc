// $Id: clover_fermact_params_w.cc,v 2.3 2006-03-01 18:58:23 bjoo Exp $
/*! \file
 *  \brief Clover fermion action parameters
 */

#include "chromabase.h"
#include "actions/ferm/fermacts/clover_fermact_params_w.h"

#include "io/param_io.h"

namespace Chroma
{

  //! Default constructor
  CloverFermActParams::CloverFermActParams()
  {
    Mass = Real(0);
    u0   = Real(1);
    clovCoeffR = clovCoeffT = Real(0);
    ext_fieldP = false;
    ext_field_strength = Real(0);
  }

  //! Read parameters
  CloverFermActParams::CloverFermActParams(XMLReader& xml, const string& path)
  {
    XMLReader paramtop(xml, path);

    // Read the stuff for the action
    if (paramtop.count("Mass") != 0) 
    {
      read(paramtop, "Mass", Mass);
      if (paramtop.count("Kappa") != 0) 
      {
	QDPIO::cerr << "Error: found both a Kappa and a Mass tag" << endl;
	QDP_abort(1);
      }
    }
    else if (paramtop.count("Kappa") != 0)
    {
      Real Kappa;
      read(paramtop, "Kappa", Kappa);
      Mass = kappaToMass(Kappa);    // Convert Kappa to Mass
      QDPIO::cout << "Kappa is " << Kappa << "Mass is " << Mass << endl << flush;
    }
    else
    {
      QDPIO::cerr << "Error: neither Mass or Kappa found" << endl;
      QDP_abort(1);
    }

    // Read optional u0
    if (paramtop.count("u0") != 0) 
      read(paramtop, "u0", u0);
    else {
      u0 = Real(1);
      QDPIO::cout << "u0 is " << u0 << endl << flush;
    }
    // Read optional anisoParam
    if (paramtop.count("AnisoParam") != 0) 
      read(paramtop, "AnisoParam", anisoParam);

    // If aniso, read all clover coeff
    if (anisoParam.anisoP)
    {
      read(paramtop, "clovCoeffR", clovCoeffR);
      read(paramtop, "clovCoeffT", clovCoeffT);
    }
    else
    {
      Real clovCoeff;
      read(paramtop, "clovCoeff", clovCoeff);
      clovCoeffR = clovCoeff;
      clovCoeffT = clovCoeff;
      QDPIO::cout << "clovCoeff=" << clovCoeff << endl;
      QDPIO::cout << "clovCoeffR=" << clovCoeffR << endl;
      QDPIO::cout << "clovCOeffT=" << clovCoeffT << endl;
    }

    if( paramtop.count("ext_field_strength") != 0 ) {
      ext_fieldP = true;
      read(paramtop, "ext_field_strength", ext_field_strength);
    }
  }

  //! Read parameters
  void read(XMLReader& xml, const string& path, CloverFermActParams& param)
  {
    CloverFermActParams tmp(xml, path);
    param = tmp;
  }


  //! Write parameters
  void write(XMLWriter& xml, const string& path, const CloverFermActParams& param)
  {
    push(xml, path);

    write(xml, "Mass", param.Mass);
    write(xml, "u0", param.u0);

    if (param.anisoParam.anisoP)
    {
      write(xml, "clovCoeffR", param.clovCoeffR);
      write(xml, "clovCoeffT", param.clovCoeffT);
    }
    else
    {
      write(xml, "clovCoeff", param.clovCoeffR);
    }
    
    if( param.ext_fieldP == true ) { 
      write(xml, "ext_field_strength", param.ext_field_strength);
    }

    pop(xml);
  }


}

