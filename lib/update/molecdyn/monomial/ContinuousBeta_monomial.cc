// -*- C++ -*-
// 
/*! \file
 */

#include "chromabase.h"

#include "update/molecdyn/monomial/ContinuousBeta_monomial.h"
#include "update/molecdyn/monomial/monomial_factory.h"
#include "actions/gauge/gaugeacts/gaugeacts_aggregate.h"
#include "actions/gauge/gaugeacts/gaugeact_factory.h"

namespace Chroma 
{ 
  namespace ContinuousBetaMonomialEnv 
  {
    namespace
    {
      //! Callback function for the factory
      Monomial< multi1d<LatticeColorMatrix>,
		multi1d<LatticeColorMatrix> >*
      createMonomial(XMLReader& xml, const string& path) 
      {
	QDPIO::cout << "Create monomial: " << name << endl;

	return new ContinuousBetaMonomial(ContinuousBetaMonomialParams(xml, path));
      }

      //! Local registration flag
      bool registered = false;
    }

    const std::string name("CONTINUOUS_BETA_GAUGE_MONOMIAL");
    
    //! Register all the factories
    bool registerAll() 
    {
      bool success = true; 
      if (! registered)
      {
	success &= GaugeActsEnv::registerAll();
	success &= TheMonomialFactory::Instance().registerObject(name, createMonomial);
	registered = true;
      }
      return success;
    }
  } //end namespace ContinuousBetaMonomialEnv



  // Read the parameters
  ContinuousBetaMonomialParams::ContinuousBetaMonomialParams(XMLReader& xml_in, const string& path)
  {
    // Get the top of the parameter XML tree
    XMLReader paramtop(xml_in, path);
    
    try {
      // Read the inverter Parameters
      XMLReader xml_tmp(paramtop, "./GaugeAction");
      std::ostringstream os;
      xml_tmp.print(os);
      gauge_act = os.str();
   
    }
    catch(const string& s) {
      QDPIO::cerr << "Caught Exception while reading parameters: " << s <<endl;
      QDP_abort(1);
    }

    QDPIO::cout << "ContinuousBetaMonomialParams: read \n" << gauge_act << endl;
  }

  //! Read Parameters
  void read(XMLReader& xml, const std::string& path,
	    ContinuousBetaMonomialParams& params)
  {
    GaugeMonomialParams tmp(xml, path);
    params = tmp;
  }

  //! Write Parameters
  void write(XMLWriter& xml, const std::string& path,
	     const ContinuousBetaMonomialParams& params)
  {
    // Not implemented
    QDPIO::cerr << ContinuousBetaMonomialEnv::name << ": write not implemented" << endl;
    QDP_abort(1);
  }


  // Constructor
  ContinuousBetaMonomial::ContinuousBetaMonomial(const ContinuousBetaMonomialParams& param_)
  {
    std::istringstream is(param_.gauge_act);
    XMLReader gaugeact_reader(is);

    // Get the name of the gauge act
    std::string gaugeact_string;
    try { 
      read(gaugeact_reader, "/GaugeAction/Name", gaugeact_string);
    }
    catch( const std::string& e) { 
      QDPIO::cerr << "Error grepping the gaugeact name: " << e<<  endl;
      QDP_abort(1);
    }

    // Throw an exception if not found
    gaugeact = TheGaugeActFactory::Instance().createObject(gaugeact_string, 
							   gaugeact_reader, 
							   "/GaugeAction");
  }
  
  void refreshInternalFields(const AbsFieldState<P,Q>& s) 
      {
	  	  ReverceMCGaugeActEnv::setbeta();
      }

} //end namespace Chroma


