// $Id: unprec_two_flavor_hasenbusch_monomial5d_w.cc,v 3.5 2006-12-28 15:39:00 bjoo Exp $
/*! @file
 * @brief Two-flavor collection of unpreconditioned 4D ferm monomials
 */

#include "chromabase.h"
#include "update/molecdyn/monomial/unprec_two_flavor_hasenbusch_monomial5d_w.h"
#include "update/molecdyn/monomial/monomial_factory.h"

#include "actions/ferm/fermacts/fermact_factory_w.h"
#include "actions/ferm/fermacts/fermacts_aggregate_w.h"

#include "update/molecdyn/predictor/chrono_predictor_factory.h"
#include "update/molecdyn/predictor/zero_guess_predictor.h"


namespace Chroma 
{ 
 
  namespace UnprecTwoFlavorHasenbuschWilsonTypeFermMonomial5DEnv 
  {
    namespace
    {
      //! Callback function for the factory
      Monomial< multi1d<LatticeColorMatrix>,
		multi1d<LatticeColorMatrix> >* createMonomial(XMLReader& xml, const string& path) 
      {
	return new UnprecTwoFlavorHasenbuschWilsonTypeFermMonomial5D(
	  TwoFlavorHasenbuschWilsonTypeFermMonomialParams(xml, path));
      }
 
      //! Local registration flag
      bool registered = false;
    }

    const std::string name("TWO_FLAVOR_UNPREC_HASENBUSCH_FERM_MONOMIAL5D");

    //! Register all the factories
    bool registerAll() 
    {
      bool success = true; 
      if (! registered)
      {
	success &= WilsonTypeFermActs5DEnv::registerAll();
	success &= TheMonomialFactory::Instance().registerObject(name, createMonomial);
	registered = true;
      }
      return success;
    }
  } //end namespace Unprec TwoFlavorHasenbuschWilsonFermMonomialEnv


  // Constructor
  UnprecTwoFlavorHasenbuschWilsonTypeFermMonomial5D::UnprecTwoFlavorHasenbuschWilsonTypeFermMonomial5D(
    const TwoFlavorHasenbuschWilsonTypeFermMonomialParams& param) 
  {
    START_CODE();

    inv_param = param.inv_param;

    std::istringstream is(param.fermact.xml);
    XMLReader fermact_reader(is);

    std::istringstream is_prec(param.fermact_prec.xml);
    XMLReader fermact_prec_reader(is_prec);

    // Check that the two 
    if( param.fermact_prec.id != param.fermact.id ) { 
      QDPIO::cerr << "For now both the numerator and the denominator fermacts mast be the same: You have asked for " 
		  << param.fermact.id 
		  << " in the denominator and " 
		  << param.fermact_prec.id << " in the numerator" << endl;
      QDP_abort(1);
    }
    

    QDPIO::cout << "UnprecTwoFlavorHasenbuschWilsonTypeFermMonomial: construct " << param.fermact.id << endl;
    WilsonTypeFermAct5D<T,P,Q>* tmp_act = 
      TheWilsonTypeFermAct5DFactory::Instance().createObject(param.fermact.id, fermact_reader, param.fermact.path);
   

    UnprecWilsonTypeFermAct5D<T,P,Q>* downcast = 
      dynamic_cast<UnprecWilsonTypeFermAct5D<T,P,Q>*>(tmp_act);

    // Check success of the downcast 
    if( downcast == 0x0 ) {
      QDPIO::cerr << "Unable to downcast FermAct to UnprecWilsonTypeFermAct in UnprecTwoFlavorWilsonTypeFermMonomial()" << endl;
      QDP_abort(1);
    }

    fermact = downcast;    

    QDPIO::cout << "UnprecTwoFlavorHasenbuschWilsonTypeFermMonomial: construct " << param.fermact_prec.id << endl;

    WilsonTypeFermAct5D<T,P,Q>* tmp_act_prec = 
      TheWilsonTypeFermAct5DFactory::Instance().createObject(param.fermact_prec.id, 
							     fermact_prec_reader, 
							     param.fermact_prec.path);

    UnprecWilsonTypeFermAct5D<T,P,Q>* downcast_prec = 
      dynamic_cast<UnprecWilsonTypeFermAct5D<T,P,Q>*>(tmp_act_prec);


    // Check success of the downcast 
    if( downcast_prec == 0x0 ) {
      QDPIO::cerr << "Unable to downcast FermAct to UnprecWilsonTypeFermAct in UnprecTwoFlavorWilsonTypeFermMonomial()" << endl;
      QDP_abort(1);
    }

    fermact_prec = downcast_prec;    

    if (fermact->size() != fermact_prec->size()) { 
      QDPIO::cerr << "Error: numerator action has to have the same length in the 5th dimension as the denominator action." << endl;
      QDPIO::cerr << "N5 in FermionAction " << fermact->size() << endl;
      QDPIO::cerr << "N5 in FermionActionPrec " << fermact_prec->size() << endl;
      QDP_abort(1);
    }
      

    //------------------------------------

    // Get Chronological predictor
    AbsChronologicalPredictor5D<LatticeFermion>* tmp = 0x0;
    if( param.predictor.xml == "" ) {
      // No predictor specified use zero guess
       tmp = new ZeroGuess5DChronoPredictor(fermact->size());
    }
    else 
    {
      try 
      { 
	std::istringstream chrono_is(param.predictor.xml);
	XMLReader chrono_xml(chrono_is);
	tmp = The5DChronologicalPredictorFactory::Instance().createObject(param.predictor.id, 
									  fermact->size(),
									  chrono_xml, 
									  param.predictor.path);
      }
      catch(const std::string& e ) { 
	QDPIO::cerr << "Caught Exception Reading XML: " << e << endl;
	QDP_abort(1);
      }


    }
     
    if( tmp == 0x0 ) { 
      QDPIO::cerr << "Failed to create the 4D ChronoPredictor" << endl;
      QDP_abort(1);
    }
    chrono_predictor = tmp;

    QDPIO::cout << "Initing PF field" << endl;
    getPhi().resize( fermact->size() );

    QDPIO::cout << "UnprecTwoFlavorWilsonTypeFermMonomial: finished " << param.fermact.id << endl;
    
    END_CODE();
  }

} //end namespace Chroma


