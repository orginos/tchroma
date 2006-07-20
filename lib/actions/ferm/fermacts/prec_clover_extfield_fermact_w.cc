// $Id: prec_clover_extfield_fermact_w.cc,v 3.1 2006-07-20 20:06:52 edwards Exp $
/*! \file
 *  \brief Even-odd preconditioned Clover fermion action with an external field
 */

#include "chromabase.h"
#include "actions/ferm/linop/prec_clover_extfield_linop_w.h"
#include "actions/ferm/fermacts/prec_clover_extfield_fermact_w.h"

#include "actions/ferm/fermacts/fermact_factory_w.h"
//#include "actions/ferm/fermacts/ferm_createstate_reader_w.h"

#include "actions/ferm/fermacts/extfield_fermstate_w.h"
#include "actions/ferm/fermacts/extfield_factory_w.h"
#include "actions/ferm/fermacts/extfield_aggregate_w.h"
#include "actions/ferm/fermbcs/fermbcs_reader_w.h"


namespace Chroma
{

  //! Hooks to register the class with the fermact factory
  namespace EvenOddPrecCloverExtFieldFermActEnv
  {
    // Helper function for the FermAction readers
    Handle< CreateFermState<LatticeFermion,
			    multi1d<LatticeColorMatrix>, 
			    multi1d<LatticeColorMatrix> > > reader(XMLReader& xml_in, 
								   const std::string& path)
    {
      XMLReader top(xml_in, path);

      {
	ostringstream os;
	top.printCurrentContext(os);
	cout << __func__ << ": top = XX" << os.str() << "XX" << endl;
      }

      std::string fermstate = "FermState";
      if (top.count(fermstate) == 0)
      {
	QDPIO::cerr << EvenOddPrecCloverExtFieldFermActEnv::name 
		    << ": did not find appropriate FermState" << endl;
	QDP_abort(1);
      }

      XMLReader paramtop(top, fermstate);

      // Read the array of function names
      multi1d< Handle< ExternalField > > ext_field(Nd);

      try
      {
	for(int mu=0; mu < ext_field.size(); ++mu)
	{
	// Create the query for the element 
	  std::ostringstream element_xpath;
	  element_xpath << "ExternalField/elem[" << (mu+1) << "]";

	  string name;
	  read(paramtop, element_xpath.str() + "/Name", name);

	  QDPIO::cout << "External field type = " << name << endl;
	  ext_field[mu] = TheExternalFieldFactory::Instance().createObject(name,
									   paramtop,
									   element_xpath.str());
	}
      }
      catch(const std::string& e) 
      {
	QDPIO::cerr << __func__ << ": caught exception reading XML: " << e << endl;
	QDP_abort(1);
      }
      
      return new CreateExtFieldFermState<LatticeFermion,
	multi1d<LatticeColorMatrix>, 
	multi1d<LatticeColorMatrix> >(WilsonTypeFermBCEnv::reader(paramtop, "FermionBC"),
				      ext_field);
    }


    //! Callback function
    WilsonTypeFermAct<LatticeFermion,
		      multi1d<LatticeColorMatrix>,
		      multi1d<LatticeColorMatrix> >* createFermAct4D(XMLReader& xml_in,
								     const std::string& path)
    {
      return new EvenOddPrecCloverExtFieldFermAct(reader(xml_in, path), 
						  CloverFermActParams(xml_in, path));
    }

    //! Callback function
    /*! Differs in return type */
    FermionAction<LatticeFermion,
		  multi1d<LatticeColorMatrix>,
		  multi1d<LatticeColorMatrix> >* createFermAct(XMLReader& xml_in,
							       const std::string& path)
    {
      return createFermAct4D(xml_in, path);
    }

    //! Name to be used
    const std::string name = "CLOVER_EXTERNAL_FIELD";

    //! Register all the factories
    bool registerAll()
    {
      bool foo = true;
      foo &= ExternalFieldEnv::registered;
      foo &= Chroma::TheFermionActionFactory::Instance().registerObject(name, createFermAct);
      foo &= Chroma::TheWilsonTypeFermActFactory::Instance().registerObject(name, createFermAct4D);
    }

    //! Register the fermact
    const bool registered = registerAll();
  }

  //! Produce a linear operator for this action
  /*!
   * The operator acts on the odd subset
   *
   * \param state	    gauge field     	       (Read)
   */
  EvenOddPrecLogDetLinearOperator<LatticeFermion,
				  multi1d<LatticeColorMatrix>,
				  multi1d<LatticeColorMatrix> >* 
  EvenOddPrecCloverExtFieldFermAct::linOp(Handle< FermState<T,P,Q> > state) const
  {
    return new EvenOddPrecCloverExtFieldLinOp(state, param);
  }

}

