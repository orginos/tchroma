// -*- C++ -*-
//
/*! \file
 *  \ContinousBeta Monomial
 */

#ifndef __ContinuousBeta_monomial_h__
#define __ContinuousBeta_monomial_h__

#include "chromabase.h"

#include "update/molecdyn/field_state.h"
#include "update/molecdyn/monomial/abs_monomial.h"
#include "update/molecdyn/monomial/force_monitors.h"
#include "io/xmllog_io.h"

namespace Chroma 
{
  /*! @ingroup monomial */
  namespace ContinuousBetaMonomialEnv 
  {
    extern const string name;
    bool registerAll();
  }


  // Parameter structure
  /*! @ingroup monomial */
  struct ContinuousBetaMonomialParams 
  {
    // Base Constructor
	  ContinuousBetaMonomialParams();

    // Read monomial from some root path
	  ContinuousBetaMonomialParams(XMLReader& in, const std::string& path);
    string gauge_act;
  };


  //! Wrapper class for  CB monomials
  /*! @ingroup monomial
   *
   * Monomial is expected to be the same for these fermacts
   */
  class ContinuousBetaMonomial :
    public ExactMonomial< multi1d<LatticeColorMatrix>,
                          multi1d<LatticeColorMatrix> >    
  {
  public: 
    typedef multi1d<LatticeColorMatrix>  P;
    typedef multi1d<LatticeColorMatrix>  Q;

    //! Construct out of a parameter struct. Check against the desired GaugeAct name
    ContinuousBetaMonomial(const ContinuousBetaMonomialParams& param_);

    //! Copy Constructor
    ContinuousBetaMonomial(const ContinuousBetaMonomial& m) : gaugeact((m.gaugeact)) {}

    //! Create a suitable state and compute F
    void dsdq(P& F, const AbsFieldState<P,Q>& s) 
    {
      START_CODE();

      XMLWriter& xml_out = TheXMLLogWriter::Instance();
      push(xml_out, "ContinuousBetaMonomial");

      // Make a gauge connect state
      Handle< GaugeState<P,Q> > g_state(getGaugeAct().createState(s.getQ()));

      getGaugeAct().deriv(F, g_state);

      monitorForces(xml_out, "Forces", F);
      pop(xml_out);

      END_CODE();
    }


    //! Gauge action value
    Double S(const AbsFieldState<P,Q>& s)  
    {
      START_CODE();

      XMLWriter& xml_out = TheXMLLogWriter::Instance();
      push(xml_out, "ContinuousBetaMonomial");

      Handle< GaugeState<P,Q> > g_state(getGaugeAct().createState(s.getQ()));
      Double action = getGaugeAct().S(g_state);

      write(xml_out, "S", action);
      pop(xml_out);

      END_CODE();

      return action;
    }
	
	
    void refreshInternalFields(const AbsFieldState<P,Q>& s) 
    {
      //No internal fields to refresh => Nop
    	ReverceMCGaugeActEnv::setbeta();
    }

    void setInternalFields(const Monomial<P,Q>& m) 
    {
      // No internal fields to refresh => Nop
    }

    protected:
      const GaugeAction<P,Q>& getGaugeAct(void) const { 
	return *gaugeact;
      }

    private:
      // Hide empty constructor and =
      ContinuousBetaMonomial();
      void operator=(const ContinuousBetaMonomial&);

    private:
      // A handle for the gaugeact
      Handle< GaugeAction<P,Q> > gaugeact;
    };


}; //end namespace chroma

#endif
