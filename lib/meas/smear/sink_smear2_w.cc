// $Id: sink_smear2_w.cc,v 1.8 2004-01-06 04:56:47 edwards Exp $
/*! \file
 *  \brief Control routine for types of propagator smearing
 */

#include "chromabase.h"
#include "meas/smear/sink_smear2_w.h"
#include "meas/smear/gaus_smear.h"

using namespace QDP;


//! "Smear" the quark propagator at the sink by a covariant Gaussian
/*!
 * \ingroup smear
 *
 * This routine is specific to Wilson fermions!
 *
 * Arguments:
 *
 *  \param u                   gauge field ( Read )
 *  \param quark_propagator    quark propagator ( Modify )
 *  \param Wvf_kind            wave function kind: Gaussian or exponential
 *                             ( Read )
 *  \param wvf_param           wvf_param of "shell" wave function ( Read )
 *  \param wvfIntPar           number of iterations to approximate Gaussian
 *                             or terminate CG inversion for Wuppertal smearing
 *                             ( Read )
 *  \param j_decay             direction of decay ( Read ) 
 */

void sink_smear2(const multi1d<LatticeColorMatrix>& u,
                 LatticePropagator& quark_propagator, 
		 WvfKind Wvf_kind, const Real& wvf_param, int wvfIntPar, 
		 int j_decay)
{
  Real RsdCG;

  START_CODE("sink_smear2");

  switch (Wvf_kind) 
  {
  case WVF_KIND_GAUGE_INV_GAUSSIAN:
    gausSmear(u, quark_propagator, wvf_param, wvfIntPar, j_decay);
    break;

#if 0
  case WVF_KIND_WUPPERTAL:
    RsdCG = fuzz;
    wupp_smear(u, quark_propagator, wvf_param, wvfIntPar, j_decay, RsdCG);
    break;

  case WVF_KIND_JACOBI:
    RsdCG = fuzz;
    jacobi_smear(u, quark_propagator, wvf_param, wvfIntPar, j_decay);
    break;
#endif

  default :
    QDPIO::cerr << "sink_smear2: Unknown gauge invariant wave function" << endl;
    QDP_abort(1);
  }
    
  END_CODE("sink_smear2");
}

