// $Id: prec_nef_linop_array_w.cc,v 1.11 2005-03-15 17:23:41 bjoo Exp $
/*! \file
 *  \brief  4D-style even-odd preconditioned NEF domain-wall linear operator
 */

#include "chromabase.h"
#include "actions/ferm/linop/prec_nef_linop_array_w.h"


namespace Chroma 
{ 

  //! Creation routine
  /*! \ingroup fermact
   *
   * \param u_            gauge field   (Read)
   * \param WilsonMass_   DWF height    (Read)
   * \param b5_           NEF parameter (Read)
   * \param c5_           NEF parameter (Read)
   * \param m_q_          quark mass    (Read)
   * \param N5_           extent of 5D  (Read)
   */
  void 
  EvenOddPrecNEFDWLinOpArray::create(const multi1d<LatticeColorMatrix>& u_, 
				     const Real& WilsonMass_, const Real &b5_, 
				     const Real &c5_, const Real& m_q_, int N5_)
  {
    START_CODE();
  
    WilsonMass = WilsonMass_;
    m_q = m_q_;
    b5  = b5_;
    c5  = c5_;
    N5  = N5_;
  
    D.create(u_);
  
  
    c5InvTwoKappa = 1.0 - c5*(Nd-WilsonMass) ;
    c5TwoKappa = 1.0 / c5InvTwoKappa ;
  
    b5InvTwoKappa = 1.0 + b5*(Nd-WilsonMass) ;
    b5TwoKappa = 1.0 / b5InvTwoKappa ;
  
    //InvTwoKappa = b5InvTwoKappa/c5InvTwoKappa ; 
    TwoKappa =  c5InvTwoKappa/b5InvTwoKappa ;
    Kappa = TwoKappa/2.0 ;
  
    invDfactor =1.0/(1.0  + m_q*pow(TwoKappa,N5)) ;


    END_CODE();
  }


  //! Apply the even-even (odd-odd) coupling piece of the domain-wall fermion operator
  /*!
   * \ingroup linop
   *
   * The operator acts on the entire lattice
   *
   * \param psi 	  Pseudofermion field     	       (Read)
   * \param isign   Flag ( PLUS | MINUS )   	       (Read)
   * \param cb      checkerboard ( 0 | 1 )               (Read)
   */
  void 
  EvenOddPrecNEFDWLinOpArray::applyDiag(multi1d<LatticeFermion>& chi, 
					const multi1d<LatticeFermion>& psi, 
					enum PlusMinus isign,
					const int cb) const
  {
    START_CODE();

    if( chi.size() != N5 ) chi.resize(N5);

    // Real c5Fact(0.5*c5InvTwoKappa) ; // The 0.5 is for the P+ and P-

    Real c5InvTwoKappamf = m_q*c5InvTwoKappa;
    switch ( isign ) {
    
    case PLUS:
    {
      for(int s(1);s<N5-1;s++) { // 1/2k psi[s] - P_- * psi[s+1] - P_+ * psi[s-1]


	//	chi[s][rb[cb]] = b5InvTwoKappa*psi[s] - 
	//  c5Fact*( psi[s+1] + psi[s-1] + GammaConst<Ns,Ns*Ns-1>()*(psi[s-1] - psi[s+1]) ) ;

	// Recoded using chiralProject and BLASology
	chi[s][rb[cb]] = b5InvTwoKappa*psi[s] - c5InvTwoKappa*chiralProjectPlus(psi[s-1]);
	chi[s][rb[cb]] -= c5InvTwoKappa*chiralProjectMinus(psi[s+1]);

      }



      //s=0 -- 1/2k psi[0] - P_- * psi[1] + mf* P_+ * psi[N5-1]
      //      chi[0][rb[cb]] = b5InvTwoKappa*psi[0] - 
      //	c5Fact*( psi[1]   - m_q*psi[N5m1] - GammaConst<Ns,Ns*Ns-1>()*(m_q*psi[N5m1] + psi[1]) ) ;

      // Recoded using chiralProject with BLAS-ology. c5Fact-s factor of 
      // 1/2 absorbed by projectors
      chi[0][rb[cb]] = b5InvTwoKappa*psi[0] - c5InvTwoKappa*chiralProjectMinus(psi[1]);
      chi[0][rb[cb]] += c5InvTwoKappamf*chiralProjectPlus(psi[N5-1]);


      //s=N5-1 -- 1/2k psi[N5-1] +mf* P_- * psi[0]  -  P_+ * psi[N5-2]
      // chi[N5m1][rb[cb]] = b5InvTwoKappa*psi[N5m1] - 
      //	c5Fact*( psi[N5m2] - m_q *psi[0] + GammaConst<Ns,Ns*Ns-1>()*(psi[N5m2] + m_q * psi[0]) );

      // Recoded with chiral projector and BLAS ology
      chi[N5-1][rb[cb]] = b5InvTwoKappa*psi[N5-1] - c5InvTwoKappa*chiralProjectPlus(psi[N5-2]);
      chi[N5-1][rb[cb]] += c5InvTwoKappamf*chiralProjectMinus(psi[0]);
    }
    break ;

    case MINUS:
    {    
      for(int s(1);s<N5-1;s++) { // 1/2k psi[s] - P_+ * psi[s+1] - P_- * psi[s-1]
	// chi[s][rb[cb]] = b5InvTwoKappa*psi[s] - 
	//  c5Fact*( psi[s+1] + psi[s-1] + GammaConst<Ns,Ns*Ns-1>()*(psi[s+1] - psi[s-1]) ) ;

	// Recoded using chiral projector and BLASology
	chi[s][rb[cb]] = b5InvTwoKappa*psi[s] - c5InvTwoKappa*chiralProjectPlus(psi[s+1]);
	chi[s][rb[cb]] -= c5InvTwoKappa*chiralProjectMinus(psi[s-1]);
      }

      //s=0 -- 1/2k psi[0] - P_+ * psi[1] + mf* P_- * psi[N5-1]
      //      chi[0][rb[cb]] = b5InvTwoKappa*psi[0] - 
      //	c5Fact*( psi[1]   - m_q*psi[N5-1] + GammaConst<Ns,Ns*Ns-1>()*( psi[1]+m_q*psi[N5-1]) ) ;
      chi[0][rb[cb]] = b5InvTwoKappa*psi[0] - c5InvTwoKappa*chiralProjectPlus(psi[1]);
      chi[0][rb[cb]] += c5InvTwoKappamf * chiralProjectMinus(psi[N5-1]);
      


      //s=N5-1 -- 1/2k psi[N5-1] + mf* P_+ * psi[0]  -  P_- * psi[N5-2]
      // chi[N5-1][rb[cb]] = b5InvTwoKappa*psi[N5-1] - 
      //	c5Fact*( psi[N5-2] - m_q *psi[0] - GammaConst<Ns,Ns*Ns-1>()*(psi[N5-2] + m_q * psi[0]) );

      // Recoded using Chiral Projector and BLAS ology
      chi[N5-1][rb[cb]] = b5InvTwoKappa*psi[N5-1] - c5InvTwoKappa*chiralProjectMinus(psi[N5-2]);
      chi[N5-1][rb[cb]] += c5InvTwoKappamf*chiralProjectPlus(psi[0]);
    }
    break ;
    }

    END_CODE();
  }


  //! Apply the inverse even-even (odd-odd) coupling piece of the domain-wall fermion operator
  /*!
   * \ingroup linop
   *
   * The operator acts on the entire lattice
   *
   * \param psi 	  Pseudofermion field     	       (Read)
   * \param isign   Flag ( PLUS | MINUS )   	       (Read)
   * \param cb      checkerboard ( 0 | 1 )               (Read)
   */
  void 
  EvenOddPrecNEFDWLinOpArray::applyDiagInv(multi1d<LatticeFermion>& chi, 
					   const multi1d<LatticeFermion>& psi, 
					   enum PlusMinus isign,
					   const int cb) const
  {
    START_CODE();
 
    if( chi.size() != N5 ) chi.resize(N5);
   
    // Copy and scale by TwoKappa (1/M0)
    for(int s(0);s<N5;s++)
      chi[s][rb[cb]] = b5TwoKappa * psi[s] ;


    switch ( isign ) {

    case PLUS:
    {
      
      // First apply the inverse of Lm 
      // Real fact(0.5*m_q*TwoKappa) ;
      
      // Recoded with chiral projectors
      // factor of 0.5 in fact absorbed into projectors
      Real fact = m_q*TwoKappa;

      for(int s(0);s<N5-1;s++){
	chi[N5-1][rb[cb]] -= fact * chiralProjectMinus(chi[s]);
	fact *= TwoKappa ;
      }
      
      //Now apply the inverse of L. Forward elimination 
      for(int s(1);s<N5;s++)
	chi[s][rb[cb]] += TwoKappa*chiralProjectPlus(chi[s-1]);
      
      //The inverse of D  now
      chi[N5-1][rb[cb]] *= invDfactor ;
      // That was easy....
      
      //The inverse of R. Back substitution...... Getting there! 
      for(int s(N5-2);s>-1;s--)
	chi[s][rb[cb]] += TwoKappa*chiralProjectMinus(chi[s+1]);

      
      //Finally the inverse of Rm 
      LatticeFermion tt;

      // Former factor of 0.5 in fact is absorbed into chiralProjector
      fact = m_q*TwoKappa;
      tt[rb[cb]] = fact*chiralProjectPlus(chi[N5-1]);
      for(int s(0);s<N5-1;s++){
	chi[s][rb[cb]] -= tt  ;
	tt[rb[cb]] *= TwoKappa ;
      }
    }
    break ;
    
    case MINUS:
    {
       
      // First apply the inverse of Lm 
      // Recoded using chiral projectors. The former factor of 0.5 in fact
      // is absorbed into the projectors
      Real fact = m_q*TwoKappa ;
      for(int s(0);s<N5-1;s++){
	chi[N5-1][rb[cb]] -= fact * chiralProjectPlus(chi[s]);
	fact *= TwoKappa ;
      }
      
      //Now apply the inverse of L. Forward elimination 
      for(int s(1);s<N5;s++)
	chi[s][rb[cb]] += TwoKappa*chiralProjectMinus(chi[s-1]);
      
      //The inverse of D  now
      chi[N5-1][rb[cb]] *= invDfactor ;
      // That was easy....
      
      //The inverse of R. Back substitution...... Getting there! 
      for(int s(N5-2);s>-1;s--)
	chi[s][rb[cb]] += TwoKappa*chiralProjectPlus(chi[s+1]);
      
      //Finally the inverse of Rm 
      LatticeFermion tt;
      fact = m_q*TwoKappa;

      tt[rb[cb]] = fact*chiralProjectMinus(chi[N5-1]);
      for(int s(0);s<N5-1;s++){
	chi[s][rb[cb]] -= tt  ;
	tt[rb[cb]] *= TwoKappa ;
      }
    }
    break ;
    }

    //Fixup the normalization. This step can probably be incorporated into
    // the above algerbra for more efficiency
    //for(int s(0);s<N5;s++)
    //  chi[s][rb[cb]] *= c5TwoKappa ;

    //Done! That was not that bad after all....
    //See, I told you so...
    END_CODE();
  }

  //! Apply the even-odd (odd-even) coupling piece of the NEF operator
  /*!
   * \ingroup linop
   *
   * The operator acts on the entire lattice
   *
   * \param psi 	  Pseudofermion field     	       (Read)
   * \param isign   Flag ( PLUS | MINUS )   	       (Read)
   * \param cb      checkerboard ( 0 | 1 )               (Read)
   */
  void 
  EvenOddPrecNEFDWLinOpArray::applyOffDiag(multi1d<LatticeFermion>& chi, 
					   const multi1d<LatticeFermion>& psi, 
					   enum PlusMinus isign,
					   const int cb) const 
  {
    START_CODE();

    Real fb5 = -Real(0.5)*b5 ;

    // Recoding with chiral projectors, a former factor of 0.5 is absorbed
    // into the projector
    Real fc5 = -Real(0.5)*c5 ;
    Real fc5mf = fc5*m_q;

    if( chi.size() != N5 ) chi.resize(N5);
  
    switch ( isign ) 
    {
    case PLUS:
    {
      multi1d<LatticeFermion> tmp(N5);
      int otherCB = (cb + 1)%2 ;
      
      for(int s = 1; s < N5-1; s++){
	//	tmp[s][rb[otherCB]] = fb5*psi[s] + 
	//   fc5*(psi[s+1] + psi[s-1] +
	//       GammaConst<Ns,Ns*Ns-1>()*(psi[s-1]-psi[s+1]));
	//
	// Recoded with chiral projectors and BLAS ology.
	tmp[s][rb[otherCB]] = fb5*psi[s] + fc5*chiralProjectPlus(psi[s-1]);
	tmp[s][rb[otherCB]]+= fc5*chiralProjectMinus(psi[s+1]);
      }
      
      
      // tmp[0][rb[otherCB]] = fb5*psi[0]  + 
      //	fc5*(psi[1] - m_q*psi[N5-1] 
      //     - GammaConst<Ns,Ns*Ns-1>()*(m_q*psi[N5-1] + psi[1]));
      //
      // Recoded with chiralProjec and BLAS ology
      tmp[0][rb[otherCB]] = fb5*psi[0] +fc5*chiralProjectMinus(psi[1]);
      tmp[0][rb[otherCB]] -= fc5mf*chiralProjectPlus(psi[N5-1]);

      // tmp[N5-1][rb[otherCB]] = fb5*psi[N5-1] + 
      //	fc5*( psi[N5-2] - m_q *psi[0] +
      //           GammaConst<Ns,Ns*Ns-1>()*(psi[N5-2] + m_q * psi[0]));
      tmp[N5-1][rb[otherCB]] = fb5*psi[N5-1] + fc5*chiralProjectPlus(psi[N5-2]);      tmp[N5-1][rb[otherCB]] -= fc5mf*chiralProjectMinus(psi[0]);

      // Replace this with a vector Dslash in time
      for(int s=0; s < N5; s++) { 
	D.apply(chi[s],tmp[s], isign, cb);
      }
      
    }
    break ;
    
    case MINUS:
    { 
      multi1d<LatticeFermion> tmp(N5) ;

      // Replace this with a vector Dslash in time
      for(int s(0);s<N5;s++){
	D.apply(tmp[s],psi[s],isign,cb);
      }

      for(int s(1);s<N5-1;s++){
	//	chi[s][rb[cb]] = fb5*tmp[s] + 
	//  fc5*(tmp[s+1] + tmp[s-1] -
	//       GammaConst<Ns,Ns*Ns-1>()*(tmp[s-1]-tmp[s+1]));
	//  
	//  Recoded using chiralProject and BLAS ology
	chi[s][rb[cb]] = fb5*tmp[s] + fc5*chiralProjectPlus(tmp[s+1]);
	chi[s][rb[cb]] += fc5*chiralProjectMinus(tmp[s-1]);
      }

      // chi[0][rb[cb]] = fb5*tmp[0]  + 
      //	fc5*(tmp[1] - m_q*tmp[N5-1] + 
      // GammaConst<Ns,Ns*Ns-1>()*(m_q*tmp[N5-1] + tmp[1]));
      //
      // Recoded using chiralProject and BLAS ology
      chi[0][rb[cb]] = fb5*tmp[0] + fc5*chiralProjectPlus(tmp[1]);
      chi[0][rb[cb]] -= fc5mf*chiralProjectMinus(tmp[N5-1]);

      
      // chi[N5-1][rb[cb]] = fb5*tmp[N5-1] + 
      //	fc5*( tmp[N5-2] - m_q *tmp[0] -
      //          GammaConst<Ns,Ns*Ns-1>()*(tmp[N5-2] + m_q * tmp[0]));
      chi[N5-1][rb[cb]] = fb5*tmp[N5-1] + fc5*chiralProjectMinus(tmp[N5-2]);
      chi[N5-1][rb[cb]] -= fc5mf*chiralProjectPlus(tmp[0]);

    }
    break ;
    }

    //Done! That was not that bad after all....
    //See, I told you so...
    END_CODE();
  }



  //! Apply the Dminus operator on a lattice fermion. See my notes ;-)
  void 
  EvenOddPrecNEFDWLinOpArray::Dminus(LatticeFermion& chi,
				     const LatticeFermion& psi,
				     enum PlusMinus isign,
				     int s5) const
  {
    LatticeFermion tt ;
    D.apply(tt,psi,isign,0);
    D.apply(tt,psi,isign,1);
    chi = c5InvTwoKappa*psi + (0.5*c5)*tt ;//It really is -(-0.5*c5)D 
  }


}; // End Namespace Chroma

