// -*- C++ -*-
/*! \file
 *  \brief Null based map object, factory registration
 */

#include <string>
#include "chromabase.h"
#include "qdp_map_obj_null.h"
#include "util/ferm/map_obj/map_obj_factory_w.h"
#include "util/ferm/key_prop_colorvec.h"

namespace Chroma 
{ 
  
  namespace MapObjectNullEnv 
  {
    namespace
    {
      //! Callback function
      QDP::MapObject<int,EVPair<LatticeColorVector> >* createMapObjIntKeyCV(XMLReader& xml_in,
									    const std::string& path) 
      {
	// Doesn't need parameters...
	return new QDP::MapObjectNull<int,EVPair<LatticeColorVector> >();
      }

      //! Callback function
      QDP::MapObject<KeyPropColorVec_t,LatticeFermion>* createMapObjKeyPropColorVecLF(XMLReader& xml_in,
										      const std::string& path) 
      {
	// Doesn't need parameters...
	return new QDP::MapObjectNull<KeyPropColorVec_t,LatticeFermion>();
      }

      //! Local registration flag
      bool registered = false;

      //! Name to be used
      const std::string name = "MAP_OBJECT_NULL";
    } // namespace anonymous

    std::string getName() {return name;}

    //! Register all the factories
    bool registerAll() 
    {
      bool success = true; 
      if (! registered)
      {
	success &= Chroma::TheMapObjIntKeyColorEigenVecFactory::Instance().registerObject(name, createMapObjIntKeyCV);
	success &= Chroma::TheMapObjKeyPropColorVecFactory::Instance().registerObject(name, createMapObjKeyPropColorVecLF);
	registered = true;
      }
      return success;
    }

  } // Namespace MapObjectNullEnv

} // Chroma

