#ifndef SLDNuCovMat_h
#define SLDNuCovMat_h 1

#include <iostream>
#include <vector>
#include <string>

#include "marlin/Processor.h"
#include "lcio.h"
#include <EVENT/MCParticle.h>
#include <EVENT/Vertex.h>
#include <EVENT/ReconstructedParticle.h>
#include "TLorentzVector.h"
#include "IMPL/ReconstructedParticleImpl.h"
#include "IMPL/LCCollectionVec.h"
#include <EVENT/LCCollection.h>
#include "UTIL/LCRelationNavigator.h"
using namespace lcio ;
using namespace marlin ;

class SLDNuCovMat : public Processor
{
public:
	virtual Processor *newProcessor()
	{
		return new SLDNuCovMat;
	}
	SLDNuCovMat();
	virtual ~SLDNuCovMat() = default;
	SLDNuCovMat( const SLDNuCovMat& ) = delete;
	SLDNuCovMat &operator = ( const SLDNuCovMat& ) = delete;
	virtual void init();
	virtual void processEvent( EVENT::LCEvent *pLCEvent );

	typedef std::vector<int>				IntVector;
	typedef std::vector<double>				DoubleVector;
	typedef std::vector<float>				FloatVector;
	typedef std::vector<EVENT::ReconstructedParticle*>	pfoVector;
	typedef std::vector<EVENT::Vertex*>			vtxVector;
	virtual void check();
	virtual void end();

private:

	std::string				m_inputJetCollection{};
	std::string				m_JetSLDLinkName{};
	std::string				m_SLDNuLinkName{};
	std::string				m_recoNumcNuLinkName{};
	int 					m_NuCovMatScenarion{};

	int						m_nRun;
	int						m_nEvt;
	int						m_nRunSum;
	int						m_nEvtSum;

};

#endif
