#include "SLDNuCovMat.h"
using namespace lcio ;
using namespace marlin ;
using namespace std ;

SLDNuCovMat aSLDNuCovMat ;

SLDNuCovMat::SLDNuCovMat() :
Processor("SLDNuCovMat"),
m_nRun(0),
m_nEvt(0),
m_nRunSum(0),
m_nEvtSum(0)
{
	_description = "SLDNuCovMat updates the covariance matrix of the neutrino from semi-leptonic decays" ;
	registerInputCollection( 	LCIO::RECONSTRUCTEDPARTICLE,
								"JetCollectionName" ,
								"Name of the Jet collection"  ,
								m_inputJetCollection ,
								std::string("Durham_2Jets")
							);

	registerInputCollection(	LCIO::LCRELATION,
								"SLDNuLinkName",
								"Name of the NeutrinoSemiLeptonicDecayLinkName output collection",
								m_SLDNuLinkName,
								std::string("SLDNuLinkName")
							);

	registerInputCollection(	LCIO::LCRELATION,
								"JetSLDLinkName",
								"Name of the JetSemiLeptonicDecayLinkName output collection",
								m_JetSLDLinkName,
								std::string("JetSLDLinkName")
							);

	registerInputCollection(	LCIO::LCRELATION,
								"recoNumcNuLinkName",
								"Name of the trueNeutrino-reconstructedNeutrino output Link collection",
								m_recoNumcNuLinkName,
								std::string("recoNumcNuLinkName")
							);

	registerProcessorParameter(	"NuCovMatScenarion",
								"Scenario for updating Covariance Matrix of neutrino: (0) do not update, (1) Fill by 0s, (2) Fill by Residuals",
								m_NuCovMatScenarion,
								int(0)
							);


}

void SLDNuCovMat::init()
{
	streamlog_out(DEBUG) << "	init called  " << std::endl;
	printParameters();
}

void SLDNuCovMat::processEvent( EVENT::LCEvent *pLCEvent )
{
	LCCollection *inputJetCollection{};
	int n_Jets = -1;
	streamlog_out(MESSAGE) << "" << std::endl;
	streamlog_out(MESSAGE) << "	////////////////////////////////////////////////////////////////////////////" << std::endl;
	streamlog_out(MESSAGE) << "	////////////////////	Processing event 	" << pLCEvent->getEventNumber() << "	////////////////////" << std::endl;
	streamlog_out(MESSAGE) << "	////////////////////////////////////////////////////////////////////////////" << std::endl;
	try
	{
		inputJetCollection = pLCEvent->getCollection( m_inputJetCollection );
		LCRelationNavigator JetSLDNav( pLCEvent->getCollection( m_JetSLDLinkName ) );
		LCRelationNavigator SLDNuNav( pLCEvent->getCollection( m_SLDNuLinkName ) );
		LCRelationNavigator recoNumcNuNav( pLCEvent->getCollection( m_recoNumcNuLinkName ) );

		n_Jets = inputJetCollection->getNumberOfElements();
		streamlog_out(DEBUG4) << "	Input Jet collection (" << m_inputJetCollection << ") has " << n_Jets << " element (Jets) " << std::endl;
		if ( n_Jets == -1 )
		{
			return;
		}
		for ( int i_jet = 0; i_jet < n_Jets ; ++i_jet )
		{
			ReconstructedParticle* jet = dynamic_cast<ReconstructedParticle*>( inputJetCollection->getElementAt( i_jet ) );
			const EVENT::LCObjectVec& SLDVertices = JetSLDNav.getRelatedToObjects( jet );
			streamlog_out(DEBUG3) << "	Jet[ " << i_jet << " ] has " << SLDVertices.size() << " solved semi-leptonic decays" << std::endl;
			for ( unsigned int i_sld = 0 ; i_sld < SLDVertices.size() ; ++i_sld )
			{
				Vertex *sldVertex = (Vertex*) SLDVertices.at( i_sld );
				const EVENT::LCObjectVec& neutrinos = SLDNuNav.getRelatedToObjects( sldVertex );
				streamlog_out(DEBUG2) << "	SLD[ " << i_sld << " ] has " << neutrinos.size() << " neutrino solutions" << std::endl;
				for ( unsigned int i_nu = 0 ; i_nu < neutrinos.size() ; ++i_nu )
				{
					streamlog_out(DEBUG1) << "	Setting Covariance Matrix of neutrino[ " << i_nu << " ]" << std::endl;
					ReconstructedParticleImpl* neutrino = dynamic_cast<ReconstructedParticleImpl*>( neutrinos.at( i_nu ) );
					std::vector<float> NuCovMat( 10 , 0.0 );
					const EVENT::LCObjectVec& linkedTrueNeutrinos = recoNumcNuNav.getRelatedToObjects( neutrinos.at( 0 ) );
					EVENT::MCParticle* linkedTrueNeutrino = (EVENT::MCParticle*) linkedTrueNeutrinos.at( 0 );
					TLorentzVector recoNeutrinoFourMomentum = TLorentzVector( neutrino->getMomentum() , neutrino->getEnergy() );
					TLorentzVector trueNeutrinoFourMomentum = TLorentzVector( linkedTrueNeutrino->getMomentum() , linkedTrueNeutrino->getEnergy() );
					TLorentzVector fourMomentumResidual = recoNeutrinoFourMomentum - trueNeutrinoFourMomentum;
					streamlog_out(DEBUG1) << "	true Four-momentum of neutrino[ " << i_nu << " ]: ( " << trueNeutrinoFourMomentum.Px() << " , " << trueNeutrinoFourMomentum.Py() << " , " << trueNeutrinoFourMomentum.Pz() << " , " << trueNeutrinoFourMomentum.E() << " )" << std::endl;
					streamlog_out(DEBUG1) << "	reco Four-momentum of neutrino[ " << i_nu << " ]: ( " << recoNeutrinoFourMomentum.Px() << " , " << recoNeutrinoFourMomentum.Py() << " , " << recoNeutrinoFourMomentum.Pz() << " , " << recoNeutrinoFourMomentum.E() << " )" << std::endl;
					streamlog_out(DEBUG1) << "	Initial Covariance Matrix of neutrino[ " << i_nu << " ]" << std::endl;
					for( int i_element = 0 ; i_element < 10 ; ++i_element ) streamlog_out(DEBUG1) << "			CovMat[ " << i_element << " ] = " << neutrino->getCovMatrix()[ i_element ] << std::endl;
					if ( i_nu == 0 )
					{
						neutrino->setCovMatrix( NuCovMat );
					}
					else if ( m_NuCovMatScenarion == 1 )
					{
						neutrino->setCovMatrix( NuCovMat );
					}
					else if ( m_NuCovMatScenarion == 2 )
					{
						NuCovMat[ 0 ] = fourMomentumResidual.Px() * fourMomentumResidual.Px();
						NuCovMat[ 1 ] = fourMomentumResidual.Py() * fourMomentumResidual.Px();
						NuCovMat[ 2 ] = fourMomentumResidual.Py() * fourMomentumResidual.Py();
						NuCovMat[ 3 ] = fourMomentumResidual.Pz() * fourMomentumResidual.Px();
						NuCovMat[ 4 ] = fourMomentumResidual.Pz() * fourMomentumResidual.Py();
						NuCovMat[ 5 ] = fourMomentumResidual.Pz() * fourMomentumResidual.Pz();
						NuCovMat[ 6 ] = fourMomentumResidual.E() * fourMomentumResidual.Px();
						NuCovMat[ 7 ] = fourMomentumResidual.E() * fourMomentumResidual.Py();
						NuCovMat[ 8 ] = fourMomentumResidual.E() * fourMomentumResidual.Pz();
						NuCovMat[ 9 ] = fourMomentumResidual.E() * fourMomentumResidual.E();
						neutrino->setCovMatrix( NuCovMat );
					}
					streamlog_out(DEBUG1) << "	Covariance Matrix of neutrino[ " << i_nu << " ]" << std::endl;
					for( int i_element = 0 ; i_element < 10 ; ++i_element ) streamlog_out(DEBUG1) << "			CovMat[ " << i_element << " ] = " << neutrino->getCovMatrix()[ i_element ] << std::endl;
					streamlog_out(DEBUG1) << "	Covariance Matrix of neutrino[ " << i_nu << " ] is set succesfully" << std::endl;
				}
			}
		}
	}
	catch(DataNotAvailableException &e)
	{
		streamlog_out(MESSAGE) << "Input collection not found in event " << pLCEvent->getEventNumber() << std::endl;
	}

}

void SLDNuCovMat::check()
{

}

void SLDNuCovMat::end()
{

}