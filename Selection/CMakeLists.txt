include_directories( $ENV{PANDORA_INC} )
include_directories( $ENV{LARPANDORACONTENT_INC} )

add_subdirectory(SelectionTools)
add_subdirectory(AnalysisTools)
add_subdirectory(job)

art_make( MODULE_LIBRARIES larcorealg_Geometry
                           larcore_Geometry_Geometry_service
                           larsim_Simulation nutools_ParticleNavigation
                           lardataobj_Simulation
			   lardataobj_RecoBase
			   lardataobj_MCBase	
                           lardata_Utilities
			   nusimdata_SimulationBase
                           larevt_Filters
                           larsim_MCCheater_BackTrackerService_service
                           larsim_MCCheater_ParticleInventoryService_service
                           lardataobj_RawData
                           larevt_CalibrationDBI_Providers
                           nutools_RandomUtils_NuRandomService_service
                           ${ART_FRAMEWORK_CORE}
                           ${ART_FRAMEWORK_PRINCIPAL}
                           ${ART_FRAMEWORK_SERVICES_REGISTRY}
                           ${ART_FRAMEWORK_SERVICES_OPTIONAL}
                           ${ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR_SERVICE}
                           ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
                           art_Persistency_Common
                           art_Persistency_Provenance
                           art_Utilities
                           canvas
                           ${MF_MESSAGELOGGER}
                           ${MF_UTILITIES}
                           ${FHICLCPP}
                           ${CETLIB}
                           cetlib_except
                           ${CLHEP}
                           ${ROOT_GEOM}
                           ${ROOT_XMLIO}
                           ${ROOT_GDML}
                           ${ROOT_BASIC_LIB_LIST}
            		   ${PANDORASDK}	
			   larpandora_LArPandoraInterface
        )

install_headers()
install_fhicl()
install_source()
