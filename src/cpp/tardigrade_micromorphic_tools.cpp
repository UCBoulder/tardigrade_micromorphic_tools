/*!
 * Micromorphic Tools.cpp
 * ====================
 *
 * A collection of tools and utilities which seek to make the 
 * implementation of micromorphic continuum mechanics based constitutive 
 * theories easier.
 *
 */

#include<tardigrade_micromorphic_tools.h>

namespace tardigradeMicromorphicTools{

    void computePsi( const variableVector &deformationGradient, const variableVector &microDeformation,
                         variableVector &Psi ){
        /*!
         * Computes the micromorphic quantity psi defined as:
         * \f$\Psi_{IJ} = F_{i I} \Chi_{i J}\f$
         *
         * \param &deformationGradient: The deformation gradient.
         * \param &microDeformation: The micro-deformation.
         * \param &Psi: The micro-displacement metric Psi.
         */

        //Assume 3d
        constexpr unsigned int dim = 3;
        constexpr unsigned int sot_dim = dim * dim;

        TARDIGRADE_ERROR_TOOLS_CHECK( deformationGradient.size() == sot_dim, "The deformation gradient doesn't have the correct size" );

        TARDIGRADE_ERROR_TOOLS_CHECK( microDeformation.size() == sot_dim, "The micro-deformation doesn't have the correct size" );

        Psi = variableVector( sot_dim, 0 );

        Eigen::Map< const Eigen::Matrix< variableType, dim, dim, Eigen::RowMajor > > map1( deformationGradient.data( ), dim, dim );
        Eigen::Map< const Eigen::Matrix< variableType, dim, dim, Eigen::RowMajor > > map2( microDeformation.data( ), dim, dim );
        Eigen::Map< Eigen::Matrix< variableType, dim, dim, Eigen::RowMajor > > map3( Psi.data( ), dim, dim );

        map3 = ( map1.transpose( ) * map2 ).eval( );

        return;
    }

    void computePsi( const variableVector &deformationGradient, const variableVector &microDeformation,
                         variableVector &Psi, variableMatrix &dPsidF, variableMatrix &dPsidChi ){
        /*!
         * Computes the micromorphic quantity psi defined as:
         * \f$\Psi_{IJ} = F_{i I} \Chi_{i J}\f$
         *
         * along with the jacobians
         *
         * \f$\frac{ \partial \Psi_{IJ} }{ \partial F_{k K} } = \delta_{I K} \Chi_{k J}\f$
         * 
         * \f$\frac{ \partial \Psi_{IJ} }{ \partial \Chi_{k K} } = F_{k I} \delta_{J K}\f$
         *
         * \param &deformationGradient: The deformation gradient.
         * \param &microDeformation: The micro-deformation.
         * \param &Psi: The micro-displacement metric Psi
         * \param &dPsidF: The jacobian of Psi w.r.t. the deformation gradient.
         * \param &dPsidChi: The jacobian of Psi w.r.t. the micro-deformation.
         */

        constexpr unsigned int dim = 3;
        constexpr unsigned int sot_dim = dim * dim;

        variableVector _dPsidF;
        variableVector _dPsidChi;

        TARDIGRADE_ERROR_TOOLS_CATCH( computePsi( deformationGradient, microDeformation, Psi, _dPsidF, _dPsidChi ) );

        dPsidF   = tardigradeVectorTools::inflate( _dPsidF, sot_dim, sot_dim );
        dPsidChi = tardigradeVectorTools::inflate( _dPsidChi, sot_dim, sot_dim );

        return;

    }

    void computePsi( const variableVector &deformationGradient, const variableVector &microDeformation,
                         variableVector &Psi, variableVector &dPsidF, variableVector &dPsidChi ){
        /*!
         * Computes the micromorphic quantity psi defined as:
         * \f$\Psi_{IJ} = F_{i I} \Chi_{i J}\f$
         *
         * along with the jacobians
         *
         * \f$\frac{ \partial \Psi_{IJ} }{ \partial F_{k K} } = \delta_{I K} \Chi_{k J}\f$
         * 
         * \f$\frac{ \partial \Psi_{IJ} }{ \partial \Chi_{k K} } = F_{k I} \delta_{J K}\f$
         *
         * \param &deformationGradient: The deformation gradient.
         * \param &microDeformation: The micro-deformation.
         * \param &Psi: The micro-displacement metric Psi
         * \param &dPsidF: The jacobian of Psi w.r.t. the deformation gradient.
         * \param &dPsidChi: The jacobian of Psi w.r.t. the micro-deformation.
         */

        //Assume 3d
        constexpr unsigned int dim = 3;
        constexpr unsigned int sot_dim = dim * dim;

        TARDIGRADE_ERROR_TOOLS_CATCH( computePsi( deformationGradient, microDeformation, Psi ) );

        dPsidF   = variableVector( sot_dim * sot_dim, 0 );
        dPsidChi = variableVector( sot_dim * sot_dim, 0 );

        for ( unsigned int I = 0; I < dim; I++ ){
            for ( unsigned int J = 0; J < dim; J++ ){
                for ( unsigned int k = 0; k < dim; k++ ){
                    dPsidF[ dim * sot_dim * I + sot_dim * J + dim * k + I ] = microDeformation[ dim * k + J ];
                    dPsidChi[ dim * sot_dim * I + sot_dim * J + dim * k + J ] = deformationGradient[ dim * k + I ];
                }
            }
        }
        return;
    }

    void computeGamma( const variableVector &deformationGradient, const variableVector &gradChi,
                           variableVector &Gamma ){
        /*!
         * Compute the deformation metric Gamma:
         *
         * \f$ \Gamma_{IJK} = F_{iI} \Chi_{iJ,K} \f$
         *
         * \param &deformationGradient: The deformation gradient.
         * \param &gradChi: The gradient of the micro-deformation tensor
         * \   w.r.t. the reference configuration.
         * \param &Gamma: The micromorphic deformation metric Gamma.
         */

        //Assume 3d
        constexpr unsigned int dim = 3;
        constexpr unsigned int sot_dim = dim * dim;
        constexpr unsigned int tot_dim = sot_dim * dim;

        TARDIGRADE_ERROR_TOOLS_CHECK( deformationGradient.size() == sot_dim, "The deformation gradient isn't the right size");

        TARDIGRADE_ERROR_TOOLS_CHECK( gradChi.size() == tot_dim, "The micro-deformation gradient isn't the right size");

        Gamma = variableVector( tot_dim, 0 );

        for ( unsigned int i = 0; i < dim; i++ ){
            for ( unsigned int I = 0; I < dim; I++ ){
                for ( unsigned int J = 0; J < dim; J++ ){
                    for ( unsigned int K = 0; K < dim; K++ ){
                        Gamma[ dim * dim * I + dim * J + K ] += deformationGradient[ dim * i + I ] * gradChi[ dim * dim * i + dim * J + K ];
                    }
                }
            }
        }

        return;
    }

    void computeGamma( const variableVector &deformationGradient, const variableVector &gradChi,
                           variableVector &Gamma, variableMatrix &dGammadF, variableMatrix &dGammadGradChi ){
        /*!
         * Compute the deformation metric Gamma:
         *
         * \f$ Gamma_{IJK} = F_{iI} \Chi_{iJ,K} \f$
         *
         * Also return the Jacobians
         * \f$ \frac{ \partial Gamma_{IJK} }{ \partial F_{lL} } = \delta_{IL} \Chi_{lJ,K} \f$
         * 
         * \f$ \frac{ \partial Gamma_{IJK} }{ \partial \Chi_{lL,M} } = F_{lI} \delta_{JL} \delta_{KM} \f$
         *
         * \param &deformationGradient: The deformation gradient.
         * \param &gradChi: The gradient of the micro-deformation tensor
         *     w.r.t. the reference configuration.
         * \param &Gamma: The micromorphic deformation metric Gamma.
         * \param &dGammadF: The gradient of Gamma w.r.t. the deformation gradient.
         * \param &dGammadGradChi: The gradient of Gamma w.r.t. the gradient of Chi in the reference 
         *     configuration.
         */

        constexpr unsigned int dim = 3;
        constexpr unsigned int sot_dim = dim * dim;
        constexpr unsigned int tot_dim = sot_dim * dim;

        variableVector _dGammadF;
        variableVector _dGammadGradChi;

        TARDIGRADE_ERROR_TOOLS_CATCH( computeGamma( deformationGradient, gradChi, Gamma, _dGammadF, _dGammadGradChi ) );

        dGammadF       = tardigradeVectorTools::inflate( _dGammadF, tot_dim, sot_dim );
        dGammadGradChi = tardigradeVectorTools::inflate( _dGammadGradChi, tot_dim, tot_dim ); 

        return;

    }

    void computeGamma( const variableVector &deformationGradient, const variableVector &gradChi,
                           variableVector &Gamma, variableVector &dGammadF, variableVector &dGammadGradChi ){
        /*!
         * Compute the deformation metric Gamma:
         *
         * \f$\Gamma_{IJK} = F_{iI} \Chi_{iJ,K}\f$
         *
         * Also return the Jacobians
         * \f$\frac{ \partial Gamma_{IJK} }{ \partial F_{lL} } = \delta_{IL} \Chi_{lJ,K}\f$
         * 
         * \f$\frac{ \partial Gamma_{IJK} }{ \partial \Chi_{lL,M} } = F_{lI} \delta_{JL} \delta_{KM}\f$
         *
         * \param &deformationGradient: The deformation gradient.
         * \param &gradChi: The gradient of the micro-deformation tensor
         *     w.r.t. the reference configuration.
         * \param &Gamma: The micromorphic deformation metric Gamma.
         * \param &dGammadF: The gradient of Gamma w.r.t. the deformation gradient.
         * \param &dGammadGradChi: The gradient of Gamma w.r.t. the gradient of Chi in the reference 
         *     configuration.
         */

        //Assume 3d
        constexpr unsigned int dim = 3;
        constexpr unsigned int sot_dim = dim * dim;
        constexpr unsigned int tot_dim = sot_dim * dim;

        TARDIGRADE_ERROR_TOOLS_CATCH( computeGamma( deformationGradient, gradChi, Gamma ) );

        dGammadF       = variableVector( tot_dim * sot_dim, 0 );
        dGammadGradChi = variableVector( tot_dim * tot_dim, 0 );

        for ( unsigned int I = 0; I < dim; I++ ){
            for ( unsigned int J = 0; J < dim; J++ ){
                for ( unsigned int K = 0; K < dim; K++ ){
                    for ( unsigned int l = 0; l < dim; l++ ){
                        dGammadF[ dim * dim * sot_dim * I + dim * sot_dim * J + sot_dim * K + dim * l + I ] = gradChi[ dim * dim * l + dim * J + K ];
                        dGammadGradChi[ dim * dim * tot_dim * I + dim * tot_dim * J + tot_dim * K + dim * dim * l + dim * J + K ] = deformationGradient[ dim * l + I ];
                    }
                }
            }
        }

        return;
    }

    void computeMicroStrain( const variableVector &Psi, variableVector &microStrain ){
        /*!
         * Compute the microstrain defined as:
         * \f$ \Epsilon_{IJ} = \Psi_{IJ} - eye_{IJ} \f$
         *
         * \param &Psi: The micro-deformation metric Psi.
         * \param &microStrain: The micro-strain.
         */

        //Assume 3d
        constexpr unsigned int dim = 3;
        constexpr unsigned int sot_dim = dim * dim;

        TARDIGRADE_ERROR_TOOLS_CHECK( Psi.size() == sot_dim, "Psi is not of the correct size" );

        microStrain = Psi;
        for ( unsigned int i = 0; i < dim; i++ ){ microStrain[ dim * i + i ] -= 1.; }

        return;
    }

    void computeMicroStrain( const variableVector &Psi, variableVector &microStrain,
                                 variableVector &dMicroStraindPsi ){
        /*!
         * Compute the microstrain defined as:
         * \f$ \Epsilon_{IJ} = \Psi_{IJ} - eye_{IJ} \f$
         *
         * and also compute the jacobian
         * \f$ \frac{ \partial \Epsilon_{IJ} }{ \partial \Psi_{KL} } = \delta_{IK} \delta_{JL} \f$
         *
         * \param &Psi: The micro-deformation metric Psi.
         * \param &microStrain: The micro-strain.
         * \param &dMicroStraindPsi: The jacobian of the micro-strain.
         */

        constexpr unsigned int dim = 3;
        constexpr unsigned int sot_dim = dim * dim;

        TARDIGRADE_ERROR_TOOLS_CATCH( computeMicroStrain( Psi, microStrain ) );

        dMicroStraindPsi = variableVector( sot_dim * sot_dim, 0 );
        for ( unsigned int i = 0; i < sot_dim; i++ ){ dMicroStraindPsi[ sot_dim * i + i ] = 1; };

        return;

    }

    void computeMicroStrain( const variableVector &Psi, variableVector &microStrain,
                                 variableMatrix &dMicroStraindPsi ){
        /*!
         * Compute the microstrain defined as:
         * \f$ \Epsilon_{IJ} = \Psi_{IJ} - eye_{IJ} \f$
         *
         * and also compute the jacobian
         * \f$ \frac{ \partial \Epsilon_{IJ} }{ \partial \Psi_{KL} } = \delta_{IK} \delta_{JL} \f$
         *
         * \param &Psi: The micro-deformation metric Psi.
         * \param &microStrain: The micro-strain.
         * \param &dMicroStraindPsi: The jacobian of the micro-strain.
         */

        constexpr unsigned int dim = 3;
        constexpr unsigned int sot_dim = dim * dim;

        variableVector _dMicroStraindPsi;

        TARDIGRADE_ERROR_TOOLS_CATCH( computeMicroStrain( Psi, microStrain, _dMicroStraindPsi ) );

        dMicroStraindPsi = tardigradeVectorTools::inflate( _dMicroStraindPsi, sot_dim, sot_dim );

        return;

    }

    void pushForwardPK2Stress( const variableVector &PK2Stress,
                                   const variableVector &deformationGradient,
                                   variableVector &cauchyStress ){
        /*!
         * Push forward the PK2-stress in the reference configuration to the 
         * configuration (the Cauchy stress) indicated by the deformation gradient.
         *
         * \f$ \sigma_{ij} = (1 / J ) F_{i I} S_{I J} F_{j J} \f$
         *
         * \param &PK2Stress: The PK2 stress in the 
         *     reference configuration.
         * \param &deformationGradient: The deformation gradient 
         *     mapping between configurations.
         * \param &cauchyStress: The Cauchy stress in the current 
         *     configuration.
         */

        variableType detF;
        TARDIGRADE_ERROR_TOOLS_CATCH( pushForwardReferenceMicroStress( PK2Stress, deformationGradient, 
                                                                       detF, cauchyStress ) );

        return;

    }

    void pushForwardPK2Stress( const variableVector &PK2Stress,
                                   const variableVector &deformationGradient,
                                   variableVector &cauchyStress,
                                   variableVector &dCauchyStressdPK2Stress,
                                   variableVector &dCauchyStressdDeformationGradient ){
        /*!
         * Push forward the PK2 stress in the reference configuration to the 
         * configuration indicated by the deformation gradient.
         *
         * \f$ \sigma_{ij} = (1 / J ) F_{iI} S_{IJ} F_{jJ} \f$
         *
         * Also computes the jacobians:
         * \f$ \frac{ \partial \cauchy_{ij} }{\partial \Sigma_{KL} } = ( 1 / J ) F_{iK} F_{jL} \f$
         *
         * \f$ \frac{ \partial \cauchy_{ij} }{\partial F_{kK} } = ( \delta_{i k} \delta_{I K} S_{I J} F_{j J}
         *                                                  + F_{i I} S_{I J} \delta_{j k} \delta_{J K}
         *                                                  - \cauchy_{i j} dDetFdF_{kK} ) / J \f$
         *
         * \param &PK2Stress: The PK2 stress in the 
         *     reference configuration.
         * \param &deformationGradient: The deformation gradient 
         *     mapping between configurations.
         * \param &cauchyStress: The Cauchy stress in the current 
         *     configuration.
         * \param &dCauchyStressdPK2Stress: The jacobian of 
         *     the Cauchy w.r.t. the PK2 tress in the reference configuration.
         * \param &dCauchyStressdDeformationGradient: The jacobian of 
         *     the Cauchy stress w.r.t. the deformation gradient.
         */

        TARDIGRADE_ERROR_TOOLS_CATCH( pushForwardReferenceMicroStress( PK2Stress, deformationGradient, cauchyStress,
                                                                       dCauchyStressdPK2Stress,
                                                                       dCauchyStressdDeformationGradient ) );

        return;

    }


    void pushForwardPK2Stress( const variableVector &PK2Stress,
                                   const variableVector &deformationGradient,
                                   variableVector &cauchyStress,
                                   variableMatrix &dCauchyStressdPK2Stress,
                                   variableMatrix &dCauchyStressdDeformationGradient ){
        /*!
         * Push forward the PK2 stress in the reference configuration to the 
         * configuration indicated by the deformation gradient.
         *
         * \f$ \sigma_{ij} = (1 / J ) F_{iI} S_{IJ} F_{jJ} \f$
         *
         * Also computes the jacobians:
         * \f$ \frac{ \partial \cauchy_{ij} }{\partial \Sigma_{KL} } = ( 1 / J ) F_{iK} F_{jL} \f$
         * 
         * \f$ \frac{ \partial \cauchy_{ij} }{\partial F_{kK} } = ( \delta_{i k} \delta_{I K} S_{I J} F_{j J}
         *                                                  + F_{i I} S_{I J} \delta_{j k} \delta_{J K}
         *                                                  - \cauchy_{i j} dDetFdF_{kK} ) / J \f$
         *
         * \param &PK2Stress: The PK2 stress in the 
         *     reference configuration.
         * \param &deformationGradient: The deformation gradient 
         *     mapping between configurations.
         * \param &cauchyStress: The Cauchy stress in the current 
         *     configuration.
         * \param &dCauchyStressdPK2Stress: The jacobian of 
         *     the Cauchy w.r.t. the PK2 tress in the reference configuration.
         * \param &dCauchyStressdDeformationGradient: The jacobian of 
         *     the Cauchy stress w.r.t. the deformation gradient.
         */

        TARDIGRADE_ERROR_TOOLS_CATCH( pushForwardReferenceMicroStress( PK2Stress, deformationGradient, cauchyStress,
                                                                       dCauchyStressdPK2Stress,
                                                                       dCauchyStressdDeformationGradient ) );
        return;

    }

    void pullBackCauchyStress( const variableVector &cauchyStress,
                                   const variableVector &deformationGradient,
                                   variableVector &PK2Stress ){
        /*!
         * Pull back the Cauchy stress in the configuration indicated by the deformation gradient
         * to the PK2 stress.
         *
         * \f$ S_{IJ} = J F_{Ii}^{-1} \sigma_{ij} F_{Jj}^{-1} \f$
         *
         * \param &cauchyStress: The Cauchy stress in the current configuration of
         *     the provided deformation gradient.
         * \param &deformationGradient: The deformation gradient mapping between the 
         *     reference configuration and the current configuration.
         * \param &PK2Stress: The PK2 stress in the reference configuration.
         */

        TARDIGRADE_ERROR_TOOLS_CATCH( pullBackMicroStress( cauchyStress, deformationGradient, PK2Stress ) );

        return;

    }

    void pullBackCauchyStress( const variableVector &cauchyStress,
                                   const variableVector &deformationGradient,
                                   variableVector &PK2Stress, variableVector &dPK2StressdCauchyStress,
                                   variableVector &dPK2StressdDeformationGradient ){
        /*!
         * Pull back the Cauchy stress in the configuration indicated by the deformation gradient
         * to the PK2 stress.
         *
         * \f$ S_{IJ} = J F_{Ii}^{-1} \sigma_{ij} F_{Jj}^{-1} \f$
         *
         * Also computes the Jacobians
         *
         * \f$ \frac{ \partial S_{IJ} }{ \partial \sigma_{kl} } = J F_{Ik}^{-1} F_{Jl}^{-1} \f$
         * 
         * \f$ \frac{ \partial S_{IJ} }{ \partial F_{kK} } = F_{Kk}^{-1} S_{IJ} - F_{Ik}^{-1} S_{KJ} - S_{IK} F_{Jk}^{-1} \f$
         *
         * \param &cauchyStress: The Cauchy stress in the current configuration of
         *     the provided deformation gradient.
         * \param &deformationGradient: The deformation gradient mapping between the 
         *     reference configuration and the current configuration.
         * \param &PK2Stress: The PK2 stress in the reference configuration.
         * \param &dPK2StressdCauchyStress: The derivative of the PK2 stress w.r.t. the Cauchy stress
         * \param &dPK2StressdDeformationGradient: The derivative of the PK2 stress w.r.t. the deformation gradient
         */

        TARDIGRADE_ERROR_TOOLS_CATCH( pullBackMicroStress( cauchyStress, deformationGradient, PK2Stress,
                                                           dPK2StressdCauchyStress, dPK2StressdDeformationGradient ) );

        return;

    }


    void pullBackCauchyStress( const variableVector &cauchyStress,
                                   const variableVector &deformationGradient,
                                   variableVector &PK2Stress, variableMatrix &dPK2StressdCauchyStress,
                                   variableMatrix &dPK2StressdDeformationGradient ){
        /*!
         * Pull back the Cauchy stress in the configuration indicated by the deformation gradient
         * to the PK2 stress.
         *
         * \f$ S_{IJ} = J F_{Ii}^{-1} \sigma_{ij} F_{Jj}^{-1} \f$
         *
         * Also computes the Jacobians
         *
         * \f$ \frac{ \partial S_{IJ} }{ \partial \sigma_{kl} } = J F_{Ik}^{-1} F_{Jl}^{-1} \f$
         * 
         * \f$ \frac{ \partial S_{IJ} }{ \partial F_{kK} } = F_{Kk}^{-1} S_{IJ} - F_{Ik}^{-1} S_{KJ} - S_{IK} F_{Jk}^{-1} \f$
         *
         * \param &cauchyStress: The Cauchy stress in the current configuration of
         *     the provided deformation gradient.
         * \param &deformationGradient: The deformation gradient mapping between the 
         *     reference configuration and the current configuration.
         * \param &PK2Stress: The PK2 stress in the reference configuration.
         * \param &dPK2StressdCauchyStress: The derivative of the PK2 stress w.r.t. the Cauchy stress
         * \param &dPK2StressdDeformationGradient: The derivative of the PK2 stress w.r.t. the deformation gradient
         */

        TARDIGRADE_ERROR_TOOLS_CATCH( pullBackMicroStress( cauchyStress, deformationGradient, PK2Stress,
                                                           dPK2StressdCauchyStress, dPK2StressdDeformationGradient ) );

        return;

    }

    void pushForwardReferenceMicroStress( const variableVector &referenceMicroStress,
                                              const variableVector &deformationGradient,
                                              variableVector &microStress ){
        /*!
         * Push forward the micro-stress in the reference configuration to the 
         * configuration indicated by the deformation gradient.
         *
         * \f$ s_{ij} = (1 / J ) F_{i I} \Sigma_{I J} F_{j J} \f$
         *
         * \param &referenceMicroStress: The micro-stress in the 
         *     reference configuration.
         * \param &deformationGradient: The deformation gradient 
         *     mapping between configurations.
         * \param &microStress: The micro-stress in the current 
         *     configuration.
         */

        variableType detF;
        return pushForwardReferenceMicroStress( referenceMicroStress, deformationGradient, 
                                                detF, microStress );

    }

    void pushForwardReferenceMicroStress( const variableVector &referenceMicroStress,
                                              const variableVector &deformationGradient,
                                              variableType &detF, variableVector &microStress ){
        /*!
         * Push forward the micro-stress in the reference configuration to the 
         * configuration indicated by the deformation gradient.
         *
         * \f$ s_{ij} = (1 / J ) F_{i I} \Sigma_{I J} F_{j J} \f$
         *
         * \param &referenceMicroStress: The micro-stress in the 
         *     reference configuration.
         * \param &deformationGradient: The deformation gradient 
         *     mapping between configurations.
         * \param &detF: The determinant of the deformation gradient.
         * \param &microStress: The micro-stress in the current 
         *     configuration.
         */

        //Assume 3d
        constexpr unsigned int dim = 3;
        constexpr unsigned int sot_dim = dim * dim;

        TARDIGRADE_ERROR_TOOLS_CHECK( referenceMicroStress.size() == sot_dim, "The reference micro-stress has an incorrect size");

        TARDIGRADE_ERROR_TOOLS_CHECK( deformationGradient.size() == sot_dim, "The deformation gradient has an incorrect size");

        microStress = variableVector( sot_dim, 0 );

        Eigen::Map< const Eigen::Matrix< variableType, dim, dim, Eigen::RowMajor > > map1( deformationGradient.data( ), dim, dim );
        Eigen::Map< const Eigen::Matrix< variableType, dim, dim, Eigen::RowMajor > > map2( referenceMicroStress.data( ), dim, dim );
        Eigen::Map< Eigen::Matrix< variableType, dim, dim, Eigen::RowMajor > > map3( microStress.data( ), dim, dim );

        detF = map1.determinant( );

        map3 = ( map1 * map2 * map1.transpose( ) / detF ).eval( );

        return;
    }

    void pushForwardReferenceMicroStress( const variableVector &referenceMicroStress,
                                              const variableVector &deformationGradient,
                                              variableVector &microStress,
                                              variableMatrix &dMicroStressdReferenceMicroStress,
                                              variableMatrix &dMicroStressdDeformationGradient ){
        /*!
         * Push forward the micro-stress in the reference configuration to the 
         * configuration indicated by the deformation gradient.
         *
         * \f$ s_{ij} = (1 / J ) F_{iI} \Sigma_{IJ} F_{jJ} \f$
         *
         * Also computes the jacobians:
         * \f$ \frac{ \partial s_{ij} }{\partial \Sigma_{KL} } = ( 1 / J ) F_{iK} F_{jL} \f$
         * 
         * \f$ \frac{ \partial s_{ij} }{\partial F_{kK} } = ( \delta_{i k} \delta_{I K} \Sigma_{I J} F_{j J}
         *                                              + F_{i I} \Sigma_{I J} \delta_{j k} \delta_{J K}
         *                                              - s_{i j} dDetFdF_{kK} ) / J \f$
         *
         * \param &referenceMicroStress: The micro-stress in the 
         *     reference configuration.
         * \param deformationGradient: The deformation gradient 
         *     mapping between configurations.
         * \param microStress: The micro-stress in the current 
         *     configuration.
         * \param dMicroStressdReferenceMicroStress: The jacobian of 
         *     the micro-stress w.r.t. the micro-stress in the reference configuration.
         * \param dMicroStressdDeformationGradient: The jacobian of 
         *     the micro-stress w.r.t. the deformation gradient.
         */

        variableVector _dMicroStressdReferenceMicroStress;
        variableVector _dMicroStressdDeformationGradient;

        TARDIGRADE_ERROR_TOOLS_CATCH( pushForwardReferenceMicroStress( referenceMicroStress, deformationGradient,
                                                                       microStress, _dMicroStressdReferenceMicroStress,
                                                                       _dMicroStressdDeformationGradient ) );

        dMicroStressdReferenceMicroStress = tardigradeVectorTools::inflate( _dMicroStressdReferenceMicroStress, 9, 9 );
        dMicroStressdDeformationGradient  = tardigradeVectorTools::inflate( _dMicroStressdDeformationGradient, 9, 9 );

        return;

    }
    void pushForwardReferenceMicroStress( const variableVector &referenceMicroStress,
                                              const variableVector &deformationGradient,
                                              variableVector &microStress,
                                              variableVector &dMicroStressdReferenceMicroStress,
                                              variableVector &dMicroStressdDeformationGradient ){
        /*!
         * Push forward the micro-stress in the reference configuration to the 
         * configuration indicated by the deformation gradient.
         *
         * \f$ s_{ij} = (1 / J ) F_{iI} \Sigma_{IJ} F_{jJ} \f$
         *
         * Also computes the jacobians:
         * \f$ \frac{ \partial s_{ij} }{\partial \Sigma_{KL} } = ( 1 / J ) F_{iK} F_{jL} \f$
         * 
         * \f$ \frac{ \partial s_{ij} }{\partial F_{kK} } = ( \delta_{i k} \delta_{I K} \Sigma_{I J} F_{j J}
         *                                              + F_{i I} \Sigma_{I J} \delta_{j k} \delta_{J K}
         *                                              - s_{i j} dDetFdF_{kK} ) / J \f$
         *
         * \param &referenceMicroStress: The micro-stress in the 
         *     reference configuration.
         * \param deformationGradient: The deformation gradient 
         *     mapping between configurations.
         * \param microStress: The micro-stress in the current 
         *     configuration.
         * \param dMicroStressdReferenceMicroStress: The jacobian of 
         *     the micro-stress w.r.t. the micro-stress in the reference configuration.
         * \param dMicroStressdDeformationGradient: The jacobian of 
         *     the micro-stress w.r.t. the deformation gradient.
         */

        //Assume 3d
        constexpr unsigned int dim = 3;
        constexpr unsigned int sot_dim = dim * dim;

        variableType detF;
        TARDIGRADE_ERROR_TOOLS_CATCH( pushForwardReferenceMicroStress( referenceMicroStress, deformationGradient,
                                                                       detF, microStress ) );

        //Assemble the inverse deformation gradient
        variableVector inverseDeformationGradient = deformationGradient;
        Eigen::Map< Eigen::Matrix< variableType, dim, dim, Eigen::RowMajor > > map( inverseDeformationGradient.data( ), dim, dim );
        map = map.inverse( ).eval( );

        //Assemble the jacobians
        dMicroStressdReferenceMicroStress = variableVector( sot_dim * sot_dim, 0 );
        dMicroStressdDeformationGradient  = variableVector( sot_dim * sot_dim, 0 );

        variableVector temp_sot1( sot_dim, 0 );
        variableVector temp_sot2( sot_dim, 0 );

        for ( unsigned int i = 0; i < dim; i++ ){
            for ( unsigned int j = 0; j < dim; j++ ){
                for ( unsigned int k = 0; k < dim; k++ ){
                    temp_sot1[ dim * i + j ] += referenceMicroStress[ dim * j + k ] * deformationGradient[ dim * i + k ] / detF;
                    temp_sot2[ dim * i + j ] += referenceMicroStress[ dim * k + j ] * deformationGradient[ dim * i + k ] / detF;
                }
            }
        }

        for ( unsigned int i = 0; i < dim; i++ ){
            for ( unsigned int j = 0; j < dim; j++ ){
                for ( unsigned int k = 0; k < dim; k++ ){
                    dMicroStressdDeformationGradient[ dim * sot_dim * i + sot_dim * j + dim * i + k] += temp_sot1[ dim * j + k ];
                    dMicroStressdDeformationGradient[ dim * sot_dim * i + sot_dim * j + dim * j + k] += temp_sot2[ dim * i + k ];

                    for ( unsigned int K = 0; K < dim; K++ ){
                        dMicroStressdReferenceMicroStress[ dim * sot_dim * i + sot_dim * j + dim * k + K ] = deformationGradient[ dim * i + k ]
                                                                                                           * deformationGradient[ dim * j + K ] / detF;
                        
                        dMicroStressdDeformationGradient[ dim * sot_dim * i + sot_dim * j + dim * k + K] -= microStress[ dim * i + j ] * inverseDeformationGradient[ dim * K + k ];

                    }

                }

            }

        }

        return;
    }

    void pullBackMicroStress( const variableVector &microStress,
                                  const variableVector &deformationGradient,
                                  variableVector &referenceMicroStress ){
        /*!
         * Push forward the micro-stress in the reference configuration to the 
         * configuration indicated by the deformation gradient.
         *
         * \f$ \Sigma_{IJ} = J F_{I i}^{-1} s_{ij} F_{J j}^{-1} \f$
         *
         * \param &microStress: The micro-stress in the current 
         *     configuration.
         * \param &deformationGradient: The deformation gradient 
         *     mapping between configurations.
         * \param &referenceMicroStress: The micro-stress in the 
         *     reference configuration.
         */

        variableType detF;
        variableVector inverseDeformationGradient;
        return pullBackMicroStress( microStress, deformationGradient, 
                                    detF, inverseDeformationGradient, referenceMicroStress );

    }

    void pullBackMicroStress( const variableVector &microStress,
                                  const variableVector &deformationGradient,
                                  variableType &detF, variableVector &inverseDeformationGradient,
                                  variableVector &referenceMicroStress ){
        /*!
         * Push forward the micro-stress in the reference configuration to the 
         * configuration indicated by the deformation gradient.
         *
         * \f$ \Sigma_{IJ} = J F_{I i}^{-1} s_{ij} F_{J j}^{-1} \f$
         *
         * \param &microStress: The micro-stress in the current 
         *     configuration.
         * \param &deformationGradient: The deformation gradient 
         *     mapping between configurations.
         * \param &detF: The determinant of the deformation gradient.
         * \param &inverseDeformationGradient: The inverse of the deformation gradient.
         * \param &referenceMicroStress: The micro-stress in the 
         *     reference configuration.
         */

        //Assume 3d
        constexpr unsigned int dim = 3;
        constexpr unsigned int sot_dim = dim * dim;

        TARDIGRADE_ERROR_TOOLS_CHECK( microStress.size() == sot_dim, "The micro-stress has an incorrect size");

        TARDIGRADE_ERROR_TOOLS_CHECK( deformationGradient.size() == sot_dim, "The deformation gradient has an incorrect size");

        referenceMicroStress = variableVector( sot_dim, 0 );

        Eigen::Map< const Eigen::Matrix< variableType, dim, dim, Eigen::RowMajor > > cmap( deformationGradient.data( ), dim, dim );
        detF = cmap.determinant( );

        inverseDeformationGradient = deformationGradient;
        Eigen::Map< Eigen::Matrix< variableType, dim, dim, Eigen::RowMajor > > map( inverseDeformationGradient.data( ), dim, dim );
        map = map.inverse( ).eval( );

        variableVector temp_sot( sot_dim, 0 );
        for ( unsigned int I = 0; I < dim; I++ ){
            for ( unsigned int i = 0; i < dim; i++ ){
                for ( unsigned int J = 0; J < dim; J++ ){
                    temp_sot[ dim * I + J ] += inverseDeformationGradient[ dim * I + i ] * microStress[ dim * i + J ];
                }
            }
        }

        for ( unsigned int I = 0; I < dim; I++ ){
            for ( unsigned int J = 0; J < dim; J++ ){
                for ( unsigned int i = 0; i < dim; i++ ){
                        referenceMicroStress[ dim * I + J ] += temp_sot[ dim * I + i ] * inverseDeformationGradient[ dim * J + i ];
                }
            }
        }

        referenceMicroStress *= detF;

        return;
    }

    void pullBackMicroStress( const variableVector &microStress,
                                  const variableVector &deformationGradient,
                                  variableVector &referenceMicroStress,
                                  variableMatrix &dReferenceMicroStressdMicroStress,
                                  variableMatrix &dReferenceMicroStressdDeformationGradient ){
        /*!
         * Push forward the micro-stress in the reference configuration to the 
         * configuration indicated by the deformation gradient.
         *
         * \f$ \Sigma_{IJ} = J F_{Ii}^{-1} s_{IJ} F_{Jj}^{-1} \f$
         *
         * Also computes the jacobians:
         * \f$ \frac{ \partial \Sigma_{IJ} }{ \partial s_{kl} } = J F_{Ik}^{-1} F_{Jl}^{-1} \f$
         * 
         * \f$\frac{ \partial \Sigma_{IJ} }{ \partial F_{kK} } = F_{Kk}^{-1} \Sigma_{IJ} - F_{Ik}^{-1} \Sigma_{KJ} - \Sigma_{IK} F_{Jk}^{-1}\f$
         *
         * \param &microStress: The micro-stress in the current 
         *     configuration.
         * \param &deformationGradient: The deformation gradient 
         *     mapping between configurations.
         * \param &referenceMicroStress: The micro-stress in the 
         *     reference configuration.
         * \param &dReferenceMicroStressdMicroStress: The jacobian of 
         *     the reference micro-stress w.r.t. the micro-stress in the reference configuration.
         * \param &dReferenceMicroStressdDeformationGradient: The jacobian of 
         *     the reference micro-stress w.r.t. the deformation gradient.
         */

        variableVector _dReferenceMicroStressdMicroStress;
        variableVector _dReferenceMicroStressdDeformationGradient;

        TARDIGRADE_ERROR_TOOLS_CATCH( pullBackMicroStress( microStress, deformationGradient,
                                                           referenceMicroStress, _dReferenceMicroStressdMicroStress,
                                                           _dReferenceMicroStressdDeformationGradient ) );

        dReferenceMicroStressdMicroStress = tardigradeVectorTools::inflate( _dReferenceMicroStressdMicroStress, 9, 9 );
        dReferenceMicroStressdDeformationGradient  = tardigradeVectorTools::inflate( _dReferenceMicroStressdDeformationGradient, 9, 9 );

        return;


    }

    void pullBackMicroStress( const variableVector &microStress,
                                  const variableVector &deformationGradient,
                                  variableVector &referenceMicroStress,
                                  variableVector &dReferenceMicroStressdMicroStress,
                                  variableVector &dReferenceMicroStressdDeformationGradient ){
        /*!
         * Push forward the micro-stress in the reference configuration to the 
         * configuration indicated by the deformation gradient.
         *
         * \f$\Sigma_{IJ} = J F_{Ii}^{-1} s_{IJ} F_{Jj}^{-1}\f$
         *
         * Also computes the jacobians:
         * \f$\frac{ \partial \Sigma_{IJ} }{ \partial s_{kl} } = J F_{Ik}^{-1} F_{Jl}^{-1}\f$
         * 
         * \f$\frac{ \partial \Sigma_{IJ} }{ \partial F_{kK} } = F_{Kk}^{-1} \Sigma_{IJ} - F_{Ik}^{-1} \Sigma_{KJ} - \Sigma_{IK} F_{Jk}^{-1}\f$
         *
         * \param &microStress: The micro-stress in the current 
         *     configuration.
         * \param &deformationGradient: The deformation gradient 
         *     mapping between configurations.
         * \param &referenceMicroStress: The micro-stress in the 
         *     reference configuration.
         * \param &dReferenceMicroStressdMicroStress: The jacobian of 
         *     the reference micro-stress w.r.t. the micro-stress in the reference configuration.
         * \param &dReferenceMicroStressdDeformationGradient: The jacobian of 
         *     the reference micro-stress w.r.t. the deformation gradient.
         */

        //Assume 3d
        constexpr unsigned int dim = 3;
        constexpr unsigned int sot_dim = dim * dim;

        variableType detF;
        variableVector inverseDeformationGradient;
        TARDIGRADE_ERROR_TOOLS_CATCH( pullBackMicroStress( microStress, deformationGradient,
                                                           detF, inverseDeformationGradient, referenceMicroStress ) );

        //Assemble the jacobians
        dReferenceMicroStressdMicroStress = variableVector( sot_dim * sot_dim, 0 );
        dReferenceMicroStressdDeformationGradient = variableVector( sot_dim * sot_dim, 0 );

        for ( unsigned int I = 0; I < dim; I++ ){
            for ( unsigned int J = 0; J < dim; J++ ){
                for ( unsigned int k = 0; k < dim; k++ ){
                    for ( unsigned int K = 0; K < dim; K++ ){
                        dReferenceMicroStressdMicroStress[ dim * sot_dim * I + sot_dim * J + dim * k + K ]
                            = detF * inverseDeformationGradient[ dim * I + k ] * inverseDeformationGradient[ dim * J + K ];
                        dReferenceMicroStressdDeformationGradient[ dim * sot_dim * I + sot_dim * J + dim * k + K ]
                            = inverseDeformationGradient[ dim * K + k ] * referenceMicroStress[ dim * I + J ]
                            - inverseDeformationGradient[ dim * I + k ] * referenceMicroStress[ dim * K + J ]
                            - inverseDeformationGradient[ dim * J + k ] * referenceMicroStress[ dim * I + K ];
                    }
                }
            }
        }

        return;
    }

    void pushForwardHigherOrderStress( const variableVector &referenceHigherOrderStress,
                                           const variableVector &deformationGradient,
                                           const variableVector &microDeformation,
                                           variableVector &higherOrderStress ){
        /*!
         * Compute the push-forward operation on the higher order stress.
         *
         * \f$m_{ijk} = \frac{1}{J} F_{iI} F_{jJ} \Chi_{kK} M_{IJK}\f$
         *
         * \param &referenceHigherOrderStress: The higher order stress in the 
         *     reference configuration.
         * \param &deformationGradient: The deformation gradient which maps 
         *     between the reference and current configurations.
         * \param &microDeformation: The micro-deformation tensor.
         * \param &higherOrderStress: The higher order stress in the current configuration.
         */

        variableType detF;
        return pushForwardHigherOrderStress( referenceHigherOrderStress, deformationGradient,
                                             microDeformation, detF, higherOrderStress );
    }

    void pushForwardHigherOrderStress( const variableVector &referenceHigherOrderStress,
                                           const variableVector &deformationGradient,
                                           const variableVector &microDeformation,
                                           variableType &detF,
                                           variableVector &higherOrderStress ){
        /*!
         * Compute the push-forward operation on the higher order stress.
         *
         * \f$m_{ijk} = \frac{1}{J} F_{iI} F_{jJ} \Chi_{kK} M_{IJK}\f$
         *
         * \param &referenceHigherOrderStress: The higher order stress in the 
         *     reference configuration.
         * \param &deformationGradient: The deformation gradient which maps 
         *     between the reference and current configurations.
         * \param &microDeformation: The micro-deformation tensor.
         * \param &detF: The determinant of the deformation gradient.
         * \param &higherOrderStress: The higher order stress in the current configuration.
         */

        //Assume 3d
        constexpr unsigned int dim = 3;
        constexpr unsigned int sot_dim = dim * dim;
        constexpr unsigned int tot_dim = sot_dim * dim;

        TARDIGRADE_ERROR_TOOLS_CHECK( referenceHigherOrderStress.size() == tot_dim, "The reference higher order stress doesn't have the correct size" );

        TARDIGRADE_ERROR_TOOLS_CHECK( deformationGradient.size() == sot_dim, "The deformation gradient doesn't have the correct size" );

        TARDIGRADE_ERROR_TOOLS_CHECK( microDeformation.size() == sot_dim, "The micro-deformation doesn't have the correct size" );

        Eigen::Map< const Eigen::Matrix< variableType, dim, dim, Eigen::RowMajor > > map( deformationGradient.data( ), dim, dim );
        detF = map.determinant( );

        higherOrderStress = variableVector( dim * dim * dim, 0 );

        variableVector temp_tot( tot_dim, 0 );
        for ( unsigned int i = 0; i < dim; i++ ){
            for ( unsigned int j = 0; j < dim; j++ ){
                for ( unsigned int k = 0; k < dim; k++ ){
                    for ( unsigned int I = 0; I < dim; I++ ){
                                higherOrderStress[ dim * dim * i + dim * j + k ] += deformationGradient[ dim * i + I ]
                                                                                  * referenceHigherOrderStress[ dim * dim * I + dim * j + k ];
                    }
                }
            }
        }
        higherOrderStress /= detF;

        for ( unsigned int i = 0; i < dim; i++ ){
            for ( unsigned int j = 0; j < dim; j++ ){
                for ( unsigned int k = 0; k < dim; k++ ){
                    for ( unsigned int I = 0; I < dim; I++ ){
                                temp_tot[ dim * dim * i + dim * j + k ] += deformationGradient[ dim * j + I ]
                                                                         * higherOrderStress[ dim * dim * i + dim * I + k ];
                    }
                }
            }
        }

        std::fill( higherOrderStress.begin( ), higherOrderStress.end( ), 0. );

        for ( unsigned int i = 0; i < dim; i++ ){
            for ( unsigned int j = 0; j < dim; j++ ){
                for ( unsigned int k = 0; k < dim; k++ ){
                    for ( unsigned int I = 0; I < dim; I++ ){
                                higherOrderStress[ dim * dim * i + dim * j + k ] += microDeformation[ dim * k + I ]
                                                                                  * temp_tot[ dim * dim * i + dim * j + I ];
                    }
                }
            }
        }

        return;
    }

    void pushForwardHigherOrderStress( const variableVector &referenceHigherOrderStress,
                                           const variableVector &deformationGradient,
                                           const variableVector &microDeformation,
                                           variableVector &higherOrderStress,
                                           variableMatrix &dHigherOrderStressdReferenceHigherOrderStress,
                                           variableMatrix &dHigherOrderStressdDeformationGradient,
                                           variableMatrix &dHigherOrderStressdMicroDeformation ){
        /*!
         * Compute the push-forward operation on the higher order stress.
         *
         * \f$ m_{ijk} = \frac{1}{J} F_{iI} F_{jJ} \Chi_{kK} M_{IJK} \f$
         *
         * Also returns the Jacobians
         *
         * \f$ \frac{ \partial m_{ijk} }{ \partial M_{LMN} } = \frac{1}{J} F_{iL} F_{jM} \Chi_{kN} \f$
         * 
         * \f$ \frac{ \partial m_{ijk} }{ \partial F_{lM} } = \left( \delta_{il} F_{jN} \Chi_{kO} M_{MNO}
         *                                                     + F_{iN} \delta_{jl} \Chi_{kO} M_{NMO}
         *                                                     - m_{ijk} dDetFdF_{lM} \right)/J\f$
         * 
         * \f$ \frac{ \partial m_{ijk} }{ \partial \Chi_{lM} } = \frac{1}{J} F_{iN} F_{jO} \delta_{kl} M_{NOM} \f$
         *
         * \param &referenceHigherOrderStress: The higher order stress in the 
         *     reference configuration.
         * \param &deformationGradient: The deformation gradient which maps 
         *     between the reference and current configurations.
         * \param &microDeformation: The micro-deformation tensor.
         * \param &higherOrderStress: The higher order stress in the current configuration.
         * \param &dHigherOrderStressdReferenceHigherOrderStress: The Jacobian of the higher order stress w.r.t. the reference higher order stress
         * \param &dHigherOrderStressdDeformationGradient: The Jacobian of the higher order stress w.r.t. the deformation gradient
         * \param &dHigherOrderStressdMicroDeformation: The Jacobian of the higher order stress w.r.t. the micro deformation
         */

        variableVector _dHigherOrderStressdReferenceHigherOrderStress;
        variableVector _dHigherOrderStressdDeformationGradient;
        variableVector _dHigherOrderStressdMicroDeformation;

        TARDIGRADE_ERROR_TOOLS_CATCH( pushForwardHigherOrderStress( referenceHigherOrderStress,
                                                                    deformationGradient,
                                                                    microDeformation,
                                                                    higherOrderStress,
                                                                    _dHigherOrderStressdReferenceHigherOrderStress,
                                                                    _dHigherOrderStressdDeformationGradient,
                                                                    _dHigherOrderStressdMicroDeformation ) );

        dHigherOrderStressdReferenceHigherOrderStress = tardigradeVectorTools::inflate( _dHigherOrderStressdReferenceHigherOrderStress , 27, 27 );
        dHigherOrderStressdDeformationGradient        = tardigradeVectorTools::inflate( _dHigherOrderStressdDeformationGradient        , 27,  9 );
        dHigherOrderStressdMicroDeformation           = tardigradeVectorTools::inflate( _dHigherOrderStressdMicroDeformation           , 27,  9 );

        return;

    }
    void pushForwardHigherOrderStress( const variableVector &referenceHigherOrderStress,
                                           const variableVector &deformationGradient,
                                           const variableVector &microDeformation,
                                           variableVector &higherOrderStress,
                                           variableVector &dHigherOrderStressdReferenceHigherOrderStress,
                                           variableVector &dHigherOrderStressdDeformationGradient,
                                           variableVector &dHigherOrderStressdMicroDeformation ){
        /*!
         * Compute the push-forward operation on the higher order stress.
         *
         * \f$ m_{ijk} = \frac{1}{J} F_{iI} F_{jJ} \Chi_{kK} M_{IJK} \f$
         *
         * Also returns the Jacobians
         *
         * \f$ \frac{ \partial m_{ijk} }{ \partial M_{LMN} } = \frac{1}{J} F_{iL} F_{jM} \Chi_{kN} \f$
         * 
         * \f$\frac{ \partial m_{ijk} }{ \partial F_{lM} } = \left( \delta_{il} F_{jN} \Chi_{kO} M_{MNO}
         *                                                     + F_{iN} \delta_{jl} \Chi_{kO} M_{NMO}
         *                                                     - m_{ijk} dDetFdF_{lM} \right)/J \f$
         * 
         * \f$ \frac{ \partial m_{ijk} }{ \partial \Chi_{lM} } = \frac{1}{J} F_{iN} F_{jO} \delta_{kl} M_{NOM} \f$
         *
         * \param &referenceHigherOrderStress: The higher order stress in the 
         *     reference configuration.
         * \param &deformationGradient: The deformation gradient which maps 
         *     between the reference and current configurations.
         * \param &microDeformation: The micro-deformation tensor.
         * \param &higherOrderStress: The higher order stress in the current configuration.
         * \param &dHigherOrderStressdReferenceHigherOrderStress: The Jacobian of the higher order stress w.r.t. the reference higher order stress
         * \param &dHigherOrderStressdDeformationGradient: The Jacobian of the higher order stress w.r.t. the deformation gradient
         * \param &dHigherOrderStressdMicroDeformation: The Jacobian of the higher order stress w.r.t. the micro deformation
         */

        //Assume 3d
        constexpr unsigned int dim = 3;
        constexpr unsigned int sot_dim = dim * dim;
        constexpr unsigned int tot_dim = sot_dim * dim;

        variableType detF;

        TARDIGRADE_ERROR_TOOLS_CATCH( pushForwardHigherOrderStress( referenceHigherOrderStress, deformationGradient, microDeformation,
                                                                    detF, higherOrderStress ) );

        //Assemble the inverse of the deformation gradient
        variableVector inverseDeformationGradient = deformationGradient;
        Eigen::Map< Eigen::Matrix< variableType, dim, dim, Eigen::RowMajor > > map( inverseDeformationGradient.data( ), dim, dim );
        map = map.inverse( ).eval( );

        dHigherOrderStressdReferenceHigherOrderStress = variableVector( tot_dim * tot_dim, 0 );
        dHigherOrderStressdDeformationGradient        = variableVector( tot_dim * sot_dim, 0 );
        dHigherOrderStressdMicroDeformation           = variableVector( tot_dim * sot_dim, 0 );

        variableVector temp_tot1a( tot_dim, 0 );
        variableVector temp_tot2a( tot_dim, 0 );
        variableVector temp_tot3a( tot_dim, 0 );

        for ( unsigned int i = 0; i < dim; i++ ){
            for ( unsigned int j = 0; j < dim; j++ ){
                for ( unsigned int k = 0; k < dim; k++ ){
                    for ( unsigned int l = 0; l < dim; l++ ){
                        temp_tot1a[ dim * dim * i + dim * j + k ]
                            += deformationGradient[ dim * i + l ] * referenceHigherOrderStress[ dim * dim * j + dim * l + k ];
                        temp_tot2a[ dim * dim * i + dim * j + k ]
                            += deformationGradient[ dim * i + l ] * referenceHigherOrderStress[ dim * dim * l + dim * j + k ];
                        temp_tot3a[ dim * dim * i + dim * j + k ]
                            += deformationGradient[ dim * i + l ] * referenceHigherOrderStress[ dim * dim * l + dim * k + j ];
                    }
                }
            }
        }

        temp_tot1a /= detF;
        temp_tot2a /= detF;
        temp_tot3a /= detF;

        for ( unsigned int i = 0; i < dim; i++ ){
            for ( unsigned int j = 0; j < dim; j++ ){
                for ( unsigned int k = 0; k < dim; k++ ){
                    for ( unsigned int l = 0; l < dim; l++ ){
                        for ( unsigned int M = 0; M < dim; M++ ){

                            dHigherOrderStressdDeformationGradient[ dim * dim * sot_dim * i + dim * sot_dim * j +sot_dim * k + dim * i + l ]
                                += temp_tot1a[ dim * dim * j + dim * l + M ] * microDeformation[ dim * k + M ];

                            dHigherOrderStressdDeformationGradient[ dim * dim * sot_dim * i + dim * sot_dim * j +sot_dim * k + dim * j + l ]
                                += temp_tot2a[ dim * dim * i + dim * l + M ] * microDeformation[ dim * k + M ]; 

                            dHigherOrderStressdMicroDeformation[ dim * dim * sot_dim * i + dim * sot_dim * j + sot_dim * k + dim * k + l ]
                                += temp_tot3a[ dim * dim * i + dim * l + M ] * deformationGradient[ dim * j + M ];

                            for ( unsigned int N = 0; N < dim; N++ ){
                                dHigherOrderStressdReferenceHigherOrderStress[ dim * dim * tot_dim * i + dim * tot_dim * j + tot_dim * k + dim * dim * l + dim * M + N ] = deformationGradient[ dim * i + l ] * deformationGradient[ dim * j + M ] * microDeformation[ dim * k + N] / detF;

                            }
                            dHigherOrderStressdDeformationGradient[ dim * dim * sot_dim * i + dim * sot_dim * j + sot_dim * k + dim * l + M ] -= higherOrderStress[ dim * dim * i + dim * j + k ] * inverseDeformationGradient[ dim * M + l ];

                        }
                    }
                }
            }
        }

        return;
    }

    void pullBackHigherOrderStress( const variableVector &higherOrderStress,
                                        const variableVector &deformationGradient,
                                        const variableVector &microDeformation,
                                        variableVector &referenceHigherOrderStress ){
        /*!
         * Compute the pull back operation on the higher order stress.
         *
         * \f$ M_{IJK} = J F_{Ii}^{-1} F_{Jj}^{-1} \chi_{Kk}^{-1} m_{ijk} \f$
         *
         * \param &higherOrderStress: The higher order stress in the current configuration.
         * \param &deformationGradient: The deformation gradient which maps 
         *     between the reference and current configurations.
         * \param &microDeformation: The micro-deformation tensor.
         * \param &referenceHigherOrderStress: The higher order stress in the reference configuration.
         */

        variableType detF;
        variableVector inverseDeformationGradient, inverseMicroDeformation;
        return pullBackHigherOrderStress( higherOrderStress, deformationGradient,
                                          microDeformation, detF, inverseDeformationGradient,
                                          inverseMicroDeformation, referenceHigherOrderStress );
    }

    void pullBackHigherOrderStress( const variableVector &higherOrderStress,
                                        const variableVector &deformationGradient,
                                        const variableVector &microDeformation,
                                        variableType &detF, variableVector &inverseDeformationGradient,
                                        variableVector &inverseMicroDeformation,
                                        variableVector &referenceHigherOrderStress ){
        /*!
         * Compute the push-forward operation on the higher order stress.
         *
         * \f$ M_{IJK} = J F_{Ii}^{-1} F_{Jj}^{-1} \chi_{Kk}^{-1} m_{ijk} \f$
         *
         * \param &higherOrderStress: The higher order stress in the current configuration.
         * \param &deformationGradient: The deformation gradient which maps 
         *     between the reference and current configurations.
         * \param &microDeformation: The micro-deformation tensor.
         * \param &detF: The determinant of the deformation gradient.
         * \param &inverseDeformationGradient: The inverse of the deformation gradient.
         * \param &inverseMicroDeformation: The inverse of the micro deformation.
         * \param &referenceHigherOrderStress: The higher order stress in the reference configuration.
         */

        //Assume 3d
        constexpr unsigned int dim = 3;
        constexpr unsigned int sot_dim = dim * dim;
        constexpr unsigned int tot_dim = sot_dim * dim;

        TARDIGRADE_ERROR_TOOLS_CHECK( higherOrderStress.size() == tot_dim, "The current higher order stress doesn't have the correct size" );

        TARDIGRADE_ERROR_TOOLS_CHECK( deformationGradient.size() == sot_dim, "The deformation gradient doesn't have the correct size" );

        TARDIGRADE_ERROR_TOOLS_CHECK( microDeformation.size() == dim * dim, "The micro-deformation doesn't have the correct size" );

        Eigen::Map< const Eigen::Matrix< variableType, dim, dim, Eigen::RowMajor > > cmap( deformationGradient.data( ), dim, dim );
        detF = cmap.determinant( );

        inverseDeformationGradient = deformationGradient;
        Eigen::Map< Eigen::Matrix< variableType, dim, dim, Eigen::RowMajor > > map( inverseDeformationGradient.data( ), dim, dim );
        map = map.inverse( ).eval( );

        inverseMicroDeformation = microDeformation;
        new (&map) Eigen::Map< Eigen::Matrix< variableType, dim, dim, Eigen::RowMajor > >( inverseMicroDeformation.data( ), dim, dim );
        map = map.inverse( ).eval( );

        referenceHigherOrderStress = variableVector( tot_dim, 0 );

        variableVector temp_tot1( tot_dim, 0 );

        for ( unsigned int I = 0; I < dim; I++ ){
            for ( unsigned int i = 0; i < dim; i++ ){
                for ( unsigned int J = 0; J < dim; J++ ){
                    for ( unsigned int K = 0; K < dim; K++ ){
                        referenceHigherOrderStress[ dim * dim * I + dim * J + K ]
                            += inverseDeformationGradient[ dim * I + i ]
                             * higherOrderStress[ dim * dim * i + dim * J + K ];
                    }
                }
            }
        }
        referenceHigherOrderStress *= detF;

        for ( unsigned int I = 0; I < dim; I++ ){
            for ( unsigned int J = 0; J < dim; J++ ){
                for ( unsigned int i = 0; i < dim; i++ ){
                    for ( unsigned int K = 0; K < dim; K++ ){
                        temp_tot1[ dim * dim * I + dim * J + K ]
                            += inverseDeformationGradient[ dim * J + i ]
                             * referenceHigherOrderStress[ dim * dim * I + dim * i + K ];
                    }
                }
            }
        }

        std::fill( referenceHigherOrderStress.begin( ),
                   referenceHigherOrderStress.end( ),
                   0. );

        for ( unsigned int I = 0; I < dim; I++ ){
            for ( unsigned int J = 0; J < dim; J++ ){
                for ( unsigned int K = 0; K < dim; K++ ){
                    for ( unsigned int i = 0; i < dim; i++ ){
                        referenceHigherOrderStress[ dim * dim * I + dim * J + K ]
                            += inverseMicroDeformation[ dim * K + i ]
                             * temp_tot1[ dim * dim * I + dim * J + i ];
                    }
                }
            }
        }

        return;
    }

    void pullBackHigherOrderStress( const variableVector &higherOrderStress,
                                        const variableVector &deformationGradient,
                                        const variableVector &microDeformation,
                                        variableVector &referenceHigherOrderStress,
                                        variableMatrix &dReferenceHigherOrderStressdHigherOrderStress,
                                        variableMatrix &dReferenceHigherOrderStressdDeformationGradient,
                                        variableMatrix &dReferenceHigherOrderStressdMicroDeformation ){
        /*!
         * Compute the pull back operation on the higher order stress.
         *
         * \f$ M_{IJK} = J F_{Ii}^{-1} F_{Jj}^{-1} \chi_{Kk}^{-1} m_{ijk} \f$
         *
         * Also returns the Jacobians
         *
         * \f$ \frac{ \partial M_{IJK} }{ \partial m_{lmn} } = J F_{Il}^{-1} F_{Jm}^{-1} \chi_{Kn}^{-1} \f$
         * 
         * \f$ \frac{ \partial M_{IJK} }{ \partial F_{lL} } = F_{Ll}^{-1} M_{IJK} - F_{Il}^{-1} M_{LJK} - F_{Jl}^{-1} M_{ILK} \f$
         * 
         * \f$ \frac{ \partial M_{IJK} }{ \partial \chi_{lL} } = -\chi_{Kl}^{-1} M_{IJL} \f$
         *
         * \param &referenceHigherOrderStress: The higher order stress in the 
         *     reference configuration.
         * \param &deformationGradient: The deformation gradient which maps 
         *     between the reference and current configurations.
         * \param &microDeformation: The micro-deformation tensor.
         * \param &higherOrderStress: The higher order stress in the current configuration.
         * \param &dReferenceHigherOrderStressdHigherOrderStress: The derivative of the reference higher order stress w.r.t. the higher order stress
         * \param &dReferenceHigherOrderStressdDeformationGradient: The derivative of the reference higher order stress w.r.t. the deformation gradient
         * \param &dReferenceHigherOrderStressdMicroDeformation: The derivative of the reference higher order stress w.r.t. the micro deformation
         */

        variableVector _dReferenceHigherOrderStressdHigherOrderStress;
        variableVector _dReferenceHigherOrderStressdDeformationGradient;
        variableVector _dReferenceHigherOrderStressdMicroDeformation;

        TARDIGRADE_ERROR_TOOLS_CATCH( pullBackHigherOrderStress( higherOrderStress, deformationGradient, microDeformation,
                                                                 referenceHigherOrderStress, _dReferenceHigherOrderStressdHigherOrderStress,
                                                                 _dReferenceHigherOrderStressdDeformationGradient,
                                                                 _dReferenceHigherOrderStressdMicroDeformation ) );

        dReferenceHigherOrderStressdHigherOrderStress   = tardigradeVectorTools::inflate( _dReferenceHigherOrderStressdHigherOrderStress  , 27, 27 );
        dReferenceHigherOrderStressdDeformationGradient = tardigradeVectorTools::inflate( _dReferenceHigherOrderStressdDeformationGradient, 27,  9 );
        dReferenceHigherOrderStressdMicroDeformation    = tardigradeVectorTools::inflate( _dReferenceHigherOrderStressdMicroDeformation   , 27,  9 );

        return;

    }

    void pullBackHigherOrderStress( const variableVector &higherOrderStress,
                                        const variableVector &deformationGradient,
                                        const variableVector &microDeformation,
                                        variableVector &referenceHigherOrderStress,
                                        variableVector &dReferenceHigherOrderStressdHigherOrderStress,
                                        variableVector &dReferenceHigherOrderStressdDeformationGradient,
                                        variableVector &dReferenceHigherOrderStressdMicroDeformation ){
        /*!
         * Compute the pull back operation on the higher order stress.
         *
         * \f$ M_{IJK} = J F_{Ii}^{-1} F_{Jj}^{-1} \chi_{Kk}^{-1} m_{ijk} \f$
         *
         * Also returns the Jacobians
         *
         * \f$ \frac{ \partial M_{IJK} }{ \partial m_{lmn} } = J F_{Il}^{-1} F_{Jm}^{-1} \chi_{Kn}^{-1} \f$
         * 
         * \f$ \frac{ \partial M_{IJK} }{ \partial F_{lL} } = F_{Ll}^{-1} M_{IJK} - F_{Il}^{-1} M_{LJK} - F_{Jl}^{-1} M_{ILK} \f$
         * 
         * \f$ \frac{ \partial M_{IJK} }{ \partial \chi_{lL} } = -\chi_{Kl}^{-1} M_{IJL} \f$
         *
         * \param &referenceHigherOrderStress: The higher order stress in the 
         *     reference configuration.
         * \param &deformationGradient: The deformation gradient which maps 
         *     between the reference and current configurations.
         * \param &microDeformation: The micro-deformation tensor.
         * \param &higherOrderStress: The higher order stress in the current configuration.
         * \param &dReferenceHigherOrderStressdHigherOrderStress: The derivative of the reference higher order stress w.r.t. the higher order stress
         * \param &dReferenceHigherOrderStressdDeformationGradient: The derivative of the reference higher order stress w.r.t. the deformation gradient
         * \param &dReferenceHigherOrderStressdMicroDeformation: The derivative of the reference higher order stress w.r.t. the micro deformation
         */

        //Assume 3d
        constexpr unsigned int dim = 3;
        constexpr unsigned int sot_dim = dim * dim;
        constexpr unsigned int tot_dim = sot_dim * dim;

        variableType detF;
        variableVector inverseDeformationGradient, inverseMicroDeformation;

        TARDIGRADE_ERROR_TOOLS_CATCH( pullBackHigherOrderStress( higherOrderStress, deformationGradient, microDeformation,
                                                                 detF, inverseDeformationGradient, inverseMicroDeformation,
                                                                 referenceHigherOrderStress ) );

        dReferenceHigherOrderStressdHigherOrderStress   = variableVector( tot_dim * tot_dim, 0 );
        dReferenceHigherOrderStressdDeformationGradient = variableVector( tot_dim * sot_dim, 0 );
        dReferenceHigherOrderStressdMicroDeformation    = variableVector( tot_dim * sot_dim, 0 );

        for ( unsigned int I = 0; I < dim; I++){
            for ( unsigned int J = 0; J < dim; J++){
                for ( unsigned int K = 0; K < dim; K++){
                    for ( unsigned int l = 0; l < dim; l++){
                        for ( unsigned int m = 0; m < dim; m++){
                            dReferenceHigherOrderStressdDeformationGradient[ dim * dim * sot_dim * I + dim * sot_dim * J + sot_dim * K + dim * l + m ]
                                += inverseDeformationGradient[ dim * m + l ] * referenceHigherOrderStress[ dim * dim * I + dim * J + K ]
                                 - inverseDeformationGradient[ dim * I + l ] * referenceHigherOrderStress[ dim * dim * m + dim * J + K ]
                                 - inverseDeformationGradient[ dim * J + l ] * referenceHigherOrderStress[ dim * dim * I + dim * m + K ];

                            dReferenceHigherOrderStressdMicroDeformation[ dim * dim * sot_dim * I + dim * sot_dim * J + sot_dim * K + dim * l + m ]
                                -= inverseMicroDeformation[ dim * K + l ] * referenceHigherOrderStress[ dim * dim * I + dim * J + m ];

                            for ( unsigned int n = 0; n < dim; n++){

                                dReferenceHigherOrderStressdHigherOrderStress[ dim * dim * tot_dim * I + dim * tot_dim * J + tot_dim * K + dim * dim * l + dim * m + n ]
                                    += detF * inverseDeformationGradient[ dim * I + l ]
                                     * inverseDeformationGradient[ dim * J + m ]
                                     * inverseMicroDeformation[ dim * K + n ];
                            }
                        }
                    }
                }
            }
        }

        return;
    }

    void computeDeviatoricHigherOrderStress( const variableVector &higherOrderStress,
                                                 variableVector &deviatoricHigherOrderStress ){
        /*!
         * Compute the deviatoric part of the higher order stress.
         *
         * \f$ \text{dev} ( m_{ijk} ) = m_{ijk} - ( 1 / 3 ) m_{llk} \delta_{ij} \f$
         *
         * \param &higherOrderStress: The higher order stress in the current configuration.
         * \param &deviatoricHigherOrderStress: The deviatoric part of the higher order stress.
         */

        //Assume 3d
        constexpr unsigned int dim = 3;
        constexpr unsigned int sot_dim = dim * dim;
        constexpr unsigned int tot_dim = sot_dim * dim;

        TARDIGRADE_ERROR_TOOLS_CHECK( higherOrderStress.size() == tot_dim, "The higher order stress has an incorrect size" );

        deviatoricHigherOrderStress = higherOrderStress;

        for ( unsigned int i = 0; i < dim; i++ ){
            for ( unsigned int k = 0; k < dim; k++ ){
                for ( unsigned int j = 0; j < dim; j++ ){
                    deviatoricHigherOrderStress[ dim * dim * i + dim * i + j ] -= higherOrderStress[ dim * dim * k + dim * k + j ] / 3;
                }
            }
        }

        return;
    }

    void computeDeviatoricHigherOrderStress( const variableVector &higherOrderStress,
                                                 variableVector &deviatoricHigherOrderStress,
                                                 variableVector &dDeviatoricHigherOrderStressdHigherOrderStress){
        /*!
         * Compute the deviatoric part of the higher order stress.
         *
         * \f$ \text{dev} ( m_{ijk} ) = m_{ijk} - ( 1 / 3 ) m_{llk} \delta_{ij} \f$
         *
         * Also compute the Jacobian
         * \f$ \frac{ \partial \text{dev} ( m_{ijk} ) }{ \partial m_{mno} } = \delta_{im} \delta_{jn} \delta_{ko} - ( 1 / 3 ) \delta_{mn} \delta_{ko} \delta_{ij} \f$
         *
         * \param &higherOrderStress: The higher order stress in the current configuration.
         * \param &deviatoricHigherOrderStress: The deviatoric part of the higher order stress.
         * \param &dDeviatoricHigherOrderStressdHigherOrderStress: The gradient of the deviatoric part of the 
         *     higher order stress w.r.t. the higher order stress.
         */

        //Assume 3d
        constexpr unsigned int dim = 3;
        constexpr unsigned int sot_dim = dim * dim;
        constexpr unsigned int tot_dim = sot_dim * dim;

        TARDIGRADE_ERROR_TOOLS_CATCH( computeDeviatoricHigherOrderStress( higherOrderStress, deviatoricHigherOrderStress ) )

        dDeviatoricHigherOrderStressdHigherOrderStress = variableVector( tot_dim * tot_dim, 0 );

        for ( unsigned int i = 0; i < dim; i++ ){
            for ( unsigned int j = 0; j < dim; j++ ){
                for ( unsigned int k = 0; k < dim; k++ ){
                    dDeviatoricHigherOrderStressdHigherOrderStress[ tot_dim * dim * dim * i + tot_dim * dim * j + tot_dim * k + dim * dim * i + dim * j + k ]
                        += 1;
                    dDeviatoricHigherOrderStressdHigherOrderStress[ tot_dim * dim * dim * i + tot_dim * dim * i + tot_dim * j + dim * dim * k + dim * k + j ]
                        -= 1. / 3;
                }
            }
        }

        return;
    }

    void computeDeviatoricHigherOrderStress( const variableVector &higherOrderStress,
                                                 variableVector &deviatoricHigherOrderStress,
                                                 variableMatrix &dDeviatoricHigherOrderStressdHigherOrderStress){
        /*!
         * Compute the deviatoric part of the higher order stress.
         *
         * \f$ \text{dev} ( m_{ijk} ) = m_{ijk} - ( 1 / 3 ) m_{llk} \delta_{ij} \f$
         *
         * Also compute the Jacobian
         * \f$ \frac{ \partial \text{dev} ( m_{ijk} ) }{ \partial m_{mno} } = \delta_{im} \delta_{jn} \delta_{ko} - ( 1 / 3 ) \delta_{mn} \delta_{ko} \delta_{ij} \f$
         *
         * \param &higherOrderStress: The higher order stress in the current configuration.
         * \param &deviatoricHigherOrderStress: The deviatoric part of the higher order stress.
         * \param &dDeviatoricHigherOrderStressdHigherOrderStress: The gradient of the deviatoric part of the 
         *     higher order stress w.r.t. the higher order stress.
         */

        //Assume 3d
        constexpr unsigned int dim = 3;
        constexpr unsigned int sot_dim = dim * dim;
        constexpr unsigned int tot_dim = sot_dim * dim;

        variableVector _dDeviatoricHigherOrderStressdHigherOrderStress;

        TARDIGRADE_ERROR_TOOLS_CATCH( computeDeviatoricHigherOrderStress( higherOrderStress, deviatoricHigherOrderStress, _dDeviatoricHigherOrderStressdHigherOrderStress ) );

        dDeviatoricHigherOrderStressdHigherOrderStress = tardigradeVectorTools::inflate( _dDeviatoricHigherOrderStressdHigherOrderStress, tot_dim, tot_dim );

        return;
    }

    void computeReferenceHigherOrderStressPressure( const variableVector &referenceHigherOrderStress,
                                                        const variableVector &rightCauchyGreenDeformation,
                                                        variableVector &referenceHigherOrderPressure ){
        /*!
         * Compute the pressure for a higher-order stress in the reference configuration.
         * \f$p_K = \frac{1}{3} C_{AB} M_{ABK}\f$
         *
         * where \f$C_{AB}\f$ is the right Cauchy-Green deformation tensor and 
         * \f$M_{ABK}\f$ is the higher order stress tensor in the reference configuration.
         *
         * \param &referenceHigherOrderStress: The higher order stress in the 
         *     reference configuration.
         * \param &rightCauchyGreenDeformation: The right Cauchy-Green deformation
         *     tensor.
         * \param &referenceHigherOrderPressure: The higher order pressure.
         */

        //Assume 3D
        constexpr unsigned int dim = 3;
        constexpr unsigned int sot_dim = dim * dim;
        constexpr unsigned int tot_dim = sot_dim * dim;

        TARDIGRADE_ERROR_TOOLS_CHECK( rightCauchyGreenDeformation.size() == sot_dim,  "The right Cauchy-Green deformation tensor must have nine terms." );

        TARDIGRADE_ERROR_TOOLS_CHECK( referenceHigherOrderStress.size() == tot_dim,  "The higher order stress tensor must have 27 terms." );

        referenceHigherOrderPressure = variableVector( dim, 0 );

        for ( unsigned int AB = 0; AB < sot_dim; AB++ ){
            for ( unsigned int K = 0; K < dim; K++ ){
                    referenceHigherOrderPressure[K] += rightCauchyGreenDeformation[ AB ]
                                                     * referenceHigherOrderStress[ dim * AB + K ];
            }
        }

        referenceHigherOrderPressure /= 3;
        return;
    }

    void computeReferenceHigherOrderStressPressure( const variableVector &referenceHigherOrderStress,
                                                        const variableVector &rightCauchyGreenDeformation,
                                                        variableVector &referenceHigherOrderPressure,
                                                        variableMatrix &dpdM, variableMatrix &dpdC ){
        /*!
         * Compute the pressure for a higher-order stress in the reference configuration.
         * \f$p_K = \frac{1}{3} C_{AB} M_{ABK}\f$
         *
         * Also compute the Jacobians
         * \f$\frac{ \partial p_K }{ \partial M_{NOP} } = \frac{1}{3} C_{NO} \delta_{KP}\f$
         * 
         * \f$\frac{ \partial p_K }{ \partial C_{NO} } = \frac{1}{3} M_{NOK}\f$
         *
         * where \f$C_{AB}\f$ is the right Cauchy-Green deformation tensor and 
         * \f$M_{ABK}\f$ is the higher order stress tensor in the reference configuration.
         *
         * \param &referenceHigherOrderStress: The higher order stress in the 
         *     reference configuration.
         * \param &rightCauchyGreenDeformation: The right Cauchy-Green deformation
         *     tensor.
         * \param &referenceHigherOrderPressure: The higher order pressure.
         * \param &dpdM: The Jacobian of the pressure w.r.t. the higher order stress.
         * \param &dpdC: The Jacobian of the pressure w.r.t. the right Cauchy-Green 
         *     deformation tensor.
         */

        variableVector _dpdM;
        variableVector _dpdC;

        TARDIGRADE_ERROR_TOOLS_CATCH( computeReferenceHigherOrderStressPressure( referenceHigherOrderStress, rightCauchyGreenDeformation,
                                                                                 referenceHigherOrderPressure,
                                                                                 _dpdM, _dpdC ) );

        dpdM = tardigradeVectorTools::inflate( _dpdM, 3, 27 );
        dpdC = tardigradeVectorTools::inflate( _dpdC, 3,  9 );

        return;

    }

    void computeReferenceHigherOrderStressPressure( const variableVector &referenceHigherOrderStress,
                                                        const variableVector &rightCauchyGreenDeformation,
                                                        variableVector &referenceHigherOrderPressure,
                                                        variableVector &dpdM, variableVector &dpdC ){
        /*!
         * Compute the pressure for a higher-order stress in the reference configuration.
         * \f$p_K = \frac{1}{3} C_{AB} M_{ABK}\f$
         *
         * Also compute the Jacobians
         * \f$\frac{ \partial p_K }{ \partial M_{NOP} } = \frac{1}{3} C_{NO} \delta_{KP}\f$
         * 
         * \f$\frac{ \partial p_K }{ \partial C_{NO} } = \frac{1}{3} M_{NOK}\f$
         *
         * where \f$C_{AB}\f$ is the right Cauchy-Green deformation tensor and 
         * \f$M_{ABK}\f$ is the higher order stress tensor in the reference configuration.
         *
         * \param &referenceHigherOrderStress: The higher order stress in the 
         *     reference configuration.
         * \param &rightCauchyGreenDeformation: The right Cauchy-Green deformation
         *     tensor.
         * \param &referenceHigherOrderPressure: The higher order pressure.
         * \param &dpdM: The Jacobian of the pressure w.r.t. the higher order stress.
         * \param &dpdC: The Jacobian of the pressure w.r.t. the right Cauchy-Green 
         *     deformation tensor.
         */

        //Assume 3D
        constexpr unsigned int dim = 3;
        constexpr unsigned int sot_dim = dim * dim;
        constexpr unsigned int tot_dim = sot_dim * dim;

        TARDIGRADE_ERROR_TOOLS_CATCH( computeReferenceHigherOrderStressPressure( referenceHigherOrderStress, rightCauchyGreenDeformation,
                                                                                 referenceHigherOrderPressure ) );

        dpdM = variableVector( dim * tot_dim, 0 );
        dpdC = variableVector( dim * sot_dim, 0 );

        for ( unsigned int K = 0; K < dim; K++ ){
            for ( unsigned int N = 0; N < dim; N++ ){
                for ( unsigned int O = 0; O < dim; O++ ){
                    dpdC[ sot_dim * K + dim * N + O ] = referenceHigherOrderStress[ dim * dim * N + dim * O + K ];
                    dpdM[ tot_dim * K + dim * dim * N + dim * O + K ] = rightCauchyGreenDeformation[ dim * N + O ];
                }
            }
        }

        dpdM /= 3;
        dpdC /= 3;

        return;
    }

    void computeReferenceHigherOrderStressPressure( const variableVector &referenceHigherOrderStress,
                                                        const variableVector &rightCauchyGreenDeformation,
                                                        variableVector &referenceHigherOrderPressure,
                                                        variableMatrix &dpdM, variableMatrix &dpdC,
                                                        variableMatrix &d2pdMdC ){
        /*!
         * Compute the pressure for a higher-order stress in the reference configuration.
         * \f$p_K = \frac{1}{3} C_{AB} M_{ABK}\f$
         *
         * Also compute the Jacobians
         * \f$\frac{ \partial p_K }{ \partial M_{NOP} } = \frac{1}{3} C_{NO} \delta_{KP}\f$
         * 
         * \f$\frac{ \partial p_K }{ \partial C_{NO} } = \frac{1}{3} M_{NOK}\f$
         * 
         * \f$\frac{ \partial^2 p_K}{ \partial M_{NOP} C_{QR} } = \frac{1}{3} \delta_{NQ} \delta_{OR} \delta_{KP}\f$
         *
         * where \f$C_{AB}\f$ is the right Cauchy-Green deformation tensor and
         * \f$M_{ABK}\f$ is the higher order stress tensor in the reference configuration.
         *
         * \param &referenceHigherOrderStress: The higher order stress in the
         *     reference configuration.
         * \param &rightCauchyGreenDeformation: The right Cauchy-Green deformation
         *     tensor.
         * \param &referenceHigherOrderPressure: The higher order pressure.
         * \param &dpdM: The Jacobian of the pressure w.r.t. the higher order stress.
         * \param &dpdC: The Jacobian of the pressure w.r.t. the right Cauchy-Green
         *     deformation tensor.
         * \param &d2pdMdC: The second order jacobian of the pressure w.r.t the 
         *     reference higher order stress and right Cauchy-Green deformation tensor.
         */

        variableVector _dpdM;
        variableVector _dpdC;
        variableVector _d2pdMdC;

        TARDIGRADE_ERROR_TOOLS_CATCH( computeReferenceHigherOrderStressPressure( referenceHigherOrderStress, rightCauchyGreenDeformation,
                                                                                 referenceHigherOrderPressure,
                                                                                 _dpdM, _dpdC, _d2pdMdC ) );

        dpdM    = tardigradeVectorTools::inflate(    _dpdM, 3,  27 );
        dpdC    = tardigradeVectorTools::inflate(    _dpdC, 3,   9 );
        d2pdMdC = tardigradeVectorTools::inflate( _d2pdMdC, 3, 243 );

        return;
    }

    void computeReferenceHigherOrderStressPressure( const variableVector &referenceHigherOrderStress,
                                                        const variableVector &rightCauchyGreenDeformation,
                                                        variableVector &referenceHigherOrderPressure,
                                                        variableVector &dpdM, variableVector &dpdC,
                                                        variableVector &d2pdMdC ){
        /*!
         * Compute the pressure for a higher-order stress in the reference configuration.
         * \f$p_K = \frac{1}{3} C_{AB} M_{ABK}\f$
         *
         * Also compute the Jacobians
         * \f$\frac{ \partial p_K }{ \partial M_{NOP} } = \frac{1}{3} C_{NO} \delta_{KP}\f$
         * 
         * \f$\frac{ \partial p_K }{ \partial C_{NO} } = \frac{1}{3} M_{NOK}\f$
         * 
         * \f$\frac{ \partial^2 p_K}{ \partial M_{NOP} C_{QR} } = \frac{1}{3} \delta_{NQ} \delta_{OR} \delta_{KP}\f$
         *
         * where \f$C_{AB}\f$ is the right Cauchy-Green deformation tensor and
         * \f$M_{ABK}\f$ is the higher order stress tensor in the reference configuration.
         *
         * \param &referenceHigherOrderStress: The higher order stress in the
         *     reference configuration.
         * \param &rightCauchyGreenDeformation: The right Cauchy-Green deformation
         *     tensor.
         * \param &referenceHigherOrderPressure: The higher order pressure.
         * \param &dpdM: The Jacobian of the pressure w.r.t. the higher order stress.
         * \param &dpdC: The Jacobian of the pressure w.r.t. the right Cauchy-Green
         *     deformation tensor.
         * \param &d2pdMdC: The second order jacobian of the pressure w.r.t the 
         *     reference higher order stress and right Cauchy-Green deformation tensor.
         */

        //Assume 3D
        constexpr unsigned int dim = 3;
        constexpr unsigned int sot_dim = dim * dim;
        constexpr unsigned int tot_dim = sot_dim * dim;

        TARDIGRADE_ERROR_TOOLS_CATCH( computeReferenceHigherOrderStressPressure( referenceHigherOrderStress, rightCauchyGreenDeformation,
                                                                                 referenceHigherOrderPressure, dpdM, dpdC ) );

        d2pdMdC = variableVector( dim * tot_dim * sot_dim, 0 );
        for ( unsigned int K = 0; K < dim; K++ ){
            for ( unsigned int N = 0; N < dim; N++ ){
                for ( unsigned int O = 0; O < dim; O++ ){
                    d2pdMdC[ tot_dim * sot_dim * K + dim * dim * dim * dim * N + dim * dim * dim * O + dim * dim * K + dim * N + O ] += 1. / 3;
                }
            }
        }

        return;
    }

    void computeDeviatoricReferenceHigherOrderStress( const variableVector &referenceHigherOrderStress,
                                                          const variableVector &rightCauchyGreenDeformation,
                                                          variableVector &deviatoricReferenceHigherOrderStress ){
        /*!
         * Compute the deviatoric part of the higher order stress in the reference configuration.
         *
         * \f$ \text{dev} ( M_{IJK} ) = M_{IJK} - ( 1 / 3 ) (C^{-1})_{IJ} C_{AB} M_{ABK} \f$
         *
         * \param &referenceHigherOrderStress: The higher order stress in the reference configuration.
         * \param &rightCauchyGreenDeformation: The right Cauchy-Green deformation tensor.
         * \param &deviatoricReferenceHigherOrderStress: The deviatoric part of the higher order tensor in the 
         *     reference configuration.
         */

        //Compute the pressure
        variableVector pressure;
        TARDIGRADE_ERROR_TOOLS_CATCH( computeReferenceHigherOrderStressPressure( referenceHigherOrderStress, rightCauchyGreenDeformation, pressure ) );

        return computeDeviatoricReferenceHigherOrderStress( referenceHigherOrderStress, rightCauchyGreenDeformation,
                                                            pressure, deviatoricReferenceHigherOrderStress );
    }

    void computeDeviatoricReferenceHigherOrderStress( const variableVector &referenceHigherOrderStress,
                                                          const variableVector &rightCauchyGreenDeformation,
                                                          const variableVector &pressure,
                                                          variableVector &deviatoricReferenceHigherOrderStress ){
        /*!
         * Compute the deviatoric part of the higher order stress in the reference configuration.
         *
         * \f$ \text{dev} ( M_{IJK} ) = M_{IJK} - ( 1 / 3 ) (C^{-1})_{IJ} C_{AB} M_{ABK} \f$
         *
         * \param &referenceHigherOrderStress: The higher order stress in the reference configuration.
         * \param &rightCauchyGreenDeformation: The right Cauchy-Green deformation tensor.
         * \param &pressure: The pressure of the higher-order stress.
         * \param &deviatoricReferenceHigherOrderStress: The deviatoric part of the higher order tensor in the 
         *     reference configuration.
         */

        //Assume 3d
        constexpr unsigned int dim = 3;

        deviatoricReferenceHigherOrderStress = referenceHigherOrderStress;

        variableVector invRCG = rightCauchyGreenDeformation;
        Eigen::Map< Eigen::Matrix< variableType, dim, dim, Eigen::RowMajor > > map( invRCG.data( ), dim, dim );
        map = map.inverse( ).eval( );
        
        for ( unsigned int I = 0; I < dim; I++ ){
            for ( unsigned int J = 0; J < dim; J++ ){
                for ( unsigned int K = 0; K < dim; K++ ){
                    deviatoricReferenceHigherOrderStress[ dim * dim * I + dim * J + K ] -= invRCG[ dim * I + J ] * pressure[ K ];
                }
            }
        }

        return;
    }

    void computeDeviatoricReferenceHigherOrderStress( const variableVector &referenceHigherOrderStress,
                                                          const variableVector &rightCauchyGreenDeformation,
                                                          variableVector &deviatoricReferenceHigherOrderStress,
                                                          variableMatrix &dDeviatoricReferenceHigherOrderStressdReferenceHigherOrderStress,
                                                          variableMatrix &dDeviatoricReferenceHigherOrderStressdRCG ){
        /*!
         * Compute the deviatoric part of the higher order stress in the reference configuration.
         *
         * \f$ \text{dev} ( M_{IJK} ) = M_{IJK} - \frac{1}{3} (C^{-1})_{IJ} C_{AB} M_{ABK} = M_{IJK} - ( C^{-1} )_{IJ} p_{K} \f$
         *
         * Also compute Jacobians:
         * \f$ \frac{ \partial dev ( M_{IJK} ) }{ \partial M_{LMN} } = \delta_{IL} \delta_{JM} \delta_{KN} - (C^{-1})_{IJ} \frac{ \partial p_K }{ \partial M_{LMN} } \f$
         *
         * \f$ \frac{ \partial dev ( M_{IJK} ) }{ \partial C_{LM} } = \left( (C^{-1})_{IL} (C^{-1})_{MJ} p_{K} - (C^{-1})_{IJ} \frac{ \partial p_K }{ \partial C_{LM} } \right) \f$
         *
         * \param &referenceHigherOrderStress: The higher order stress in the reference configuration.
         * \param &rightCauchyGreenDeformation: The right Cauchy-Green deformation tensor.
         * \param &deviatoricReferenceHigherOrderStress: The deviatoric part of the higher order tensor in the 
         *     reference configuration.
         * \param &dDeviatoricReferenceHigherOrderStressdReferenceHigherOrderStress: The Jacobian of the 
         *     deviatoric part of the reference higher order stress w.r.t. the reference higher order stress.
         * \param &dDeviatoricReferenceHigherOrderStressdRCG: The Jacobian of the 
         *     deviatoric part of the reference higher order stress w.r.t. the right Cauchy-Green deformation tensor.
         */

        variableVector _dDeviatoricReferenceHigherOrderStressdReferenceHigherOrderStress;
        variableVector _dDeviatoricReferenceHigherOrderStressdRCG;

        TARDIGRADE_ERROR_TOOLS_CATCH( computeDeviatoricReferenceHigherOrderStress( referenceHigherOrderStress, rightCauchyGreenDeformation,
                                                                                   deviatoricReferenceHigherOrderStress,
                                                                                   _dDeviatoricReferenceHigherOrderStressdReferenceHigherOrderStress,
                                                                                   _dDeviatoricReferenceHigherOrderStressdRCG ) );

        dDeviatoricReferenceHigherOrderStressdReferenceHigherOrderStress = tardigradeVectorTools::inflate( _dDeviatoricReferenceHigherOrderStressdReferenceHigherOrderStress, 27, 27 );
        dDeviatoricReferenceHigherOrderStressdRCG                        = tardigradeVectorTools::inflate( _dDeviatoricReferenceHigherOrderStressdRCG                       , 27,  9 );

        return;

    }

    void computeDeviatoricReferenceHigherOrderStress( const variableVector &referenceHigherOrderStress,
                                                          const variableVector &rightCauchyGreenDeformation,
                                                          variableVector &deviatoricReferenceHigherOrderStress,
                                                          variableVector &dDeviatoricReferenceHigherOrderStressdReferenceHigherOrderStress,
                                                          variableVector &dDeviatoricReferenceHigherOrderStressdRCG ){
        /*!
         * Compute the deviatoric part of the higher order stress in the reference configuration.
         *
         * \f$ \text{dev} ( M_{IJK} ) = M_{IJK} - \frac{1}{3} (C^{-1})_{IJ} C_{AB} M_{ABK} = M_{IJK} - ( C^{-1} )_{IJ} p_{K} \f$
         *
         * Also compute Jacobians:
         * \f$ \frac{ \partial \text{ dev } ( M_{IJK} ) }{ \partial M_{LMN} } = \delta_{IL} \delta_{JM} \delta_{KN} - (C^{-1})_{IJ} \frac{ \partial p_K }{ \partial M_{LMN} } \f$
         * 
         * \f$ \frac{ \partial \text{ dev } ( M_{IJK} ) }{ \partial C_{LM} } = \left( (C^{-1})_{IL} (C^{-1})_{MJ} p_{K} - (C^{-1})_{IJ} \frac{ \partial p_K }{ \partial C_{LM} } \right) \f$
         *
         * \param &referenceHigherOrderStress: The higher order stress in the reference configuration.
         * \param &rightCauchyGreenDeformation: The right Cauchy-Green deformation tensor.
         * \param &deviatoricReferenceHigherOrderStress: The deviatoric part of the higher order tensor in the 
         *     reference configuration.
         * \param &dDeviatoricReferenceHigherOrderStressdReferenceHigherOrderStress: The Jacobian of the 
         *     deviatoric part of the reference higher order stress w.r.t. the reference higher order stress.
         * \param &dDeviatoricReferenceHigherOrderStressdRCG: The Jacobian of the 
         *     deviatoric part of the reference higher order stress w.r.t. the right Cauchy-Green deformation tensor.
         */

        variableVector invRCG, pressure;
        variableVector dpdM, dpdC;

        //Compute the pressure
        TARDIGRADE_ERROR_TOOLS_CATCH( computeReferenceHigherOrderStressPressure( referenceHigherOrderStress, rightCauchyGreenDeformation, pressure,
                                                                                 dpdM, dpdC ) );

        return computeDeviatoricReferenceHigherOrderStress( referenceHigherOrderStress, rightCauchyGreenDeformation,
                                                            pressure, dpdM, dpdC,
                                                            deviatoricReferenceHigherOrderStress,
                                                            dDeviatoricReferenceHigherOrderStressdReferenceHigherOrderStress,
                                                            dDeviatoricReferenceHigherOrderStressdRCG );

    }

    void computeDeviatoricReferenceHigherOrderStress( const variableVector &referenceHigherOrderStress,
                                                          const variableVector &rightCauchyGreenDeformation,
                                                          const variableVector &pressure,
                                                          const variableVector &dPressuredStress,
                                                          const variableVector &dPressuredRCG,
                                                          variableVector &deviatoricReferenceHigherOrderStress,
                                                          variableVector &dDeviatoricReferenceHigherOrderStressdReferenceHigherOrderStress,
                                                          variableVector &dDeviatoricReferenceHigherOrderStressdRCG ){
        /*!
         * Compute the deviatoric part of the higher order stress in the reference configuration.
         *
         * \f$ \text{dev} ( M_{IJK} ) = M_{IJK} - \frac{1}{3} (C^{-1})_{IJ} C_{AB} M_{ABK} = M_{IJK} - ( C^{-1} )_{IJ} p_{K} \f$
         *
         * Also compute Jacobians:
         * \f$\frac{ \partial \text{dev} ( M_{IJK} ) }{ \partial M_{LMN} } = \delta_{IL} \delta_{JM} \delta_{KN} - (C^{-1})_{IJ} \frac{ \partial p_K }{ \partial M_{LMN} }\f$
         * 
         * \f$\frac{ \partial \text{dev} ( M_{IJK} ) }{ \partial C_{LM} } = \left( (C^{-1})_{IL} (C^{-1})_{MJ} p_{K} - (C^{-1})_{IJ} \frac{ \partial p_K }{ \partial C_{LM} } \right)\f$
         *
         * \param &referenceHigherOrderStress: The higher order stress in the reference configuration.
         * \param &rightCauchyGreenDeformation: The right Cauchy-Green deformation tensor.
         * \param &pressure: The higher order pressure.
         * \param &dPressuredStress: The Jacobian of the pressure w.r.t. the stress.
         * \param &dPressuredRCG: The Jacobian of the pressure w.r.t. the right Cauchy Green 
         *     deformation tensor.
         * \param &deviatoricReferenceHigherOrderStress: The deviatoric part of the higher order tensor in the 
         *     reference configuration.
         * \param &dDeviatoricReferenceHigherOrderStressdReferenceHigherOrderStress: The Jacobian of the 
         *     deviatoric part of the reference higher order stress w.r.t. the reference higher order stress.
         * \param &dDeviatoricReferenceHigherOrderStressdRCG: The Jacobian of the 
         *     deviatoric part of the reference higher order stress w.r.t. the right Cauchy-Green deformation tensor.
         */

        //Assume 3d
        constexpr unsigned int dim = 3;
        constexpr unsigned int sot_dim = dim * dim;
        constexpr unsigned int tot_dim = sot_dim * dim;

        deviatoricReferenceHigherOrderStress = referenceHigherOrderStress;

        variableVector invRCG = rightCauchyGreenDeformation;
        Eigen::Map< Eigen::Matrix< variableType, dim, dim, Eigen::RowMajor > > map( invRCG.data( ), dim, dim );
        map = map.inverse( ).eval( );

        //Compute the deviatoric higher order stress        
        for ( unsigned int I = 0; I < dim; I++ ){
            for ( unsigned int J = 0; J < dim; J++ ){
                for ( unsigned int K = 0; K < dim; K++ ){
                    deviatoricReferenceHigherOrderStress[ dim * dim * I + dim * J + K ] -= invRCG[ dim * I + J ] * pressure[ K ];
                }
            }
        }

        //Compute the jacobians
        dDeviatoricReferenceHigherOrderStressdReferenceHigherOrderStress = variableVector( tot_dim * tot_dim, 0 );
        dDeviatoricReferenceHigherOrderStressdRCG = variableVector( tot_dim * sot_dim, 0 );

        for ( unsigned int I = 0; I < tot_dim; I++ ){
            dDeviatoricReferenceHigherOrderStressdReferenceHigherOrderStress[ tot_dim * I + I ] = 1.;
        }

        for ( unsigned int I = 0; I < dim; I++ ){
            for ( unsigned int J = 0; J < dim; J++ ){
                for ( unsigned int K = 0; K < dim; K++ ){
                    for ( unsigned int L = 0; L < dim; L++ ){
                        for ( unsigned int M = 0; M < dim; M++ ){
                            for ( unsigned int N = 0; N < dim; N++ ){
                                dDeviatoricReferenceHigherOrderStressdReferenceHigherOrderStress[ dim * dim * tot_dim * I + dim * tot_dim * J + tot_dim * K + dim * dim * L + dim * M + N ]
                                    -= invRCG[ dim * I + J ] * dPressuredStress[ tot_dim * K + dim * dim * L + dim * M + N ];
                            }

                            dDeviatoricReferenceHigherOrderStressdRCG[ dim * dim * sot_dim * I + dim * sot_dim * J + sot_dim * K + dim * L + M ]
                                = invRCG[ dim * I + L ] * invRCG[ dim * M + J ] * pressure[ K ] - invRCG[ dim * I + J ] * dPressuredRCG[ sot_dim * K + dim * L + M ];
                        }
                    }
                }
            }
        }

        return;
    }

    void computeDeviatoricReferenceHigherOrderStress( const variableVector &referenceHigherOrderStress,
                                                          const variableVector &rightCauchyGreenDeformation,
                                                          variableVector &deviatoricReferenceHigherOrderStress,
                                                          variableMatrix &dDeviatoricReferenceHigherOrderStressdReferenceHigherOrderStress,
                                                          variableMatrix &dDeviatoricReferenceHigherOrderStressdRCG,
                                                          variableMatrix &d2MdMdRCG ){
        /*!
         * Compute the deviatoric part of the higher order stress in the reference configuration.
         *
         * \f$ \text{dev} ( M_{IJK} ) = M_{IJK} - \frac{1}{3} (C^{-1})_{IJ} C_{AB} M_{ABK} = M_{IJK} - ( C^{-1} )_{IJ} p_{K} \f$
         *
         * Also compute Jacobians:
         * \f$ \frac{ \partial \text{dev} ( M_{IJK} ) }{ \partial M_{LMN} } = \delta_{IL} \delta_{JM} \delta_{KN} - (C^{-1})_{IJ} \frac{ \partial p_K }{ \partial M_{LMN} } \f$
         * 
         * \f$ \frac{ \partial \text{dev} ( M_{IJK} ) }{ \partial C_{LM} } = \left( (C^{-1})_{IL} (C^{-1})_{MJ} p_{K} - (C^{-1})_{IJ} \frac{ \partial p_K }{ \partial C_{LM} } \right) \f$
         *
         * \f$ \frac{ \partial \text{dev} ( M_{IJK} ) }{ \partial M_{LMN} \partial C_{OP} } = (C^{-1})_{IO} (C^{-1})_{PJ} \frac{ \partial p_{k }{ \partial M_{LMN} } - (C^{-1})_{IJ} \frac{ \partial^2 p_K}{ \partial M_{LMN} \partial C_{OP} } \f$
         *
         * \param &referenceHigherOrderStress: The higher order stress in the reference configuration.
         * \param &rightCauchyGreenDeformation: The right Cauchy-Green deformation tensor.
         * \param &deviatoricReferenceHigherOrderStress: The deviatoric part of the higher order tensor in the 
         *     reference configuration.
         * \param &dDeviatoricReferenceHigherOrderStressdReferenceHigherOrderStress: The Jacobian of the 
         *     deviatoric part of the reference higher order stress w.r.t. the reference higher order stress.
         * \param &dDeviatoricReferenceHigherOrderStressdRCG: The Jacobian of the 
         *     deviatoric part of the reference higher order stress w.r.t. the right Cauchy-Green deformation tensor.
         * \param &d2MdMdRCG: The mixed second derivative of the deviatoric part of the reference higher order 
         *     stress tensor with respect to the reference higher order stress tensor and the right Cauchy-Green deformation 
         *     tensor.
         */

        variableVector _dDeviatoricReferenceHigherOrderStressdReferenceHigherOrderStress;
        variableVector _dDeviatoricReferenceHigherOrderStressdRCG;
        variableVector _d2MdMdRCG;

        TARDIGRADE_ERROR_TOOLS_CATCH( computeDeviatoricReferenceHigherOrderStress( referenceHigherOrderStress, rightCauchyGreenDeformation,
                                                                                   deviatoricReferenceHigherOrderStress,
                                                                                   _dDeviatoricReferenceHigherOrderStressdReferenceHigherOrderStress,
                                                                                   _dDeviatoricReferenceHigherOrderStressdRCG,
                                                                                   _d2MdMdRCG ) );

        dDeviatoricReferenceHigherOrderStressdReferenceHigherOrderStress = tardigradeVectorTools::inflate( _dDeviatoricReferenceHigherOrderStressdReferenceHigherOrderStress, 27,  27 );
        dDeviatoricReferenceHigherOrderStressdRCG                        = tardigradeVectorTools::inflate( _dDeviatoricReferenceHigherOrderStressdRCG                       , 27,   9 );
        d2MdMdRCG                                                        = tardigradeVectorTools::inflate( _d2MdMdRCG                                                       , 27, 243 );

        return;
    }

    void computeDeviatoricReferenceHigherOrderStress( const variableVector &referenceHigherOrderStress,
                                                          const variableVector &rightCauchyGreenDeformation,
                                                          variableVector &deviatoricReferenceHigherOrderStress,
                                                          variableVector &dDeviatoricReferenceHigherOrderStressdReferenceHigherOrderStress,
                                                          variableVector &dDeviatoricReferenceHigherOrderStressdRCG,
                                                          variableVector &d2MdMdRCG ){
        /*!
         * Compute the deviatoric part of the higher order stress in the reference configuration.
         *
         * \f$ \text{dev} ( M_{IJK} ) = M_{IJK} - \frac{1}{3} (C^{-1})_{IJ} C_{AB} M_{ABK} = M_{IJK} - ( C^{-1} )_{IJ} p_{K} \f$
         *
         * Also compute Jacobians:
         * \f$ \frac{ \partial dev ( M_{IJK} ) }{ \partial M_{LMN} } = \delta_{IL} \delta_{JM} \delta_{KN} - (C^{-1})_{IJ} \frac{ \partial p_K }{ \partial M_{LMN} } \f$
         * 
         * \f$ \frac{ \partial dev ( M_{IJK} ) }{ \partial C_{LM} } = \left( (C^{-1})_{IL} (C^{-1})_{MJ} p_{K} - (C^{-1})_{IJ} \frac{ \partial p_K }{ \partial C_{LM} } \right) \f$
         *
         * \f$ \frac{ \partial dev( M_{IJK} ) }{ \partial M_{LMN} \partial C_{OP} } = (C^{-1})_{IO} (C^{-1})_{PJ} \frac{ \partial p_{k }{ \partial M_{LMN} } - (C^{-1})_{IJ} \frac{ \partial^2 p_K}{ \partial M_{LMN} \partial C_{OP} } \f$
         *
         * \param &referenceHigherOrderStress: The higher order stress in the reference configuration.
         * \param &rightCauchyGreenDeformation: The right Cauchy-Green deformation tensor.
         * \param &deviatoricReferenceHigherOrderStress: The deviatoric part of the higher order tensor in the 
         * \   reference configuration.
         * \param &dDeviatoricReferenceHigherOrderStressdReferenceHigherOrderStress: The Jacobian of the 
         * \   deviatoric part of the reference higher order stress w.r.t. the reference higher order stress.
         * \param &dDeviatoricReferenceHigherOrderStressdRCG: The Jacobian of the 
         * \   deviatoric part of the reference higher order stress w.r.t. the right Cauchy-Green deformation tensor.
         * \param &d2MdMdRCG: The mixed second derivative of the deviatoric part of the reference higher order 
         *     stress tensor with respect to the reference higher order stress tensor and the right Cauchy-Green deformation 
         *     tensor.
         */
        variableVector invRCG, pressure;
        variableVector dpdM, dpdC, d2pdMdRCG;

        //Compute the pressure
        TARDIGRADE_ERROR_TOOLS_CATCH( computeReferenceHigherOrderStressPressure( referenceHigherOrderStress, rightCauchyGreenDeformation, pressure,
                                                                                 dpdM, dpdC, d2pdMdRCG )  );

        return computeDeviatoricReferenceHigherOrderStress( referenceHigherOrderStress, rightCauchyGreenDeformation,
                                                            pressure, dpdM, dpdC, d2pdMdRCG,
                                                            deviatoricReferenceHigherOrderStress,
                                                            dDeviatoricReferenceHigherOrderStressdReferenceHigherOrderStress,
                                                            dDeviatoricReferenceHigherOrderStressdRCG,
                                                            d2MdMdRCG );
    }

    void computeDeviatoricReferenceHigherOrderStress( const variableVector &referenceHigherOrderStress,
                                                          const variableVector &rightCauchyGreenDeformation,
                                                          const variableVector &pressure,
                                                          const variableVector &dPressuredStress,
                                                          const variableVector &dPressuredRCG,
                                                          const variableVector &d2PressuredStressdRCG,
                                                          variableVector &deviatoricReferenceHigherOrderStress,
                                                          variableVector &dDeviatoricReferenceHigherOrderStressdReferenceHigherOrderStress,
                                                          variableVector &dDeviatoricReferenceHigherOrderStressdRCG,
                                                          variableVector &d2MdMdRCG ){
        /*!
         * Compute the deviatoric part of the higher order stress in the reference configuration.
         *
         * \f$ \text{dev} ( M_{IJK} ) = M_{IJK} - \frac{1}{3} (C^{-1})_{IJ} C_{AB} M_{ABK} = M_{IJK} - ( C^{-1} )_{IJ} p_{K} \f$
         *
         * Also compute Jacobians:
         * \f$ \frac{ \partial \text{ dev } ( M_{IJK} ) }{ \partial M_{LMN} } = \delta_{IL} \delta_{JM} \delta_{KN} - (C^{-1})_{IJ} \frac{ \partial p_K }{ \partial M_{LMN} } \f$
         * 
         * \f$ \frac{ \partial \text{ dev } ( M_{IJK} ) }{ \partial C_{LM} } = \left( (C^{-1})_{IL} (C^{-1})_{MJ} p_{K} - (C^{-1})_{IJ} \frac{ \partial p_K }{ \partial C_{LM} } \right) \f$
         *
         * \f$ \frac{ \partial \text{ dev } ( M_{IJK} ) }{ \partial M_{LMN} \partial C_{OP} } = (C^{-1})_{IO} (C^{-1})_{PJ} \frac{ \partial p_{k }{ \partial M_{LMN} } - (C^{-1})_{IJ} \frac{ \partial^2 p_K}{ \partial M_{LMN} \partial C_{OP} } \f$
         *
         * \param &referenceHigherOrderStress: The higher order stress in the reference configuration.
         * \param &rightCauchyGreenDeformation: The right Cauchy-Green deformation tensor.
         * \param &pressure: The higher order pressure.
         * \param &dPressuredStress: The Jacobian of the pressure w.r.t. the stress.
         * \param &dPressuredRCG: The Jacobian of the pressure w.r.t. the right Cauchy-Green 
         * \   deformation tensor.
         * \param &d2PressuredStressdRCG: The second order Jacobian of the pressure w.r.t. the 
         * \   stress and the right Cauchy-Green deformation tensor.
         * \param &deviatoricReferenceHigherOrderStress: The deviatoric part of the higher order tensor in the 
         * \   reference configuration.
         * \param &dDeviatoricReferenceHigherOrderStressdReferenceHigherOrderStress: The Jacobian of the 
         * \   deviatoric part of the reference higher order stress w.r.t. the reference higher order stress.
         * \param &dDeviatoricReferenceHigherOrderStressdRCG: The Jacobian of the 
         * \   deviatoric part of the reference higher order stress w.r.t. the right Cauchy-Green deformation tensor.
         * \param &d2MdMdRCG: The mixed second derivative of the deviatoric part of the reference higher order 
         *     stress tensor with respect to the reference higher order stress tensor and the right Cauchy-Green deformation 
         *     tensor.
         */
        
        //Assume 3d
        constexpr unsigned int dim = 3;
        constexpr unsigned int sot_dim = dim * dim;
        constexpr unsigned int tot_dim = sot_dim * dim;

        deviatoricReferenceHigherOrderStress = referenceHigherOrderStress;

        variableVector invRCG = rightCauchyGreenDeformation;
        Eigen::Map< Eigen::Matrix< variableType, dim, dim, Eigen::RowMajor > > map( invRCG.data( ), dim, dim );
        map = map.inverse( ).eval( );

        //Compute the deviatoric higher order stress        
        for ( unsigned int I = 0; I < dim; I++ ){
            for ( unsigned int J = 0; J < dim; J++ ){
                for ( unsigned int K = 0; K < dim; K++ ){
                    deviatoricReferenceHigherOrderStress[ dim * dim * I + dim * J + K ] -= invRCG[ dim * I + J ] * pressure[ K ];
                }
            }
        }

        //Compute the first order Jacobians
        dDeviatoricReferenceHigherOrderStressdReferenceHigherOrderStress = variableVector( tot_dim * tot_dim, 0 );
        dDeviatoricReferenceHigherOrderStressdRCG = variableVector( tot_dim * sot_dim, 0 );

        for ( unsigned int I = 0; I < tot_dim; I++ ){
            dDeviatoricReferenceHigherOrderStressdReferenceHigherOrderStress[ tot_dim * I + I ] = 1.;
        }

        for ( unsigned int I = 0; I < dim; I++ ){
            for ( unsigned int J = 0; J < dim; J++ ){
                for ( unsigned int K = 0; K < dim; K++ ){
                    for ( unsigned int L = 0; L < dim; L++ ){
                        for ( unsigned int M = 0; M < dim; M++ ){
                            for ( unsigned int N = 0; N < dim; N++ ){
                                dDeviatoricReferenceHigherOrderStressdReferenceHigherOrderStress[ dim * dim * tot_dim * I + dim * tot_dim * J + tot_dim * K + dim * dim * L + dim * M + N ]
                                    -= invRCG[ dim * I + J ] * dPressuredStress[ tot_dim * K + dim * dim * L + dim * M + N ];
                            }

                            dDeviatoricReferenceHigherOrderStressdRCG[ dim * dim * sot_dim * I + dim * sot_dim * J + sot_dim * K + dim * L + M ]
                                = invRCG[ dim * I + L ] * invRCG[ dim * M + J ] * pressure[ K ] - invRCG[ dim * I + J ] * dPressuredRCG[ sot_dim * K + dim * L + M ];
                        }
                    }
                }
            }
        }

        //Compute the second order Jacobian

        d2MdMdRCG = variableVector( tot_dim * tot_dim * sot_dim, 0 );

        for ( unsigned int I = 0; I < dim; I++ ){
            for ( unsigned int J = 0; J < dim; J++ ){
                for ( unsigned int K = 0; K < dim; K++ ){
                    for ( unsigned int L = 0; L < dim; L++ ){
                        for ( unsigned int M = 0; M < dim; M++ ){
                            for ( unsigned int N = 0; N < dim; N++ ){
                                for ( unsigned int O = 0; O < dim; O++ ){
                                    for ( unsigned int P = 0; P < dim; P++ ){
                                        d2MdMdRCG[ dim * dim * tot_dim * sot_dim * I + dim * tot_dim * sot_dim * J + tot_dim * sot_dim * K + dim * dim * dim * dim * L + dim * dim * dim * M + dim * dim * N + dim * O + P ]
                                            = invRCG[ dim * I + O ] * invRCG[ dim * P + J ] * dPressuredStress[ tot_dim * K + dim * dim * L + dim * M + N ]
                                            - invRCG[ dim * I + J ] * d2PressuredStressdRCG[ tot_dim * sot_dim * K + dim * dim * dim * dim * L + dim * dim * dim * M + dim * dim * N + dim * O + P ];
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        return;
    }

    void computeDeviatoricSecondOrderStress( const variableVector &secondOrderStress,
                                                 variableVector &deviatoricSecondOrderStress ){
        /*!
         * Compute the deviatoric part of a second order stress measure in the current configuration.
         * \f$ \hat{s}_{ij} = s_{ij} - \frac{1}{3} s_{mm} \delta_{ij} \f$
         *
         * \param &secondOrderStress: The stress measure in the current configuration.
         * \param &deviatoricSecondOrderStress: The deviatoric part of the second order stress 
         *     measure.
         */

        //Assume 3d
        constexpr unsigned int dim = 3;

        variableType trace = tardigradeVectorTools::trace( secondOrderStress );

        deviatoricSecondOrderStress = secondOrderStress;
        for ( unsigned int i = 0; i < dim; i++ ){ deviatoricSecondOrderStress[ dim * i + i ] -= trace / 3.; }

        return;
    }

    void computeDeviatoricSecondOrderStress( const variableVector &secondOrderStress,
                                                 variableVector &deviatoricSecondOrderStress,
                                                 variableVector &dDeviatoricStressdStress ){
        /*!
         * Compute the deviatoric part of a second order stress measure in the current configuration.
         * \f$ \hat{s}_{ij} = s_{ij} - \frac{1}{3} s_{mm} \delta_{ij} \f$
         *
         * Also return the jacobian
         * \f$ \frac{\partial \hat{s}_{ij}}{\partial s_{kl} } \delta_{ik} \delta_{jl} - \frac{1}{3} \delta_{ij} \delta_{kl} \f$
         *
         * \param &secondOrderStress: The stress measure in the current configuration.
         * \param &deviatoricSecondOrderStress: The deviatoric part of the second order stress 
         *     measure.
         * \param &dDeviatoricStressdStress: The jacobian of the deviatoric stress.
         */

        //Assume 3d
        constexpr unsigned int dim = 3;
        constexpr unsigned int sot_dim = dim * dim;

        TARDIGRADE_ERROR_TOOLS_CATCH( computeDeviatoricSecondOrderStress( secondOrderStress, deviatoricSecondOrderStress ) );

        dDeviatoricStressdStress = variableVector( sot_dim * sot_dim, 0 );
        for ( unsigned int i = 0; i < sot_dim; i++ ){
            dDeviatoricStressdStress[ sot_dim * i + i ] = 1;
        }

        for ( unsigned int i = 0; i < dim; i++ ){
            for ( unsigned int j = 0; j < dim; j++ ){
                dDeviatoricStressdStress[ sot_dim * dim * i + sot_dim * i + dim * j + j ] -= 1. / 3;
            }
        }

        return;
    }

    void computeDeviatoricSecondOrderStress( const variableVector &secondOrderStress,
                                                 variableVector &deviatoricSecondOrderStress,
                                                 variableMatrix &dDeviatoricStressdStress ){
        /*!
         * Compute the deviatoric part of a second order stress measure in the current configuration.
         * \f$ \hat{s}_{ij} = s_{ij} - \frac{1}{3} s_{mm} \delta_{ij} \f$
         *
         * Also return the jacobian
         * \f$ \frac{\partial \hat{s}_{ij}}{\partial s_{kl} } \delta_{ik} \delta_{jl} - \frac{1}{3} \delta_{ij} \delta_{kl} \f$
         *
         * \param &secondOrderStress: The stress measure in the current configuration.
         * \param &deviatoricSecondOrderStress: The deviatoric part of the second order stress 
         *     measure.
         * \param &dDeviatoricStressdStress: The jacobian of the deviatoric stress.
         */

        //Assume 3d
        constexpr unsigned int dim = 3;
        constexpr unsigned int sot_dim = dim * dim;

        variableVector _dDeviatoricStressdStress;

        TARDIGRADE_ERROR_TOOLS_CATCH( computeDeviatoricSecondOrderStress( secondOrderStress, deviatoricSecondOrderStress, _dDeviatoricStressdStress ) );

        dDeviatoricStressdStress = tardigradeVectorTools::inflate( _dDeviatoricStressdStress, sot_dim, sot_dim );

        return;
    }

    void computeReferenceSecondOrderStressPressure( const variableVector &referenceStressMeasure,
                                                        const variableVector &rightCauchyGreen,  variableType &pressure ){
        /*!
         * Compute the pressure part of a second order stress measure in the reference configuration.
         * \f$ p = \frac{1}{3} C_{IJ} S_{IJ} \f$
         *
         * where \f$ C_{IJ} \f$ is the right Cauchy-Green deformation tensor and \f$ S_{IJ} \f$ is the stress measure.
         *
         * \param &referenceStressMeasure: The stress measure in the reference configuration.
         * \param &rightCauchyGreen: The right Cauchy-Green deformation tensor between the 
         *     current configuration and the reference configuration of the stress tensor.
         * \param &pressure: The computed pressure.
         */

        TARDIGRADE_ERROR_TOOLS_CHECK( referenceStressMeasure.size() == rightCauchyGreen.size(), "The stress measure and right Cauchy-Green deformation tensors aren't the same size" );

        pressure = tardigradeVectorTools::dot( referenceStressMeasure, rightCauchyGreen ) / 3;

        return;
    }

    void computeReferenceSecondOrderStressPressure( const variableVector &referenceStressMeasure,
                                                        const variableVector &rightCauchyGreen,  variableType &pressure,
                                                        variableVector &dpdStress, variableVector &dpdRCG ){
        /*!
         * Compute the pressure part of a second order stress measure in the reference configuration.
         * \f$ p = \frac{1}{3} C_{IJ} S_{IJ} \f$
         *
         * where \f$C_{IJ}\f$ is the right Cauchy-Green deformation tensor and \f$ S_{IJ} \f$ is the stress measure.
         *
         * Also compute the Jacobians
         * \f$ \frac{ \partial p }{ \partial C_{IJ} } = \frac{1}{3} S_{IJ} \f$
         * 
         * \f$ \frac{ \partial p }{ \partial S_{IJ} } = \frac{1}{3} C_{IJ} \f$
         *
         * \param &referenceStressMeasure: The stress measure in the reference configuration.
         * \param &rightCauchyGreen: The right Cauchy-Green deformation tensor between the 
         *     current configuration and the reference configuration of the stress tensor.
         * \param &pressure: The computed pressure.
         * \param &dpdStress: The Jacobian of the pressure w.r.t. the stress.
         * \param &dpdRCG: The Jacobian of the pressure w.r.t. the right Cauchy-Green deformation tensor.
         */

        TARDIGRADE_ERROR_TOOLS_CATCH( computeReferenceSecondOrderStressPressure( referenceStressMeasure, rightCauchyGreen, pressure ) );

        dpdStress = rightCauchyGreen / 3;
        dpdRCG    = referenceStressMeasure / 3;

        return;
    }

    void computeReferenceSecondOrderStressPressure( const variableVector &referenceStressMeasure,
                                                        const variableVector &rightCauchyGreen,  variableType &pressure,
                                                        variableVector &dpdStress, variableVector &dpdRCG,
                                                        variableMatrix &d2pdStressdRCG ){
        /*!
         * Compute the pressure part of a second order stress measure in the reference configuration.
         * \f$ p = \frac{1}{3} C_{IJ} S_{IJ} \f$
         *
         * where C_{IJ} is the right Cauchy-Green deformation tensor and S_{IJ} is the stress measure.
         *
         * Also compute the Jacobians
         * \f$ \frac{ \partial p }{ \partial C_{IJ} } = \frac{1}{3} S_{IJ} \f$
         * 
         * \f$ \frac{ \partial p }{ \partial S_{IJ} } = \frac{1}{3} C_{IJ} \f$
         *
         * \f$ \frac{ \partial^2 p}{ \partial S_{IJ} \partial C_{KL} } = \frac{1}{3} \delta_{IK} \delta_{JL} \f$
         *
         * \param &referenceStressMeasure: The stress measure in the reference configuration.
         * \param &rightCauchyGreen: The right Cauchy-Green deformation tensor between the 
         *     current configuration and the reference configuration of the stress tensor.
         * \param &pressure: The computed pressure.
         * \param &dpdStress: The Jacobian of the pressure w.r.t. the stress.
         * \param &dpdRCG: The Jacobian of the pressure w.r.t. the right Cauchy-Green deformation tensor.
         * \param &d2pdStressdRCG: The second-order Jacobian of the pressure w.r.t. the stress and the 
         *     right Cauchy-Green deformation tensor.
         */

        variableVector _d2pdStressdRCG;

        TARDIGRADE_ERROR_TOOLS_CATCH( computeReferenceSecondOrderStressPressure( referenceStressMeasure, rightCauchyGreen, pressure,
                                                                                 dpdStress, dpdRCG, _d2pdStressdRCG ) );

        d2pdStressdRCG = tardigradeVectorTools::inflate( _d2pdStressdRCG, 9, 9 );

        return;

    }

    void computeReferenceSecondOrderStressPressure( const variableVector &referenceStressMeasure,
                                                        const variableVector &rightCauchyGreen,  variableType &pressure,
                                                        variableVector &dpdStress, variableVector &dpdRCG,
                                                        variableVector &d2pdStressdRCG ){
        /*!
         * Compute the pressure part of a second order stress measure in the reference configuration.
         * \f$ p = \frac{1}{3} C_{IJ} S_{IJ} \f$
         *
         * where \f$C_{IJ}\f$ is the right Cauchy-Green deformation tensor and \f$S_{IJ}\f$ is the stress measure.
         *
         * Also compute the Jacobians
         * \f$\frac{ \partial p }{ \partial C_{IJ} } = \frac{1}{3} S_{IJ}\f$
         * \f$\frac{ \partial p }{ \partial S_{IJ} } = \frac{1}{3} C_{IJ}\f$
         *
         * \f$\frac{ \partial^2 p}{ \partial S_{IJ} \partial C_{KL} } = \frac{1}{3} \delta_{IK} \delta_{JL}\f$
         *
         * \param &referenceStressMeasure: The stress measure in the reference configuration.
         * \param &rightCauchyGreen: The right Cauchy-Green deformation tensor between the 
         *     current configuration and the reference configuration of the stress tensor.
         * \param &pressure: The computed pressure.
         * \param &dpdStress: The Jacobian of the pressure w.r.t. the stress.
         * \param &dpdRCG: The Jacobian of the pressure w.r.t. the right Cauchy-Green deformation tensor.
         * \param &d2pdStressdRCG: The second-order Jacobian of the pressure w.r.t. the stress and the 
         *     right Cauchy-Green deformation tensor.
         */

        constexpr unsigned int dim = 3;
        constexpr unsigned int sot_dim = dim * dim;

        TARDIGRADE_ERROR_TOOLS_CATCH( computeReferenceSecondOrderStressPressure( referenceStressMeasure, rightCauchyGreen, pressure,
                                                                                 dpdStress, dpdRCG ) );

        d2pdStressdRCG = variableVector( sot_dim * sot_dim, 0 );
        for ( unsigned int i = 0; i < sot_dim; i++ ){ d2pdStressdRCG[ sot_dim * i + i ] = 1./3; }

        return;
    }

    void computeDeviatoricReferenceSecondOrderStress( const variableVector &secondOrderReferenceStress,
                                                          const variableVector &rightCauchyGreenDeformation,
                                                          variableVector &deviatoricSecondOrderReferenceStress ){
        /*!
         * Compute the deviatoric part of a second order stress measure in the reference configuration.
         * \f$\hat{S}_{IJ} = S_{IJ} - \frac{1}{3} C_{AB} S_{AB} (C^{-1})_{IJ}\f$
         *
         * \param &secondOrderReferenceStress: The stress measure in the reference configuration.
         * \param &rightCauchyGreenDeformation: The right Cauchy Green Deformation tensor of the 
         *     deformation between configurations.
         * \param &deviatoricSecondOrderReferenceStress: The deviatoric part of the second order 
         *     stress in the reference configuration.
         */

        variableType pressure;
        TARDIGRADE_ERROR_TOOLS_CATCH( computeReferenceSecondOrderStressPressure( secondOrderReferenceStress, rightCauchyGreenDeformation, pressure ) );

        return computeDeviatoricReferenceSecondOrderStress( secondOrderReferenceStress, rightCauchyGreenDeformation,
                                                            pressure, deviatoricSecondOrderReferenceStress );
    }

    void computeDeviatoricReferenceSecondOrderStress( const variableVector &secondOrderReferenceStress,
                                                          const variableVector &rightCauchyGreenDeformation,
                                                          const variableType &pressure,
                                                          variableVector &deviatoricSecondOrderReferenceStress ){
        /*!
         * Compute the deviatoric part of a second order stress measure in the reference configuration.
         * \f$\hat{S}_{IJ} = S_{IJ} - \frac{1}{3} C_{AB} S_{AB} (C^{-1})_{IJ}\f$
         *
         * \param &secondOrderReferenceStress: The stress measure in the reference configuration.
         * \param &rightCauchyGreenDeformation: The right Cauchy Green Deformation tensor of the 
         *     deformation between configurations.
         * \param &pressure: The pressure computed from the reference stress.
         * \param &deviatoricSecondOrderReferenceStress: The deviatoric part of the second order 
         *     stress in the reference configuration.
         */

        //Assume 3d
        constexpr unsigned int dim = 3;

        variableVector invRCG = rightCauchyGreenDeformation;
        Eigen::Map< Eigen::Matrix< variableType, dim, dim, Eigen::RowMajor > > map( invRCG.data( ), dim, dim );
        map = map.inverse( ).eval( );

        deviatoricSecondOrderReferenceStress = secondOrderReferenceStress - pressure * invRCG;

        return;
    }

    void computeDeviatoricReferenceSecondOrderStress( const variableVector &secondOrderReferenceStress,
                                                          const variableVector &rightCauchyGreenDeformation,
                                                          variableVector &deviatoricSecondOrderReferenceStress,
                                                          variableMatrix &dDeviatoricReferenceStressdReferenceStress,
                                                          variableMatrix &dDeviatoricReferenceStressdRCG ){
        /*!
         * Compute the deviatoric part of a second order stress measure in the reference configuration.
         * \f$\hat{S}_{IJ} = S_{IJ} - \frac{1}{3} C_{AB} S_{AB} (C^{-1})_{IJ}\f$
         * 
         * Also compute the Jacobians.
         * \f$\frac{ \partial \hat{S}_{IJ} }{ \partial S_{KL} } = \delta_{IK} \delta_{LJ} - \frac{1}{3} C_{KL} (C^{-1})_{IJ}\f$
         * 
         * \f$\frac{ \partial \hat{S}_{IJ} }{ \partial C_{KL} } = \frac{1}{3} \left( C_{AB} S_{AB} (C^{-1})_{IK} (C^{-1})_{LJ} - S_{KL} (C^{-1}_{IJ}) \right)\f$
         *
         * \param &secondOrderReferenceStress: The stress measure in the reference configuration.
         * \param &rightCauchyGreenDeformation: The right Cauchy Green Deformation tensor of the 
         *     deformation between configurations.
         * \param &deviatoricSecondOrderReferenceStress: The deviatoric part of the second order 
         *     stress in the reference configuration.
         * \param &dDeviatoricReferenceStressdReferenceStress: The jacobian w.r.t. the reference stress.
         * \param &dDeviatoricReferenceStressdRCG: The jacobian w.r.t. the right Cauchy Green deformation
         *     tensor.
         */

        variableVector _dDeviatoricReferenceStressdReferenceStress;
        variableVector _dDeviatoricReferenceStressdRCG;

        TARDIGRADE_ERROR_TOOLS_CATCH( computeDeviatoricReferenceSecondOrderStress( secondOrderReferenceStress,
                                                                                   rightCauchyGreenDeformation,
                                                                                   deviatoricSecondOrderReferenceStress,
                                                                                   _dDeviatoricReferenceStressdReferenceStress,
                                                                                   _dDeviatoricReferenceStressdRCG ) );

        dDeviatoricReferenceStressdReferenceStress = tardigradeVectorTools::inflate( _dDeviatoricReferenceStressdReferenceStress, 9, 9 );
        dDeviatoricReferenceStressdRCG             = tardigradeVectorTools::inflate( _dDeviatoricReferenceStressdRCG            , 9, 9 );

        return;

    }

    void computeDeviatoricReferenceSecondOrderStress( const variableVector &secondOrderReferenceStress,
                                                          const variableVector &rightCauchyGreenDeformation,
                                                          variableVector &deviatoricSecondOrderReferenceStress,
                                                          variableVector &dDeviatoricReferenceStressdReferenceStress,
                                                          variableVector &dDeviatoricReferenceStressdRCG ){
        /*!
         * Compute the deviatoric part of a second order stress measure in the reference configuration.
         * \f$\hat{S}_{IJ} = S_{IJ} - \frac{1}{3} C_{AB} S_{AB} (C^{-1})_{IJ}\f$
         * 
         * Also compute the Jacobians.
         * \f$\frac{ \partial \hat{S}_{IJ} }{ \partial S_{KL} } = \delta_{IK} \delta_{LJ} - \frac{1}{3} C_{KL} (C^{-1})_{IJ}\f$
         * 
         * \f$\frac{ \partial \hat{S}_{IJ} }{ \partial C_{KL} } = \frac{1}{3} \left( C_{AB} S_{AB} (C^{-1})_{IK} (C^{-1})_{LJ} - S_{KL} (C^{-1}_{IJ}) \right)\f$
         *
         * \param &secondOrderReferenceStress: The stress measure in the reference configuration.
         * \param &rightCauchyGreenDeformation: The right Cauchy Green Deformation tensor of the 
         *     deformation between configurations.
         * \param &deviatoricSecondOrderReferenceStress: The deviatoric part of the second order 
         *     stress in the reference configuration.
         * \param &dDeviatoricReferenceStressdReferenceStress: The jacobian w.r.t. the reference stress.
         * \param &dDeviatoricReferenceStressdRCG: The jacobian w.r.t. the right Cauchy Green deformation
         *     tensor.
         */

        variableType pressure;
        variableVector dpdS, dpdC;
        TARDIGRADE_ERROR_TOOLS_CATCH( computeReferenceSecondOrderStressPressure( secondOrderReferenceStress, rightCauchyGreenDeformation, pressure,
                                                                                 dpdS, dpdC ) );

        return computeDeviatoricReferenceSecondOrderStress( secondOrderReferenceStress, rightCauchyGreenDeformation,
                                                            pressure, dpdS, dpdC,
                                                            deviatoricSecondOrderReferenceStress,
                                                            dDeviatoricReferenceStressdReferenceStress,
                                                            dDeviatoricReferenceStressdRCG );
    }

    void computeDeviatoricReferenceSecondOrderStress( const variableVector &secondOrderReferenceStress,
                                                          const variableVector &rightCauchyGreenDeformation,
                                                          const variableType &pressure,
                                                          const variableVector &dPressuredStress,
                                                          const variableVector &dPressuredRCG,
                                                          variableVector &deviatoricSecondOrderReferenceStress,
                                                          variableVector &dDeviatoricReferenceStressdReferenceStress,
                                                          variableVector &dDeviatoricReferenceStressdRCG ){
        /*!
         * Compute the deviatoric part of a second order stress measure in the reference configuration.
         * \f$\hat{S}_{IJ} = S_{IJ} - \frac{1}{3} C_{AB} S_{AB} (C^{-1})_{IJ}\f$
         * 
         * Also compute the Jacobians.
         * \f$\frac{ \partial \hat{S}_{IJ} }{ \partial S_{KL} } = \delta_{IK} \delta_{LJ} - \frac{1}{3} C_{KL} (C^{-1})_{IJ}\f$
         * 
         * \f$\frac{ \partial \hat{S}_{IJ} }{ \partial C_{KL} } = \frac{1}{3} \left( C_{AB} S_{AB} (C^{-1})_{IK} (C^{-1})_{LJ} - S_{KL} (C^{-1}_{IJ}) \right)\f$
         *
         * \param &secondOrderReferenceStress: The stress measure in the reference configuration.
         * \param &rightCauchyGreenDeformation: The right Cauchy Green Deformation tensor of the 
         *     deformation between configurations.
         * \param &pressure: The pressure of the reference stress measure.
         * \param &dPressuredStress: The Jacobian of the pressure w.r.t. the reference stress.
         * \param &dPressuredRCG: The Jacobian of the pressure w.r.t. the right Cauchy-Green
         *     deformation measure.
         * \param &deviatoricSecondOrderReferenceStress: The deviatoric part of the second order 
         *     stress in the reference configuration.
         * \param &dDeviatoricReferenceStressdReferenceStress: The jacobian w.r.t. the reference stress.
         * \param &dDeviatoricReferenceStressdRCG: The jacobian w.r.t. the right Cauchy Green deformation
         *     tensor.
         */
        //Assume 3d
        constexpr unsigned int dim = 3;
        constexpr unsigned int sot_dim = dim * dim;

        variableVector invRCG = rightCauchyGreenDeformation;
        Eigen::Map< Eigen::Matrix< variableType, dim, dim, Eigen::RowMajor > > map1( invRCG.data( ), dim, dim );
        map1 = map1.inverse( ).eval( );

        deviatoricSecondOrderReferenceStress = secondOrderReferenceStress - pressure * invRCG;

        //Compute the first order jacobians
        dDeviatoricReferenceStressdReferenceStress = variableVector( sot_dim * sot_dim, 0 );
        dDeviatoricReferenceStressdRCG = variableVector( sot_dim * sot_dim, 0 );
        for ( unsigned int i = 0; i < sot_dim; i++ ){ dDeviatoricReferenceStressdReferenceStress[ sot_dim * i + i ] = 1.; }

        Eigen::Map< const Eigen::Matrix< variableType, 1, sot_dim, Eigen::RowMajor > > map2( dPressuredStress.data( ), 1, sot_dim );
        Eigen::Map< Eigen::Matrix< variableType, sot_dim, sot_dim, Eigen::RowMajor > > map3( dDeviatoricReferenceStressdReferenceStress.data( ), sot_dim, sot_dim );
        Eigen::Map< const Eigen::Matrix< variableType, 1, sot_dim, Eigen::RowMajor > > map4( dPressuredRCG.data( ), 1, sot_dim );
        Eigen::Map< Eigen::Matrix< variableType, sot_dim, sot_dim, Eigen::RowMajor > > map5( dDeviatoricReferenceStressdRCG.data( ), sot_dim, sot_dim );

        map3 = ( map3 - map1.reshaped<Eigen::RowMajor>( ) * map2 ).eval( );

        map5 = ( -map1.reshaped<Eigen::RowMajor>( ) * map4 ).eval( );

        for ( unsigned int I = 0; I < dim; I++ ){
            for ( unsigned int J = 0; J < dim; J++ ){
                for ( unsigned int K = 0; K < dim; K++ ){
                    for ( unsigned int L = 0; L < dim; L++ ){
                        dDeviatoricReferenceStressdRCG[ dim * sot_dim * I + sot_dim * J + dim * K + L ] += pressure * invRCG[ dim * I + K ] * invRCG[ dim * L + J ];
                    }
                }
            }
        }

        return;
    }

    void computeDeviatoricReferenceSecondOrderStress( const variableVector &secondOrderReferenceStress,
                                                          const variableVector &rightCauchyGreenDeformation,
                                                          const variableType &pressure,
                                                          const variableVector &dPressuredStress,
                                                          const variableVector &dPressuredRCG,
                                                          const variableVector &d2PressuredStressdRCG,
                                                          variableVector &deviatoricSecondOrderReferenceStress,
                                                          variableVector &dDeviatoricReferenceStressdReferenceStress,
                                                          variableVector &dDeviatoricReferenceStressdRCG,
                                                          variableVector &d2DevSdSdRCG ){
        /*!
         * Compute the deviatoric part of a second order stress measure in the reference configuration.
         * \f$\hat{S}_{IJ} = S_{IJ} - \frac{1}{3} C_{AB} S_{AB} (C^{-1})_{IJ}\f$
         * 
         * Also compute the Jacobians.
         * \f$\frac{ \partial \hat{S}_{IJ} }{ \partial S_{KL} } = \delta_{IK} \delta_{LJ} - \frac{1}{3} C_{KL} (C^{-1})_{IJ}\f$
         * 
         * \f$\frac{ \partial \hat{S}_{IJ} }{ \partial C_{KL} } = \frac{1}{3} \left( C_{AB} S_{AB} (C^{-1})_{IK} (C^{-1})_{LJ} - S_{KL} (C^{-1}_{IJ}) \right)\f$
         *
         * \f$\frac{ \partial \hat{S}_{IJ} }{ \partial S_{KL} \partial C_{MN} } = - \frac{ \partial p^2 }{ \partial S_{KL} \partial C_{MN}} (C^{-1})_{IJ} + \frac{ \partial p }{ \partial S_{KL} } (C^{-1})_{IM} (C^{-1})_{NJ}\f$
         *
         * \param &secondOrderReferenceStress: The stress measure in the reference configuration.
         * \param &rightCauchyGreenDeformation: The right Cauchy-Green deformation tensor of the 
         *     deformation between configurations.
         * \param &pressure: The pressure of the reference stress.
         * \param &dPressuredStress: The Jacobian of the pressure w.r.t. the stress.
         * \param &dPressuredRCG: The Jacobian of the pressure w.r.t. the right Cauchy-Green
         *     deformation tensor.
         * \param &d2PressuredStressdRCG: The Jacobian of the pressure w.r.t. the stress and 
         *     the right Cauchy-Green deformation tensor.
         * \param &deviatoricSecondOrderReferenceStress: The deviatoric part of the second order 
         *     stress in the reference configuration.
         * \param &dDeviatoricReferenceStressdReferenceStress: The Jacobian w.r.t. the reference stress.
         * \param &dDeviatoricReferenceStressdRCG: The Jacobian w.r.t. the right Cauchy-Green deformation
         *     tensor.
         * \param &d2DevSdSdRCG: The second order mixed Jacobian w.r.t. the stress and right Cauchy-Green 
         *     deformation tensor. Stored [IJ][KLMN]
         */
        //Assume 3d
        constexpr unsigned int dim = 3;
        constexpr unsigned int sot_dim = dim * dim;

        variableVector invRCG = rightCauchyGreenDeformation;
        Eigen::Map< Eigen::Matrix< variableType, dim, dim, Eigen::RowMajor > > map1( invRCG.data( ), dim, dim );
        map1 = map1.inverse( ).eval( );

        deviatoricSecondOrderReferenceStress = secondOrderReferenceStress - pressure * invRCG;

        //Compute the first order jacobians
        dDeviatoricReferenceStressdReferenceStress = variableVector( sot_dim * sot_dim, 0 );
        dDeviatoricReferenceStressdRCG = variableVector( sot_dim * sot_dim, 0 );
        for ( unsigned int i = 0; i < sot_dim; i++ ){ dDeviatoricReferenceStressdReferenceStress[ sot_dim * i + i ] = 1.; }

        Eigen::Map< const Eigen::Matrix< variableType, 1, sot_dim, Eigen::RowMajor > > map2( dPressuredStress.data( ), 1, sot_dim );
        Eigen::Map< Eigen::Matrix< variableType, sot_dim, sot_dim, Eigen::RowMajor > > map3( dDeviatoricReferenceStressdReferenceStress.data( ), sot_dim, sot_dim );
        Eigen::Map< const Eigen::Matrix< variableType, 1, sot_dim, Eigen::RowMajor > > map4( dPressuredRCG.data( ), 1, sot_dim );
        Eigen::Map< Eigen::Matrix< variableType, sot_dim, sot_dim, Eigen::RowMajor > > map5( dDeviatoricReferenceStressdRCG.data( ), sot_dim, sot_dim );

        map3 = ( map3 - map1.reshaped<Eigen::RowMajor>( ) * map2 ).eval( );

        map5 = ( -map1.reshaped<Eigen::RowMajor>( ) * map4 ).eval( );

        for ( unsigned int I = 0; I < dim; I++ ){
            for ( unsigned int J = 0; J < dim; J++ ){
                for ( unsigned int K = 0; K < dim; K++ ){
                    for ( unsigned int L = 0; L < dim; L++ ){
                        dDeviatoricReferenceStressdRCG[ dim * sot_dim * I + sot_dim * J + dim * K + L ] += pressure * invRCG[ dim * I + K ] * invRCG[ dim * L + J ];
                    }
                }
            }
        }

        //Compute the second order jacobians
        d2DevSdSdRCG = variableVector( sot_dim * sot_dim * sot_dim, 0 );
        for ( unsigned int I = 0; I < dim; I++ ){
            for ( unsigned int J = 0; J < dim; J++ ){
                for ( unsigned int K = 0; K < dim; K++ ){
                    for ( unsigned int L = 0; L < dim; L++ ){
                        for ( unsigned int M = 0; M < dim; M++ ){
                            for ( unsigned int N = 0; N < dim; N++ ){
                                d2DevSdSdRCG[ dim * sot_dim * sot_dim * I + sot_dim * sot_dim * J + dim * dim * dim * K + dim * dim * L + dim * M + N ] = 
                                    -d2PressuredStressdRCG[ dim * sot_dim * K + sot_dim * L + dim * M + N ] * invRCG[ dim * I + J ]
                                    +dPressuredStress[ dim * K + L ] * invRCG[ dim * I + M ] * invRCG[ dim * N + J ];
                            }
                        }
                    }
                }
            }
        }

        return;
    }

    void computeDeviatoricReferenceSecondOrderStress( const variableVector &secondOrderReferenceStress,
                                                          const variableVector &rightCauchyGreenDeformation,
                                                          variableVector &deviatoricSecondOrderReferenceStress,
                                                          variableMatrix &dDeviatoricReferenceStressdReferenceStress,
                                                          variableMatrix &dDeviatoricReferenceStressdRCG,
                                                          variableMatrix &d2DevSdSdRCG ){
        /*!
         * Compute the deviatoric part of a second order stress measure in the reference configuration.
         * \f$\hat{S}_{IJ} = S_{IJ} - \frac{1}{3} C_{AB} S_{AB} (C^{-1})_{IJ}\f$
         * 
         * Also compute the Jacobians.
         * \f$\frac{ \partial \hat{S}_{IJ} }{ \partial S_{KL} } = \delta_{IK} \delta_{LJ} - \frac{1}{3} C_{KL} (C^{-1})_{IJ}\f$
         * 
         * \f$\frac{ \partial \hat{S}_{IJ} }{ \partial C_{KL} } = \frac{1}{3} \left( C_{AB} S_{AB} (C^{-1})_{IK} (C^{-1})_{LJ} - S_{KL} (C^{-1}_{IJ}) \right)\f$
         *
         * \f$\frac{ \partial \hat{S}_{IJ} }{ \partial S_{KL} \partial C_{MN} } = - \frac{ \partial p^2 }{ \partial S_{KL} \partial C_{MN}} (C^{-1})_{IJ} + \frac{ \partial p }{ \partial S_{KL} } (C^{-1})_{IM} (C^{-1})_{NJ}\f$
         *
         * \param &secondOrderReferenceStress: The stress measure in the reference configuration.
         * \param &rightCauchyGreenDeformation: The right Cauchy-Green deformation tensor of the 
         *     deformation between configurations.
         * \param &deviatoricSecondOrderReferenceStress: The deviatoric part of the second order 
         *     stress in the reference configuration.
         * \param &dDeviatoricReferenceStressdReferenceStress: The Jacobian w.r.t. the reference stress.
         * \param &dDeviatoricReferenceStressdRCG: The Jacobian w.r.t. the right Cauchy-Green deformation
         *     tensor.
         * \param &d2DevSdSdRCG: The second order mixed Jacobian w.r.t. the stress and right Cauchy-Green 
         *     deformation tensor.
         */

        variableVector _dDeviatoricReferenceStressdReferenceStress;
        variableVector _dDeviatoricReferenceStressdRCG;
        variableVector _d2DevSdSdRCG;

        TARDIGRADE_ERROR_TOOLS_CATCH( computeDeviatoricReferenceSecondOrderStress( secondOrderReferenceStress,
                                                                                   rightCauchyGreenDeformation,
                                                                                   deviatoricSecondOrderReferenceStress,
                                                                                   _dDeviatoricReferenceStressdReferenceStress,
                                                                                   _dDeviatoricReferenceStressdRCG,
                                                                                   _d2DevSdSdRCG ) );

        dDeviatoricReferenceStressdReferenceStress = tardigradeVectorTools::inflate( _dDeviatoricReferenceStressdReferenceStress , 9,  9 );
        dDeviatoricReferenceStressdRCG             = tardigradeVectorTools::inflate( _dDeviatoricReferenceStressdRCG             , 9,  9 );
        d2DevSdSdRCG                               = tardigradeVectorTools::inflate( _d2DevSdSdRCG                               , 9, 81 );

        return;

    }
    void computeDeviatoricReferenceSecondOrderStress( const variableVector &secondOrderReferenceStress,
                                                          const variableVector &rightCauchyGreenDeformation,
                                                          variableVector &deviatoricSecondOrderReferenceStress,
                                                          variableVector &dDeviatoricReferenceStressdReferenceStress,
                                                          variableVector &dDeviatoricReferenceStressdRCG,
                                                          variableVector &d2DevSdSdRCG ){
        /*!
         * Compute the deviatoric part of a second order stress measure in the reference configuration.
         * \f$\hat{S}_{IJ} = S_{IJ} - \frac{1}{3} C_{AB} S_{AB} (C^{-1})_{IJ}\f$
         * 
         * Also compute the Jacobians.
         * \f$\frac{ \partial \hat{S}_{IJ} }{ \partial S_{KL} } = \delta_{IK} \delta_{LJ} - \frac{1}{3} C_{KL} (C^{-1})_{IJ}\f$
         * 
         * \f$\frac{ \partial \hat{S}_{IJ} }{ \partial C_{KL} } = \frac{1}{3} \left( C_{AB} S_{AB} (C^{-1})_{IK} (C^{-1})_{LJ} - S_{KL} (C^{-1}_{IJ}) \right)\f$
         *
         * \f$\frac{ \partial \hat{S}_{IJ} }{ \partial S_{KL} \partial C_{MN} } = - \frac{ \partial p^2 }{ \partial S_{KL} \partial C_{MN}} (C^{-1})_{IJ} + \frac{ \partial p }{ \partial S_{KL} } (C^{-1})_{IM} (C^{-1})_{NJ}\f$
         *
         * \param &secondOrderReferenceStress: The stress measure in the reference configuration.
         * \param &rightCauchyGreenDeformation: The right Cauchy-Green deformation tensor of the 
         *     deformation between configurations.
         * \param &deviatoricSecondOrderReferenceStress: The deviatoric part of the second order 
         *     stress in the reference configuration.
         * \param &dDeviatoricReferenceStressdReferenceStress: The Jacobian w.r.t. the reference stress.
         * \param &dDeviatoricReferenceStressdRCG: The Jacobian w.r.t. the right Cauchy-Green deformation
         *     tensor.
         * \param &d2DevSdSdRCG: The second order mixed Jacobian w.r.t. the stress and right Cauchy-Green 
         *     deformation tensor.
         */


        variableType pressure;
        variableVector dpdS, dpdC;
        variableVector d2pdSdC;
        TARDIGRADE_ERROR_TOOLS_CATCH( computeReferenceSecondOrderStressPressure( secondOrderReferenceStress, rightCauchyGreenDeformation, pressure,
                                                                                 dpdS, dpdC, d2pdSdC ) );

        return computeDeviatoricReferenceSecondOrderStress( secondOrderReferenceStress, rightCauchyGreenDeformation,
                                                            pressure, dpdS, dpdC, d2pdSdC,
                                                            deviatoricSecondOrderReferenceStress,
                                                            dDeviatoricReferenceStressdReferenceStress,
                                                            dDeviatoricReferenceStressdRCG,
                                                            d2DevSdSdRCG );
    }

    void computeSecondOrderReferenceStressDecomposition( const variableVector &secondOrderReferenceStress,
                                                             const variableVector &rightCauchyGreenDeformation,
                                                             variableVector &deviatoricSecondOrderReferenceStress,
                                                             variableType &pressure ){
        /*!
         * Compute the decomposition of a second-order stress measure into pressure and deviatoric parts.
         *
         * \param &secondOrderReferenceStress: The second-order stress in the reference
         *     configuration.
         * \param &rightCauchyGreenDeformation: The right Cauchy-Green deformation tensor 
         *     between the current configuration and the configuration the stress is located in.
         * \param &deviatoricSecondOrderReferenceStress: The deviatoric part of the second order 
         *     stress tensor.
         * \param &pressure: The pressure of the stress tensor.
         */

        TARDIGRADE_ERROR_TOOLS_CATCH( computeReferenceSecondOrderStressPressure( secondOrderReferenceStress,
                                                                                 rightCauchyGreenDeformation,
                                                                                 pressure ) );

        TARDIGRADE_ERROR_TOOLS_CATCH( computeDeviatoricReferenceSecondOrderStress( secondOrderReferenceStress,
                                                                                   rightCauchyGreenDeformation,
                                                                                   pressure,
                                                                                   deviatoricSecondOrderReferenceStress ) );

        return;
    }

    void computeSecondOrderReferenceStressDecomposition( const variableVector &secondOrderReferenceStress,
                                                             const variableVector &rightCauchyGreenDeformation,
                                                             variableVector &deviatoricSecondOrderReferenceStress,
                                                             variableType &pressure, variableMatrix &dDevStressdStress,
                                                             variableMatrix &dDevStressdRCG, variableVector &dPressuredStress,
                                                             variableVector &dPressuredRCG ){
        /*!
         * Compute the decomposition of a second-order stress measure into pressure and deviatoric parts.
         *
         * Includes the first order Jacobians
         *
         * \param &secondOrderReferenceStress: The second-order stress in the reference
         *     configuration.
         * \param &rightCauchyGreenDeformation: The right Cauchy-Green deformation tensor 
         *     between the current configuration and the configuration the stress is located in.
         * \param &deviatoricSecondOrderReferenceStress: The deviatoric part of the second order 
         *     stress tensor.
         * \param &pressure: The pressure of the stress tensor.
         * \param &dDevStressdStress: The Jacobian of the deviatoric part of the reference 
         *     stress w.r.t. the reference stress.
         * \param &dDevStressdRCG: The Jacobian of the deviatoric part of the reference stress
         *     w.r.t. the right Cauchy-Green deformation tensor.
         * \param &dPressuredStress: The Jacobian of the pressure w.r.t. the reference stress.
         * \param &dPressuredRCG: The Jacobian of the pressure w.r.t. the right Cauchy-Green 
         *     deformation tensor.
         */

        variableVector _dDevStressdStress;
        variableVector _dDevStressdRCG;

        TARDIGRADE_ERROR_TOOLS_CATCH( computeSecondOrderReferenceStressDecomposition( secondOrderReferenceStress,
                                                                                      rightCauchyGreenDeformation,
                                                                                      deviatoricSecondOrderReferenceStress,
                                                                                      pressure, _dDevStressdStress,
                                                                                      _dDevStressdRCG, dPressuredStress,
                                                                                      dPressuredRCG ) );

        dDevStressdStress = tardigradeVectorTools::inflate( _dDevStressdStress, 9, 9 );
        dDevStressdRCG    = tardigradeVectorTools::inflate( _dDevStressdRCG,    9, 9 );

        return;

    }

    void computeSecondOrderReferenceStressDecomposition( const variableVector &secondOrderReferenceStress,
                                                             const variableVector &rightCauchyGreenDeformation,
                                                             variableVector &deviatoricSecondOrderReferenceStress,
                                                             variableType &pressure, variableVector &dDevStressdStress,
                                                             variableVector &dDevStressdRCG, variableVector &dPressuredStress,
                                                             variableVector &dPressuredRCG ){
        /*!
         * Compute the decomposition of a second-order stress measure into pressure and deviatoric parts.
         *
         * Includes the first order Jacobians
         *
         * \param &secondOrderReferenceStress: The second-order stress in the reference
         *     configuration.
         * \param &rightCauchyGreenDeformation: The right Cauchy-Green deformation tensor 
         *     between the current configuration and the configuration the stress is located in.
         * \param &deviatoricSecondOrderReferenceStress: The deviatoric part of the second order 
         *     stress tensor.
         * \param &pressure: The pressure of the stress tensor.
         * \param &dDevStressdStress: The Jacobian of the deviatoric part of the reference 
         *     stress w.r.t. the reference stress.
         * \param &dDevStressdRCG: The Jacobian of the deviatoric part of the reference stress
         *     w.r.t. the right Cauchy-Green deformation tensor.
         * \param &dPressuredStress: The Jacobian of the pressure w.r.t. the reference stress.
         * \param &dPressuredRCG: The Jacobian of the pressure w.r.t. the right Cauchy-Green 
         *     deformation tensor.
         */

        TARDIGRADE_ERROR_TOOLS_CATCH( computeReferenceSecondOrderStressPressure( secondOrderReferenceStress,
                                                                                 rightCauchyGreenDeformation,
                                                                                 pressure, dPressuredStress, dPressuredRCG ) );

        TARDIGRADE_ERROR_TOOLS_CATCH( computeDeviatoricReferenceSecondOrderStress( secondOrderReferenceStress,
                                                                                   rightCauchyGreenDeformation,
                                                                                   pressure, dPressuredStress, dPressuredRCG,
                                                                                   deviatoricSecondOrderReferenceStress,
                                                                                   dDevStressdStress, dDevStressdRCG ) );

        return;
    }

    void computeSecondOrderReferenceStressDecomposition( const variableVector &secondOrderReferenceStress,
                                                             const variableVector &rightCauchyGreenDeformation,
                                                             variableVector &deviatoricSecondOrderReferenceStress,
                                                             variableType &pressure, variableMatrix &dDevStressdStress,
                                                             variableMatrix &dDevStressdRCG, variableVector &dPressuredStress,
                                                             variableVector &dPressuredRCG, variableMatrix &d2DevStressdStressdRCG,
                                                             variableMatrix &d2PressuredStressdRCG ){
        /*!
         * Compute the decomposition of a second-order stress measure into pressure and deviatoric parts.
         *
         * Includes the first and some of the second order Jacobians
         *
         * \param &secondOrderReferenceStress: The second-order stress in the reference
         *     configuration.
         * \param &rightCauchyGreenDeformation: The right Cauchy-Green deformation tensor 
         *     between the current configuration and the configuration the stress is located in.
         * \param &deviatoricSecondOrderReferenceStress: The deviatoric part of the second order 
         *     stress tensor.
         * \param &pressure: The pressure of the stress tensor.
         * \param &dDevStressdStress: The Jacobian of the deviatoric part of the reference 
         *     stress w.r.t. the reference stress.
         * \param &dDevStressdRCG: The Jacobian of the deviatoric part of the reference stress
         *     w.r.t. the right Cauchy-Green deformation tensor.
         * \param &dPressuredStress: The Jacobian of the pressure w.r.t. the reference stress.
         * \param &dPressuredRCG: The Jacobian of the pressure w.r.t. the right Cauchy-Green 
         *     deformation tensor.
         * \param &d2DevStressdStressdRCG: The second order Jacobian of the deviatoric part of the 
         *     reference stress w.r.t. the reference stress and the right Cauchy-Green deformation tensor. 
         *     \f$\frac{ \partial^2 dev( S )_{IJ} }{ \partial S_{KL} \partial C_{MN} }\f$
         * \param &d2PressuredStressdRCG: The second order Jacobian of the deviatoric part of the 
         *     reference stress w.r.t. the reference stress and the right Cauchy-Green deformation tensor.
         */

        variableVector _dDevStressdStress;
        variableVector _dDevStressdRCG;
        variableVector _d2DevStressdStressdRCG;
        variableVector _d2PressuredStressdRCG;

        TARDIGRADE_ERROR_TOOLS_CATCH( computeSecondOrderReferenceStressDecomposition( secondOrderReferenceStress,
                                                                                      rightCauchyGreenDeformation,
                                                                                      deviatoricSecondOrderReferenceStress,
                                                                                      pressure, _dDevStressdStress,
                                                                                      _dDevStressdRCG,dPressuredStress,
                                                                                      dPressuredRCG, _d2DevStressdStressdRCG,
                                                                                      _d2PressuredStressdRCG ) );

        dDevStressdStress      = tardigradeVectorTools::inflate( _dDevStressdStress     , 9,  9 );
        dDevStressdRCG         = tardigradeVectorTools::inflate( _dDevStressdRCG        , 9,  9 );
        d2DevStressdStressdRCG = tardigradeVectorTools::inflate( _d2DevStressdStressdRCG, 9, 81 );
        d2PressuredStressdRCG  = tardigradeVectorTools::inflate( _d2PressuredStressdRCG , 9,  9 );

        return;

    }

    void computeSecondOrderReferenceStressDecomposition( const variableVector &secondOrderReferenceStress,
                                                             const variableVector &rightCauchyGreenDeformation,
                                                             variableVector &deviatoricSecondOrderReferenceStress,
                                                             variableType &pressure, variableVector &dDevStressdStress,
                                                             variableVector &dDevStressdRCG, variableVector &dPressuredStress,
                                                             variableVector &dPressuredRCG, variableVector &d2DevStressdStressdRCG,
                                                             variableVector &d2PressuredStressdRCG ){
        /*!
         * Compute the decomposition of a second-order stress measure into pressure and deviatoric parts.
         *
         * Includes the first and some of the second order Jacobians
         *
         * \param &secondOrderReferenceStress: The second-order stress in the reference
         *     configuration.
         * \param &rightCauchyGreenDeformation: The right Cauchy-Green deformation tensor 
         *     between the current configuration and the configuration the stress is located in.
         * \param &deviatoricSecondOrderReferenceStress: The deviatoric part of the second order 
         *     stress tensor.
         * \param &pressure: The pressure of the stress tensor.
         * \param &dDevStressdStress: The Jacobian of the deviatoric part of the reference 
         *     stress w.r.t. the reference stress.
         * \param &dDevStressdRCG: The Jacobian of the deviatoric part of the reference stress
         *     w.r.t. the right Cauchy-Green deformation tensor.
         * \param &dPressuredStress: The Jacobian of the pressure w.r.t. the reference stress.
         * \param &dPressuredRCG: The Jacobian of the pressure w.r.t. the right Cauchy-Green 
         *     deformation tensor.
         * \param &d2DevStressdStressdRCG: The second order Jacobian of the deviatoric part of the 
         *     reference stress w.r.t. the reference stress and the right Cauchy-Green deformation tensor. 
         *     \f$\frac{ \partial^2 dev( S )_{IJ} }{ \partial S_{KL} \partial C_{MN} }\f$
         * \param &d2PressuredStressdRCG: The second order Jacobian of the deviatoric part of the 
         *     reference stress w.r.t. the reference stress and the right Cauchy-Green deformation tensor.
         */

        TARDIGRADE_ERROR_TOOLS_CATCH( computeReferenceSecondOrderStressPressure( secondOrderReferenceStress,
                                                                                 rightCauchyGreenDeformation,
                                                                                 pressure, dPressuredStress, dPressuredRCG,
                                                                                 d2PressuredStressdRCG ) );

        TARDIGRADE_ERROR_TOOLS_CATCH( computeDeviatoricReferenceSecondOrderStress( secondOrderReferenceStress,
                                                                                   rightCauchyGreenDeformation,
                                                                                   pressure, dPressuredStress, dPressuredRCG,
                                                                                   d2PressuredStressdRCG,
                                                                                   deviatoricSecondOrderReferenceStress,
                                                                                   dDevStressdStress, dDevStressdRCG,
                                                                                   d2DevStressdStressdRCG ) );

        return;
    }

    void computeHigherOrderReferenceStressDecomposition( const variableVector &higherOrderReferenceStress,
                                                             const variableVector &rightCauchyGreenDeformation,
                                                             variableVector &deviatoricHigherOrderReferenceStress,
                                                             variableVector &pressure ){
        /*!
         * Compute the decomposition of the higher-order stress measure into pressure and deviatoric parts.
         *
         * \param &higherOrderReferenceStress: The higher order stress in the reference
         *     configuration.
         * \param &rightCauchyGreenDeformation: The right Cauchy-Green deformation tensor 
         *     between the current configuration and the configuration the stress is located in.
         * \param &deviatoricHigherOrderReferenceStress: The deviatoric part of the higher order 
         *     stress tensor.
         * \param &pressure: The pressure of the higher order stress tensor.
         */

        TARDIGRADE_ERROR_TOOLS_CATCH( computeReferenceHigherOrderStressPressure( higherOrderReferenceStress,
                                                                                 rightCauchyGreenDeformation,
                                                                                 pressure ) );

        TARDIGRADE_ERROR_TOOLS_CATCH( computeDeviatoricReferenceHigherOrderStress( higherOrderReferenceStress,
                                                                                   rightCauchyGreenDeformation,
                                                                                   pressure,
                                                                                   deviatoricHigherOrderReferenceStress ) );

        return;
    }

    void computeHigherOrderReferenceStressDecomposition( const variableVector &higherOrderReferenceStress,
                                                             const variableVector &rightCauchyGreenDeformation,
                                                             variableVector &deviatoricHigherOrderReferenceStress,
                                                             variableVector &pressure, variableMatrix &dDevStressdStress,
                                                             variableMatrix &dDevStressdRCG, variableMatrix &dPressuredStress,
                                                             variableMatrix &dPressuredRCG ){
        /*!
         * Compute the decomposition of the higher-order stress measure into pressure and deviatoric parts.
         *
         * Also return the Jacobians
         *
         * \param &higherOrderReferenceStress: The higher order stress in the reference
         *     configuration.
         * \param &rightCauchyGreenDeformation: The right Cauchy-Green deformation tensor 
         *     between the current configuration and the configuration the stress is located in.
         * \param &deviatoricHigherOrderReferenceStress: The deviatoric part of the higher order 
         *     stress tensor.
         * \param &pressure: The pressure of the higher order stress tensor.
         * \param &dDevStressdStress: The Jacobian of the deviatoric stress w.r.t. the stress.
         * \param &dDevStressdRCG: The Jacobian of the deviatoric stress w.r.t. the right 
         *     Cauchy-Green deformation tensor.
         * \param &dPressuredStress: The Jacobian of the pressure w.r.t. the stress.
         * \param &dPressuredRCG: The Jacobian of the pressure w.r.t. the right 
         *     Cauchy-Green deformation tensor.
         */

        variableVector _dDevStressdStress;
        variableVector _dDevStressdRCG;
        variableVector _dPressuredStress;
        variableVector _dPressuredRCG;

        TARDIGRADE_ERROR_TOOLS_CATCH( computeHigherOrderReferenceStressDecomposition( higherOrderReferenceStress,
                                                                                      rightCauchyGreenDeformation,
                                                                                      deviatoricHigherOrderReferenceStress,
                                                                                      pressure, _dDevStressdStress,
                                                                                      _dDevStressdRCG, _dPressuredStress,
                                                                                      _dPressuredRCG ) );

        dDevStressdStress = tardigradeVectorTools::inflate( _dDevStressdStress , 27, 27 );
        dDevStressdRCG    = tardigradeVectorTools::inflate( _dDevStressdRCG    , 27,  9 );
        dPressuredStress  = tardigradeVectorTools::inflate( _dPressuredStress  ,  3, 27 );
        dPressuredRCG     = tardigradeVectorTools::inflate( _dPressuredRCG     ,  3,  9 );

        return;

    }

    void computeHigherOrderReferenceStressDecomposition( const variableVector &higherOrderReferenceStress,
                                                             const variableVector &rightCauchyGreenDeformation,
                                                             variableVector &deviatoricHigherOrderReferenceStress,
                                                             variableVector &pressure, variableVector &dDevStressdStress,
                                                             variableVector &dDevStressdRCG, variableVector &dPressuredStress,
                                                             variableVector &dPressuredRCG ){
        /*!
         * Compute the decomposition of the higher-order stress measure into pressure and deviatoric parts.
         *
         * Also return the Jacobians
         *
         * \param &higherOrderReferenceStress: The higher order stress in the reference
         *     configuration.
         * \param &rightCauchyGreenDeformation: The right Cauchy-Green deformation tensor 
         *     between the current configuration and the configuration the stress is located in.
         * \param &deviatoricHigherOrderReferenceStress: The deviatoric part of the higher order 
         *     stress tensor.
         * \param &pressure: The pressure of the higher order stress tensor.
         * \param &dDevStressdStress: The Jacobian of the deviatoric stress w.r.t. the stress.
         * \param &dDevStressdRCG: The Jacobian of the deviatoric stress w.r.t. the right 
         *     Cauchy-Green deformation tensor.
         * \param &dPressuredStress: The Jacobian of the pressure w.r.t. the stress.
         * \param &dPressuredRCG: The Jacobian of the pressure w.r.t. the right 
         *     Cauchy-Green deformation tensor.
         */

        TARDIGRADE_ERROR_TOOLS_CATCH( computeReferenceHigherOrderStressPressure( higherOrderReferenceStress,
                                                                                 rightCauchyGreenDeformation,
                                                                                 pressure, dPressuredStress,
                                                                                 dPressuredRCG ) );

        TARDIGRADE_ERROR_TOOLS_CATCH( computeDeviatoricReferenceHigherOrderStress( higherOrderReferenceStress,
                                                                                   rightCauchyGreenDeformation,
                                                                                   pressure, dPressuredStress,
                                                                                   dPressuredRCG,
                                                                                   deviatoricHigherOrderReferenceStress,
                                                                                   dDevStressdStress, dDevStressdRCG ) );

        return;
    }

    void computeHigherOrderReferenceStressDecomposition( const variableVector &higherOrderReferenceStress,
                                                             const variableVector &rightCauchyGreenDeformation,
                                                             variableVector &deviatoricHigherOrderReferenceStress,
                                                             variableVector &pressure, variableMatrix &dDevStressdStress,
                                                             variableMatrix &dDevStressdRCG, variableMatrix &dPressuredStress,
                                                             variableMatrix &dPressuredRCG, variableMatrix &d2DevStressdStressdRCG,
                                                             variableMatrix &d2PressuredStressdRCG ){
        /*!
         * Compute the decomposition of the higher-order stress measure into pressure and deviatoric parts.
         *
         * Also return the Jacobians
         *
         * \param &higherOrderReferenceStress: The higher order stress in the reference
         *     configuration.
         * \param &rightCauchyGreenDeformation: The right Cauchy-Green deformation tensor 
         *     between the current configuration and the configuration the stress is located in.
         * \param &deviatoricHigherOrderReferenceStress: The deviatoric part of the higher order 
         *     stress tensor.
         * \param &pressure: The pressure of the higher order stress tensor.
         * \param &dDevStressdStress: The Jacobian of the deviatoric stress w.r.t. the stress.
         * \param &dDevStressdRCG: The Jacobian of the deviatoric stress w.r.t. the right 
         *     Cauchy-Green deformation tensor.
         * \param &dPressuredStress: The Jacobian of the pressure w.r.t. the stress.
         * \param &dPressuredRCG: The Jacobian of the pressure w.r.t. the right 
         *     Cauchy-Green deformation tensor.
         * \param &d2DevStressdStressdRCG: The second order jacobian of the deviatoric stress 
         *     w.r.t. the stress and the right Cauchy-Green deformation tensor.
         * \param &d2PressuredStressdRCG: The second order jacobian of the pressure
         *     w.r.t. the stress and the right Cauchy-Green deformation tensor.
         */

        variableVector _dDevStressdStress;
        variableVector _dDevStressdRCG;
        variableVector _dPressuredStress;
        variableVector _dPressuredRCG;
        variableVector _d2DevStressdStressdRCG;
        variableVector _d2PressuredStressdRCG;

        TARDIGRADE_ERROR_TOOLS_CATCH( computeHigherOrderReferenceStressDecomposition( higherOrderReferenceStress,
                                                                                      rightCauchyGreenDeformation,
                                                                                      deviatoricHigherOrderReferenceStress,
                                                                                      pressure, _dDevStressdStress,
                                                                                      _dDevStressdRCG, _dPressuredStress,
                                                                                      _dPressuredRCG, _d2DevStressdStressdRCG,
                                                                                      _d2PressuredStressdRCG ) );

        dDevStressdStress      = tardigradeVectorTools::inflate( _dDevStressdStress     , 27,  27 );
        dDevStressdRCG         = tardigradeVectorTools::inflate( _dDevStressdRCG        , 27,   9 );
        dPressuredStress       = tardigradeVectorTools::inflate( _dPressuredStress      ,  3,  27 );
        dPressuredRCG          = tardigradeVectorTools::inflate( _dPressuredRCG         ,  3,   9 );
        d2DevStressdStressdRCG = tardigradeVectorTools::inflate( _d2DevStressdStressdRCG, 27, 243 );
        d2PressuredStressdRCG  = tardigradeVectorTools::inflate( _d2PressuredStressdRCG ,  3, 243 );

        return;

    }

    void computeHigherOrderReferenceStressDecomposition( const variableVector &higherOrderReferenceStress,
                                                             const variableVector &rightCauchyGreenDeformation,
                                                             variableVector &deviatoricHigherOrderReferenceStress,
                                                             variableVector &pressure, variableVector &dDevStressdStress,
                                                             variableVector &dDevStressdRCG, variableVector &dPressuredStress,
                                                             variableVector &dPressuredRCG, variableVector &d2DevStressdStressdRCG,
                                                             variableVector &d2PressuredStressdRCG ){
        /*!
         * Compute the decomposition of the higher-order stress measure into pressure and deviatoric parts.
         *
         * Also return the Jacobians
         *
         * \param &higherOrderReferenceStress: The higher order stress in the reference
         *     configuration.
         * \param &rightCauchyGreenDeformation: The right Cauchy-Green deformation tensor 
         *     between the current configuration and the configuration the stress is located in.
         * \param &deviatoricHigherOrderReferenceStress: The deviatoric part of the higher order 
         *     stress tensor.
         * \param &pressure: The pressure of the higher order stress tensor.
         * \param &dDevStressdStress: The Jacobian of the deviatoric stress w.r.t. the stress.
         * \param &dDevStressdRCG: The Jacobian of the deviatoric stress w.r.t. the right 
         *     Cauchy-Green deformation tensor.
         * \param &dPressuredStress: The Jacobian of the pressure w.r.t. the stress.
         * \param &dPressuredRCG: The Jacobian of the pressure w.r.t. the right 
         *     Cauchy-Green deformation tensor.
         * \param &d2DevStressdStressdRCG: The second order jacobian of the deviatoric stress 
         *     w.r.t. the stress and the right Cauchy-Green deformation tensor.
         * \param &d2PressuredStressdRCG: The second order jacobian of the pressure
         *     w.r.t. the stress and the right Cauchy-Green deformation tensor.
         */

        TARDIGRADE_ERROR_TOOLS_CATCH( computeReferenceHigherOrderStressPressure( higherOrderReferenceStress,
                                                                                 rightCauchyGreenDeformation,
                                                                                 pressure, dPressuredStress,
                                                                                 dPressuredRCG, d2PressuredStressdRCG ) );

        TARDIGRADE_ERROR_TOOLS_CATCH( computeDeviatoricReferenceHigherOrderStress( higherOrderReferenceStress,
                                                                                   rightCauchyGreenDeformation,
                                                                                   pressure, dPressuredStress,
                                                                                   dPressuredRCG, d2PressuredStressdRCG,
                                                                                   deviatoricHigherOrderReferenceStress,
                                                                                   dDevStressdStress, dDevStressdRCG,
                                                                                   d2DevStressdStressdRCG ) );

        return;
    }

    void computeHigherOrderStressNorm( const variableVector &higherOrderStress, variableVector &higherOrderStressNorm ){
        /*!
         * Compute the norm of the higher order stress which is defined as
         * \f$|| M ||_K = \sqrt{ M_{IJK} M_{IJK} }\f$
         *
         * where K is not summed over.
         *
         * \param &higherOrderStress: The higher order stress tensor.
         * \param &higherOrderStressNorm: The norm of the higher order stress.
         */

        //Assume 3D
        constexpr unsigned int dim = 3;
        constexpr unsigned int sot_dim = dim * dim;
        constexpr unsigned int tot_dim = sot_dim * dim;

        TARDIGRADE_ERROR_TOOLS_CHECK( higherOrderStress.size() == tot_dim,  "The higher order stress does not have the expected dimension" );
        
        higherOrderStressNorm = variableVector( dim, 0 );
        for ( unsigned int K = 0; K < 3; K++ ){
            for ( unsigned int I = 0; I < 3; I++ ){
                for ( unsigned int J = 0; J < 3; J++ ){
                     higherOrderStressNorm[ K ] += higherOrderStress[ dim * dim * I + dim * J + K ]
                                                 * higherOrderStress[ dim * dim * I + dim * J + K ];
                }
            }
            higherOrderStressNorm[ K ] = std::sqrt( higherOrderStressNorm[ K ] );
        }

        return;
    }

    void computeHigherOrderStressNorm( const variableVector &higherOrderStress, variableVector &higherOrderStressNorm,
                                           variableMatrix &dHigherOrderStressNormdHigherOrderStress,
                                           double tol ){
        /*!
         * Compute the norm of the higher order stress which is defined as
         * \f$|| M ||_K = \sqrt{ M_{IJK} M_{IJK} }\f$
         *
         * where K is not summed over.
         *
         * Also computes the Jacobians
         *
         * \f$\frac{ \partial || M ||_K }{ \partial M_{LMN} } = \frac{ M_{LMK} \delta_{KN} }{ || M ||_K }\f$
         *
         * where K is not summed over
         *
         * \param &higherOrderStress: The higher order stress tensor.
         * \param &higherOrderStressNorm: The norm of the higher order stress.
         * \param &dHigherOrderStressNormdHigherOrderStress: The Jacobian of the higher order
         *     stress norm w.r.t. the higher order stress.
         * \param tol: The tolerance of the higher order stress norm when computing the derivative.
         *     Prevents nans.
         */

        variableVector _dHigherOrderStressNormdHigherOrderStress;

        TARDIGRADE_ERROR_TOOLS_CATCH( computeHigherOrderStressNorm( higherOrderStress, higherOrderStressNorm,
                                                                    _dHigherOrderStressNormdHigherOrderStress,
                                                                    tol ) );

        dHigherOrderStressNormdHigherOrderStress = tardigradeVectorTools::inflate( _dHigherOrderStressNormdHigherOrderStress, 3, 27 );

        return;

    }
    void computeHigherOrderStressNorm( const variableVector &higherOrderStress, variableVector &higherOrderStressNorm,
                                           variableVector &dHigherOrderStressNormdHigherOrderStress,
                                           double tol ){
        /*!
         * Compute the norm of the higher order stress which is defined as
         * \f$|| M ||_K = \sqrt{ M_{IJK} M_{IJK} }\f$
         *
         * where K is not summed over.
         *
         * Also computes the Jacobians
         *
         * \f$\frac{ \partial || M ||_K }{ \partial M_{LMN} } = \frac{ M_{LMK} \delta_{KN} }{ || M ||_K }\f$
         *
         * where K is not summed over
         *
         * \param &higherOrderStress: The higher order stress tensor.
         * \param &higherOrderStressNorm: The norm of the higher order stress.
         * \param &dHigherOrderStressNormdHigherOrderStress: The Jacobian of the higher order
         *     stress norm w.r.t. the higher order stress.
         * \param tol: The tolerance of the higher order stress norm when computing the derivative.
         *     Prevents nans.
         */

        //Assume 3D
        constexpr unsigned int dim = 3;
        constexpr unsigned int sot_dim = dim * dim;
        constexpr unsigned int tot_dim = sot_dim * dim;

        TARDIGRADE_ERROR_TOOLS_CATCH( computeHigherOrderStressNorm( higherOrderStress, higherOrderStressNorm ) );

        dHigherOrderStressNormdHigherOrderStress = variableVector( dim * tot_dim, 0 );
        for ( unsigned int K = 0; K < 3; K++ ){
            for ( unsigned int L = 0; L < 3; L++ ){
                for ( unsigned int M = 0; M < 3; M++ ){
                        dHigherOrderStressNormdHigherOrderStress[ tot_dim * K + dim * dim * L + dim * M + K ] = higherOrderStress[ dim * dim * L + dim * M + K ] / ( higherOrderStressNorm[ K ] + tol );
                }
            }
        }

        return;
    }

    void computeHigherOrderStressNorm( const variableVector &higherOrderStress, variableVector &higherOrderStressNorm,
                                           variableMatrix &dHigherOrderStressNormdHigherOrderStress,
                                           variableMatrix &d2HigherOrderStressNormdHigherOrderStress2,
                                           double tol ){
        /*!
         * Compute the norm of the higher order stress which is defined as
         * \f$|| M ||_K = \sqrt{ M_{IJK} M_{IJK} }\f$
         *
         * where K is not summed over.
         *
         * Also computes the Jacobians
         *
         * \f$\frac{ \partial || M ||_K }{ \partial M_{LMN} } = \frac{ M_{LMK} \delta_{KN} }{ || M ||_K }\f$
         * 
         * \f$\frac{ \partial^2 || M ||_K }{ \partial M_{LMN} \partial M_{OPQ} } = \frac{1}{ || M ||_K } \left[ \delta_{LO} \delta_{MP} \delta_{KQ} \delta_{KN} - \frac{ M_{LMK} \delta_{KN} }{ || M ||_K } \frac{ M_{OPK} \delta_{KQ} }{ || M ||_K } \right]\f$
         *
         * where K is not summed over
         *
         * \param &higherOrderStress: The higher order stress tensor.
         * \param &higherOrderStressNorm: The norm of the higher order stress.
         * \param &dHigherOrderStressNormdHigherOrderStress: The Jacobian of the higher order
         *     stress norm w.r.t. the higher order stress.
         * \param &d2HigherOrderStressNormdHigherOrderStress2: The second order Jacobian of the 
         *     higher order stress norm w.r.t. the higher order stress.
         * \param tol: The tolerance of the higher order stress norm when computing the derivatives.
         *     Prevents nans.
         */

        variableVector _dHigherOrderStressNormdHigherOrderStress;
        variableVector _d2HigherOrderStressNormdHigherOrderStress2;

        TARDIGRADE_ERROR_TOOLS_CATCH( computeHigherOrderStressNorm( higherOrderStress, higherOrderStressNorm,
                                                                    _dHigherOrderStressNormdHigherOrderStress,
                                                                    _d2HigherOrderStressNormdHigherOrderStress2,
                                                                    tol ) );

        dHigherOrderStressNormdHigherOrderStress   = tardigradeVectorTools::inflate( _dHigherOrderStressNormdHigherOrderStress, 3, 27 );
        d2HigherOrderStressNormdHigherOrderStress2 = tardigradeVectorTools::inflate( _d2HigherOrderStressNormdHigherOrderStress2, 3, 27 * 27 );

        return;

    }

    void computeHigherOrderStressNorm( const variableVector &higherOrderStress, variableVector &higherOrderStressNorm,
                                           variableVector &dHigherOrderStressNormdHigherOrderStress,
                                           variableVector &d2HigherOrderStressNormdHigherOrderStress2,
                                           double tol ){
        /*!
         * Compute the norm of the higher order stress which is defined as
         * \f$|| M ||_K = \sqrt{ M_{IJK} M_{IJK} }\f$
         *
         * where K is not summed over.
         *
         * Also computes the Jacobians
         *
         * \f$\frac{ \partial || M ||_K }{ \partial M_{LMN} } = \frac{ M_{LMK} \delta_{KN} }{ || M ||_K }\f$
         * 
         * \f$\frac{ \partial^2 || M ||_K }{ \partial M_{LMN} \partial M_{OPQ} } = \frac{1}{ || M ||_K } \left[ \delta_{LO} \delta_{MP} \delta_{KQ} \delta_{KN} - \frac{ M_{LMK} \delta_{KN} }{ || M ||_K } \frac{ M_{OPK} \delta_{KQ} }{ || M ||_K } \right]\f$
         *
         * where K is not summed over
         *
         * \param &higherOrderStress: The higher order stress tensor.
         * \param &higherOrderStressNorm: The norm of the higher order stress.
         * \param &dHigherOrderStressNormdHigherOrderStress: The Jacobian of the higher order
         *     stress norm w.r.t. the higher order stress.
         * \param &d2HigherOrderStressNormdHigherOrderStress2: The second order Jacobian of the 
         *     higher order stress norm w.r.t. the higher order stress.
         * \param tol: The tolerance of the higher order stress norm when computing the derivatives.
         *     Prevents nans.
         */

        //Assume 3D
        constexpr unsigned int dim = 3;
        constexpr unsigned int sot_dim = dim * dim;
        constexpr unsigned int tot_dim = sot_dim * dim;

        TARDIGRADE_ERROR_TOOLS_CATCH( computeHigherOrderStressNorm( higherOrderStress, higherOrderStressNorm, dHigherOrderStressNormdHigherOrderStress, tol ) );

        d2HigherOrderStressNormdHigherOrderStress2 = variableVector( dim * tot_dim * tot_dim, 0 );

        for ( unsigned int K = 0; K < dim; K++ ){
            for ( unsigned int L = 0; L < dim; L++ ){
                for ( unsigned int M = 0; M < dim; M++ ){
                    d2HigherOrderStressNormdHigherOrderStress2[ tot_dim * tot_dim * K + dim * dim * dim * dim * dim * L + dim * dim * dim * dim * M + dim * dim * dim * K + dim * dim * L + dim * M + K ]
                        += 1. / ( higherOrderStressNorm[ K ] + tol );
                    for ( unsigned int N = 0; N < dim; N++ ){
                        for ( unsigned int O = 0; O < dim; O++ ){
                            for ( unsigned int P = 0; P < dim; P++ ){
                                for ( unsigned int Q = 0; Q < dim; Q++ ){
                                    d2HigherOrderStressNormdHigherOrderStress2[ tot_dim * tot_dim * K + dim * dim * dim * dim * dim * L + dim * dim * dim * dim * M + dim * dim * dim * N + dim * dim * O + dim * P + Q ]
                                        -=   dHigherOrderStressNormdHigherOrderStress[ tot_dim * K + dim * dim * L + dim * M + N ]
                                         *   dHigherOrderStressNormdHigherOrderStress[ tot_dim * K + dim * dim * O + dim * P + Q ]
                                         / ( higherOrderStressNorm[ K ] + tol );
                                }
                            }
                        }
                    }
                }
            }
        }

        return;
    }

    void assembleDeformationGradient( const variableVector &displacementGradient, variableVector &deformationGradient ){
        /*!
         * Assemble the deformation gradient from the gradient of the displacement.
         *
         * \param &displacementGradient: The gradient of the displacement.
         * \param &deformationGradient: The deformation gradient.
         */

        //Assume 3D
        constexpr unsigned int dim = 3;
        constexpr unsigned int sot_dim = dim * dim;

        TARDIGRADE_ERROR_TOOLS_CHECK( displacementGradient.size() == sot_dim, "The gradient of the deformation is not 3D" );

        deformationGradient = displacementGradient;
        for ( unsigned int i = 0; i < dim; i++ ){ deformationGradient[ dim * i + i ] += 1; }

        return;
    }

    void assembleDeformationGradient( const variableVector &displacementGradient, variableVector &deformationGradient,
                                          variableVector &dFdGradU ){
        /*!
         * Assemble the deformation gradient from the gradient of the displacement.
         *
         * \param &displacementGradient: The gradient of the displacement.
         * \param &deformationGradient: The deformation gradient.
         * \param &dFdGradU: The Jacobian of the deformation gradient w.r.t. the displacement gradient.
         */

        //Assume 3D
        constexpr unsigned int dim = 3;
        constexpr unsigned int sot_dim = dim * dim;

        TARDIGRADE_ERROR_TOOLS_CATCH( assembleDeformationGradient( displacementGradient, deformationGradient ) );

        dFdGradU = variableVector( sot_dim * sot_dim, 0 );
        for ( unsigned int i = 0; i < sot_dim; i++ ){ dFdGradU[ sot_dim * i + i ] = 1; }

        return;
    }

    void assembleMicroDeformation( const variableVector &microDisplacement, variableVector &microDeformation ){
        /*!
         * Assemble the micro deformation from the micro displacement
         *
         * \param &microDisplacement: The micro degrees of freedom.
         * \param &microDeformation: The micro deformation.
         */

        //Assume 3D
        constexpr unsigned int dim = 3;
        constexpr unsigned int sot_dim = dim * dim;

        TARDIGRADE_ERROR_TOOLS_CHECK( microDisplacement.size() == sot_dim, "The micro degrees of freedom must be 3D" );

        microDeformation = microDisplacement;
        for ( unsigned int i = 0; i < dim; i++ ){ microDeformation[ dim * i + i ] += 1; }

        return;
    }

    void assembleMicroDeformation( const variableVector &microDisplacement, variableVector &microDeformation,
                                       variableVector &dChidPhi ){
        /*!
         * Assemble the micro deformation from the micro displacement
         *
         * \param &microDisplacement: The micro degrees of freedom.
         * \param &microDeformation: The micro deformation.
         * \param &dChidPhi: The Jacobian of the micro deformation w.r.t. the micro displacement
         */

        //Assume 3D
        constexpr unsigned int dim = 3;
        constexpr unsigned int sot_dim = dim * dim;

        TARDIGRADE_ERROR_TOOLS_CATCH( assembleMicroDeformation( microDisplacement, microDeformation ) );

        dChidPhi = variableVector( sot_dim * sot_dim, 0 );
        for ( unsigned int i = 0; i < sot_dim; i++ ){ dChidPhi[ sot_dim * i + i ] = 1; }

        return;
    }

    void assembleGradientMicroDeformation( const variableVector &gradientMicroDisplacement,
                                               variableVector &gradientMicroDeformation ){
        /*!
         * Assemble the gradient of the micro deformation from the gradient of the micro displacement
         * in the reference configuration
         *
         * \param &gradientMicroDisplacement: The gradient of the micro displacement
         * \param &gradientMicroDeformation: The gradient of the micro deformation.
         */

        //Assume 3D
        constexpr unsigned int dim = 3;
        constexpr unsigned int sot_dim = dim * dim;

        TARDIGRADE_ERROR_TOOLS_CHECK( gradientMicroDisplacement.size() == sot_dim * dim, "The gradient of the micro displacement must be 3D" );

        gradientMicroDeformation = gradientMicroDisplacement;

        return;
    } 

    void assembleGradientMicroDeformation( const variableVector &gradientMicroDisplacement,
                                               variableVector &gradientMicroDeformation,
                                               variableVector &dGradChidGradPhi ){
        /*!
         * Assemble the gradient of the micro deformation from the gradient of the micro displacement
         * in the reference configuration
         *
         * \param &gradientMicroDisplacement: The gradient of the micro displacement
         * \param &gradientMicroDeformation: The gradient of the micro deformation.
         * \param &dGradChidGradPhi: The gradient of the micro deformation gradient w.r.t.
         *     the micro displacement gradient.
         */

        //Assume 3D
        constexpr unsigned int dim = 3;
        constexpr unsigned int sot_dim = dim * dim;
        constexpr unsigned int tot_dim = sot_dim * dim;

        TARDIGRADE_ERROR_TOOLS_CATCH( assembleGradientMicroDeformation( gradientMicroDisplacement, gradientMicroDeformation ) );
        
        dGradChidGradPhi = variableVector( tot_dim * tot_dim, 0 );
        for ( unsigned int i = 0; i < tot_dim; i++ ){ dGradChidGradPhi[ tot_dim * i + i ] = 1; }

        return;
    }

    void dCauchyStressdPK2Stress( const variableVector &deformationGradient,
                                  variableVector &dCauchyStressdPK2Stress ){
        /*!
         * Compute the derivative of the Cauchy stress w.r.t. the PK2 stress
         * 
         * \param &deformationGradient: The deformation gradient from the reference configuration to the current configuration
         * \param &dCauchyStressdPK2Stress: The derivative of the Cauchy stress w.r.t. the PK2 stress
         */

        //Assume 3d
        constexpr unsigned int dim = 3;
        constexpr unsigned int sot_dim = dim * dim;

        variableType detF;

        Eigen::Map< const Eigen::Matrix< variableType, dim, dim, Eigen::RowMajor > > map( deformationGradient.data( ), dim, dim );
        detF = map.determinant( );

        //Assemble the derivative
        dCauchyStressdPK2Stress = variableVector( sot_dim * sot_dim, 0 );

        for ( unsigned int i = 0; i < dim; i++ ){
            for ( unsigned int j = 0; j < dim; j++ ){
                for ( unsigned int k = 0; k < dim; k++ ){
                    for ( unsigned int K = 0; K < dim; K++ ){
                        dCauchyStressdPK2Stress[ dim * sot_dim * i + sot_dim * j + dim * k + K ] = deformationGradient[ dim * i + k ]
                                                                                                 * deformationGradient[ dim * j + K ] / detF;
                    }
                }
            }
        }

    }

    void dCauchyStressdPK2Stress( const variableVector &deformationGradient,
                                  variableVector &dCauchyStressdPK2Stress,
                                  variableVector &dRdF ){
        /*!
         * Compute the derivative of the Cauchy stress w.r.t. the PK2 stress
         * 
         * \param &deformationGradient: The deformation gradient from the reference configuration to the current configuration
         * \param &dCauchyStressdPK2Stress: The derivative of the Cauchy stress w.r.t. the PK2 stress
         * \param &dRdF: The derivative of the result w.r.t. the deformation gradient
         */

        //Assume 3d
        constexpr unsigned int dim = 3;
        constexpr unsigned int sot_dim = dim * dim;

        variableType detF;

        variableVector inverseDeformationGradient = deformationGradient;
        Eigen::Map< const Eigen::Matrix< variableType, dim, dim, Eigen::RowMajor > > map1( deformationGradient.data( ), dim, dim );
        Eigen::Map< Eigen::Matrix< variableType, dim, dim, Eigen::RowMajor > > map2( inverseDeformationGradient.data( ), dim, dim );
        detF = map1.determinant( );
        map2 = map2.inverse( ).transpose( ).eval( );

        //Assemble the derivative
        dCauchyStressdPK2Stress = variableVector( sot_dim * sot_dim, 0 );
        dRdF = variableVector( sot_dim * sot_dim * sot_dim, 0 );

        for ( unsigned int i = 0; i < dim; i++ ){
            for ( unsigned int j = 0; j < dim; j++ ){
                for ( unsigned int k = 0; k < dim; k++ ){
                    for ( unsigned int K = 0; K < dim; K++ ){
                        dCauchyStressdPK2Stress[ dim * sot_dim * i + sot_dim * j + dim * k + K ] += deformationGradient[ dim * i + k ]
                                                                                                  * deformationGradient[ dim * j + K ] / detF;

                        dRdF[ dim * sot_dim * sot_dim * i + sot_dim * sot_dim * j + dim * sot_dim * k + sot_dim * K + dim * i + k ] += deformationGradient[ dim * j + K ] / detF;

                        dRdF[ dim * sot_dim * sot_dim * i + sot_dim * sot_dim * j + dim * sot_dim * k + sot_dim * K + dim * j + K ] += deformationGradient[ dim * i + k ] / detF;

                        for ( unsigned int aA = 0; aA < sot_dim; ++aA ){
                            dRdF[ dim * sot_dim * sot_dim * i + sot_dim * sot_dim * j + dim * sot_dim * k + sot_dim * K + aA ] +=
                                                                                                 - deformationGradient[ dim * i + k ]
                                                                                                 * deformationGradient[ dim * j + K ]
                                                                                                 * inverseDeformationGradient[ aA ] / detF;
                        }
                    }
                }
            }
        }

    }

    void dSymmetricMicroStressdReferenceSymmetricMicroStress( const variableVector &deformationGradient,
                                                              variableVector &dSymmetricMicroStressdReferenceSymmetricMicroStress ){
        /*!
         * Compute the derivative of the symmetric micro stress w.r.t. the reference symmetric micro stress
         * 
         * \param &deformationGradient: The deformation gradient from the reference configuration to the current configuration
         * \param &dSymmetricMicroStressdReferenceSymmetricMicroStress: The derivative of the symmetric micro stress w.r.t. the reference symmetric micro stress
         */

        dCauchyStressdPK2Stress( deformationGradient, dSymmetricMicroStressdReferenceSymmetricMicroStress );

    }

    void dSymmetricMicroStressdReferenceSymmetricMicroStress( const variableVector &deformationGradient,
                                                              variableVector &dSymmetricMicroStressdReferenceSymmetricMicroStress,
                                                              variableVector &dRdF ){
        /*!
         * Compute the derivative of the symmetric micro stress w.r.t. the reference symmetric micro stress
         * 
         * \param &deformationGradient: The deformation gradient from the reference configuration to the current configuration
         * \param &dSymmetricMicroStressdReferenceSymmetricMicroStress: The derivative of the symmetric micro stress w.r.t. the reference symmetric micro stress
         * \param &dRdF: The derivative of the result w.r.t. the deformation gradient
         */

        dCauchyStressdPK2Stress( deformationGradient, dSymmetricMicroStressdReferenceSymmetricMicroStress );

    }

    void dHigherOrderStressdReferenceHigherOrderStress( const variableVector &deformationGradient,
                                                        const variableVector &microDeformation,
                                                        variableVector &dHigherOrderStressdReferenceHigherOrderStress ){
        /*!
         * Compute the derivative of the higher order stress w.r.t. the reference higher order stress
         * 
         * \param &deformationGradient: The deformation gradient from the reference configuration to the current configuration
         * \param &microDeformation: The micro deformationt from the reference configuration to the current configuration
         * \param &dHigherOrderStressdReferenceHigherOrderStress: The derivative of the higher order stress w.r.t. the reference higher order stress
         */

        //Assume 3d
        constexpr unsigned int dim = 3;
        constexpr unsigned int tot_dim = dim * dim * dim;

        variableType detF;

        Eigen::Map< const Eigen::Matrix< variableType, dim, dim, Eigen::RowMajor > > map( deformationGradient.data( ), dim, dim );
        detF = map.determinant( );

        dHigherOrderStressdReferenceHigherOrderStress = variableVector( tot_dim * tot_dim, 0 );

        for ( unsigned int i = 0; i < dim; i++ ){
            for ( unsigned int j = 0; j < dim; j++ ){
                for ( unsigned int k = 0; k < dim; k++ ){
                    for ( unsigned int l = 0; l < dim; l++ ){
                        for ( unsigned int M = 0; M < dim; M++ ){
                            for ( unsigned int N = 0; N < dim; N++ ){
                                dHigherOrderStressdReferenceHigherOrderStress[ dim * dim * tot_dim * i + dim * tot_dim * j + tot_dim * k + dim * dim * l + dim * M + N ] = deformationGradient[ dim * i + l ] * deformationGradient[ dim * j + M ] * microDeformation[ dim * k + N] / detF;

                            }
                        }
                    }
                }
            }
        }
    }

    void dHigherOrderStressdReferenceHigherOrderStress( const variableVector &deformationGradient,
                                                        const variableVector &microDeformation,
                                                        variableVector &dHigherOrderStressdReferenceHigherOrderStress,
                                                        variableVector &dRdF,
                                                        variableVector &dRdChi ){
        /*!
         * Compute the derivative of the higher order stress w.r.t. the reference higher order stress
         * 
         * \param &deformationGradient: The deformation gradient from the reference configuration to the current configuration
         * \param &microDeformation: The micro deformationt from the reference configuration to the current configuration
         * \param &dHigherOrderStressdReferenceHigherOrderStress: The derivative of the higher order stress w.r.t. the reference higher order stress
         * \param &dRdF: The derivative of the result w.r.t. the deformation gradient
         * \param &dRdF: The derivative of the result w.r.t. the micro deformation
         */

        //Assume 3d
        constexpr unsigned int dim = 3;
        constexpr unsigned int sot_dim = dim * dim;
        constexpr unsigned int tot_dim = dim * dim * dim;

        variableType detF;

        variableVector inverseDeformationGradient = deformationGradient;
        Eigen::Map< const Eigen::Matrix< variableType, dim, dim, Eigen::RowMajor > > map1( deformationGradient.data( ), dim, dim );
        Eigen::Map< Eigen::Matrix< variableType, dim, dim, Eigen::RowMajor > > map2( inverseDeformationGradient.data( ), dim, dim );
        detF = map1.determinant( );
        map2 = map2.inverse( ).transpose( ).eval( );

        dHigherOrderStressdReferenceHigherOrderStress = variableVector( tot_dim * tot_dim, 0 );
        dRdF =  variableVector( tot_dim * tot_dim * sot_dim, 0 );
        dRdChi = variableVector( tot_dim * tot_dim * sot_dim, 0 );

        for ( unsigned int i = 0; i < dim; i++ ){
            for ( unsigned int j = 0; j < dim; j++ ){
                for ( unsigned int k = 0; k < dim; k++ ){
                    for ( unsigned int l = 0; l < dim; l++ ){
                        for ( unsigned int M = 0; M < dim; M++ ){
                            for ( unsigned int N = 0; N < dim; N++ ){
                                dHigherOrderStressdReferenceHigherOrderStress[ dim * dim * tot_dim * i + dim * tot_dim * j + tot_dim * k + dim * dim * l + dim * M + N ] += deformationGradient[ dim * i + l ] * deformationGradient[ dim * j + M ] * microDeformation[ dim * k + N] / detF;

                                dRdF[ dim * dim * tot_dim * sot_dim * i + dim * tot_dim * sot_dim * j + tot_dim * sot_dim * k + dim * dim * sot_dim * l + dim * sot_dim * M + sot_dim * N + dim * i + l ] += deformationGradient[ dim * j + M ] * microDeformation[ dim * k + N] / detF;

                                dRdF[ dim * dim * tot_dim * sot_dim * i + dim * tot_dim * sot_dim * j + tot_dim * sot_dim * k + dim * dim * sot_dim * l + dim * sot_dim * M + sot_dim * N + dim * j + M ] += deformationGradient[ dim * i + l ] * microDeformation[ dim * k + N] / detF;

                                dRdChi[ dim * dim * tot_dim * sot_dim * i + dim * tot_dim * sot_dim * j + tot_dim * sot_dim * k + dim * dim * sot_dim * l + dim * sot_dim * M + sot_dim * N + dim * k + N ] += deformationGradient[ dim * i + l ] * deformationGradient[ dim * j + M ] / detF;

                                for ( unsigned int aA = 0; aA < sot_dim; ++aA ){

                                    dRdF[ dim * dim * tot_dim * sot_dim * i + dim * tot_dim * sot_dim * j + tot_dim * sot_dim * k + dim * dim * sot_dim * l + dim * sot_dim * M + sot_dim * N + aA ] -= deformationGradient[ dim * i + l ] * deformationGradient[ dim * j + M ] * microDeformation[ dim * k + N]  * inverseDeformationGradient[ aA ] / detF;
                                }
                            }
                        }
                    }
                }
            }
        }
    }

}
