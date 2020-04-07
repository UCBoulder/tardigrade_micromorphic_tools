/*!
 * Micromorphic Tools.cpp
 * ====================
 *
 * A collection of tools and utilities which seek to make the 
 * implementation of micromorphic continuum mechanics based constitutive 
 * theories easier.
 *
 */

#include<micromorphic_tools.h>

namespace micromorphicTools{

    errorOut computePsi( const variableVector &deformationGradient, const variableVector &microDeformation,
                         variableVector &Psi ){
        /*!
         * Computes the micromorphic quantity psi defined as:
         * \Psi_{IJ} = F_{i I} \Chi_{i J}
         *
         * :param const variableVector &deformationGradient: The deformation gradient.
         * :param const variableVector &microDeformation: The micro-deformation.
         * :param variableVector &Psi: The micro-displacement metric Psi.
         */

        //Assume 3d
        unsigned int dim = 3;

        if ( deformationGradient.size() != dim * dim ){
            return new errorNode( "computePsi", "The deformation gradient doesn't have the correct size" );
        }

        if ( microDeformation.size() != dim * dim ){
            return new errorNode( "computePsi", "The micro-deformation doesn't have the correct size" );
        }

        Psi = vectorTools::matrixMultiply( deformationGradient, microDeformation,
                                           dim, dim, dim, dim, 1, 0 );

        return NULL;
    }

    errorOut computePsi( const variableVector &deformationGradient, const variableVector &microDeformation,
                         variableVector &Psi, variableMatrix &dPsidF, variableMatrix &dPsidChi ){
        /*!
         * Computes the micromorphic quantity psi defined as:
         * \Psi_{IJ} = F_{i I} \Chi_{i J}
         *
         * along with the jacobians
         *
         * \frac{ \partial \Psi_{IJ} }{ \partial F_{k K} } = \delta_{I K} \Chi_{k J}
         * \frac{ \partial \Psi_{IJ} }{ \partial \Chi_{k K} } = F_{k I} \delta_{J K}
         *
         * :param const variableVector &deformationGradient: The deformation gradient.
         * :param const variableVector &microDeformation: The micro-deformation.
         * :param variableVector &Psi: The micro-displacement metric Psi
         * :param variableMatrix &dPsidF: The jacobian of Psi w.r.t. the deformation gradient.
         * :param variableMatrix &dPsidChi: The jacobian of Psi w.r.t. the micro-deformation.
         */

        //Assume 3d
        unsigned int dim = 3;

        errorOut error = computePsi( deformationGradient, microDeformation, Psi );

        if (error){
            errorOut result = new errorNode( "computePsi (jacobian)", "Error in the computation of Psi" );
            result->addNext( error );
            return result;
        }

        constantVector eye( dim * dim, 0 );
        vectorTools::eye( eye );

        dPsidF  = variableMatrix( Psi.size(), variableVector( deformationGradient.size(), 0 ) );
        dPsidChi = variableMatrix( Psi.size(), variableVector( microDeformation.size(), 0 ) );

        for ( unsigned int I = 0; I < dim; I++ ){
            for ( unsigned int J = 0; J < dim; J++ ){
                for ( unsigned int k = 0; k < dim; k++ ){
                    for ( unsigned int K = 0; K < dim; K++ ){
                        dPsidF[ dim * I + J ][ dim * k + K ] = eye[ dim * I + K ] * microDeformation[ dim * k + J ];
                        dPsidChi[ dim * I + J ][ dim * k + K ] = deformationGradient[ dim * k + I ] * eye[ dim * J + K ];
                    }
                }
            }
        }
        return NULL;
    }

    errorOut computeGamma( const variableVector &deformationGradient, const variableVector &gradChi,
                           variableVector &Gamma ){
        /*!
         * Compute the deformation metric Gamma:
         *
         * Gamma_{IJK} = F_{iI} \Chi_{iJ,K}
         *
         * :param const variableVector &deformationGradient: The deformation gradient.
         * :param const variableVector &gradChi: The gradient of the micro-deformation tensor
         *     w.r.t. the reference configuration.
         * :param variableVector &Gamma: The micromorphic deformation metric Gamma.
         */

        //Assume 3d
        unsigned int dim = 3;

        if ( deformationGradient.size() != dim * dim ){
            return new errorNode("computeGamma", "The deformation gradient isn't the right size");
        }

        if ( gradChi.size() != dim * dim * dim ){
            return new errorNode("computeGamma", "The micro-deformation gradient isn't the right size");
        }

        Gamma = variableVector( dim * dim * dim, 0 );

        for ( unsigned int I = 0; I < dim; I++ ){
            for ( unsigned int J = 0; J < dim; J++ ){
                for ( unsigned int K = 0; K < dim; K++ ){
                    for ( unsigned int i = 0; i < dim; i++ ){
                        Gamma[ dim * dim * I + dim * J + K ] += deformationGradient[ dim * i + I ] * gradChi[ dim * dim * i + dim * J + K ];
                    }
                }
            }
        }

        return NULL;
    }

    errorOut computeGamma( const variableVector &deformationGradient, const variableVector &gradChi,
                           variableVector &Gamma, variableMatrix &dGammadF, variableMatrix &dGammadGradChi ){
        /*!
         * Compute the deformation metric Gamma:
         *
         * Gamma_{IJK} = F_{iI} \Chi_{iJ,K}
         *
         * Also return the Jacobians
         * \frac{ \partial Gamma_{IJK} }{ \partial F_{lL} } = \delta_{IL} \Chi_{lJ,K}
         * \frac{ \partial Gamma_{IJK} }{ \partial \Chi_{lL,M} } = F_{lI} \delta_{JL} \delta_{KM}
         *
         * :param const variableVector &deformationGradient: The deformation gradient.
         * :param const variableVector &gradChi: The gradient of the micro-deformation tensor
         *     w.r.t. the reference configuration.
         * :param variableVector &Gamma: The micromorphic deformation metric Gamma.
         * :param variableMatrix &dGammadF: The gradient of Gamma w.r.t. the deformation gradient.
         * :param variableMatrix &dGammadGradChi: The gradient of Gamma w.r.t. the gradient of Chi in the reference 
         *     configuration.
         */

        //Assume 3d
        unsigned int dim = 3;

        errorOut error = computeGamma( deformationGradient, gradChi, Gamma );

        if ( error ){
            errorOut result = new errorNode("computeGamma (jacobian)", "Error in computation of Gamma");
            result->addNext(error);
            return result;
        }

        dGammadF      = variableMatrix( dim * dim * dim, variableVector( dim * dim, 0 ) );
        dGammadGradChi = variableMatrix( dim * dim * dim, variableVector( dim * dim * dim, 0 ) );

        constantVector eye( dim * dim, 0 );
        vectorTools::eye( eye );

        for ( unsigned int I = 0; I < dim; I++ ){
            for ( unsigned int J = 0; J < dim; J++ ){
                for ( unsigned int K = 0; K < dim; K++ ){
                    for ( unsigned int l = 0; l < dim; l++ ){
                        for ( unsigned int L = 0; L < dim; L++ ){
                            dGammadF[ dim * dim * I + dim * J + K ][ dim * l + L ] = eye[ dim * I + L ]
                                                                                   * gradChi[ dim * dim * l + dim * J + K ];
                            for ( unsigned int M = 0; M < dim; M++ ){
                                dGammadGradChi[ dim * dim * I + dim * J + K ][ dim * dim * l + dim * L + M ] = deformationGradient[ dim * l + I ]
                                                                                                            * eye[ dim * J + L ] 
                                                                                                            * eye[ dim * K + M ];
                            }
                        }
                    }
                }
            }
        }

        return NULL;
    }

    errorOut computeMicroStrain( const variableVector &Psi, variableVector &microStrain ){
        /*!
         * Compute the microstrain defined as:
         * \Epsilon_{IJ} = \Psi_{IJ} - eye_{IJ}
         *
         * :param const variableVector &Psi: The micro-deformation metric Psi.
         * :param variableVector &microStrain: The micro-strain.
         */

        //Assume 3d
        unsigned int dim = 3;

        if ( Psi.size() != dim * dim ){
            return new errorNode( "computeMicroStrain", "Psi is not of the correct size" );
        }

        constantVector eye( dim * dim, 0 );
        vectorTools::eye( eye );

        microStrain = Psi - eye;

        return NULL;
    }

    errorOut computeMicroStrain( const variableVector &Psi, variableVector &microStrain,
                                 variableMatrix &dMicroStraindPsi ){
        /*!
         * Compute the microstrain defined as:
         * \Epsilon_{IJ} = \Psi_{IJ} - eye_{IJ}
         *
         * and also compute the jacobian
         * \frac{ \partial \Epsilon_{IJ} }{ \partial \Psi_{KL} } = \delta_{IK} \delta_{JL}
         *
         * :param const variableVector &Psi: The micro-deformation metric Psi.
         * :param variableVector &microStrain: The micro-strain.
         * :param variableMatrix &dMicroStraindPsi: The jacobian of the micro-strain.
         */

        errorOut error = computeMicroStrain( Psi, microStrain );

        if (error){
            errorOut result = new errorNode( "computeMicroStrain (jacobian)", "Error in computation of micro-strain" );
            result->addNext( error );
            return result;
        }

        dMicroStraindPsi = vectorTools::eye< variableType >( Psi.size() );

        return NULL;
    }

    errorOut pushForwardPK2Stress( const variableVector &PK2Stress,
                                   const variableVector &deformationGradient,
                                   variableVector &cauchyStress ){
        /*!
         * Push forward the PK2-stress in the reference configuration to the 
         * configuration (the Cauchy stress) indicated by the deformation gradient.
         *
         * \cauchy_{ij} = (1 / J ) F_{i I} S_{I J} F_{j J}
         *
         * :param const variableVector &PK2Stress: The PK2 stress in the 
         *     reference configuration.
         * :param const variableVector &deformationGradient: The deformation gradient 
         *     mapping between configurations.
         * :param variableVector &cauchyStress: The Cauchy stress in the current 
         *     configuration.
         */

        variableType detF;
        errorOut error = pushForwardReferenceMicroStress( PK2Stress, deformationGradient, 
                                                          detF, cauchyStress );

        if ( error ){
            errorOut result = new errorNode( "pushForwardPK2Stress",
                                             "Error in push-forward operation (micro-stress and PK2 are identical)" );
            result->addNext( error );
            return result;
        }
        return NULL;
    }

    errorOut pushForwardPK2Stress( const variableVector &PK2Stress,
                                   const variableVector &deformationGradient,
                                   variableVector &cauchyStress,
                                   variableMatrix &dCauchyStressdPK2Stress,
                                   variableMatrix &dCauchyStressdDeformationGradient ){
        /*!
         * Push forward the PK2 stress in the reference configuration to the 
         * configuration indicated by the deformation gradient.
         *
         * \sigma_{ij} = (1 / J ) F_{iI} S_{IJ} F_{jJ}
         *
         * Also computes the jacobians:
         * \frac{ \partial \cauchy_{ij} }{\partial \Sigma_{KL} } = ( 1 / J ) F_{iK} F_{jL}
         * \frac{ \partial \cauchy_{ij} }{\partial F_{kK} } = ( \delta_{i k} \delta_{I K} S_{I J} F_{j J}
         *                                                  + F_{i I} S_{I J} \delta_{j k} \delta_{J K}
         *                                                  - \cauchy_{i j} dDetFdF_{kK} ) / J
         *
         * :param const variableVector &referenceMicroStress: The PK2 stress in the 
         *     reference configuration.
         * :param const variableVector &deformationGradient: The deformation gradient 
         *     mapping between configurations.
         * :param variableVector &cauchyStress: The Cauchy stress in the current 
         *     configuration.
         * :param variableMatrix &dCauchyStressdReferenceMicroStress: The jacobian of 
         *     the Cauchy w.r.t. the PK2 tress in the reference configuration.
         * :param variableMatrix &dCauchyStressdDeformationGradient: The jacobian of 
         *     the Cauchy stress w.r.t. the deformation gradient.
         */

        errorOut error = pushForwardReferenceMicroStress( PK2Stress, deformationGradient, cauchyStress,
                                                          dCauchyStressdPK2Stress,
                                                          dCauchyStressdDeformationGradient );

        if ( error ){
            errorOut result = new errorNode( "pushForwardPK2Stress (jacobian)",
                                             "Error in push-forward operation (micro-stress and PK2 are identical)" );
            result->addNext( error );
            return result;
        }
        return NULL;
    }

    errorOut pullBackCauchyStress( const variableVector &cauchyStress,
                                   const variableVector &deformationGradient,
                                   variableVector &PK2Stress ){
        /*!
         * Pull back the Cauchy stress in the configuration indicated by the deformation gradient
         * to the PK2 stress.
         *
         * S_{IJ} = J F_{Ii}^{-1} \sigma_{ij} F_{Jj}^{-1}
         *
         * :param const variableVector &cauchyStress: The Cauchy stress in the current configuration of
         *     the provided deformation gradient.
         * :param const variableVector &deformationGradient: The deformation gradient mapping between the 
         *     reference configuration and the current configuration.
         * :param variableVector &PK2Stress: The PK2 stress in the reference configuration.
         */

        errorOut error = pullBackMicroStress( cauchyStress, deformationGradient, PK2Stress );

        if ( error ){
            errorOut result = new errorNode( "pullBackCauchyStress",
                                             "Error in pull-back operation (micro-stress and Cauchy stress are identical)" );
            result->addNext( error );
            return result;
        }
        return NULL;
    }

    errorOut pullBackCauchyStress( const variableVector &cauchyStress,
                                   const variableVector &deformationGradient,
                                   variableVector &PK2Stress, variableMatrix &dPK2StressdCauchyStress,
                                   variableMatrix &dPK2StressdDeformationGradient ){
        /*!
         * Pull back the Cauchy stress in the configuration indicated by the deformation gradient
         * to the PK2 stress.
         *
         * S_{IJ} = J F_{Ii}^{-1} \sigma_{ij} F_{Jj}^{-1}
         *
         * Also computes the Jacobians
         *
         * \frac{ \partial S_{IJ} }{ \partial \sigma_{kl} } = J F_{Ik}^{-1} F_{Jl}^{-1}
         * \frac{ \partial S_{IJ} }{ \partial F_{kK} } = F_{Kk}^{-1} S_{IJ} - F_{Ik}^{-1} S_{KJ} - S_{IK} F_{Jk}^{-1}
         *
         * :param const variableVector &cauchyStress: The Cauchy stress in the current configuration of
         *     the provided deformation gradient.
         * :param const variableVector &deformationGradient: The deformation gradient mapping between the 
         *     reference configuration and the current configuration.
         * :param variableVector &PK2Stress: The PK2 stress in the reference configuration.
         */

        errorOut error = pullBackMicroStress( cauchyStress, deformationGradient, PK2Stress,
                                              dPK2StressdCauchyStress, dPK2StressdDeformationGradient );

        if ( error ){
            errorOut result = new errorNode( "pullBackCauchyStress (jacobian)",
                                             "Error in pull-back operation (micro-stress and Cauchy stress are identical)" );
            result->addNext( error );
            return result;
        }
        return NULL;
    }

    errorOut pushForwardReferenceMicroStress( const variableVector &referenceMicroStress,
                                              const variableVector &deformationGradient,
                                              variableVector &microStress ){
        /*!
         * Push forward the micro-stress in the reference configuration to the 
         * configuration indicated by the deformation gradient.
         *
         * s_{ij} = (1 / J ) F_{i I} \Sigma_{I J} F_{j J}
         *
         * :param const variableVector &referenceMicroStress: The micro-stress in the 
         *     reference configuration.
         * :param const variableVector &deformationGradient: The deformation gradient 
         *     mapping between configurations.
         * :param variableVector &microStress: The micro-stress in the current 
         *     configuration.
         */

        variableType detF;
        return pushForwardReferenceMicroStress( referenceMicroStress, deformationGradient, 
                                                detF, microStress );

    }

    errorOut pushForwardReferenceMicroStress( const variableVector &referenceMicroStress,
                                              const variableVector &deformationGradient,
                                              variableType &detF, variableVector &microStress ){
        /*!
         * Push forward the micro-stress in the reference configuration to the 
         * configuration indicated by the deformation gradient.
         *
         * s_{ij} = (1 / J ) F_{i I} \Sigma_{I J} F_{j J}
         *
         * :param const variableVector &referenceMicroStress: The micro-stress in the 
         *     reference configuration.
         * :param const variableVector &deformationGradient: The deformation gradient 
         *     mapping between configurations.
         * :param variableType &detF: The determinant of the deformation gradient.
         * :param variableVector &microStress: The micro-stress in the current 
         *     configuration.
         */

        //Assume 3d
        unsigned int dim = 3;

        if ( referenceMicroStress.size() != dim * dim ){
            return new errorNode("pushForwardReferenceMicroStress", "The reference micro-stress has an incorrect size");
        }

        if ( deformationGradient.size() != dim * dim ){
            return new errorNode("pushForwardReferenceMicroStress", "The deformation gradient has an incorrect size");
        }

        microStress = variableVector( dim * dim, 0 );

        detF = vectorTools::determinant( deformationGradient, dim, dim );

        for ( unsigned int i = 0; i < dim; i++ ){
            for ( unsigned int j = 0; j < dim; j++ ){
                for ( unsigned int I = 0; I < dim; I++ ){
                    for ( unsigned int J = 0; J < dim; J++ ){
                        microStress[ dim * i + j ] += deformationGradient[ dim * i + I ]
                                                    * referenceMicroStress[ dim * I + J ]
                                                    * deformationGradient[ dim * j + J ];
                    }
                }
                microStress[ dim * i + j ] /= detF;
            }
        }

        return NULL;
    }

    errorOut pushForwardReferenceMicroStress( const variableVector &referenceMicroStress,
                                              const variableVector &deformationGradient,
                                              variableVector &microStress,
                                              variableMatrix &dMicroStressdReferenceMicroStress,
                                              variableMatrix &dMicroStressdDeformationGradient ){
        /*!
         * Push forward the micro-stress in the reference configuration to the 
         * configuration indicated by the deformation gradient.
         *
         * s_{ij} = (1 / J ) F_{iI} \Sigma_{IJ} F_{jJ}
         *
         * Also computes the jacobians:
         * \frac{ \partial s_{ij} }{\partial \Sigma_{KL} } = ( 1 / J ) F_{iK} F_{jL}
         * \frac{ \partial s_{ij} }{\partial F_{kK} } = ( \delta_{i k} \delta_{I K} \Sigma_{I J} F_{j J}
         *                                              + F_{i I} \Sigma_{I J} \delta_{j k} \delta_{J K}
         *                                              - s_{i j} dDetFdF_{kK} ) / J
         *
         * :param const variableVector &referenceMicroStress: The micro-stress in the 
         *     reference configuration.
         * :param const variableVector &deformationGradient: The deformation gradient 
         *     mapping between configurations.
         * :param variableVector &microStress: The micro-stress in the current 
         *     configuration.
         * :param variableMatrix &dmicroStressdReferenceMicroStress: The jacobian of 
         *     the micro-stress w.r.t. the micro-stress in the reference configuration.
         * :param variableMatrix &dmicroSTressdDeformationGradient: The jacobian of 
         *     the micro-stress w.r.t. the deformation gradient.
         */

        //Assume 3d
        unsigned int dim = 3;

        variableType detF;
        errorOut error = pushForwardReferenceMicroStress( referenceMicroStress, deformationGradient,
                                                          detF, microStress );

        if (error){
            errorOut result = new errorNode( "pushForwardReferenceMicroStress (jacobian)", "Error in computation of push forward of micro-stress" );
            result->addNext(error);
            return result;
        }

        //Assemble the jacobian of the determinant of the deformation gradient
        variableVector inverseDeformationGradient = vectorTools::inverse( deformationGradient, dim, dim );

        variableVector dDetFdF( dim * dim, 0 );

        for (unsigned int i = 0; i < dim; i++ ){
            for (unsigned int I = 0; I < dim; I++ ){
                dDetFdF[ dim * i + I ] = inverseDeformationGradient[ dim * I + i ] * detF;
            }
        }

        //Assemble the jacobians
        dMicroStressdReferenceMicroStress = variableMatrix( microStress.size(), variableVector( referenceMicroStress.size(), 0 ) );
        dMicroStressdDeformationGradient = variableMatrix( microStress.size(), variableVector( deformationGradient.size(), 0 ) );

        constantVector eye( dim * dim );
        vectorTools::eye( eye );

        for ( unsigned int i = 0; i < dim; i++ ){
            for ( unsigned int j = 0; j < dim; j++ ){
                for ( unsigned int k = 0; k < dim; k++ ){
                    for ( unsigned int K = 0; K < dim; K++ ){
                        dMicroStressdReferenceMicroStress[ dim * i + j ][ dim * k + K ] = deformationGradient[ dim * i + k ]
                                                                                        * deformationGradient[ dim * j + K ] / detF;
                        
                        for ( unsigned int I = 0; I < dim; I++ ){

                            dMicroStressdDeformationGradient[ dim * i + j ][ dim * k + K] += eye[ dim * i + k ] * referenceMicroStress[ dim * K + I ] * deformationGradient[ dim * j + I ]
                                                                                          + deformationGradient[ dim * i + I ] * referenceMicroStress[ dim * I + K ] * eye[ dim * j + k ];
                        }

                        dMicroStressdDeformationGradient[ dim * i + j ][ dim * k + K] -= microStress[ dim * i + j ] * dDetFdF[ dim * k + K ];
                        dMicroStressdDeformationGradient[ dim * i + j ][ dim * k + K] /= detF; 
                    }
                }
            }
        }

        return NULL;
    }

    errorOut pullBackMicroStress( const variableVector &microStress,
                                  const variableVector &deformationGradient,
                                  variableVector &referenceMicroStress ){
        /*!
         * Push forward the micro-stress in the reference configuration to the 
         * configuration indicated by the deformation gradient.
         *
         * \Sigma_{IJ} = J F_{I i}^{-1} s_{ij} F_{J j}^{-1}
         *
         * :param const variableVector &microStress: The micro-stress in the current 
         *     configuration.
         * :param const variableVector &deformationGradient: The deformation gradient 
         *     mapping between configurations.
         * :param variableVector &referenceMicroStress: The micro-stress in the 
         *     reference configuration.
         */

        variableType detF;
        variableVector inverseDeformationGradient;
        return pullBackMicroStress( microStress, deformationGradient, 
                                    detF, inverseDeformationGradient, referenceMicroStress );

    }

    errorOut pullBackMicroStress( const variableVector &microStress,
                                  const variableVector &deformationGradient,
                                  variableType &detF, variableVector &inverseDeformationGradient,
                                  variableVector &referenceMicroStress ){
        /*!
         * Push forward the micro-stress in the reference configuration to the 
         * configuration indicated by the deformation gradient.
         *
         * \Sigma_{IJ} = J F_{I i}^{-1} s_{ij} F_{J j}^{-1}
         *
         * :param const variableVector &microStress: The micro-stress in the current 
         *     configuration.
         * :param const variableVector &deformationGradient: The deformation gradient 
         *     mapping between configurations.
         * :param variableType &detF: The determinant of the deformation gradient.
         * :param variableVector &inverseDeformationGradient: The inverse of the deformation gradient.
         * :param variableVector &referenceMicroStress: The micro-stress in the 
         *     reference configuration.
         */

        //Assume 3d
        unsigned int dim = 3;

        if ( microStress.size() != dim * dim ){
            return new errorNode("pullBackdMicroStress", "The micro-stress has an incorrect size");
        }

        if ( deformationGradient.size() != dim * dim ){
            return new errorNode("pullBackdMicroStress", "The deformation gradient has an incorrect size");
        }

        referenceMicroStress = variableVector( dim * dim, 0 );

        detF = vectorTools::determinant( deformationGradient, dim, dim );
        inverseDeformationGradient = vectorTools::inverse( deformationGradient, dim, dim );

        for ( unsigned int I = 0; I < dim; I++ ){
            for ( unsigned int J = 0; J < dim; J++ ){
                for ( unsigned int i = 0; i < dim; i++ ){
                    for ( unsigned int j = 0; j < dim; j++ ){
                        referenceMicroStress[ dim * I + J ] += inverseDeformationGradient[ dim * I + i ]
                                                             * microStress[ dim * i + j ]
                                                             * inverseDeformationGradient[ dim * J + j ];
                    }
                }
                referenceMicroStress[ dim * I + J ] *= detF;
            }
        }

        return NULL;
    }

    errorOut pullBackMicroStress( const variableVector &microStress,
                                  const variableVector &deformationGradient,
                                  variableVector &referenceMicroStress,
                                  variableMatrix &dReferenceMicroStressdMicroStress,
                                  variableMatrix &dReferenceMicroStressdDeformationGradient ){
        /*!
         * Push forward the micro-stress in the reference configuration to the 
         * configuration indicated by the deformation gradient.
         *
         * \Sigma_{IJ} = J F_{Ii}^{-1} s_{IJ} F_{Jj}^{-1}
         *
         * Also computes the jacobians:
         * \frac{ \partial \Sigma_{IJ} }{ \partial s_{kl} } = J F_{Ik}^{-1} F_{Jl}^{-1}
         * \frac{ \partial \Sigma_{IJ} }{ \partial F_{kK} } = F_{Kk}^{-1} \Sigma_{IJ} - F_{Ik}^{-1} \Sigma_{KJ} - \Sigma_{IK} F_{Jk}^{-1}
         *
         * :param const variableVector &microStress: The micro-stress in the current 
         *     configuration.
         * :param const variableVector &deformationGradient: The deformation gradient 
         *     mapping between configurations.
         * :param variableVector &referenceMicroStress: The micro-stress in the 
         *     reference configuration.
         * :param variableMatrix &dReferenceMicroStressdMicroStress: The jacobian of 
         *     the reference micro-stress w.r.t. the micro-stress in the reference configuration.
         * :param variableMatrix &dReferenceMicroStressdDeformationGradient: The jacobian of 
         *     the reference micro-stress w.r.t. the deformation gradient.
         */

        //Assume 3d
        unsigned int dim = 3;

        variableType detF;
        variableVector inverseDeformationGradient;
        errorOut error = pullBackMicroStress( microStress, deformationGradient,
                                              detF, inverseDeformationGradient, referenceMicroStress );

        if (error){
            errorOut result = new errorNode( "pullBackMicroStress (jacobian)", "Error in computation of pull back of micro-stress" );
            result->addNext(error);
            return result;
        }

        //Assemble the jacobian of the determinant of the deformation gradient
        variableVector dDetFdF( dim * dim, 0 );

        for (unsigned int i = 0; i < dim; i++ ){
            for (unsigned int I = 0; I < dim; I++ ){
                dDetFdF[ dim * i + I ] = inverseDeformationGradient[ dim * I + i ] * detF;
            }
        }

        //Assemble the jacobians
        dReferenceMicroStressdMicroStress = variableMatrix( referenceMicroStress.size(), variableVector( microStress.size(), 0 ) );
        dReferenceMicroStressdDeformationGradient = variableMatrix( referenceMicroStress.size(), variableVector( deformationGradient.size(), 0 ) );

        constantVector eye( dim * dim );
        vectorTools::eye( eye );

        for ( unsigned int I = 0; I < dim; I++ ){
            for ( unsigned int J = 0; J < dim; J++ ){
                for ( unsigned int k = 0; k < dim; k++ ){
                    for ( unsigned int K = 0; K < dim; K++ ){
                        dReferenceMicroStressdMicroStress[ dim * I + J ][ dim * k + K ]
                            = detF * inverseDeformationGradient[ dim * I + k ] * inverseDeformationGradient[ dim * J + K ];
                        dReferenceMicroStressdDeformationGradient[ dim * I + J ][ dim * k + K ]
                            = inverseDeformationGradient[ dim * K + k ] * referenceMicroStress[ dim * I + J ]
                            - inverseDeformationGradient[ dim * I + k ] * referenceMicroStress[ dim * K + J ]
                            - inverseDeformationGradient[ dim * J + k ] * referenceMicroStress[ dim * I + K ];
                    }
                }
            }
        }

        return NULL;
    }

    errorOut pushForwardHigherOrderStress( const variableVector &referenceHigherOrderStress,
                                           const variableVector &deformationGradient,
                                           const variableVector &microDeformation,
                                           variableVector &higherOrderStress ){
        /*!
         * Compute the push-forward operation on the higher order stress.
         *
         * m_{ijk} = \frac{1}{J} F_{iI} F_{jJ} \Chi_{kK} M_{IJK}
         *
         * :param const variableVector &referenceHigherOrderStress: The higher order stress in the 
         *     reference configuration.
         * :param const variableVector &deformationGradient: The deformation gradient which maps 
         *     between the reference and current configurations.
         * :param const variableVector &microDeformation: The micro-deformation tensor.
         * :param variableVector &higherOrderStress: The higher order stress in the current configuration.
         */

        variableType detF;
        return pushForwardHigherOrderStress( referenceHigherOrderStress, deformationGradient,
                                             microDeformation, detF, higherOrderStress );
    }

    errorOut pushForwardHigherOrderStress( const variableVector &referenceHigherOrderStress,
                                           const variableVector &deformationGradient,
                                           const variableVector &microDeformation,
                                           variableType &detF,
                                           variableVector &higherOrderStress ){
        /*!
         * Compute the push-forward operation on the higher order stress.
         *
         * m_{ijk} = \frac{1}{J} F_{iI} F_{jJ} \Chi_{kK} M_{IJK}
         *
         * :param const variableVector &referenceHigherOrderStress: The higher order stress in the 
         *     reference configuration.
         * :param const variableVector &deformationGradient: The deformation gradient which maps 
         *     between the reference and current configurations.
         * :param const variableVector &microDeformation: The micro-deformation tensor.
         * :param variableType &detF: The determinant of the deformation gradient.
         * :param variableVector &higherOrderStress: The higher order stress in the current configuration.
         */

        //Assume 3d
        unsigned int dim = 3;

        if ( referenceHigherOrderStress.size() != dim * dim * dim ){
            return new errorNode( "pushForwardHigherOrderStress", "The reference higher order stress doesn't have the correct size" );
        }

        if ( deformationGradient.size() != dim * dim ){
            return new errorNode( "pushForwardHigherOrderStress", "The deformation gradient doesn't have the correct size" );
        }

        if ( microDeformation.size() != dim * dim ){
            return new errorNode( "pushForwardHigherOrderStress", "The micro-deformation doesn't have the correct size" );
        }

        detF = vectorTools::determinant( deformationGradient, dim, dim );

        higherOrderStress = variableVector( dim * dim * dim, 0 );

        for ( unsigned int i = 0; i < dim; i++ ){
            for ( unsigned int j = 0; j < dim; j++ ){
                for ( unsigned int k = 0; k < dim; k++ ){
                    for ( unsigned int I = 0; I < dim; I++ ){
                        for ( unsigned int J = 0; J < dim; J++ ){
                            for ( unsigned int K = 0; K < dim; K++ ){
                                higherOrderStress[ dim * dim * i + dim * j + k ] += deformationGradient[ dim * i + I ]
                                                                                  * deformationGradient[ dim * j + J ]
                                                                                  * microDeformation[ dim * k + K ]
                                                                                  * referenceHigherOrderStress[ dim * dim * I + dim * J + K ];
                            }
                        }
                    }
                    higherOrderStress[ dim * dim * i + dim * j + k ] /= detF;
                }
            }
        }

        return NULL;
    }

    errorOut pushForwardHigherOrderStress( const variableVector &referenceHigherOrderStress,
                                           const variableVector &deformationGradient,
                                           const variableVector &microDeformation,
                                           variableVector &higherOrderStress,
                                           variableMatrix &dHigherOrderStressdReferenceHigherOrderStress,
                                           variableMatrix &dHigherOrderStressdDeformationGradient,
                                           variableMatrix &dHigherOrderStressdMicroDeformation ){
        /*!
         * Compute the push-forward operation on the higher order stress.
         *
         * m_{ijk} = \frac{1}{J} F_{iI} F_{jJ} \Chi_{kK} M_{IJK}
         *
         * Also returns the Jacobians
         *
         * \frac{ \partial m_{ijk} }{ \partial M_{LMN} } = \frac{1}{J} F_{iL} F_{jM} \Chi_{kN}
         * \frac{ \partial m_{ijk} }{ \partial F_{lM} } = \left( \delta_{il} F_{jN} \Chi_{kO} M_{MNO}
         *                                                     + F_{iN} \delta_{jl} \Chi_{kO} M_{NMO}
         *                                                     - m_{ijk} dDetFdF_{lM} \right)/J
         * \frac{ \partial m_{ijk} }{ \partial \Chi_{lM} } = \frac{1}{J} F_{iN} F_{jO} \delta_{kl} M_{NOM}
         *
         * :param const variableVector &referenceHigherOrderStress: The higher order stress in the 
         *     reference configuration.
         * :param const variableVector &deformationGradient: The deformation gradient which maps 
         *     between the reference and current configurations.
         * :param const variableVector &microDeformation: The micro-deformation tensor.
         * :param variableVector &higherOrderStress: The higher order stress in the current configuration.
         */

        //Assume 3d
        unsigned int dim = 3;

        variableType detF;

        errorOut error = pushForwardHigherOrderStress( referenceHigherOrderStress, deformationGradient, microDeformation,
                                                       detF, higherOrderStress );
        
        if (error){
            errorOut result = new errorNode( "pushForwardHigherOrderStress (jacobian)", "Error in computation of push forward of the higher order stress" );
            result->addNext(error);
            return result;
        }

        //Assemble the jacobian of the determinant of the deformation gradient
        variableVector inverseDeformationGradient = vectorTools::inverse( deformationGradient, dim, dim );

        variableVector dDetFdF( dim * dim, 0 );

        for (unsigned int i = 0; i < dim; i++ ){
            for (unsigned int I = 0; I < dim; I++ ){
                dDetFdF[ dim * i + I ] = inverseDeformationGradient[ dim * I + i ] * detF;
            }
        }

        constantVector eye( dim * dim );
        vectorTools::eye( eye );

        dHigherOrderStressdReferenceHigherOrderStress = variableMatrix( dim * dim * dim, variableVector( dim * dim * dim, 0 ) );
        dHigherOrderStressdDeformationGradient = variableMatrix( dim * dim * dim, variableVector( dim * dim, 0 ) );
        dHigherOrderStressdMicroDeformation = variableMatrix( dim * dim * dim, variableVector( dim * dim, 0 ) );

        for ( unsigned int i = 0; i < dim; i++ ){
            for ( unsigned int j = 0; j < dim; j++ ){
                for ( unsigned int k = 0; k < dim; k++ ){
                    for ( unsigned int l = 0; l < dim; l++ ){
                        for ( unsigned int M = 0; M < dim; M++ ){
                            for ( unsigned int N = 0; N < dim; N++ ){
                                dHigherOrderStressdReferenceHigherOrderStress[ dim * dim * i + dim * j + k ][ dim * dim * l + dim * M + N ] = deformationGradient[ dim * i + l ] * deformationGradient[ dim * j + M ] * microDeformation[ dim * k + N] / detF;
                                for ( unsigned int O = 0; O < dim; O++ ){
                                    dHigherOrderStressdDeformationGradient[ dim * dim * i + dim * j + k ][ dim * l + M ] += 
                                        eye[ dim * i + l ] * deformationGradient[ dim * j + N ] * microDeformation[ dim * k + O ] * referenceHigherOrderStress[ dim * dim * M + dim * N + O ]
                                      + deformationGradient[ dim * i + N ] * eye[ dim * j + l ] * microDeformation[ dim * k + O ] * referenceHigherOrderStress[ dim * dim * N + dim * M + O ];
                                    dHigherOrderStressdMicroDeformation[ dim * dim * i + dim * j + k ][ dim * l + M ] += deformationGradient[ dim * i + N ] * deformationGradient[ dim * j + O ] * eye[ dim * k + l ] * referenceHigherOrderStress[ dim * dim * N + dim * O + M ];
                                }
                            }
                            dHigherOrderStressdDeformationGradient[ dim * dim * i + dim * j + k ][ dim * l + M ] -= higherOrderStress[ dim * dim * i + dim * j + k ] * dDetFdF[ dim * l + M ];
                            dHigherOrderStressdDeformationGradient[ dim * dim * i + dim * j + k ][ dim * l + M ] /= detF;
                            dHigherOrderStressdMicroDeformation[ dim * dim * i + dim * j + k ][ dim * l + M ] /= detF;
                        }
                    }
                }
            }
        }

        return NULL;
    }

    errorOut pullBackHigherOrderStress( const variableVector &higherOrderStress,
                                        const variableVector &deformationGradient,
                                        const variableVector &microDeformation,
                                        variableVector &referenceHigherOrderStress ){
        /*!
         * Compute the pull back operation on the higher order stress.
         *
         * M_{IJK} = J F_{Ii}^{-1} F_{Jj}^{-1} \chi_{Kk}^{-1} m_{ijk}
         *
         * :param const variableVector &higherOrderStress: The higher order stress in the current configuration.
         * :param const variableVector &deformationGradient: The deformation gradient which maps 
         *     between the reference and current configurations.
         * :param const variableVector &microDeformation: The micro-deformation tensor.
         * :param variableVector &referenceHigherOrderStress: The higher order stress in the reference configuration.
         */

        variableType detF;
        variableVector inverseDeformationGradient, inverseMicroDeformation;
        return pullBackHigherOrderStress( higherOrderStress, deformationGradient,
                                          microDeformation, detF, inverseDeformationGradient,
                                          inverseMicroDeformation, referenceHigherOrderStress );
    }

    errorOut pullBackHigherOrderStress( const variableVector &higherOrderStress,
                                        const variableVector &deformationGradient,
                                        const variableVector &microDeformation,
                                        variableType &detF, variableVector &inverseDeformationGradient,
                                        variableVector &inverseMicroDeformation,
                                        variableVector &referenceHigherOrderStress ){
        /*!
         * Compute the push-forward operation on the higher order stress.
         *
         * M_{IJK} = J F_{Ii}^{-1} F_{Jj}^{-1} \chi_{Kk}^{-1} m_{ijk}
         *
         * :param const variableVector &higherOrderStress: The higher order stress in the current configuration.
         * :param const variableVector &deformationGradient: The deformation gradient which maps 
         *     between the reference and current configurations.
         * :param const variableVector &microDeformation: The micro-deformation tensor.
         * :param variableType &detF: The determinant of the deformation gradient.
         * :param variableVector &inverseDeformationGradient: The inverse of the deformation gradient.
         * :param variableVector &inverseMicroDeformation: The inverse of the micro deformation.
         * :param variableVector &referenceHigherOrderStress: The higher order stress in the reference configuration.
         */

        //Assume 3d
        unsigned int dim = 3;

        if ( higherOrderStress.size() != dim * dim * dim ){
            return new errorNode( "pullBackHigherOrderStress", "The current higher order stress doesn't have the correct size" );
        }

        if ( deformationGradient.size() != dim * dim ){
            return new errorNode( "pushForwardHigherOrderStress", "The deformation gradient doesn't have the correct size" );
        }

        if ( microDeformation.size() != dim * dim ){
            return new errorNode( "pushForwardHigherOrderStress", "The micro-deformation doesn't have the correct size" );
        }

        detF = vectorTools::determinant( deformationGradient, dim, dim );
        inverseDeformationGradient = vectorTools::inverse( deformationGradient, dim, dim );
        inverseMicroDeformation = vectorTools::inverse( microDeformation, dim, dim );

        referenceHigherOrderStress = variableVector( dim * dim * dim, 0 );

        for ( unsigned int I = 0; I < dim; I++ ){
            for ( unsigned int J = 0; J < dim; J++ ){
                for ( unsigned int K = 0; K < dim; K++ ){
                    for ( unsigned int i = 0; i < dim; i++ ){
                        for ( unsigned int j = 0; j < dim; j++ ){
                            for ( unsigned int k = 0; k < dim; k++ ){
                                referenceHigherOrderStress[ dim * dim * I + dim * J + K ]
                                    += inverseDeformationGradient[ dim * I + i ]
                                     * inverseDeformationGradient[ dim * J + j ]
                                     * inverseMicroDeformation[ dim * K + k ]
                                     * higherOrderStress[ dim * dim * i + dim * j + k ];
                            }
                        }
                    }
                    referenceHigherOrderStress[ dim * dim * I + dim * J + K ] *= detF;
                }
            }
        }

        return NULL;
    }

    errorOut pullBackHigherOrderStress( const variableVector &higherOrderStress,
                                        const variableVector &deformationGradient,
                                        const variableVector &microDeformation,
                                        variableVector &referenceHigherOrderStress,
                                        variableMatrix &dReferenceHigherOrderStressdHigherOrderStress,
                                        variableMatrix &dReferenceHigherOrderStressdDeformationGradient,
                                        variableMatrix &dReferenceHigherOrderStressdMicroDeformation ){
        /*!
         * Compute the pull back operation on the higher order stress.
         *
         * M_{IJK} = J F_{Ii}^{-1} F_{Jj}^{-1} \chi_{Kk}^{-1} m_{ijk}
         *
         * Also returns the Jacobians
         *
         * \frac{ \partial M_{IJK} }{ \partial m_{lmn} } = J F_{Il}^{-1} F_{Jm}^{-1} \chi_{Kn}^{-1}
         * \frac{ \partial M_{IJK} }{ \partial F_{lL} } = F_{Ll}^{-1} M_{IJK} - F_{Il}^{-1} M_{LJK} - F_{Jl}^{-1} M_{ILK}
         * \frac{ \partial M_{IJK} }{ \partial \chi_{lL} } = -\chi_{Kl}^{-1} M_{IJL}
         *
         * :param const variableVector &referenceHigherOrderStress: The higher order stress in the 
         *     reference configuration.
         * :param const variableVector &deformationGradient: The deformation gradient which maps 
         *     between the reference and current configurations.
         * :param const variableVector &microDeformation: The micro-deformation tensor.
         * :param variableVector &higherOrderStress: The higher order stress in the current configuration.
         */

        //Assume 3d
        unsigned int dim = 3;

        variableType detF;
        variableVector inverseDeformationGradient, inverseMicroDeformation;

        errorOut error = pullBackHigherOrderStress( higherOrderStress, deformationGradient, microDeformation,
                                                    detF, inverseDeformationGradient, inverseMicroDeformation,
                                                    referenceHigherOrderStress );
        
        if (error){
            errorOut result = new errorNode( "pullBackHigherOrderStress (jacobian)", "Error in computation of pull back of the higher order stress" );
            result->addNext(error);
            return result;
        }

        dReferenceHigherOrderStressdHigherOrderStress   = variableMatrix( dim * dim * dim, variableVector( dim * dim * dim, 0 ) );
        dReferenceHigherOrderStressdDeformationGradient = variableMatrix( dim * dim * dim, variableVector( dim * dim, 0 ) );
        dReferenceHigherOrderStressdMicroDeformation    = variableMatrix( dim * dim * dim, variableVector( dim * dim, 0 ) );

        for ( unsigned int I = 0; I < dim; I++){
            for ( unsigned int J = 0; J < dim; J++){
                for ( unsigned int K = 0; K < dim; K++){
                    for ( unsigned int l = 0; l < dim; l++){
                        for ( unsigned int m = 0; m < dim; m++){
                            dReferenceHigherOrderStressdDeformationGradient[ dim * dim * I + dim * J + K ][ dim * l + m ]
                                += inverseDeformationGradient[ dim * m + l ] * referenceHigherOrderStress[ dim * dim * I + dim * J + K ]
                                 - inverseDeformationGradient[ dim * I + l ] * referenceHigherOrderStress[ dim * dim * m + dim * J + K ]
                                 - inverseDeformationGradient[ dim * J + l ] * referenceHigherOrderStress[ dim * dim * I + dim * m + K ];

                            dReferenceHigherOrderStressdMicroDeformation[ dim * dim * I + dim * J + K ][ dim * l + m ]
                                -= inverseMicroDeformation[ dim * K + l ] * referenceHigherOrderStress[ dim * dim * I + dim * J + m ];

                            for ( unsigned int n = 0; n < dim; n++){

                                dReferenceHigherOrderStressdHigherOrderStress[ dim * dim * I + dim * J + K ][ dim * dim * l + dim * m + n ]
                                    += detF * inverseDeformationGradient[ dim * I + l ]
                                     * inverseDeformationGradient[ dim * J + m ]
                                     * inverseMicroDeformation[ dim * K + n ];
                            }
                        }
                    }
                }
            }
        }

        return NULL;
    }

    errorOut computeDeviatoricHigherOrderStress( const variableVector &higherOrderStress,
                                                 variableVector &deviatoricHigherOrderStress ){
        /*!
         * Compute the deviatoric part of the higher order stress.
         *
         * dev ( m_{ijk} ) = m_{ijk} - ( 1 / 3 ) m_{llk} \delta_{ij}
         *
         * :param const variableVector &higherOrderStress: The higher order stress in the current configuration.
         * :param variableVector &deviatoricHigherOrderStress: The deviatoric part of the higher order stress.
         */

        //Assume 3d
        unsigned int dim = 3;

        if ( higherOrderStress.size() != dim * dim * dim ){
            return new errorNode( "computeDeviatoricHigherOrderStress", "The higher order stress has an incorrect size" );
        }

        deviatoricHigherOrderStress = higherOrderStress;

        constantVector eye( dim * dim );
        vectorTools::eye( eye );

        for ( unsigned int i = 0; i < dim; i++ ){
            for ( unsigned int j = 0; j < dim; j++ ){
                for ( unsigned int k = 0; k < dim; k++ ){
                    for ( unsigned int l = 0; l < dim; l++ ){
                        deviatoricHigherOrderStress[ dim * dim * i + dim * j + k ] -= higherOrderStress[ dim * dim * l + dim * l + k ] 
                                                                                    * eye[ dim * i + j ] / 3;
                    }
                }
            }
        }

        return NULL;
    }

    errorOut computeDeviatoricHigherOrderStress( const variableVector &higherOrderStress,
                                                 variableVector &deviatoricHigherOrderStress,
                                                 variableMatrix &dDeviatoricHigherOrderStressdHigherOrderStress){
        /*!
         * Compute the deviatoric part of the higher order stress.
         *
         * dev ( m_{ijk} ) = m_{ijk} - ( 1 / 3 ) m_{llk} \delta_{ij}
         *
         * Also compute the Jacobian
         * \frac{ \partial dev ( m_{ijk} ) }{ \partial m_{mno} } = \delta_{im} \delta_{jn} \delta_{ko} - ( 1 / 3 ) \delta_{mn} \delta_{ko} \delta_{ij}
         *
         * :param const variableVector &higherOrderStress: The higher order stress in the current configuration.
         * :param variableVector &deviatoricHigherOrderStress: The deviatoric part of the higher order stress.
         * :param variableMatrix &dDeviatoricHigherOrderStressdHigherOrderStress: The gradient of the deviatoric part of the 
         *     higher order stress w.r.t. the higher order stress.
         */

        //Assume 3d
        unsigned int dim = 3;

        errorOut error = computeDeviatoricHigherOrderStress( higherOrderStress, deviatoricHigherOrderStress );

        if ( error ){
            errorOut result = new errorNode( "computeDeviatoricHigherOrderStress (jacobian)",
                                             "Error in the computation of the deviatoric part of the higher order stress" );
            result->addNext(error);
            return result;
        }

        constantVector eye( dim * dim );
        vectorTools::eye( eye );

        dDeviatoricHigherOrderStressdHigherOrderStress = variableMatrix( dim * dim * dim, variableVector( dim * dim * dim, 0 ) );

        for ( unsigned int i = 0; i < dim; i++ ){
            for ( unsigned int j = 0; j < dim; j++ ){
                for ( unsigned int k = 0; k < dim; k++ ){
                    for ( unsigned int m = 0; m < dim; m++ ){
                        for ( unsigned int n = 0; n < dim; n++ ){
                            for ( unsigned int o = 0; o < dim; o++ ){
                                dDeviatoricHigherOrderStressdHigherOrderStress[ dim * dim * i + dim * j + k ][ dim * dim * m + dim * n + o ] = eye[ dim * i + m ] * eye[ dim * j + n ] * eye[ dim * k + o ] - eye[ dim * m + n ] * eye[ dim * k + o ] * eye[ dim * i + j ] / 3;
                            }
                        }
                    }
                }
            }
        }

        return NULL;
    }

    errorOut computeReferenceHigherOrderStressPressure( const variableVector &referenceHigherOrderStress,
                                                        const variableVector &rightCauchyGreenDeformation,
                                                        variableVector &referenceHigherOrderPressure ){
        /*!
         * Compute the pressure for a higher-order stress in the reference configuration.
         * $p_K = \frac{1}{3} C_{AB} M_{ABK}$
         *
         * where $C_{AB}$ is the right Cauchy-Green deformation tensor and 
         * M_{ABK} is the higher order stress tensor in the reference configuration.
         *
         * :param const variableVector &referenceHigherOrderStress: The higher order stress in the 
         *     reference configuration.
         * :param const variableVector &rightCauchyGreenDeformation: The right Cauchy-Green deformation
         *     tensor.
         * :param variableVector &referenceHigherOrderPressure: The higher order pressure.
         */

        //Assume 3D
        unsigned int dim = 3;

        if ( rightCauchyGreenDeformation.size() != dim * dim ){
            return new errorNode( "computeReferenceHigherOrderStressPressure",
                                  "The right Cauchy-Green deformation tensor must have nine terms." );
        }

        if ( referenceHigherOrderStress.size() != dim * dim * dim ){
            return new errorNode( "computeReferenceHigherOrderStressPressure",
                                  "The higher order stress tensor must have 27 terms." );
        }

        referenceHigherOrderPressure = variableVector( dim, 0 );

        for ( unsigned int K = 0; K < dim; K++ ){
            for ( unsigned int A = 0; A < dim; A++ ){
                for ( unsigned int B = 0; B < dim; B++ ){
                    referenceHigherOrderPressure[K] += rightCauchyGreenDeformation[ dim * A + B ]
                                                     * referenceHigherOrderStress[ dim * dim * A + dim * B + K ];
                }
            }
        }

        referenceHigherOrderPressure /= 3;
        return NULL;
    }

    errorOut computeReferenceHigherOrderStressPressure( const variableVector &referenceHigherOrderStress,
                                                        const variableVector &rightCauchyGreenDeformation,
                                                        variableVector &referenceHigherOrderPressure,
                                                        variableMatrix &dpdM, variableMatrix &dpdC ){
        /*!
         * Compute the pressure for a higher-order stress in the reference configuration.
         * $p_K = \frac{1}{3} C_{AB} M_{ABK}$
         *
         * Also compute the Jacobians
         * $\frac{ \partial p_K }{ \partial M_{NOP} } = \frac{1}{3} C_{NO} \delta_{KP}$
         * $\frac{ \partial p_K }{ \partial C_{NO} } = \frac{1}{3} M_{NOK}$
         *
         * where $C_{AB}$ is the right Cauchy-Green deformation tensor and 
         * M_{ABK} is the higher order stress tensor in the reference configuration.
         *
         * :param const variableVector &referenceHigherOrderStress: The higher order stress in the 
         *     reference configuration.
         * :param const variableVector &rightCauchyGreenDeformation: The right Cauchy-Green deformation
         *     tensor.
         * :param variableVector &referenceHigherOrderPressure: The higher order pressure.
         * :param variableMatrix &dpdM: The Jacobian of the pressure w.r.t. the higher order stress.
         * :param variableMatrix &dpdC: The Jacobian of the pressure w.r.t. the right Cauchy-Green 
         *     deformation tensor.
         */

        //Assume 3D
        unsigned int dim = 3;

        errorOut error = computeReferenceHigherOrderStressPressure( referenceHigherOrderStress, rightCauchyGreenDeformation,
                                                                    referenceHigherOrderPressure );

        if ( error ){
            errorOut result = new errorNode( "computeReferenceHigherOrderStressPressure (jacobian)",
                                             "Error in computation of reference higher order pressure" );
            result->addNext( error );
            return result;
        }

        constantVector eye( dim * dim );
        vectorTools::eye( eye );

        dpdM = variableMatrix( referenceHigherOrderPressure.size(), variableVector( referenceHigherOrderStress.size(), 0 ) );
        dpdC = variableMatrix( referenceHigherOrderPressure.size(), variableVector( rightCauchyGreenDeformation.size(), 0 ) );

        for ( unsigned int K = 0; K < dim; K++ ){
            for ( unsigned int N = 0; N < dim; N++ ){
                for ( unsigned int O = 0; O < dim; O++ ){
                    dpdC[ K ][ dim * N + O ] = referenceHigherOrderStress[ dim * dim * N + dim * O + K ];
                    for ( unsigned int P = 0; P < dim; P++ ){
                        dpdM[ K ][ dim * dim * N + dim * O + P ] = rightCauchyGreenDeformation[ dim * N + O ] * eye[ dim * K + P ];
                    }
                }
            }
        }

        dpdM /= 3;
        dpdC /= 3;

        return NULL;
    }

    errorOut computeReferenceHigherOrderStressPressure( const variableVector &referenceHigherOrderStress,
                                                        const variableVector &rightCauchyGreenDeformation,
                                                        variableVector &referenceHigherOrderPressure,
                                                        variableMatrix &dpdM, variableMatrix &dpdC,
                                                        variableMatrix &d2pdMdC ){
        /*!
         * Compute the pressure for a higher-order stress in the reference configuration.
         * $p_K = \frac{1}{3} C_{AB} M_{ABK}$
         *
         * Also compute the Jacobians
         * $\frac{ \partial p_K }{ \partial M_{NOP} } = \frac{1}{3} C_{NO} \delta_{KP}$
         * $\frac{ \partial p_K }{ \partial C_{NO} } = \frac{1}{3} M_{NOK}$
         * $\frac{ \partial^2 p_K}{ \partial M_{NOP} C_{QR} } = \frac{1}{3} \delta_{NQ} \delta_{OR} \delta_{KP}
         *
         * where $C_{AB}$ is the right Cauchy-Green deformation tensor and
         * M_{ABK} is the higher order stress tensor in the reference configuration.
         *
         * :param const variableVector &referenceHigherOrderStress: The higher order stress in the
         *     reference configuration.
         * :param const variableVector &rightCauchyGreenDeformation: The right Cauchy-Green deformation
         *     tensor.
         * :param variableVector &referenceHigherOrderPressure: The higher order pressure.
         * :param variableMatrix &dpdM: The Jacobian of the pressure w.r.t. the higher order stress.
         * :param variableMatrix &dpdC: The Jacobian of the pressure w.r.t. the right Cauchy-Green
         *     deformation tensor.
         * :param variableMatrix &d2pdMdC: The second order jacobian of the pressure w.r.t the 
         *     reference higher order stress and right Cauchy-Green deformation tensor. This Jacobian is organized
         *     [ K ][ NOPQR ]
         */

        //Assume 3D
        unsigned int dim = 3;

        errorOut error = computeReferenceHigherOrderStressPressure( referenceHigherOrderStress, rightCauchyGreenDeformation,
                                                                    referenceHigherOrderPressure, dpdM, dpdC );

        if ( error ){
            errorOut result = new errorNode( "computeReferenceHigherOrderStressPressure (second order jacobian)",
                                             "Error in computation of reference higher order pressure" );
            result->addNext( error );
            return result;
        }

        variableVector eye( dim * dim );
        vectorTools::eye( eye );

        d2pdMdC = variableMatrix( dim, variableVector( dim * dim * dim * dim * dim, 0 ) );

        for ( unsigned int K = 0; K < dim; K++ ){
            for ( unsigned int N = 0; N < dim; N++ ){
                for ( unsigned int O = 0; O < dim; O++ ){
                    for ( unsigned int P = 0; P < dim; P++ ){
                        for ( unsigned int Q = 0; Q < dim; Q++ ){
                            for ( unsigned int R = 0; R < dim; R++ ){
                                d2pdMdC[ K ][ dim * dim * dim * dim * N + dim * dim * dim * O + dim * dim * P + dim * Q + R ] = 
                                    eye[ dim * N + Q ] * eye[ dim * O + R ] * eye[ dim * K + P ] / 3;
                            }
                        }
                    }
                }
            }
        }

        return NULL;
    }

    errorOut computeDeviatoricReferenceHigherOrderStress( const variableVector &referenceHigherOrderStress,
                                                          const variableVector &rightCauchyGreenDeformation,
                                                          variableVector &deviatoricReferenceHigherOrderStress ){
        /*!
         * Compute the deviatoric part of the higher order stress in the reference configuration.
         *
         * dev ( M_{IJK} ) = M_{IJK} - ( 1 / 3 ) (C^{-1})_{IJ} C_{AB} M_{ABK}
         *
         * :param const variableVector &referenceHigherOrderStress: The higher order stress in the reference configuration.
         * :param const variableVector &rightCauchyGreenDeformation: The right Cauchy-Green deformation tensor.
         * :param variableVector &deviatoricReferenceHigherOrderStress: The deviatoric part of the higher order tensor in the 
         *     reference configuration.
         */

        //Compute the pressure
        variableVector pressure;
        errorOut error = computeReferenceHigherOrderStressPressure( referenceHigherOrderStress, rightCauchyGreenDeformation, pressure );

        if ( error ){
            errorOut result = new errorNode( "computeDeviatoricReferenceHigherOrderStress",
                                             "Error in computation of higher order pressure" );
            result->addNext( error );
            return result;
        }

        return computeDeviatoricReferenceHigherOrderStress( referenceHigherOrderStress, rightCauchyGreenDeformation,
                                                            pressure, deviatoricReferenceHigherOrderStress );
    }

    errorOut computeDeviatoricReferenceHigherOrderStress( const variableVector &referenceHigherOrderStress,
                                                          const variableVector &rightCauchyGreenDeformation,
                                                          const variableVector &pressure,
                                                          variableVector &deviatoricReferenceHigherOrderStress ){
        /*!
         * Compute the deviatoric part of the higher order stress in the reference configuration.
         *
         * dev ( M_{IJK} ) = M_{IJK} - ( 1 / 3 ) (C^{-1})_{IJ} C_{AB} M_{ABK}
         *
         * :param const variableVector &referenceHigherOrderStress: The higher order stress in the reference configuration.
         * :param const variableVector &rightCauchyGreenDeformation: The right Cauchy-Green deformation tensor.
         * :param const variableType &pressure: The pressure of the higher-order stress.
         * :param variableVector &deviatoricReferenceHigherOrderStress: The deviatoric part of the higher order tensor in the 
         *     reference configuration.
         */

        //Assume 3d
        unsigned int dim = 3;

        deviatoricReferenceHigherOrderStress = referenceHigherOrderStress;

        variableVector invRCG = vectorTools::inverse( rightCauchyGreenDeformation, dim, dim );
        
        for ( unsigned int I = 0; I < dim; I++ ){
            for ( unsigned int J = 0; J < dim; J++ ){
                for ( unsigned int K = 0; K < dim; K++ ){
                    deviatoricReferenceHigherOrderStress[ dim * dim * I + dim * J + K ] -= invRCG[ dim * I + J ] * pressure[ K ];
                }
            }
        }

        return NULL;
    }

    errorOut computeDeviatoricReferenceHigherOrderStress( const variableVector &referenceHigherOrderStress,
                                                          const variableVector &rightCauchyGreenDeformation,
                                                          variableVector &deviatoricReferenceHigherOrderStress,
                                                          variableMatrix &dDeviatoricReferenceHigherOrderStressdReferenceHigherOrderStress,
                                                          variableMatrix &dDeviatoricReferenceHigherOrderStressdRCG ){
        /*!
         * Compute the deviatoric part of the higher order stress in the reference configuration.
         *
         * dev ( M_{IJK} ) = M_{IJK} - \frac{1}{3} (C^{-1})_{IJ} C_{AB} M_{ABK} = M_{IJK} - ( C^{-1} )_{IJ} p_{K}
         *
         * Also compute Jacobians:
         * \frac{ \partial dev ( M_{IJK} ) }{ \partial M_{LMN} } = \delta_{IL} \delta_{JM} \delta_{KN} - (C^{-1})_{IJ} \frac{ \partial p_K }{ \partial M_{LMN} }
         * \frac{ \partial dev ( M_{IJK} ) }{ \partial C_{LM} } = \left( (C^{-1})_{IL} (C^{-1})_{MJ} p_{K} - (C^{-1})_{IJ} \frac{ \partial p_K }{ \partial C_{LM} } \right)
         *
         * :param const variableVector &referenceHigherOrderStress: The higher order stress in the reference configuration.
         * :param const variableVector &rightCauchyGreenDeformation: The right Cauchy-Green deformation tensor.
         * :param variableVector &deviatoricReferenceHigherOrderStress: The deviatoric part of the higher order tensor in the 
         *     reference configuration.
         * :param variableVector &dDeviatoricReferenceHigherOrdeterStressdReferenceHigherOrderStress: The Jacobian of the 
         *     deviatoric part of the reference higher order stress w.r.t. the reference higher order stress.
         * :param variableVector &dDeviatoricReferenceHigherOrdeterStressdRCG: The Jacobian of the 
         *     deviatoric part of the reference higher order stress w.r.t. the right Cauchy-Green deformation tensor.
         */

        variableVector invRCG, pressure;
        variableMatrix dpdM, dpdC;

        //Compute the pressure
        errorOut error = computeReferenceHigherOrderStressPressure( referenceHigherOrderStress, rightCauchyGreenDeformation, pressure,
                                                                    dpdM, dpdC );

        if ( error ){
            errorOut result = new errorNode( "computeDeviatoricReferenceHigherOrderStress (jacobian)",
                                             "Error in computation of higher order pressure" );
            result->addNext( error );
            return result;
        }

        return computeDeviatoricReferenceHigherOrderStress( referenceHigherOrderStress, rightCauchyGreenDeformation,
                                                            pressure, dpdM, dpdC, deviatoricReferenceHigherOrderStress,
                                                            dDeviatoricReferenceHigherOrderStressdReferenceHigherOrderStress,
                                                            dDeviatoricReferenceHigherOrderStressdRCG );
    }

    errorOut computeDeviatoricReferenceHigherOrderStress( const variableVector &referenceHigherOrderStress,
                                                          const variableVector &rightCauchyGreenDeformation,
                                                          const variableVector &pressure,
                                                          const variableMatrix &dPressuredStress,
                                                          const variableMatrix &dPressuredRCG,
                                                          variableVector &deviatoricReferenceHigherOrderStress,
                                                          variableMatrix &dDeviatoricReferenceHigherOrderStressdReferenceHigherOrderStress,
                                                          variableMatrix &dDeviatoricReferenceHigherOrderStressdRCG ){
        /*!
         * Compute the deviatoric part of the higher order stress in the reference configuration.
         *
         * dev ( M_{IJK} ) = M_{IJK} - \frac{1}{3} (C^{-1})_{IJ} C_{AB} M_{ABK} = M_{IJK} - ( C^{-1} )_{IJ} p_{K}
         *
         * Also compute Jacobians:
         * \frac{ \partial dev ( M_{IJK} ) }{ \partial M_{LMN} } = \delta_{IL} \delta_{JM} \delta_{KN} - (C^{-1})_{IJ} \frac{ \partial p_K }{ \partial M_{LMN} }
         * \frac{ \partial dev ( M_{IJK} ) }{ \partial C_{LM} } = \left( (C^{-1})_{IL} (C^{-1})_{MJ} p_{K} - (C^{-1})_{IJ} \frac{ \partial p_K }{ \partial C_{LM} } \right)
         *
         * :param const variableVector &referenceHigherOrderStress: The higher order stress in the reference configuration.
         * :param const variableVector &rightCauchyGreenDeformation: The right Cauchy-Green deformation tensor.
         * :param const variableVector &pressure: The higher order pressure.
         * :param const variableMatrix &dPressuredStress: The Jacobian of the pressure w.r.t. the stress.
         * :param const variableMatrix &dPressuredStress: The Jacobian of the pressure w.r.t. the right Cauchy Green 
         *     deformation tensor.
         * :param variableVector &deviatoricReferenceHigherOrderStress: The deviatoric part of the higher order tensor in the 
         *     reference configuration.
         * :param variableVector &dDeviatoricReferenceHigherOrdeterStressdReferenceHigherOrderStress: The Jacobian of the 
         *     deviatoric part of the reference higher order stress w.r.t. the reference higher order stress.
         * :param variableVector &dDeviatoricReferenceHigherOrdeterStressdRCG: The Jacobian of the 
         *     deviatoric part of the reference higher order stress w.r.t. the right Cauchy-Green deformation tensor.
         */
        
        //Assume 3d
        unsigned int dim = 3;

        deviatoricReferenceHigherOrderStress = referenceHigherOrderStress;

        variableVector invRCG = vectorTools::inverse( rightCauchyGreenDeformation, dim, dim );

        //Compute the deviatoric higher order stress        
        for ( unsigned int I = 0; I < dim; I++ ){
            for ( unsigned int J = 0; J < dim; J++ ){
                for ( unsigned int K = 0; K < dim; K++ ){
                    deviatoricReferenceHigherOrderStress[ dim * dim * I + dim * J + K ] -= invRCG[ dim * I + J ] * pressure[ K ];
                }
            }
        }

        //Compute the jacobians
        constantVector eye( dim * dim );
        vectorTools::eye( eye );

        dDeviatoricReferenceHigherOrderStressdReferenceHigherOrderStress = variableMatrix( dim * dim * dim, variableVector( dim * dim * dim, 0 ) );
        dDeviatoricReferenceHigherOrderStressdRCG = variableMatrix( dim * dim * dim, variableVector( dim * dim, 0 ) );

        for ( unsigned int I = 0; I < dim; I++ ){
            for ( unsigned int J = 0; J < dim; J++ ){
                for ( unsigned int K = 0; K < dim; K++ ){
                    for ( unsigned int L = 0; L < dim; L++ ){
                        for ( unsigned int M = 0; M < dim; M++ ){
                            for ( unsigned int N = 0; N < dim; N++ ){
                                dDeviatoricReferenceHigherOrderStressdReferenceHigherOrderStress[ dim * dim * I + dim * J + K ][ dim * dim * L + dim * M + N ]
                                    = eye[ dim * I + L ] * eye[ dim * J + M ] * eye[ dim * K + N ] - invRCG[ dim * I + J ] * dPressuredStress[ K ][ dim * dim * L + dim * M + N ];
                            }

                            dDeviatoricReferenceHigherOrderStressdRCG[ dim * dim * I + dim * J + K ][ dim * L + M ]
                                = invRCG[ dim * I + L ] * invRCG[ dim * M + J ] * pressure[ K ] - invRCG[ dim * I + J ] * dPressuredRCG[ K ][ dim * L + M ];
                        }
                    }
                }
            }
        }

        return NULL;
    }

    errorOut computeDeviatoricReferenceHigherOrderStress( const variableVector &referenceHigherOrderStress,
                                                          const variableVector &rightCauchyGreenDeformation,
                                                          variableVector &deviatoricReferenceHigherOrderStress,
                                                          variableMatrix &dDeviatoricReferenceHigherOrderStressdReferenceHigherOrderStress,
                                                          variableMatrix &dDeviatoricReferenceHigherOrderStressdRCG,
                                                          variableMatrix &d2MdMdRCG ){
        /*!
         * Compute the deviatoric part of the higher order stress in the reference configuration.
         *
         * dev ( M_{IJK} ) = M_{IJK} - \frac{1}{3} (C^{-1})_{IJ} C_{AB} M_{ABK} = M_{IJK} - ( C^{-1} )_{IJ} p_{K}
         *
         * Also compute Jacobians:
         * \frac{ \partial dev ( M_{IJK} ) }{ \partial M_{LMN} } = \delta_{IL} \delta_{JM} \delta_{KN} - (C^{-1})_{IJ} \frac{ \partial p_K }{ \partial M_{LMN} }
         * \frac{ \partial dev ( M_{IJK} ) }{ \partial C_{LM} } = \left( (C^{-1})_{IL} (C^{-1})_{MJ} p_{K} - (C^{-1})_{IJ} \frac{ \partial p_K }{ \partial C_{LM} } \right)
         *
         * \frac{ \partial dev( M_{IJK} ) }{ \partial M_{LMN} \partial C_{OP} } = (C^{-1})_{IO} (C^{-1})_{PJ} \frac{ \partial p_{k }{ \partial M_{LMN} } - (C^{-1})_{IJ} \frac{ \partial^2 p_K}{ \partial M_{LMN} \partial C_{OP} } 
         *
         * :param const variableVector &referenceHigherOrderStress: The higher order stress in the reference configuration.
         * :param const variableVector &rightCauchyGreenDeformation: The right Cauchy-Green deformation tensor.
         * :param variableVector &deviatoricReferenceHigherOrderStress: The deviatoric part of the higher order tensor in the 
         *     reference configuration.
         * :param variableVector &dDeviatoricReferenceHigherOrdeterStressdReferenceHigherOrderStress: The Jacobian of the 
         *     deviatoric part of the reference higher order stress w.r.t. the reference higher order stress.
         * :param variableVector &dDeviatoricReferenceHigherOrdeterStressdRCG: The Jacobian of the 
         *     deviatoric part of the reference higher order stress w.r.t. the right Cauchy-Green deformation tensor.
         * :param variableMatrix &d2DevMdMdRCG: The mixed second derivative of the deviatoric part of the reference higher order 
         *     stress tensor with respect to the reference higher order stress tensor and the right Cauchy-Green deformation 
         *     tensor.
         */

        variableVector invRCG, pressure;
        variableMatrix dpdM, dpdC, d2pdMdRCG;

        //Compute the pressure
        errorOut error = computeReferenceHigherOrderStressPressure( referenceHigherOrderStress, rightCauchyGreenDeformation, pressure,
                                                                    dpdM, dpdC, d2pdMdRCG );

        if ( error ){
            errorOut result = new errorNode( "computeDeviatoricReferenceHigherOrderStress (second order jacobian)",
                                             "Error in computation of higher order pressure" );
            result->addNext( error );
            return result;
        }

        return computeDeviatoricReferenceHigherOrderStress( referenceHigherOrderStress, rightCauchyGreenDeformation,
                                                            pressure, dpdM, dpdC, d2pdMdRCG,
                                                            deviatoricReferenceHigherOrderStress,
                                                            dDeviatoricReferenceHigherOrderStressdReferenceHigherOrderStress,
                                                            dDeviatoricReferenceHigherOrderStressdRCG,
                                                            d2MdMdRCG );
    }

    errorOut computeDeviatoricReferenceHigherOrderStress( const variableVector &referenceHigherOrderStress,
                                                          const variableVector &rightCauchyGreenDeformation,
                                                          const variableVector &pressure,
                                                          const variableMatrix &dPressuredStress,
                                                          const variableMatrix &dPressuredRCG,
                                                          const variableMatrix &d2PressuredStressdRCG,
                                                          variableVector &deviatoricReferenceHigherOrderStress,
                                                          variableMatrix &dDeviatoricReferenceHigherOrderStressdReferenceHigherOrderStress,
                                                          variableMatrix &dDeviatoricReferenceHigherOrderStressdRCG,
                                                          variableMatrix &d2MdMdRCG ){
        /*!
         * Compute the deviatoric part of the higher order stress in the reference configuration.
         *
         * dev ( M_{IJK} ) = M_{IJK} - \frac{1}{3} (C^{-1})_{IJ} C_{AB} M_{ABK} = M_{IJK} - ( C^{-1} )_{IJ} p_{K}
         *
         * Also compute Jacobians:
         * \frac{ \partial dev ( M_{IJK} ) }{ \partial M_{LMN} } = \delta_{IL} \delta_{JM} \delta_{KN} - (C^{-1})_{IJ} \frac{ \partial p_K }{ \partial M_{LMN} }
         * \frac{ \partial dev ( M_{IJK} ) }{ \partial C_{LM} } = \left( (C^{-1})_{IL} (C^{-1})_{MJ} p_{K} - (C^{-1})_{IJ} \frac{ \partial p_K }{ \partial C_{LM} } \right)
         *
         * \frac{ \partial dev( M_{IJK} ) }{ \partial M_{LMN} \partial C_{OP} } = (C^{-1})_{IO} (C^{-1})_{PJ} \frac{ \partial p_{k }{ \partial M_{LMN} } - (C^{-1})_{IJ} \frac{ \partial^2 p_K}{ \partial M_{LMN} \partial C_{OP} } 
         *
         * :param const variableVector &referenceHigherOrderStress: The higher order stress in the reference configuration.
         * :param const variableVector &rightCauchyGreenDeformation: The right Cauchy-Green deformation tensor.
         * :param const variableVector &pressure: The higher order pressure.
         * :param const variableMatrix &dPressuredStress: The Jacobian of the pressure w.r.t. the stress.
         * :param const variableMatrix &dPressuredStress: The Jacobian of the pressure w.r.t. the right Cauchy-Green 
         *     deformation tensor.
         * :param const variableMatrix &d2PressuredStressdRCG: The second order Jacobian of the pressure w.r.t. the 
         *     stress and the right Cauchy-Green deformation tensor.
         * :param variableVector &deviatoricReferenceHigherOrderStress: The deviatoric part of the higher order tensor in the 
         *     reference configuration.
         * :param variableVector &dDeviatoricReferenceHigherOrdeterStressdReferenceHigherOrderStress: The Jacobian of the 
         *     deviatoric part of the reference higher order stress w.r.t. the reference higher order stress.
         * :param variableVector &dDeviatoricReferenceHigherOrdeterStressdRCG: The Jacobian of the 
         *     deviatoric part of the reference higher order stress w.r.t. the right Cauchy-Green deformation tensor.
         * :param variableMatrix &d2DevMdMdRCG: The mixed second derivative of the deviatoric part of the reference higher order 
         *     stress tensor with respect to the reference higher order stress tensor and the right Cauchy-Green deformation 
         *     tensor.
         */
        
        //Assume 3d
        unsigned int dim = 3;

        deviatoricReferenceHigherOrderStress = referenceHigherOrderStress;

        variableVector invRCG = vectorTools::inverse( rightCauchyGreenDeformation, dim, dim );

        //Compute the deviatoric higher order stress        
        for ( unsigned int I = 0; I < dim; I++ ){
            for ( unsigned int J = 0; J < dim; J++ ){
                for ( unsigned int K = 0; K < dim; K++ ){
                    deviatoricReferenceHigherOrderStress[ dim * dim * I + dim * J + K ] -= invRCG[ dim * I + J ] * pressure[ K ];
                }
            }
        }

        //Compute the first order Jacobians
        constantVector eye( dim * dim );
        vectorTools::eye( eye );

        dDeviatoricReferenceHigherOrderStressdReferenceHigherOrderStress = variableMatrix( dim * dim * dim, variableVector( dim * dim * dim, 0 ) );
        dDeviatoricReferenceHigherOrderStressdRCG = variableMatrix( dim * dim * dim, variableVector( dim * dim, 0 ) );

        for ( unsigned int I = 0; I < dim; I++ ){
            for ( unsigned int J = 0; J < dim; J++ ){
                for ( unsigned int K = 0; K < dim; K++ ){
                    for ( unsigned int L = 0; L < dim; L++ ){
                        for ( unsigned int M = 0; M < dim; M++ ){
                            for ( unsigned int N = 0; N < dim; N++ ){
                                dDeviatoricReferenceHigherOrderStressdReferenceHigherOrderStress[ dim * dim * I + dim * J + K ][ dim * dim * L + dim * M + N ]
                                    = eye[ dim * I + L ] * eye[ dim * J + M ] * eye[ dim * K + N ] - invRCG[ dim * I + J ] * dPressuredStress[ K ][ dim * dim * L + dim * M + N ];
                            }

                            dDeviatoricReferenceHigherOrderStressdRCG[ dim * dim * I + dim * J + K ][ dim * L + M ]
                                = invRCG[ dim * I + L ] * invRCG[ dim * M + J ] * pressure[ K ] - invRCG[ dim * I + J ] * dPressuredRCG[ K ][ dim * L + M ];
                        }
                    }
                }
            }
        }

        //Compute the second order Jacobian

        d2MdMdRCG = variableMatrix( dim * dim * dim, variableVector( dim * dim * dim * dim * dim, 0 ) );

        for ( unsigned int I = 0; I < dim; I++ ){
            for ( unsigned int J = 0; J < dim; J++ ){
                for ( unsigned int K = 0; K < dim; K++ ){
                    for ( unsigned int L = 0; L < dim; L++ ){
                        for ( unsigned int M = 0; M < dim; M++ ){
                            for ( unsigned int N = 0; N < dim; N++ ){
                                for ( unsigned int O = 0; O < dim; O++ ){
                                    for ( unsigned int P = 0; P < dim; P++ ){
                                        d2MdMdRCG[ dim * dim * I + dim * J + K ][ dim * dim * dim * dim * L + dim * dim * dim * M + dim * dim * N + dim * O + P ]
                                            = invRCG[ dim * I + O ] * invRCG[ dim * P + J ] * dPressuredStress[ K ][ dim * dim * L + dim * M + N ]
                                            - invRCG[ dim * I + J ] * d2PressuredStressdRCG[ K ][ dim * dim * dim * dim * L + dim * dim * dim * M + dim * dim * N + dim * O + P ];
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        return NULL;
    }

    errorOut computeDeviatoricSecondOrderStress( const variableVector &secondOrderStress,
                                                 variableVector &deviatoricSecondOrderStress ){
        /*!
         * Compute the deviatoric part of a second order stress measure in the current configuration.
         * \hat{s}_{ij} = s_{ij} - \frac{1}{3} s_{mm} \delta_{ij}
         *
         * :param const variableVector &secondOrderStress: The stress measure in the current configuration.
         * :param variableVector &deviatoricSecondOrderStress: The deviatoric part of the second order stress 
         *     measure.
         */

        //Assume 3d
        unsigned int dim = 3;

        variableVector eye( dim * dim );
        vectorTools::eye( eye );

        variableType trace = vectorTools::trace( secondOrderStress );

        deviatoricSecondOrderStress = secondOrderStress - trace * eye / 3.;

        return NULL;
    }

    errorOut computeDeviatoricSecondOrderStress( const variableVector &secondOrderStress,
                                                 variableVector &deviatoricSecondOrderStress,
                                                 variableMatrix &dDeviatoricStressdStress ){
        /*!
         * Compute the deviatoric part of a second order stress measure in the current configuration.
         * \hat{s}_{ij} = s_{ij} - \frac{1}{3} s_{mm} \delta_{ij}
         *
         * Also return the jacobian
         * \frac{\partial \hat{s}_{ij}}{\partial s_{kl} } \delta_{ik} \delta_{jl} - \frac{1}{3} \delta_{ij} \delta_{kl}
         *
         * :param const variableVector &secondOrderStress: The stress measure in the current configuration.
         * :param variableVector &deviatoricSecondOrderStress: The deviatoric part of the second order stress 
         *     measure.
         * :param variableMatrix &dDeviatoricStressdStress: The jacobian of the deviatoric stress.
         */

        //Assume 3d
        unsigned int dim = 3;

        variableVector eye( dim * dim );
        vectorTools::eye( eye );

        errorOut error = computeDeviatoricSecondOrderStress( secondOrderStress, deviatoricSecondOrderStress );

        if ( error ){
            errorOut result = new errorNode( "computeDeviatoricSecondOrderStress (jacobian)",
                                             "Error in the computation of the deviatoric stress" );
            result->addNext( error );
            return result;
        }

        dDeviatoricStressdStress = variableMatrix( deviatoricSecondOrderStress.size(), variableVector( secondOrderStress.size(), 0 ) );

        for ( unsigned int i = 0; i < dim; i++ ){
            for ( unsigned int j = 0; j < dim; j++ ){
                for ( unsigned int k = 0; k < dim; k++ ){
                    for ( unsigned int l = 0; l < dim; l++ ){
                        dDeviatoricStressdStress[ dim * i + j ][ dim * k + l ] = eye[ dim * i + k ] * eye[ dim * j + l ]
                                                                               - eye[ dim * i + j ] * eye[ dim * k + l ] / 3;
                    }
                }
            }
        }

        return NULL;
    }

    errorOut computeReferenceSecondOrderStressPressure( const variableVector &referenceStressMeasure,
                                                        const variableVector &rightCauchyGreen,  variableType &pressure ){
        /*!
         * Compute the pressure part of a second order stress measure in the reference configuration.
         * p = \frac{1}{3} C_{IJ} S_{IJ}
         *
         * where C_{IJ} is the right Cauchy-Green deformation tensor and S_{IJ} is the stress measure.
         *
         * :param const variableVector &referenceStressMeasure: The stress measure in the reference configuration.
         * :param const variableVector &rightCauchyGreen: The right Cauchy-Green deformation tensor between the 
         *     current configuration and the reference configuration of the stress tensor.
         * :param variableType &pressure: The computed pressure.
         */

        if ( referenceStressMeasure.size() != rightCauchyGreen.size() ){
            return new errorNode( "computeReferenceSecondOrderStressPressure",
                                  "The stress measure and right Cauchy-Green deformation tensors aren't the same size" );
        }

        pressure = vectorTools::dot( referenceStressMeasure, rightCauchyGreen ) / 3;

        return NULL;
    }

    errorOut computeReferenceSecondOrderStressPressure( const variableVector &referenceStressMeasure,
                                                        const variableVector &rightCauchyGreen,  variableType &pressure,
                                                        variableVector &dpdStress, variableVector &dpdRCG ){
        /*!
         * Compute the pressure part of a second order stress measure in the reference configuration.
         * p = \frac{1}{3} C_{IJ} S_{IJ}
         *
         * where C_{IJ} is the right Cauchy-Green deformation tensor and S_{IJ} is the stress measure.
         *
         * Also compute the Jacobians
         * \frac{ \partial p }{ \partial C_{IJ} } = \frac{1}{3} S_{IJ}
         * \frac{ \partial p }{ \partial S_{IJ} } = \frac{1}{3} C_{IJ}
         *
         * :param const variableVector &referenceStressMeasure: The stress measure in the reference configuration.
         * :param const variableVector &rightCauchyGreen: The right Cauchy-Green deformation tensor between the 
         *     current configuration and the reference configuration of the stress tensor.
         * :param variableType &pressure: The computed pressure.
         * :param variableVector &dpdStress: The Jacobian of the pressure w.r.t. the stress.
         * :param variableVector &dpdRCG: The Jacobian of the pressure w.r.t. the right Cauchy-Green deformation tensor.
         */

        errorOut error = computeReferenceSecondOrderStressPressure( referenceStressMeasure, rightCauchyGreen, pressure );

        if ( error ){
            errorOut result = new errorNode( "computeReferenceSecondOrderStressPressure (jacobian)",
                                             "Error in computation of pressure in the reference configuration" );
            result->addNext( error );
            return result;
        }

        dpdStress = rightCauchyGreen / 3;
        dpdRCG    = referenceStressMeasure / 3;

        return NULL;
    }

    errorOut computeReferenceSecondOrderStressPressure( const variableVector &referenceStressMeasure,
                                                        const variableVector &rightCauchyGreen,  variableType &pressure,
                                                        variableVector &dpdStress, variableVector &dpdRCG,
                                                        variableMatrix &d2pdStressdRCG ){
        /*!
         * Compute the pressure part of a second order stress measure in the reference configuration.
         * p = \frac{1}{3} C_{IJ} S_{IJ}
         *
         * where C_{IJ} is the right Cauchy-Green deformation tensor and S_{IJ} is the stress measure.
         *
         * Also compute the Jacobians
         * \frac{ \partial p }{ \partial C_{IJ} } = \frac{1}{3} S_{IJ}
         * \frac{ \partial p }{ \partial S_{IJ} } = \frac{1}{3} C_{IJ}
         *
         * \frac{ \partial^2 p}{ \partial S_{IJ} \partial C_{KL} } = \frac{1}{3} \delta_{IK} \delta_{JL}
         *
         * :param const variableVector &referenceStressMeasure: The stress measure in the reference configuration.
         * :param const variableVector &rightCauchyGreen: The right Cauchy-Green deformation tensor between the 
         *     current configuration and the reference configuration of the stress tensor.
         * :param variableType &pressure: The computed pressure.
         * :param variableVector &dpdStress: The Jacobian of the pressure w.r.t. the stress.
         * :param variableVector &dpdRCG: The Jacobian of the pressure w.r.t. the right Cauchy-Green deformation tensor.
         * :param variableMatrix $d2pdStressdRCG: The second-order Jacobian of the pressure w.r.t. the stress and the 
         *     right Cauchy-Green deformation tensor.
         */

        errorOut error = computeReferenceSecondOrderStressPressure( referenceStressMeasure, rightCauchyGreen, pressure,
                                                                    dpdStress, dpdRCG );

        if ( error ){
            errorOut result = new errorNode( "computeReferenceSecondOrderStressPressure (second order jacobian)",
                                             "Error in computation of pressure in the reference configuration" );
            result->addNext( error );
            return result;
        }

        d2pdStressdRCG = vectorTools::eye< constantType >( referenceStressMeasure.size() ) / 3;

        return NULL;
    }

    errorOut computeDeviatoricReferenceSecondOrderStress( const variableVector &secondOrderReferenceStress,
                                                          const variableVector &rightCauchyGreenDeformation,
                                                          variableVector &deviatoricSecondOrderReferenceStress ){
        /*!
         * Compute the deviatoric part of a second order stress measure in the reference configuration.
         * \hat{S}_{IJ} = S_{IJ} - \frac{1}{3} C_{AB} S_{AB} (C^{-1})_{IJ}
         *
         * :param const variableVector &secondOrderReferenceStress: The stress measure in the reference configuration.
         * :param const variableVector &rightCauchyGreenDeformation: The right Cauchy Green Deformation tensor of the 
         *     deformation between configurations.
         * :param variableVector &deviatoricSecondOrderReferenceStrain: The deviatoric part of the second order 
         *     stress in the reference configuration.
         */

        variableType pressure;
        errorOut error = computeReferenceSecondOrderStressPressure( secondOrderReferenceStress, rightCauchyGreenDeformation, pressure );

        if ( error ){
            errorOut result = new errorNode( "computeDeviatoricReferenceSecondOrderStress",
                                             "Error in computation of the reference pressure" );
            result->addNext( error );
            return result;
        }

        return computeDeviatoricReferenceSecondOrderStress( secondOrderReferenceStress, rightCauchyGreenDeformation,
                                                            pressure, deviatoricSecondOrderReferenceStress );
    }

    errorOut computeDeviatoricReferenceSecondOrderStress( const variableVector &secondOrderReferenceStress,
                                                          const variableVector &rightCauchyGreenDeformation,
                                                          const variableType &pressure,
                                                          variableVector &deviatoricSecondOrderReferenceStress ){
        /*!
         * Compute the deviatoric part of a second order stress measure in the reference configuration.
         * \hat{S}_{IJ} = S_{IJ} - \frac{1}{3} C_{AB} S_{AB} (C^{-1})_{IJ}
         *
         * :param const variableVector &secondOrderReferenceStress: The stress measure in the reference configuration.
         * :param const variableVector &rightCauchyGreenDeformation: The right Cauchy Green Deformation tensor of the 
         *     deformation between configurations.
         * :param const variableType &pressure: The pressure computed from the reference stress.
         * :param variableVector &deviatoricSecondOrderReferenceStrain: The deviatoric part of the second order 
         *     stress in the reference configuration.
         */

        //Assume 3d
        unsigned int dim = 3;

        variableVector invRCG = vectorTools::inverse( rightCauchyGreenDeformation, dim, dim );

        deviatoricSecondOrderReferenceStress = secondOrderReferenceStress - pressure * invRCG;

        return NULL;
    }

    errorOut computeDeviatoricReferenceSecondOrderStress( const variableVector &secondOrderReferenceStress,
                                                          const variableVector &rightCauchyGreenDeformation,
                                                          variableVector &deviatoricSecondOrderReferenceStress,
                                                          variableMatrix &dDeviatoricReferenceStressdReferenceStress,
                                                          variableMatrix &dDeviatoricReferenceStressdRCG ){
        /*!
         * Compute the deviatoric part of a second order stress measure in the reference configuration.
         * \hat{S}_{IJ} = S_{IJ} - \frac{1}{3} C_{AB} S_{AB} (C^{-1})_{IJ}
         * 
         * Also compute the Jacobians.
         * \frac{ \partial \hat{S}_{IJ} }{ \partial S_{KL} } = \delta_{IK} \delta_{LJ} - \frac{1}{3} C_{KL} (C^{-1})_{IJ}
         * \frac{ \partial \hat{S}_{IJ} }{ \partial C_{KL} } = \frac{1}{3} \left( C_{AB} S_{AB} (C^{-1})_{IK} (C^{-1})_{LJ} - S_{KL} (C^{-1}_{IJ}) \right)
         *
         * :param const variableVector &secondOrderReferenceStress: The stress measure in the reference configuration.
         * :param const variableVector &rightCauchyGreenDeformation: The right Cauchy Green Deformation tensor of the 
         *     deformation between configurations.
         * :param variableVector &deviatoricSecondOrderReferenceStrain: The deviatoric part of the second order 
         *     stress in the reference configuration.
         * :param variableMatrix &dDeviatoricReferenceStressdReferenceStress: The jacobian w.r.t. the reference stress.
         * :param variableMatrix &dDeviatoricreferenceStressdRCG: The jacobian w.r.t. the right Cauchy Green deformation
         *     tensor.
         */

        variableType pressure;
        variableVector dpdS, dpdC;
        errorOut error = computeReferenceSecondOrderStressPressure( secondOrderReferenceStress, rightCauchyGreenDeformation, pressure,
                                                                    dpdS, dpdC );

        if ( error ){
            errorOut result = new errorNode( "computeDeviatoricReferenceSecondOrderStress",
                                             "Error in computation of the reference pressure" );
            result->addNext( error );
            return result;
        }

        return computeDeviatoricReferenceSecondOrderStress( secondOrderReferenceStress, rightCauchyGreenDeformation,
                                                            pressure, dpdS, dpdC,
                                                            deviatoricSecondOrderReferenceStress,
                                                            dDeviatoricReferenceStressdReferenceStress,
                                                            dDeviatoricReferenceStressdRCG );
    }

    errorOut computeDeviatoricReferenceSecondOrderStress( const variableVector &secondOrderReferenceStress,
                                                          const variableVector &rightCauchyGreenDeformation,
                                                          const variableType &pressure,
                                                          const variableVector &dPressuredStress,
                                                          const variableVector &dPressuredRCG,
                                                          variableVector &deviatoricSecondOrderReferenceStress,
                                                          variableMatrix &dDeviatoricReferenceStressdReferenceStress,
                                                          variableMatrix &dDeviatoricReferenceStressdRCG ){
        /*!
         * Compute the deviatoric part of a second order stress measure in the reference configuration.
         * \hat{S}_{IJ} = S_{IJ} - \frac{1}{3} C_{AB} S_{AB} (C^{-1})_{IJ}
         * 
         * Also compute the Jacobians.
         * \frac{ \partial \hat{S}_{IJ} }{ \partial S_{KL} } = \delta_{IK} \delta_{LJ} - \frac{1}{3} C_{KL} (C^{-1})_{IJ}
         * \frac{ \partial \hat{S}_{IJ} }{ \partial C_{KL} } = \frac{1}{3} \left( C_{AB} S_{AB} (C^{-1})_{IK} (C^{-1})_{LJ} - S_{KL} (C^{-1}_{IJ}) \right)
         *
         * :param const variableVector &secondOrderReferenceStress: The stress measure in the reference configuration.
         * :param const variableVector &rightCauchyGreenDeformation: The right Cauchy Green Deformation tensor of the 
         *     deformation between configurations.
         * :param const variableType &pressure: The pressure of the reference stress measure.
         * :param const variableVector &dPressuredStress: The Jacobian of the pressure w.r.t. the reference stress.
         * :param const variableVector &dPressuredRCG: The Jacobian of the pressure w.r.t. the right Cauchy-Green
         *     deformation measure.
         * :param variableVector &deviatoricSecondOrderReferenceStrain: The deviatoric part of the second order 
         *     stress in the reference configuration.
         * :param variableMatrix &dDeviatoricReferenceStressdReferenceStress: The jacobian w.r.t. the reference stress.
         * :param variableMatrix &dDeviatoricreferenceStressdRCG: The jacobian w.r.t. the right Cauchy Green deformation
         *     tensor.
         */
        //Assume 3d
        unsigned int dim = 3;

        variableVector invRCG = vectorTools::inverse( rightCauchyGreenDeformation, dim, dim );

        deviatoricSecondOrderReferenceStress = secondOrderReferenceStress - pressure * invRCG;

        //Compute the first order jacobians
        dDeviatoricReferenceStressdReferenceStress = vectorTools::eye< constantType >( dim * dim );
        dDeviatoricReferenceStressdReferenceStress -= vectorTools::dyadic( invRCG, dPressuredStress );

        dDeviatoricReferenceStressdRCG = - vectorTools::dyadic( invRCG, dPressuredRCG );
        for ( unsigned int I = 0; I < dim; I++ ){
            for ( unsigned int J = 0; J < dim; J++ ){
                for ( unsigned int K = 0; K < dim; K++ ){
                    for ( unsigned int L = 0; L < dim; L++ ){
                        dDeviatoricReferenceStressdRCG[ dim * I + J ][ dim * K + L ] += pressure * invRCG[ dim * I + K ] * invRCG[ dim * L + J ];
                    }
                }
            }
        }

        return NULL;
    }

    errorOut computeDeviatoricReferenceSecondOrderStress( const variableVector &secondOrderReferenceStress,
                                                          const variableVector &rightCauchyGreenDeformation,
                                                          const variableType &pressure,
                                                          const variableVector &dPressuredStress,
                                                          const variableVector &dPressuredRCG,
                                                          const variableMatrix &d2PressuredStressdRCG,
                                                          variableVector &deviatoricSecondOrderReferenceStress,
                                                          variableMatrix &dDeviatoricReferenceStressdReferenceStress,
                                                          variableMatrix &dDeviatoricReferenceStressdRCG,
                                                          variableMatrix &d2DevSdSdRCG ){
        /*!
         * Compute the deviatoric part of a second order stress measure in the reference configuration.
         * \hat{S}_{IJ} = S_{IJ} - \frac{1}{3} C_{AB} S_{AB} (C^{-1})_{IJ}
         * 
         * Also compute the Jacobians.
         * \frac{ \partial \hat{S}_{IJ} }{ \partial S_{KL} } = \delta_{IK} \delta_{LJ} - \frac{1}{3} C_{KL} (C^{-1})_{IJ}
         * \frac{ \partial \hat{S}_{IJ} }{ \partial C_{KL} } = \frac{1}{3} \left( C_{AB} S_{AB} (C^{-1})_{IK} (C^{-1})_{LJ} - S_{KL} (C^{-1}_{IJ}) \right)
         *
         * \frac{ \partial \hat{S}_{IJ} }{ \partial S_{KL} \partial C_{MN} } = - \frac{ \partial p^2 }{ \partial S_{KL} \partial C_{MN}} (C^{-1})_{IJ} + \frac{ \partial p }{ \partial S_{KL} } (C^{-1})_{IM} (C^{-1})_{NJ}
         *
         * :param const variableVector &secondOrderReferenceStress: The stress measure in the reference configuration.
         * :param const variableVector &rightCauchyGreenDeformation: The right Cauchy-Green deformation tensor of the 
         *     deformation between configurations.
         * :param const variableType &pressure: The pressure of the reference stress.
         * :param const variableVector &dPressuredStress: The Jacobian of the pressure w.r.t. the stress.
         * :param const variableVector &dPressuredRCG: The Jacobian of the pressure w.r.t. the right Cauchy-Green
         *     deformation tensor.
         * :param const variableMatrix &dPressuredStressdRCG: The Jacobian of the pressure w.r.t. the stress and 
         *     the right Cauchy-Green deformation tensor.
         * :param variableVector &deviatoricSecondOrderReferenceStrain: The deviatoric part of the second order 
         *     stress in the reference configuration.
         * :param variableMatrix &dDeviatoricReferenceStressdReferenceStress: The Jacobian w.r.t. the reference stress.
         * :param variableMatrix &dDeviatoricreferenceStressdRCG: The Jacobian w.r.t. the right Cauchy-Green deformation
         *     tensor.
         * :param variableMatrix &d2DevSdSdRCG: The second order mixed Jacobian w.r.t. the stress and right Cauchy-Green 
         *     deformation tensor. Stored [IJ][KLMN]
         */
        //Assume 3d
        unsigned int dim = 3;

        variableVector invRCG = vectorTools::inverse( rightCauchyGreenDeformation, dim, dim );

        deviatoricSecondOrderReferenceStress = secondOrderReferenceStress - pressure * invRCG;

        //Compute the first order jacobians
        dDeviatoricReferenceStressdReferenceStress = vectorTools::eye< constantType >( dim * dim );
        dDeviatoricReferenceStressdReferenceStress -= vectorTools::dyadic( invRCG, dPressuredStress );

        dDeviatoricReferenceStressdRCG = - vectorTools::dyadic( invRCG, dPressuredRCG );
        for ( unsigned int I = 0; I < dim; I++ ){
            for ( unsigned int J = 0; J < dim; J++ ){
                for ( unsigned int K = 0; K < dim; K++ ){
                    for ( unsigned int L = 0; L < dim; L++ ){
                        dDeviatoricReferenceStressdRCG[ dim * I + J ][ dim * K + L ] += pressure * invRCG[ dim * I + K ] * invRCG[ dim * L + J ];
                    }
                }
            }
        }

        //Compute the second order jacobians
        d2DevSdSdRCG = variableMatrix( dim * dim, variableVector( dim * dim * dim * dim, 0 ) );
        for ( unsigned int I = 0; I < dim; I++ ){
            for ( unsigned int J = 0; J < dim; J++ ){
                for ( unsigned int K = 0; K < dim; K++ ){
                    for ( unsigned int L = 0; L < dim; L++ ){
                        for ( unsigned int M = 0; M < dim; M++ ){
                            for ( unsigned int N = 0; N < dim; N++ ){
                                d2DevSdSdRCG[ dim * I + J ][ dim * dim * dim * K + dim * dim * L + dim * M + N ] = 
                                    -d2PressuredStressdRCG[ dim * K + L ][ dim * M + N ] * invRCG[ dim * I + J ]
                                    +dPressuredStress[ dim * K + L ] * invRCG[ dim * I + M ] * invRCG[ dim * N + J ];
                            }
                        }
                    }
                }
            }
        }

        return NULL;
    }

    errorOut computeDeviatoricReferenceSecondOrderStress( const variableVector &secondOrderReferenceStress,
                                                          const variableVector &rightCauchyGreenDeformation,
                                                          variableVector &deviatoricSecondOrderReferenceStress,
                                                          variableMatrix &dDeviatoricReferenceStressdReferenceStress,
                                                          variableMatrix &dDeviatoricReferenceStressdRCG,
                                                          variableMatrix &d2DevSdSdRCG ){
        /*!
         * Compute the deviatoric part of a second order stress measure in the reference configuration.
         * \hat{S}_{IJ} = S_{IJ} - \frac{1}{3} C_{AB} S_{AB} (C^{-1})_{IJ}
         * 
         * Also compute the Jacobians.
         * \frac{ \partial \hat{S}_{IJ} }{ \partial S_{KL} } = \delta_{IK} \delta_{LJ} - \frac{1}{3} C_{KL} (C^{-1})_{IJ}
         * \frac{ \partial \hat{S}_{IJ} }{ \partial C_{KL} } = \frac{1}{3} \left( C_{AB} S_{AB} (C^{-1})_{IK} (C^{-1})_{LJ} - S_{KL} (C^{-1}_{IJ}) \right)
         *
         * \frac{ \partial \hat{S}_{IJ} }{ \partial S_{KL} \partial C_{MN} } = - \frac{ \partial p^2 }{ \partial S_{KL} \partial C_{MN}} (C^{-1})_{IJ} + \frac{ \partial p }{ \partial S_{KL} } (C^{-1})_{IM} (C^{-1})_{NJ}
         *
         * :param const variableVector &secondOrderReferenceStress: The stress measure in the reference configuration.
         * :param const variableVector &rightCauchyGreenDeformation: The right Cauchy-Green deformation tensor of the 
         *     deformation between configurations.
         * :param variableVector &deviatoricSecondOrderReferenceStrain: The deviatoric part of the second order 
         *     stress in the reference configuration.
         * :param variableMatrix &dDeviatoricReferenceStressdReferenceStress: The Jacobian w.r.t. the reference stress.
         * :param variableMatrix &dDeviatoricreferenceStressdRCG: The Jacobian w.r.t. the right Cauchy-Green deformation
         *     tensor.
         * :param variableMatrix &d2DevSdSdRCG: The second order mixed Jacobian w.r.t. the stress and right Cauchy-Green 
         *     deformation tensor. Stored [IJ][KLMN]
         */

        variableType pressure;
        variableVector dpdS, dpdC;
        variableMatrix d2pdSdC;
        errorOut error = computeReferenceSecondOrderStressPressure( secondOrderReferenceStress, rightCauchyGreenDeformation, pressure,
                                                                    dpdS, dpdC, d2pdSdC );

        if ( error ){
            errorOut result = new errorNode( "computeDeviatoricReferenceSecondOrderStress",
                                             "Error in computation of the reference pressure" );
            result->addNext( error );
            return result;
        }

        return computeDeviatoricReferenceSecondOrderStress( secondOrderReferenceStress, rightCauchyGreenDeformation,
                                                            pressure, dpdS, dpdC, d2pdSdC,
                                                            deviatoricSecondOrderReferenceStress,
                                                            dDeviatoricReferenceStressdReferenceStress,
                                                            dDeviatoricReferenceStressdRCG,
                                                            d2DevSdSdRCG );
    }

    errorOut computeSecondOrderReferenceStressDecomposition( const variableVector &secondOrderReferenceStress,
                                                             const variableVector &rightCauchyGreenDeformation,
                                                             variableVector &deviatoricSecondOrderReferenceStress,
                                                             variableType &pressure ){
        /*!
         * Compute the decomposition of a second-order stress measure into pressure and deviatoric parts.
         *
         * :param const variableVector &secondOrderReferenceStress: The second-order stress in the reference
         *     configuration.
         * :param const variableVector &rightCauchyGreenDeformation: The right Cauchy-Green deformation tensor 
         *     between the current configuration and the configuration the stress is located in.
         * :param variableVector &deviatoricSecondOrderReferenceStress: The deviatoric part of the second order 
         *     stress tensor.
         * :param variableType &pressure: The pressure of the stress tensor.
         */

        errorOut error = computeReferenceSecondOrderStressPressure( secondOrderReferenceStress,
                                                                    rightCauchyGreenDeformation,
                                                                    pressure );

        if ( error ){
            errorOut result = new errorNode( "computeSecondOrderReferenceStressDecomposition",
                                             "Error in computation of pressure" );
            result->addNext( error );
            return result;
        }

        error = computeDeviatoricReferenceSecondOrderStress( secondOrderReferenceStress,
                                                             rightCauchyGreenDeformation,
                                                             pressure,
                                                             deviatoricSecondOrderReferenceStress );

        if ( error ){
            errorOut result = new errorNode( "computeSecondOrderReferenceStressDecomposition",
                                             "Error in computation of the deviatoric part of the stress" );
            result->addNext( error );
            return result;
        }

        return NULL;
    }

    errorOut computeSecondOrderReferenceStressDecomposition( const variableVector &secondOrderReferenceStress,
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
         * :param const variableVector &secondOrderReferenceStress: The second-order stress in the reference
         *     configuration.
         * :param const variableVector &rightCauchyGreenDeformation: The right Cauchy-Green deformation tensor 
         *     between the current configuration and the configuration the stress is located in.
         * :param variableVector &deviatoricSecondOrderReferenceStress: The deviatoric part of the second order 
         *     stress tensor.
         * :param variableType &pressure: The pressure of the stress tensor.
         * :param variableMatrix &dDevStressdStress: The Jacobian of the deviatoric part of the reference 
         *     stress w.r.t. the reference stress.
         * :param variableMatrix &dDevStressdRCG: The Jacobian of the deviatoric part of the reference stress
         *     w.r.t. the right Cauchy-Green deformation tensor.
         * :param variableMatrix &dPressuredStress: The Jacobian of the pressure w.r.t. the reference stress.
         * :param variableMatrix &dPressuredRCG: The Jacobian of the pressure w.r.t. the right Cauchy-Green 
         *     deformation tensor.
         */

        errorOut error = computeReferenceSecondOrderStressPressure( secondOrderReferenceStress,
                                                                    rightCauchyGreenDeformation,
                                                                    pressure, dPressuredStress, dPressuredRCG );

        if ( error ){
            errorOut result = new errorNode( "computeSecondOrderReferenceStressDecomposition (jacobian)",
                                             "Error in computation of pressure" );
            result->addNext( error );
            return result;
        }

        error = computeDeviatoricReferenceSecondOrderStress( secondOrderReferenceStress,
                                                             rightCauchyGreenDeformation,
                                                             pressure, dPressuredStress, dPressuredRCG,
                                                             deviatoricSecondOrderReferenceStress,
                                                             dDevStressdStress, dDevStressdRCG );

        if ( error ){
            errorOut result = new errorNode( "computeSecondOrderReferenceStressDecomposition (jacobian)",
                                             "Error in computation of the deviatoric part of the stress" );
            result->addNext( error );
            return result;
        }

        return NULL;
    }

    errorOut computeSecondOrderReferenceStressDecomposition( const variableVector &secondOrderReferenceStress,
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
         * :param const variableVector &secondOrderReferenceStress: The second-order stress in the reference
         *     configuration.
         * :param const variableVector &rightCauchyGreenDeformation: The right Cauchy-Green deformation tensor 
         *     between the current configuration and the configuration the stress is located in.
         * :param variableVector &deviatoricSecondOrderReferenceStress: The deviatoric part of the second order 
         *     stress tensor.
         * :param variableType &pressure: The pressure of the stress tensor.
         * :param variableMatrix &dDevStressdStress: The Jacobian of the deviatoric part of the reference 
         *     stress w.r.t. the reference stress.
         * :param variableMatrix &dDevStressdRCG: The Jacobian of the deviatoric part of the reference stress
         *     w.r.t. the right Cauchy-Green deformation tensor.
         * :param variableMatrix &dPressuredStress: The Jacobian of the pressure w.r.t. the reference stress.
         * :param variableMatrix &dPressuredRCG: The Jacobian of the pressure w.r.t. the right Cauchy-Green 
         *     deformation tensor.
         * :param variableMatrix &d2DevStressdStressdRCG: The second order Jacobian of the deviatoric part of the 
         *     reference stress w.r.t. the reference stress and the right Cauchy-Green deformation tensor. 
         *     $\frac{ \partial^2 dev( S )_{IJ} }{ \partial S_{KL} \partial C_{MN} }$
         *     stored as [IJ][KLMN]
         * :param variableMatrix &d2PressuredStressdRCG: THe second order Jacobian of the deviatoric part of the 
         *     reference stress w.r.t. the reference stress and the right Cauchy-Green deformation tensor.
         */

        errorOut error = computeReferenceSecondOrderStressPressure( secondOrderReferenceStress,
                                                                    rightCauchyGreenDeformation,
                                                                    pressure, dPressuredStress, dPressuredRCG,
                                                                    d2PressuredStressdRCG );

        if ( error ){
            errorOut result = new errorNode( "computeSecondOrderReferenceStressDecomposition (second order jacobian)",
                                             "Error in computation of pressure" );
            result->addNext( error );
            return result;
        }

        error = computeDeviatoricReferenceSecondOrderStress( secondOrderReferenceStress,
                                                             rightCauchyGreenDeformation,
                                                             pressure, dPressuredStress, dPressuredRCG,
                                                             d2PressuredStressdRCG,
                                                             deviatoricSecondOrderReferenceStress,
                                                             dDevStressdStress, dDevStressdRCG,
                                                             d2DevStressdStressdRCG );

        if ( error ){
            errorOut result = new errorNode( "computeSecondOrderReferenceStressDecomposition ( second order jacobian)",
                                             "Error in computation of the deviatoric part of the stress" );
            result->addNext( error );
            return result;
        }

        return NULL;
    }

    errorOut computeHigherOrderReferenceStressDecomposition( const variableVector &higherOrderReferenceStress,
                                                             const variableVector &rightCauchyGreenDeformation,
                                                             variableVector &deviatoricHigherOrderReferenceStress,
                                                             variableVector &pressure ){
        /*!
         * Compute the decomposition of the higher-order stress measure into pressure and deviatoric parts.
         *
         * :param const variableVector &higherOrderReferenceStress: The higher order stress in the reference
         *     configuration.
         * :param const variableVector &rightCauchyGreenDeformation: The right Cauchy-Green deformation tensor 
         *     between the current configuration and the configuration the stress is located in.
         * :param variableVector &deviatoricHigherOrderReferenceStress: The deviatoric part of the higher order 
         *     stress tensor.
         * :param variableVector &pressure: The pressure of the higher order stress tensor.
         */

        errorOut error = computeReferenceHigherOrderStressPressure( higherOrderReferenceStress,
                                                                    rightCauchyGreenDeformation,
                                                                    pressure );

        if ( error ){
            errorOut result = new errorNode( "computeHigherOrderReferenceStressDecomposition",
                                             "Error in computation of pressure" );
            result->addNext( error );
            return result;
        }

        error = computeDeviatoricReferenceHigherOrderStress( higherOrderReferenceStress,
                                                             rightCauchyGreenDeformation,
                                                             pressure,
                                                             deviatoricHigherOrderReferenceStress );

        if ( error ){
            errorOut result = new errorNode( "computeHigherOrderReferenceStressDecomposition",
                                             "Error in computation of the deviatoric part of the stress" );
            result->addNext( error );
            return result;
        }

        return NULL;
    }

    errorOut computeHigherOrderReferenceStressDecomposition( const variableVector &higherOrderReferenceStress,
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
         * :param const variableVector &higherOrderReferenceStress: The higher order stress in the reference
         *     configuration.
         * :param const variableVector &rightCauchyGreenDeformation: The right Cauchy-Green deformation tensor 
         *     between the current configuration and the configuration the stress is located in.
         * :param variableVector &deviatoricHigherOrderReferenceStress: The deviatoric part of the higher order 
         *     stress tensor.
         * :param variableVector &pressure: The pressure of the higher order stress tensor.
         * :param variableMatrix &dDevStressdStress: The Jacobian of the deviatoric stress w.r.t. the stress.
         * :param variableMatrix &dDevStressdRCG: The Jacobian of the deviatoric stress w.r.t. the right 
         *     Cauchy-Green deformation tensor.
         * :param variableMatrix &dPressuredStress: The Jacobian of the pressure w.r.t. the stress.
         * :param variableMatrix &dPressuredRCG: The Jacobian of the pressure w.r.t. the right 
         *     Cauchy-Green deformation tensor.
         */

        errorOut error = computeReferenceHigherOrderStressPressure( higherOrderReferenceStress,
                                                                    rightCauchyGreenDeformation,
                                                                    pressure, dPressuredStress,
                                                                    dPressuredRCG );

        if ( error ){
            errorOut result = new errorNode( "computeHigherOrderReferenceStressDecomposition (jacobian)",
                                             "Error in computation of pressure" );
            result->addNext( error );
            return result;
        }

        error = computeDeviatoricReferenceHigherOrderStress( higherOrderReferenceStress,
                                                             rightCauchyGreenDeformation,
                                                             pressure, dPressuredStress,
                                                             dPressuredRCG,
                                                             deviatoricHigherOrderReferenceStress,
                                                             dDevStressdStress, dDevStressdRCG );

        if ( error ){
            errorOut result = new errorNode( "computeHigherOrderReferenceStressDecomposition (jacobian)",
                                             "Error in computation of the deviatoric part of the stress" );
            result->addNext( error );
            return result;
        }

        return NULL;
    }

    errorOut computeHigherOrderReferenceStressDecomposition( const variableVector &higherOrderReferenceStress,
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
         * :param const variableVector &higherOrderReferenceStress: The higher order stress in the reference
         *     configuration.
         * :param const variableVector &rightCauchyGreenDeformation: The right Cauchy-Green deformation tensor 
         *     between the current configuration and the configuration the stress is located in.
         * :param variableVector &deviatoricHigherOrderReferenceStress: The deviatoric part of the higher order 
         *     stress tensor.
         * :param variableVector &pressure: The pressure of the higher order stress tensor.
         * :param variableMatrix &dDevStressdStress: The Jacobian of the deviatoric stress w.r.t. the stress.
         * :param variableMatrix &dDevStressdRCG: The Jacobian of the deviatoric stress w.r.t. the right 
         *     Cauchy-Green deformation tensor.
         * :param variableMatrix &dPressuredStress: The Jacobian of the pressure w.r.t. the stress.
         * :param variableMatrix &dPressuredRCG: The Jacobian of the pressure w.r.t. the right 
         *     Cauchy-Green deformation tensor.
         * :param variableMatrix &d2DevStressdStressdRCG: The second order jacobian of the deviatoric stress 
         *     w.r.t. the stress and the right Cauchy-Green deformation tensor.
         * :param variableMatrix &d2PressuredStressdRCG: The second order jacobian of the pressure
         *     w.r.t. the stress and the right Cauchy-Green deformation tensor.
         */

        errorOut error = computeReferenceHigherOrderStressPressure( higherOrderReferenceStress,
                                                                    rightCauchyGreenDeformation,
                                                                    pressure, dPressuredStress,
                                                                    dPressuredRCG, d2PressuredStressdRCG );

        if ( error ){
            errorOut result = new errorNode( "computeHigherOrderReferenceStressDecomposition (second order jacobian)",
                                             "Error in computation of pressure" );
            result->addNext( error );
            return result;
        }

        error = computeDeviatoricReferenceHigherOrderStress( higherOrderReferenceStress,
                                                             rightCauchyGreenDeformation,
                                                             pressure, dPressuredStress,
                                                             dPressuredRCG, d2PressuredStressdRCG,
                                                             deviatoricHigherOrderReferenceStress,
                                                             dDevStressdStress, dDevStressdRCG,
                                                             d2DevStressdStressdRCG );

        if ( error ){
            errorOut result = new errorNode( "computeHigherOrderReferenceStressDecomposition (second order jacobian)",
                                             "Error in computation of the deviatoric part of the stress" );
            result->addNext( error );
            return result;
        }

        return NULL;
    }

    errorOut computeHigherOrderStressNorm( const variableVector &higherOrderStress, variableVector &higherOrderStressNorm ){
        /*!
         * Compute the norm of the higher order stress which is defined as
         * || M ||_K = \sqrt{ M_{IJK} M_{IJK} }
         *
         * where K is not summed over.
         *
         * :param const variableVector &higherOrderStress: The higher order stress tensor.
         * :param variableVector &higherOrderStressNorm: The norm of the higher order stress.
         */

        //Assume 3D
        unsigned int dim = 3;

        if ( higherOrderStress.size() != dim * dim * dim ){
            return new errorNode( "computeHigherOrderStressNorm",
                                  "The higher order stress does not have the expected dimension" );
        }
        
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

        return NULL;
    }

    errorOut computeHigherOrderStressNorm( const variableVector &higherOrderStress, variableVector &higherOrderStressNorm,
                                           variableMatrix &dHigherOrderStressNormdHigherOrderStress ){
        /*!
         * Compute the norm of the higher order stress which is defined as
         * || M ||_K = \sqrt{ M_{IJK} M_{IJK} }
         *
         * where K is not summed over.
         *
         * Also computes the Jacobians
         *
         * \frac{ \partial || M ||_K }{ \partial M_{LMN} } = \frac{ M_{LMK} \delta_{KN} }{ || M ||_K }
         *
         * where K is not summed over
         *
         * :param const variableVector &higherOrderStress: The higher order stress tensor.
         * :param variableVector &higherOrderStressNorm: The norm of the higher order stress.
         * :param variableMatrix &dHigherOrderSTressNormdHigherOrderStress: The Jacobian of the higher order
         *     stress norm w.r.t. the higher order stress.
         */

        //Assume 3D
        unsigned int dim = 3;

        errorOut error = computeHigherOrderStressNorm( higherOrderStress, higherOrderStressNorm );

        if ( error ){
            errorOut result = new errorNode( "computeHigherOrderStressNorm (jacobian)",
                                             "Error in computation of higher order stress norm" );
            result->addNext( error );
            return result;
        }

        constantVector eye( dim * dim );
        vectorTools::eye( eye );
        dHigherOrderStressNormdHigherOrderStress = variableMatrix( dim, variableVector( dim * dim * dim, 0 ) );
        for ( unsigned int K = 0; K < 3; K++ ){
            for ( unsigned int L = 0; L < 3; L++ ){
                for ( unsigned int M = 0; M < 3; M++ ){
                    for ( unsigned int N = 0; N < 3; N++ ){
                        dHigherOrderStressNormdHigherOrderStress[ K ][ dim * dim * L + dim * M + N ] = higherOrderStress[ dim * dim * L + dim * M + K ]  * eye[ dim * K + N ] / higherOrderStressNorm[ K ];
                    }
                }
            }
        }

        return NULL;
    }

    errorOut computeHigherOrderStressNorm( const variableVector &higherOrderStress, variableVector &higherOrderStressNorm,
                                           variableMatrix &dHigherOrderStressNormdHigherOrderStress,
                                           variableMatrix &d2HigherOrderStressNormdHigherOrderStress2 ){
        /*!
         * Compute the norm of the higher order stress which is defined as
         * || M ||_K = \sqrt{ M_{IJK} M_{IJK} }
         *
         * where K is not summed over.
         *
         * Also computes the Jacobians
         *
         * \frac{ \partial || M ||_K }{ \partial M_{LMN} } = \frac{ M_{LMK} \delta_{KN} }{ || M ||_K }
         * \frac{ \partial^2 || M ||_K }{ \partial M_{LMN} \partial M_{OPQ} } = \frac{1}{ || M ||_K } \left[ \delta_{LO} \delta_{MP} \delta_{KQ} \delta_{KN} - \frac{ M_{LMK} \delta_{KN} }{ || M ||_K } \frac{ M_{OPK} \delta_{KQ} }{ || M ||_K } \right]
         *
         * where K is not summed over
         *
         * :param const variableVector &higherOrderStress: The higher order stress tensor.
         * :param variableVector &higherOrderStressNorm: The norm of the higher order stress.
         * :param variableMatrix &dHigherOrderSTressNormdHigherOrderStress: The Jacobian of the higher order
         *     stress norm w.r.t. the higher order stress.
         * :param variableMatrix &d2HigherOrderStressNormdHigherOrderStress2: The second order Jacobian of the 
         *     higher order stress norm w.r.t. the higher order stress.
         */

        //Assume 3D
        unsigned int dim = 3;

        errorOut error = computeHigherOrderStressNorm( higherOrderStress, higherOrderStressNorm, dHigherOrderStressNormdHigherOrderStress );

        if ( error ){
            errorOut result = new errorNode( "computeHigherOrderStressNorm (second order jacobian)",
                                             "Error in computation of higher order stress norm" );
            result->addNext( error );
            return result;
        }

        constantVector eye( dim * dim );
        vectorTools::eye( eye );

        d2HigherOrderStressNormdHigherOrderStress2 = variableMatrix( dim, variableVector( dim * dim * dim * dim * dim * dim, 0 ) );

        for ( unsigned int K = 0; K < dim; K++ ){
            for ( unsigned int L = 0; L < dim; L++ ){
                for ( unsigned int M = 0; M < dim; M++ ){
                    for ( unsigned int N = 0; N < dim; N++ ){
                        for ( unsigned int O = 0; O < dim; O++ ){
                            for ( unsigned int P = 0; P < dim; P++ ){
                                for ( unsigned int Q = 0; Q < dim; Q++ ){
                                    d2HigherOrderStressNormdHigherOrderStress2[ K ][ dim * dim * dim * dim * dim * L + dim * dim * dim * dim * M + dim * dim * dim * N + dim * dim * O + dim * P + Q ]
                                        = ( eye[ dim * L + O ] * eye[ dim * M + P ] * eye[ dim * K + Q ] * eye[ dim * K + N ]
                                        -   dHigherOrderStressNormdHigherOrderStress[ K ][ dim * dim * L + dim * M + N ]
                                        *   dHigherOrderStressNormdHigherOrderStress[ K ][ dim * dim * O + dim * P + Q ] )
                                        / higherOrderStressNorm[ K ];
                                }
                            }
                        }
                    }
                }
            }
        }

        return NULL;
    }

    errorOut assembleDeformationGradient( const variableMatrix &displacementGradient, variableVector &deformationGradient ){
        /*!
         * Assemble the deformation gradient from the gradient of the displacement.
         *
         * :param const variableMatrix &displacementGradient: The gradient of the displacement.
         * :param variableVector &deformationGradient: The deformation gradient.
         */

        //Assume 3D
        unsigned int dim = 3;

        if ( displacementGradient.size() != dim ){
            return new errorNode( "assembleDeformationGradient",
                                  "The gradient of the deformation is not 3D" );
        }

        for ( unsigned int i = 0; i < dim; i++ ){
            if ( displacementGradient[ i ].size() != dim ){
                return new errorNode( "assembleDeformationGradient",
                                      "The gradient of the deformation is not a regular matrix" );
            }
        }

        deformationGradient = vectorTools::appendVectors( displacementGradient + vectorTools::eye< variableType >( dim ) );

        return NULL;
    }

    errorOut assembleDeformationGradient( const variableMatrix &displacementGradient, variableVector &deformationGradient,
                                          variableMatrix &dFdGradU ){
        /*!
         * Assemble the deformation gradient from the gradient of the displacement.
         *
         * :param const variableMatrix &displacementGradient: The gradient of the displacement.
         * :param variableVector &deformationGradient: The deformation gradient.
         * :param variableMatrix &dFdGradU: The Jacobian of the deformation gradient w.r.t. the displacement gradient.
         */

        //Assume 3D
        unsigned int dim = 3;

        errorOut error = assembleDeformationGradient( displacementGradient, deformationGradient );

        if ( error ){
            errorOut result = new errorNode( "assembleDeformationGradient (jacobian)",
                                             "Error in the computation of the deformation gradient" );
            result->addNext( error );
            return result;
        }

        dFdGradU = vectorTools::eye< variableType >( dim * dim );

        return NULL;
    }

    errorOut assembleMicroDeformation( const variableVector &microDisplacement, variableVector &microDeformation ){
        /*!
         * Assemble the micro deformation from the micro displacement
         *
         * :param const variableVector &microDisplacement: The micro degrees of freedom.
         * :param variableVector &microDeformation: The micro deformation.
         */

        //Assume 3D
        unsigned int dim = 3;

        if ( microDisplacement.size() != dim * dim ){
            return new errorNode( "assembleMicroDeformation",
                                  "The micro degrees of freedom must be 3D" );
        }

        constantVector eye( dim * dim );
        vectorTools::eye( eye );

        microDeformation = microDisplacement + eye;

        return NULL;
    }

    errorOut assembleMicroDeformation( const variableVector &microDisplacement, variableVector &microDeformation,
                                       variableMatrix &dChidPhi ){
        /*!
         * Assemble the micro deformation from the micro displacement
         *
         * :param const variableVector &microDisplacement: The micro degrees of freedom.
         * :param variableVector &microDeformation: The micro deformation.
         * :param variableMatrix &dChidPhi: The Jacobian of the micro deformation w.r.t. the micro displacement
         */

        //Assume 3D
        unsigned int dim = 3;

        errorOut error = assembleMicroDeformation( microDisplacement, microDeformation );

        if ( error ){
            errorOut result = new errorNode( "assembleMicroDeformation (jacobian)",
                                             "Errir in the computation of the micro deformation" );
            result->addNext( error );
            return result;
        }

        dChidPhi = vectorTools::eye< variableType >( dim * dim );

        return NULL;
    }

    errorOut assembleGradientMicroDeformation( const variableMatrix &gradientMicroDisplacement,
                                               variableVector &gradientMicroDeformation ){
        /*!
         * Assemble the gradient of the micro deformation from the gradient of the micro displacement
         * in the reference configuration
         *
         * :param const variableVector &gradientMicroDisplacement: The gradient of the micro displacement
         * :param variableVector &gradientMicroDeformation: The gradient of the micro deformation.
         */

        //Assume 3D
        unsigned int dim = 3;

        if ( gradientMicroDisplacement.size() != dim * dim ){
            return new errorNode( "assembleGradientMicroDeformation",
                                  "The gradient of the micro displacement must be 3D" );
        }

        for ( unsigned int i = 0; i < dim; i++ ){
            if ( gradientMicroDisplacement[ i ].size() != dim ){
                return new errorNode( "assembleGradientMicroDeformation",
                                      "The gradient of the micro displacement must be 3D" );
            }
        }

        gradientMicroDeformation = variableVector( dim * dim * dim, 0 );

        for ( unsigned int i = 0; i < dim; i++ ){
            for ( unsigned int I = 0; I < dim; I++ ){
                for ( unsigned int J = 0; J < dim; J++ ){
                    gradientMicroDeformation[ dim * dim * i + dim * I + J ] = gradientMicroDisplacement[ dim * i + I ][ J ];
                }
            }
        }

        return NULL;
    } 

    errorOut assembleGradientMicroDeformation( const variableMatrix &gradientMicroDisplacement,
                                               variableVector &gradientMicroDeformation,
                                               variableMatrix &dGradChidGradPhi ){
        /*!
         * Assemble the gradient of the micro deformation from the gradient of the micro displacement
         * in the reference configuration
         *
         * :param const variableVector &gradientMicroDisplacement: The gradient of the micro displacement
         * :param variableVector &gradientMicroDeformation: The gradient of the micro deformation.
         * :param variableMatrix &dGradChidGradPhi: The gradient of the micro deformation gradient w.r.t.
         *     the micro displacement gradient.
         */

        //Assume 3D
        unsigned int dim = 3;

        errorOut error = assembleGradientMicroDeformation( gradientMicroDisplacement, gradientMicroDeformation );

        if ( error ){
            errorOut result = new errorNode( "assembleGradientMicroDeformation (jacobian)",
                                             "Error in the computation of the micro deformation gradient" );
            result->addNext( error );
            return result;
        }
        
        dGradChidGradPhi = vectorTools::eye< variableType >( dim * dim * dim );

        return NULL;
    }
}
