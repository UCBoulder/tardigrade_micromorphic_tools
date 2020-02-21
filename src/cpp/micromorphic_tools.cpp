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
         * \Psi_{IJ} = F_{i I} \Xi_{i J}
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
                         variableVector &Psi, variableMatrix &dPsidF, variableMatrix &dPsidXi ){
        /*!
         * Computes the micromorphic quantity psi defined as:
         * \Psi_{IJ} = F_{i I} \Xi_{i J}
         *
         * along with the jacobians
         *
         * \frac{ \partial \Psi_{IJ} }{ \partial F_{k K} } = \delta_{I K} \Xi_{k J}
         * \frac{ \partial \Psi_{IJ} }{ \partial \Xi_{k K} } = F_{k I} \delta_{J K}
         *
         * :param const variableVector &deformationGradient: The deformation gradient.
         * :param const variableVector &microDeformation: The micro-deformation.
         * :param variableVector &Psi: The micro-displacement metric Psi
         * :param variableMatrix &dPsidF: The jacobian of Psi w.r.t. the deformation gradient.
         * :param variableMatrix &dPsidXi: The jacobian of Psi w.r.t. the micro-deformation.
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
        dPsidXi = variableMatrix( Psi.size(), variableVector( microDeformation.size(), 0 ) );

        for ( unsigned int I = 0; I < dim; I++ ){
            for ( unsigned int J = 0; J < dim; J++ ){
                for ( unsigned int k = 0; k < dim; k++ ){
                    for ( unsigned int K = 0; K < dim; K++ ){
                        dPsidF[ dim * I + J ][ dim * k + K ] = eye[ dim * I + K ] * microDeformation[ dim * k + J ];
                        dPsidXi[ dim * I + J ][ dim * k + K ] = deformationGradient[ dim * k + I ] * eye[ dim * J + K ];
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

        error = computeMicroStrain( Psi, microStrain );

        if (error){
            errorOut result = new errorNode( "computeMicroStrain (jacobian)", "Error in computation of micro-strain" );
            result->addNext( error );
            return result;
        }

        dMicroStraindPsi = vectorTools::eye< variableType >( Psi.size() );

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
         * :param const variableType &detF: The determinant of the deformation gradient.
         * :param variableVector &microStress: The micro-stress in the current 
         *     configuration.
         */

        //Assume 3d
        unsigned int dim = 3;

        microStress = variableVector( dim * dim, 0 );

        detF = vectorTools::determinant( deformationGradient, dim, dim );

        for ( unsigned int i = 0; i < dim; i++ ){
            for ( unsigned int j = 0; j < dim; j++ ){
                for ( unsigned int I = 0; I < dim; I++ ){
                    for ( unsigned int J = 0; J < dim; J++ ){
                        microStress[ dim * i + j ] += deformationGradient[ dim * i + I ]
                                                    * referenceMicroStress[ dim * I + J ]
                                                    * deformationGradient[ dim * j + J ] / detF;
                    }
                }
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
         * \frac{ \partial s_{ij} }{\partial F_{kK} } = - ( 1 / J ) \frac{\partial J}{\partial F_{k K} } s_{ij} 
         *     + ( 1 / J ) \delta_{ik} \Sigma_{KJ} F_{jJ} + ( 1 / J ) F_{iI} \Sigma_{IK} \delta_{jk}
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
        error = pushForwardReferenceMicroStress( referenceMicroStress, deformationGradient,
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

        for ( unsigned int i = 0; i < dim; i++ ){
            for ( unsigned int j = 0; j < dim; j++ ){
                for ( unsigned int k = 0; k < dim; k++ ){
                    for ( unsigned int K = 0; K < dim; K++ ){
                        dMicroStressdReferenceMicroStress[ dim * i + j ][ dim * k + K ] = deformationGradient[ dim * i + k ]
                                                                                        * deformationGradient[ dim * j + K ] / detF;
                        
                        for ( unsigned int I = 0; I < dim; I++ ){

                            dMicroStressdDeformationGradient[ dim * i + j ][ dim * k + K] = ( eye[ dim * i + k ] * referenceMicroStress[ dim * K + I ] * deformationGradient[ dim * j + I ]
                                                                                          + deformationGradient[ dim * i + I ] * referenceMicroStress[ dim * I + K ] * eye[ dim * j + k ]
                                                                                          - microStress[ dim * i + j ] * dDetFdF[ dim * k + K ] ) / detF;
                        }
                    }
                }
            }
        }

        return NULL;
    }

}
