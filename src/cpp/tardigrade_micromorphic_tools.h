/*!
 * Micromorphic Tools.h
 * ====================
 *
 * A collection of tools and utilities which seek to make the 
 * implementation of micromorphic continuum mechanics based constitutive 
 * theories easier.
 *
 */

#ifndef TARDIGRADE_MICROMORPHIC_TOOLS_H
#define TARDIGRADE_MICROMORPHIC_TOOLS_H

#include<tardigrade_error_tools.h>
#define USE_EIGEN
#include<tardigrade_vector_tools.h>
#include<tardigrade_constitutive_tools.h>

namespace tardigradeMicromorphicTools{

    typedef double variableType;
    typedef std::vector< variableType > variableVector;
    typedef std::vector< variableVector > variableMatrix;

    typedef double parameterType;
    typedef std::vector< parameterType > parameterVector;
    typedef std::vector< parameterVector > parameterMatrix;

    typedef double constantType;
    typedef std::vector< constantType > constantVector;
    typedef std::vector< constantVector > constantMatrix;

    typedef tardigradeErrorTools::Node errorNode;
    typedef errorNode* errorOut;

    errorOut computePsi( const variableVector &deformationGradient, const variableVector &microDeformation,
                         variableVector &Psi );

    errorOut computePsi( const variableVector &deformationGradient, const variableVector &microDeformation,
                         variableVector &Psi, variableMatrix &dPsidF, variableMatrix &dPsidChi );

    errorOut computeGamma( const variableVector &deformationGradient, const variableVector &gradChi,
                           variableVector &Gamma );

    errorOut computeGamma( const variableVector &deformationGradient, const variableVector &gradChi,
                           variableVector &Gamma, variableMatrix &dGammadF, variableMatrix &dGammadGradChi );

    errorOut computeMicroStrain( const variableVector &Psi, variableVector &microStrain );

    errorOut computeMicroStrain( const variableVector &Psi, variableVector &microStrain,
                                 variableMatrix &dMicroStraindPsi );


    errorOut pushForwardPK2Stress( const variableVector &PK2Stress,
                                   const variableVector &deformationGradient,
                                   variableVector &cauchyStress );

    errorOut pushForwardPK2Stress( const variableVector &PK2Stress,
                                   const variableVector &deformationGradient,
                                   variableVector &cauchyStress, 
                                   variableVector &dcauchyStressdPK2Stress,
                                   variableVector &dcauchyStressdDeformationGradient );

    errorOut pushForwardPK2Stress( const variableVector &PK2Stress,
                                   const variableVector &deformationGradient,
                                   variableVector &cauchyStress, 
                                   variableMatrix &dcauchyStressdPK2Stress,
                                   variableMatrix &dcauchyStressdDeformationGradient );

    errorOut pullBackCauchyStress( const variableVector &cauchyStress,
                                   const variableVector &deformationGradient,
                                   variableVector &PK2Stress );

    errorOut pullBackCauchyStress( const variableVector &cauchyStress,
                                   const variableVector &deformationGradient,
                                   variableVector &PK2Stress, 
                                   variableMatrix &dcauchyStressdPK2Stress,
                                   variableMatrix &dcauchyStressdDeformationGradient );

    errorOut pushForwardReferenceMicroStress( const variableVector &referenceMicroStress,
                                              const variableVector &deformationGradient,
                                              variableType &detF, variableVector &microStress );

    errorOut pushForwardReferenceMicroStress( const variableVector &referenceMicroStress,
                                              const variableVector &deformationGradient,
                                              variableVector &microStress );

    errorOut pushForwardReferenceMicroStress( const variableVector &referenceMicroStress,
                                              const variableVector &deformationGradient,
                                              variableVector &microStress, 
                                              variableVector &dMicroStressdReferenceMicroStress,
                                              variableVector &dMicroStressdDeformationGradient );

    errorOut pushForwardReferenceMicroStress( const variableVector &referenceMicroStress,
                                              const variableVector &deformationGradient,
                                              variableVector &microStress, 
                                              variableMatrix &dMicroStressdReferenceMicroStress,
                                              variableMatrix &dMicroStressdDeformationGradient );

    errorOut pullBackMicroStress( const variableVector &microStress,
                                  const variableVector &deformationGradient,
                                  variableType &detF, variableVector &inverseDeformationGradient,
                                  variableVector &referenceMicroStress );

    errorOut pullBackMicroStress( const variableVector &microStress,
                                  const variableVector &deformationGradient,
                                  variableVector &referenceMicroStress );

    errorOut pullBackMicroStress( const variableVector &microStress,
                                  const variableVector &deformationGradient,
                                  variableVector &referenceMicroStress, 
                                  variableMatrix &dReferenceMicroStressdMicroStress,
                                  variableMatrix &dReferenceMicroStressdDeformationGradient );

    errorOut pushForwardHigherOrderStress( const variableVector &referenceHigherOrderStress,
                                           const variableVector &deformationGradient,
                                           const variableVector &microDeformation,
                                           variableVector &higherOrderStress );

    errorOut pushForwardHigherOrderStress( const variableVector &referenceHigherOrderStress,
                                           const variableVector &deformationGradient,
                                           const variableVector &microDeformation,
                                           variableType &detF,
                                           variableVector &higherOrderStress );

    errorOut pushForwardHigherOrderStress( const variableVector &referenceHigherOrderStress,
                                           const variableVector &deformationGradient,
                                           const variableVector &microDeformation,
                                           variableVector &higherOrderStress,
                                           variableMatrix &dHigherOrderStressdReferenceHigherOrderStress,
                                           variableMatrix &dHigherOrderStressdDeformationGradient,
                                           variableMatrix &dHigherOrderStressdMicroDeformation );

    errorOut pullBackHigherOrderStress( const variableVector &higherOrderStress,
                                        const variableVector &deformationGradient,
                                        const variableVector &microDeformation,
                                        variableVector &referenceHigherOrderStress );

    errorOut pullBackHigherOrderStress( const variableVector &higherOrderStress,
                                        const variableVector &deformationGradient,
                                        const variableVector &microDeformation,
                                        variableType &detF, variableVector &inverseDeformationGradient,
                                        variableVector &inverseMicroDeformation,
                                        variableVector &referenceHigherOrderStress );

    errorOut pullBackHigherOrderStress( const variableVector &higherOrderStress,
                                        const variableVector &deformationGradient,
                                        const variableVector &microDeformation,
                                        variableVector &referenceHigherOrderStress,
                                        variableMatrix &dHigherOrderStressdReferenceHigherOrderStress,
                                        variableMatrix &dHigherOrderStressdDeformationGradient,
                                        variableMatrix &dHigherOrderStressdMicroDeformation );

    errorOut computeDeviatoricSecondOrderStress( const variableVector &secondOrderStress,
                                                 variableVector &deviatoricSecondOrderStress );

    errorOut computeDeviatoricSecondOrderStress( const variableVector &secondOrderStress,
                                                 variableVector &deviatoricSecondOrderStress,
                                                 variableMatrix &dDeviatoricStressdStress );

    errorOut computeReferenceSecondOrderStressPressure( const variableVector &referenceStressMeasure,
                                                        const variableVector &rightCauchyGreen,  variableType &pressure );

    errorOut computeReferenceSecondOrderStressPressure( const variableVector &referenceStressMeasure,
                                                        const variableVector &rightCauchyGreen, variableType &pressure,
                                                        variableVector &dpdStress, variableVector &dpdRCG );

    errorOut computeReferenceSecondOrderStressPressure( const variableVector &referenceStressMeasure,
                                                        const variableVector &rightCauchyGreen, variableType &pressure,
                                                        variableVector &dpdStress, variableVector &dpdRCG,
                                                        variableMatrix &d2pdStressdRCG );

    errorOut computeDeviatoricReferenceSecondOrderStress( const variableVector &secondOrderReferenceStress,
                                                          const variableVector &rightCauchyGreenDeformation,
                                                          variableVector &deviatoricSecondOrderReferenceStress );

    errorOut computeDeviatoricReferenceSecondOrderStress( const variableVector &secondOrderReferenceStress,
                                                          const variableVector &rightCauchyGreenDeformation,
                                                          const variableType &pressure,
                                                          variableVector &deviatoricSecondOrderReferenceStress );

    errorOut computeDeviatoricReferenceSecondOrderStress( const variableVector &secondOrderReferenceStress,
                                                          const variableVector &rightCauchyGreenDeformation,
                                                          variableVector &deviatoricSecondOrderReferenceStress,
                                                          variableMatrix &dDeviatoricReferenceStressdReferenceStress,
                                                          variableMatrix &dDeviatoricReferenceStressdRCG );

    errorOut computeDeviatoricReferenceSecondOrderStress( const variableVector &secondOrderReferenceStress,
                                                          const variableVector &rightCauchyGreenDeformation,
                                                          const variableType &pressure,
                                                          const variableVector &dPressuredStress,
                                                          const variableVector &dPressuredRCG,
                                                          variableVector &deviatoricSecondOrderReferenceStress,
                                                          variableMatrix &dDeviatoricReferenceStressdReferenceStress,
                                                          variableMatrix &dDeviatoricReferenceStressdRCG );

    errorOut computeDeviatoricReferenceSecondOrderStress( const variableVector &secondOrderReferenceStress,
                                                          const variableVector &rightCauchyGreenDeformation,
                                                          variableVector &deviatoricSecondOrderReferenceStress,
                                                          variableMatrix &dDeviatoricReferenceStressdReferenceStress,
                                                          variableMatrix &dDeviatoricReferenceStressdRCG,
                                                          variableMatrix &d2DevSdSdRCG);

    errorOut computeDeviatoricReferenceSecondOrderStress( const variableVector &secondOrderReferenceStress,
                                                          const variableVector &rightCauchyGreenDeformation,
                                                          const variableType &pressure,
                                                          const variableVector &dPressuredStress,
                                                          const variableVector &dPressuredRCG,
                                                          const variableMatrix &d2PressuredStressdRCG,
                                                          variableVector &deviatoricSecondOrderReferenceStress,
                                                          variableMatrix &dDeviatoricReferenceStressdReferenceStress,
                                                          variableMatrix &dDeviatoricReferenceStressdRCG,
                                                          variableMatrix &d2DevSdSdRCG);

    errorOut computeDeviatoricHigherOrderStress( const variableVector &higherOrderStress,
                                                 variableVector &deviatoricHigherOrderStress );
    
    errorOut computeDeviatoricHigherOrderStress( const variableVector &higherOrderStress,
                                                 variableVector &deviatoricHigherOrderStress,
                                                 variableMatrix &dDeviatoricHigherOrderStressdHigherOrderStress);

    errorOut computeReferenceHigherOrderStressPressure( const variableVector &referenceHigherOrderStress, 
                                                        const variableVector &rightCauchyGreenDeformation,
                                                        variableVector &referenceHigherOrderPressure );

    errorOut computeReferenceHigherOrderStressPressure( const variableVector &referenceHigherOrderStress, 
                                                        const variableVector &rightCauchyGreenDeformation,
                                                        variableVector &referenceHigherOrderPressure,
                                                        variableMatrix &dpdM, variableMatrix &dpdC );

    errorOut computeReferenceHigherOrderStressPressure( const variableVector &referenceHigherOrderStress,
                                                        const variableVector &rightCauchyGreenDeformation,
                                                        variableVector &referenceHigherOrderPressure,
                                                        variableMatrix &dpdM, variableMatrix &dpdC,
                                                        variableMatrix &d2pdMdC );

    errorOut computeDeviatoricReferenceHigherOrderStress( const variableVector &referenceHigherOrderStress,
                                                          const variableVector &rightCauchyGreenDeformation,
                                                          variableVector &deviatoricReferenceHigherOrderStress );

    errorOut computeDeviatoricReferenceHigherOrderStress( const variableVector &referenceHigherOrderStress,
                                                          const variableVector &rightCauchyGreenDeformation,
                                                          const variableVector &pressure,
                                                          variableVector &deviatoricReferenceHigherOrderStress );
    
    errorOut computeDeviatoricReferenceHigherOrderStress( const variableVector &referenceHigherOrderStress,
                                                          const variableVector &rightCauchyGreenDeformation,
                                                          variableVector &deviatoricReferenceHigherOrderStress,
                                                          variableMatrix &dDeviatoricReferenceHigherOrderStressdReferenceHigherOrderStress,
                                                          variableMatrix &dDeviatoricReferenceHigherOrderStressdRCG );

    errorOut computeDeviatoricReferenceHigherOrderStress( const variableVector &referenceHigherOrderStress,
                                                          const variableVector &rightCauchyGreenDeformation,
                                                          const variableVector &pressure,
                                                          const variableMatrix &dPressuredStress,
                                                          const variableMatrix &dPressuredRCG,
                                                          variableVector &deviatoricReferenceHigherOrderStress,
                                                          variableMatrix &dDeviatoricReferenceHigherOrderStressdReferenceHigherOrderStress,
                                                          variableMatrix &dDeviatoricReferenceHigherOrderStressdRCG );

    errorOut computeDeviatoricReferenceHigherOrderStress( const variableVector &referenceHigherOrderStress,
                                                          const variableVector &rightCauchyGreenDeformation,
                                                          variableVector &deviatoricReferenceHigherOrderStress,
                                                          variableMatrix &dDeviatoricReferenceHigherOrderStressdReferenceHigherOrderStress,
                                                          variableMatrix &dDeviatoricReferenceHigherOrderStressdRCG,
                                                          variableMatrix &d2DevMdMdRCG );

    errorOut computeDeviatoricReferenceHigherOrderStress( const variableVector &referenceHigherOrderStress,
                                                          const variableVector &rightCauchyGreenDeformation,
                                                          const variableVector &pressure,
                                                          const variableMatrix &dPressuredStress,
                                                          const variableMatrix &dPressuredRCG,
                                                          const variableMatrix &d2PressuredStressdRCG,
                                                          variableVector &deviatoricReferenceHigherOrderStress,
                                                          variableMatrix &dDeviatoricReferenceHigherOrderStressdReferenceHigherOrderStress,
                                                          variableMatrix &dDeviatoricReferenceHigherOrderStressdRCG,
                                                          variableMatrix &d2DevMdMdRCG );

    errorOut computeSecondOrderReferenceStressDecomposition( const variableVector &secondOrderReferenceStress,
                                                             const variableVector &rightCauchyGreenDeformation,
                                                             variableVector &deviatoricSecondOrderReferenceStress,
                                                             variableType &pressure );

    errorOut computeSecondOrderReferenceStressDecomposition( const variableVector &secondOrderReferenceStress,
                                                             const variableVector &rightCauchyGreenDeformation,
                                                             variableVector &deviatoricSecondOrderReferenceStress,
                                                             variableType &pressure, variableMatrix &dDevStressdStress,
                                                             variableMatrix &dDevStressdRCG, variableVector &dPressuredStress,
                                                             variableVector &dPressuredRCG );

    errorOut computeSecondOrderReferenceStressDecomposition( const variableVector &secondOrderReferenceStress,
                                                             const variableVector &rightCauchyGreenDeformation,
                                                             variableVector &deviatoricSecondOrderReferenceStress,
                                                             variableType &pressure, variableMatrix &dDevStressdStress,
                                                             variableMatrix &dDevStressdRCG, variableVector &dPressuredStress,
                                                             variableVector &dPressuredRCG, variableMatrix &d2DevStressdStressdRCG,
                                                             variableMatrix &d2PressuredStressdRCG );

    errorOut computeHigherOrderReferenceStressDecomposition( const variableVector &higherOrderReferenceStress,
                                                             const variableVector &rightCauchyGreenDeformation,
                                                             variableVector &deviatoricHigherOrderReferenceStress,
                                                             variableVector &pressure );

    errorOut computeHigherOrderReferenceStressDecomposition( const variableVector &higherOrderReferenceStress,
                                                             const variableVector &rightCauchyGreenDeformation,
                                                             variableVector &deviatoricHigherOrderReferenceStress,
                                                             variableVector &pressure, variableMatrix &dDevStressdStress,
                                                             variableMatrix &dDevStressdRCG, variableMatrix &dPressuredStress,
                                                             variableMatrix &dPressuredRCG );

    errorOut computeHigherOrderReferenceStressDecomposition( const variableVector &higherOrderReferenceStress,
                                                             const variableVector &rightCauchyGreenDeformation,
                                                             variableVector &deviatoricHigherOrderReferenceStress,
                                                             variableVector &pressure, variableMatrix &dDevStressdStress,
                                                             variableMatrix &dDevStressdRCG, variableMatrix &dPressuredStress,
                                                             variableMatrix &dPressuredRCG, variableMatrix &d2DevStressdStressdRCG,
                                                             variableMatrix &d2PressuredStressdRCG );

    errorOut computeHigherOrderStressNorm( const variableVector &higherOrderStress, variableVector &higherOrderStressNorm );

    errorOut computeHigherOrderStressNorm( const variableVector &higherOrderStress, variableVector &higherOrderStressNorm,
                                           variableMatrix &dHigherOrderStressNormdHigherOrderStress,
                                           double tol = 1e-9 );

    errorOut computeHigherOrderStressNorm( const variableVector &higherOrderStress, variableVector &higherOrderStressNorm,
                                           variableMatrix &dHigherOrderStressNormdHigherOrderStress,
                                           variableMatrix &d2HigherOrderStressNormdHigherOrderStress2,
                                           double tol = 1e-9 );

    errorOut assembleDeformationGradient( const variableMatrix &displacementGradient, variableVector &deformationGradient );

    errorOut assembleDeformationGradient( const variableMatrix &displacementGradient, variableVector &deformationGradient,
                                          variableMatrix &dFdGradU );

    errorOut assembleMicroDeformation( const variableVector &microDisplacement, variableVector &microDeformation );

    errorOut assembleMicroDeformation( const variableVector &microDisplacement, variableVector &microDeformation,
                                       variableMatrix &dChidPhi );

    errorOut assembleGradientMicroDeformation( const variableMatrix &gradientMicroDisplacement,
                                               variableVector &gradientMicroDeformation );

    errorOut assembleGradientMicroDeformation( const variableMatrix &gradientMicroDisplacement,
                                               variableVector &gradientMicroDeformation, variableMatrix &dGradChidGradPhi );
}

#endif
