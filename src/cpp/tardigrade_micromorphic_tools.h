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

    typedef double variableType; //!< Type definition for variable values
    typedef std::vector< variableType > variableVector; //!< Type definition for vectors of variable values
    typedef std::vector< variableVector > variableMatrix; //!< Type definition for matrices of variable values

    typedef double parameterType; //!< Type definition for parameter values
    typedef std::vector< parameterType > parameterVector; //!< Type definition for vectors of parameters
    typedef std::vector< parameterVector > parameterMatrix; //!< Type definition for matrices of parameters

    typedef double constantType; //!< Type definition for constants
    typedef std::vector< constantType > constantVector; //!< Type definition for vectors of constants
    typedef std::vector< constantVector > constantMatrix; //!< Type definition for matrices of constants

    void computePsi( const variableVector &deformationGradient, const variableVector &microDeformation,
                         variableVector &Psi );

    void computePsi( const variableVector &deformationGradient, const variableVector &microDeformation,
                         variableVector &Psi, variableVector &dPsidF, variableVector &dPsidChi );

    void computePsi( const variableVector &deformationGradient, const variableVector &microDeformation,
                         variableVector &Psi, variableMatrix &dPsidF, variableMatrix &dPsidChi );

    void computeGamma( const variableVector &deformationGradient, const variableVector &gradChi,
                           variableVector &Gamma );

    void computeGamma( const variableVector &deformationGradient, const variableVector &gradChi,
                           variableVector &Gamma, variableVector &dGammadF, variableVector &dGammadGradChi );

    void computeGamma( const variableVector &deformationGradient, const variableVector &gradChi,
                           variableVector &Gamma, variableMatrix &dGammadF, variableMatrix &dGammadGradChi );

    void computeMicroStrain( const variableVector &Psi, variableVector &microStrain );

    void computeMicroStrain( const variableVector &Psi, variableVector &microStrain,
                                 variableVector &dMicroStraindPsi );

    void computeMicroStrain( const variableVector &Psi, variableVector &microStrain,
                                 variableMatrix &dMicroStraindPsi );

    void pushForwardPK2Stress( const variableVector &PK2Stress,
                                   const variableVector &deformationGradient,
                                   variableVector &cauchyStress );

    void pushForwardPK2Stress( const variableVector &PK2Stress,
                                   const variableVector &deformationGradient,
                                   variableVector &cauchyStress, 
                                   variableVector &dcauchyStressdPK2Stress,
                                   variableVector &dcauchyStressdDeformationGradient );

    void pushForwardPK2Stress( const variableVector &PK2Stress,
                               const variableVector &deformationGradient,
                               variableVector &cauchyStress, 
                               variableMatrix &dcauchyStressdPK2Stress,
                               variableMatrix &dcauchyStressdDeformationGradient );

    void dCauchyStressdPK2Stress( const variableVector &deformationGradient,
                                  variableVector &dCauchyStressdPK2Stress );

    void pullBackCauchyStress( const variableVector &cauchyStress,
                                   const variableVector &deformationGradient,
                                   variableVector &PK2Stress );

    void pullBackCauchyStress( const variableVector &cauchyStress,
                                   const variableVector &deformationGradient,
                                   variableVector &PK2Stress, 
                                   variableVector &dcauchyStressdPK2Stress,
                                   variableVector &dcauchyStressdDeformationGradient );

    void pullBackCauchyStress( const variableVector &cauchyStress,
                                   const variableVector &deformationGradient,
                                   variableVector &PK2Stress, 
                                   variableMatrix &dcauchyStressdPK2Stress,
                                   variableMatrix &dcauchyStressdDeformationGradient );

    void pushForwardReferenceMicroStress( const variableVector &referenceMicroStress,
                                              const variableVector &deformationGradient,
                                              variableType &detF, variableVector &microStress );

    void pushForwardReferenceMicroStress( const variableVector &referenceMicroStress,
                                              const variableVector &deformationGradient,
                                              variableVector &microStress );

    void pushForwardReferenceMicroStress( const variableVector &referenceMicroStress,
                                              const variableVector &deformationGradient,
                                              variableVector &microStress, 
                                              variableVector &dMicroStressdReferenceMicroStress,
                                              variableVector &dMicroStressdDeformationGradient );

    void pushForwardReferenceMicroStress( const variableVector &referenceMicroStress,
                                              const variableVector &deformationGradient,
                                              variableVector &microStress, 
                                              variableMatrix &dMicroStressdReferenceMicroStress,
                                              variableMatrix &dMicroStressdDeformationGradient );

    void dSymmetricMicroStressdReferenceSymmetricMicroStress( const variableVector &deformationGradient,
                                                              variableVector &dSymmetricMicroStressdReferenceSymmetricMicroStress );

    void pullBackMicroStress( const variableVector &microStress,
                                  const variableVector &deformationGradient,
                                  variableType &detF, variableVector &inverseDeformationGradient,
                                  variableVector &referenceMicroStress );

    void pullBackMicroStress( const variableVector &microStress,
                                  const variableVector &deformationGradient,
                                  variableVector &referenceMicroStress );

    void pullBackMicroStress( const variableVector &microStress,
                                  const variableVector &deformationGradient,
                                  variableVector &referenceMicroStress, 
                                  variableVector &dReferenceMicroStressdMicroStress,
                                  variableVector &dReferenceMicroStressdDeformationGradient );

    void pullBackMicroStress( const variableVector &microStress,
                                  const variableVector &deformationGradient,
                                  variableVector &referenceMicroStress, 
                                  variableMatrix &dReferenceMicroStressdMicroStress,
                                  variableMatrix &dReferenceMicroStressdDeformationGradient );

    void pushForwardHigherOrderStress( const variableVector &referenceHigherOrderStress,
                                           const variableVector &deformationGradient,
                                           const variableVector &microDeformation,
                                           variableVector &higherOrderStress );

    void pushForwardHigherOrderStress( const variableVector &referenceHigherOrderStress,
                                           const variableVector &deformationGradient,
                                           const variableVector &microDeformation,
                                           variableType &detF,
                                           variableVector &higherOrderStress );

    void pushForwardHigherOrderStress( const variableVector &referenceHigherOrderStress,
                                           const variableVector &deformationGradient,
                                           const variableVector &microDeformation,
                                           variableVector &higherOrderStress,
                                           variableVector &dHigherOrderStressdReferenceHigherOrderStress,
                                           variableVector &dHigherOrderStressdDeformationGradient,
                                           variableVector &dHigherOrderStressdMicroDeformation );

    void pushForwardHigherOrderStress( const variableVector &referenceHigherOrderStress,
                                           const variableVector &deformationGradient,
                                           const variableVector &microDeformation,
                                           variableVector &higherOrderStress,
                                           variableMatrix &dHigherOrderStressdReferenceHigherOrderStress,
                                           variableMatrix &dHigherOrderStressdDeformationGradient,
                                           variableMatrix &dHigherOrderStressdMicroDeformation );

    void dHigherOrderStressdReferenceHigherOrderStress( const variableVector &deformationGradient,
                                                        const variableVector &microDeformation,
                                                        variableVector &dHigherOrderStressdReferenceHigherOrderStress );

    void pullBackHigherOrderStress( const variableVector &higherOrderStress,
                                        const variableVector &deformationGradient,
                                        const variableVector &microDeformation,
                                        variableVector &referenceHigherOrderStress );

    void pullBackHigherOrderStress( const variableVector &higherOrderStress,
                                        const variableVector &deformationGradient,
                                        const variableVector &microDeformation,
                                        variableType &detF, variableVector &inverseDeformationGradient,
                                        variableVector &inverseMicroDeformation,
                                        variableVector &referenceHigherOrderStress );

    void pullBackHigherOrderStress( const variableVector &higherOrderStress,
                                        const variableVector &deformationGradient,
                                        const variableVector &microDeformation,
                                        variableVector &referenceHigherOrderStress,
                                        variableVector &dHigherOrderStressdReferenceHigherOrderStress,
                                        variableVector &dHigherOrderStressdDeformationGradient,
                                        variableVector &dHigherOrderStressdMicroDeformation );

    void pullBackHigherOrderStress( const variableVector &higherOrderStress,
                                        const variableVector &deformationGradient,
                                        const variableVector &microDeformation,
                                        variableVector &referenceHigherOrderStress,
                                        variableMatrix &dHigherOrderStressdReferenceHigherOrderStress,
                                        variableMatrix &dHigherOrderStressdDeformationGradient,
                                        variableMatrix &dHigherOrderStressdMicroDeformation );

    void computeDeviatoricSecondOrderStress( const variableVector &secondOrderStress,
                                                 variableVector &deviatoricSecondOrderStress );

    void computeDeviatoricSecondOrderStress( const variableVector &secondOrderStress,
                                                 variableVector &deviatoricSecondOrderStress,
                                                 variableVector &dDeviatoricStressdStress );

    void computeDeviatoricSecondOrderStress( const variableVector &secondOrderStress,
                                                 variableVector &deviatoricSecondOrderStress,
                                                 variableMatrix &dDeviatoricStressdStress );

    void computeReferenceSecondOrderStressPressure( const variableVector &referenceStressMeasure,
                                                        const variableVector &rightCauchyGreen,  variableType &pressure );

    void computeReferenceSecondOrderStressPressure( const variableVector &referenceStressMeasure,
                                                        const variableVector &rightCauchyGreen, variableType &pressure,
                                                        variableVector &dpdStress, variableVector &dpdRCG );

    void computeReferenceSecondOrderStressPressure( const variableVector &referenceStressMeasure,
                                                        const variableVector &rightCauchyGreen, variableType &pressure,
                                                        variableVector &dpdStress, variableVector &dpdRCG,
                                                        variableVector &d2pdStressdRCG );

    void computeReferenceSecondOrderStressPressure( const variableVector &referenceStressMeasure,
                                                        const variableVector &rightCauchyGreen, variableType &pressure,
                                                        variableVector &dpdStress, variableVector &dpdRCG,
                                                        variableMatrix &d2pdStressdRCG );

    void computeDeviatoricReferenceSecondOrderStress( const variableVector &secondOrderReferenceStress,
                                                          const variableVector &rightCauchyGreenDeformation,
                                                          variableVector &deviatoricSecondOrderReferenceStress );

    void computeDeviatoricReferenceSecondOrderStress( const variableVector &secondOrderReferenceStress,
                                                          const variableVector &rightCauchyGreenDeformation,
                                                          const variableType &pressure,
                                                          variableVector &deviatoricSecondOrderReferenceStress );

    void computeDeviatoricReferenceSecondOrderStress( const variableVector &secondOrderReferenceStress,
                                                          const variableVector &rightCauchyGreenDeformation,
                                                          variableVector &deviatoricSecondOrderReferenceStress,
                                                          variableMatrix &dDeviatoricReferenceStressdReferenceStress,
                                                          variableMatrix &dDeviatoricReferenceStressdRCG );

    void computeDeviatoricReferenceSecondOrderStress( const variableVector &secondOrderReferenceStress,
                                                          const variableVector &rightCauchyGreenDeformation,
                                                          variableVector &deviatoricSecondOrderReferenceStress,
                                                          variableVector &dDeviatoricReferenceStressdReferenceStress,
                                                          variableVector &dDeviatoricReferenceStressdRCG );

    void computeDeviatoricReferenceSecondOrderStress( const variableVector &secondOrderReferenceStress,
                                                          const variableVector &rightCauchyGreenDeformation,
                                                          const variableType &pressure,
                                                          const variableVector &dPressuredStress,
                                                          const variableVector &dPressuredRCG,
                                                          variableVector &deviatoricSecondOrderReferenceStress,
                                                          variableVector &dDeviatoricReferenceStressdReferenceStress,
                                                          variableVector &dDeviatoricReferenceStressdRCG );

    void computeDeviatoricReferenceSecondOrderStress( const variableVector &secondOrderReferenceStress,
                                                          const variableVector &rightCauchyGreenDeformation,
                                                          variableVector &deviatoricSecondOrderReferenceStress,
                                                          variableMatrix &dDeviatoricReferenceStressdReferenceStress,
                                                          variableMatrix &dDeviatoricReferenceStressdRCG,
                                                          variableMatrix &d2DevSdSdRCG);

    void computeDeviatoricReferenceSecondOrderStress( const variableVector &secondOrderReferenceStress,
                                                          const variableVector &rightCauchyGreenDeformation,
                                                          variableVector &deviatoricSecondOrderReferenceStress,
                                                          variableVector &dDeviatoricReferenceStressdReferenceStress,
                                                          variableVector &dDeviatoricReferenceStressdRCG,
                                                          variableVector &d2DevSdSdRCG);

    void computeDeviatoricReferenceSecondOrderStress( const variableVector &secondOrderReferenceStress,
                                                          const variableVector &rightCauchyGreenDeformation,
                                                          const variableType &pressure,
                                                          const variableVector &dPressuredStress,
                                                          const variableVector &dPressuredRCG,
                                                          const variableVector &d2PressuredStressdRCG,
                                                          variableVector &deviatoricSecondOrderReferenceStress,
                                                          variableVector &dDeviatoricReferenceStressdReferenceStress,
                                                          variableVector &dDeviatoricReferenceStressdRCG,
                                                          variableVector &d2DevSdSdRCG);

    void computeDeviatoricHigherOrderStress( const variableVector &higherOrderStress,
                                                 variableVector &deviatoricHigherOrderStress );
    
    void computeDeviatoricHigherOrderStress( const variableVector &higherOrderStress,
                                                 variableVector &deviatoricHigherOrderStress,
                                                 variableVector &dDeviatoricHigherOrderStressdHigherOrderStress);

    void computeDeviatoricHigherOrderStress( const variableVector &higherOrderStress,
                                                 variableVector &deviatoricHigherOrderStress,
                                                 variableMatrix &dDeviatoricHigherOrderStressdHigherOrderStress);

    void computeReferenceHigherOrderStressPressure( const variableVector &referenceHigherOrderStress, 
                                                        const variableVector &rightCauchyGreenDeformation,
                                                        variableVector &referenceHigherOrderPressure );

    void computeReferenceHigherOrderStressPressure( const variableVector &referenceHigherOrderStress, 
                                                        const variableVector &rightCauchyGreenDeformation,
                                                        variableVector &referenceHigherOrderPressure,
                                                        variableVector &dpdM, variableVector &dpdC );

    void computeReferenceHigherOrderStressPressure( const variableVector &referenceHigherOrderStress, 
                                                        const variableVector &rightCauchyGreenDeformation,
                                                        variableVector &referenceHigherOrderPressure,
                                                        variableMatrix &dpdM, variableMatrix &dpdC );

    void computeReferenceHigherOrderStressPressure( const variableVector &referenceHigherOrderStress,
                                                        const variableVector &rightCauchyGreenDeformation,
                                                        variableVector &referenceHigherOrderPressure,
                                                        variableVector &dpdM, variableVector &dpdC,
                                                        variableVector &d2pdMdC );

    void computeReferenceHigherOrderStressPressure( const variableVector &referenceHigherOrderStress,
                                                        const variableVector &rightCauchyGreenDeformation,
                                                        variableVector &referenceHigherOrderPressure,
                                                        variableMatrix &dpdM, variableMatrix &dpdC,
                                                        variableMatrix &d2pdMdC );

    void computeDeviatoricReferenceHigherOrderStress( const variableVector &referenceHigherOrderStress,
                                                          const variableVector &rightCauchyGreenDeformation,
                                                          variableVector &deviatoricReferenceHigherOrderStress );

    void computeDeviatoricReferenceHigherOrderStress( const variableVector &referenceHigherOrderStress,
                                                          const variableVector &rightCauchyGreenDeformation,
                                                          const variableVector &pressure,
                                                          variableVector &deviatoricReferenceHigherOrderStress );
    
    void computeDeviatoricReferenceHigherOrderStress( const variableVector &referenceHigherOrderStress,
                                                          const variableVector &rightCauchyGreenDeformation,
                                                          variableVector &deviatoricReferenceHigherOrderStress,
                                                          variableVector &dDeviatoricReferenceHigherOrderStressdReferenceHigherOrderStress,
                                                          variableVector &dDeviatoricReferenceHigherOrderStressdRCG );

    void computeDeviatoricReferenceHigherOrderStress( const variableVector &referenceHigherOrderStress,
                                                          const variableVector &rightCauchyGreenDeformation,
                                                          variableVector &deviatoricReferenceHigherOrderStress,
                                                          variableMatrix &dDeviatoricReferenceHigherOrderStressdReferenceHigherOrderStress,
                                                          variableMatrix &dDeviatoricReferenceHigherOrderStressdRCG );

    void computeDeviatoricReferenceHigherOrderStress( const variableVector &referenceHigherOrderStress,
                                                          const variableVector &rightCauchyGreenDeformation,
                                                          const variableVector &pressure,
                                                          const variableVector &dPressuredStress,
                                                          const variableVector &dPressuredRCG,
                                                          variableVector &deviatoricReferenceHigherOrderStress,
                                                          variableVector &dDeviatoricReferenceHigherOrderStressdReferenceHigherOrderStress,
                                                          variableVector &dDeviatoricReferenceHigherOrderStressdRCG );

    void computeDeviatoricReferenceHigherOrderStress( const variableVector &referenceHigherOrderStress,
                                                          const variableVector &rightCauchyGreenDeformation,
                                                          variableVector &deviatoricReferenceHigherOrderStress,
                                                          variableVector &dDeviatoricReferenceHigherOrderStressdReferenceHigherOrderStress,
                                                          variableVector &dDeviatoricReferenceHigherOrderStressdRCG,
                                                          variableVector &d2DevMdMdRCG );

    void computeDeviatoricReferenceHigherOrderStress( const variableVector &referenceHigherOrderStress,
                                                          const variableVector &rightCauchyGreenDeformation,
                                                          variableVector &deviatoricReferenceHigherOrderStress,
                                                          variableMatrix &dDeviatoricReferenceHigherOrderStressdReferenceHigherOrderStress,
                                                          variableMatrix &dDeviatoricReferenceHigherOrderStressdRCG,
                                                          variableMatrix &d2DevMdMdRCG );

    void computeDeviatoricReferenceHigherOrderStress( const variableVector &referenceHigherOrderStress,
                                                          const variableVector &rightCauchyGreenDeformation,
                                                          const variableVector &pressure,
                                                          const variableVector &dPressuredStress,
                                                          const variableVector &dPressuredRCG,
                                                          const variableVector &d2PressuredStressdRCG,
                                                          variableVector &deviatoricReferenceHigherOrderStress,
                                                          variableVector &dDeviatoricReferenceHigherOrderStressdReferenceHigherOrderStress,
                                                          variableVector &dDeviatoricReferenceHigherOrderStressdRCG,
                                                          variableVector &d2DevMdMdRCG );

    void computeSecondOrderReferenceStressDecomposition( const variableVector &secondOrderReferenceStress,
                                                             const variableVector &rightCauchyGreenDeformation,
                                                             variableVector &deviatoricSecondOrderReferenceStress,
                                                             variableType &pressure );

    void computeSecondOrderReferenceStressDecomposition( const variableVector &secondOrderReferenceStress,
                                                             const variableVector &rightCauchyGreenDeformation,
                                                             variableVector &deviatoricSecondOrderReferenceStress,
                                                             variableType &pressure, variableVector &dDevStressdStress,
                                                             variableVector &dDevStressdRCG, variableVector &dPressuredStress,
                                                             variableVector &dPressuredRCG );

    void computeSecondOrderReferenceStressDecomposition( const variableVector &secondOrderReferenceStress,
                                                             const variableVector &rightCauchyGreenDeformation,
                                                             variableVector &deviatoricSecondOrderReferenceStress,
                                                             variableType &pressure, variableMatrix &dDevStressdStress,
                                                             variableMatrix &dDevStressdRCG, variableVector &dPressuredStress,
                                                             variableVector &dPressuredRCG );

    void computeSecondOrderReferenceStressDecomposition( const variableVector &secondOrderReferenceStress,
                                                             const variableVector &rightCauchyGreenDeformation,
                                                             variableVector &deviatoricSecondOrderReferenceStress,
                                                             variableType &pressure, variableMatrix &dDevStressdStress,
                                                             variableMatrix &dDevStressdRCG, variableVector &dPressuredStress,
                                                             variableVector &dPressuredRCG, variableMatrix &d2DevStressdStressdRCG,
                                                             variableMatrix &d2PressuredStressdRCG );

    void computeSecondOrderReferenceStressDecomposition( const variableVector &secondOrderReferenceStress,
                                                             const variableVector &rightCauchyGreenDeformation,
                                                             variableVector &deviatoricSecondOrderReferenceStress,
                                                             variableType &pressure, variableVector &dDevStressdStress,
                                                             variableVector &dDevStressdRCG, variableVector &dPressuredStress,
                                                             variableVector &dPressuredRCG, variableVector &d2DevStressdStressdRCG,
                                                             variableVector &d2PressuredStressdRCG );

    void computeHigherOrderReferenceStressDecomposition( const variableVector &higherOrderReferenceStress,
                                                             const variableVector &rightCauchyGreenDeformation,
                                                             variableVector &deviatoricHigherOrderReferenceStress,
                                                             variableVector &pressure );

    void computeHigherOrderReferenceStressDecomposition( const variableVector &higherOrderReferenceStress,
                                                             const variableVector &rightCauchyGreenDeformation,
                                                             variableVector &deviatoricHigherOrderReferenceStress,
                                                             variableVector &pressure, variableVector &dDevStressdStress,
                                                             variableVector &dDevStressdRCG, variableVector &dPressuredStress,
                                                             variableVector &dPressuredRCG );

    void computeHigherOrderReferenceStressDecomposition( const variableVector &higherOrderReferenceStress,
                                                             const variableVector &rightCauchyGreenDeformation,
                                                             variableVector &deviatoricHigherOrderReferenceStress,
                                                             variableVector &pressure, variableMatrix &dDevStressdStress,
                                                             variableMatrix &dDevStressdRCG, variableMatrix &dPressuredStress,
                                                             variableMatrix &dPressuredRCG );

    void computeHigherOrderReferenceStressDecomposition( const variableVector &higherOrderReferenceStress,
                                                             const variableVector &rightCauchyGreenDeformation,
                                                             variableVector &deviatoricHigherOrderReferenceStress,
                                                             variableVector &pressure, variableVector &dDevStressdStress,
                                                             variableVector &dDevStressdRCG, variableVector &dPressuredStress,
                                                             variableVector &dPressuredRCG, variableVector &d2DevStressdStressdRCG,
                                                             variableVector &d2PressuredStressdRCG );

    void computeHigherOrderReferenceStressDecomposition( const variableVector &higherOrderReferenceStress,
                                                             const variableVector &rightCauchyGreenDeformation,
                                                             variableVector &deviatoricHigherOrderReferenceStress,
                                                             variableVector &pressure, variableMatrix &dDevStressdStress,
                                                             variableMatrix &dDevStressdRCG, variableMatrix &dPressuredStress,
                                                             variableMatrix &dPressuredRCG, variableMatrix &d2DevStressdStressdRCG,
                                                             variableMatrix &d2PressuredStressdRCG );

    void computeHigherOrderStressNorm( const variableVector &higherOrderStress, variableVector &higherOrderStressNorm );

    void computeHigherOrderStressNorm( const variableVector &higherOrderStress, variableVector &higherOrderStressNorm,
                                           variableVector &dHigherOrderStressNormdHigherOrderStress,
                                           double tol = 1e-9 );

    void computeHigherOrderStressNorm( const variableVector &higherOrderStress, variableVector &higherOrderStressNorm,
                                           variableMatrix &dHigherOrderStressNormdHigherOrderStress,
                                           double tol = 1e-9 );

    void computeHigherOrderStressNorm( const variableVector &higherOrderStress, variableVector &higherOrderStressNorm,
                                           variableVector &dHigherOrderStressNormdHigherOrderStress,
                                           variableVector &d2HigherOrderStressNormdHigherOrderStress2,
                                           double tol = 1e-9 );

    void computeHigherOrderStressNorm( const variableVector &higherOrderStress, variableVector &higherOrderStressNorm,
                                           variableMatrix &dHigherOrderStressNormdHigherOrderStress,
                                           variableMatrix &d2HigherOrderStressNormdHigherOrderStress2,
                                           double tol = 1e-9 );

    void assembleDeformationGradient( const variableVector &displacementGradient, variableVector &deformationGradient );

    void assembleDeformationGradient( const variableVector &displacementGradient, variableVector &deformationGradient,
                                          variableVector &dFdGradU );

    void assembleMicroDeformation( const variableVector &microDisplacement, variableVector &microDeformation );

    void assembleMicroDeformation( const variableVector &microDisplacement, variableVector &microDeformation,
                                       variableVector &dChidPhi );

    void assembleGradientMicroDeformation( const variableVector &gradientMicroDisplacement,
                                               variableVector &gradientMicroDeformation );

    void assembleGradientMicroDeformation( const variableVector &gradientMicroDisplacement,
                                               variableVector &gradientMicroDeformation, variableVector &dGradChidGradPhi );
}

#endif
