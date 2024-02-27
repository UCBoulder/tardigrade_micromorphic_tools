//Tests for tardigrade_constitutive_tools

#include<tardigrade_micromorphic_tools.h>
#include<sstream>
#include<fstream>
#include<iostream>

#define BOOST_TEST_MODULE test_tardigrade_micromorphic_tools
#include <boost/test/included/unit_test.hpp>

typedef tardigradeMicromorphicTools::constantType constantType;
typedef tardigradeMicromorphicTools::constantVector constantVector;
typedef tardigradeMicromorphicTools::constantMatrix constantMatrix;

typedef tardigradeMicromorphicTools::parameterType parameterType;
typedef tardigradeMicromorphicTools::parameterVector parameterVector;
typedef tardigradeMicromorphicTools::parameterMatrix parameterMatrix;

typedef tardigradeMicromorphicTools::variableType variableType;
typedef tardigradeMicromorphicTools::variableVector variableVector;
typedef tardigradeMicromorphicTools::variableMatrix variableMatrix;

typedef tardigradeMicromorphicTools::errorNode errorNode;
typedef tardigradeMicromorphicTools::errorOut errorOut;

struct cout_redirect{
    cout_redirect( std::streambuf * new_buffer)
        : old( std::cout.rdbuf( new_buffer ) )
    { }

    ~cout_redirect( ) {
        std::cout.rdbuf( old );
    }

    private:
        std::streambuf * old;
};

struct cerr_redirect{
    cerr_redirect( std::streambuf * new_buffer)
        : old( std::cerr.rdbuf( new_buffer ) )
    { }

    ~cerr_redirect( ) {
        std::cerr.rdbuf( old );
    }

    private:
        std::streambuf * old;
};

BOOST_AUTO_TEST_CASE( test_computePsi ){
    /*!
     * Tests of the compute Psi function.
     *
     * :param std::ofstream &results: The output file.
     */

    variableVector F  = { 0, 1, 2, 3, 4, 5, 6, 7, 8 };
    variableVector Chi = { 9, 10, 11, 12, 13, 14, 15, 16, 17 };

    variableVector answer = { 126, 135, 144, 162, 174, 186, 198, 213, 228 };

    variableVector result;

    errorOut error = tardigradeMicromorphicTools::computePsi( F, Chi, result );

    BOOST_CHECK( !error );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( result, answer ) );

    //Test Jacobians
    variableVector resultJ;
    variableMatrix dPsidF, dPsidChi;

    error = tardigradeMicromorphicTools::computePsi( F, Chi, resultJ, dPsidF, dPsidChi );

    BOOST_CHECK( !error );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( resultJ, answer ) );

    //Test dPsidF
    
    constantType eps = 1e-6;
    for ( unsigned int i = 0; i < F.size(); i++ ){
        constantVector delta( F.size(), 0 );
        delta[i] = eps * fabs( F[i] ) + eps;

        error = tardigradeMicromorphicTools::computePsi( F + delta, Chi, resultJ );

        BOOST_CHECK( !error );

        constantVector gradCol = ( resultJ - result ) / delta[i];

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){

            BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( gradCol[ j ], dPsidF[ j ][ i ] ) );

        }

    }

    for ( unsigned int i = 0; i < Chi.size(); i++ ){
        constantVector delta( Chi.size(), 0 );
        delta[i] = eps * fabs( Chi[i] ) + eps;

        error = tardigradeMicromorphicTools::computePsi( F, Chi + delta, resultJ );

        BOOST_CHECK( !error );

        constantVector gradCol = ( resultJ - result ) / delta[i];

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){

            BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( gradCol[ j ], dPsidChi[ j ][ i ] ) );
        }
    }

    return;

}

BOOST_AUTO_TEST_CASE( test_computeMicroStrain ){
    /*!
     * Test the computation of the micro-strain
     *
     * :param std::ofstream &results: The output file
     */

    variableVector Psi = { 1, 2, 3, 4, 5, 6, 7, 8, 9 };

    variableVector answer = { 0, 2, 3, 4, 4, 6, 7, 8, 8 };

    variableVector result;

    errorOut error = tardigradeMicromorphicTools::computeMicroStrain( Psi, result );

    BOOST_CHECK( !error );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( result, answer ) );

    //Test Jacobians 
    
    variableVector resultJ;
    variableMatrix dMicroStraindPsi;
    error = tardigradeMicromorphicTools::computeMicroStrain( Psi, resultJ, dMicroStraindPsi );

    BOOST_CHECK( !error );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( resultJ, answer ) );

    constantType eps = 1e-6;
    for ( unsigned int i = 0; i < Psi.size(); i++ ){
        constantVector delta( Psi.size(), 0 );
        delta[i] = eps * fabs( Psi[i] ) + eps;

        error = tardigradeMicromorphicTools::computeMicroStrain( Psi + delta, resultJ );

        BOOST_CHECK( !error );

        constantVector gradCol = ( resultJ - result ) / delta[i];

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){

            BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( gradCol[ j ], dMicroStraindPsi[ j ][ i ] ) );

        }

    }

    return;
}

BOOST_AUTO_TEST_CASE( test_pushForwardPK2Stress ){
    /*!
     * Test the computation of the push-foward operation on the PK2 Stress.
     *
     * :param std::ofstream &results: The output file.
     */

    variableVector PK2Stress = { 1, 2, 3, 4, 5, 6, 7, 8, 9 };
    variableVector deformationGradient  = { -1.08831037, -0.66333427, -0.48239487,
                                            -0.904554  ,  1.28942848, -0.02156112,
                                            -0.08464824, -0.07730218,  0.86415668 };

    variableVector answer = { -10.75427056,   2.95576352,   4.7810659 ,
                                5.36943821,  -1.06947471,  -1.91553073,
                                7.58260457,  -1.61246489,  -2.82366599 };

    variableVector result;

    errorOut error = tardigradeMicromorphicTools::pushForwardPK2Stress( PK2Stress, deformationGradient, result );

    BOOST_CHECK( !error );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( result, answer ) );

    //Test Jacobian
    variableVector resultJ;
    variableMatrix dCauchydPK2, dCauchydF;

    error = tardigradeMicromorphicTools::pushForwardPK2Stress( PK2Stress, deformationGradient, resultJ,
                                                     dCauchydPK2, dCauchydF );

    BOOST_CHECK( !error );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( resultJ, answer ) );

    //Test dCauchydPK2
    constantType eps = 1e-6;
    for ( unsigned int i = 0; i < PK2Stress.size(); i++ ){
        constantVector delta( PK2Stress.size(), 0 );
        delta[i] = eps * fabs( PK2Stress[i] ) + eps;

        error = tardigradeMicromorphicTools::pushForwardPK2Stress( PK2Stress + delta, deformationGradient, resultJ );

        BOOST_CHECK( !error );

        constantVector gradCol = ( resultJ - result ) / delta[i];

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){

            BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( gradCol[ j ], dCauchydPK2[ j ][ i ] ) );
        }
    }

    for ( unsigned int i = 0; i < deformationGradient.size(); i++ ){
        constantVector delta( deformationGradient.size(), 0 );
        delta[i] = eps * fabs( deformationGradient[i] ) + eps;

        variableVector resultp, resultm;

        error = tardigradeMicromorphicTools::pushForwardPK2Stress( PK2Stress, deformationGradient + delta, resultp );

        BOOST_CHECK( !error );

        error = tardigradeMicromorphicTools::pushForwardPK2Stress( PK2Stress, deformationGradient - delta, resultm );

        BOOST_CHECK( !error );

        constantVector gradCol = ( resultp - resultm ) / ( 2 * delta[i] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){

            BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( gradCol[ j ], dCauchydF[ j ][ i ] ) );
        }
    }

    return;
}

BOOST_AUTO_TEST_CASE( test_pushForwardReferenceMicroStress ){
    /*!
     * Test the computation of the push-foward operation on the reference micro-stress.
     *
     * :param std::ofstream &results: The output file.
     */

    variableVector referenceMicroStress = { 1, 2, 3, 4, 5, 6, 7, 8, 9 };
    variableVector deformationGradient  = { -1.08831037, -0.66333427, -0.48239487,
                                            -0.904554  ,  1.28942848, -0.02156112,
                                            -0.08464824, -0.07730218,  0.86415668 };

    variableVector answer = { -10.75427056,   2.95576352,   4.7810659 ,
                                5.36943821,  -1.06947471,  -1.91553073,
                                7.58260457,  -1.61246489,  -2.82366599 };

    variableVector result;

    errorOut error = tardigradeMicromorphicTools::pushForwardReferenceMicroStress( referenceMicroStress, deformationGradient, result );

    BOOST_CHECK( !error );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( result, answer ) );

    //Test Jacobian
    variableVector resultJ;
    variableMatrix dsdS, dsdF;

    error = tardigradeMicromorphicTools::pushForwardReferenceMicroStress( referenceMicroStress, deformationGradient, resultJ,
                                                                dsdS, dsdF );

    BOOST_CHECK( !error );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( resultJ, answer ) );

    //Test dsdS
    constantType eps = 1e-6;
    for ( unsigned int i = 0; i < referenceMicroStress.size(); i++ ){
        constantVector delta( referenceMicroStress.size(), 0 );
        delta[i] = eps * fabs( referenceMicroStress[i] ) + eps;

        error = tardigradeMicromorphicTools::pushForwardReferenceMicroStress( referenceMicroStress + delta, deformationGradient, resultJ );

        BOOST_CHECK( !error );

        constantVector gradCol = ( resultJ - result ) / delta[i];

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){

            BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( gradCol[ j ], dsdS[ j ][ i ] ) );

        }
    }

    for ( unsigned int i = 0; i < deformationGradient.size(); i++ ){
        constantVector delta( deformationGradient.size(), 0 );
        delta[i] = eps * fabs( deformationGradient[i] ) + eps;

        variableVector resultp, resultm;

        error = tardigradeMicromorphicTools::pushForwardReferenceMicroStress( referenceMicroStress, deformationGradient + delta, resultp );

        BOOST_CHECK( !error );

        error = tardigradeMicromorphicTools::pushForwardReferenceMicroStress( referenceMicroStress, deformationGradient - delta, resultm );

        BOOST_CHECK( !error );

        constantVector gradCol = ( resultp - resultm ) / ( 2 * delta[i] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){

            BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( gradCol[ j ], dsdF[ j ][ i ] ) );
        }
    }

    return;
}

BOOST_AUTO_TEST_CASE( test_computeGamma ){
    /*!
     * Test the computation of the deformation gradient Gamma
     *
     * :param std::ofstream &results: The output file.
     */

    variableVector deformationGradient = { -1, -2, -3, -4, -5, -6, -7, -8, -9 };
    variableVector gradChi = { 1,  2,  3,  4,  5,  6,  7,  8,  9,
                             10, 11, 12, 13, 14, 15, 16, 17, 18,
                             19, 20, 21, 22, 23, 24, 25, 26, 27 };

    variableVector answer = { -174, -186, -198, -210, -222, -234, -246, -258, -270,
                              -204, -219, -234, -249, -264, -279, -294, -309, -324,
                              -234, -252, -270, -288, -306, -324, -342, -360, -378 };

    variableVector result;

    errorOut error = tardigradeMicromorphicTools::computeGamma( deformationGradient, gradChi, result );

    BOOST_CHECK( !error );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( result, answer ) );

    //Test Jacobian
    variableVector resultJ;
    variableMatrix dGammadF, dGammadGradChi;

    error = tardigradeMicromorphicTools::computeGamma( deformationGradient, gradChi, resultJ, 
                                             dGammadF, dGammadGradChi );

    BOOST_CHECK( !error );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( resultJ, answer ) );

    //Test dGammadF
    constantType eps = 1e-6;
    for ( unsigned int i = 0; i < deformationGradient.size(); i++ ){
        constantVector delta( deformationGradient.size(), 0 );
        delta[i] = eps * fabs( deformationGradient[i] ) + eps;

        variableVector resultp, resultm;

        error = tardigradeMicromorphicTools::computeGamma( deformationGradient + delta, gradChi, resultp );

        BOOST_CHECK( !error );

        error = tardigradeMicromorphicTools::computeGamma( deformationGradient - delta, gradChi, resultm );

        BOOST_CHECK( !error );

        constantVector gradCol = ( resultp - resultm ) / ( 2 * delta[i] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){

            BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( gradCol[ j ], dGammadF[ j ][ i ] ) );

        }
    }

    //Test dGammadGradChi
    for ( unsigned int i = 0; i < gradChi.size(); i++ ){
        constantVector delta( gradChi.size(), 0 );
        delta[i] = eps * fabs( gradChi[i] ) + eps;

        variableVector resultp, resultm;

        error = tardigradeMicromorphicTools::computeGamma( deformationGradient, gradChi + delta, resultp );

        BOOST_CHECK( !error );

        error = tardigradeMicromorphicTools::computeGamma( deformationGradient, gradChi - delta, resultm );

        BOOST_CHECK( !error );

        constantVector gradCol = ( resultp - resultm ) / ( 2 * delta[i] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){

            BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( gradCol[ j ], dGammadGradChi[ j ][ i ] ) );
        }
    }

    return;

}

BOOST_AUTO_TEST_CASE( test_pushForwardHigherOrderStress ){
    /*!
     * Tests for the push-forward operation of the higher order stress.
     *
     * :param std::ofstream &results: The output file.
     */

    variableVector referenceHigherOrderStress = {  1,  2,  3,  4,  5,  6,  7,  8,  9,
                                                  10, 11, 12, 13, 14, 15, 16, 17, 18,
                                                  19, 20, 21, 22, 23, 24, 25, 26, 27 };

    variableVector deformationGradient = { 0.29524861, -0.25221581, -1.60534711,
                                          -0.08817703,  0.28447808,  0.06451703,
                                           0.39849454,  0.38681512,  0.87840084 };

    variableVector microDeformation = { -0.25781969, -0.39826899, -0.79493259,
                                         0.38104724, -0.00830511, -0.51985409,
                                        -0.36415661, -0.6871168 ,  0.54018665 };

    variableVector answer = { -370.21000924,  -44.9887908 , -120.76625915,   57.76488049,
                                 7.10106323,   18.73840733,  356.3462823 ,   44.06705478,
                               115.25802557,   49.68640146,    6.28202573,   15.89296064,
                                -7.62049271,   -0.98037632,   -2.41571071,  -46.58548828,
                                -6.04841312,  -14.69638769,  280.56462547,   36.38392302,
                                88.5657901 ,  -42.53702297,   -5.63795901,  -13.27041478,
                              -258.4236358 ,  -34.6543976 ,  -80.10152498 };

    variableVector result;

    errorOut error = tardigradeMicromorphicTools::pushForwardHigherOrderStress( referenceHigherOrderStress, deformationGradient,
                                                                      microDeformation, result );

    BOOST_CHECK( !error );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( result, answer ) );

    //Test Jacobian
    variableVector resultJ;
    variableMatrix dHigherOrderStressdReferenceHigherOrderStress;
    variableMatrix dHigherOrderStressdDeformationGradient;
    variableMatrix dHigherOrderStressdMicroDeformation;

    error = tardigradeMicromorphicTools::pushForwardHigherOrderStress( referenceHigherOrderStress,
                                                             deformationGradient, microDeformation,
                                                             resultJ,
                                                             dHigherOrderStressdReferenceHigherOrderStress,
                                                             dHigherOrderStressdDeformationGradient,
                                                             dHigherOrderStressdMicroDeformation );

    BOOST_CHECK( !error );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( resultJ, answer ) );

    //Test dHigherOrderStressdReferenceHigherOrderStress
    constantType eps = 1e-6;
    for ( unsigned int i = 0; i < referenceHigherOrderStress.size(); i++ ){
        constantVector delta( referenceHigherOrderStress.size(), 0 );
        delta[i] = eps * fabs( referenceHigherOrderStress[i] ) + eps;

        variableVector resultp, resultm;

        error = tardigradeMicromorphicTools::pushForwardHigherOrderStress( referenceHigherOrderStress + delta,
                                                                 deformationGradient, microDeformation,
                                                                 resultp );

        BOOST_CHECK( !error );

        error = tardigradeMicromorphicTools::pushForwardHigherOrderStress( referenceHigherOrderStress - delta,
                                                                 deformationGradient, microDeformation,
                                                                 resultm );

        BOOST_CHECK( !error );

        constantVector gradCol = ( resultp - resultm ) / ( 2 * delta[i] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){

            BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( gradCol[ j ],  dHigherOrderStressdReferenceHigherOrderStress[ j ][ i ] ) );

        }

    }

    //Test dHigherOrderStressdDeformationGradient
    for ( unsigned int i = 0; i < deformationGradient.size(); i++ ){
        constantVector delta( deformationGradient.size(), 0 );
        delta[i] = eps * fabs( deformationGradient[i] ) + eps;

        variableVector resultp, resultm;

        error = tardigradeMicromorphicTools::pushForwardHigherOrderStress( referenceHigherOrderStress,
                                                                 deformationGradient + delta, microDeformation,
                                                                 resultp );

        BOOST_CHECK( !error );

        error = tardigradeMicromorphicTools::pushForwardHigherOrderStress( referenceHigherOrderStress,
                                                                 deformationGradient - delta, microDeformation,
                                                                 resultm );

        BOOST_CHECK( !error );

        constantVector gradCol = ( resultp - resultm ) / ( 2 * delta[i] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){

            BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( gradCol[ j ], dHigherOrderStressdDeformationGradient[ j ][ i ] ) );

        }

    }

    //Test dHigherOrderStressdMicroDeformation
    for ( unsigned int i = 0; i < microDeformation.size(); i++ ){
        constantVector delta( microDeformation.size(), 0 );
        delta[i] = eps * fabs( microDeformation[i] ) + eps;

        variableVector resultp, resultm;

        error = tardigradeMicromorphicTools::pushForwardHigherOrderStress( referenceHigherOrderStress,
                                                                 deformationGradient, microDeformation + delta,
                                                                 resultp );

        BOOST_CHECK( !error );

        error = tardigradeMicromorphicTools::pushForwardHigherOrderStress( referenceHigherOrderStress,
                                                                 deformationGradient, microDeformation - delta,
                                                                 resultm );

        BOOST_CHECK( !error );

        constantVector gradCol = ( resultp - resultm ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){

            BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( gradCol[ j ], dHigherOrderStressdMicroDeformation[ j ][ i ] ) );

        }

    }

    return;

}

BOOST_AUTO_TEST_CASE( test_computeDeviatoricHigherOrderStress ){
    /*!
     * Test the computation of the deviatoric higher order stress.
     *
     * :param std::ofstream &results: The output file.
     */

    variableVector higherOrderStress = {  1,  2,  3,  4,  5,  6,  7,  8,  9,
                                         10, 11, 12, 13, 14, 15, 16, 17, 18,
                                         19, 20, 21, 22, 23, 24, 25, 26, 27 };

    variableVector answer = { -12., -12., -12.,   4.,   5.,   6.,   7.,   8.,   9.,
                               10.,  11.,  12.,   0.,   0.,   0.,  16.,  17.,  18.,
                               19.,  20.,  21.,  22.,  23.,  24.,  12.,  12.,  12. };

    variableVector result;

    errorOut error = tardigradeMicromorphicTools::computeDeviatoricHigherOrderStress( higherOrderStress, result );

    BOOST_CHECK( !error );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( answer, result ) );

    //Test Jacobians
    
    variableVector resultJ;
    variableMatrix dDeviatoricHigherOrderStressdHigherOrderStress;

    error = tardigradeMicromorphicTools::computeDeviatoricHigherOrderStress( higherOrderStress, resultJ,
                                                                   dDeviatoricHigherOrderStressdHigherOrderStress );

    BOOST_CHECK( !error );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( answer, resultJ ) );

    constantType eps = 1e-6;
    for ( unsigned int i = 0; i < higherOrderStress.size(); i++ ){
        constantVector delta( higherOrderStress.size(), 0 );
        delta[i] = eps * fabs( higherOrderStress[i] ) + eps;

        variableVector resultp, resultm;

        error = tardigradeMicromorphicTools::computeDeviatoricHigherOrderStress( higherOrderStress + delta, resultp );

        BOOST_CHECK( !error );

        error = tardigradeMicromorphicTools::computeDeviatoricHigherOrderStress( higherOrderStress - delta, resultm );

        BOOST_CHECK( !error );

        constantVector gradCol = ( resultp - resultm ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){

            BOOST_CHECK( tardigradeVectorTools::fuzzyEquals(  gradCol[ j ], dDeviatoricHigherOrderStressdHigherOrderStress[ j ][ i ] ) );

        }

    }

    return;

}

BOOST_AUTO_TEST_CASE( test_computeDeviatoricReferenceHigherOrderStress ){
    /*!
     * Test the computation of the deviatoric part of the reference higher order stress.
     *
     * :param std::ofstream &results: The output file
     */

    variableVector M = { 0.80732114,  0.79202055,  0.17990022,  0.97454675,  0.703207  ,
                        -0.58236697,  0.53324571, -0.93438873, -0.40650796,  0.14071918,
                         0.66933708, -0.67854069, -0.30317772, -0.93821882,  0.97270622,
                         0.00295302, -0.12441126,  0.30539971, -0.0580227 ,  0.89696105,
                         0.17567709, -0.9592962 ,  0.63535407,  0.95437804, -0.64531877,
                         0.69978907,  0.81327586 };

    variableVector C = { 0.3991656 ,  0.43459435,  0.15398811,
                        -0.20239202, -0.50763359,  0.04756988,
                        -0.0573016 , -0.95939895,  0.2693173 };

    variableVector answer = { 5.27633771,   1.68477127,  -3.21111472,  13.96704974,
                              3.29864744, -10.4408607 ,  -4.31689678,  -1.90327486,
                              3.27369895,  -2.4001685 ,   0.161758  ,   1.24944234,
                             -6.01118653,  -2.07847616,   5.30384773,   2.46397506,
                              0.3672135 ,  -1.56198261,  -8.15866529,  -0.72125951,
                              6.32230906, -18.52878201,  -2.87440491,  14.2858097 ,
                              4.98150383,   1.82382829,  -3.45626291 };

    variableVector pressure;
    variableMatrix dpdM, dpdC, d2pdMdC;
    tardigradeMicromorphicTools::computeReferenceHigherOrderStressPressure( M, C, pressure, dpdM, dpdC, d2pdMdC );

    variableVector result;

    errorOut error = tardigradeMicromorphicTools::computeDeviatoricReferenceHigherOrderStress( M, C, result );

    BOOST_CHECK( !error );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( result, answer ) );

    error = tardigradeMicromorphicTools::computeDeviatoricReferenceHigherOrderStress( M, C, pressure, result );

    BOOST_CHECK( !error );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( result, answer ) );

    variableVector resultJ;
    variableMatrix dDevMdM, dDevMdC;

    error = tardigradeMicromorphicTools::computeDeviatoricReferenceHigherOrderStress( M, C, resultJ, dDevMdM, dDevMdC );

    BOOST_CHECK( !error );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( resultJ, answer ) );

    variableVector resultJP;
    variableMatrix dDevMdMP, dDevMdCP;

    error = tardigradeMicromorphicTools::computeDeviatoricReferenceHigherOrderStress( M, C, pressure, dpdM, dpdC, resultJP, dDevMdMP, dDevMdCP );

    BOOST_CHECK( !error );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( resultJP, answer ) );

    variableVector resultJ2;
    variableMatrix dDevMdMJ2, dDevMdCJ2, d2DevMdMdC;

    error = tardigradeMicromorphicTools::computeDeviatoricReferenceHigherOrderStress( M, C, resultJ2, dDevMdMJ2, dDevMdCJ2, d2DevMdMdC );

    BOOST_CHECK( !error );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( resultJ, answer ) );

    variableVector resultJ2P;
    variableMatrix dDevMdMJ2P, dDevMdCJ2P, d2DevMdMdCP;

    error = tardigradeMicromorphicTools::computeDeviatoricReferenceHigherOrderStress( M, C, pressure, dpdM, dpdC, d2pdMdC,
                                                                            resultJ2P, dDevMdMJ2P, dDevMdCJ2P, d2DevMdMdCP );

    BOOST_CHECK( !error );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( resultJ2P, answer ) );

    //Test dDevMdM
    constantType eps = 1e-6;
    for ( unsigned int i = 0; i < M.size(); i++ ){
        constantVector delta( M.size(), 0 );
        delta[i] = eps * fabs( M[i] ) + eps;

        variableVector resultp, resultm;

        error = tardigradeMicromorphicTools::computeDeviatoricReferenceHigherOrderStress( M + delta, C, resultp );

        BOOST_CHECK( !error );

        error = tardigradeMicromorphicTools::computeDeviatoricReferenceHigherOrderStress( M - delta, C, resultm );

        BOOST_CHECK( !error );

        constantVector gradCol = ( resultp - resultm ) / ( 2 * delta[i] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){

            BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( gradCol[ j ], dDevMdM[ j ][ i ] ) );

        }

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){

            BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( gradCol[ j ], dDevMdMP[ j ][ i ] ) );

        }

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){

            BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( gradCol[ j ], dDevMdMJ2[ j ][ i ] ) );

        }

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){

            BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( gradCol[ j ], dDevMdMJ2P[ j ][ i ] ) );

        }

    }

    //Test dDevMdC
    for ( unsigned int i = 0; i < C.size(); i++ ){
        constantVector delta( C.size(), 0 );
        delta[i] = eps * fabs( C[i] ) + eps;

        variableVector resultp, resultm;

        error = tardigradeMicromorphicTools::computeDeviatoricReferenceHigherOrderStress( M, C + delta, resultp );

        BOOST_CHECK( !error );

        error = tardigradeMicromorphicTools::computeDeviatoricReferenceHigherOrderStress( M, C - delta, resultm );

        BOOST_CHECK( !error );

        constantVector gradCol = ( resultp - resultm ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){

            BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( gradCol[ j ], dDevMdC[ j ][ i ] ) );

        }

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){

            BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( gradCol[ j ], dDevMdCP[ j ][ i ] ) );

        }

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){

            BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( gradCol[ j ], dDevMdCJ2[ j ][ i ] ) );

        }

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){

            BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( gradCol[ j ], dDevMdCJ2P[ j ][ i ] ) );

        }

    }

    //Test d2DevMdMdC
    for ( unsigned int i = 0; i < C.size(); i++ ){
        constantVector delta( C.size(), 0 );
        delta[i] = eps * fabs( C[i] ) + eps;

        variableMatrix resultp, resultm;

        error = tardigradeMicromorphicTools::computeDeviatoricReferenceHigherOrderStress( M, C + delta, resultJ2, resultp, dDevMdCJ2 );

        BOOST_CHECK( !error );

        error = tardigradeMicromorphicTools::computeDeviatoricReferenceHigherOrderStress( M, C - delta, resultJ2, resultm, dDevMdCJ2 );

        BOOST_CHECK( !error );

        constantMatrix gradCol = ( resultp - resultm ) / ( 2 * delta[i] );

        for ( unsigned int j = 0; j < 3; j++ ){
            for ( unsigned int k = 0; k < 3; k++ ){
                for ( unsigned int l = 0; l < 3; l++ ){
                    for ( unsigned int m = 0; m < 3; m++ ){
                        for ( unsigned int n = 0; n < 3; n++ ){
                            for ( unsigned int o = 0; o < 3; o++ ){
                                BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( gradCol[ 9 * j + 3 * k + l ][ 9 * m + 3 * n + o ],
                                                                       d2DevMdMdC[ 9 * j + 3 * k + l ][ 81 * m + 27 * n + 9 * o + 3 * (int)( i / 3 ) + i % 3 ] ) );

                                BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( gradCol[ 9 * j + 3 * k + l ][ 9 * m + 3 * n + o ],
                                                                       d2DevMdMdCP[ 9 * j + 3 * k + l ][ 81 * m + 27 * n + 9 * o + 3 * (int)( i / 3 ) + i % 3 ] ) );
                            }
                        }
                    }
                }
            }
        }
    }

    return;

}

BOOST_AUTO_TEST_CASE( test_computeDeviatoricSecondOrderStress ){
    /*!
     * Test the computation the deviatoric part of a second order 
     * stress tensor in the current configuration.
     *
     * :param std::ofstream &results: The output file.
     */

    variableVector stress = { 1, 2, 3, 4, 5, 6, 7, 8, 9 };

    variableVector answer = { -4.,  2.,  3.,  4.,  0.,  6.,  7.,  8.,  4. };

    variableVector result;

    errorOut error = tardigradeMicromorphicTools::computeDeviatoricSecondOrderStress( stress, result );

    BOOST_CHECK( !error );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( answer, result ) );

    variableVector resultJ;
    variableMatrix jacobian;

    error = tardigradeMicromorphicTools::computeDeviatoricSecondOrderStress( stress, resultJ, jacobian );

    BOOST_CHECK( !error );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( answer, resultJ ) );

    constantType eps = 1e-6;
    for ( unsigned int i = 0; i < stress.size(); i++ ){
        constantVector delta( stress.size(), 0 );
        delta[i] = eps * fabs( stress[i] ) + eps;

        variableVector resultp, resultm;

        error = tardigradeMicromorphicTools::computeDeviatoricSecondOrderStress( stress + delta, resultp );

        BOOST_CHECK( !error );

        error = tardigradeMicromorphicTools::computeDeviatoricSecondOrderStress( stress - delta, resultm );

        BOOST_CHECK( !error );

        constantVector gradCol = ( resultp - resultm ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){

            BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( gradCol[ j ], jacobian[ j ][ i ] ) );

        }

    }

    return;

}

BOOST_AUTO_TEST_CASE( test_computeReferenceSecondOrderStressPressure ){
    /*!
     * Test the computation of the pressure of a second order stress in the 
     * reference configuration.
     *
     * :param std::ofstream &rsults: The output file.
     */

    variableVector S = { 0.77189588, -0.84417528,  0.95929231,
                        -0.50465708,  0.50576944,  0.05335127,
                         0.81510751,  0.76814059, -0.82146208 };

    variableVector C = { 0.03468919, -0.31275742, -0.57541261,
                        -0.27865312, -0.45844965,  0.52325004,
                        -0.0439162 , -0.80201065, -0.44921044 };

    variableType answer = -0.20245462701026676;

    variableType result;

    errorOut error = tardigradeMicromorphicTools::computeReferenceSecondOrderStressPressure( S, C, result );

    BOOST_CHECK( !error );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( answer, result ) );

    //Test the Jacobians

    variableType resultJ;
    variableVector dpdS, dpdC;

    error = tardigradeMicromorphicTools::computeReferenceSecondOrderStressPressure( S, C, resultJ, dpdS, dpdC );

    BOOST_CHECK( !error );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( resultJ, answer ) );

    variableType resultJ2;
    variableVector dpdSJ2, dpdCJ2;
    variableMatrix d2pdSdC;

    error = tardigradeMicromorphicTools::computeReferenceSecondOrderStressPressure( S, C, resultJ2, dpdSJ2, dpdCJ2, d2pdSdC );

    BOOST_CHECK( !error );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( resultJ2, answer ) );

    //Test dpdS
    constantType eps = 1e-6;
    for ( unsigned int i = 0; i < S.size(); i++ ){
        constantVector delta( S.size(), 0 );
        delta[i] = eps * fabs( S[i] ) + eps;

        variableType resultp, resultm;

        error = tardigradeMicromorphicTools::computeReferenceSecondOrderStressPressure( S + delta, C, resultp );

        BOOST_CHECK( !error );

        error = tardigradeMicromorphicTools::computeReferenceSecondOrderStressPressure( S - delta, C, resultm );

        BOOST_CHECK( !error );

        constantType gradCol = ( resultp - resultm ) / ( 2 * delta[i] );

        BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( gradCol, dpdS[ i ] ) );

        BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( gradCol, dpdSJ2[ i ] ) );

    }

    //Test dpdC
    for ( unsigned int i = 0; i < C.size(); i++ ){
        constantVector delta( C.size(), 0 );
        delta[i] = eps * fabs( C[i] ) + eps;

        variableType resultp, resultm;

        error = tardigradeMicromorphicTools::computeReferenceSecondOrderStressPressure( S, C + delta, resultp );

        BOOST_CHECK( !error );

        error = tardigradeMicromorphicTools::computeReferenceSecondOrderStressPressure( S, C - delta, resultm );

        BOOST_CHECK( !error );

        constantType gradCol = ( resultp - resultm ) / ( 2 * delta[ i ] );

        BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( gradCol, dpdC[ i ] ) );

        BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( gradCol, dpdCJ2[ i ] ) );

    }

    //Test d2pdSdC
    for ( unsigned int i = 0; i < C.size(); i++ ){
        constantVector delta( C.size(), 0 );
        delta[i] = eps * fabs( C[i] ) + eps;

        variableVector resultp, resultm;

        error = tardigradeMicromorphicTools::computeReferenceSecondOrderStressPressure( S, C + delta, resultJ, resultp, dpdCJ2 );

        BOOST_CHECK( !error );

        error = tardigradeMicromorphicTools::computeReferenceSecondOrderStressPressure( S, C - delta, resultJ, resultm, dpdCJ2 );

        BOOST_CHECK( !error );

        constantVector gradCol = ( resultp - resultm ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){

            BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( gradCol[ j ], d2pdSdC[ j ][ i ] ) );

        }

    }

    return;

}

BOOST_AUTO_TEST_CASE( test_computeDeviatoricReferenceSecondOrderStress ){
    /*!
     * Test the computation of the deviatoric part of the second higher order stress.
     *
     * :param std::ofstream &results: The output file
     */

    variableVector S = { 0.77189588, -0.84417528,  0.95929231,
                        -0.50465708,  0.50576944,  0.05335127,
                         0.81510751,  0.76814059, -0.82146208 };

    variableVector C = { 0.03468919, -0.31275742, -0.57541261,
                        -0.27865312, -0.45844965,  0.52325004,
                        -0.0439162 , -0.80201065, -0.44921044 };

    variableVector answer = { -1.81433044, -2.17117481,  2.726381  ,
                               0.10781401,  0.6746562 , -0.53446584,
                              -0.0255477 ,  0.59634539, -0.39543207 };

    variableType pressure;
    variableVector dpdS, dpdC;
    variableMatrix d2pdSdC;

    errorOut error = tardigradeMicromorphicTools::computeReferenceSecondOrderStressPressure( S, C, pressure, dpdS, dpdC, d2pdSdC );

    BOOST_CHECK( !error );

    variableVector result;

    error = tardigradeMicromorphicTools::computeDeviatoricReferenceSecondOrderStress( S, C, result );

    BOOST_CHECK( !error );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( result, answer ) );

    error = tardigradeMicromorphicTools::computeDeviatoricReferenceSecondOrderStress( S, C, pressure, result );

    BOOST_CHECK( !error );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( result, answer ) );

    //Test the jacobians
    variableVector resultJ;
    variableMatrix dDevSdS, dDevSdC;

    error = tardigradeMicromorphicTools::computeDeviatoricReferenceSecondOrderStress( S, C, resultJ, dDevSdS, dDevSdC );

    BOOST_CHECK( !error );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( resultJ, answer ) );

    variableVector resultJP;
    variableVector _dDevSdSJP, _dDevSdCJP;

    error = tardigradeMicromorphicTools::computeDeviatoricReferenceSecondOrderStress( S, C, pressure, dpdS, dpdC, 
                                                                            resultJP, _dDevSdSJP, _dDevSdCJP );

    variableMatrix dDevSdSJP = tardigradeVectorTools::inflate( _dDevSdSJP, 9, 9 );
    variableMatrix dDevSdCJP = tardigradeVectorTools::inflate( _dDevSdCJP, 9, 9 );

    BOOST_CHECK( !error );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( resultJP, answer ) );

    variableVector resultJ2;
    variableMatrix dDevSdSJ2, dDevSdCJ2, d2DevSdSdC;

    error = tardigradeMicromorphicTools::computeDeviatoricReferenceSecondOrderStress( S, C, resultJ2, dDevSdSJ2, dDevSdCJ2, d2DevSdSdC );

    BOOST_CHECK( !error );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( resultJ2, answer ) );

    variableVector resultJ2P;
    variableVector _dDevSdSJ2P, _dDevSdCJ2P, _d2DevSdSdCP;

    error = tardigradeMicromorphicTools::computeDeviatoricReferenceSecondOrderStress( S, C, pressure, dpdS, dpdC, tardigradeVectorTools::appendVectors( d2pdSdC ),
                                                                            resultJ2P, _dDevSdSJ2P, _dDevSdCJ2P, _d2DevSdSdCP );

    variableMatrix dDevSdSJ2P  = tardigradeVectorTools::inflate( _dDevSdSJ2P, 9,   9 );
    variableMatrix dDevSdCJ2P  = tardigradeVectorTools::inflate( _dDevSdCJ2P, 9,   9 );
    variableMatrix d2DevSdSdCP = tardigradeVectorTools::inflate( _d2DevSdSdCP, 9, 81 );

    BOOST_CHECK( !error );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( resultJ2P, answer ) );

    //Test of dpdS
    constantType eps = 1e-6;
    for ( unsigned int i = 0; i < S.size(); i++ ){
        constantVector delta( S.size(), 0 );
        delta[i] = eps * fabs( S[i] ) + eps;

        variableVector resultp, resultm;

        error = tardigradeMicromorphicTools::computeDeviatoricReferenceSecondOrderStress( S + delta, C, resultp );

        BOOST_CHECK( !error );

        error = tardigradeMicromorphicTools::computeDeviatoricReferenceSecondOrderStress( S - delta, C, resultm );

        BOOST_CHECK( !error );

        constantVector gradCol = ( resultp - resultm ) / ( 2 * delta[i] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){

            BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( gradCol[ j ], dDevSdS[ j ][ i ] ) );

        }

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){

            BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( gradCol[ j ], dDevSdSJP[ j ][ i ] ) );

        }

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){

            BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( gradCol[ j ], dDevSdSJ2[ j ][ i ] ) );

        }

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){

            BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( gradCol[ j ], dDevSdSJ2P[ j ][ i ] ) );

        }

    }

    //Test of dpdC
    for ( unsigned int i = 0; i < C.size(); i++ ){
        constantVector delta( C.size(), 0 );
        delta[i] = eps * fabs( C[i] ) + eps;

        variableVector resultp, resultm;

        error = tardigradeMicromorphicTools::computeDeviatoricReferenceSecondOrderStress( S, C + delta, resultp );

        BOOST_CHECK( !error );

        error = tardigradeMicromorphicTools::computeDeviatoricReferenceSecondOrderStress( S, C - delta, resultm );

        BOOST_CHECK( !error );

        constantVector gradCol = ( resultp - resultm ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){

            BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( gradCol[ j ], dDevSdC[ j ][ i ] ) );

        }

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){

            BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( gradCol[ j ], dDevSdCJP[ j ][ i ] ) );

        }

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){

            BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( gradCol[ j ], dDevSdCJ2[ j ][ i ] ) );

        }

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){

            BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( gradCol[ j ], dDevSdCJ2P[ j ][ i ] ) );

        }

    }

    //Test of d2pdSdC
    for ( unsigned int i = 0; i < C.size(); i++ ){
        constantVector delta( C.size(), 0 );
        delta[i] = eps * fabs( C[i] ) + eps;

        variableMatrix resultp, resultm;

        error = tardigradeMicromorphicTools::computeDeviatoricReferenceSecondOrderStress( S, C + delta, resultJ, resultp, dDevSdCJ2 );

        BOOST_CHECK( !error );

        error = tardigradeMicromorphicTools::computeDeviatoricReferenceSecondOrderStress( S, C - delta, resultJ, resultm, dDevSdCJ2 );

        BOOST_CHECK( !error );

        constantMatrix gradCol = ( resultp - resultm ) / ( 2 * delta[i] );

        for ( unsigned int j = 0; j < 3; j++ ){
            for ( unsigned int k = 0; k < 3; k++ ){
                for ( unsigned int l = 0; l < 3; l++ ){
                    for ( unsigned int m = 0; m < 3; m++ ){

                        BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( gradCol[ 3 * j + k ][ 3 * l + m ],
                                                               d2DevSdSdC[ 3 * j + k ][ 27 * l + 9 * m + 3 * ( int )( i / 3 ) + i % 3 ] ) );

                        BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( gradCol[ 3 * j + k ][ 3 * l + m ],
                                                               d2DevSdSdCP[ 3 * j + k ][ 27 * l + 9 * m + 3 * ( int )( i / 3 ) + i % 3 ] ) );

                    }

                }

            }

        }

    }

    return;

}

BOOST_AUTO_TEST_CASE( test_computeSecondOrderReferenceStressDecomposition ){
    /*!
     * Test the computation of the deviatoric - volumetric (pressure) stress 
     * decomposition.
     *
     * :param std::ofstream &results: The output file.
     */

    variableVector S = { 0.77189588, -0.84417528,  0.95929231,
                        -0.50465708,  0.50576944,  0.05335127,
                         0.81510751,  0.76814059, -0.82146208 };

    variableVector C = { 0.03468919, -0.31275742, -0.57541261,
                        -0.27865312, -0.45844965,  0.52325004,
                        -0.0439162 , -0.80201065, -0.44921044 };

    variableVector deviatoricAnswer = { -1.81433044, -2.17117481,  2.726381  ,
                                         0.10781401,  0.6746562 , -0.53446584,
                                        -0.0255477 ,  0.59634539, -0.39543207 };

    variableType pressureAnswer = -0.20245462701026676;

    variableType pressureResult;
    variableVector deviatoricResult;

    errorOut error = tardigradeMicromorphicTools::computeSecondOrderReferenceStressDecomposition( S, C, deviatoricResult, pressureResult );

    BOOST_CHECK( !error );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( pressureResult, pressureAnswer ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( deviatoricResult, deviatoricAnswer ) );

    //Test the Jacobians
    variableType pressureResultJ;
    variableVector deviatoricResultJ;

    variableVector dPressuredStress, dPressuredRCG;
    variableMatrix dDevStressdStress, dDevStressdRCG;

    error = tardigradeMicromorphicTools::computeSecondOrderReferenceStressDecomposition( S, C, deviatoricResultJ, pressureResultJ,
                                                                               dDevStressdStress, dDevStressdRCG,
                                                                               dPressuredStress, dPressuredRCG );

    BOOST_CHECK( !error );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( pressureResultJ, pressureAnswer ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( deviatoricResultJ, deviatoricAnswer ) );

    variableType pressureResultJ2;
    variableVector deviatoricResultJ2;

    variableVector dPressuredStressJ2, dPressuredRCGJ2;
    variableMatrix dDevStressdStressJ2, dDevStressdRCGJ2;

    variableMatrix d2DevStressdStressdRCGJ2, d2PressuredStressdRCGJ2;

    error = tardigradeMicromorphicTools::computeSecondOrderReferenceStressDecomposition( S, C, deviatoricResultJ2, pressureResultJ2,
                                                                               dDevStressdStressJ2, dDevStressdRCGJ2,
                                                                               dPressuredStressJ2, dPressuredRCGJ2,
                                                                               d2DevStressdStressdRCGJ2, d2PressuredStressdRCGJ2 );

    BOOST_CHECK( !error );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( pressureResultJ2, pressureAnswer ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( deviatoricResultJ2, deviatoricAnswer ) );

    //Test deriatives w.r.t. the stress.
    constantType eps = 1e-6;
    for ( unsigned int i = 0; i < S.size(); i++ ){
        constantVector delta( S.size(), 0 );
        delta[i] = eps * fabs( S[i] ) + eps;

        variableVector vec_resultp, vec_resultm;
        variableType   sca_resultp, sca_resultm;

        error = tardigradeMicromorphicTools::computeSecondOrderReferenceStressDecomposition( S + delta, C, vec_resultp, sca_resultp );

        BOOST_CHECK( !error );

        error = tardigradeMicromorphicTools::computeSecondOrderReferenceStressDecomposition( S - delta, C, vec_resultm, sca_resultm );

        BOOST_CHECK( !error );

        constantVector gradCol = ( vec_resultp - vec_resultm ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){

            BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( gradCol[ j ], dDevStressdStress[ j ][ i ] ) );

            BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( gradCol[ j ], dDevStressdStressJ2[ j ][ i ] ) );

        }

        constantType gradScalar = ( sca_resultp - sca_resultm ) / ( 2 * delta[ i ] );

        BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( gradScalar, dPressuredStress[ i ] ) );

        BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( gradScalar, dPressuredStressJ2[ i ] ) );

    }

    //Test deriatives w.r.t. the right Cauchy-Green deformation tensor
    variableMatrix dDevStressdStressP, dDevStressdRCGP;
    variableVector dPressuredStressP, dPressuredRCGP;
    for ( unsigned int i = 0; i < C.size(); i++ ){
        constantVector delta( C.size(), 0 );
        delta[i] = eps * fabs( C[i] ) + eps;

        variableVector vec_resultp, vec_resultm;
        variableType   sca_resultp, sca_resultm;

        error = tardigradeMicromorphicTools::computeSecondOrderReferenceStressDecomposition( S, C + delta, vec_resultp, sca_resultp );

        BOOST_CHECK( !error );

        error = tardigradeMicromorphicTools::computeSecondOrderReferenceStressDecomposition( S, C - delta, vec_resultm, sca_resultm );

        BOOST_CHECK( !error );

        constantVector gradCol = ( vec_resultp - vec_resultm ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){

            BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( gradCol[ j ], dDevStressdRCG[ j ][ i ] ) );

            BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( gradCol[ j ], dDevStressdRCGJ2[ j ][ i ] ) );

        }

        constantType gradScalar = ( sca_resultp - sca_resultm ) / ( 2 * delta[i] );

        BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( gradScalar, dPressuredRCG[ i ] ) );

        BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( gradScalar, dPressuredRCGJ2[ i ] ) );

        variableMatrix mat_resultp, mat_resultm;

        error = tardigradeMicromorphicTools::computeSecondOrderReferenceStressDecomposition( S, C + delta, deviatoricResultJ, pressureResultJ,
                                                                                   mat_resultp, dDevStressdRCGP, 
                                                                                   vec_resultp, dPressuredRCGP );

        BOOST_CHECK( !error );

        error = tardigradeMicromorphicTools::computeSecondOrderReferenceStressDecomposition( S, C - delta, deviatoricResultJ, pressureResultJ,
                                                                                   mat_resultm, dDevStressdRCGP, 
                                                                                   vec_resultm, dPressuredRCGP );

        BOOST_CHECK( !error );

        constantMatrix gradMat = ( mat_resultp - mat_resultm ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < 3; j++ ){
            for ( unsigned int k = 0; k < 3; k++ ){
                for ( unsigned int l = 0; l < 3; l++ ){
                    for ( unsigned int m = 0; m < 3; m++ ){
                        
                        BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( gradMat[ 3 * j + k ][ 3 * l + m ],
                                                               d2DevStressdStressdRCGJ2[ 3 * j + k ][ 27 * l + 9 * m + 3 * ( int )( i / 3 ) + i % 3 ] ) );

                    }
                }
            }
        }

        gradCol = ( vec_resultp - vec_resultm ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){

            BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( gradCol[ j ], d2PressuredStressdRCGJ2[ j ][ i ] ) );

        }

    }

    return;

}

BOOST_AUTO_TEST_CASE( test_computeReferenceHigherOrderStressPressure ){
    /*!
     * Test the computation of the higher order stress pressure.
     *
     * :param std::ofstream &results: The output file.
     */

    variableVector M = {  1,  2,  3,  4,  5,  6,  7,  8,  9,
                         10, 11, 12, 13, 14, 15, 16, 17, 18,
                         19, 20, 21, 22, 23, 24, 25, 26, 27 };

    variableVector C = { 0.34852835, 0.47540122, 1.11252634,
                         0.47540122, 1.49184663, 1.57435946,
                         1.11252634, 1.57435946, 3.68235756 };

    variableVector answer = { 69.06947837, 73.01858056, 76.96768276 };

    variableVector result;

    errorOut error = tardigradeMicromorphicTools::computeReferenceHigherOrderStressPressure( M, C, result );

    BOOST_CHECK( !error );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( answer, result ) );

    //Test the Jacobians
    
    variableVector resultJ;
    variableMatrix dpdM, dpdC;
    error = tardigradeMicromorphicTools::computeReferenceHigherOrderStressPressure( M, C, resultJ, dpdM, dpdC );

    BOOST_CHECK( !error );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( answer, resultJ ) );

    variableVector resultJ2;
    variableMatrix dpdMJ2, dpdCJ2, d2pdMdC;
    error = tardigradeMicromorphicTools::computeReferenceHigherOrderStressPressure( M, C, resultJ2, dpdMJ2, dpdCJ2, d2pdMdC );

    BOOST_CHECK( !error );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( answer, resultJ2 ) );

    //Test dpdM
    constantType eps = 1e-6;
    for ( unsigned int i = 0; i < M.size(); i++ ){
        constantVector delta( M.size(), 0 );
        delta[i] = eps * fabs( M[i] ) + eps;

        variableVector resultp, resultm;

        error = tardigradeMicromorphicTools::computeReferenceHigherOrderStressPressure( M + delta, C, resultp );

        BOOST_CHECK( !error );

        error = tardigradeMicromorphicTools::computeReferenceHigherOrderStressPressure( M - delta, C, resultm );

        BOOST_CHECK( !error );

        constantVector gradCol = ( resultp - resultm ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){

            BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( gradCol[ j ], dpdM[ j ][ i ] ) );

        }

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){

            BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( gradCol[ j ], dpdMJ2[ j ][ i ] ) );

        }

    }

    //Test dpdC
    for ( unsigned int i = 0; i < C.size(); i++ ){
        constantVector delta( C.size(), 0 );
        delta[i] = eps * fabs( C[i] ) + eps;

        variableVector resultp, resultm;

        error = tardigradeMicromorphicTools::computeReferenceHigherOrderStressPressure( M, C + delta, resultp );

        BOOST_CHECK( !error );

        error = tardigradeMicromorphicTools::computeReferenceHigherOrderStressPressure( M, C - delta, resultm );

        BOOST_CHECK( !error );

        constantVector gradCol = ( resultp - resultm ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){

            BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( gradCol[ j ], dpdC[ j ][ i ] ) );

        }

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){

            BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( gradCol[ j ], dpdCJ2[ j ][ i ] ) );

        }

    }

    //Test d2pdMdC
    for ( unsigned int i = 0; i < C.size(); i++ ){
        constantVector delta( C.size(), 0 );
        delta[i] = eps * fabs( C[i] ) + eps;

        variableMatrix resultp, resultm;

        error = tardigradeMicromorphicTools::computeReferenceHigherOrderStressPressure( M, C + delta, resultJ, resultp, dpdCJ2 );

        BOOST_CHECK( !error );

        error = tardigradeMicromorphicTools::computeReferenceHigherOrderStressPressure( M, C - delta, resultJ, resultm, dpdCJ2 );

        BOOST_CHECK( !error );

        constantMatrix gradCol = ( resultp - resultm ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < 3; j++ ){
            for ( unsigned int k = 0; k < 3; k++ ){
                for ( unsigned int l = 0; l < 3; l++ ){
                    for ( unsigned int m = 0; m < 3; m++ ){
                        BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( gradCol[ j ][ 9 * k + 3 * l + m ], 
                                                               d2pdMdC[ j ][ 81 * k + 27 * l + 9 * m + 3 * (int)( i / 3 ) + i % 3 ] ) );
                    }
                }
            }
        }
    }

    return;

}

BOOST_AUTO_TEST_CASE( test_computeHigherOrderReferenceStressDecomposition ){
    /*!
     * Test the computation of the decomposition of the higher order stress
     * into deviatoric and volumetric (pressure) parts.
     * 
     * :param std::ofstream &results: The output file
     */

    variableVector M = { 0.80732114,  0.79202055,  0.17990022,  0.97454675,  0.703207  ,
                        -0.58236697,  0.53324571, -0.93438873, -0.40650796,  0.14071918,
                         0.66933708, -0.67854069, -0.30317772, -0.93821882,  0.97270622,
                         0.00295302, -0.12441126,  0.30539971, -0.0580227 ,  0.89696105,
                         0.17567709, -0.9592962 ,  0.63535407,  0.95437804, -0.64531877,
                         0.69978907,  0.81327586 };

    variableVector C = { 0.3991656 ,  0.43459435,  0.15398811,
                        -0.20239202, -0.50763359,  0.04756988,
                        -0.0573016 , -0.95939895,  0.2693173 };

    variableVector deviatoricAnswer = { 5.27633771,   1.68477127,  -3.21111472,  13.96704974,
                                        3.29864744, -10.4408607 ,  -4.31689678,  -1.90327486,
                                        3.27369895,  -2.4001685 ,   0.161758  ,   1.24944234,
                                       -6.01118653,  -2.07847616,   5.30384773,   2.46397506,
                                        0.3672135 ,  -1.56198261,  -8.15866529,  -0.72125951,
                                        6.32230906, -18.52878201,  -2.87440491,  14.2858097 ,
                                        4.98150383,   1.82382829,  -3.45626291 };

    variableVector pressureAnswer = { 0.56778037,  0.11342234, -0.43082224 };

    variableVector deviatoricResult, pressureResult;
    errorOut error = tardigradeMicromorphicTools::computeHigherOrderReferenceStressDecomposition( M, C, deviatoricResult, pressureResult );

    BOOST_CHECK( !error );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( deviatoricResult, deviatoricAnswer ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( pressureResult, pressureAnswer ) );

    //Test the Jacobians
    variableVector pressureResultJ;
    variableVector deviatoricResultJ;

    variableMatrix dPressuredStress, dPressuredRCG;
    variableMatrix dDevStressdStress, dDevStressdRCG;

    error = tardigradeMicromorphicTools::computeHigherOrderReferenceStressDecomposition( M, C, deviatoricResultJ, pressureResultJ,
                                                                               dDevStressdStress, dDevStressdRCG,
                                                                               dPressuredStress, dPressuredRCG );

    BOOST_CHECK( !error );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( pressureResultJ, pressureAnswer ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( deviatoricResultJ, deviatoricAnswer ) );

    variableVector pressureResultJ2;
    variableVector deviatoricResultJ2;

    variableMatrix dPressuredStressJ2, dPressuredRCGJ2;
    variableMatrix dDevStressdStressJ2, dDevStressdRCGJ2;

    variableMatrix d2DevStressdStressdRCGJ2, d2PressuredStressdRCGJ2;

    error = tardigradeMicromorphicTools::computeHigherOrderReferenceStressDecomposition( M, C, deviatoricResultJ2, pressureResultJ2,
                                                                               dDevStressdStressJ2, dDevStressdRCGJ2,
                                                                               dPressuredStressJ2, dPressuredRCGJ2,
                                                                               d2DevStressdStressdRCGJ2, d2PressuredStressdRCGJ2 );

    BOOST_CHECK( !error );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( pressureResultJ2, pressureAnswer ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( deviatoricResultJ2, deviatoricAnswer ) );

    //Test deriatives w.r.t. the stress.
    constantType eps = 1e-6;
    for ( unsigned int i = 0; i < M.size(); i++ ){
        constantVector delta( M.size(), 0 );
        delta[i] = eps * fabs( M[i] ) + eps;

        variableVector dev_resultp, dev_resultm;
        variableVector pre_resultp, pre_resultm;

        error = tardigradeMicromorphicTools::computeHigherOrderReferenceStressDecomposition( M + delta, C, dev_resultp, pre_resultp );

        BOOST_CHECK( !error );

        error = tardigradeMicromorphicTools::computeHigherOrderReferenceStressDecomposition( M - delta, C, dev_resultm, pre_resultm );

        BOOST_CHECK( !error );

        constantVector gradCol = ( dev_resultp - dev_resultm ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){

            BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( gradCol[ j ], dDevStressdStress[ j ][ i ] ) );

            BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( gradCol[ j ], dDevStressdStressJ2[ j ][ i ] ) );

        }

        gradCol = ( pre_resultp - pre_resultm ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){

            BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( gradCol[ j ], dPressuredStress[ j ][ i ] ) );

            BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( gradCol[ j ], dPressuredStressJ2[ j ][ i ] ) );

        }

    }

    //Test deriatives w.r.t. the right Cauchy-Green deformation tensor
    variableMatrix dDevStressdStressP, dDevStressdRCGP;
    variableMatrix dPressuredStressP, dPressuredRCGP;
    for ( unsigned int i = 0; i < C.size(); i++ ){
        constantVector delta( C.size(), 0 );
        delta[i] = eps * fabs( C[i] ) + eps;

        variableVector dev_resultp, pre_resultp;
        variableVector dev_resultm, pre_resultm;

        error = tardigradeMicromorphicTools::computeHigherOrderReferenceStressDecomposition( M, C + delta, dev_resultp, pre_resultp );

        BOOST_CHECK( !error );

        error = tardigradeMicromorphicTools::computeHigherOrderReferenceStressDecomposition( M, C - delta, dev_resultm, pre_resultm );

        BOOST_CHECK( !error );

        constantVector gradCol = ( dev_resultp - dev_resultm ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){

            BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( gradCol[ j ], dDevStressdRCG[ j ][ i ] ) );

            BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( gradCol[ j ], dDevStressdRCGJ2[ j ][ i ] ) );

        }

        gradCol = ( pre_resultp - pre_resultm ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){

            BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( gradCol[ j ], dPressuredRCG[ j ][ i ] ) );

            BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( gradCol[ j ], dPressuredRCGJ2[ j ][ i ] ) );

        }

        variableMatrix dev_mat_resultp, pre_mat_resultp;
        variableMatrix dev_mat_resultm, pre_mat_resultm;

        error = tardigradeMicromorphicTools::computeHigherOrderReferenceStressDecomposition( M, C + delta, deviatoricResultJ, pressureResultJ,
                                                                                   dev_mat_resultp, dDevStressdRCGP, 
                                                                                   pre_mat_resultp, dPressuredRCGP );

        BOOST_CHECK( !error );

        error = tardigradeMicromorphicTools::computeHigherOrderReferenceStressDecomposition( M, C - delta, deviatoricResultJ, pressureResultJ,
                                                                                   dev_mat_resultm, dDevStressdRCGP, 
                                                                                   pre_mat_resultm, dPressuredRCGP );

        BOOST_CHECK( !error );

        constantMatrix gradMat = ( dev_mat_resultp - dev_mat_resultm ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < 3; j++ ){
            for ( unsigned int k = 0; k < 3; k++ ){
                for ( unsigned int l = 0; l < 3; l++ ){
                    for ( unsigned int m = 0; m < 3; m++ ){
                        for ( unsigned int n = 0; n < 3; n++ ){
                            for ( unsigned int o = 0; o < 3; o++ ){
                                BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( gradMat[ 9 * j + 3 * k + l ][ 9 * m + 3 * n + o ],
                                                                       d2DevStressdStressdRCGJ2[ 9 * j + 3 * k + l ][ 81 * m + 27 * n + + 9 * o + 3 * ( int )( i / 3 ) + i % 3 ] ) );
                            }
                        }
                    }
                }
            }
        }

        gradMat = ( pre_mat_resultp - pre_mat_resultm ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < 3; j++ ){
            for ( unsigned int k = 0; k < 3; k++ ){
                for ( unsigned int l = 0; l < 3; l++ ){
                    for ( unsigned int m = 0; m < 3; m++ ){
                        BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( gradMat[ j ][ 9 * k + 3 * l + m ],
                                                               d2PressuredStressdRCGJ2[ j ][ 81 * k + 27 * l + 9 * m + 3 * ( int )( i / 3 ) + i % 3 ] ) );
                    }
                }
            }
        }
    }

    return;
    
}

BOOST_AUTO_TEST_CASE( test_computeHigherOrderStressNorm ){
    /*!
     * Test the computation of the special higher order stress norm.
     *
     * :param std::ofstream &results: The output file.
     */

    variableVector M = { 0.80732114,  0.79202055,  0.17990022,  0.97454675,  0.703207  ,
                        -0.58236697,  0.53324571, -0.93438873, -0.40650796,  0.14071918,
                         0.66933708, -0.67854069, -0.30317772, -0.93821882,  0.97270622,
                         0.00295302, -0.12441126,  0.30539971, -0.0580227 ,  0.89696105,
                         0.17567709, -0.9592962 ,  0.63535407,  0.95437804, -0.64531877,
                         0.69978907,  0.81327586 };

    variableVector answer = { 1.82692071, 2.24422424, 1.90780645 };

    variableVector result;

    errorOut error = tardigradeMicromorphicTools::computeHigherOrderStressNorm( M, result );

    BOOST_CHECK( !error );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( answer, result ) );

    //Test the Jacobians
    variableVector resultJ;
    variableMatrix dNormMdMJ;

    error = tardigradeMicromorphicTools::computeHigherOrderStressNorm( M, resultJ, dNormMdMJ );

    BOOST_CHECK( !error );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( answer, resultJ ) );

    variableVector resultJ2;
    variableMatrix dNormMdMJ2, d2NormMdM2J2;

    error = tardigradeMicromorphicTools::computeHigherOrderStressNorm( M, resultJ2, dNormMdMJ2, d2NormMdM2J2 );

    BOOST_CHECK( !error );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( answer, resultJ2 ) );

    //Test dNormMdM
    constantType eps = 1e-6;
    for ( unsigned int i = 0; i < M.size(); i++ ){
        constantVector delta( M.size(), 0 );
        delta[i] = eps * fabs( M[i] ) + eps;

        variableVector resultP, resultM;

        error = tardigradeMicromorphicTools::computeHigherOrderStressNorm( M + delta, resultP );

        BOOST_CHECK( !error );

        error = tardigradeMicromorphicTools::computeHigherOrderStressNorm( M - delta, resultM );

        BOOST_CHECK( !error );

        variableVector gradCol = ( resultP - resultM ) / ( 2 * delta[i] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( gradCol[ j ], dNormMdMJ[ j ][ i ] ) );
        }

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( gradCol[ j ], dNormMdMJ2[ j ][ i ] ) );
        }

        variableMatrix derP, derM;

        error = tardigradeMicromorphicTools::computeHigherOrderStressNorm( M + delta, resultP, derP );

        BOOST_CHECK( !error );

        error = tardigradeMicromorphicTools::computeHigherOrderStressNorm( M - delta, resultM, derM );

        BOOST_CHECK( !error );

        variableMatrix gradMat = ( derP - derM ) / ( 2 * delta[i] );

        unsigned int n = ( int )( i / 9 );
        unsigned int o = ( int )( ( i - 9 * n ) / 3 );
        unsigned int p = ( i - 9 * n - 3 * o ) % 3;

        for ( unsigned int j = 0; j < 3; j++ ){
            for ( unsigned int k = 0; k < 3; k++ ){
                for ( unsigned int l = 0; l < 3; l++ ){
                    for ( unsigned int m = 0; m < 3; m++ ){
                        BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( gradMat[ j ][ 9 * k + 3 * l + m ],
                                                               d2NormMdM2J2[ j ][ 243 * k + 81 * l + 27 * m + 9 * n + 3 * o + p ] ) );
                    }
                }
            }
        }
    }

    return;

}

BOOST_AUTO_TEST_CASE( test_assembleDeformationGradient ){
    /*!
     * Test the assembly of the deformation gradient from the gradient of the displacement.
     *
     * :param std::ofstream &results: The output file.
     */

    variableMatrix displacementGradient = { { 1, 2, 3 },
                                            { 4, 5, 6 },
                                            { 7, 8, 9 } };

    variableVector answer = { 2, 2, 3, 4, 6, 6, 7, 8, 10 };

    variableVector result;

    errorOut error = tardigradeMicromorphicTools::assembleDeformationGradient( displacementGradient, result );

    BOOST_CHECK( !error );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( answer, result ) );

    //Tests of the Jacobians
    variableVector resultJ;
    variableMatrix dFdGradU;
    error = tardigradeMicromorphicTools::assembleDeformationGradient( displacementGradient, resultJ, dFdGradU );

    BOOST_CHECK( !error );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( answer, resultJ ) );

    //Test the Jacobian w.r.t. the displacement gradient
    constantType eps = 1e-6;
    for ( unsigned int i = 0; i < 9; i++ ){
        constantMatrix delta = constantMatrix( 3, constantVector( 3, 0 ) );
        delta[ ( int )( i / 3) ][ i % 3 ] = eps * fabs( displacementGradient[ ( int )( i / 3 ) ][ i % 3 ] ) + eps;

        variableVector P, M;

        error = tardigradeMicromorphicTools::assembleDeformationGradient( displacementGradient + delta, P );

        BOOST_CHECK( !error );

        error = tardigradeMicromorphicTools::assembleDeformationGradient( displacementGradient - delta, M );

        BOOST_CHECK( !error );

        variableVector gradCol = ( P - M ) / ( 2 * delta[ ( int )( i / 3 ) ][ i % 3 ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( gradCol[ j ], dFdGradU[ j ][ i ] ) );
        }
    }

    return;

}

BOOST_AUTO_TEST_CASE( test_assembleMicroDeformation ){
    /*!
     * Test the assembly of the micro deformation from the micro displacement
     *
     * :param std::ofstream &results: The output file.
     */

    variableVector microDisplacement = { 1, 2, 3, 4, 5, 6, 7, 8, 9 };

    variableVector answer = { 2, 2, 3, 4, 6, 6, 7, 8, 10 };

    variableVector result;

    errorOut error = tardigradeMicromorphicTools::assembleMicroDeformation( microDisplacement, result );

    BOOST_CHECK( !error );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( answer, result ) );

    //Test the Jacobians
    variableVector resultJ;
    variableMatrix dChidPhi;

    error = tardigradeMicromorphicTools::assembleMicroDeformation( microDisplacement, resultJ, dChidPhi );

    BOOST_CHECK( !error );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( answer, resultJ ) );

    //Test the jacobians w.r.t. phi
    constantType eps = 1e-6;
    for ( unsigned int i = 0; i < microDisplacement.size(); i++ ){
        constantVector delta( microDisplacement.size(), 0 );
        delta[i] = eps * fabs( microDisplacement[ i ] ) + eps;

        variableVector P, M;

        error = tardigradeMicromorphicTools::assembleMicroDeformation( microDisplacement + delta, P );

        BOOST_CHECK( !error );

        error = tardigradeMicromorphicTools::assembleMicroDeformation( microDisplacement - delta, M );

        BOOST_CHECK( !error );

        variableVector gradCol = ( P - M ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( gradCol[ j ], dChidPhi[ j ][ i ] ) );
        }
    }

    return;

}

BOOST_AUTO_TEST_CASE( test_assembleGradientMicroDeformation ){
    /*!
     * Test the assembly of the micro deformation gradient
     *
     * :param std::ofstream &results: The output file.
     */

    variableMatrix gradientMicroDisplacement = { {  1,  2,  3 },
                                                 {  4,  5,  6 },
                                                 {  7,  8,  9 },
                                                 { 10, 11, 12 },
                                                 { 13, 14, 15 },
                                                 { 16, 17, 18 },
                                                 { 19, 20, 21 },
                                                 { 22, 23, 24 },
                                                 { 25, 26, 27 } };

    variableVector answer = { 1,  2,  3,  4,  5,  6,  7,  8,  9,
                             10, 11, 12, 13, 14, 15, 16, 17, 18,
                             19, 20, 21, 22, 23, 24, 25, 26, 27 };

    variableVector result;

    errorOut error = tardigradeMicromorphicTools::assembleGradientMicroDeformation( gradientMicroDisplacement, result );

    BOOST_CHECK( !error );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( answer, result ) );

    //Test the Jacobians
    variableVector resultJ;
    variableMatrix dGradChidGradPhi;

    error = tardigradeMicromorphicTools::assembleGradientMicroDeformation( gradientMicroDisplacement, resultJ, dGradChidGradPhi );

    BOOST_CHECK( !error );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( answer, resultJ ) );

    //Test the Jacobian w.r.t. the micro displacement gradient
    constantType eps = 1e-6;
    for ( unsigned int i = 0; i < 27; i++ ){
        constantMatrix delta( 9, constantVector( 3, 0 ) );
        delta[ ( int )( i / 3 ) ][ i % 3 ] = eps * fabs( gradientMicroDisplacement[ ( int )( i / 3 ) ][ i % 3 ] ) + eps;

        variableVector P, M;

        error = tardigradeMicromorphicTools::assembleGradientMicroDeformation( gradientMicroDisplacement + delta, P );

        BOOST_CHECK( !error );

        error = tardigradeMicromorphicTools::assembleGradientMicroDeformation( gradientMicroDisplacement - delta, M );

        BOOST_CHECK( !error );

        variableVector gradCol = ( P - M ) / ( 2 * delta[ ( int )( i / 3 ) ][ i % 3 ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( gradCol[ j ], dGradChidGradPhi[ j ][ i ] ) );
        }
    }

    return;

}

BOOST_AUTO_TEST_CASE( test_pullBackCauchyStress ){
    /*!
     * Test the pull back operation on the Cauchy stress
     *
     * :param std::ofstream &results: The output file.
     */

    variableVector cauchyStress = { 1, 2, 3, 4, 5, 6, 7, 8, 9 };
    variableVector deformationGradient  = { -1.08831037, -0.66333427, -0.48239487,
                                            -0.904554  ,  1.28942848, -0.02156112,
                                            -0.08464824, -0.07730218,  0.86415668 };

    variableVector answer = { -12.84536712,  -0.74608931,  13.40350133,
                                0.10367879,   0.04818336,   0.22121131,
                               20.04248815,   1.49351009, -18.33987341 };

    variableVector result, resultPF;

    errorOut error = tardigradeMicromorphicTools::pullBackCauchyStress( cauchyStress, deformationGradient, result );

    BOOST_CHECK( !error );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( result, answer ) );

    error = tardigradeMicromorphicTools::pushForwardPK2Stress( result, deformationGradient, resultPF );

    BOOST_CHECK( !error );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( resultPF, cauchyStress ) );

    //Test the Jacobians
    variableVector resultJ;
    variableMatrix dPK2dCauchy, dPK2dF;

    error = tardigradeMicromorphicTools::pullBackCauchyStress( cauchyStress, deformationGradient, resultJ, dPK2dCauchy, dPK2dF );

    BOOST_CHECK( !error );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( resultJ, answer ) );

    //Test the Jacobian w.r.t. the Cauchy stress
    constantType eps = 1e-6;
    for ( unsigned int i = 0; i < cauchyStress.size(); i++ ){
        constantVector delta( cauchyStress.size(), 0 );
        delta[i] = eps * fabs( cauchyStress[ i ] ) + eps;

        variableVector P, M;

        error = tardigradeMicromorphicTools::pullBackCauchyStress( cauchyStress + delta, deformationGradient, P );

        BOOST_CHECK( !error );

        error = tardigradeMicromorphicTools::pullBackCauchyStress( cauchyStress - delta, deformationGradient, M );

        BOOST_CHECK( !error );

        variableVector gradCol = ( P - M ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( gradCol[ j ], dPK2dCauchy[ j ][ i ] ) );
        }
    }

    //Test the Jacobian w.r.t. the deformation gradient
    for ( unsigned int i = 0; i < deformationGradient.size(); i++ ){
        constantVector delta( deformationGradient.size(), 0 );
        delta[i] = eps * fabs( deformationGradient[ i ] ) + eps;

        variableVector P, M;

        error = tardigradeMicromorphicTools::pullBackCauchyStress( cauchyStress, deformationGradient + delta, P );

        BOOST_CHECK( !error );

        error = tardigradeMicromorphicTools::pullBackCauchyStress( cauchyStress, deformationGradient - delta, M );

        BOOST_CHECK( !error );

        variableVector gradCol = ( P - M ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( gradCol[ j ], dPK2dF[ j ][ i ] ) );
        }
    }

    return;

}

BOOST_AUTO_TEST_CASE( test_pullBackMicroStress ){
    /*!
     * Test the pull back operation on the micro stress
     *
     * :param std::ofstream &results: The output file.
     */

    variableVector microStress = { 1, 2, 3, 4, 5, 6, 7, 8, 9 };
    variableVector deformationGradient  = { -1.08831037, -0.66333427, -0.48239487,
                                            -0.904554  ,  1.28942848, -0.02156112,
                                            -0.08464824, -0.07730218,  0.86415668 };

    variableVector answer = { -12.84536712,  -0.74608931,  13.40350133,
                                0.10367879,   0.04818336,   0.22121131,
                               20.04248815,   1.49351009, -18.33987341 };

    variableVector result, resultPF;

    errorOut error = tardigradeMicromorphicTools::pullBackMicroStress( microStress, deformationGradient, result );

    BOOST_CHECK( !error );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( result, answer ) );

    error = tardigradeMicromorphicTools::pushForwardReferenceMicroStress( result, deformationGradient, resultPF );

    BOOST_CHECK( !error );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( resultPF, microStress ) );

    //Test the Jacobians
    variableVector resultJ;
    variableMatrix dSigmads, dSigmadF;

    error = tardigradeMicromorphicTools::pullBackMicroStress( microStress, deformationGradient, resultJ, dSigmads, dSigmadF );

    BOOST_CHECK( !error );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( resultJ, answer ) );

    //Test the Jacobian w.r.t. the micro stress
    constantType eps = 1e-6;
    for ( unsigned int i = 0; i < microStress.size(); i++ ){
        constantVector delta( microStress.size(), 0 );
        delta[i] = eps * fabs( microStress[ i ] ) + eps;

        variableVector P, M;

        error = tardigradeMicromorphicTools::pullBackMicroStress( microStress + delta, deformationGradient, P );

        BOOST_CHECK( !error );

        error = tardigradeMicromorphicTools::pullBackMicroStress( microStress - delta, deformationGradient, M );

        BOOST_CHECK( !error );

        variableVector gradCol = ( P - M ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( gradCol[ j ], dSigmads[ j ][ i ] ) );
        }
    }

    //Test the Jacobian w.r.t. the deformation gradient
    for ( unsigned int i = 0; i < deformationGradient.size(); i++ ){
        constantVector delta( deformationGradient.size(), 0 );
        delta[i] = eps * fabs( deformationGradient[ i ] ) + eps;

        variableVector P, M;

        error = tardigradeMicromorphicTools::pullBackMicroStress( microStress, deformationGradient + delta, P );

        BOOST_CHECK( !error );

        error = tardigradeMicromorphicTools::pullBackMicroStress( microStress, deformationGradient - delta, M );

        BOOST_CHECK( !error );

        variableVector gradCol = ( P - M ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( gradCol[ j ], dSigmadF[ j ][ i ] ) );
        }
    }

    return;

}

BOOST_AUTO_TEST_CASE( test_pullBackHigherOrderStress ){
    /*!
     * Test the push forward operation on the higher order stress tensor.
     *
     * :param std::ofstream &results: The output file.
     */

    variableVector higherOrderStress = {  1,  2,  3,  4,  5,  6,  7,  8,  9,
                                         10, 11, 12, 13, 14, 15, 16, 17, 18,
                                         19, 20, 21, 22, 23, 24, 25, 26, 27 };

    variableVector deformationGradient = { 0.29524861, -0.25221581, -1.60534711,
                                          -0.08817703,  0.28447808,  0.06451703,
                                           0.39849454,  0.38681512,  0.87840084 };

    variableVector microDeformation = { -0.25781969, -0.39826899, -0.79493259,
                                         0.38104724, -0.00830511, -0.51985409,
                                        -0.36415661, -0.6871168 ,  0.54018665 };

    variableVector answer = {   12.16637846,  -17.55206802,   -2.44609546,   46.49442604,
                               -66.93990198,   -9.1624221 ,  -11.22600619,   16.12630462,
                                 2.16292349,   36.56411317,  -52.34132435,   -6.7949869 ,
                               133.6526827 , -190.51226706,  -23.73373276,  -30.65330414,
                                43.46860157,    5.13634497,   -6.18698921,    8.71843323,
                                 0.96159717,  -20.55921994,   28.62926331,    2.72986588,
                                 4.14350171,   -5.66535471,   -0.40778148 };

    variableVector result;

    errorOut error = tardigradeMicromorphicTools::pullBackHigherOrderStress( higherOrderStress, deformationGradient,
                                                                   microDeformation, result );

    BOOST_CHECK( !error );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( answer, result ) );

    variableVector resultPF;

    error = tardigradeMicromorphicTools::pushForwardHigherOrderStress( result, deformationGradient, microDeformation, resultPF );

    BOOST_CHECK( !error );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( higherOrderStress, resultPF ) );

    //Test the Jacobian
    variableVector resultJ;
    variableMatrix dMdm, dMdF, dMdChi;

    error = tardigradeMicromorphicTools::pullBackHigherOrderStress( higherOrderStress, deformationGradient,
                                                          microDeformation, resultJ, dMdm, dMdF, dMdChi );

    BOOST_CHECK( !error );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( answer, resultJ ) );

    //Test the jacobian w.r.t. the higher order stress
    constantType eps = 1e-6;
    for ( unsigned int i = 0; i < higherOrderStress.size(); i++ ){
        constantVector delta( higherOrderStress.size(), 0 );
        delta[ i ] = eps * fabs( higherOrderStress[ i ] ) + eps;

        variableVector P, M;

        error = tardigradeMicromorphicTools::pullBackHigherOrderStress( higherOrderStress + delta, deformationGradient,
                                                              microDeformation, P );

        BOOST_CHECK( !error );    

        error = tardigradeMicromorphicTools::pullBackHigherOrderStress( higherOrderStress - delta, deformationGradient,
                                                              microDeformation, M );

        BOOST_CHECK( !error );    

        variableVector gradCol = ( P - M ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( gradCol[ j ], dMdm[ j ][ i ] ) );
        }
    }

    //Test the jacobian w.r.t. the deformation gradient
    for ( unsigned int i = 0; i < deformationGradient.size(); i++ ){
        constantVector delta( deformationGradient.size(), 0 );
        delta[ i ] = eps * fabs( deformationGradient[ i ] ) + eps;

        variableVector P, M;

        error = tardigradeMicromorphicTools::pullBackHigherOrderStress( higherOrderStress, deformationGradient + delta,
                                                              microDeformation, P );

        BOOST_CHECK( !error );
    
        error = tardigradeMicromorphicTools::pullBackHigherOrderStress( higherOrderStress, deformationGradient - delta,
                                                              microDeformation, M );
   
        BOOST_CHECK( !error ); 

        variableVector gradCol = ( P - M ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( gradCol[ j ], dMdF[ j ][ i ] ) );
        }
    }

    //Test the jacobian w.r.t. the micro deformation
    for ( unsigned int i = 0; i < microDeformation.size(); i++ ){
        constantVector delta( microDeformation.size(), 0 );
        delta[ i ] = eps * fabs( microDeformation[ i ] ) + eps;

        variableVector P, M;

        error = tardigradeMicromorphicTools::pullBackHigherOrderStress( higherOrderStress, deformationGradient,
                                                              microDeformation + delta, P );

        BOOST_CHECK( !error );    

        error = tardigradeMicromorphicTools::pullBackHigherOrderStress( higherOrderStress, deformationGradient,
                                                              microDeformation - delta, M );

        BOOST_CHECK( !error );    

        variableVector gradCol = ( P - M ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( gradCol[ j ], dMdChi[ j ][ i ] ) );
        }
    }

    return;

}
