/*
  This target aims to be an example of usage of the current PZ version
*/
#include <TPZGeoMeshTools.h>
#include "TPZAnalyticSolution.h"
#include <TPZGmshReader.h>
#include "TPZCompMeshTools.h"
#include "pzlog.h"
#include "TPZLinearAnalysis.h"
#include "TPZMultiphysicsCompMesh.h"
#include "pzbuildmultiphysicsmesh.h"
#include "TPZHDivApproxCreator.h"

#include "pzshapequad.h"
#include "pzshapetriang.h"
#include "pzshapecube.h"
#include "pzshapetetra.h"
#include "TPZTimer.h"
#include "TPZSimpleTimer.h"
#include "TPZVTKGenerator.h"
#include "TPZRefPattern.h"
#include "tpzgeoelrefpattern.h"
#include "TPZRefPatternDataBase.h"
#include "DarcyFlow/TPZMixedDarcyFlow.h"
#include "Elasticity/TPZMixedElasticityND.h"
#include "Elasticity/TPZElasticity2D.h"
#include <TPZSSpStructMatrix.h> //symmetric sparse matrix storage
#include <pzstepsolver.h> //for TPZStepSolver
#include "pzblockdiag.h"
#include "pzbdstrmatrix.h"
#include "TPZVTKGeoMesh.h"
#include "TPZH1ApproxCreator.h"
#include <iostream>
#include <fstream>

std::ofstream printerrors("results_errors.txt",std::ios::app);

enum EMatid  {ENone, EDomain, EDomain2, ELeft, EBottom, ERight, ETop, ECircle};

/**
   @brief Reads the test mesh from gmsh
   @param[in] file_name the .msh mesh file.
*/
TPZGeoMesh* ReadMeshFromGmsh(std::string file_name);

void InsertMaterials(TPZCompMesh* cmesh);

int main() {

    #ifdef PZ_LOG 
    TPZLogger::InitializePZLOG();
#endif

    TPZGeoMesh* gmesh = ReadMeshFromGmsh("../Atividades Felipe/Geograde/Geograde.msh");

    {
        // Prints gmesh mesh properties
        std::string vtk_name = "geoMesh.vtk";
        std::ofstream vtkfile(vtk_name.c_str());
        TPZVTKGeoMesh::PrintGMeshVTK(gmesh, vtkfile, true);
    
    }
    TPZCompMesh *cmesh = new TPZCompMesh(gmesh); 
    cmesh->SetDefaultOrder(1);

    //Insert Materials
    InsertMaterials(cmesh);

    cmesh->SetAllCreateFunctionsContinuous();
    cmesh->AutoBuild();

    //Prints gmesh mesh properties
    std::string vtk_name = "compMesh.vtk";
    std::ofstream vtkfile(vtk_name.c_str());

    TPZVTKGeoMesh::PrintCMeshVTK(cmesh, vtkfile, true);

    //Create the analysis environment
    TPZLinearAnalysis an(cmesh,RenumType::ENone);

    // Solve problem
    constexpr int nThreads{20};
    TPZSSpStructMatrix<STATE,TPZStructMatrixOR<STATE>> stiffness(cmesh);   
    stiffness.SetNumThreads(nThreads);
    an.SetStructuralMatrix(stiffness);

    ///Setting a direct solver
    TPZStepSolver<STATE> step;
    step.SetDirect(ELDLt);//ELU //ECholesky // ELDLt
    an.SetSolver(step);
    
    //Assemble and solve the problem.
    an.Run();
  
    //Printing results in a vtk file
    std::string plotfile = "result";
    //Fields to be printed
    TPZVec<std::string> fields={"Displacement","SigmaX","SigmaY","TauXY"};
    
    int res=0;
    auto vtk = TPZVTKGenerator(an.Mesh(), fields, plotfile, res);
    vtk.Do(); 

    return 0;
}



void InsertMaterials(TPZCompMesh *cmesh){

    //Material do Solo
    double E1 = 1.;
    double nu1 = .1;

    TPZElasticity2D* matElast = new TPZElasticity2D(EDomain);
    matElast -> SetElasticity(E1,nu1);

    cmesh->InsertMaterialObject(matElast);

    //Material da Geograde.
    double E2 = 1000;
    double nu2 =.49;

    TPZElasticity2D* matElast2 = new TPZElasticity2D(EDomain2);
    matElast2 -> SetElasticity(E2,nu2);

    cmesh->InsertMaterialObject(matElast2);

    //Problema vetorial:

    /*
    -> Domain = mat solo; 
    -> Domain2 = mat geograde; 
    -> Left = laterais e fundo fixos; 
    -> Top = solo sob tensão; 
    -> Right = geograde
    */

    TPZFMatrix<STATE> val1(2,2,0.);
    TPZManVector<STATE> val2(2,0.);
    TPZManVector<STATE> val3(2,0.);

    val2[0]=8;
    val2[1]=0;

    TPZBndCondT<STATE> *BCond1 = matElast->CreateBC(matElast, ELeft, 3, val1, val2);
    cmesh->InsertMaterialObject(BCond1);

    val2[0]=0;
    val2[1]=0;

    TPZBndCondT<STATE> *BCond3 = matElast->CreateBC(matElast, ERight, 0, val1, val2);
    cmesh->InsertMaterialObject(BCond3);

    val2[0]=0;
    val2[1]=-2;

    TPZBndCondT<STATE> *BCond2 = matElast->CreateBC(matElast, ETop, 1, val1, val2);
    cmesh->InsertMaterialObject(BCond2);

}



TPZGeoMesh* ReadMeshFromGmsh(std::string file_name)
{
    //read mesh from gmsh
    TPZGeoMesh *gmesh;
    gmesh = new TPZGeoMesh();
    {
        TPZGmshReader reader;
        // essa interface permite voce mapear os nomes dos physical groups para
        // o matid que voce mesmo escolher
        TPZManVector<std::map<std::string,int>,4> stringtoint(4);
        stringtoint[2]["Domain"] = EDomain;
        stringtoint[2]["Domain2"] = EDomain2;
        stringtoint[1]["Left"] = ELeft;
        stringtoint[1]["Right"] = ERight;
        stringtoint[1]["Top"] = ETop;
        stringtoint[1]["Bottom"] = EBottom;
        stringtoint[1]["Circle"] = ECircle;

        reader.SetDimNamePhysical(stringtoint);
        reader.GeometricGmshMesh(file_name,gmesh);

        //Prints gmesh mesh properties
        std::string vtk_name = "geoMesh.vtk";
        std::ofstream vtkfile(vtk_name.c_str());

        TPZVTKGeoMesh::PrintGMeshVTK(gmesh, vtkfile, true);
    }

    return gmesh;
}