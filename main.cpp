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
#include <cmath>
#include "Poisson/TPZMatPoisson.h"
#include "TPZGeoMeshTools.h"
#include "Projection/TPZL2Projection.h"

auto forcefunction = [](const TPZVec<REAL> &loc,
    TPZVec<STATE>&u){
    const auto &x=loc[0];
    const auto &y=loc[1];
    const auto &z=loc[2];

    u[0] = 1.;
};


std::ofstream printerrors("results_errors.txt",std::ios::app);

enum EMatid  {ENone, EDomain, ELeft, EBottom, ERight, ETop, ECircle};

/**
   @brief Reads the test mesh from gmsh
   @param[in] file_name the .msh mesh file.
*/
TPZGeoMesh* ReadMeshFromGmsh(std::string file_name);

void InsertMaterials(TPZCompMesh* cmesh);
void InsertMaterialsL2Projection(TPZCompMesh* cmesh);

int main() {

    #ifdef PZ_LOG 
    TPZLogger::InitializePZLOG();
#endif

    // TPZGeoMesh* gmesh = ReadMeshFromGmsh("../carro.msh");


        TPZVec<int> matids(3);
        matids[0] = EDomain; 
        matids[1] = ELeft; 
        matids[2] = ERight; 

        auto gmeshx = TPZGeoMeshTools::CreateGeoMesh1D(0., 2., 5,
                    matids, true);
        auto gmeshxm = TPZGeoMeshTools::CreateGeoMesh1D(0., 2., 5,
                    matids, true);
        auto gmeshy = TPZGeoMeshTools::CreateGeoMesh1D(0., 2., 10,
                    matids, true);
        auto gmeshym = TPZGeoMeshTools::CreateGeoMesh1D(0., 2., 10,
                    matids, true);

    {
        // Prints gmesh mesh properties
        std::string vtk_namex = "geoMeshx.vtk";
        std::ofstream vtkfile(vtk_namex.c_str());
        TPZVTKGeoMesh::PrintGMeshVTK(gmeshx, vtkfile, true);
    }

    {
        std::string vtk_name = "geoMeshy.vtk";
        std::ofstream vtkfile(vtk_name.c_str());
        TPZVTKGeoMesh::PrintGMeshVTK(gmeshy, vtkfile, true);
    }
    

    {
        // Prints gmesh mesh properties
        // std::string vtk_name = "geoMesh.vtk";
        // std::ofstream vtkfile(vtk_name.c_str());
        // TPZVTKGeoMesh::PrintGMeshVTK(gmesh, vtkfile, true);
    
    }

    TPZCompMesh *cmeshxk = new TPZCompMesh(gmeshx); 
    cmeshxk->SetDefaultOrder(1);
    cmeshxk->SetDimModel(1);
    TPZCompMesh *cmeshxm = new TPZCompMesh(gmeshxm);
    cmeshxm->SetDefaultOrder(1);
    cmeshxm->SetDimModel(1);

    TPZCompMesh *cmeshyk = new TPZCompMesh(gmeshy); 
    cmeshyk->SetDefaultOrder(1);
    cmeshyk->SetDimModel(1);
    TPZCompMesh *cmeshym = new TPZCompMesh(gmeshym);
    cmeshym->SetDefaultOrder(1);
    cmeshym->SetDimModel(1);

    //Insert Materials
    InsertMaterials(cmeshxk);
    InsertMaterialsL2Projection(cmeshxm);

    cmeshxk->SetAllCreateFunctionsContinuous();
    cmeshxk->AutoBuild();
    cmeshxm->SetAllCreateFunctionsContinuous();
    cmeshxm->AutoBuild();

    InsertMaterials(cmeshyk);
    InsertMaterialsL2Projection(cmeshym);

    cmeshyk->SetAllCreateFunctionsContinuous();
    cmeshyk->AutoBuild();
    cmeshym->SetAllCreateFunctionsContinuous();
    cmeshym->AutoBuild();

    TPZSSpStructMatrix<STATE,TPZStructMatrixOR<STATE>> Kx(cmeshxk), Mx(cmeshxm);
    int nequationx = cmeshxk->NEquations(); 
    TPZFMatrix<double> Fx(nequationx,1,0.) , Zx(nequationx,1,0.);
    TPZFMatrix<STATE> matAuxx(nequationx,nequationx,0.);
    TPZFMatrix<STATE> matAuxxM(nequationx,nequationx,0.);
    TPZAutoPointer<TPZGuiInterface> guiInterfacex;
    Kx.Assemble(matAuxx,Fx,guiInterfacex);
    matAuxx.Print("Kx = ", std::cout,EMathematicaInput);
    Fx.Print("Fx = ", std::cout,EMathematicaInput);
    Mx.Assemble(matAuxxM,Fx,guiInterfacex);
    matAuxxM.Print("Mx = ", std::cout,EMathematicaInput);
    

    TPZSSpStructMatrix<STATE,TPZStructMatrixOR<STATE>> Ky(cmeshyk), My(cmeshym);
    int nequationy = cmeshyk->NEquations(); 
    TPZFMatrix<double> Fy(nequationy,1,0.) , Zy(nequationy,1,0.);  

    TPZFMatrix<STATE> matAuxy(nequationy,nequationy,0.);
    TPZFMatrix<STATE> matAuxyM(nequationy,nequationy,0.);
    Ky.Assemble(matAuxy,Fy,guiInterfacex);
    matAuxy.Print("Ky = ", std::cout,EMathematicaInput);
    Fy.Print("Fy = ", std::cout,EMathematicaInput);
    My.Assemble(matAuxyM,Fy,guiInterfacex);
    matAuxyM.Print("My = ", std::cout,EMathematicaInput);

    int modos = 5;
    double e = 10;
    int limite = std::pow(10,-10);

    // for (int i = 0; i < modos; i++){
    //     while( e > limite ){
    //         // Resolve βxi
    //         // Resolve δxi
    //         // Resolve Zx 
    //         // Resolve (Kx + Mx)Φpn = Fx − Z (problema 12)
    //         // Resolve βyi
    //         // Resolve δyi
    //         // Resolve Zy 
    //         // Resolve (Ky + My)Ψpn = Fy − Z (problema 22)
    //         // Calcula e
    //         e = std::pow(10,-11);
    //     }
    // }


    //Prints gmesh mesh properties
    // {
    //     std::string vtk_name = "compMesh.vtk";
    //     std::ofstream vtkfile(vtk_name.c_str());
    //     TPZVTKGeoMesh::PrintCMeshVTK(cmeshxk, vtkfile, true);
    // }

    // {
    //     std::string vtk_name = "compMeshy.vtk";
    //     std::ofstream vtkfile(vtk_name.c_str());
    //     TPZVTKGeoMesh::PrintCMeshVTK(cmeshyk, vtkfile, true);
    // }
    // //Create the analysis environment
    // TPZLinearAnalysis anx(cmeshxk,RenumType::ENone);
    // TPZLinearAnalysis any(cmeshyk,RenumType::ENone);

    // // Solve problem
    // constexpr int nThreads{20};
    // TPZSSpStructMatrix<STATE,TPZStructMatrixOR<STATE>> stiffnessx(cmeshxk);   
    // stiffnessx.SetNumThreads(nThreads);
    // anx.SetStructuralMatrix(stiffnessx);

    // TPZSSpStructMatrix<STATE,TPZStructMatrixOR<STATE>> stiffnessy(cmeshyk);   
    // any.SetStructuralMatrix(stiffnessy);

    // ///Setting a direct solver
    // TPZStepSolver<STATE> step;
    // step.SetDirect(ELDLt);//ELU //ECholesky // ELDLt
    // anx.SetSolver(step);
    // any.SetSolver(step);
    
    // //Assemble and solve the problem.
    // anx.Run();
    // auto solx = anx.Solution();
    // any.Run();
    // auto soly = any.Solution();
    // solx.Print("Solution in x = ",std::cout,EFormatted);
    // soly.Print("Solution in y = ",std::cout,EFormatted);
    
    // //Printing results in a vtk file
    // std::string plotfilex = "resultx";
    // std::string plotfiley = "resulty";
    // //Fields to be printed
    // TPZVec<std::string> fields={"Solution","Derivative"};
    
    // int res=1;
    // auto vtkx = TPZVTKGenerator(anx.Mesh(), fields, plotfilex, res);
    // auto vtky = TPZVTKGenerator(any.Mesh(), fields, plotfiley, res);
    // vtkx.Do(); 
    // vtky.Do();

    return 0;
}



void InsertMaterials(TPZCompMesh *cmeshx){

    // double E = 1.;
    // double nu = .1;

    // TPZElasticity2D* matElast = new TPZElasticity2D(EDomain);
    // matElast -> SetElasticity(E,nu);

    // cmesh->InsertMaterialObject(matElast);

    // TPZFMatrix<STATE> val1(2,2,0.);
    // TPZManVector<STATE> val2(2,0.);
    
    // TPZBndCondT<STATE> *BCond1 = matElast->CreateBC(matElast, ECircle, 0, val1, val2);
    // cmesh->InsertMaterialObject(BCond1);

    // //val2[0]=1.;
    // val2[1] = -10;
    // TPZBndCondT<STATE> *BCond2 = matElast->CreateBC(matElast, ERight, 1, val1, val2);
    // cmesh->InsertMaterialObject(BCond2);

    // val2[1] = -20;
    // TPZBndCondT<STATE> *BCond3 = matElast->CreateBC(matElast, ELeft, 1, val1, val2);
    // cmesh->InsertMaterialObject(BCond3);


    TPZMatPoisson<STATE>* matpoisson = new TPZMatPoisson(EDomain,1);
    cmeshx->InsertMaterialObject(matpoisson);
    // matpoisson -> SetScaleFactor(2);
    matpoisson->SetForcingFunction(forcefunction,1);

    TPZFMatrix<STATE> val1(1,1,0.);
    TPZManVector<STATE> val2(1,0.);
    
    TPZBndCondT<STATE> *BCond1 = matpoisson -> CreateBC(matpoisson , ELeft, 0, val1, val2);
    cmeshx->InsertMaterialObject(BCond1);

    // val2[0] = 0;
    TPZBndCondT<STATE> *BCond2 = matpoisson -> CreateBC(matpoisson , ERight, 0, val1, val2);
    cmeshx->InsertMaterialObject(BCond2);

}


void InsertMaterialsL2Projection(TPZCompMesh *cmeshxk){

    TPZL2Projection<STATE>* matpoisson = new TPZL2Projection(EDomain,1);
    cmeshxk->InsertMaterialObject(matpoisson);
    // matpoisson -> SetScaleFactor(2);

    TPZFMatrix<STATE> val1(1,1,0.);
    TPZManVector<STATE> val2(1,0.);
    
    TPZBndCondT<STATE> *BCond1 = matpoisson -> CreateBC(matpoisson , ELeft, 0, val1, val2);
    cmeshxk->InsertMaterialObject(BCond1);
    matpoisson->SetForcingFunction(forcefunction,1);


    val2[0] = 0;
    TPZBndCondT<STATE> *BCond2 = matpoisson -> CreateBC(matpoisson , ERight, 1, val1, val2);
    cmeshxk->InsertMaterialObject(BCond2);

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