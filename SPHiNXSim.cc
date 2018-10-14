/*
 * ==========================================================================================================
 *                                       SPHiNXSim - SPHiNX Simulation
 * ==========================================================================================================
 * 
 * Versions -
 *      2013-10 - Maxime Chauvin
 * 
 * Command line -
 *      Program name - SPHiNXSim
 *      Arg[1] - Macro file name
 *      Arg[2] - Output file name with no extension
 *
 *      ./SPHiNXSim run.mac output
 *
 * Output -
 *
 *     output.root
 *       G4data (Tree)
 *         Primary   (Branch)
 *           int Id              (Leaf)
 *           float Energy        (Leaf)
 *           float Px            (Leaf)
 *           float Py            (Leaf)
 *           float Pz            (Leaf)
 *           float Dx            (Leaf)
 *           float Dy            (Leaf)
 *           float Dz            (Leaf)
 *         Slow      (Branch)
 *           int nhits           (Leaf)
 *           float EdepTot       (Leaf)
 *           float Edep[nhits]   (Leaf)
 *           float Px[nhits]     (Leaf)
 *           float Py[nhits]     (Leaf)
 *           float Pz[nhits]     (Leaf)
 *           float UnitId[nhits] (Leaf)
 *           int nUnits          (Leaf)
 *           float EdepTot[61]   (Leaf)
 *         Fast      (Branch)
 *           ...
 *         HighZ     (Branch)
 *            ...
 * 
 *     output.txt
 *       8.40211 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 
 *       89.4844 0 0 0 0 0 0 349.865 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 
 *       5.04741 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 
 *       ...
 *
 * ==========================================================================================================
 */

#include "G4RunManager.hh"
#include "G4UImanager.hh"

#include "DetectorConstruction.hh"
#include "PhysicsList.hh"
#include "PrimaryGeneratorAction.hh"
#include "RunAction.hh"
#include "EventAction.hh"

#ifdef G4VIS_USE
 #include "G4VisExecutive.hh"
 #include "G4UIExecutive.hh"
#endif

int main(int argc,char** argv) {

  //Construct the default run manager
  G4RunManager* runManager = new G4RunManager;

  //set mandatory initialization classes
  DetectorConstruction* detector = new DetectorConstruction();
  runManager->SetUserInitialization(detector);
    
  PhysicsList* physics = new PhysicsList();
  runManager->SetUserInitialization(physics);
    
  //set user action classes
  PrimaryGeneratorAction* gen_action = new PrimaryGeneratorAction();
  runManager->SetUserAction(gen_action);
  
  RunAction* run_action = new RunAction();
  runManager->SetUserAction(run_action);

  EventAction* event_action = new EventAction();
  runManager->SetUserAction(event_action);

  //Initialize Geant4 kernel
  runManager->Initialize();

  //get the pointer to the User Interface manager 
  G4UImanager* UIm = G4UImanager::GetUIpointer();  
  
  if (argc==1) {// interactive mode
  #ifdef G4VIS_USE
    event_action->CreateROOTFile("vis.root");
    G4VisManager* visManager = new G4VisExecutive;
    visManager->Initialize();
    G4UIExecutive* UIi = new G4UIExecutive(argc, argv);
    UIm->ApplyCommand("/control/execute ../vis0.mac");
    //UIm->ApplyCommand("/control/execute autoOffaxis.mac");
    UIi->SessionStart();
    delete UIi;
    delete visManager;
  #endif
  }
  else if (argc==3) {// batch mode
	G4String ext_ascii = ".txt";
    event_action->CreateASCIIFile(argv[2] + ext_ascii);
    G4String ext_root = ".root";
    event_action->CreateROOTFile(argv[2] + ext_root);
    G4String fileName = argv[1];
  	G4String command = "/control/execute ";
    UIm->ApplyCommand(command + fileName);
  }
  else {
	G4cout << "Not the right number of arguments. Usage:" << G4endl;
	G4cout << " -interactive mode ./SPHiNXSim" << G4endl;
	G4cout << " -batch mode ./SPHiNXSim macro.mac output" << G4endl;
	exit(1);
  }

  //job termination
  delete runManager;
  return 0;
}
