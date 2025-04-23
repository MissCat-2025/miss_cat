#include "MissCatApp.h"
#include "Moose.h"
#include "AppFactory.h"
#include "ModulesApp.h"
#include "MooseSyntax.h"

InputParameters
MissCatApp::validParams()
{
  InputParameters params = MooseApp::validParams();
  params.set<bool>("use_legacy_material_output") = false;
  params.set<bool>("use_legacy_initial_residual_evaluation_behavior") = false;
  return params;
}

MissCatApp::MissCatApp(InputParameters parameters) : MooseApp(parameters)
{
  MissCatApp::registerAll(_factory, _action_factory, _syntax);
}

MissCatApp::~MissCatApp() {}

void
MissCatApp::registerAll(Factory & f, ActionFactory & af, Syntax & syntax)
{
  ModulesApp::registerAllObjects<MissCatApp>(f, af, syntax);
  Registry::registerObjectsTo(f, {"MissCatApp"});
  Registry::registerActionsTo(af, {"MissCatApp"});

  /* register custom execute flags, action syntax, etc. here */
}

void
MissCatApp::registerApps()
{
  registerApp(MissCatApp);
}

/***************************************************************************************************
 *********************** Dynamic Library Entry Points - DO NOT MODIFY ******************************
 **************************************************************************************************/
extern "C" void
MissCatApp__registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  MissCatApp::registerAll(f, af, s);
}
extern "C" void
MissCatApp__registerApps()
{
  MissCatApp::registerApps();
}
