//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html
#include "MissCatTestApp.h"
#include "MissCatApp.h"
#include "Moose.h"
#include "AppFactory.h"
#include "MooseSyntax.h"

InputParameters
MissCatTestApp::validParams()
{
  InputParameters params = MissCatApp::validParams();
  params.set<bool>("use_legacy_material_output") = false;
  params.set<bool>("use_legacy_initial_residual_evaluation_behavior") = false;
  return params;
}

MissCatTestApp::MissCatTestApp(InputParameters parameters) : MooseApp(parameters)
{
  MissCatTestApp::registerAll(
      _factory, _action_factory, _syntax, getParam<bool>("allow_test_objects"));
}

MissCatTestApp::~MissCatTestApp() {}

void
MissCatTestApp::registerAll(Factory & f, ActionFactory & af, Syntax & s, bool use_test_objs)
{
  MissCatApp::registerAll(f, af, s);
  if (use_test_objs)
  {
    Registry::registerObjectsTo(f, {"MissCatTestApp"});
    Registry::registerActionsTo(af, {"MissCatTestApp"});
  }
}

void
MissCatTestApp::registerApps()
{
  registerApp(MissCatApp);
  registerApp(MissCatTestApp);
}

/***************************************************************************************************
 *********************** Dynamic Library Entry Points - DO NOT MODIFY ******************************
 **************************************************************************************************/
// External entry point for dynamic application loading
extern "C" void
MissCatTestApp__registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  MissCatTestApp::registerAll(f, af, s);
}
extern "C" void
MissCatTestApp__registerApps()
{
  MissCatTestApp::registerApps();
}
