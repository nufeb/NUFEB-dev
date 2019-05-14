/////////////////////////////////////////////////////////////////////////////
// Name:        minimal.cpp
// Purpose:     Minimal wxWidgets sample
// Author:      Julian Smart
// Modified by:
// Created:     04/01/98
// RCS-ID:      $Id$
// Copyright:   (c) Julian Smart
// Licence:     wxWindows licence
/////////////////////////////////////////////////////////////////////////////

// ============================================================================
// declarations
// ============================================================================

// ----------------------------------------------------------------------------
// headers
// ----------------------------------------------------------------------------

// For compilers that support precompilation, includes "wx/wx.h".
#include "wx/wxprec.h"

#ifdef __BORLANDC__
    #pragma hdrstop
#endif

// for all others, include the necessary headers (this file is usually all you
// need because it includes almost all "standard" wxWidgets headers)
#ifndef WX_PRECOMP
    #include "wx/wx.h"
#endif

#include <wx/grid.h>
#include <iostream>
#include <fstream>

// ----------------------------------------------------------------------------
// resources
// ----------------------------------------------------------------------------

// the application icon (under Windows and OS/2 it is in resources and even
// though we could still include the XPM here it would be unused)
#ifndef wxHAS_IMAGES_IN_RESOURCES
    #include "sample.xpm"
#endif

// ----------------------------------------------------------------------------
// private classes
// ----------------------------------------------------------------------------

// Define a new application type, each program should derive a class from wxApp
class MyApp : public wxApp
{
public:
    // override base class virtuals
    // ----------------------------

    // this one is called on application startup and is a good place for the app
    // initialization (doing it here and not in the ctor allows to have an error
    // return: if OnInit() returns false, the application terminates)
    virtual bool OnInit();
};

// Define a new frame type: this is going to be our main frame
class MyFrame : public wxFrame
{
public:
    // ctor(s)
    MyFrame(const wxString& title);

    // event handlers (these functions should _not_ be virtual)
    void OnQuit(wxCommandEvent& event);
    void OnAbout(wxCommandEvent& event);
    void InputBrowse(wxCommandEvent& event);
    void OutputBrowse(wxCommandEvent& event);
    void Generate(wxCommandEvent& event);
    void Run(wxCommandEvent& event);
    void ClickGrowth(wxCommandEvent& event);
    void ClickDeath(wxCommandEvent& event);
    void ClickDivision(wxCommandEvent& event);
    void ClickPhysics(wxCommandEvent& event);
    void ClickGravity(wxCommandEvent& event);

private:
    // any class wishing to process wxWidgets events must use this macro
    wxDECLARE_EVENT_TABLE();

    wxRadioButton *periodicX;
    wxRadioButton *fixedX;
    wxRadioButton *periodicY;
    wxRadioButton *fixedY;
    wxRadioButton *periodicZ;
    wxRadioButton *fixedZ;

    wxTextCtrl *neighbor;
    wxTextCtrl *timestep;

    wxTextCtrl *inputFile;
    wxTextCtrl *outputFile;

    wxCheckBox *particlePhysics;

    wxTextCtrl *kn;
    wxTextCtrl *kt;
    wxTextCtrl *gamman;
    wxTextCtrl *gammat;
    wxTextCtrl *mu;
    wxCheckBox *flag;

    wxCheckBox *gravity;

    wxTextCtrl *magnitude;
    wxComboBox *direction;

    wxCheckBox *hetGrowth;
    wxCheckBox *hetDeath;
    wxCheckBox *hetDivision;
    wxCheckBox *hetEPSExcretion;

    wxCheckBox *aobGrowth;
    wxCheckBox *aobDeath;
    wxCheckBox *aobDivision;

    wxCheckBox *nobGrowth;
    wxCheckBox *nobDeath;
    wxCheckBox *nobDivision;

    wxCheckBox *allGrowth;
    wxCheckBox *allDeath;
    wxCheckBox *allDivision;

    wxTextCtrl *hetGrowthRate;
    wxTextCtrl *hetDecayRate;
    wxTextCtrl *hetYieldRate;

    wxTextCtrl *aobGrowthRate;
    wxTextCtrl *aobDecayRate;
    wxTextCtrl *aobYieldRate;

    wxTextCtrl *nobGrowthRate;
    wxTextCtrl *nobDecayRate;
    wxTextCtrl *nobYieldRate;

    wxTextCtrl *epsDecayRate;
    wxTextCtrl *epsYieldRate;

    wxTextCtrl *hetOxygen;
    wxTextCtrl *aobOxygen;
    wxTextCtrl *nobOxygen;

    wxTextCtrl *hetNitrite;
    wxTextCtrl *nobNitrite;

    wxTextCtrl *hetCarbon;
    wxTextCtrl *hetNitrate;
    wxTextCtrl *aobAmmonia;
    wxTextCtrl *hetReduction;
    wxTextCtrl *deathFactor;

    wxTextCtrl *epsDensity;
    wxTextCtrl *epsRatio;
};

// ----------------------------------------------------------------------------
// constants
// ----------------------------------------------------------------------------

// IDs for the controls and the menu commands
enum
{
    // menu items
    Minimal_Quit = wxID_EXIT,

    // it is important for the id corresponding to the "About" command to have
    // this standard value as otherwise it won't be handled properly under Mac
    // (where it is special and put into the "Apple" menu)
    Minimal_About = wxID_ABOUT,
    BUTTON_InputBrowse = wxID_HIGHEST + 1,
    BUTTON_OutputBrowse = wxID_HIGHEST + 2,
    BUTTON_Generate = wxID_HIGHEST + 3,
    BUTTON_Run = wxID_HIGHEST + 4,
    CHECKBOX_Growth = wxID_HIGHEST + 5,
    CHECKBOX_Death = wxID_HIGHEST + 6,
    CHECKBOX_Division = wxID_HIGHEST + 7,
    CHECKBOX_Physics = wxID_HIGHEST + 8,
    CHECKBOX_Gravity = wxID_HIGHEST + 9
};

// ----------------------------------------------------------------------------
// event tables and other macros for wxWidgets
// ----------------------------------------------------------------------------

// the event tables connect the wxWidgets events with the functions (event
// handlers) which process them. It can be also done at run-time, but for the
// simple menu events like this the static method is much simpler.
wxBEGIN_EVENT_TABLE(MyFrame, wxFrame)
    EVT_MENU(Minimal_Quit,  MyFrame::OnQuit)
    EVT_MENU(Minimal_About, MyFrame::OnAbout)
    EVT_BUTTON ( BUTTON_InputBrowse, MyFrame::InputBrowse )
    EVT_BUTTON ( BUTTON_OutputBrowse, MyFrame::OutputBrowse )
    EVT_BUTTON ( BUTTON_Generate, MyFrame::Generate )
    EVT_BUTTON ( BUTTON_Run, MyFrame::Run )
    EVT_CHECKBOX ( CHECKBOX_Growth, MyFrame::ClickGrowth )
    EVT_CHECKBOX ( CHECKBOX_Death, MyFrame::ClickDeath )
    EVT_CHECKBOX ( CHECKBOX_Division, MyFrame::ClickDivision )
    EVT_CHECKBOX ( CHECKBOX_Physics, MyFrame::ClickPhysics )
    EVT_CHECKBOX ( CHECKBOX_Gravity, MyFrame::ClickGravity )
wxEND_EVENT_TABLE()

// Create a new application object: this macro will allow wxWidgets to create
// the application object during program execution (it's better than using a
// static object for many reasons) and also implements the accessor function
// wxGetApp() which will return the reference of the right type (i.e. MyApp and
// not wxApp)
IMPLEMENT_APP(MyApp)

// ============================================================================
// implementation
// ============================================================================

// ----------------------------------------------------------------------------
// the application class
// ----------------------------------------------------------------------------

// 'Main program' equivalent: the program execution "starts" here
bool MyApp::OnInit()
{
    // call the base class initialization method, currently it only parses a
    // few common command-line options but it could be do more in the future
    if ( !wxApp::OnInit() )
        return false;

    // create the main application window
    MyFrame *frame = new MyFrame("NUFEB Model");

    // and show it (the frames, unlike simple controls, are not shown when
    // created initially)
    frame->Show(true);

    // success: wxApp::OnRun() will be called which will enter the main message
    // loop and the application will run. If we returned false here, the
    // application would exit immediately.
    return true;
}

// ----------------------------------------------------------------------------
// main frame
// ----------------------------------------------------------------------------

// frame constructor
MyFrame::MyFrame(const wxString& title)
       : wxFrame(NULL, wxID_ANY, title, wxDefaultPosition, wxSize(1250,1250))
{
    // set the frame icon
    SetIcon(wxICON(sample));

#if wxUSE_MENUS
    // create a menu bar
    wxMenu *fileMenu = new wxMenu;

    // the "About" item should be in the help menu
    wxMenu *helpMenu = new wxMenu;
    helpMenu->Append(Minimal_About, "&About\tF1", "Show about dialog");

    fileMenu->Append(Minimal_Quit, "E&xit\tAlt-X", "Quit this program");

    // now append the freshly created menu to the menu bar...
    wxMenuBar *menuBar = new wxMenuBar();
    menuBar->Append(fileMenu, "&File");
    menuBar->Append(helpMenu, "&Help");

    // ... and attach this menu bar to the frame
    SetMenuBar(menuBar);
#endif // wxUSE_MENUS

#if wxUSE_STATUSBAR
    // create a status bar just for fun (by default with 1 pane only)
    CreateStatusBar(2);
    SetStatusText("Enter parameters");
#endif // wxUSE_STATUSBAR


    wxPanel *panel = new wxPanel(this, -1);

    wxBoxSizer *main = new wxBoxSizer(wxVERTICAL);

    wxBoxSizer *sizer = new wxBoxSizer(wxVERTICAL);


    sizer->Add(new wxStaticText(panel, -1, "All units should be specified in SI"));

    sizer->Add(new wxStaticText(panel, -1, "", wxDefaultPosition));

    wxBoxSizer *vbox1 = new wxBoxSizer(wxVERTICAL);
    wxBoxSizer *vbox2 = new wxBoxSizer(wxVERTICAL);
    wxBoxSizer *vbox3 = new wxBoxSizer(wxVERTICAL);
    wxBoxSizer *vbox4 = new wxBoxSizer(wxVERTICAL);

    wxBoxSizer *hboxBoundary = new wxBoxSizer(wxHORIZONTAL);


    periodicX = new wxRadioButton(panel, -1, "", wxDefaultPosition, wxDefaultSize, wxRB_GROUP);
    fixedX = new wxRadioButton(panel, -1, "");
    periodicY = new wxRadioButton(panel, -1, "", wxDefaultPosition, wxDefaultSize, wxRB_GROUP);
    fixedY = new wxRadioButton(panel, -1, "");
    periodicZ = new wxRadioButton(panel, -1, "", wxDefaultPosition, wxDefaultSize, wxRB_GROUP);
    fixedZ = new wxRadioButton(panel, -1, "");

    vbox1->Add(new wxStaticText(panel, -1, "Boundary"));
    vbox1->Add(new wxStaticText(panel, -1, "Periodic:", wxDefaultPosition, wxSize(-1, 26)));
    vbox1->Add(new wxStaticText(panel, -1, "Fixed:", wxDefaultPosition, wxSize(-1, 26)));

    vbox2->Add(new wxStaticText(panel, -1, "x"), 0, wxALIGN_CENTER);
    vbox2->Add(periodicX, 0, wxALIGN_CENTER);
    vbox2->Add(fixedX, 0, wxALIGN_CENTER);

    vbox3->Add(new wxStaticText(panel, -1, "y"), 0, wxALIGN_CENTER);
    vbox3->Add(periodicY, 0, wxALIGN_CENTER);
    vbox3->Add(fixedY, 0, wxALIGN_CENTER);

    vbox4->Add(new wxStaticText(panel, -1, "z"), 0, wxALIGN_CENTER);
    vbox4->Add(periodicZ, 0, wxALIGN_CENTER);
    vbox4->Add(fixedZ, 0, wxALIGN_CENTER);

    hboxBoundary->Add(vbox1);
    hboxBoundary->Add(vbox2);
    hboxBoundary->Add(vbox3);
    hboxBoundary->Add(vbox4);

    sizer->Add(hboxBoundary);

    sizer->Add(new wxStaticText(panel, -1, "", wxDefaultPosition));


    wxBoxSizer *hboxInput = new wxBoxSizer(wxHORIZONTAL);


    inputFile = new wxTextCtrl(panel, -1, "", wxDefaultPosition, wxSize(400,-1), wxTE_LEFT | wxTE_READONLY);

    wxButton *inputFileBrowse = new wxButton(panel, BUTTON_InputBrowse, _T("Browse"));

    hboxInput->Add(new wxStaticText(panel, -1, "Input File:"));
    hboxInput->Add(inputFile);
    hboxInput->Add(inputFileBrowse);

    sizer->Add(hboxInput);


    wxBoxSizer *hboxNeighbor = new wxBoxSizer(wxHORIZONTAL);


    neighbor = new wxTextCtrl(panel, -1, "1.0e-5");

    hboxNeighbor->Add(new wxStaticText(panel, -1, "Neighbor Skin Distance (m):"));
    hboxNeighbor->Add(neighbor);

    sizer->Add(hboxNeighbor);


    wxBoxSizer *hboxTimestep = new wxBoxSizer(wxHORIZONTAL);


    timestep = new wxTextCtrl(panel, -1, "1.0");

    hboxTimestep->Add(new wxStaticText(panel, -1, "Time Step (s):"));
    hboxTimestep->Add(timestep);

    sizer->Add(hboxTimestep);

    sizer->Add(new wxStaticText(panel, -1, "", wxDefaultPosition));

    particlePhysics = new wxCheckBox(panel, CHECKBOX_Physics, "Particle Physics");
    particlePhysics->SetValue(true);

    sizer->Add(particlePhysics);

    wxBoxSizer *hboxPhysics = new wxBoxSizer(wxHORIZONTAL);

    kn = new wxTextCtrl(panel, -1, "200000000");
    kt = new wxTextCtrl(panel, -1, "");
    gamman = new wxTextCtrl(panel, -1, "15000000");
    gammat = new wxTextCtrl(panel, -1, "");
    mu = new wxTextCtrl(panel, -1, "0.5");
    flag = new wxCheckBox(panel, -1, "Flag");
    flag->SetValue(true);

    hboxPhysics->Add(new wxStaticText(panel, -1, "k_n:"));
    hboxPhysics->Add(kn);
    hboxPhysics->Add(new wxStaticText(panel, -1, "k_t:"));
    hboxPhysics->Add(kt);
    hboxPhysics->Add(new wxStaticText(panel, -1, "gamma_n:"));
    hboxPhysics->Add(gamman);
    hboxPhysics->Add(new wxStaticText(panel, -1, "gamma_t:"));
    hboxPhysics->Add(gammat);
    hboxPhysics->Add(new wxStaticText(panel, -1, "mu:"));
    hboxPhysics->Add(mu);
    hboxPhysics->Add(flag);

    sizer->Add(hboxPhysics);

    gravity = new wxCheckBox(panel, CHECKBOX_Gravity, "Gravity");
    gravity->SetValue(true);

    sizer->Add(gravity);

    wxBoxSizer *hboxGravity = new wxBoxSizer(wxHORIZONTAL);

    wxString dir[6];

    magnitude = new wxTextCtrl(panel, -1, "9.8");
    direction = new wxComboBox(panel, -1, "down", wxDefaultPosition, wxDefaultSize, 6, dir);
    direction->SetString(0,"up");
    direction->SetString(1,"down");
    direction->SetString(2,"left");
    direction->SetString(3,"right");
    direction->SetString(4,"forward");
    direction->SetString(5,"backward");

    hboxGravity->Add(new wxStaticText(panel, -1, "magnitude (m/s^2):"));
    hboxGravity->Add(magnitude);
    hboxGravity->Add(new wxStaticText(panel, -1, "direction:"));
    hboxGravity->Add(direction);

    sizer->Add(hboxGravity);

    sizer->Add(new wxStaticText(panel, -1, "", wxDefaultPosition));

    wxBoxSizer *vbox12 = new wxBoxSizer(wxVERTICAL);
    wxBoxSizer *vbox13 = new wxBoxSizer(wxVERTICAL);
    wxBoxSizer *vbox14 = new wxBoxSizer(wxVERTICAL);
    wxBoxSizer *vbox15 = new wxBoxSizer(wxVERTICAL);
    wxBoxSizer *vbox16 = new wxBoxSizer(wxVERTICAL);

    wxBoxSizer *hboxFixes = new wxBoxSizer(wxHORIZONTAL);


    vbox12->Add(new wxStaticText(panel, -1, ""));
    vbox12->Add(new wxStaticText(panel, -1, "Growth:", wxDefaultPosition, wxSize(-1, 26)));
    vbox12->Add(new wxStaticText(panel, -1, "Death:", wxDefaultPosition, wxSize(-1, 26)));
    vbox12->Add(new wxStaticText(panel, -1, "Division:", wxDefaultPosition, wxSize(-1, 26)));
    vbox12->Add(new wxStaticText(panel, -1, "EPS Excretion:", wxDefaultPosition, wxSize(-1, 26)));


    hetGrowth = new wxCheckBox(panel, -1, "");
    hetDeath = new wxCheckBox(panel, -1, "");
    hetDivision = new wxCheckBox(panel, -1, "");
    hetEPSExcretion = new wxCheckBox(panel, -1, "");

    vbox13->Add(new wxStaticText(panel, -1, "HET"), 0, wxALIGN_CENTER);
    vbox13->Add(hetGrowth, 0, wxALIGN_CENTER);
    vbox13->Add(hetDeath, 0, wxALIGN_CENTER);
    vbox13->Add(hetDivision, 0, wxALIGN_CENTER);
    vbox13->Add(hetEPSExcretion, 0, wxALIGN_CENTER);

    aobGrowth = new wxCheckBox(panel, -1, "");
    aobDeath = new wxCheckBox(panel, -1, "");
    aobDivision = new wxCheckBox(panel, -1, "");

    vbox14->Add(new wxStaticText(panel, -1, "AOB"), 0, wxALIGN_CENTER);
    vbox14->Add(aobGrowth, 0, wxALIGN_CENTER);
    vbox14->Add(aobDeath, 0, wxALIGN_CENTER);
    vbox14->Add(aobDivision, 0, wxALIGN_CENTER);
    vbox14->Add(new wxStaticText(panel, -1, ""), 0, wxALIGN_CENTER);

    nobGrowth = new wxCheckBox(panel, -1, "");
    nobDeath = new wxCheckBox(panel, -1, "");
    nobDivision = new wxCheckBox(panel, -1, "");

    vbox15->Add(new wxStaticText(panel, -1, "NOB"), 0, wxALIGN_CENTER);
    vbox15->Add(nobGrowth, 0, wxALIGN_CENTER);
    vbox15->Add(nobDeath, 0, wxALIGN_CENTER);
    vbox15->Add(nobDivision, 0, wxALIGN_CENTER);
    vbox15->Add(new wxStaticText(panel, -1, ""), 0, wxALIGN_CENTER);

    allGrowth = new wxCheckBox(panel, CHECKBOX_Growth, "");
    allDeath = new wxCheckBox(panel, CHECKBOX_Death, "");
    allDivision = new wxCheckBox(panel, CHECKBOX_Division, "");

    allGrowth->SetValue(true);
    hetGrowth->Enable(false);
	aobGrowth->Enable(false);
	nobGrowth->Enable(false);
	allDeath->SetValue(true);
    hetDeath->Enable(false);
	aobDeath->Enable(false);
	nobDeath->Enable(false);
	allDivision->SetValue(true);
    hetDivision->Enable(false);
	aobDivision->Enable(false);
	nobDivision->Enable(false);
	hetEPSExcretion->SetValue(true);

    vbox16->Add(new wxStaticText(panel, -1, "All"), 0, wxALIGN_CENTER);
    vbox16->Add(allGrowth, 0, wxALIGN_CENTER);
    vbox16->Add(allDeath, 0, wxALIGN_CENTER);
    vbox16->Add(allDivision, 0, wxALIGN_CENTER);
    vbox16->Add(new wxStaticText(panel, -1, ""), 0, wxALIGN_CENTER);

    hboxFixes->Add(vbox12);
    hboxFixes->Add(vbox13);
    hboxFixes->Add(new wxStaticText(panel, -1, "  "));
    hboxFixes->Add(vbox14);
    hboxFixes->Add(new wxStaticText(panel, -1, "  "));
    hboxFixes->Add(vbox15);
    hboxFixes->Add(new wxStaticText(panel, -1, "  "));
    hboxFixes->Add(vbox16);

    sizer->Add(hboxFixes);

    sizer->Add(new wxStaticText(panel, -1, "", wxDefaultPosition));


    wxBoxSizer *vbox5 = new wxBoxSizer(wxVERTICAL);
    wxBoxSizer *vbox6 = new wxBoxSizer(wxVERTICAL);
    wxBoxSizer *vbox7 = new wxBoxSizer(wxVERTICAL);
    wxBoxSizer *vbox8 = new wxBoxSizer(wxVERTICAL);
    wxBoxSizer *vbox9 = new wxBoxSizer(wxVERTICAL);

    wxBoxSizer *hboxRates = new wxBoxSizer(wxHORIZONTAL);

    vbox5->Add(new wxStaticText(panel, -1, ""));
    vbox5->Add(new wxStaticText(panel, -1, "Maximum Growth Rate (1/s):", wxDefaultPosition, wxSize(-1, 34)));
    vbox5->Add(new wxStaticText(panel, -1, "Decay Rate (1/s):", wxDefaultPosition, wxSize(-1, 34)));
    vbox5->Add(new wxStaticText(panel, -1, "Yield Coefficient (-):", wxDefaultPosition, wxSize(-1, 34)));


    hetGrowthRate = new wxTextCtrl(panel, -1, "6.944e-5");
    hetDecayRate = new wxTextCtrl(panel, -1, "4.6296e-6");
    hetYieldRate = new wxTextCtrl(panel, -1, "0.61");

    vbox6->Add(new wxStaticText(panel, -1, "HET"), 0, wxALIGN_CENTER);
    vbox6->Add(hetGrowthRate, 0, wxALIGN_CENTER);
    vbox6->Add(hetDecayRate, 0, wxALIGN_CENTER);
    vbox6->Add(hetYieldRate, 0, wxALIGN_CENTER);


    aobGrowthRate = new wxTextCtrl(panel, -1, "3.4722e-5");
    aobDecayRate = new wxTextCtrl(panel, -1, "1.2731e-6");
    aobYieldRate = new wxTextCtrl(panel, -1, "0.33");


    vbox7->Add(new wxStaticText(panel, -1, "AOB"), 0, wxALIGN_CENTER);
    vbox7->Add(aobGrowthRate, 0, wxALIGN_CENTER);
    vbox7->Add(aobDecayRate, 0, wxALIGN_CENTER);
    vbox7->Add(aobYieldRate, 0, wxALIGN_CENTER);


    nobGrowthRate = new wxTextCtrl(panel, -1, "3.4722e-5");
    nobDecayRate = new wxTextCtrl(panel, -1, "1.2731e-6");
    nobYieldRate = new wxTextCtrl(panel, -1, "0.083");


    vbox8->Add(new wxStaticText(panel, -1, "NOB"), 0, wxALIGN_CENTER);
    vbox8->Add(nobGrowthRate, 0, wxALIGN_CENTER);
    vbox8->Add(nobDecayRate, 0, wxALIGN_CENTER);
    vbox8->Add(nobYieldRate, 0, wxALIGN_CENTER);


    epsDecayRate = new wxTextCtrl(panel, -1, "1.9676e-6");
    epsYieldRate = new wxTextCtrl(panel, -1, "0.18");


    vbox9->Add(new wxStaticText(panel, -1, "EPS"), 0, wxALIGN_CENTER);
    vbox9->Add(new wxStaticText(panel, -1, "", wxDefaultPosition, wxSize(-1, 34)), 0, wxALIGN_CENTER);
    vbox9->Add(epsDecayRate, 0, wxALIGN_CENTER);
    vbox9->Add(epsYieldRate, 0, wxALIGN_CENTER);

    hboxRates->Add(vbox5);
    hboxRates->Add(vbox6);
    hboxRates->Add(vbox7);
    hboxRates->Add(vbox8);
    hboxRates->Add(vbox9);


    sizer->Add(hboxRates);

    sizer->Add(new wxStaticText(panel, -1, "", wxDefaultPosition));


    wxBoxSizer *vbox10 = new wxBoxSizer(wxVERTICAL);
    wxBoxSizer *vbox11 = new wxBoxSizer(wxVERTICAL);

    wxBoxSizer *hboxAffinities = new wxBoxSizer(wxHORIZONTAL);

    vbox10->Add(new wxStaticText(panel, -1, "Dissolved Oxygen Affinity (kg/m^3):", wxDefaultPosition, wxSize(-1, 34)));
    vbox10->Add(new wxStaticText(panel, -1, "Nitrite Affinity (kg/m^3):", wxDefaultPosition, wxSize(-1, 34)));
    vbox10->Add(new wxStaticText(panel, -1, "HET Carbon Source Affinity (kg/m^3):", wxDefaultPosition, wxSize(-1, 34)));
    vbox10->Add(new wxStaticText(panel, -1, "HET Nitrate Affinity (kg/m^3):", wxDefaultPosition, wxSize(-1, 34)));
    vbox10->Add(new wxStaticText(panel, -1, "AOB Ammonia Affinity (kg/m^3):", wxDefaultPosition, wxSize(-1, 34)));
    vbox10->Add(new wxStaticText(panel, -1, "HET Reduction Factor (-):", wxDefaultPosition, wxSize(-1, 34)));
    vbox10->Add(new wxStaticText(panel, -1, "Death Factor (-):", wxDefaultPosition, wxSize(-1, 34)));


    wxBoxSizer *hboxOxygen = new wxBoxSizer(wxHORIZONTAL);

    hetOxygen = new wxTextCtrl(panel, -1, "8.1e-4");
    aobOxygen = new wxTextCtrl(panel, -1, "5e-4");
    nobOxygen = new wxTextCtrl(panel, -1, "6.8e-4");

    hboxOxygen->Add(new wxStaticText(panel, -1, "HET:"));
    hboxOxygen->Add(hetOxygen);
    hboxOxygen->Add(new wxStaticText(panel, -1, "AOB:"));
    hboxOxygen->Add(aobOxygen);
    hboxOxygen->Add(new wxStaticText(panel, -1, "NOB:"));
    hboxOxygen->Add(nobOxygen);


    wxBoxSizer *hboxNitrite = new wxBoxSizer(wxHORIZONTAL);

    hetNitrite = new wxTextCtrl(panel, -1, "3e-4");
    nobNitrite = new wxTextCtrl(panel, -1, "1.3e-3");

    hboxNitrite->Add(new wxStaticText(panel, -1, "HET:"));
    hboxNitrite->Add(hetNitrite);
    hboxNitrite->Add(new wxStaticText(panel, -1, "NOB:"));
    hboxNitrite->Add(nobNitrite);

    hetCarbon = new wxTextCtrl(panel, -1, "1e-2");
    hetNitrate = new wxTextCtrl(panel, -1, "3e-4");
    aobAmmonia = new wxTextCtrl(panel, -1, "1e-3");
    hetReduction = new wxTextCtrl(panel, -1, "0.6");
    deathFactor = new wxTextCtrl(panel, -1, "2.0");


    vbox11->Add(hboxOxygen);
    vbox11->Add(hboxNitrite);
    vbox11->Add(hetCarbon);
    vbox11->Add(hetNitrate);
    vbox11->Add(aobAmmonia);
    vbox11->Add(hetReduction);
    vbox11->Add(deathFactor);

    hboxAffinities->Add(vbox10);
    hboxAffinities->Add(vbox11);

    sizer->Add(hboxAffinities);


    wxBoxSizer *hboxEPS = new wxBoxSizer(wxHORIZONTAL);

    epsDensity = new wxTextCtrl(panel, -1, "30");
    epsRatio = new wxTextCtrl(panel, -1, "1.25");

    hboxEPS->Add(new wxStaticText(panel, -1, "EPS:"));
    hboxEPS->Add(new wxStaticText(panel, -1, "Density (kg/m^3):"));
    hboxEPS->Add(epsDensity);
    hboxEPS->Add(new wxStaticText(panel, -1, "Extraction Ratio (-):"));
    hboxEPS->Add(epsRatio);


    sizer->Add(hboxEPS);

    sizer->Add(new wxStaticText(panel, -1, "", wxDefaultPosition));

    wxBoxSizer *hboxOutput = new wxBoxSizer(wxHORIZONTAL);

    outputFile = new wxTextCtrl(panel, -1, "", wxDefaultPosition, wxSize(400,-1), wxTE_LEFT | wxTE_READONLY);

    wxButton *outputFileBrowse = new wxButton(panel, BUTTON_OutputBrowse, _T("Browse"));

    hboxOutput->Add(new wxStaticText(panel, -1, "Output File:"));
    hboxOutput->Add(outputFile);
    hboxOutput->Add(outputFileBrowse);

    // wxListView *outputSpecies = new wxListView(panel, -1);

    // hboxOutput->Add(outputSpecies);

    wxString str[20];

    wxComboBox *outputOptions = new wxComboBox(panel, -1, "id", wxDefaultPosition, wxDefaultSize, 20, str);

    outputOptions->SetString(0,"id");
    outputOptions->SetString(1,"type");
    outputOptions->SetString(2,"mass");
    outputOptions->SetString(3,"x");
    outputOptions->SetString(4,"y");
    outputOptions->SetString(5,"z");
    outputOptions->SetString(6,"vx");
    outputOptions->SetString(7,"vy");
    outputOptions->SetString(8,"vz");
    outputOptions->SetString(9,"fx");
    outputOptions->SetString(10,"fy");
    outputOptions->SetString(11,"fz");
    outputOptions->SetString(12,"radius");
    outputOptions->SetString(13,"diameter");
    outputOptions->SetString(14,"outer diameter");
    outputOptions->SetString(15,"substrate");
    outputOptions->SetString(16,"O2");
    outputOptions->SetString(17,"NO2");
    outputOptions->SetString(18,"NO3");
    outputOptions->SetString(19,"NH4");

    hboxOutput->Add(outputOptions);

    sizer->Add(hboxOutput);

    sizer->Add(new wxStaticText(panel, -1, "", wxDefaultPosition));

    wxBoxSizer *hboxExecute = new wxBoxSizer(wxHORIZONTAL);

    wxButton *generateButton = new wxButton(panel, BUTTON_Generate, _T("Generate Input Script"));
    wxButton *runBrowse = new wxButton(panel, BUTTON_Run, _T("Run Model"));

    hboxExecute->Add(generateButton);
    hboxExecute->Add(runBrowse);

    sizer->Add(hboxExecute);

    main->Add(sizer, 1, wxEXPAND | wxALL, 20);

    panel->SetSizer(main);

    Centre();







//     wxGrid *grid = new wxGrid( this,
//                     -1,
//                     wxPoint( 0, 240 ),
//                     wxSize( 400, 300 ) );

//     grid->CreateGrid( 100, 10 );
// // We can set the sizes of individual rows and columns
// // in pixels
// grid->SetRowSize( 0, 60 );
// grid->SetColSize( 0, 120 );
// // And set grid cell contents as strings
// grid->SetCellValue( 0, 0, "wxGrid is good" );
// // We can specify that some cells are read->only
// grid->SetCellValue( 0, 3, "This is read->only" );
// grid->SetReadOnly( 0, 3 );
// // Colours can be specified for grid cell contents
// grid->SetCellValue(3, 3, "green on grey");
// grid->SetCellTextColour(3, 3, *wxGREEN);
// grid->SetCellBackgroundColour(3, 3, *wxLIGHT_GREY);
// // We can specify the some cells will store numeric
// // values rather than strings. Here we set grid column 5
// // to hold floating point values displayed with width of 6
// // and precision of 2
// grid->SetColFormatFloat(5, 6, 2);
// grid->SetCellValue(0, 6, "3.1415");





}


// event handlers

void MyFrame::OnQuit(wxCommandEvent& WXUNUSED(event))
{
	delete periodicX;
    delete fixedX;
    delete periodicY;
    delete fixedY;
    delete periodicZ;
    delete fixedZ;

    delete neighbor;
    delete timestep;

	delete inputFile;
	delete outputFile;

	delete particlePhysics;

    delete kn;
    delete kt;
    delete gamman;
    delete gammat;
    delete mu;
    delete flag;

    delete gravity;

    delete magnitude;
    delete direction;

	delete hetGrowth;
    delete hetDeath;
    delete hetDivision;
    delete hetEPSExcretion;

    delete aobGrowth;
    delete aobDeath;
    delete aobDivision;

    delete nobGrowth;
    delete nobDeath;
    delete nobDivision;

    delete allGrowth;
    delete allDeath;
    delete allDivision;

    delete hetGrowthRate;
    delete hetDecayRate;
    delete hetYieldRate;

    delete aobGrowthRate;
    delete aobDecayRate;
    delete aobYieldRate;

    delete nobGrowthRate;
    delete nobDecayRate;
    delete nobYieldRate;

    delete epsDecayRate;
    delete epsYieldRate;

    delete hetOxygen;
    delete aobOxygen;
    delete nobOxygen;

    delete hetNitrite;
    delete nobNitrite;

    delete hetCarbon;
    delete hetNitrate;
    delete aobAmmonia;
    delete hetReduction;
    delete deathFactor;

    delete epsDensity;
    delete epsRatio;

    // true is to force the frame to close
    Close(true);
}

void MyFrame::OnAbout(wxCommandEvent& WXUNUSED(event))
{
    wxMessageBox(wxString::Format
                 (
                    "Welcome to %s!\n"
                    "\n"
                    "This is the minimal wxWidgets sample\n"
                    "running under %s.",
                    wxVERSION_STRING,
                    wxGetOsDescription()
                 ),
                 "About wxWidgets minimal sample",
                 wxOK | wxICON_INFORMATION,
                 this);
}

void MyFrame::InputBrowse(wxCommandEvent& WXUNUSED(event))
{
    wxFileDialog* OpenDialog = new wxFileDialog(
		this, "Choose a file to open", wxEmptyString, wxEmptyString, 
		"Input Files (*.in)|*.in|Text Files (*.txt)|*.txt",
		wxFD_OPEN, wxDefaultPosition);
 
	// Creates a "open file" dialog with 4 file types
	if (OpenDialog->ShowModal() == wxID_OK) // if the user click "Open" instead of "Cancel"
	{
		inputFile->SetValue(OpenDialog->GetPath());
	}
 
	// Clean up after ourselves
	OpenDialog->Destroy();
}

void MyFrame::OutputBrowse(wxCommandEvent& WXUNUSED(event))
{
    wxFileDialog* OpenDialog = new wxFileDialog(
		this, "Choose a file to save", wxEmptyString, wxEmptyString, 
		"Output Files (*.out)|*.out|Text Files (*.txt)|*.txt|Bubble MD (*.bubblemd)|*.bubblemd",
		wxFD_SAVE, wxDefaultPosition);
 
	// Creates a "open file" dialog with 4 file types
	if (OpenDialog->ShowModal() == wxID_OK) // if the user click "Open" instead of "Cancel"
	{
		outputFile->SetValue(OpenDialog->GetPath());
	}
 
	// Clean up after ourselves
	OpenDialog->Destroy();
}

void MyFrame::Generate(wxCommandEvent& WXUNUSED(event))
{
	wxFileDialog* OpenDialog = new wxFileDialog(
		this, "Choose a file to save", wxEmptyString, wxEmptyString, 
		"LAMMPS Input Scripts (*.lammps)|*.lammps",
		wxFD_SAVE, wxDefaultPosition);
 
	// Creates a "open file" dialog with 4 file types
	if (OpenDialog->ShowModal() == wxID_OK) // if the user click "Open" instead of "Cancel"
	{
		std::ofstream myfile (OpenDialog->GetPath());
		myfile << "# NUFEB simulation\n\n";
		myfile << "atom_style	bio\n";
		myfile << "atom_modify	map array sort 10000 1e-5\n";
		myfile << "boundary	";
		if (periodicX->GetValue()) {
			myfile << "pp	";
		}
		else {
			myfile << "ff	";
		}
		if (periodicY->GetValue()) {
			myfile << "pp	";
		}
		else {
			myfile << "ff	";
		}
		if (periodicZ->GetValue()) {
			myfile << "pp\n";
		}
		else {
			myfile << "ff\n";
		}
		myfile << "newton		off\n\n";
		myfile << "communicate single vel yes\n";
		if (inputFile->GetValue().compare("") != 0) {
			myfile << "read_data 	";
			myfile << inputFile->GetValue();
			myfile << "\n\n";
		}
		else {
			myfile << "\n";
		}
		myfile << "group HET type 1\n";
		myfile << "group AOB type 2\n";
		myfile << "group NOB type 3\n";
		myfile << "group EPS type 4\n";
		myfile << "group inert type 5\n\n";
		myfile << "neighbor	";
		myfile << neighbor->GetValue();
		myfile << " bin\n";
		myfile << "neigh_modify	delay 0\n\n";


		if (particlePhysics->IsChecked()) {
			myfile << "pair_style  gran/hooke/history ";
			if (kn->GetValue().compare("") == 0) {
				myfile << "NULL ";
			}
			else {
				myfile << kn->GetValue();
				myfile << " ";
			}
			if (kt->GetValue().compare("") == 0) {
				myfile << "NULL ";
			}
			else {
				myfile << kt->GetValue();
				myfile << " ";
			}
			if (gamman->GetValue().compare("") == 0) {
				myfile << "NULL ";
			}
			else {
				myfile << gamman->GetValue();
				myfile << " ";
			}
			if (gammat->GetValue().compare("") == 0) {
				myfile << "NULL ";
			}
			else {
				myfile << gammat->GetValue();
				myfile << " ";
			}
			myfile << mu->GetValue();
			myfile << " ";
			if (flag->IsChecked()) {
				myfile << "1";
			}
			else {
				myfile << "0";
			}
			myfile << "\npair_coeff  * *\n";
		}
		myfile << "timestep	";
		myfile << timestep->GetValue();
		myfile << "\n\n";
		myfile << "velocity    all set 0.0 0.0 0.0 units box\n\n";
		myfile << "fix		1 all nve/sphere\n";
		if (particlePhysics->IsChecked() && gravity->IsChecked()) {
			myfile << "fix		2 all gravity ";
			myfile << magnitude->GetValue();
			myfile << " vector ";
			if (direction->GetValue().compare("down") == 0) {
				myfile << "0 -1 0\n";
			}
			else if (direction->GetValue().compare("up") == 0) {
				myfile << "0 1 0\n";
			}
			else if (direction->GetValue().compare("left") == 0) {
				myfile << "-1 0 0\n";				
			}
			else if (direction->GetValue().compare("right") == 0) {
				myfile << "1 0 0\n";				
			}
			else if (direction->GetValue().compare("forward") == 0) {
				myfile << "0 0 1\n";
			}
			else {
				myfile << "0 0 -1\n";				
			}
		}
		myfile << "\n";	
		myfile << "variable KsHET equal ";
		myfile << hetCarbon->GetValue();
		myfile << "\n";
		myfile << "variable Ko2HET equal ";
		myfile << hetOxygen->GetValue();
		myfile << "\n";
		myfile << "variable Kno2HET equal ";
		myfile << hetNitrite->GetValue();
		myfile << "\n";
		myfile << "variable Kno3HET equal ";
		myfile << hetNitrate->GetValue();
		myfile << "\n";
		myfile << "variable Knh4AOB equal ";
		myfile << aobAmmonia->GetValue();
		myfile << "\n";
		myfile << "variable Ko2AOB equal ";
		myfile << aobOxygen->GetValue();
		myfile << "\n";
		myfile << "variable Kno2NOB equal ";
		myfile << nobNitrite->GetValue();
		myfile << "\n";
		myfile << "variable Ko2NOB equal ";
		myfile << nobOxygen->GetValue();
		myfile << "\n";
		myfile << "variable MumHET equal ";
		myfile << hetGrowthRate->GetValue();
		myfile << "\n";
		myfile << "variable MumAOB equal ";
		myfile << aobGrowthRate->GetValue();
		myfile << "\n";
		myfile << "variable MumNOB equal ";
		myfile << nobGrowthRate->GetValue();
		myfile << "\n";
		myfile << "variable etaHET equal ";
		myfile << hetReduction->GetValue();
		myfile << "\n";
		myfile << "variable bHET equal ";
		myfile << hetDecayRate->GetValue();
		myfile << "\n";
		myfile << "variable bAOB equal ";
		myfile << aobDecayRate->GetValue();
		myfile << "\n";
		myfile << "variable bNOB equal ";
		myfile << nobDecayRate->GetValue();
		myfile << "\n";
		myfile << "variable bEPS equal ";
		myfile << epsDecayRate->GetValue();
		myfile << "\n";
		myfile << "variable YEPS equal ";
		myfile << epsYieldRate->GetValue();
		myfile << "\n";
		myfile << "variable YHET equal ";
		myfile << hetYieldRate->GetValue();
		myfile << "\n";
		myfile << "variable EPSdens equal ";
		myfile << epsDensity->GetValue();
		myfile << "\n";
		myfile << "variable EPSratio equal ";
		myfile << epsRatio->GetValue();
		myfile << "\n";
		myfile << "variable factor equal ";
		myfile << deathFactor->GetValue();
		myfile << "\n\n";
		if (allGrowth->IsChecked()) {
			myfile << "fix g1 all nugrowth 1 v_KsHET v_Ko2HET v_Kno2HET v_Kno3HET v_Knh4AOB v_Ko2AOB v_Kno2NOB v_Ko2NOB v_MumHET v_MumAOB v_MumNOB v_etaHET v_bEPS v_YEPS v_YHET v_EPSdens\n";
		}
		else {
			if (hetGrowth->IsChecked()) {
				myfile << "fix g1 HET nugrowth 1 v_KsHET v_Ko2HET v_Kno2HET v_Kno3HET v_Knh4AOB v_Ko2AOB v_Kno2NOB v_Ko2NOB v_MumHET v_MumAOB v_MumNOB v_etaHET v_bEPS v_YEPS v_YHET v_EPSdens\n";
			}
			if (aobGrowth->IsChecked()) {
				myfile << "fix g2 AOB nugrowth 1 v_KsHET v_Ko2HET v_Kno2HET v_Kno3HET v_Knh4AOB v_Ko2AOB v_Kno2NOB v_Ko2NOB v_MumHET v_MumAOB v_MumNOB v_etaHET v_bEPS v_YEPS v_YHET v_EPSdens\n";
			}
			if (nobGrowth->IsChecked()) {
				myfile << "fix g3 NOB nugrowth 1 v_KsHET v_Ko2HET v_Kno2HET v_Kno3HET v_Knh4AOB v_Ko2AOB v_Kno2NOB v_Ko2NOB v_MumHET v_MumAOB v_MumNOB v_etaHET v_bEPS v_YEPS v_YHET v_EPSdens\n";
			}
		}
		if (allDeath->IsChecked()) {
			myfile << "fix dt1 HET death 1 v_bHET v_factor ";
			myfile << rand();
			myfile << "\n";
			myfile << "fix dt2 AOB death 1 v_bAOB v_factor ";
			myfile << rand();
			myfile << "\n";
			myfile << "fix dt3 NOB death 1 v_bNOB v_factor ";
			myfile << rand();
			myfile << "\n";
		}
		else {
			if (hetDeath->IsChecked()) {
				myfile << "fix dt1 HET death 1 v_bHET v_factor ";
				myfile << rand();
				myfile << "\n";
			}
			if (aobDeath->IsChecked()) {
				myfile << "fix dt2 AOB death 1 v_bAOB v_factor ";
				myfile << rand();
				myfile << "\n";
			}
			if (nobDeath->IsChecked()) {
				myfile << "fix dt3 NOB death 1 v_bNOB v_factor ";
				myfile << rand();
				myfile << "\n";
			}
		}
		if (allDivision->IsChecked()) {
			myfile << "fix d1 all divide 1 v_EPSdens 2.0 ";
			myfile << rand();
			myfile << "\n";
		}
		else {
			if (hetDivision->IsChecked()) {
				myfile << "fix d1 HET divide 1 v_EPSdens 2.0 ";
				myfile << rand();
				myfile << "\n";
			}
			if (aobDivision->IsChecked()) {
				myfile << "fix d2 AOB divide 1 v_EPSdens 2.0 ";
				myfile << rand();
				myfile << "\n";
			}
			if (nobDivision->IsChecked()) {
				myfile << "fix d3 NOB divide 1 v_EPSdens 2.0 ";
				myfile << rand();
				myfile << "\n";
			}
		}
		if (hetEPSExcretion->IsChecked()) {
			myfile << "fix e1 HET eps_extract 1 v_EPSratio v_EPSdens ";
			myfile << rand();
			myfile << "\n";
		}
		myfile << "\n";
		myfile << "dump		id all custom 2000 snapshot.bubblemd id type radius vx vy vz x y z outerradius outermass\n\n";
		myfile << "thermo_style    custom step atoms ke vol\n\n";
		myfile << "thermo		1\n";
		myfile << "thermo_modify	lost error\n\n";
		myfile << "run 172800\n";
  		myfile.close();
	}
 
	// Clean up after ourselves
	OpenDialog->Destroy();

}

void MyFrame::Run(wxCommandEvent& WXUNUSED(event))
{
	
}

void MyFrame::ClickGrowth(wxCommandEvent& event)
{
	if (allGrowth->IsChecked()) {
		hetGrowth->Enable(false);
		aobGrowth->Enable(false);
		nobGrowth->Enable(false);
	}
	else {
		hetGrowth->Enable(true);
		aobGrowth->Enable(true);
		nobGrowth->Enable(true);
	}
}

void MyFrame::ClickDeath(wxCommandEvent& event)
{
	if (allDeath->IsChecked()) {
		hetDeath->Enable(false);
		aobDeath->Enable(false);
		nobDeath->Enable(false);
	}
	else {
		hetDeath->Enable(true);
		aobDeath->Enable(true);
		nobDeath->Enable(true);
	}
}

void MyFrame::ClickDivision(wxCommandEvent& event)
{
	if (allDivision->IsChecked()) {
		hetDivision->Enable(false);
		aobDivision->Enable(false);
		nobDivision->Enable(false);
	}
	else {
		hetDivision->Enable(true);
		aobDivision->Enable(true);
		nobDivision->Enable(true);
	}
}

void MyFrame::ClickPhysics(wxCommandEvent& event)
{
	if (particlePhysics->IsChecked()) {
		kn->Enable(true);
    	kt->Enable(true);
    	gamman->Enable(true);
    	gammat->Enable(true);
    	mu->Enable(true);
    	flag->Enable(true);
    	gravity->Enable(true);
    	ClickGravity(event);
	}
	else {
		kn->Enable(false);
    	kt->Enable(false);
    	gamman->Enable(false);
    	gammat->Enable(false);
    	mu->Enable(false);
    	flag->Enable(false);
    	gravity->Enable(false);
    	magnitude->Enable(false);
    	direction->Enable(false);
	}
}

void MyFrame::ClickGravity(wxCommandEvent& event)
{
	if (gravity->IsChecked()) {
		magnitude->Enable(true);
    	direction->Enable(true);
	}
	else {
		magnitude->Enable(false);
    	direction->Enable(false);
	}
}
