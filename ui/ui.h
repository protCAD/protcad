# include <string>

// Computational Constants
# define MAX_SOLUTION_CONDITIONS 100
# define MAX_SOLUTES 8
# define MAX_SOLVENTS 8
# define CONC_SIZE 110
# define LINE_INPUT_SIZE 100
# define NO_VALUE "----"

using namespace std;

# ifndef UI_H
# define UI_H

# include <QWidget>

struct fromUIToZPRED
{string delimiter;
string* id;
string files;
string outputTitle;
string* pH;
string* proteinConc;
string* solvent;
string* solventConc;
string* solventConcType;
string* solute;
string* soluteConc;
string* soluteConcType;
string* temperature;
string apbs;
string multivalue;
string hydropro;
string msms;
string msmsAtmTypeNumbers;
string msmsPdbToXyzr;
string msmsPdbToXyzrn;
string pdb2pqr;
string zpred;
string maxThreads;
// Forced Solution Parameters (-1 if not specified)
string* fdensity;
string* fviscosity;
string* fdielectric;
string* fXsp;
// Generate Electric Potential Profile
string* genPotProfile;
// Use HYDROPRO
string* useHydropro;
// Optional Parameters Read File Reference for When Optional Parameters are set (makes .optionalParameters.txt)
// APBS Output Write Types (; delimited)
string apbsWriteTypes;
string dim_x;
string dim_y;
string dim_z;
string pdie;
// HYDROPRO Calculation Type
string calcType;
// PDB2PQR Force Field (PROPKA)
string forceField;
// MSMS Parameters
string lowPntDensity;
string highPntDensity;
};

QT_BEGIN_NAMESPACE
class QLineEdit;
QT_END_NAMESPACE
class Button;
# include <QDialog>

QT_BEGIN_NAMESPACE
class QCheckBox;
class QComboBox;
class QDateTimeEdit;
class QDial;
class QGroupBox;
class QLabel;
class QLineEdit;
class QProgressBar;
class QPushButton;
class QRadioButton;
class QScrollBar;
class QSlider;
class QSpinBox;
class QTabWidget;
class QTableWidget;
class QTextEdit;
QT_END_NAMESPACE

/*class apbsInstallInfo : public QWidget
{	
	Q_OBJECT
public:
	apbsInstallInfo(QWidget *parent = 0);
private:
};

class hydroproInstallInfo : public QWidget
{	
	Q_OBJECT
public:
	hydroproInstallInfo(QWidget *parent = 0);
private:
};

class msmsInstallInfo : public QWidget
{	
	Q_OBJECT
public:
	msmsInstallInfo(QWidget *parent = 0);
private:
};

class pInstallInfo : public QWidget
{	
	Q_OBJECT
public:
	pInstallInfo(QWidget *parent = 0);
private:
};

class zpredInstallInfo : public QWidget
{	
	Q_OBJECT
public:
	zpredInstallInfo(QWidget *parent = 0);
	void writeZPREDCFile(string oFile);
	
private:
};

class zpredProgressDisplay : public QWidget
{	
	Q_OBJECT
public:
	zpredProgressDisplay(QWidget *parent = 0);	
	QLabel* displayLabel;
	string displayFldr;
	int PID;
public slots:
	void killZPRED();
	void updateDisplay();
private:
};*/

class msgBox : public QWidget
{	
	Q_OBJECT
public:
	msgBox(QWidget *parent = 0);	
	QLabel* displayLabel;
private:
};

/*class apbsOptions : public QWidget
{
    Q_OBJECT
public:
	apbsOptions(QWidget *parent = 0);
	QLabel *dx;
	QLabel *dy;
	QLabel *dz;
	QLabel *proteinDielectric;
	QLineEdit *proteinDielectricInput;
	// APBS Write Type: pot Check Box
	QCheckBox *pot;
	// APBS Write Type: atompot Check Box
	QCheckBox *atompot;
	QCheckBox *dielx;
	QCheckBox *diely;
	QCheckBox *dielz;
	QCheckBox *charge;
	QCheckBox *edens;
	QCheckBox *ivdw;
	QCheckBox *lap;
	QCheckBox *ndens;
	QCheckBox *qdens;
	QCheckBox *smol;
	QCheckBox *sspl;
	QCheckBox *vdw;
	QCheckBox *kappa;
	//
	void defineParameters(string iFile);
public slots:
	void atompotChecked();
	void chargeChecked();
	void dielxChecked();
	void dielyChecked();
	void dielzChecked();
	void edensChecked();
	void ivdwChecked();
	void kappaChecked();
	void lapChecked();
	void ndensChecked();
	void qdensChecked();
	void smolChecked();
	void ssplChecked();
	void vdwChecked();
	void labelDX(const QString s);
	void labelDY(const QString s);
	void labelDZ(const QString s);
	void proteinDielectricDefined();
private:
};

class hydroproOptions : public QWidget
{
    Q_OBJECT
public:
	hydroproOptions(QWidget *parent = 0);
	// Calculation Type
	QLabel *calcType;
	void defineParameters(string iFile);
public slots:
	void labelCalcType(const QString s);
private:
};

class msmsOptions : public QWidget
{
    Q_OBJECT
public:
	msmsOptions(QWidget *parent = 0);
	// msms Low Surface Point Density
	QLabel *lowPntDensity;
	QLineEdit *lowPntDensityInput;
	// msms High Surface Point Density
	QLabel *highPntDensity;
	QLineEdit *highPntDensityInput;
	void defineParameters(string iFile);
public slots:
	void lowPntDensityDefined();
	void highPntDensityDefined();
private:
};

class pdb2pqrOptions : public QWidget
{
    Q_OBJECT
public:
	pdb2pqrOptions(QWidget *parent = 0);
	// Optional PROPKA Input: Force Field (AMBER, CHARMM, PARSE, PEOEPB, SWANSON, TYL06)
	QLabel *forceField;
	void defineParameters(string iFile);
public slots:
	void labelForceField(const QString s);
private:
};

class methodsInfo : public QWidget
{
    Q_OBJECT
public:
    methodsInfo(QWidget *parent = 0);
public slots:
private:
};

class solventInfo : public QWidget
{
    Q_OBJECT
public:
    solventInfo(QWidget *parent = 0);
public slots:
private:
};

class zpredInfo : public QWidget
{
    Q_OBJECT
public:
    zpredInfo(QWidget *parent = 0);
public slots:
private:
};

class zUI : public QWidget
{
    Q_OBJECT

public:
   zUI(QWidget *parent = 0);
	QGroupBox* create_oneSoluteBox(int tabIndex);
	QGroupBox* create_twoSoluteBox(int tabIndex);
	QGroupBox* create_threeSoluteBox(int tabIndex);
	QGroupBox* create_fourSoluteBox(int tabIndex);
	QGroupBox* create_fiveSoluteBox(int tabIndex);
	QGroupBox* create_sixSoluteBox(int tabIndex);
	QGroupBox* create_sevenSoluteBox(int tabIndex);
	QGroupBox* create_eightSoluteBox(int tabIndex);
	//
	QGroupBox* create_oneSolventBox(int tabIndex);
	QGroupBox* create_twoSolventBox(int tabIndex);
	QGroupBox* create_threeSolventBox(int tabIndex);
	QGroupBox* create_fourSolventBox(int tabIndex);
	QGroupBox* create_fiveSolventBox(int tabIndex);
	QGroupBox* create_sixSolventBox(int tabIndex);
	QGroupBox* create_sevenSolventBox(int tabIndex);
	QGroupBox* create_eightSolventBox(int tabIndex);
	//
	void writeZPREDInputFile(string oFile);
	void defineDownloadFolders(string iFile);
	void updateHistoryFile(string prgmNm,string Fldr);
	QGroupBox* create_fBox(int tabIndex);
	// Message Box Launcher
	void launchMsgBox(string errTxt);
	// Input Data Structure Array
	fromUIToZPRED uiInput;
	// Number of PDB Files Input by USer
	int numPDBs;
	// Array of PDB File Names
	string* pdbFiles;
	//
	QLabel* fileLabel;
	// Numbe of Solution Conditions Input by User for Zeta Potential Computation	
	int numSolCond=1;
	// Number of Solvents specified in each Solution Condition
	int* numSolvent;//=new int[numSolCond];
	// Number of Solutes specified in each Solution Condition
	int* numSolute;//=new int[numSolCond];
	// Solution Condition Tab
	QTabWidget* scTab;
	// Force User Specified Parameters instead of Computation
	QCheckBox* forceParmBox[MAX_SOLUTION_CONDITIONS];
	// Solvent Selection Group
	QGroupBox* solventBox;
	// Solute Selection Group
	QGroupBox* soluteBox;
	QPushButton* remSolCondButton;
	// Important Labels/Text Input
	QLineEdit* jobName;
	QLineEdit* pH[MAX_SOLUTION_CONDITIONS];	// pH
	QLineEdit* T[MAX_SOLUTION_CONDITIONS];		// temperature K
	QLineEdit* pC[MAX_SOLUTION_CONDITIONS];	// protein concentration g/L
	QLineEdit* fVisc[MAX_SOLUTION_CONDITIONS];// User-Specified Solution Viscosity Pa s	
	QLineEdit* fDens[MAX_SOLUTION_CONDITIONS];// User-Specified Solution Density kg/L	
	QLineEdit* fDiel[MAX_SOLUTION_CONDITIONS];// User-Specified relative Dielectric
	QLineEdit* fXsp[MAX_SOLUTION_CONDITIONS];	// User-Specified Hydration Layer Thickness [Angstroms]
	QCheckBox* potProfileBox[MAX_SOLUTION_CONDITIONS]; // CheckBox to generate electric potential profile	

	QLabel* jobNameLabel;
	QLabel* pHLabel[MAX_SOLUTION_CONDITIONS];
	QLabel* TLabel[MAX_SOLUTION_CONDITIONS];
	QLabel* pCLabel[MAX_SOLUTION_CONDITIONS];
	QLabel* seleSolvLabel[MAX_SOLUTION_CONDITIONS];				// Select Solvent Label
	QLabel* seleSolvConcLabel[MAX_SOLUTION_CONDITIONS];		// Select Solvent Concentration Label
	QLabel* seleSolvConcTypeLabel[MAX_SOLUTION_CONDITIONS]; // Select Solvent Concentration Type Label
	QLabel* seleSoluLabel[MAX_SOLUTION_CONDITIONS];				// Select Solute Label
	QLabel* seleSoluConcLabel[MAX_SOLUTION_CONDITIONS];		// Select Solute Concentration Label
	QLabel* seleSoluConcTypeLabel[MAX_SOLUTION_CONDITIONS]; // Select Solute Concentration Type Label
	QLabel* apbsLabel;
	QLabel* mvLabel;
	QLabel* hLabel;
	QLabel* msLabel;
	QLabel* msaLabel;
	QLabel* msxLabel;
	QLabel* msxnLabel;
	QLabel* pLabel;
	QLabel* zLabel;
	// Solute(s) Label Array
	QLabel* soluteLabel[MAX_SOLUTION_CONDITIONS][MAX_SOLUTES];
	// Formatted Solute(s) Label Array
	QLabel* fSoluteLabel[MAX_SOLUTION_CONDITIONS][MAX_SOLUTES];
	// Solvent(s) Label Array
	QLabel* solventLabel[MAX_SOLUTION_CONDITIONS][MAX_SOLVENTS];
	// Solute(s) Concentration(s) Line Value Array
	QLineEdit* soluteConcLine[MAX_SOLUTION_CONDITIONS][MAX_SOLUTES];
	// Solvent(s) Concentration(s) Line Value Array
	QLineEdit* solventConcLine[MAX_SOLUTION_CONDITIONS][MAX_SOLVENTS];
	// Solute(s) Concentration(s) Label Array
	QLabel* soluteConcLabel[MAX_SOLUTION_CONDITIONS][MAX_SOLUTES];
	// Solvent(s) Concentration(s) Label Array
	QLabel* solventConcLabel[MAX_SOLUTION_CONDITIONS][MAX_SOLVENTS];
	// Solute(s) Concentration(s) Type Label Array
	QLabel* soluteConcTypeLabel[MAX_SOLUTION_CONDITIONS];
	// Solvent(s) Concentration(s) Type Label Array
	QLabel* solventConcTypeLabel[MAX_SOLUTION_CONDITIONS];
	// File Reference for When Downloads already Performed (makes .history.txt)
	string refFile;
	// File Reference for When Optional Parameters are set (makes .optionalParameters.txt)
	string optRefFile;
	// APBS
	QLineEdit* apbsFilePath;
	string apbsFldr;
	string hydroproFldr;
	string msmsFldr;
	string pdb2pqrFldr;
	string zpredFldr;
	string filesFldr;
	// MULTIVALUE
	QLineEdit* mvFilePath;
	// HYDROPRO
	QLineEdit* hFilePath;
	// MSMS
	QLineEdit* msFilePath;
	// MSMS atmTypeNumbers File
	QLineEdit* msaFilePath;
	// MSMS XYZR File
	QLineEdit* msxFilePath;
	// MSMS XYZRN File
	QLineEdit* msxnFilePath;
	// PDB2PQR main.py File
	QLineEdit* pFilePath;
	// ZPRED Executable File Path
	QLineEdit* zFilePath;
	//
	static const int numModeledSolutes=59;
	//string modeledSolutes[numModeledSolutes]={"AlCl₃","Al₂(SO₄)₃","BaCl₂","CaCl₂","CdCl₂","CdSO₄","CoCl₂","CoSO₄","CrCl₃","Cr₂(SO₄)₃","CuCl₂","CuSO₄","FeCl₂","FeSO₄","FeCl₃","Fe₂(SO₄)₃","HCl","HCN","HNO₃","H₃PO₄","H₂SO₄","KCl","K₂CO₃","KNO₃","KOH","K₂SO₄","LiCl","Li₂SO₄","MgCl₂","MgSO₄","MnCl₂","MnSO₄","NaBr","NaCl","NaClO₃","Na₂CO₃","NaF","NaHCO₃","NaH₂PO₄","Na₂HPO₄","NaHSO₃","NaI","Na₂MoO₄","NaNO₂","NaNO₃","NaOH","Na₃PO₄","Na₂SO₃","Na₂S₂O₃","Na₂SO₄","NH₃","NH₄Cl","NH₄NO₃","(NH₄)₂SO₄","NiCl₂","NiSO₄","SrCl₂","ZnCl₂","ZnSO₄"};
	//string modeledSolutes[numModeledSolutes]={"AlCl3","Al2(SO4)3","BaCl2","CaCl2","CdCl2","CdSO4","CoCl2","CoSO4","CrCl3","Cr2(SO4)3","CuCl2","CuSO4","FeCl2","FeSO4","FeCl3","Fe2(SO4)3","HCl","HCN","HNO3","H3PO4","H2SO4","KCl","K2CO3","KNO3","KOH","K2SO4","LiCl","Li2SO4","MgCl2","MgSO4","MnCl2","MnSO4","NaBr","NaCl","NaClO3","Na2CO3","NaF","NaHCO3","NaH2PO4","Na3HPO4","NaHSO3","NaI","Na2MoO4","NaNO2","NaNO3","NaOH","Na3PO4","Na2SO3","Na2S2O3","Na2SO4","NH3","NH4Cl","NH4NO3","(NH4)2SO4","NiCl2","NiSO4","SrCl2","ZnCl2","ZnSO4"};
	//string modeledSolutes[numModeledSolutes]={"AlCl3","Al2(SO4)3","BaCl2","CaCl2","CdCl2","CdSO4","CoCl2","CoSO4","CrCl3","Cr2(SO4)3","CuCl2","CuSO4","FeCl2","FeSO4","FeCl3","Fe2(SO4)3","HCl","HCN","HNO3","H3PO4","H2SO4","KCl","KClO4","K2CO3","KH2PO4","K2HPO4","KNO3","KOH","K2SO4","LiCl","Li2SO4","MgCl2","MgSO4","MnCl2","MnSO4","NaBr","NaCl","NaClO3","Na2CO3","NaF","NaHCO3","NaH2PO4","Na2HPO4","NaHSO3","NaI","Na2MoO4","NaNO2","NaNO3","NaOH","Na3PO4","Na2SO3","Na2S2O3","Na2SO4","NH3","NH4Cl","NH4NO3","(NH4)2SO4","NiCl2","NiSO4","SrCl2","ZnCl2","ZnSO4"};
	string modeledSolutes[numModeledSolutes]={"AlCl3","Al2(SO4)3","BaCl2","CaCl2","CdCl2","CdSO4","CoCl2","CoSO4","CuCl2","CuSO4","FeCl2","FeSO4","FeCl3","HCl","HNO3","H3PO4","H2SO4","KCl","KClO4","K2CO3","KH2PO4","K2HPO4","KNO3","KOH","K2SO4","LiCl","Li2SO4","MgCl2","MgSO4","MnCl2","MnSO4","NaBr","NaCl","NaClO3","Na2CO3","NaF","NaHCO3","NaH2PO4","Na2HPO4","NaHSO3","NaI","Na2MoO4","NaNO2","NaNO3","NaOH","Na3PO4","Na2SO3","Na2S2O3","Na2SO4","NH3","NH4Cl","NH4NO3","(NH4)2SO4","NiCl2","NiSO4","SrCl2","ZnCl2","ZnSO4","H3Citrate"};
	//
	static const int numModeledSolvents=10;
	string modeledSolvents[numModeledSolvents]={"acetic acid","acetone","butanol","dimethyl sulfoxide","ethanol","glycerol","isopropanol","methanol","propanol","water"};	
	//
	static const int numModeledSolventConcTypes=3;
	string modeledSolventConcTypes[numModeledSolventConcTypes]={"Mass Fraction","Mole Fraction","Volume Fraction"};
	// saveTemps Check Box (boolean for saving all files used in computation)
	QCheckBox *saveTempsBox;
	QCheckBox *hydroproBox;	// Specify whether to use HYDROPRO or HULLRAD
	// Number of Parallel Computations
	QLineEdit* maxThreadsLine;
	QLabel* maxThreadsLineLabel;
	apbsOptions *AO;
	apbsInstallInfo* AII;
	hydroproOptions *HO;
	hydroproInstallInfo* HII;
	msmsOptions *MO;
	msmsInstallInfo* MII;
	pdb2pqrOptions *PO;
	pInstallInfo* PII;
	zpredInstallInfo* ZII;
	zpredProgressDisplay* ZPD;
	msgBox* MB;
	// ZPRED Process ID
	int zpredPID;
public slots:
	void calcZetaPotential();
	void openAPBS();
	void openHYDROPRO();
	void openMSMS();
	void openMSMS_A();
	void openMSMS_X();
	void openMSMS_XN();
	void openMULTIVALUE();
	void openPDB2PQR();
	void openZPRED();
	void openPdbFiles();
	void addSolCondTab();
	void remSolCondTab();
	void addSolute();
	void addSolvent();
	void remSolute();
	void remSolvent();
	// GUI Update
	void maxThreads_defined();
	void jobName_defined();
	void pH_defined();
	void T_defined();
	void pC_defined();
	void solvent_defined(const QString&);
	void solventConc_defined();
	void solventConcType_defined(const QString&);
	void solute_defined(const QString&);
	void soluteConc_defined();
	void soluteConcType_defined(const QString&);
	// For Specifiying Parameters
	void show_fBox();
	// APBS
	void apbsFilePath_defined();
	// MULTIVALUE
	void mvFilePath_defined();
	// HYDROPRO
	void hFilePath_defined();
	// MSMS
	void msFilePath_defined();
	// MSMS atmTypeNumbers File
	void msaFilePath_defined();
	// MSMS XYZR File
	void msxFilePath_defined();
	// MSMS XYZRN File
	void msxnFilePath_defined();
	// PDB2PQR main.py File
	void pFilePath_defined();
	// ZPRED
	void zFilePath_defined();
	// AutoFill FilePaths
	void autoFillFilePaths();
	// Experimental Methods Launcher
	void launchMethodsInfo();	
	// Solvent Info Launcher
	void launchSolventInfo();
	// Valildation Info Launcher
	void launchValidationInfo();
	// ZPRED Info Launcher
	void launchZpredInfo();
	void launchZpredInstallInfo();
	void launchZpredProgressDisplay();
	// APBS Options Launcher
	void launchApbsOptions();
	void launchApbsInstallInfo();
	// HYDROPRO Options Launcher
	void launchHydroproOptions();
	void launchHydroproInstallInfo();
	// MSMS options Launcher
	void launchMsmsOptions();
	void launchMsmsInstallInfo();
	// PDB2PQR options Launcher
	void launchPdb2pqrOptions();
	void launchPdb2pqrInstallInfo();
	// Solute Name
	void label1solute(int index);
	void label2solute(int index);
	void label3solute(int index);
	void label4solute(int index);
	void label5solute(int index);
	void label6solute(int index);
	void label7solute(int index);
	void label8solute(int index);
	// Solute Concentration
	void label1soluteConc();
	void label2soluteConc();
	void label3soluteConc();
	void label4soluteConc();
	void label5soluteConc();
	void label6soluteConc();
	void label7soluteConc();
	void label8soluteConc();
	// Solute Concentration Type
	void label1soluteConcType(const QString& s);
	// Solvent Name
	void label1solvent(const QString& s);
	void label2solvent(const QString& s);
	void label3solvent(const QString& s);
	void label4solvent(const QString& s);
	void label5solvent(const QString& s);
	void label6solvent(const QString& s);
	void label7solvent(const QString& s);
	void label8solvent(const QString& s);
	// Solvent Concentration
	void label1solventConc();
	void label2solventConc();
	void label3solventConc();
	void label4solventConc();
	void label5solventConc();
	void label6solventConc();
	void label7solventConc();
	void label8solventConc();
	// Solvent Concentration Type
	void label1solventConcType(const QString& s);

private:
	methodsInfo *MI;
	solventInfo *SI;
	zpredInfo *ZI;	
};*/

// makeGif UI
/*class gUI : public QWidget
{
    Q_OBJECT
public:
   gUI(QWidget *parent = 0);
	void write_pymolFunction_File(bool xROT,bool yROT,bool zROT,string pyFnNm,string outFile);
	QLabel* pdbLabel;
	int numPDBs;
	string* pdbFiles;
	QCheckBox* xRotationBox;
	QCheckBox* yRotationBox;
	QCheckBox* zRotationBox;
public slots:
	void makeGif();
	void open_pdbFiles();
private:

};*/

// protCAD UI
class pUI : public QWidget
{
    Q_OBJECT

public:
   pUI(QWidget *parent = 0);

	void write_protEvolver_pymolFunction_File(string pyFnNm,string outFile);
	// protAlign
	QLabel* protAlignPDBLabel1;
	QLabel* protAlignPDBLabel2;
	string protAlign_pdbFile1;
	string protAlign_pdbFile2;
	// protDielectric

	// protEvolver
	QLabel* protEvolverPDBLabel;
	QLabel* protEvolverActiveChainLabel;
	QLabel* protEvolverActivePositionLabel;
	QLabel* protEvolverRandomPositionLabel;
	QLabel* protEvolverFrozenPositionLabel;
	QLabel* protEvolverAminoAcidLabel;
	QLabel* maxThreadsLabel;
	QLineEdit* protEvolverActiveChainInput;
	QLineEdit* protEvolverActivePositionInput;
	QLineEdit* protEvolverRandomPositionInput;
	QLineEdit* protEvolverFrozenPositionInput;
	QLineEdit* protEvolverAminoAcidInput;
	QLineEdit* maxThreadsLine;
	QCheckBox* protEvolverRelaxationBox;
	string protEvolver_pdbFile;
	// makeGif UI
//	gUI* GIF;
	// ZPRED UI
//	zUI* ZPRED;
public slots:
	void open_protAlignPDBFile1();
	void open_protAlignPDBFile2();
	void open_protEvolverPDBFile();
	void protEvolverActiveChainInput_defined();
	void protEvolverActivePositionInput_defined();
	void protEvolverRandomPositionInput_defined();
	void protEvolverFrozenPositionInput_defined();
	void protEvolverAminoAcidInput_defined();
//	void runMakeGif();
//	void runZPRED();
	void runProtEvolver();
	void view();
private:

};

string cnvrtNumToStrng(float Num,int numberAfterDecimalpoint);
string cnvrtNumToStrng(double Num,int numberAfterDecimalpoint);
string cnvrtNumToStrng(long double Num,int numberAfterDecimalpoint);
string cnvrtNumToStrng(int Num,int numberAfterDecimalpoint);
string cnvrtNumToStrng(unsigned int Num,int numberAfterDecimalpoint);
string* fill_string_array(string Data,int numPnts,string delimiter);
int* fill_int_array(string Data,int numPnts,string delimiter);
double* fill_double_array(string Data,int numPnts,string delimiter);
string fixSubscriptNumbers(string s);
string getBaseFolder(string f);
string makeUpperCase(string X);
string checkFinalBackSlash(string s);
string setStringWidth(string In,int width);
void updateParameterFile(string Name,string Value);
double interpolate(double x,double x1,double x2,double y1,double y2);
double getSoluteAqueousSolubility(string solute,double T);
string base64_decode(string const& encoded_string);
void decryptString(string encStr,const string salt,const string pass,string outFile);
void handleErrors(void);

# endif

//⁰ ¹ ² ³ ⁴ ⁵ ⁶ ⁷ ⁸ ⁹ ⁺ ⁻ ⁼ ⁽ ⁾
//ᵃ ᵇ ᶜ ᵈ ᵉ ᶠ ᵍ ʰ ⁱ ʲ ᵏ ˡ ᵐ ⁿ ᵒ ᵖ ʳ ˢ ᵗ ᵘ ᵛ ʷ ˣ ʸ ᶻ 
//ᴬ ᴮ ᴰ ᴱ ᴳ ᴴ ᴵ ᴶ ᴷ ᴸ ᴹ ᴺ ᴼ ᴾ ᴿ ᵀ ᵁ ⱽ ᵂ 
//ₐ ₑ ₕ ᵢ ⱼ ₖ ₗ ₘ ₙ ₒ ₚ ᵣ ₛ ₜ ᵤ ᵥ ₓ
//ᵅ ᵝ ᵞ ᵟ ᵋ ᶿ ᶥ ᶲ ᵠ ᵡ 
//ᵦ ᵧ ᵨ ᵩ ᵪ  
//1₀ ₁ ₂ ₃ ₄ ₅ ₆ ₇ ₈ ₉ ₊ ₋ ₌ ₍ ₎
