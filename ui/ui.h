# include <string>

// Computational Constants
# define LINE_INPUT_SIZE 100
# define NO_VALUE "----"

using namespace std;

# ifndef UI_H
# define UI_H

# include <QWidget>

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


class msgBox : public QWidget
{	
	Q_OBJECT
public:
	msgBox(QWidget *parent = 0);	
	QLabel* displayLabel;
private:
};

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
	string protEvolver_path;
	QPushButton* protEvolverPDBButton;
	QPushButton* xButton3;
public slots:
	void open_protAlignPDBFile1();
	void open_protAlignPDBFile2();
	void open_protEvolverPDBFile();
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
vector<string> split (const string &s, char delim);
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
