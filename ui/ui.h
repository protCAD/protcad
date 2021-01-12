# include <fstream>
# include <unistd.h>
# include <iostream>
# include <sstream>
# include <string>
# include <vector>
# include <thread>
# include <QWidget>
# include <QProcess>
# include <QtWidgets>
# include <QDialog>

// Computational Constants
# define LINE_INPUT_SIZE 100
# define NO_VALUE "----"

using namespace std;

# ifndef UI_H
# define UI_H

QT_BEGIN_NAMESPACE
class QLineEdit;
QT_END_NAMESPACE
class Button;


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
	QCheckBox* protEvolverPolarityBox;
	string protEvolver_pdbFile;
	string protEvolver_path;
	QPushButton* protEvolverPDBButton;
	QPushButton* xButton2;
	QProcess * process = new QProcess;
public slots:
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
string makeUpperCase(string X);
string checkFinalBackSlash(string s);
string setStringWidth(string In,int width);

# endif
