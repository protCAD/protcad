# include <cmath>
# include <cstdio>
# include <cstring>
# include <fstream>
# include <iostream>
# include <sstream>
# include <string>
# include <thread>
# include <vector>
# include <unistd.h>
# include <QtWidgets>
# include <QProcess>

# include "ui.h"

using namespace std;

// protCAD UI
pUI::pUI(QWidget *parent) : QWidget(parent)
{
	// Initialize Variables to Track when User Changes their value
	protAlign_pdbFile1="";
	protAlign_pdbFile2="";
	
	QFont font;
	font.setBold(true);
	QTabWidget* theTab=new QTabWidget;

	//-----protALIGN (tab1)
	QWidget *tPg1 = new QWidget;
	// Labels
	QLabel* pg1Label=new QLabel(tr("Align PDB file #1 to PDB file #2\n")); pg1Label->setAlignment(Qt::AlignLeft);
	pg1Label->setWordWrap(true);
	protAlignPDBLabel1=new QLabel;
	protAlignPDBLabel1->setText("<b><span style=\"color:red;\">Select PDB file #1</span>:</b>");
	protAlignPDBLabel1->setAlignment(Qt::AlignLeft);
	protAlignPDBLabel2=new QLabel;
	protAlignPDBLabel2->setText("<b><span style=\"color:red;\">Select PDB file #2</span>:</b>");
	protAlignPDBLabel2->setAlignment(Qt::AlignLeft);
	// Buttons
	QPushButton* protAlignPDBButton1=new QPushButton(tr("..."));	
	protAlignPDBButton1->setToolTip(tr("Tooltip"));
	protAlignPDBButton1->setFixedWidth(40);
	protAlignPDBButton1->setFont(font);
	connect(protAlignPDBButton1,SIGNAL(clicked()),this,SLOT(open_protAlignPDBFile1()));
	QPushButton* protAlignPDBButton2=new QPushButton(tr("..."));
	protAlignPDBButton2->setToolTip(tr("Tooltip"));
	protAlignPDBButton2->setFixedWidth(40);
	protAlignPDBButton2->setFont(font);
	connect(protAlignPDBButton2,SIGNAL(clicked()),this,SLOT(open_protAlignPDBFile2()));
	QPushButton* xButton1=new QPushButton(tr("EXECUTE"));
	xButton1->setFixedWidth(200);
	xButton1->setFont(font);
	// Design Tab Page
	QGridLayout *pg1box = new QGridLayout;
	pg1box->setSizeConstraint(QLayout::SetFixedSize);
	pg1box->addWidget(pg1Label,0,0,1,5);
	pg1box->addWidget(protAlignPDBLabel1,1,0,1,1);pg1box->addWidget(protAlignPDBButton1,1,1,1,1);pg1box->addWidget(protAlignPDBLabel2,1,3,1,1);pg1box->addWidget(protAlignPDBButton2,1,4,1,1);
	pg1box->addWidget(xButton1,2,0,1,5);
	// Set Layout
	tPg1->setLayout(pg1box);
	// Update Tab Widget
	theTab->addTab(tPg1,tr("protAlign"));

	//----protDIELECTRIC (tab2)
	QWidget *tPg2 = new QWidget;
	// Labels
	QLabel* pg2Label=new QLabel(tr("Description Label\n")); pg2Label->setAlignment(Qt::AlignLeft);
	pg2Label->setWordWrap(true);
	// Buttons
	QPushButton* xButton2=new QPushButton(tr("EXECUTE"));
	xButton2->setFixedWidth(200);
	xButton2->setFont(font);
	// Design Tab Page Layout
	QGridLayout *pg2box = new QGridLayout;
	pg2box->setSizeConstraint(QLayout::SetFixedSize);
	pg2box->addWidget(pg2Label,0,0,1,2);
	pg2box->addWidget(xButton2,1,0,1,2);
	// Set Layout
	tPg2->setLayout(pg2box);
	// Update Tab Widget
	theTab->addTab(tPg2,tr("protDielectric"));

	//-----protEVOLVER (tab3)
	QWidget *tPg3 = new QWidget;
	// Labels
	QLabel* pg3Label=new QLabel(tr("Evolve a sequence for a fold\n")); pg3Label->setAlignment(Qt::AlignLeft);
	pg3Label->setWordWrap(true);
	protEvolverPDBLabel=new QLabel;
	protEvolverPDBLabel->setText("<b><span style=\"color:black;\">Select PDB file</span>:</b>");
	protEvolverPDBLabel->setAlignment(Qt::AlignLeft);
	protEvolverActiveChainLabel=new QLabel;
	protEvolverActiveChainLabel->setText("<b><span style=\"color:black;\">Specify active chain(s)</span>:</b>");
	protEvolverActiveChainLabel->setAlignment(Qt::AlignLeft);
	protEvolverActivePositionLabel=new QLabel;
	protEvolverActivePositionLabel->setText("<b><span style=\"color:black;\">Specify active position(s)</span>:</b>");
	protEvolverActivePositionLabel->setAlignment(Qt::AlignLeft);
	protEvolverRandomPositionLabel=new QLabel;
	protEvolverRandomPositionLabel->setText("<b><span style=\"color:black;\">Specify random position(s)</span>:</b>");
	protEvolverRandomPositionLabel->setAlignment(Qt::AlignLeft);
	protEvolverFrozenPositionLabel=new QLabel;
	protEvolverFrozenPositionLabel->setText("<b><span style=\"color:black;\">Specify frozen position(s)</span>:</b>");
	protEvolverFrozenPositionLabel->setAlignment(Qt::AlignLeft);
	protEvolverAminoAcidLabel=new QLabel;
	protEvolverAminoAcidLabel->setText("<b><span style=\"color:black;\">Specify Amino Acids</span>:</b>");
	protEvolverAminoAcidLabel->setAlignment(Qt::AlignLeft);
	maxThreadsLabel=new QLabel;
	maxThreadsLabel->setText("<b><span style=\"color:black;\">Specify Number of Threads</span>:</b>");
	maxThreadsLabel->setAlignment(Qt::AlignLeft);
	// Line Input
	protEvolverActiveChainInput=new QLineEdit;
	protEvolverActiveChainInput->setAlignment(Qt::AlignCenter);
	//protEvolverActiveChainInput->setText(tr(NO_VALUE));
	protEvolverActiveChainInput->setFixedWidth(250);
	//connect(protEvolverActiveChainInput,SIGNAL(returnPressed()),this,SLOT(protEvolverActiveChainInput_defined()));
	protEvolverActivePositionInput=new QLineEdit;
	protEvolverActivePositionInput->setAlignment(Qt::AlignCenter);
	//protEvolverActivePositionInput->setText(tr(NO_VALUE));
	protEvolverActivePositionInput->setFixedWidth(250);
	//connect(protEvolverActivePositionInput,SIGNAL(returnPressed()),this,SLOT(protEvolverActivePositionInput_defined()));
	protEvolverRandomPositionInput=new QLineEdit;
	protEvolverRandomPositionInput->setAlignment(Qt::AlignCenter);
	//protEvolverRandomPositionInput->setText(tr(NO_VALUE));
	protEvolverRandomPositionInput->setFixedWidth(250);
	//connect(protEvolverRandomPositionInput,SIGNAL(returnPressed()),this,SLOT(protEvolverRandomPositionInput_defined()));
	protEvolverFrozenPositionInput=new QLineEdit;
	protEvolverFrozenPositionInput->setAlignment(Qt::AlignCenter);
	//protEvolverFrozenPositionInput->setText(tr(NO_VALUE));
	protEvolverFrozenPositionInput->setFixedWidth(250);
	//connect(protEvolverFrozenPositionInput,SIGNAL(returnPressed()),this,SLOT(protEvolverFrozenPositionInput_defined()));
	protEvolverAminoAcidInput=new QLineEdit;
	protEvolverAminoAcidInput->setAlignment(Qt::AlignCenter);
	//protEvolverAminoAcidInput->setText(tr(NO_VALUE));
	protEvolverAminoAcidInput->setFixedWidth(250);
	//connect(protEvolverAminoAcidInput,SIGNAL(returnPressed()),this,SLOT(protEvolverAminoAcidInput_defined()));
	// Check Boxes
	protEvolverRelaxationBox=new QCheckBox(tr("Allow Backbone Relaxation?")); protEvolverRelaxationBox->setChecked(false); protEvolverRelaxationBox->setFont(font);
	// Number of THreads LineEdit/INput
	unsigned int number_of_threads = thread::hardware_concurrency();
	//std::cout << n << " concurrent threads are supported.\n";
	string numThreadValue;
	if(number_of_threads!=0)
		{numThreadValue=cnvrtNumToStrng(number_of_threads,0);}
	else
		{numThreadValue="4";}
	maxThreadsLine=new QLineEdit();
	maxThreadsLine->setAlignment(Qt::AlignCenter);
	maxThreadsLine->setText(tr(numThreadValue.c_str()));
	maxThreadsLine->setMaxLength(3);
	maxThreadsLine->setFixedWidth(80);
	// Buttons
	QPushButton* viewButton=new QPushButton(tr("View"));
	connect(viewButton,SIGNAL(clicked()),this,SLOT(view()));
	protEvolverPDBButton=new QPushButton(tr("..."));	
	protEvolverPDBButton->setToolTip(tr("Tooltip"));
	protEvolverPDBButton->setFixedWidth(80);
	protEvolverPDBButton->setFont(font);
	connect(protEvolverPDBButton,SIGNAL(clicked()),this,SLOT(open_protEvolverPDBFile()));
	xButton3=new QPushButton(tr("RUN"));	
	xButton3->setFixedWidth(200);
	xButton3->setFont(font);
	xButton3->setCheckable(true);
	xButton3->setChecked(false);
	connect(xButton3,SIGNAL(clicked()),this,SLOT(runProtEvolver()));
	// Design Tab Page Layout
	QGridLayout *pg3box = new QGridLayout;
	pg3box->setSizeConstraint(QLayout::SetFixedSize);
	pg3box->addWidget(pg3Label,0,0,1,3);
	pg3box->addWidget(protEvolverPDBLabel,1,0,1,1);pg3box->addWidget(protEvolverPDBButton,1,1,1,1);pg3box->addWidget(viewButton,1,2,1,1);
	pg3box->addWidget(protEvolverActiveChainLabel,2,0,1,1);pg3box->addWidget(protEvolverActiveChainInput,2,1,1,1);
	pg3box->addWidget(protEvolverActivePositionLabel,3,0,1,1);pg3box->addWidget(protEvolverActivePositionInput,3,1,1,1);
	pg3box->addWidget(protEvolverRandomPositionLabel,4,0,1,1);pg3box->addWidget(protEvolverRandomPositionInput,4,1,1,1);
	pg3box->addWidget(protEvolverFrozenPositionLabel,5,0,1,1);pg3box->addWidget(protEvolverFrozenPositionInput,5,1,1,1);
	pg3box->addWidget(protEvolverAminoAcidLabel,6,0,1,1);pg3box->addWidget(protEvolverAminoAcidInput,6,1,1,1);
	pg3box->addWidget(maxThreadsLabel,7,0,1,1);pg3box->addWidget(maxThreadsLine,7,1,1,1);
	pg3box->addWidget(protEvolverRelaxationBox,8,0,1,3);
	pg3box->addWidget(xButton3,9,0,1,3);
	// Set Layout
	tPg3->setLayout(pg3box);
	// Update Tab Widget
	theTab->addTab(tPg3,tr("protEvolver"));

	//---main layout
	QGridLayout *mainLayout = new QGridLayout;
	mainLayout->setSizeConstraint(QLayout::SetFixedSize);
	mainLayout->addWidget(theTab,0,0);
	setLayout(mainLayout);
	setWindowTitle(tr("protCAD"));
	resize(QDesktopWidget().availableGeometry(this).size());
}

void pUI::open_protAlignPDBFile1()
{
	QString Fldr=QDir::currentPath();
	string sFldr=Fldr.toStdString()+"/";
	QStringList filenames = QFileDialog::getOpenFileNames(this,tr("files"),tr(sFldr.c_str()),tr("All files (*)") );
	string tmp;
	if(filenames.count()!=0)
		{tmp=filenames.at(0).toLocal8Bit().constData();
		protAlign_pdbFile1=tmp;
		protAlignPDBLabel1->setText("<b><span style=\"color:black;\">Select PDB file #1</span>:</b>");}
}

void pUI::open_protAlignPDBFile2()
{
	QString Fldr=QDir::currentPath();
	string sFldr=Fldr.toStdString()+"/";
	QStringList filenames = QFileDialog::getOpenFileNames(this,tr("files"),tr(sFldr.c_str()),tr("All files (*)") );
	string tmp;
	if(filenames.count()!=0)
		{tmp=filenames.at(0).toLocal8Bit().constData();
		protAlign_pdbFile2=tmp;
		protAlignPDBLabel2->setText("<b><span style=\"color:black;\">Select PDB file #2</span>:</b>");}
}

void pUI::open_protEvolverPDBFile()
{
	QString Fldr=QDir::currentPath();
	string sFldr=Fldr.toStdString()+"/";
	QStringList filenames = QFileDialog::getOpenFileNames(this,tr("files"),tr(sFldr.c_str()),tr("All files (*)") );
	string tmp;
	if(filenames.count()!=0)
	{
		tmp=filenames.at(0).toLocal8Bit().constData();
		protEvolver_pdbFile=tmp;
		string name, path="";
		vector<string> v = split (tmp, '/');
		for (unsigned int i = 0; i < v.size(); i++) 
		{
			name=v[i]; if (i< v.size()-1){path+=v[i];path+="/";}
		}
		protEvolver_path = path;
		protEvolverPDBButton->setText(QString::fromStdString(name));
		protEvolverPDBLabel->setText("<b><span style=\"color:black;\">Select PDB file</span>:</b>");
	}
}

void pUI::runProtEvolver()
{	
	if (xButton3->isChecked())
	{
		xButton3->setChecked(true);
		xButton3->setText("STOP");
		//write input file from UI feilds
		string data="";
		data+="PDB file,";
		data+=protEvolver_pdbFile;
		data+=",\n";
		data+="Active Chains,";
		data+=protEvolverActiveChainInput->text().toStdString();
		data+="\n";
		data+="Active Positions,";
		data+=protEvolverActivePositionInput->text().toStdString();
		data+="\n";
		data+="Random Positions,";
		data+=protEvolverRandomPositionInput->text().toStdString();
		data+="\n";
		data+="Frozen Positions,";
		data+=protEvolverFrozenPositionInput->text().toStdString();
		data+="\n";
		data+="Amino Acids,";
		data+=protEvolverAminoAcidInput->text().toStdString();
		data+="\n";
		data+="Backbone Relaxation,";
		if(protEvolverRelaxationBox->isChecked()){data+="true";}
		else{data+="false";}
		data+=",\n";
		
		string sFldr=protEvolver_path;
		string inputFile=sFldr+"evolver.in";
		ofstream fOut;
		fOut.open(inputFile.c_str());
		if(fOut.fail()){}
		fOut<<data;fOut.close();

		string tmp=maxThreadsLine->text().toStdString();
		int nT=atoi(tmp.c_str());
		string cmd="cd "+protEvolver_path+" && protEvolver "+inputFile+" &> evolver.log &";	
		for(int i=0;i<nT;i++)
		{
			if (i == 1){usleep(500000);}
			int statusCode=system(cmd.c_str());
			if (statusCode == -1)
			{fprintf(stderr, "program failed to run, errno = %d\n", errno);}
		}
	}
	else
	{
		xButton3->setChecked(false);xButton3->setText("RUN");
		string cmd="killall protEvolver";
		int statusCode=system(cmd.c_str());
		if (statusCode == -1){fprintf(stderr, "program failed to run, errno = %d\n", errno);}
	}
}

void pUI::write_protEvolver_pymolFunction_File(string pyFnNm,string outFile)
{
	QString theCurFldr=QDir::currentPath();
	string sFldr=theCurFldr.toStdString()+"/PyMol_Images/";
	string Output="";
	// Write PyMol Function File
	Output+="#!/usr/bin/python\nfrom pymol import stored\nfrom time import sleep\n\n# Function Definition\n\n";
	Output+="def "+pyFnNm+"():\n\n";
	// Define PDB File(s)	
	Output+="\tpFile=\""+protEvolver_pdbFile+"\"\n";
	Output+="\n";
	// Load PDB File(s)
	Output+="\tcmd.load(pFile)\n";
	Output+="\n";
	// More Functions
	Output+="\tcmd.do(\"util.cbag\")\n";
	Output+="\tcmd.zoom(\"all\",10)\n";

	Output+="\n";
	// End File
	Output+="\treturn\n\ncmd.extend(\""+pyFnNm+"\","+pyFnNm+")";

	ofstream fOut;
	fOut.open(outFile.c_str());
	if(fOut.fail()){cerr<<"Error in write_pymolFunction_File!\nCould not open file ("<<outFile<<")\n";}//exit(EXIT_FAILURE);}
	else
		{fOut<<Output;fOut.close();}
}

void pUI::view()
{
	QString theCurFldr=QDir::currentPath();
	string sFldr=theCurFldr.toStdString()+"/";

	string pythonFunctionName="Test";
	string pythonFunctionFile=sFldr+"Test.py";
	string cmd="pymol -d run "+pythonFunctionFile+" -d "+pythonFunctionName;
	
	write_protEvolver_pymolFunction_File(pythonFunctionName,pythonFunctionFile);		
	
	int statusCode=system(cmd.c_str());
	remove(pythonFunctionFile.c_str());
	if (statusCode == -1)
		{fprintf(stderr, "pymol failed to run, errno = %d\n", errno);}
}

string cnvrtNumToStrng(int Num,int numberAfterDecimalpoint)
{	
	stringstream ss;
	ss.setf(ios::fixed);
	if(numberAfterDecimalpoint>0)
		{ss.setf(ios::showpoint);}
	ss.precision(numberAfterDecimalpoint);
	ss<<Num;
	return ss.str();
}

string cnvrtNumToStrng(unsigned int Num,int numberAfterDecimalpoint)
{	
	stringstream ss;
	ss.setf(ios::fixed);
	if(numberAfterDecimalpoint>0)
		{ss.setf(ios::showpoint);}
	ss.precision(numberAfterDecimalpoint);
	ss<<Num;
	return ss.str();
}

string cnvrtNumToStrng(double Num,int numberAfterDecimalpoint){stringstream ss;ss.setf(ios::fixed);if(numberAfterDecimalpoint>0){ss.setf(ios::showpoint);}ss.precision(numberAfterDecimalpoint);ss<<Num;return ss.str();}

string cnvrtNumToStrng(long double Num,int numberAfterDecimalpoint){stringstream ss;ss.setf(ios::fixed);if(numberAfterDecimalpoint>0){ss.setf(ios::showpoint);}ss.precision(numberAfterDecimalpoint);ss<<Num;return ss.str();}

string cnvrtNumToStrng(float Num,int numberAfterDecimalpoint){stringstream ss;ss.setf(ios::fixed);if(numberAfterDecimalpoint>0){ss.setf(ios::showpoint);}ss.precision(numberAfterDecimalpoint);ss<<Num;return ss.str();}

double* fill_double_array(string Data,int numPnts,string delimiter){double* Output=new double[numPnts];string bld="",tmp="";int Counter=0;for(uint i=0;i<Data.length();i++){tmp=Data[i];if(tmp.compare(delimiter)==0 && Counter<numPnts){Output[Counter]=strtod(bld.c_str(),NULL);Counter++;bld="";}else{bld+=Data[i];}}return Output;}

int* fill_int_array(string Data,int numPnts,string delimiter){int* Output=new int[numPnts];string bld="",tmp="";int Counter=0;for(uint i=0;i<Data.length();i++){tmp=Data[i];if(tmp.compare(delimiter)==0 && Counter<numPnts){Output[Counter]=atoi(bld.c_str());Counter++;bld="";}else{bld+=Data[i];}}return Output;}

string* fill_string_array(string Data,int numPnts,string delimiter){string* Output=new string[numPnts];string bld="",tmp="";int Counter=0;for(uint i=0;i<Data.length();i++){tmp=Data[i];if(tmp.compare(delimiter)==0 && Counter<numPnts){Output[Counter]=bld;Counter++;bld="";}else{bld+=Data[i];}}return Output;}

string fixSubscriptNumbers(string s)
{
	string Out;
	if(s.compare("0")==0){Out="₀";}	
	else if(s.compare("1")==0){Out="₁";}
	else if(s.compare("2")==0){Out="₂";}
	else if(s.compare("3")==0){Out="₃";}
	else if(s.compare("4")==0){Out="₄";}
	else if(s.compare("5")==0){Out="₅";}
	else if(s.compare("6")==0){Out="₆";}
	else if(s.compare("7")==0){Out="₇";}
	else if(s.compare("8")==0){Out="₈";}
	else if(s.compare("9")==0){Out="₉";}
	else{Out=s;}
	return Out;
}

string getBaseFolder(string f)
{
	int pos=f.rfind("/",f.length()-1);
	return f.substr(0,pos)+"/";
}

string makeUpperCase(string X)
{
	string Output="";char letter;
	for(uint i=0;i<X.length();i++)
	{
		letter=X[i];Output+=toupper(letter);
	}
	return Output;
}

string checkFinalBackSlash(string s)
{
	string tmp=s.substr(s.length()-1,1);
	if(tmp.compare("/")!=0){tmp=s+"/";}
	else{tmp=s;}
	return tmp;
}

string setStringWidth(string In,int width)
{
	int Counter=0;
	string Output="";
	for(uint i=0;i<In.length();i++)
		{if(Counter>=width)
			{Output+="\n";
			Output+=In[i];
			Counter=0;
			}
		else
			{Output+=In[i];
			Counter++;}
		}
	return Output;
}

vector<string> split (const string &s, char delim) 
{
    vector<string> result;
    stringstream ss (s);
    string item;

    while (getline (ss, item, delim)) {
        result.push_back (item);
    }

    return result;
}

static inline bool is_base64(unsigned char c){return (isalnum(c) || (c == '+') || (c == '/'));}