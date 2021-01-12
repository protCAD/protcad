# include "ui.h"

using namespace std;

// protCAD UI
pUI::pUI(QWidget *parent) : QWidget(parent)
{	
	QFont font;
	font.setBold(true);
	QTabWidget* theTab=new QTabWidget;

	//----Getting Started (tab1)
	QWidget *tPg1 = new QWidget;
	// Labels
	QLabel* pg1Label=new QLabel(tr("Description Label\n")); pg1Label->setAlignment(Qt::AlignLeft);
	pg1Label->setWordWrap(true);
	// Design Tab Page Layout
	QGridLayout *pg1box = new QGridLayout;
	pg1box->setSizeConstraint(QLayout::SetFixedSize);
	pg1box->addWidget(pg1Label,0,0,1,2);
	// Set Layout
	tPg1->setLayout(pg1box);
	// Update Tab Widget
	theTab->addTab(tPg1,tr("Getting Started"));

	//-----protEVOLVER (tab2)
	QWidget *tPg2 = new QWidget;
	// Labels
	QLabel* pg2Label=new QLabel(tr("Evolve a sequence for a fold\n")); pg2Label->setAlignment(Qt::AlignLeft);
	pg2Label->setWordWrap(true);
	protEvolverPDBLabel=new QLabel;
	protEvolverPDBLabel->setText("Select PDB file:");
	protEvolverPDBLabel->setAlignment(Qt::AlignLeft);
	protEvolverActiveChainLabel=new QLabel;
	protEvolverActiveChainLabel->setText("Specify active chain(s):");
	protEvolverActiveChainLabel->setAlignment(Qt::AlignLeft);
	protEvolverActivePositionLabel=new QLabel;
	protEvolverActivePositionLabel->setText("Specify active position(s):");
	protEvolverActivePositionLabel->setAlignment(Qt::AlignLeft);
	protEvolverRandomPositionLabel=new QLabel;
	protEvolverRandomPositionLabel->setText("Specify random position(s):");
	protEvolverRandomPositionLabel->setAlignment(Qt::AlignLeft);
	protEvolverFrozenPositionLabel=new QLabel;
	protEvolverFrozenPositionLabel->setText("Specify frozen position(s):");
	protEvolverFrozenPositionLabel->setAlignment(Qt::AlignLeft);
	protEvolverAminoAcidLabel=new QLabel;
	protEvolverAminoAcidLabel->setText("Specify Amino Acids:");
	protEvolverAminoAcidLabel->setAlignment(Qt::AlignLeft);
	maxThreadsLabel=new QLabel;
	maxThreadsLabel->setText("Specify Number of Threads:");
	maxThreadsLabel->setAlignment(Qt::AlignLeft);
	// Line Input
	protEvolverActiveChainInput=new QLineEdit;
	protEvolverActiveChainInput->setAlignment(Qt::AlignCenter);
	protEvolverActiveChainInput->setFixedWidth(250);
	protEvolverActivePositionInput=new QLineEdit;
	protEvolverActivePositionInput->setAlignment(Qt::AlignCenter);
	protEvolverActivePositionInput->setFixedWidth(250);
	protEvolverRandomPositionInput=new QLineEdit;
	protEvolverRandomPositionInput->setAlignment(Qt::AlignCenter);
	protEvolverRandomPositionInput->setFixedWidth(250);
	protEvolverFrozenPositionInput=new QLineEdit;
	protEvolverFrozenPositionInput->setAlignment(Qt::AlignCenter);
	protEvolverFrozenPositionInput->setFixedWidth(250);
	protEvolverAminoAcidInput=new QLineEdit;
	protEvolverAminoAcidInput->setAlignment(Qt::AlignCenter);
	protEvolverAminoAcidInput->setFixedWidth(250);
	// Check Boxes
	protEvolverRelaxationBox=new QCheckBox(tr("Backbone relaxation?")); protEvolverRelaxationBox->setChecked(false); protEvolverRelaxationBox->setFont(font);
	protEvolverPolarityBox=new QCheckBox(tr("Auto assign polarity?")); protEvolverPolarityBox->setChecked(true); protEvolverPolarityBox->setFont(font);
	unsigned int number_of_threads = thread::hardware_concurrency();
	string numThreadValue;
	if(number_of_threads!=0)
		{numThreadValue=cnvrtNumToStrng(number_of_threads,0);}
	else
		{numThreadValue="1";}
	maxThreadsLine=new QLineEdit();
	maxThreadsLine->setAlignment(Qt::AlignCenter);
	maxThreadsLine->setText(tr(numThreadValue.c_str()));
	maxThreadsLine->setMaxLength(3);
	maxThreadsLine->setFixedWidth(80);
	// Buttons
	QPushButton* viewButton=new QPushButton(tr("View"));
	connect(viewButton,SIGNAL(clicked()),this,SLOT(view()));
	protEvolverPDBButton=new QPushButton(tr("..."));	
	protEvolverPDBButton->setToolTip(tr("Load start PDB to evolve a sequence for."));
	protEvolverPDBButton->setFixedWidth(80);
	protEvolverPDBButton->setFont(font);
	connect(protEvolverPDBButton,SIGNAL(clicked()),this,SLOT(open_protEvolverPDBFile()));
	xButton2=new QPushButton(tr("RUN"));	
	xButton2->setFixedWidth(200);
	xButton2->setFont(font);
	xButton2->setCheckable(true);
	xButton2->setChecked(false);
	connect(xButton2,SIGNAL(clicked()),this,SLOT(runProtEvolver()));
	// Design Tab Page Layout
	QGridLayout *pg2box = new QGridLayout;
	pg2box->setSizeConstraint(QLayout::SetFixedSize);
	pg2box->addWidget(pg2Label,0,0,1,3);
	pg2box->addWidget(protEvolverPDBLabel,1,0,1,1);pg2box->addWidget(protEvolverPDBButton,1,1,1,1);pg2box->addWidget(viewButton,1,2,1,1);
	pg2box->addWidget(protEvolverActiveChainLabel,2,0,1,1);pg2box->addWidget(protEvolverActiveChainInput,2,1,1,1);
	pg2box->addWidget(protEvolverActivePositionLabel,3,0,1,1);pg2box->addWidget(protEvolverActivePositionInput,3,1,1,1);
	pg2box->addWidget(protEvolverRandomPositionLabel,4,0,1,1);pg2box->addWidget(protEvolverRandomPositionInput,4,1,1,1);
	pg2box->addWidget(protEvolverFrozenPositionLabel,5,0,1,1);pg2box->addWidget(protEvolverFrozenPositionInput,5,1,1,1);
	pg2box->addWidget(protEvolverAminoAcidLabel,6,0,1,1);pg2box->addWidget(protEvolverAminoAcidInput,6,1,1,1);
	pg2box->addWidget(maxThreadsLabel,7,0,1,1);pg2box->addWidget(maxThreadsLine,7,1,1,1);
	pg2box->addWidget(protEvolverRelaxationBox,8,0,1,3);
	pg2box->addWidget(protEvolverPolarityBox,8,1,1,3);
	pg2box->addWidget(xButton2,9,0,1,3);
	// Set Layout
	tPg2->setLayout(pg2box);
	// Update Tab Widget
	theTab->addTab(tPg2,tr("protEvolver"));

	//---main layout
	QGridLayout *mainLayout = new QGridLayout;
	mainLayout->setSizeConstraint(QLayout::SetFixedSize);
	mainLayout->addWidget(theTab,0,0);
	setLayout(mainLayout);
	setWindowTitle(tr("protCAD"));
	resize(QDesktopWidget().availableGeometry(this).size());
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
		protEvolverPDBLabel->setText("Select PDB file:");
	}
}

void pUI::runProtEvolver()
{	
	if (xButton2->isChecked())
	{
		// set button as active
		xButton2->setChecked(true);
		xButton2->setText("STOP");

		//write input file from UI feilds
		string data=""; data+="PDB file,"; data+=protEvolver_pdbFile; data+=",\n";
		data+="Active Chains,"; data+=protEvolverActiveChainInput->text().toStdString(); data+="\n";
		data+="Active Positions,"; data+=protEvolverActivePositionInput->text().toStdString(); data+="\n";
		data+="Random Positions,"; data+=protEvolverRandomPositionInput->text().toStdString(); data+="\n";
		data+="Frozen Positions,"; data+=protEvolverFrozenPositionInput->text().toStdString(); data+="\n";
		data+="Amino Acids,"; data+=protEvolverAminoAcidInput->text().toStdString(); data+="\n";
		data+="Backbone Relaxation,"; if(protEvolverRelaxationBox->isChecked()){data+="true";} else{data+="false";} data+="\n";
		data+="Polarity Assignment,"; if(protEvolverPolarityBox->isChecked()){data+="true";} else{data+="false";} data+="\n";
		string sFldr=protEvolver_path;string inputFile=sFldr+"evolver.in";
		ofstream fOut; fOut.open(inputFile.c_str()); if(fOut.fail()){} fOut<<data;fOut.close();

		// define arguments and path for run
		QString cmd = QCoreApplication::applicationDirPath() +"/"+QString("protEvolver");
		QString workdir = QString::fromStdString(sFldr);
		QStringList args = {"evolver.in"};	
		string protcadpath="";
		vector<string> v = split (cmd.toStdString(), '/');
		for (unsigned int i = 0; i < v.size(); i++) 
		{
			if (i< v.size()-2){protcadpath+=v[i];protcadpath+="/";}
		}

		// run job using number of threads chosen
		string tmp=maxThreadsLine->text().toStdString();
		int nT=atoi(tmp.c_str());
		for(int i=0;i<nT;i++)
		{
			if (i == 1){usleep(500000);}
			qint64 pid;
			QProcessEnvironment env = QProcessEnvironment::systemEnvironment();
			env.insert("PROTCADDIR", QString::fromStdString(protcadpath));
			process->setProcessEnvironment(env);
			process->startDetached(cmd, args, workdir, &pid);
		}
	}
	else
	{
		xButton2->setChecked(false);xButton2->setText("RUN");
		QString stop = "killall protEvolver";
		process->start(stop);
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