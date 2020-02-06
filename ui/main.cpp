#include <QApplication>

#include "ui.h"

int main(int argc, char *argv[])
{
	QApplication app(argc, argv);
	pUI G;
	G.show();
	return app.exec();
}
