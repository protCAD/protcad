#include <QApplication>
#include <QSplashScreen>
#include <QTimer>
#include "ui.h"

int main(int argc, char *argv[])
{
	QApplication app(argc, argv);
	QPixmap pixmap(":/splash");
	QSplashScreen splash(pixmap, Qt::WindowStaysOnTopHint);
	splash.show();
	QTimer::singleShot(2500, &splash, SLOT(close()));
	pUI G;
	QTimer::singleShot(2500, &G, SLOT(show()));
	return app.exec();
}
