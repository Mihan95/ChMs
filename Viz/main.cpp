#include <QApplication>
#include <QMainWindow>
#include <QVBoxLayout>
#include <QAction>
#include <QMenuBar>
#include <QMessageBox>
#include <QIcon>
#include <QLocale>

#include "window.h"

static
void find_min_max (FILE *fp, double &min, double &max, int &npoints, int &nlayers, double &delta_t, double &delta_x)
{
  double tmp;
  if (fscanf (fp, "%d %d %lf %lf", &npoints, &nlayers, &delta_x, &delta_t) != 4)
    printf ("xxxxx\n");

  if (fscanf (fp, "%lf", &max) != 1)
    return;
  min = max;
  for (int i = 0; i < npoints * nlayers - 1; i++)
    {
      if (fscanf (fp, "%lf", &tmp) != 1)
        return;
      if (tmp > max)
        max = tmp;
      if (tmp < min)
        min = tmp;
    }
}

int main (int argc, char *argv[])
{
  QApplication app (argc, argv);

  QLocale curLocale(QLocale("en_US"));
  QLocale::setDefault(curLocale);

  std::string filename_v = "v.txt";
  std::string filename_g = "h.txt";
  FILE *fp_v = fopen (filename_v.c_str(), "r");
  FILE *fp_g = fopen (filename_g.c_str(), "r");
  int npoints, nlayers;
  double delta_t, delta_x, max_v, max_g, min_v, min_g;
  find_min_max (fp_v, min_v, max_v, npoints, nlayers, delta_t, delta_x);
  find_min_max (fp_g, min_g, max_g, npoints, nlayers, delta_t, delta_x);
  fclose (fp_v);
  fclose (fp_g);


  QMainWindow *window = new QMainWindow;
  QMenuBar *tool_bar = new QMenuBar (window);
  Window *graph_area = new Window (window, filename_g, filename_v, npoints, nlayers, delta_x, delta_t, min_g, min_v, max_g, max_v);
  QAction *action;

  action = tool_bar->addAction ("G/G_origin", graph_area, SLOT (v_g ()));
  action->setShortcut (QString ("1"));

  action = tool_bar->addAction ("Start", graph_area, SLOT (start ()));
  action->setShortcut (QString ("2"));

  action = tool_bar->addAction ("Faster", graph_area, SLOT (faster ()));
  action->setShortcut (QString ("3"));

  action = tool_bar->addAction ("Slower", graph_area, SLOT (slower ()));
  action->setShortcut (QString ("4"));

  tool_bar->setMaximumHeight (30);

  window->setMenuBar (tool_bar);
  window->setCentralWidget (graph_area);
  window->setWindowTitle ("Graph");
  window->show ();
  int ret = app.exec ();
  delete tool_bar;
  delete graph_area;
  delete window;
  return ret;
}
