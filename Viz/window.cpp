
#include <QPainter>
#include <QMessageBox>
#include <stdio.h>
#include <math.h>
#include "window.h"
#define DEFAULT_A -1
#define DEFAULT_B 1
#define DEFAULT_N 5
#define INDENT 20
#define MAX_N_for_Chebyshev 50
#define MINIMUM_DISTANCE 5.e-15
#define LAGRANGE_CONST 5
#define MESSAGE_LEN 100
double function (double x);
double function (double x)
{
  return x;
}



Window::Window (QWidget *parent, std::string file_name_g, std::string file_name_v, int npoints, int nlayers, double delta_x, double delta_t,
                double min_g, double min_v, double max_g, double max_v)
  : QWidget (parent)
{
  m_file_name_g = file_name_g;
  m_file_name_v = file_name_v;

  m_fp_v = fopen (m_file_name_v.c_str(), "r");
  m_fp_g = fopen (m_file_name_g.c_str(), "r");

  m_npoints = npoints;
  m_nlayers = nlayers;
  m_points_array.resize (m_npoints, 0);
  m_min_g = min_g;
  m_min_v = min_v;
  m_max_g = max_g;
  m_max_v = max_v;

  m_delta_x = delta_x;
  m_delta_t = delta_t;
  m_current_function_is_g = true;
  m_current_layer = 0;

  m_timer = new QTimer (this);
  m_timer->setInterval (10);
  connect (m_timer, SIGNAL(timeout()), this, SLOT(update_layer ()));
//  v_g ();
}

void Window::start ()
{
  restart_file();
  restart ();
}

void Window::v_g ()
{
  m_current_function_is_g = !m_current_function_is_g;
  start ();
}

void Window::slower ()
{
  if (m_timer->interval () == 0)
    m_timer->setInterval (10);
  else
    m_timer->setInterval (m_timer->interval () * 2);
}

void Window::faster ()
{
  m_timer->setInterval (m_timer->interval () / 2);
}

void Window::update_layer ()
{
  if (m_current_layer >= m_nlayers)
    {
      m_timer->stop ();
      return;
    }

  if (m_current_function_is_g)
    read_layer (m_fp_g);
  else
    read_layer (m_fp_v);
  m_current_layer++;
  update ();
}

Window::~Window ()
{

}

QSize Window::minimumSizeHint () const
{
  return QSize (100, 100);
}

QSize Window::sizeHint () const
{
  return QSize (1000, 1000);
}

void Window::restart ()
{
  m_current_layer = 0;
  m_timer->start ();
}

void Window::restart_file ()
{
  int tmp1, tmp2;
  double tmp3, tmp4;
  if (m_current_function_is_g)
    {
      fclose (m_fp_g);
      m_fp_g = fopen (m_file_name_g.c_str(), "r");
      if (fscanf (m_fp_g, "%d%d%lf%lf", &tmp1, &tmp2, &tmp3, &tmp4) != 4)
        return;
    }
  else
    {
      fclose (m_fp_v);
      m_fp_v = fopen (m_file_name_v.c_str(), "r");
      if (fscanf (m_fp_v, "%d%d%lf%lf", &tmp1, &tmp2, &tmp3, &tmp4) != 4)
        return;
    }
}

void Window::read_layer (FILE *fp)
{
  for (int i = 0; i < m_npoints; i++)
    if (fscanf (fp, "%lf", &m_points_array[i]) != 1)
      return;
}

void Window::paintEvent (QPaintEvent * /* event */)
{
  QPainter painter (this);
  double max_y;
  double min_y;
  painter.setRenderHint(QPainter::Qt4CompatiblePainting);
  if (m_current_function_is_g)
    {
      max_y = m_max_g;
      min_y = m_min_g;
      painter.setPen ("green");
    }
  else
    {
      max_y = m_max_v;
      min_y = m_min_v;
      painter.setPen ("red");
    }

  // save current Coordinate System
  painter.save ();

  // make Coordinate Transformations
  painter.translate (0.5 * width (), 0.5 * height ());
  painter.scale (width () / ((m_npoints - 1) * m_delta_x), -height () / (max_y - min_y));
  painter.translate (-0.5 * ((m_npoints - 1) * m_delta_x), -0.5 * (min_y + max_y));

  for (int i = 0; i < m_npoints - 1; i++)
    {
      painter.drawLine (QPointF (i * m_delta_x, m_points_array[i]), QPointF ((i + 1) * m_delta_x, m_points_array[i + 1]));
    }

  // draw axis
  painter.setPen ("black");
  painter.drawLine (QPointF (0, 0), QPointF ((m_npoints - 1) * m_delta_x, 0));
  painter.drawLine (QPointF (0, max_y), QPointF (0, min_y));

  // restore previously saved Coordinate System
  painter.restore ();
  char str[MESSAGE_LEN];
  painter.setPen ("black");
  if (m_current_function_is_g)
    sprintf (str, "g (x, %1.3g)", m_current_layer * m_delta_t);
  else
    sprintf (str, "v (x, %1.3g)", m_current_layer * m_delta_t);
  painter.drawText (0, INDENT, str);

}


























