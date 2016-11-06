
#ifndef WINDOW_H
#define WINDOW_H

#include <QTimer>

//#define NUMBER_OF_POINTS  2000
#include <QWidget>
class Window : public QWidget
{
  Q_OBJECT

private:
  int m_npoints;
  int m_nlayers;
  int m_current_layer;
  FILE *m_fp_g;
  FILE *m_fp_v;
  std::string m_file_name_g;
  std::string m_file_name_v;
  double m_max_x;
  double m_max_g;
  double m_min_g;
  double m_max_v;
  double m_min_v;
  std::vector<double> m_points_array;
  bool m_current_function_is_g;
  double m_delta_x;
  double m_delta_t;
  QTimer *m_timer = nullptr;

public:
  Window (QWidget *parent, std::string file_name_g, std::string file_name_v, int npoints, int nlayers, double delta_x,
          double delta_t, double min_g, double min_v, double max_g, double max_v);
  ~Window ();
  QSize minimumSizeHint () const;
  QSize sizeHint () const;
  void restart ();
  void restart_file ();
  void read_layer (FILE *fp);

public slots:
 void start ();
 void v_g ();
 void slower ();
 void faster ();
 void update_layer ();
protected:
  void paintEvent (QPaintEvent *event);
};

#endif
