using System;
using System.Diagnostics;
using System.Reflection;
using System.Runtime.InteropServices;
using System.Windows.Forms;
using DiffUrSolver;
using OxyPlot;
using OxyPlot.Annotations;
using OxyPlot.Axes;
using OxyPlot.Legends;
using OxyPlot.Series;
using OxyPlot.WindowsForms;
using HorizontalAlignment = OxyPlot.HorizontalAlignment;

namespace WinFormsApp
{
    public partial class Form2 : Form
    {
        public Form2()
        {
            InitializeComponent();

            example1(Methods.RadoIIA,100);

            //example1(Methods.BackEuler, 100);

            //example1(Methods.RK4, 0.000001); // С таким шагом - считает.

            //example2();

            //example3(Methods.RK4);

            //example5(0.024, Methods.BackEuler, 0.0001, Methods.RK4);

            //example4(0.025, Methods.Euler);

            //example4(0.0010, Methods.RK4);

            //example4(0.0005, Methods.RK4);

            //example5(0.0009, Methods.Euler, 0.0001, Methods.RK4);

            //example5(0.3, Methods.BackEuler, 0.0001, Methods.RK4);

            //example6(0.001, Methods.RK4);

            //example6(0.0001, Methods.Euler, 0.001, Methods.RK4);

            //example6(0.01, Methods.Euler, 0.01, Methods.RK4);
        }


        public void example1(Methods method, double h)
        {
            DiffUrSolver.Task task = new DiffUrSolver.Task();

            task.x0 = 0.0;
            task.y0 = RealSampleFunction6(task.x0);
            task.x1 = 100.0;


            task.f = DiffUrSolver.DiffUrSolver.F1ToIntPtr(SampleFunction6);


            task.h = h;
            task.method = method;

            //task.h = 2;
            //task.method = Methods.BackEuler;

            //string title = "y'  = -50(y-cos(x))";

            string title = "y'  = -1000000 (y - sin(x)) + cos(x)";

            DrawPlot(task, RealSampleFunction6, title);
        }
        
        public void example2()
        {
            DiffUrSolver.SystemTask systemTask = new DiffUrSolver.SystemTask();

            systemTask.t0 = 0.0;
            systemTask.t1 = 4.0;
            systemTask.h = 0.1;
            systemTask.n = 2;
            systemTask.method = Methods.RK4;

            double[] y0 = new double[2];
            y0[0] = 1;
            y0[1] = 1;

            double SystemFunction1(double x, IntPtr y)
        {

            double[] array = new double[2];
            Marshal.Copy(y, array, 0, 2);


            return array[0] - 5 * array[1];
        }

            double SystemFunction2(double x, IntPtr y)
            {
                double[] array = new double[2];
                Marshal.Copy(y, array, 0, 2);

                return 5 * array[0] + array[1];
            }
            
            double RealSystemFunction1(double x)
            {
                return Math.Pow(Math.E, x) * (Math.Sin(-5 * x) + Math.Cos(-5 * x));
            }

            double RealSystemFunction2(double x)
        {
            return Math.Pow(Math.E, x) * (Math.Sin(5 * x) + Math.Cos(5 * x));
        }

            systemTask.y0 = DiffUrSolver.DiffUrSolver.DoubleArrayToIntPtr(y0);
            systemTask.f = DiffUrSolver.DiffUrSolver.F1SystemToIntPtr(new F1System[] { SystemFunction1, SystemFunction2 });

            DrawSystemPlot(systemTask, new F[] { RealSystemFunction1, RealSystemFunction2 });
        }
        
        public void example3(Methods method)
        {
            double c = 4.47;
            double m = 10.0;
            double k = 50.0;

            SecondTask task = new SecondTask();
            task.t0 = 0;
            task.t1 = 10;
            task.y_0 = 0.01;
            task.y1_0 = 0;
            task.h = 0.1;
            task.A = 1;
            task.B = c / m;
            task.C = k / m;
            task.D = 0;
            task.method = method;
            double[] result = new double[10000];

            double dzeta = c / (2 * Math.Sqrt(k * m));

            double w = Math.Sqrt(k / m);

            double wd = w * Math.Sqrt(1 - dzeta * dzeta);

            double C1 = task.y_0;

            double C2 = (task.y1_0 + C1 * dzeta * w) / wd;

            double RealFunction(double t)
            {
                return Math.Exp(-dzeta * w * t) * (C1 * Math.Cos(wd * t) + C2 * Math.Sin(wd*t));
            }

            string title = task.A + "*y'' + " + $"{c}/{m}" + "*y' + " + $"{k}/{m}" + "*y = " + task.D; 

            DrawSecondPlot(task,RealFunction,  title);
        }

        public void example4(double h, Methods method)
        {
            DiffUrSolver.SystemTask systemTask = new DiffUrSolver.SystemTask();

            systemTask.t0 = 0.0;
            systemTask.t1 = 0.3;
            systemTask.h = h;
            systemTask.n = 3;
            systemTask.method = method;

            double[] y0 = new double[3];
            y0[0] = 1;
            y0[1] = 0;
            y0[1] = 0;

            double SystemFunction1(double x, IntPtr y)
            {

                double[] array = new double[3];
                Marshal.Copy(y, array, 0, 3);


                return -0.04 * array[0] + 10000 * array[1] * array[2];
            }

            double SystemFunction2(double x, IntPtr y)
            {
                double[] array = new double[3];
                Marshal.Copy(y, array, 0, 3);

                return 0.04 * array[0] - 10000 * array[1] * array[2] - 30000000 * array[1] * array[1];
            }

            double SystemFunction3(double x, IntPtr y)
            {
                double[] array = new double[3];
                Marshal.Copy(y, array, 0, 3);

                return 30000000 * array[1] * array[1];
            }

            systemTask.y0 = DiffUrSolver.DiffUrSolver.DoubleArrayToIntPtr(y0);
            systemTask.f = DiffUrSolver.DiffUrSolver.F1SystemToIntPtr(new F1System[] { SystemFunction1, SystemFunction2, SystemFunction3 });

            DrawSystemPlot(systemTask);
        }

        public void example5(double h, Methods method, double h2, Methods method2)
        {
            DiffUrSolver.DiffUrSolver.SetThreadCount(2);
            DiffUrSolver.SystemTask systemTask = new DiffUrSolver.SystemTask();
            DiffUrSolver.SystemTask systemTask2 = new DiffUrSolver.SystemTask();

            systemTask.t0 = 0.0;
            systemTask.t1 = 0.3;
            systemTask.h = h;
            systemTask.n = 3;
            systemTask.method = method;


            systemTask2.t0 = 0.0;
            systemTask2.t1 = 0.3;
            systemTask2.h = h2;
            systemTask2.n = 3;
            systemTask2.method = method2;

            double[] y0 = new double[3];
            y0[0] = 1;
            y0[1] = 0;
            y0[1] = 0;

            double SystemFunction1(double x, IntPtr y)
            {

                double[] array = new double[3];
                Marshal.Copy(y, array, 0, 3);


                return -0.04 * array[0] + 10000 * array[1] * array[2];
            }

            double SystemFunction2(double x, IntPtr y)
            {
                double[] array = new double[3];
                Marshal.Copy(y, array, 0, 3);

                return 0.04 * array[0] - 10000 * array[1] * array[2] - 30000000 * array[1] * array[1];
            }

            double SystemFunction3(double x, IntPtr y)
            {
                double[] array = new double[3];
                Marshal.Copy(y, array, 0, 3);

                return 30000000 * array[1] * array[1];
            }

            systemTask.y0 = DiffUrSolver.DiffUrSolver.DoubleArrayToIntPtr(y0);
            systemTask.f = DiffUrSolver.DiffUrSolver.F1SystemToIntPtr(new F1System[] { SystemFunction1, SystemFunction2, SystemFunction3 });

            systemTask2.y0 = DiffUrSolver.DiffUrSolver.DoubleArrayToIntPtr(y0);
            systemTask2.f = DiffUrSolver.DiffUrSolver.F1SystemToIntPtr(new F1System[] { SystemFunction1, SystemFunction2, SystemFunction3 });

            DrawSystemPlotTwo(systemTask,systemTask2);
        }


        struct BrusselatorSystem
        {
            public IntPtr f;
            public IntPtr y0;
            public double t0;
            public int status;
        }
        static BrusselatorSystem CreateBrusselatorSystem(int Nx,int Ny, double a)
        {

            //BrusselatorSystem res = new BrusselatorSystem();

            //if (n < 1)
            //{

            //    res.status = -1;
            //    return res;
            //}

            //F1System[] system = new F1System[n * 2];
            //double[] y = new double[2 * n];

            //for (int i = 0; i < n; i++)
            //{

            //    int i1 = (n + i - 1) % n;
            //    int i2 = i;
            //    int i3 = (n + i + 1) % n;

            //    F1System u = (double x, IntPtr y) =>
            //    {
            //        double[] array = new double[2 * n];
            //        Marshal.Copy(y, array, 0, 2 * n);

            //        return 1 + array[i2 * 2] * array[i2 * 2] * array[i2 * 2 + 1] -
            //        4 * array[i2 * 2] + a * (n + 1) * (n + 1) *
            //        (array[i1 * 2] - 2 * array[i2 * 2] + array[i3 * 2]);
            //    };

            //    F1System v = (double x, IntPtr y) =>
            //    {
            //        double[] array = new double[2 * n];
            //        Marshal.Copy(y, array, 0, 2 * n);

            //        return 3 * array[i2 * 2] - array[i2 * 2] * array[i2 * 2] * array[i2 * 2 + 1] +
            //        a * (n + 1) * (n + 1) *
            //        (array[i1 * 2 + 1] - 2 * array[i2 * 2 + 1] + array[i3 * 2 + 1]);
            //    };

            //    system[i * 2] = u;
            //    system[i * 2 + 1] = v;

            //    y[i * 2] = 1 + Math.Sin(2 * Math.PI * i / n);
            //    y[i * 2 + 1] = 3;

            //}

            //res.status = 0;
            //res.f = DiffUrSolver.DiffUrSolver.F1SystemToIntPtr(system);
            //res.y0 = DiffUrSolver.DiffUrSolver.DoubleArrayToIntPtr(y);
            //res.t0 = 0;


            //return res;

            BrusselatorSystem res = new BrusselatorSystem();

            

            int N = Nx * Ny;
            int TotalSize = 2 * N;

            F1System[] system = new F1System[TotalSize];
            double[] y = new double[TotalSize];

            double A = 1.0;
            double B = 3.0;

            double Du = 1e-5;
            double Dv = 1e-1;

            double hx = 1.0 / (Nx - 1);
            double hy = 1.0 / (Ny - 1);

            double hx2 = hx * hx;
            double hy2 = hy * hy;

            int U(int i, int j)
            {
                return i + j * Nx;
            }

            int V(int i, int j)
            {
                return N + i + j * Ny;
            }

            //---------------------------------------------------------
            // Границы Неймана
            //---------------------------------------------------------

            int ClampX(int i)
            {
                if (i < 0) return 0;
                if (i >= Nx) return Nx - 1;
                return i;
            }

            int ClampY(int j)
            {
                if (j < 0) return 0;
                if (j >= Ny) return Ny - 1;
                return j;
            }

            //---------------------------------------------------------
            // Лапласиан
            //---------------------------------------------------------

            double Laplace(double[] y, int idx, int i, int j)
            {
                int il = ClampX(i - 1);
                int ir = ClampX(i + 1);

                int jb = ClampY(j - 1);
                int jt = ClampY(j + 1);

                double center = y[idx];

                double left =
                    y[idx - (i - il)];

                double right =
                    y[idx + (ir - i)];

                double down =
                    y[idx - (j - jb) * Nx];

                double up =
                    y[idx + (jt - j) * Nx];

                return
                    (left - 2.0 * center + right) / hx2 +
                    (down - 2.0 * center + up) / hy2;
            }

            //---------------------------------------------------------
            // Генератор массива функций
            //---------------------------------------------------------

            Random rnd = new Random();

            //---------------------------------------------
            // Уравнения для u
            //---------------------------------------------

            for (int j = 0; j < Ny; j++)
            {
                for (int i = 0; i < Nx; i++)
                {
                    int ii = i;
                    int jj = j;

                    int idx = U(ii, jj);

                    //---------------------------------------------
                    // Начальные условия
                    //---------------------------------------------

                    y[idx] =
                        1.0 +
                        0.01 * (rnd.NextDouble() - 0.5);

                    //---------------------------------------------
                    // Функция
                    //---------------------------------------------

                    

                    system[idx] = (double t, IntPtr y) =>
                    {

                        double[] array = new double[TotalSize];
                        Marshal.Copy(y, array, 0, TotalSize);

                        int iu = U(ii, jj);
                        int iv = V(ii, jj);

                        double u = array[iu];
                        double v = array[iv];

                        double Lu = Laplace(array, iu, ii, jj);

                        return
                            Du * Lu
                            + A
                            - (B + 1.0) * u
                            + u * u * v;
                    };
                }
            }

            //---------------------------------------------
            // Уравнения для v
            //---------------------------------------------

            for (int j = 0; j < Ny; j++)
            {
                for (int i = 0; i < Nx; i++)
                {
                    int ii = i;
                    int jj = j;

                    int idx = V(ii, jj);

                    //---------------------------------------------
                    // Начальные условия
                    //---------------------------------------------

                    y[idx] =
                        3.0 +
                        0.01 * (rnd.NextDouble() - 0.5);

                    //---------------------------------------------
                    // Функция
                    //---------------------------------------------
                    
                    system[idx] = (double t, IntPtr y) =>
                    {

                        double[] array = new double[TotalSize];
                        Marshal.Copy(y, array, 0, TotalSize);

                        int iu = U(ii, jj);
                        int iv = V(ii, jj);

                        double u = array[iu];
                        double v = array[iv];

                        double Lv = Laplace(array, iv, ii, jj);

                        return
                            Dv * Lv
                            + B * u
                            - u * u * v;
                    };
                }
            }

            res.status = 0;
            res.f = DiffUrSolver.DiffUrSolver.F1SystemToIntPtr(system);
            res.y0 = DiffUrSolver.DiffUrSolver.DoubleArrayToIntPtr(y);
            res.t0 = 0;


            return res;



        }



        public void example6(double h, Methods method)
        {
            DiffUrSolver.SystemTask systemTask = new DiffUrSolver.SystemTask();

            uint n = 30;

            BrusselatorSystem b = CreateBrusselatorSystem(((int)n), ((int)n), 0.02);


            systemTask.y0 = b.y0;
            systemTask.f = b.f;
            systemTask.t0 = b.t0;
            systemTask.t1 = 10;
            systemTask.h = h;
            systemTask.n = 2 * n * n;
            systemTask.method = method;

            DrawSystemPlot(systemTask);
        }

        public void example6(double h, Methods method, double h2, Methods method2)
        {
            uint n = 30;

            DiffUrSolver.DiffUrSolver.SetThreadCount(4);
            DiffUrSolver.DiffUrSolver.InitMemory(2 * n * n);

            DiffUrSolver.SystemTask systemTask1 = new DiffUrSolver.SystemTask();
            DiffUrSolver.SystemTask systemTask2 = new DiffUrSolver.SystemTask();

            

            BrusselatorSystem b = CreateBrusselatorSystem(((int)n), ((int)n), 0.02);


            systemTask1.y0 = b.y0;
            systemTask1.f = b.f;
            systemTask1.t0 = b.t0;
            systemTask1.t1 = 6;
            systemTask1.h = h;
            systemTask1.n = 2 * n * n;
            systemTask1.method = method;

            systemTask2.y0 = b.y0;
            systemTask2.f = b.f;
            systemTask2.t0 = b.t0;
            systemTask2.t1 = 6;
            systemTask2.h = h2;
            systemTask2.n = 2 * n * n;
            systemTask2.method = method2;

            DrawSystemPlotTwo(systemTask1, systemTask2);

        }

        public static double SampleFunction(double x, double y)
        {
            return y;
        }

        public static double RealSampleFunction(double x)
        {
            return Math.Pow(Math.E,x);
        }

        public static double SampleFunction2(double x, double y)
        {
            return Math.Cos(x);
        }

        public static double RealSampleFunction2(double x)
        {
            return Math.Sin(x);
        }

        public static double SampleFunction3(double x, double y)
        {
            return Math.Cos(x);
        }

        public static double RealSampleFunction3(double x)
        {
            return Math.Sin(x);
        }

        public static double SampleFunction4(double x, double y)
        {
            return -50 *(y -Math.Cos(x));
        }

        public static double RealSampleFunction4(double x)
        {
            return Math.Pow(Math.E, -50 *x) + 50 * Math.Sin(x)/ 2501 + 2500 * Math.Cos(x) / 2501;
        }

        public static double SampleFunction5(double x, double y)
        {
            return 5 * Math.Cos(6*x + 30);
        }

        public static double RealSampleFunction5(double x)
        {
            return 5 * Math.Sin(6 * x + 30)/6;
        }

        public static double SampleFunction6(double x, double y)
        {
            return -1000000 * (y - Math.Sin(x)) + Math.Cos(x);
        }

        public static double RealSampleFunction6(double x)
        {
            return Math.Sin(x);
        }

        private const AnnotationLayer HiddenLayer = (AnnotationLayer)100;
        private const AnnotationLayer VisibleLayer = AnnotationLayer.AboveSeries;

        private class EpsAnnotation : TextAnnotation { };

        private class MethodAnnotation : TextAnnotation { };

        private void DrawPlot(DiffUrSolver.Task task, DiffUrSolver.F realFunc,string title = "")
        {
            var plotModel = new PlotModel { Title = title };

            var lineSeries = new LineSeries
            {
                Title = "Method solution",
                MarkerType = MarkerType.Circle,
                MarkerSize = 3,
                MarkerStroke = OxyColors.Red,
                MarkerFill = OxyColors.Red,
                Color = OxyColors.Red,

            };

            var lineSeries2 = new LineSeries
            {
                Title = "Exact solution",
                MarkerType = MarkerType.Circle,
                MarkerSize = 3,
                MarkerStroke = OxyColors.Green,
                MarkerFill = OxyColors.Green,
                Color = OxyColors.Green,
            };

            int Size = (int)((task.x1 - task.x0) / task.h) + 10;

            double[] resultX = new double[Size];
            double[] resultY = new double[Size];

            int status = DiffUrSolver.DiffUrSolver.SolveDiffUrArr(task, resultY);

            resultX[0] = task.x0;
            int i = 1;

            double v = task.x0 + i * task.h;

            while ( v < task.x1)
            {
                resultX[i] = v;
                i++;
                v = task.x0 + i * task.h;
            }

            resultX[i] = task.x1;

            int resultSize = i + 1;

            lineSeries.Points.Add(new DataPoint(resultX[0], resultY[0]));
            lineSeries2.Points.Add(new DataPoint(resultX[0], realFunc(resultX[0])));
            double diff = Math.Abs(resultY[0] - realFunc(resultX[0]));
            plotModel.Annotations.Add(new EpsAnnotation
            {
                Text = diff.ToString("F6"),
                TextPosition = new DataPoint(resultX[0], (resultY[0] + realFunc(resultX[0])) / 2),
                Background = OxyColors.White,
                TextColor = OxyColors.Black,
                Layer = HiddenLayer

            });

            for (i = 1; i < resultSize; i++)
            {
                lineSeries.Points.Add(new DataPoint(resultX[i], resultY[i]));
                lineSeries2.Points.Add(new DataPoint(resultX[i], realFunc(resultX[i])));
                diff = Math.Abs(resultY[i] - realFunc(resultX[i]));
                plotModel.Annotations.Add(new EpsAnnotation
                {
                    Text = diff.ToString("F6"),
                    TextPosition = new DataPoint(resultX[i], (resultY[i] + realFunc(resultX[i])) / 2),
                    Background = OxyColors.White,
                    TextColor = OxyColors.Black,
                    Layer = HiddenLayer

                });

                plotModel.Annotations.Add(new MethodAnnotation
                {
                    Text = (resultX[i] - resultX[i - 1]).ToString(),
                    TextPosition = new DataPoint((resultX[i] + resultX[i - 1]) / 2, (resultY[i] + resultY[i - 1]) / 2),
                    Background = OxyColors.White,
                    TextColor = OxyColors.Black,
                    Layer = HiddenLayer

                });

            }


            plotModel.Series.Add(lineSeries);
            plotModel.Series.Add(lineSeries2);

            Legend l = new Legend();

            l.LegendPosition = LegendPosition.RightTop;
            l.LegendPlacement = LegendPlacement.Outside;
            l.LegendOrientation = LegendOrientation.Vertical;

            plotModel.Legends.Add(l);

            var plotView = new PlotView
            {
                Dock = DockStyle.Fill,
                Model = plotModel
            };

            this.Controls.Add(plotView);

            this.KeyPreview = true;
            this.KeyDown += (sender, e) =>
            {
                if (e.Control && e.KeyCode == Keys.E)
                {
                    foreach (var annotation in plotModel.Annotations.OfType<EpsAnnotation>())
                    {
                        annotation.Layer = annotation.Layer == HiddenLayer ? VisibleLayer : HiddenLayer;
                    }
                    plotModel.InvalidatePlot(true);
                }

                if (e.Control && e.KeyCode == Keys.M)
                {
                    foreach (var annotation in plotModel.Annotations.OfType<MethodAnnotation>())
                    {
                        annotation.Layer = annotation.Layer == HiddenLayer ? VisibleLayer : HiddenLayer;
                    }
                    plotModel.InvalidatePlot(true);
                }
            };

        }

        private void DrawSecondPlot(DiffUrSolver.SecondTask task, DiffUrSolver.F realFunc, string Title = "")
        {
            var plotModel = new PlotModel { Title = Title };

            var lineSeries = new LineSeries
            {
                Title = "Method solution",
                MarkerType = MarkerType.Circle,
                MarkerSize = 3,
                MarkerStroke = OxyColors.Red,
                MarkerFill = OxyColors.Red,
                Color = OxyColors.Red,

            };

            var lineSeries2 = new LineSeries
            {
                Title = "Exact solution",
                MarkerType = MarkerType.Circle,
                MarkerSize = 3,
                MarkerStroke = OxyColors.Green,
                MarkerFill = OxyColors.Green,
                Color = OxyColors.Green,
            };

            double[] resultX = new double[10000];
            double[] resultY = new double[10000];


            int status =  DiffUrSolver.DiffUrSolver.SolveSecondDiffUrArr(task, ref resultY);
           
            resultX[0] = task.t0;
            int i = 1;

            double v = task.t0 + i * task.h;

            while (v < task.t1)
            {
                resultX[i] = v;
                i++;
                v = task.t0 + i * task.h;
            }

            resultX[i] = task.t1;

            int resultSize = i + 1;

            lineSeries.Points.Add(new DataPoint(resultX[0], resultY[0]));
            lineSeries2.Points.Add(new DataPoint(resultX[0], realFunc(resultX[0])));
            double diff = Math.Abs(resultY[0] - realFunc(resultX[0]));
            plotModel.Annotations.Add(new EpsAnnotation
            {
                Text = diff.ToString("F6"),
                TextPosition = new DataPoint(resultX[0], (resultY[0] + realFunc(resultX[0])) / 2),
                Background = OxyColors.White,
                TextColor = OxyColors.Black,
                Layer = HiddenLayer

            });

            for (i = 1; i < resultSize; i++)
            {
                lineSeries.Points.Add(new DataPoint(resultX[i], resultY[i*2]));
                lineSeries2.Points.Add(new DataPoint(resultX[i], realFunc(resultX[i])));
                diff = Math.Abs(resultY[i] - realFunc(resultX[i]));
                plotModel.Annotations.Add(new EpsAnnotation
                {
                    Text = diff.ToString("F6"),
                    TextPosition = new DataPoint(resultX[i], (resultY[i*2] + realFunc(resultX[i])) / 2),
                    Background = OxyColors.White,
                    TextColor = OxyColors.Black,
                    Layer = HiddenLayer

                });

                plotModel.Annotations.Add(new MethodAnnotation
                {
                    Text = (resultX[i] - resultX[i - 1]).ToString(),
                    TextPosition = new DataPoint((resultX[i] + resultX[i - 1]) / 2, (resultY[i*2] + resultY[i*2 - 2]) / 2),
                    Background = OxyColors.White,
                    TextColor = OxyColors.Black,
                    Layer = HiddenLayer

                });

            }


            plotModel.Series.Add(lineSeries);
            plotModel.Series.Add(lineSeries2);

            

            plotModel.Axes.Add(new LinearAxis
            {
                Position = AxisPosition.Bottom,
                Title = "t",
                //MajorGridlineStyle = LineStyle.Solid,

                ExtraGridlines = new[] { 0.0 },
                ExtraGridlineStyle = LineStyle.Solid,
                ExtraGridlineThickness = 2,
                ExtraGridlineColor = OxyColors.Black
            });

            plotModel.Axes.Add(new LinearAxis
            {
                Position = AxisPosition.Left,
                Title = "y",
                //MajorGridlineStyle = LineStyle.Solid,

                ExtraGridlines = new[] { 0.0 },
                ExtraGridlineStyle = LineStyle.Solid,
                ExtraGridlineThickness = 2,
                ExtraGridlineColor = OxyColors.Black
            });

            plotModel.IsLegendVisible = true;

            Legend l = new Legend();

            l.LegendPosition = LegendPosition.RightTop;
            l.LegendPlacement = LegendPlacement.Outside;
            l.LegendOrientation = LegendOrientation.Vertical;

            plotModel.Legends.Add(l);

            var plotView = new PlotView
            {
                Dock = DockStyle.Fill,
                Model = plotModel
            };

            this.Controls.Add(plotView);

            this.KeyPreview = true;
            this.KeyDown += (sender, e) =>
            {
                if (e.Control && e.KeyCode == Keys.E)
                {
                    foreach (var annotation in plotModel.Annotations.OfType<EpsAnnotation>())
                    {
                        annotation.Layer = annotation.Layer == HiddenLayer ? VisibleLayer : HiddenLayer;
                    }
                    plotModel.InvalidatePlot(true);
                }

                if (e.Control && e.KeyCode == Keys.M)
                {
                    foreach (var annotation in plotModel.Annotations.OfType<MethodAnnotation>())
                    {
                        annotation.Layer = annotation.Layer == HiddenLayer ? VisibleLayer : HiddenLayer;
                    }
                    plotModel.InvalidatePlot(true);
                }
            };

        }
        private void DrawSystemPlot(DiffUrSolver.SystemTask systemTask, F[] RealSystem)
        {

            TabControl tabControl = new TabControl();
            tabControl.Dock = DockStyle.Fill;
            this.Controls.Add(tabControl);

            int size = 10000000;
            double[] resultX = new double[size];
            double[] resultY = new double[size * systemTask.n];

            int status = DiffUrSolver.DiffUrSolver.SolveDiffUrSystemArr(systemTask, resultY);

            int k = 0;
            while (systemTask.t0 + k * systemTask.h < systemTask.t1)
            {
                resultX[k] = systemTask.t0 + k * systemTask.h;
                k++;
            }

            resultX[k] = systemTask.t1;

            int resultSize = k;



            for (int i = 0; i < systemTask.n; i++)
            {
                var plotModel = new PlotModel { Title = "Пример простого графика" };

                var MylineSeries = new LineSeries
                {
                    Title = "My",
                    MarkerType = MarkerType.Circle,
                    MarkerSize = 3,
                    MarkerStroke = OxyColors.Red,
                    MarkerFill = OxyColors.Red,
                    Color = OxyColors.Red,

                };

                var ReallineSeries = new LineSeries
                {
                    Title = "Real",
                    MarkerType = MarkerType.Circle,
                    MarkerSize = 3,
                    MarkerStroke = OxyColors.Green,
                    MarkerFill = OxyColors.Green,
                    Color = OxyColors.Green,
                };

                MylineSeries.Points.Add(new DataPoint(resultX[0], resultY[i]));
                ReallineSeries.Points.Add(new DataPoint(resultX[0], RealSystem[i](resultX[0])));
                double diff = Math.Abs(resultY[i] - RealSystem[i](resultX[0]));

                plotModel.Annotations.Add(new EpsAnnotation
                {
                    Text = diff.ToString("F6"),
                    TextPosition = new DataPoint(resultX[0], (resultY[i] + RealSystem[i](resultX[0])) / 2),
                    Background = OxyColors.White,
                    TextColor = OxyColors.Black,
                    Layer = HiddenLayer

                });

                for (int j = 1; j < resultSize; j++)
                {
                    MylineSeries.Points.Add(new DataPoint(resultX[j], resultY[i + j* systemTask.n]));
                    ReallineSeries.Points.Add(new DataPoint(resultX[j], RealSystem[i](resultX[j])));
                    diff = Math.Abs(resultY[i + j * systemTask.n] - RealSystem[i](resultX[j]));

                    plotModel.Annotations.Add(new EpsAnnotation
                    {
                        Text = diff.ToString("F6"),
                        TextPosition = new DataPoint(resultX[j], (resultY[i + j * systemTask.n] + RealSystem[i](resultX[j])) / 2),
                        Background = OxyColors.White,
                        TextColor = OxyColors.Black,
                        Layer = HiddenLayer

                    });

                    plotModel.Annotations.Add(new MethodAnnotation
                    {
                        Text = (resultX[j] - resultX[j - 1]).ToString(),
                        TextPosition = new DataPoint((resultX[j] + resultX[j - 1]) / 2, (resultY[i + j * systemTask.n] + resultY[i + j * (systemTask.n-1)]) / 2),
                        Background = OxyColors.White,
                        TextColor = OxyColors.Black,
                        Layer = HiddenLayer

                    });
                }


                plotModel.Series.Add(MylineSeries);
                plotModel.Series.Add(ReallineSeries);

                var plotView = new PlotView
                {
                    Dock = DockStyle.Fill,
                    Model = plotModel
                };

                var tabPage = new TabPage("Ctranica " + i );

                tabPage.Controls.Add(plotView);
                tabControl.TabPages.Add(tabPage);

                this.KeyPreview = true;
                this.KeyDown += (sender, e) =>
                {
                    if (e.Control && e.KeyCode == Keys.E)
                    {
                        foreach (var annotation in plotModel.Annotations.OfType<EpsAnnotation>())
                        {
                            annotation.Layer = annotation.Layer == HiddenLayer ? VisibleLayer : HiddenLayer;
                        }
                        plotModel.InvalidatePlot(true);
                    }

                    if (e.Control && e.KeyCode == Keys.M)
                    {
                        foreach (var annotation in plotModel.Annotations.OfType<MethodAnnotation>())
                        {
                            annotation.Layer = annotation.Layer == HiddenLayer ? VisibleLayer : HiddenLayer;
                        }
                        plotModel.InvalidatePlot(true);
                    }
                };

            }

           /* double xMax = lineSeries.Points.Max(p => Math.Abs(p.X));
            double yMax = lineSeries.Points.Max(p => Math.Abs(p.Y));
            double maxRange = Math.Max(xMax, yMax) * 1.1; // +10% для отступов


            plotModel.Axes.Add(new LinearAxis
            {
                Position = AxisPosition.Bottom,
                Minimum = -maxRange,  // Минимальное значение
                Maximum = maxRange,   // Максимальное значение
                MajorGridlineStyle = LineStyle.Solid,
                MinorGridlineStyle = LineStyle.Dot,
                ExtraGridlines = new double[] { 0 },  // Линия X=0
                ExtraGridlineColor = OxyColors.Black,
                ExtraGridlineThickness = 2,
                Title = "Y0",
                TitlePosition = 1,
            });

            // Ось Y (вертикальная)
            plotModel.Axes.Add(new LinearAxis
            {
                Position = AxisPosition.Left,
                Minimum = -maxRange,  // Минимальное значение (должно совпадать с осью X!)
                Maximum = maxRange,   // Максимальное значение (должно совпадать с осью X!)
                MajorGridlineStyle = LineStyle.Solid,
                MinorGridlineStyle = LineStyle.Dot,
                ExtraGridlines = new double[] { 0 },  // Линия Y=0
                ExtraGridlineColor = OxyColors.Black,
                ExtraGridlineThickness = 2,
                Title = "Y1",
                TitlePosition = 1,

            });

            var plotView = new PlotView
            {
                Dock = DockStyle.Fill,
                Model = plotModel
            };

            this.Controls.Add(plotView); */

        }

        private void DrawSystemPlot(DiffUrSolver.SystemTask systemTask)
        {

            TabControl tabControl = new TabControl();
            tabControl.Dock = DockStyle.Fill;
            this.Controls.Add(tabControl);

            int size = (int)Math.Round((systemTask.t1 - systemTask.t0) / systemTask.h) + 10;

            double[] resultX = new double[size];
            double[] resultY = new double[size * systemTask.n];

            int status = DiffUrSolver.DiffUrSolver.SolveDiffUrSystemArr(systemTask, resultY);

            int k = 0;
            while (systemTask.t0 + k * systemTask.h < systemTask.t1)
            {
                resultX[k] = systemTask.t0 + k * systemTask.h;
                k++;
            }

            resultX[k] = systemTask.t1;

            int resultSize = k;



            for (int i = 0; i < systemTask.n; i++)
            {
                var plotModel = new PlotModel { Title = "" };

                var MylineSeries = new LineSeries
                {
                    Title = "My",
                    //MarkerType = MarkerType.Circle,
                    MarkerSize = 3,
                    MarkerStroke = OxyColors.Red,
                    MarkerFill = OxyColors.Red,
                    Color = OxyColors.Red,

                };

                MylineSeries.Points.Add(new DataPoint(resultX[0], resultY[i]));

                for (int j = 1; j < resultSize; j++)
                {
                    MylineSeries.Points.Add(new DataPoint(resultX[j], resultY[i + j * systemTask.n]));

                    plotModel.Annotations.Add(new MethodAnnotation
                    {
                        Text = (resultX[j] - resultX[j - 1]).ToString(),
                        TextPosition = new DataPoint((resultX[j] + resultX[j - 1]) / 2, (resultY[i + j * systemTask.n] + resultY[i + j * (systemTask.n - 1)]) / 2),
                        Background = OxyColors.White,
                        TextColor = OxyColors.Black,
                        Layer = HiddenLayer

                    });
                }




                plotModel.Series.Add(MylineSeries);

                var plotView = new PlotView
                {
                    Dock = DockStyle.Fill,
                    Model = plotModel
                };

                var tabPage = new TabPage("Ctranica " + i);

                tabPage.Controls.Add(plotView);
                tabControl.TabPages.Add(tabPage);

                //this.KeyPreview = true;
                //this.KeyDown += (sender, e) =>
                //{
                //    if (e.Control && e.KeyCode == Keys.E)
                //    {
                //        foreach (var annotation in plotModel.Annotations.OfType<EpsAnnotation>())
                //        {
                //            annotation.Layer = annotation.Layer == HiddenLayer ? VisibleLayer : HiddenLayer;
                //        }
                //        plotModel.InvalidatePlot(true);
                //    }

                //    if (e.Control && e.KeyCode == Keys.M)
                //    {
                //        foreach (var annotation in plotModel.Annotations.OfType<MethodAnnotation>())
                //        {
                //            annotation.Layer = annotation.Layer == HiddenLayer ? VisibleLayer : HiddenLayer;
                //        }
                //        plotModel.InvalidatePlot(true);
                //    }
                //};

            }


        }

        private void DrawSystemPlotTwo(DiffUrSolver.SystemTask systemTask, DiffUrSolver.SystemTask systemTask2)
        {

            TabControl tabControl = new TabControl();
            tabControl.Dock = DockStyle.Fill;
            this.Controls.Add(tabControl);

            int size1 = (int)Math.Round((systemTask.t1 - systemTask.t0) / systemTask.h) + 10;
            double[] resultX = new double[size1];
            double[] resultY = new double[size1 * systemTask.n];

            Stopwatch stopwatch = new Stopwatch();

            

            stopwatch.Start();

            int status = DiffUrSolver.DiffUrSolver.SolveDiffUrSystemArr(systemTask, resultY);

            stopwatch.Stop();

            int size2 = (int)Math.Round((systemTask2.t1 - systemTask2.t0) / systemTask2.h) + 10;

            double[] resultX2 = new double[size2];
            double[] resultY2 = new double[size2 * systemTask.n];

            int status2 = DiffUrSolver.DiffUrSolver.SolveDiffUrSystemArr(systemTask2, resultY2);

            int k = 0;
            while (systemTask.t0 + k * systemTask.h < systemTask.t1)
            {
                resultX[k] = systemTask.t0 + k * systemTask.h;
                k++;
            }

            resultX[k] = systemTask.t1;

            int resultSize = k+1;

            k = 0;
            while (systemTask2.t0 + k * systemTask2.h < systemTask2.t1)
            {
                resultX2[k] = systemTask2.t0 + k * systemTask2.h;
                k++;
            }

            resultX2[k] = systemTask2.t1;

            int resultSize2 = k+1;



            for (int i = 0; i < systemTask.n; i++)
            {
                var plotModel = new PlotModel { Title = "" };

                var MylineSeries = new LineSeries
                {
                    Title = "Method Solution",
                    //MarkerType = MarkerType.Circle,
                    MarkerSize = 3,
                    MarkerStroke = OxyColors.Red,
                    MarkerFill = OxyColors.Red,
                    Color = OxyColors.Red,

                };

                var ReallineSeries = new LineSeries
                {
                    Title = "nearly exact solution",
                    //MarkerType = MarkerType.Circle,
                    MarkerSize = 3,
                    MarkerStroke = OxyColors.Green,
                    MarkerFill = OxyColors.Green,
                    Color = OxyColors.Green,

                };

                MylineSeries.Points.Add(new DataPoint(resultX[0], resultY[i]));

                for (int j = 1; j < resultSize; j++)
                {
                    MylineSeries.Points.Add(new DataPoint(resultX[j], resultY[i + j * systemTask.n]));

                    plotModel.Annotations.Add(new MethodAnnotation
                    {
                        Text = (resultX[j] - resultX[j - 1]).ToString(),
                        TextPosition = new DataPoint((resultX[j] + resultX[j - 1]) / 2, (resultY[i + j * systemTask.n] + resultY[i + j * (systemTask.n - 1)]) / 2),
                        Background = OxyColors.White,
                        TextColor = OxyColors.Black,
                        Layer = HiddenLayer

                    });
                }

                ReallineSeries.Points.Add(new DataPoint(resultX2[0], resultY2[i]));

                for (int j = 1; j < resultSize2; j++)
                {
                    ReallineSeries.Points.Add(new DataPoint(resultX2[j], resultY2[i + j * systemTask2.n]));

                    plotModel.Annotations.Add(new MethodAnnotation
                    {
                        Text = (resultX2[j] - resultX2[j - 1]).ToString(),
                        TextPosition = new DataPoint((resultX2[j] + resultX2[j - 1]) / 2, (resultY2[i + j * systemTask2.n] + resultY2[i + j * (systemTask2.n - 1)]) / 2),
                        Background = OxyColors.White,
                        TextColor = OxyColors.Black,
                        Layer = HiddenLayer

                    });
                }

                plotModel.IsLegendVisible = true;

                Legend l = new Legend();

                l.LegendPosition = LegendPosition.RightTop;
                l.LegendPlacement = LegendPlacement.Outside;
                l.LegendOrientation = LegendOrientation.Vertical;
                l.LegendTitle = "time: " + stopwatch.ElapsedMilliseconds + " ms";

                plotModel.Legends.Add(l);

                


                plotModel.Series.Add(MylineSeries);
                plotModel.Series.Add(ReallineSeries);

                var plotView = new PlotView
                {
                    Dock = DockStyle.Fill,
                    Model = plotModel
                };

                var tabPage = new TabPage("Ctranica " + i);

                tabPage.Controls.Add(plotView);
                tabControl.TabPages.Add(tabPage);

                //this.KeyPreview = true;
                //this.KeyDown += (sender, e) =>
                //{
                //    if (e.Control && e.KeyCode == Keys.E)
                //    {
                //        foreach (var annotation in plotModel.Annotations.OfType<EpsAnnotation>())
                //        {
                //            annotation.Layer = annotation.Layer == HiddenLayer ? VisibleLayer : HiddenLayer;
                //        }
                //        plotModel.InvalidatePlot(true);
                //    }

                //    if (e.Control && e.KeyCode == Keys.M)
                //    {
                //        foreach (var annotation in plotModel.Annotations.OfType<MethodAnnotation>())
                //        {
                //            annotation.Layer = annotation.Layer == HiddenLayer ? VisibleLayer : HiddenLayer;
                //        }
                //        plotModel.InvalidatePlot(true);
                //    }
                //};

            }


        }

        //private void DrawAutoHPlot(DiffUrSolver.F1 func, DiffUrSolver.F realFunc, DiffUrSolver.Methods method)
        //{
        //    var plotModel = new PlotModel { Title = "Пример простого графика" };

        //    var lineSeries = new LineSeries
        //    {
        //        Title = "My",
        //        MarkerType = MarkerType.Circle,
        //        MarkerSize = 3,
        //        MarkerStroke = OxyColors.Red,
        //        MarkerFill = OxyColors.Red,
        //        Color = OxyColors.Red,

        //    };

        //    var lineSeries2 = new LineSeries
        //    {
        //        Title = "Real",
        //        MarkerType = MarkerType.Circle,
        //        MarkerSize = 3,
        //        MarkerStroke = OxyColors.Green,
        //        MarkerFill = OxyColors.Green,
        //        Color = OxyColors.Green,
        //    };

        //    double x0 = 0.0;
        //    double y0 = 1.0;
        //    double x1 = 4.0;
        //    double h = 0.0001;
        //    double eps = 0.0001;
        //    double[] resultX = new double[10000];
        //    double[] resultY = new double[10000];
        //    int resultSize = 0;



        //    int status = DiffUrSolver.DiffUrSolver.SolveDiffUrAutoHArr(x0, y0, x1, h, eps, func,
        //        resultX, resultY, ref resultSize, method);


        //    lineSeries.Points.Add(new DataPoint(resultX[0], resultY[0]));
        //    lineSeries2.Points.Add(new DataPoint(resultX[0], realFunc(resultX[0])));
        //    double diff = Math.Abs(resultY[0] - realFunc(resultX[0]));
        //    plotModel.Annotations.Add(new EpsAnnotation
        //    {
        //        Text = diff.ToString("F6"),
        //        TextPosition = new DataPoint(resultX[0], (resultY[0] + realFunc(resultX[0])) / 2),
        //        Background = OxyColors.White,
        //        TextColor = OxyColors.Black,
        //        Layer = HiddenLayer

        //    });

        //    for (int i = 1; i < resultSize; i++)
        //    {
        //        lineSeries.Points.Add(new DataPoint(resultX[i], resultY[i]));
        //        lineSeries2.Points.Add(new DataPoint(resultX[i], realFunc(resultX[i])));
        //        diff = Math.Abs(resultY[i] - realFunc(resultX[i]));
        //        plotModel.Annotations.Add(new EpsAnnotation
        //        {
        //            Text = diff.ToString("F6"),
        //            TextPosition = new DataPoint(resultX[i], (resultY[i] + realFunc(resultX[i])) / 2),
        //            Background = OxyColors.White,
        //            TextColor = OxyColors.Black,
        //            Layer = HiddenLayer

        //        });

        //        plotModel.Annotations.Add(new MethodAnnotation
        //        {
        //            Text = (resultX[i] - resultX[i - 1]).ToString(),
        //            TextPosition = new DataPoint((resultX[i] + resultX[i - 1]) / 2, (resultY[i] + resultY[i - 1]) / 2),
        //            Background = OxyColors.White,
        //            TextColor = OxyColors.Black,
        //            Layer = HiddenLayer

        //        });

        //    }


        //    plotModel.Series.Add(lineSeries);
        //    plotModel.Series.Add(lineSeries2);

        //    var plotView = new PlotView
        //    {
        //        Dock = DockStyle.Fill,
        //        Model = plotModel
        //    };

        //    this.Controls.Add(plotView);

        //    this.KeyPreview = true;
        //    this.KeyDown += (sender, e) =>
        //    {
        //        if (e.Control && e.KeyCode == Keys.E)
        //        {
        //            foreach (var annotation in plotModel.Annotations.OfType<EpsAnnotation>())
        //            {
        //                annotation.Layer = annotation.Layer == HiddenLayer ? VisibleLayer : HiddenLayer;
        //            }
        //            plotModel.InvalidatePlot(true);
        //        }

        //        if (e.Control && e.KeyCode == Keys.M)
        //        {
        //            foreach (var annotation in plotModel.Annotations.OfType<MethodAnnotation>())
        //            {
        //                annotation.Layer = annotation.Layer == HiddenLayer ? VisibleLayer : HiddenLayer;
        //            }
        //            plotModel.InvalidatePlot(true);
        //        }
        //    };

        //}

        //private void DrawAutoMethodPlot(DiffUrSolver.F1 func, DiffUrSolver.F realFunc, DiffUrSolver.Methods method)
        //{
        //    var plotModel = new PlotModel { Title = "Пример простого графика" };

        //    var lineSeries = new LineSeries
        //    {
        //        Title = "My",
        //        MarkerType = MarkerType.Circle,
        //        MarkerSize = 3,
        //        MarkerStroke = OxyColors.Red,
        //        MarkerFill = OxyColors.Red,

        //    };

        //    var lineSeries2 = new LineSeries
        //    {
        //        Title = "Real",
        //        MarkerType = MarkerType.Circle,
        //        MarkerSize = 3,
        //        MarkerStroke = OxyColors.Green,
        //        MarkerFill = OxyColors.Green,
        //    };

        //    double x0 = 0.0;
        //    double y0 = 0.0;
        //    double x1 = 6.0;
        //    double h = 0.05;
        //    double eps = 0.001;
        //    double[] resultX = new double[1000];
        //    double[] resultY = new double[1000];
        //    Methods[] methodsArr = new Methods[1000];
        //    int resultSize = 0;



        //    int status = DiffUrSolver.DiffUrSolver.SolveDiffUrAutoMethodArr2(x0, y0, x1, h, eps, func, resultY, ref resultSize, methodsArr, method);


        //    lineSeries.Points.Add(new DataPoint(x0, resultY[0]));
        //    lineSeries2.Points.Add(new DataPoint(x0, realFunc(x0)));
        //    double diff = Math.Abs(resultY[0] - realFunc(x0));
        //    plotModel.Annotations.Add(new EpsAnnotation
        //    {
        //        Text = diff.ToString("F6"),
        //        TextPosition = new DataPoint(x0, (resultY[0] + realFunc(x0)) / 2),
        //        Background = OxyColors.White,
        //        TextColor = OxyColors.Black,
        //        Layer = HiddenLayer

        //    });

        //    for (int i = 1; i < resultSize; i++)
        //    {
        //        lineSeries.Points.Add(new DataPoint(x0 + i * h, resultY[i]));
        //        lineSeries2.Points.Add(new DataPoint(x0 + i * h, realFunc(x0 + i * h)));
        //        diff = Math.Abs(resultY[i] - realFunc(x0 + i * h));
        //        plotModel.Annotations.Add(new EpsAnnotation
        //        {
        //            Text = diff.ToString("F6"),
        //            TextPosition = new DataPoint(x0 + i * h, (resultY[i] + realFunc(x0 + i * h)) / 2),
        //            Background = OxyColors.White,
        //            TextColor = OxyColors.Black,
        //            Layer = HiddenLayer

        //        });



        //        plotModel.Annotations.Add(new MethodAnnotation
        //        {
        //            Text = methodsArr[i].ToString(),
        //            TextPosition = new DataPoint((x0 + i * h + x0 + (i - 1) * h) / 2, (resultY[i] + resultY[i - 1]) / 2),
        //            Background = OxyColors.White,
        //            TextColor = OxyColors.Black,
        //            Layer = HiddenLayer

        //        });

        //    }


        //    plotModel.Series.Add(lineSeries);
        //    plotModel.Series.Add(lineSeries2);

        //    var plotView = new PlotView
        //    {
        //        Dock = DockStyle.Fill,
        //        Model = plotModel
        //    };

        //    this.Controls.Add(plotView);

        //    this.KeyPreview = true;
        //    this.KeyDown += (sender, e) =>
        //    {
        //        if (e.Control && e.KeyCode == Keys.E)
        //        {
        //            foreach (var annotation in plotModel.Annotations.OfType<EpsAnnotation>())
        //            {
        //                annotation.Layer = annotation.Layer == HiddenLayer ? VisibleLayer : HiddenLayer;
        //            }
        //            plotModel.InvalidatePlot(true);
        //        }

        //        if (e.Control && e.KeyCode == Keys.M)
        //        {
        //            foreach (var annotation in plotModel.Annotations.OfType<MethodAnnotation>())
        //            {
        //                annotation.Layer = annotation.Layer == HiddenLayer ? VisibleLayer : HiddenLayer;
        //            }
        //            plotModel.InvalidatePlot(true);
        //        }
        //    };

        //}
    }
}