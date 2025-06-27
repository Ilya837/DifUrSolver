using System.Runtime.InteropServices;
using System.Windows.Forms;
using DiffUrSolver;
using OxyPlot;
using OxyPlot.Annotations;
using OxyPlot.Axes;
using OxyPlot.Series;
using OxyPlot.WindowsForms;

namespace WinFormsApp
{
    public partial class Form2 : Form
    {
        public Form2()
        {
            InitializeComponent();
            //DrawPlot(SampleFunction, RealSampleFunction, DiffUrSolver.Methods.RK4);
            //DrawAutoHPlot(SampleFunction, RealSampleFunction,DiffUrSolver.Methods.RK2);
            //DrawAutoMethodPlot(SampleFunction2, RealSampleFunction2, DiffUrSolver.Methods.Euler);

            F1System[] system = new F1System[2];

            system[0] = SystemFunction1;
            system[1] = SystemFunction2;

            F[] Realsystem = new F[2];

            Realsystem[0] = RealSystemFunction1;
            Realsystem[1] = RealSystemFunction2;

            DrawSystemPlot(system,Realsystem,Methods.RK2);
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

        public static double SystemFunction1(double x, IntPtr y)
        {

            double[] array = new double[2];
            Marshal.Copy(y, array, 0, 2);


            return array[0] - 5 * array[1];
        }

        public static double SystemFunction2(double x, IntPtr y)
        {
            double[] array = new double[2];
            Marshal.Copy(y, array, 0, 2);

            return 5 *array[0] + array[1];
        }
        public static double RealSystemFunction1(double x)
        {
            return Math.Pow(Math.E, x) * (Math.Sin(-5 * x) + Math.Cos(-5 * x));
        }

        public static double RealSystemFunction2(double x)
        {
            return Math.Pow(Math.E, x) * (Math.Sin(5 * x) + Math.Cos(5 * x));
        }

        private const AnnotationLayer HiddenLayer = (AnnotationLayer)100;
        private const AnnotationLayer VisibleLayer = AnnotationLayer.AboveSeries;

        private class EpsAnnotation : TextAnnotation { };

        private class MethodAnnotation : TextAnnotation { };

        private void DrawPlot(DiffUrSolver.F1 func, DiffUrSolver.F realFunc, DiffUrSolver.Methods method)
        {
            var plotModel = new PlotModel { Title = "Пример простого графика" };

            var lineSeries = new LineSeries
            {
                Title = "My",
                MarkerType = MarkerType.Circle,
                MarkerSize = 3,
                MarkerStroke = OxyColors.Red,
                MarkerFill = OxyColors.Red,
                Color = OxyColors.Red,

            };

            var lineSeries2 = new LineSeries
            {
                Title = "Real",
                MarkerType = MarkerType.Circle,
                MarkerSize = 3,
                MarkerStroke = OxyColors.Green,
                MarkerFill = OxyColors.Green,
                Color = OxyColors.Green,
            };

            double x0 = 0.0;
            double y0 = 1.0;
            double x1 = 4.0;
            double h = 0.1;
            double[] resultX = new double[10000];
            double[] resultY = new double[10000];
            int resultSize = 0;



            int status = DiffUrSolver.DiffUrSolver.SolveDiffUrArr(x0, y0, x1, h, func,
                 resultY, ref resultSize, method);

            for (int i = 0; i < resultSize - 1; i++)
            {
                resultX[i] = x0 + i * h;
            }

            resultX[resultSize - 1] = x1;


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

            for (int i = 1; i < resultSize; i++)
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
        private void DrawSystemPlot(F1System[] system, F[] RealSystem, DiffUrSolver.Methods method)
        {
            var plotModel = new PlotModel { Title = "Пример простого графика" };

            

            var lineSeries = new LineSeries
            {
                Title = "My",
                MarkerType = MarkerType.Circle,
                MarkerSize = 3,
                MarkerStroke = OxyColors.Red,
                MarkerFill = OxyColors.Red,
                Color = OxyColors.Red,

            };

            var lineSeries2 = new LineSeries
            {
                Title = "Real",
                MarkerType = MarkerType.Circle,
                MarkerSize = 3,
                MarkerStroke = OxyColors.Green,
                MarkerFill = OxyColors.Green,
                Color = OxyColors.Green,
            };

            double x0 = 0.0;
            double[] y0 = new double[6];
            y0[0] = 1;
            y0[1] = 1;
            double x1 = 4.0;
            double h = 0.001;
            int size = 10000000;
            uint n = 2;
            double[] resultX = new double[size];
            double[] resultY = new double[size*n];



            int resultSize = 0;

            // Создаём массив IntPtr для хранения указателей на делегаты
            IntPtr[] functionPointers = new IntPtr[system.Length];

            // Фиксируем делегаты в памяти, чтобы GC не перемещал их
            GCHandle[] handles = new GCHandle[system.Length];

            for (int i = 0; i < system.Length; i++)
            {
                handles[i] = GCHandle.Alloc(system[i]);
                functionPointers[i] = Marshal.GetFunctionPointerForDelegate(system[i]);
            }

            GCHandle functionsArrayHandle = GCHandle.Alloc(functionPointers, GCHandleType.Pinned);

            int status = DiffUrSolver.DiffUrSolver.SolveDiffUrSystemArr(x0, y0, x1,n, h, functionsArrayHandle.AddrOfPinnedObject(),
                 resultY,ref resultSize, method);

            int k = 0;
            while(x0 +k*h + 0.0000001 < x1)
            {
                resultX[k] = x0 + k * h;
                k++;
            }

            resultX[k] = x1;


            //lineSeries.Points.Add(new DataPoint(resultY[0], resultY[1]));
            //lineSeries2.Points.Add(new DataPoint(RealSystemFunction1(resultX[0]), RealSystemFunction2(resultX[0])));
            //double diff = Math.Abs(resultY[0] - realFunc(resultX[0]));
            //plotModel.Annotations.Add(new EpsAnnotation
            //{
            //    Text = diff.ToString("F6"),
            //    TextPosition = new DataPoint(resultX[0], (resultY[0] + realFunc(resultX[0])) / 2),
            //    Background = OxyColors.White,
            //    TextColor = OxyColors.Black,
            //    Layer = HiddenLayer

            //});

            double diff = 0;

            for (int i = 0; i < resultSize; i++)
            {
                lineSeries.Points.Add(new DataPoint(resultY[i*n], resultY[i*n + 1]));
                lineSeries2.Points.Add(new DataPoint(RealSystem[0](resultX[i]), RealSystem[1](resultX[i])));

                diff = Math.Sqrt( Math.Pow(resultY[i * n] - RealSystem[0](resultX[i]),2) + Math.Pow(resultY[i * n + 1] - RealSystem[1](resultX[i]), 2));
                plotModel.Annotations.Add(new EpsAnnotation
                {
                    Text = diff.ToString("F6"),
                    TextPosition = new DataPoint( (resultY[i * n] + RealSystem[0](resultX[i])) / 2, (resultY[i * n + 1] + RealSystem[1](resultX[i])) / 2),
                    Background = OxyColors.White,
                    TextColor = OxyColors.Black,
                    Layer = HiddenLayer

                });

                //plotModel.Annotations.Add(new MethodAnnotation
                //{
                //    Text = (resultX[i] - resultX[i - 1]).ToString(),
                //    TextPosition = new DataPoint((resultX[i] + resultX[i - 1]) / 2, (resultY[i] + resultY[i - 1]) / 2),
                //    Background = OxyColors.White,
                //    TextColor = OxyColors.Black,
                //    Layer = HiddenLayer

                //});

            }

            

            plotModel.Series.Add(lineSeries);
            plotModel.Series.Add(lineSeries2);

            double xMax = lineSeries.Points.Max(p => Math.Abs(p.X));
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
        private void DrawAutoHPlot(DiffUrSolver.F1 func, DiffUrSolver.F realFunc,DiffUrSolver.Methods method)
        {
            var plotModel = new PlotModel { Title = "Пример простого графика" };

            var lineSeries = new LineSeries
            {
                Title = "My",
                MarkerType = MarkerType.Circle,
                MarkerSize = 3,
                MarkerStroke = OxyColors.Red,
                MarkerFill = OxyColors.Red,
                Color = OxyColors.Red,
                
            };

            var lineSeries2 = new LineSeries
            {
                Title = "Real",
                MarkerType = MarkerType.Circle,
                MarkerSize = 3,
                MarkerStroke = OxyColors.Green,
                MarkerFill= OxyColors.Green,
                Color = OxyColors.Green,
            };

            double x0 = 0.0;
            double y0 = 1.0;
            double x1 = 4.0;
            double h = 0.0001;
            double eps = 0.0001;
            double[] resultX = new double[10000];
            double[] resultY = new double[10000];
            int resultSize = 0;



            int status = DiffUrSolver.DiffUrSolver.SolveDiffUrAutoHArr(x0, y0, x1, h, eps, func,
                resultX, resultY, ref resultSize, method);


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

            for (int i =1; i< resultSize; i++)
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
                    Text = (resultX[i] - resultX[i-1]).ToString(),
                    TextPosition = new DataPoint((resultX[i] + resultX[i - 1])/2, (resultY[i] + resultY[i - 1]) / 2),
                    Background = OxyColors.White,
                    TextColor = OxyColors.Black,
                    Layer = HiddenLayer

                });

            }


            plotModel.Series.Add(lineSeries);
            plotModel.Series.Add(lineSeries2);

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

        private void DrawAutoMethodPlot(DiffUrSolver.F1 func, DiffUrSolver.F realFunc, DiffUrSolver.Methods method)
        {
            var plotModel = new PlotModel { Title = "Пример простого графика" };

            var lineSeries = new LineSeries
            {
                Title = "My",
                MarkerType = MarkerType.Circle,
                MarkerSize = 3,
                MarkerStroke = OxyColors.Red,
                MarkerFill = OxyColors.Red,

            };

            var lineSeries2 = new LineSeries
            {
                Title = "Real",
                MarkerType = MarkerType.Circle,
                MarkerSize = 3,
                MarkerStroke = OxyColors.Green,
                MarkerFill = OxyColors.Green,
            };

            double x0 = 0.0;
            double y0 = 0.0;
            double x1 = 6.0;
            double h = 0.05;
            double eps = 0.001;
            double[] resultX = new double[1000];
            double[] resultY = new double[1000];
            Methods[] methodsArr = new Methods[1000];
            int resultSize = 0;



            int status = DiffUrSolver.DiffUrSolver.SolveDiffUrAutoMethodArr2(x0, y0, x1, h, eps, func,resultY, ref resultSize, methodsArr, method);


            lineSeries.Points.Add(new DataPoint(x0 , resultY[0]));
            lineSeries2.Points.Add(new DataPoint(x0, realFunc(x0)));
            double diff = Math.Abs(resultY[0] - realFunc(x0));
            plotModel.Annotations.Add(new EpsAnnotation
            {
                Text = diff.ToString("F6"),
                TextPosition = new DataPoint(x0, (resultY[0] + realFunc(x0)) / 2),
                Background = OxyColors.White,
                TextColor = OxyColors.Black,
                Layer = HiddenLayer

            });

            for (int i = 1; i < resultSize; i++)
            {
                lineSeries.Points.Add(new DataPoint(x0 + i * h, resultY[i]));
                lineSeries2.Points.Add(new DataPoint(x0 + i * h, realFunc(x0 + i * h)));
                diff = Math.Abs(resultY[i] - realFunc(x0 + i * h));
                plotModel.Annotations.Add(new EpsAnnotation
                {
                    Text = diff.ToString("F6"),
                    TextPosition = new DataPoint(x0 + i * h, (resultY[i] + realFunc(x0 + i * h)) / 2),
                    Background = OxyColors.White,
                    TextColor = OxyColors.Black,
                    Layer = HiddenLayer

                });

                

                plotModel.Annotations.Add(new MethodAnnotation
                {
                    Text = methodsArr[i].ToString(),
                    TextPosition = new DataPoint((x0 + i * h + x0 + (i-1) * h)/2, (resultY[i] + resultY[i-1]) / 2),
                    Background = OxyColors.White,
                    TextColor = OxyColors.Black,
                    Layer = HiddenLayer

                });
                
            }


            plotModel.Series.Add(lineSeries);
            plotModel.Series.Add(lineSeries2);

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
    }
}